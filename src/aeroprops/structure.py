# SPDX-License-Identifier: MIT
"""Structure input layer: RDKit loading, normalisation, and composition.

This module is the only place that touches RDKit. RDKit is a required
dependency of aeroprops, but it is still imported lazily (attempted once
and cached) so that ``import aeroprops`` itself stays lightweight; public
functions additionally return ``None`` rather than raise if RDKit is
missing from a broken installation.
"""

_rdkit_modules = None
_rdkit_loaded = False
_normalizer = None


def load_rdkit():
    """Load the RDKit modules used by aeroprops.

    Returns
    -------
    modules : tuple or None
        ``(Chem, Descriptors)`` if RDKit is importable, otherwise ``None``.
        The import is attempted once per process and cached.

    Examples
    --------
    >>> load_rdkit() is not None
    True
    """
    global _rdkit_modules, _rdkit_loaded
    if not _rdkit_loaded:
        _rdkit_loaded = True
        try:
            from rdkit import Chem as _Chem
            from rdkit.Chem import Descriptors as _Descriptors
        except ImportError:
            _rdkit_modules = None
        else:
            _rdkit_modules = (_Chem, _Descriptors)
    return _rdkit_modules


def _normalizer_transform():
    """Return (and cache) an RDKit MolStandardize ``Normalizer`` instance."""
    global _normalizer
    if _normalizer is None:
        from rdkit.Chem.MolStandardize import rdMolStandardize

        _normalizer = rdMolStandardize.Normalizer()
    return _normalizer


def normalize_mol(mol):
    """Return `mol` passed through RDKit's ``Normalizer`` transform.

    Every structure entering aeroprops passes through this charge-form
    correction (see :func:`mol_from_smiles`): RDKit renders nitrate,
    nitro, and peroxyacyl nitrate (PAN) groups in the charge-separated
    form (e.g. ``C[N+](=O)[O-]``), and the aeroprops SMARTS catalogues
    are written against that form. Callers handing in a pre-built
    ``Mol`` object are responsible for normalising it with this function
    before fragmenting.

    Parameters
    ----------
    mol : Mol
        RDKit molecule (already parsed).

    Returns
    -------
    mol : Mol or None
        Normalised molecule, or ``None`` if `mol` is ``None`` or RDKit
        is unavailable.

    Examples
    --------
    >>> from rdkit import Chem
    >>> m = normalize_mol(Chem.MolFromSmiles("CN(=O)=O"))
    >>> Chem.MolToSmiles(m)
    'C[N+](=O)[O-]'
    >>> normalize_mol(None) is None
    True
    """
    if mol is None or load_rdkit() is None:
        return None
    return _normalizer_transform().normalize(mol)


def mol_from_smiles(smi):
    """Parse a SMILES string into a normalised RDKit molecule.

    Equivalent to parsing with ``Chem.MolFromSmiles`` followed by
    :func:`normalize_mol`.

    Parameters
    ----------
    smi : str
        SMILES string; the valence (neutral) form need not be normalised
        on input.

    Returns
    -------
    mol : Mol or None
        Normalised molecule, or ``None`` if RDKit is unavailable or the
        SMILES cannot be parsed.

    Examples
    --------
    >>> mol = mol_from_smiles("CCO")
    >>> mol.GetNumAtoms()
    3
    >>> mol_from_smiles("not a smiles") is None
    True
    """
    loaded = load_rdkit()
    if loaded is None:
        return None
    mol = loaded[0].MolFromSmiles(smi)
    if mol is None:
        return None
    return _normalizer_transform().normalize(mol)


def canonical_smiles(smi):
    """Return the canonical SMILES of the normalised molecule for `smi`.

    The two valence writings of a nitrate/nitro group map to the same
    canonical (charge-separated) form:

    Examples
    --------
    >>> canonical_smiles("CN(=O)=O") == canonical_smiles("C[N+](=O)[O-]")
    True
    """
    loaded = load_rdkit()
    mol = mol_from_smiles(smi)
    if mol is None or loaded is None:
        return None
    return loaded[0].MolToSmiles(mol)


def element_counts(mol):
    """Count atoms of each element, including implicit hydrogens.

    Parameters
    ----------
    mol : Mol
        Molecule as produced by :func:`mol_from_smiles`.

    Returns
    -------
    counts : dict or None
        Mapping of element symbol to atom count (hydrogens included via
        implicit-H bookkeeping), or ``None`` if `mol` is ``None``.

    Examples
    --------
    >>> element_counts(mol_from_smiles("CCO"))
    {'C': 2, 'H': 6, 'O': 1}
    """
    if mol is None:
        return None
    counts = {}
    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        counts[symbol] = counts.get(symbol, 0) + 1
        num_h = atom.GetTotalNumHs()
        if num_h:
            counts["H"] = counts.get("H", 0) + num_h
    return counts


def molecular_weight(mol):
    """Return the average molecular weight of `mol` in kg/mol.

    Uses RDKit's average-atomic-mass molecular weight (the same convention
    as ``thermo``).

    Parameters
    ----------
    mol : Mol
        Molecule as produced by :func:`mol_from_smiles`.

    Returns
    -------
    mw : float or None
        Average molecular weight [kg/mol], or ``None`` if `mol` is
        ``None``.

    Examples
    --------
    >>> round(molecular_weight(mol_from_smiles("CCO")), 5)
    0.04607
    """
    loaded = load_rdkit()
    if mol is None or loaded is None:
        return None
    return loaded[1].MolWt(mol) / 1000.0
