# SPDX-License-Identifier: MIT
"""Molecular-graph topology primitives.

Thin wrappers over RDKit built-ins (connected components via
``Chem.GetMolFrags``, ring data via ``mol.GetRingInfo()``, neighbour
queries via ``atom.GetNeighbors()``). aeroprops deliberately contains no
home-grown graph algorithms: what is method-specific is the translation of
paper semantics (e.g. the SIMPOL.1 "acid-side carbons" of an amide, the
EVAPORATION "in-chain oxygen") onto the primitives defined here, and those
translations live in the estimator modules.

All functions take molecules as produced by
:func:`aeroprops.structure.mol_from_smiles` and return ``None`` when RDKit
is unavailable or their preconditions are not met.
"""

from aeroprops.structure import load_rdkit


def fragments_excluding_bond(mol, begin, end):
    """Connected components of the molecular graph with one bond deleted.

    Atom indices are unchanged with respect to `mol`. Used e.g. for the
    SIMPOL.1 group-2 "acid-side carbon count": the carbons on the acyl side
    of an amide are exactly the component containing the carbonyl carbon
    after deleting the C-N bond.

    Parameters
    ----------
    mol : Mol
        Molecule.
    begin, end : int
        Atom indices of the bond to delete.

    Returns
    -------
    fragments : tuple of tuple of int or None
        One tuple of atom indices per connected component.

    Examples
    --------
    >>> from aeroprops.structure import mol_from_smiles
    >>> propane = mol_from_smiles("CCC")
    >>> fragments_excluding_bond(propane, 0, 1)
    ((0,), (1, 2))
    """
    loaded = load_rdkit()
    if mol is None or loaded is None:
        return None
    Chem = loaded[0]
    rw = Chem.RWMol(mol)
    rw.RemoveBond(begin, end)
    return tuple(tuple(frag) for frag in Chem.GetMolFrags(rw))


def separate_without_atom(mol, idx):
    """Connected components of the molecular graph with one atom deleted.

    Fragment members are reported with their **original** atom indices
    (the deleted index does not appear). Used e.g. for the EVAPORATION
    "in-chain oxygen" test: an oxygen whose removal splits the molecule
    into fragments that both contain carbon.

    Parameters
    ----------
    mol : Mol
        Molecule.
    idx : int
        Atom index to delete.

    Returns
    -------
    fragments : tuple of tuple of int or None
        One tuple of (original) atom indices per connected component.

    Examples
    --------
    Removing the bridging oxygen of diethyl ether splits the carbon
    skeleton in two; removing ethanol's terminal oxygen does not:

    >>> from aeroprops.structure import mol_from_smiles
    >>> ether = mol_from_smiles("CCOCC")
    >>> separate_without_atom(ether, 2)
    ((0, 1), (3, 4))
    >>> separate_without_atom(mol_from_smiles("CCO"), 2)
    ((0, 1),)
    """
    loaded = load_rdkit()
    if mol is None or loaded is None:
        return None
    Chem = loaded[0]
    rw = Chem.RWMol(mol)
    rw.RemoveAtom(idx)
    # Atom indices above `idx` shift down by one after deletion; map them
    # back to the original numbering.
    fragments = tuple(
        tuple(i if i < idx else i + 1 for i in frag)
        for frag in Chem.GetMolFrags(rw)
    )
    return fragments


def ring_arc_distance(mol, idx1, idx2):
    """Shorter topological distance between two ring atoms along a shared
    ring.

    The distance is the number of ring bonds on the shorter arc of the
    smallest shared ring (1 = ortho, 2 = meta, 3 = para on a six-ring).
    If the atoms share several rings, the minimum over those rings is
    returned.

    Parameters
    ----------
    mol : Mol
        Molecule.
    idx1, idx2 : int
        Atom indices.

    Returns
    -------
    distance : int or None
        Number of ring bonds on the shorter arc, or ``None`` if the two
        atoms share no ring.

    Examples
    --------
    >>> from aeroprops.structure import mol_from_smiles
    >>> benzene = mol_from_smiles("c1ccccc1")
    >>> [ring_arc_distance(benzene, 0, j) for j in (1, 2, 3, 4, 5)]
    [1, 2, 3, 2, 1]
    >>> ring_arc_distance(mol_from_smiles("CCc1ccccc1"), 0, 4) is None
    True
    """
    loaded = load_rdkit()
    if mol is None or loaded is None:
        return None
    rings = mol.GetRingInfo().AtomRings()
    best = None
    for ring in rings:
        if idx1 in ring and idx2 in ring:
            n = len(ring)
            p1, p2 = ring.index(idx1), ring.index(idx2)
            d = min(abs(p1 - p2), n - abs(p1 - p2))
            if best is None or d < best:
                best = d
    return best
