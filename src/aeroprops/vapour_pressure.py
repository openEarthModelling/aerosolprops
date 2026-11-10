# SPDX-License-Identifier: MIT
"""Vapour-pressure estimation methods (object-oriented wrapper layer).

This module owns everything that needs RDKit: molecule construction,
SMARTS matching, and graph algorithms. The coefficient catalogues and
the pure estimators live in :mod:`aeroprops.methods.simpol1`; the
graph primitives live in :mod:`aeroprops.topology`.

Why not ``thermo.group_contribution``? thermo's fragmentation
machinery (``smarts_fragment_priority`` / ``BaseGroupContribution``)
partitions a molecule into *mutually exclusive* SMARTS groups: every
heavy atom must be matched by exactly one pattern and an atom-balance
check must close. SIMPOL.1 group counts are *overlapping and additive
by design* — the carbon count ``[#6]`` covers every carbon that other
groups also touch, the amide acid-side carbons (group 2) overlap the
carbon count (group 1), and the nitrophenol/nitroester groups (29/30)
*replace* aromatic-hydroxyl / ester counts (17/11) instead of
partitioning them. No priority-ordered exclusive partition can express
that, so aeroprops ships its own minimal fragmenter here, keeping the
ergonomics of the thermo wrappers (one-time fragmentation in the
constructor, ``counts``/``success``/``status`` attributes, estimators
returning ``None`` on failure).

Examples
--------
>>> sim = SIMPOL1(smi="C1=CCCCC1")  # cyclohexene
>>> sim.success, sim.status
(True, 'OK')
>>> sim.counts  # zeroth group + 6 carbons + 1 alicyclic ring + 1 C=C
{0: 1, 1: 6, 4: 1, 5: 1}
>>> round(sim.log10_vapour_pressure(293.15), 4)  # log10 Pa (paper: -0.94 atm)
4.0567
>>> round(sim.vapour_pressure(293.15), 1)  # Pa
11395.4
>>> round(sim.vapourisation_enthalpy(333.15) / 1000.0, 1)  # kJ/mol
33.4

Fragmentation happens once; molecules with elements the method does
not parameterise fail with ``success == False`` and all estimators
return ``None``:

>>> bad = SIMPOL1(smi="CCCl")  # SIMPOL.1 covers C, H, N, O only
>>> bad.success
False
>>> bad.counts, bad.vapour_pressure() is None
({}, True)
"""

from .methods.simpol1 import (
    SIMPOL1_GROUPS,
    SIMPOL1_SMARTS,
    simpol1_bk,
    simpol1_log10_vapour_pressure,
    simpol1_vapourisation_enthalpy,
)
from .methods.myrdal_yalkowsky import (
    myrdal_yalkowsky_hbn,
    myrdal_yalkowsky_tau,
)
from .structure import (
    load_rdkit,
    mol_from_smiles,
    molecular_weight,
    normalize_mol,
)
from .topology import fragments_excluding_bond

__all__ = [
    "SIMPOL1",
    "MyrdalYalkowsky",
    "myrdal_yalkowsky_descriptors",
    "simpol1_fragmentation",
]

#: Elements covered by the SIMPOL.1 group catalogue (paper Table 5:
#: hydrocarbons plus O/N functionality only).
_SIMPOL1_ELEMENTS = frozenset(("C", "H", "N", "O"))

#: Pattern-atom index that identifies a counted instance of each
#: SMARTS-countable group (the "anchor"). Counts are the number of
#: distinct anchors, which deduplicates symmetric/scaled matches the
#: way OpenBabel ``findall`` + set() did in the paper's reference
#: implementations. Groups 2/3/4 are algorithmic (no SMARTS); groups
#: 5/26/27 are bond-anchored (see _BOND_ANCHOR).
_SINGLE_ATOM_ANCHOR = {
    1: 0, 6: 1, 7: 0, 8: 0, 9: 1, 10: 0, 11: 0, 12: 1, 13: 0, 14: 1,
    15: 0, 16: 0, 17: 0, 18: 1, 19: 1, 20: 1, 21: 1, 22: 0, 23: 0,
    24: 0, 25: 0, 28: 0, 29: 0, 30: 0,
}

#: Bond-anchored groups: (begin, end) pattern-atom indices of the
#: counted bond (5: the C=C double bond; 26/27: the O-O bond).
_BOND_ANCHOR = {5: (0, 1), 26: (1, 2), 27: (1, 2)}

_pattern_cache = None


def _compiled_patterns():
    """Compile (once per process) the SIMPOL.1 SMARTS catalogue.

    Returns
    -------
    patterns : dict or None
        Group index -> compiled pattern, or ``None`` if RDKit is not
        importable. Algorithmic groups (value ``None`` in
        :data:`SIMPOL1_SMARTS`) are not included.
    """
    global _pattern_cache
    if _pattern_cache is None:
        loaded = load_rdkit()
        if loaded is None:
            return None
        cache = {}
        for k, smarts in SIMPOL1_SMARTS.items():
            if smarts is None:
                continue
            pattern = loaded[0].MolFromSmarts(smarts)
            if pattern is None:
                # A catalogue bug (should have been caught in review);
                # fail loudly rather than silently under-count.
                raise ValueError(
                    "SIMPOL.1 SMARTS for group %d does not parse: %r" % (k, smarts)
                )
            cache[k] = pattern
        _pattern_cache = cache
    return _pattern_cache


def simpol1_fragmentation(mol):
    """Count the 31 SIMPOL.1 groups of a normalised molecule.

    Implements the counting semantics of paper Table 5 with its
    footnotes:

    * every match set is deduplicated on its counting anchor (atom or
      bond), so each group instance counts exactly once;
    * groups overlap additively: the carbon count (group 1) covers
      *all* carbons including those in any functional group, and the
      amide acid-side carbons (group 2) count on top of the carbon
      count;
    * carbonyl-attached oxygens are "claimed" by the carbonyl class
      groups (10 carboxylic acid, 11 ester, 25 PAN, 28 peroxy acid) and
      are then excluded from the ether/hydroxyl groups (7, 12, 13, 14);
      the O-O bond of peroxy acids is likewise excluded from the
      peroxide/hydroperoxide groups (26, 27);
    * nitrophenol (29) replaces aromatic hydroxyl (17) and nitroester
      (30) replaces ester (11) -- replacement, not addition;
    * an aromatic amine anchor (21) that is itself aromatic (e.g.
      pyridine N), amidic, or O-bonded does not count.

    Groups 2 (amide acid-side carbons), 3 (aromatic rings) and 4
    (non-aromatic rings) are computed graphically: group 2 by cutting
    each amide C-N bond and counting carbons in the acyl component
    (:func:`aeroprops.topology.fragments_excluding_bond`); groups 3/4
    from SSSR rings classified by aromaticity.

    Parameters
    ----------
    mol : Mol
        Normalised RDKit molecule (as produced by
        :func:`aeroprops.structure.mol_from_smiles`; pass user-supplied
        ``Mol`` objects through
        :func:`aeroprops.structure.normalize_mol` first). Must contain
        only C, H, N and O.

    Returns
    -------
    counts : dict
        Non-zero group counts keyed by group index 0-30. The zeroth
        group is always present (``{0: 1}``) on success. Empty on
        failure.
    success : bool
    status : str
        ``"OK"`` on success, otherwise a reason.

    Examples
    --------
    >>> counts, success, status = simpol1_fragmentation(
    ...     mol_from_smiles("O=C(O)CCCCC(=O)O"))
    >>> success, counts  # adipic acid: 6 C + 2 COOH
    (True, {0: 1, 1: 6, 10: 2})
    >>> simpol1_fragmentation(mol_from_smiles("CCCl"))[1]
    False
    """
    counts = {0: 1}
    loaded = load_rdkit()
    if loaded is None:
        return {}, False, "RDKit unavailable"
    Chem = loaded[0]
    patterns = _compiled_patterns()
    if mol is None or patterns is None:
        return {}, False, "Failed to construct mol"
    unsupported = sorted({a.GetSymbol() for a in mol.GetAtoms()} - _SIMPOL1_ELEMENTS)
    if unsupported:
        return {}, False, "unsupported element(s): " + ", ".join(unsupported)

    # --- raw SMARTS matches and anchor sets -------------------------
    hits = {k: mol.GetSubstructMatches(patterns[k], uniquify=True) for k in patterns}
    anchors = {k: {mt[pos] for mt in hits[k]} for k, pos in _SINGLE_ATOM_ANCHOR.items()}
    for k, (i, j) in _BOND_ANCHOR.items():
        anchors[k] = {frozenset((mt[i], mt[j])) for mt in hits[k]}

    # --- rings (groups 3/4) ----------------------------------------
    for ring in Chem.GetSSSR(mol):
        if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
            counts[3] = counts.get(3, 0) + 1
        else:
            counts[4] = counts.get(4, 0) + 1

    # --- amide acid-side carbons (group 2, graph cut) ---------------
    amide_carbonyls = anchors[22] | anchors[23] | anchors[24]
    amide_nitrogens = set()
    if amide_carbonyls:
        n_acid_side = 0
        for c in amide_carbonyls:
            carbonyl_c = mol.GetAtomWithIdx(c)
            for neighbour in carbonyl_c.GetNeighbors():
                if neighbour.GetSymbol() == "N":
                    amide_nitrogens.add(neighbour.GetIdx())
                    fragments = fragments_excluding_bond(
                        mol, c, neighbour.GetIdx()
                    )
                    acyl = next(f for f in fragments if c in f)
                    n_acid_side += sum(
                        1 for i in acyl if mol.GetAtomWithIdx(i).GetSymbol() == "C"
                    )
        counts[2] = n_acid_side

    # --- oxygen claims of carbonyl-class groups ----------------------
    claimed_oxygen = set()
    for k in (10, 11, 25, 28):
        for c in anchors[k]:
            for neighbour in mol.GetAtomWithIdx(c).GetNeighbors():
                if neighbour.GetSymbol() == "O":
                    claimed_oxygen.add(neighbour.GetIdx())
    # O-O bonds inside peroxy acids (group 28) are claimed too, so the
    # peroxide (26) / hydroperoxide (27) groups do not double-count them.
    claimed_oo_bonds = set()
    for c in anchors[28]:
        carbonyl_c = mol.GetAtomWithIdx(c)
        peroxide_o = next(
            nb for nb in carbonyl_c.GetNeighbors()
            if nb.GetSymbol() == "O"
            and any(nb2.GetSymbol() == "O" for nb2 in nb.GetNeighbors())
        )
        for neighbour in peroxide_o.GetNeighbors():
            if neighbour.GetSymbol() == "O":
                claimed_oo_bonds.add(frozenset((peroxide_o.GetIdx(), neighbour.GetIdx())))

    # --- nitroester replacement of plain ester (30 replaces 11) -----
    # Walk N (group 30 anchor) -> alpha carbon -> carbonyl carbon; if
    # that carbonyl is a plain ester anchor it is replaced (PAN-type
    # nitrate on the alpha carbon re-classifies the ester).
    replaced_esters = set()
    for n in anchors[30]:
        nitro_n = mol.GetAtomWithIdx(n)
        for alpha_c in nitro_n.GetNeighbors():
            if alpha_c.GetSymbol() == "C":
                for carbonyl_c in alpha_c.GetNeighbors():
                    if (
                        carbonyl_c.GetIdx() != n
                        and carbonyl_c.GetIdx() in anchors[11]
                    ):
                        replaced_esters.add(carbonyl_c.GetIdx())

    # --- final counts ------------------------------------------------
    counts[7] = len(anchors[7] - claimed_oxygen)
    counts[11] = len(anchors[11] - replaced_esters)
    counts[17] = len(anchors[17] - anchors[29])  # nitrophenol replaces phenol
    for k in (12, 13, 14):
        counts[k] = len(anchors[k] - claimed_oxygen)
    counts[26] = len(anchors[26] - claimed_oo_bonds)
    counts[27] = len(anchors[27] - claimed_oo_bonds)
    aromatic_amines = set()
    for n in anchors[21]:
        nitrogen = mol.GetAtomWithIdx(n)
        if nitrogen.GetIsAromatic():
            continue  # pyridine-type N is not an "aromatic amine" group
        if n in amide_nitrogens:
            continue
        if any(nb.GetSymbol() == "O" for nb in nitrogen.GetNeighbors()):
            continue
        aromatic_amines.add(n)
    counts[21] = len(aromatic_amines)
    for k in (1, 5, 6, 8, 9, 10, 15, 16, 18, 19, 20, 22, 23, 24, 25, 28, 29, 30):
        if hits[k]:
            counts[k] = len(anchors[k])

    # Sorted by group index for reproducible display; dict equality is
    # key-based, so callers comparing to literal dicts are unaffected.
    return dict(sorted((k, v) for k, v in counts.items() if v)), True, "OK"


class SIMPOL1:
    """SIMPOL.1 sub-cooled liquid vapour pressure (Pankow & Asher, 2008).

    Wraps :func:`simpol1_fragmentation` and the pure estimators of
    :mod:`aeroprops.methods.simpol1` in the thermo ``Joback`` idiom:
    fragmentation runs once in the constructor, results are stored on
    the instance, and every estimation method returns ``None`` if
    fragmentation failed.

    Parameters
    ----------
    mol : Mol, optional
        Normalisation is applied internally; pass any parsed RDKit
        molecule.
    smi : str, optional
        SMILES string, used when `mol` is not given.

    Attributes
    ----------
    counts : dict
        Non-zero SIMPOL.1 group counts (empty on failure); group
        indices 0-30 follow paper Table 5.
    success : bool
    status : str
        ``"OK"`` on success, otherwise the failure reason.

    Examples
    --------
    >>> from aeroprops.methods.simpol1 import LOG10_ATM_TO_PA
    >>> sim = SIMPOL1(smi="C1=CCCCC1")
    >>> round(sim.log10_vapour_pressure(293.15) - LOG10_ATM_TO_PA, 2)
    -0.95
    """

    def __init__(self, mol=None, smi=None):
        if mol is None:
            mol = None if smi is None else mol_from_smiles(smi)
            if mol is None:
                self.counts, self.success, self.status = {}, False, "Failed to construct mol"
                return
        else:
            mol = normalize_mol(mol)
            if mol is None:
                self.counts, self.success, self.status = {}, False, "Failed to construct mol"
                return
        self.counts, self.success, self.status = simpol1_fragmentation(mol)

    def log10_vapour_pressure(self, T=298.15):
        """Base-10 logarithm of the vapour pressure [log10 Pa].

        >>> sim = SIMPOL1(smi="OC(=O)CCCCC(=O)O")
        >>> round(sim.log10_vapour_pressure(293.15), 4)  # adipic acid
        -2.9843
        """
        if not self.success:
            return None
        return simpol1_log10_vapour_pressure(self.counts, T)

    def vapour_pressure(self, T=298.15):
        """Vapour pressure [Pa].

        >>> sim = SIMPOL1(smi="OC(=O)CCCCC(=O)O")
        >>> round(sim.vapour_pressure(293.15), 3)  # adipic acid
        0.001
        """
        log10p = self.log10_vapour_pressure(T)
        return None if log10p is None else 10.0 ** log10p

    def vapourisation_enthalpy(self, T=298.15):
        """Enthalpy of vaporisation [J/mol] (paper Eq. 10-11).

        >>> sim = SIMPOL1(smi="C1=CCCCC1")
        >>> round(sim.vapourisation_enthalpy(300.0) / 1000.0, 1)
        38.5
        """
        if not self.success:
            return None
        return simpol1_vapourisation_enthalpy(self.counts, T)


# =========================================================================
# Myrdal-Yalkowsky (Myrdal & Yalkowsky, IECR 36, 2494, 1997)
# =========================================================================

#: Elements in the method's domain: the paper's training/test compounds
#: are CHNO organics plus halogenated and sulfur species (Table 1 has a
#: thiol and a chloroalkane); anything else fails cleanly.
MYRDAL_YALKOWSKY_ELEMENTS = frozenset(
    ("C", "H", "N", "O", "S", "F", "Cl", "Br", "I")
)

#: Hydrogen-bond group SMARTS, authored clean-room from the word
#: definitions of paper Eq. 9 ("the number of alcohols, carboxylic
#: acids, or primary amines") plus the published extension of
#: Compernolle et al. (ACP 10, 6271, 2010, section 2.6: hydroperoxides
#: and peracids join the alcohol/acid root). Key = descriptor key.
#: ``hydroperoxide`` excludes the peracid carbon so the two patterns
#: never double-count; ``carboxyl`` excludes the O-O-bonded acid
#: hydroxyl; ``primary_amine`` excludes amide nitrogens.
MYRDAL_YALKOWSKY_HBOND_SMARTS = {
    "hydroxyl": "[OX2H;$([OX2H]-[c]),$([OX2H]-[CX4])]",
    "carboxyl": "[CX3](=[OX1])-[OX2H;!$([OX2H]-[OX2])]",
    "hydroperoxide": "[#6;!$([CX3]=[OX1])]-[OX2]-[OX2H]",
    "peracid": "[CX3](=[OX1])-[OX2]-[OX2H]",
    "primary_amine": "[NX3;H2;!$([NX3]-[CX3]=[OX1])]",
}

#: Functional groups outside the published parameterisation that this
#: module refuses rather than silently mis-estimating (the adjudicated
#: scope ruling): UManSysProp extrapolates nitrates/PAN into the
#: flexibility and hydrogen-bond counts with no published basis, and
#: its enol handling likewise; aeroprops reports the molecule as
#: unsupported so the caller falls through to SIMPOL.1/EVAPORATION,
#: whose parameterisations cover these groups.
MYRDAL_YALKOWSKY_REFUSED_SMARTS = {
    "organic nitrate": "[#6]-[OX2]-[N+](=[OX1])[O-]",
    "peroxy acyl nitrate": "[CX3](=[OX1])-[OX2]-[OX2]-[N+](=[OX1])[O-]",
    "enol": "[OX2H]-[CX3]=[CX3]",
}

_my_pattern_cache = None


def _my_compiled_patterns():
    """Compile (once) the Myrdal-Yalkowsky SMARTS catalogues."""
    global _my_pattern_cache
    if _my_pattern_cache is None:
        loaded = load_rdkit()
        if loaded is None:
            return None
        compiled = {}
        for catalogue in (
            MYRDAL_YALKOWSKY_HBOND_SMARTS,
            MYRDAL_YALKOWSKY_REFUSED_SMARTS,
        ):
            for key, smarts in catalogue.items():
                pattern = loaded[0].MolFromSmarts(smarts)
                if pattern is None:
                    raise ValueError(
                        "Myrdal-Yalkowsky SMARTS %r does not parse: %r"
                        % (key, smarts)
                    )
                compiled[key] = pattern
        _my_pattern_cache = compiled
    return _my_pattern_cache


def myrdal_yalkowsky_descriptors(mol):
    """Compute the Myrdal-Yalkowsky descriptors of a normalised molecule.

    Implements the paper's word definitions for Eq. 8 and Eq. 9 as the
    conventions locked against paper Table 1 (all 19 rows reproduce at
    the printed two decimals):

    * **SP3** - nonring, nonterminal sp3 atoms among C/N/O/S (paper:
      "CH2 and CH and also includes NH, N, O, and S"), where
      *nonterminal* means heavy-atom degree >= 2 and *tert-butyl
      centres* (a carbon with three terminal-carbon neighbours) are
      excluded ("tert-butyl groups do not contribute to flexibility");
      atoms with a triple bond are skipped ("triple bonds ... are not
      counted");
    * **SP2** - nonring sp2 *carbons* with heavy-atom degree >= 2 (the
      "-C=" of the paper; valeric acid and methyl salicylate lock the
      carbonyl carbon at the half weight);
    * **RING** - the number of *independent ring systems*: connected
      components of the ring atoms joined by ring bonds (naphthalene
      and anthracene are one system each, biphenyl two - deliberately
      not the SSSR ring count, see the tests);
    * **tau** - :func:`aeroprops.methods.myrdal_yalkowsky.myrdal_yalkowsky_tau`;
    * **HBN groups** - :data:`MYRDAL_YALKOWSKY_HBOND_SMARTS` counts,
      combined by :func:`aeroprops.methods.myrdal_yalkowsky.myrdal_yalkowsky_hbn`
      with the all-element molecular weight.

    Parameters
    ----------
    mol : Mol
        Normalised RDKit molecule (see
        :func:`aeroprops.structure.mol_from_smiles`).

    Returns
    -------
    descriptors : dict
        Keys ``sp3``, ``sp2``, ``ring_systems``, ``tau``, the five
        hydrogen-bond group counts (``hydroxyl``, ``carboxyl``,
        ``hydroperoxide``, ``peracid``, ``primary_amine``),
        ``hbond_root`` (the first-root sum), ``hbn`` and ``mw`` [g/mol];
        empty on failure.
    success : bool
    status : str

    Examples
    --------
    Paper Table 1 anchors: 3,5-diethylphenol (Appendix A tau = 1.5),
    2,4,4-trimethylhexane (tert-butyl centre excluded, SP3 = 4,
    tau = 3) and 1-hexanethiol (terminal -SH not counted, tau = 4):

    >>> d, ok, status = myrdal_yalkowsky_descriptors(
    ...     mol_from_smiles("CCc1cc(O)cc(CC)c1"))
    >>> (d["sp3"], d["ring_systems"], d["tau"], d["hydroxyl"])
    (2, 1, 1.5, 1)
    >>> d, _, _ = myrdal_yalkowsky_descriptors(
    ...     mol_from_smiles("CC(C)CC(C)(C)CC"))
    >>> (d["sp3"], d["tau"])
    (4, 3.0)
    >>> d, _, _ = myrdal_yalkowsky_descriptors(mol_from_smiles("CCCCCCS"))
    >>> (d["sp3"], d["tau"])
    (5, 4.0)

    The Compernolle (2010) extension and the ring-system convention:

    >>> d, _, _ = myrdal_yalkowsky_descriptors(mol_from_smiles("CC(=O)OO"))
    >>> (d["peracid"], d["hydroperoxide"], d["carboxyl"])
    (1, 0, 0)
    >>> d, _, _ = myrdal_yalkowsky_descriptors(
    ...     mol_from_smiles("c1ccc2cc3ccccc3cc2c1"))
    >>> d["ring_systems"]  # anthracene: one system (SSSR would say 3)
    1

    Outside the published parameterisation the molecule is refused:

    >>> myrdal_yalkowsky_descriptors(mol_from_smiles("CCO[N+](=O)[O-]"))[1:]
    (False, 'unsupported functional group: organic nitrate')
    """
    loaded = load_rdkit()
    patterns = _my_compiled_patterns()
    if mol is None or patterns is None:
        return {}, False, "Failed to construct mol"
    Chem = loaded[0]
    unsupported = sorted(
        {a.GetSymbol() for a in mol.GetAtoms()} - MYRDAL_YALKOWSKY_ELEMENTS
    )
    if unsupported:
        return {}, False, "unsupported element(s): " + ", ".join(unsupported)
    for key in MYRDAL_YALKOWSKY_REFUSED_SMARTS:
        if mol.GetSubstructMatches(patterns[key], uniquify=True):
            return {}, False, "unsupported functional group: " + key

    # --- SP3 / SP2 walk (paper Eq. 8 word definitions) ----------------
    ring_info = mol.GetRingInfo()
    sp3 = sp2 = 0
    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in ("C", "N", "O", "S"):
            continue  # halogens are terminal in this method's domain
        if atom.GetIsAromatic() or ring_info.NumAtomRings(atom.GetIdx()):
            continue  # "nonring"
        heavy = [nb for nb in atom.GetNeighbors() if nb.GetSymbol() != "H"]
        types = {bond.GetBondTypeAsDouble() for bond in atom.GetBonds()}
        if 3.0 in types:
            continue  # "triple bonds ... are not counted"
        if 2.0 in types:
            if atom.GetSymbol() == "C" and len(heavy) >= 2:
                sp2 += 1
            continue
        if len(heavy) < 2:
            continue  # nonterminal
        if atom.GetSymbol() == "C":
            n_methyl = sum(
                1
                for nb in heavy
                if nb.GetSymbol() == "C"
                and len(
                    [x for x in nb.GetNeighbors() if x.GetSymbol() != "H"]
                ) == 1
            )
            if n_methyl >= 3:
                continue  # tert-butyl centre
        sp3 += 1

    # --- ring systems: ring atoms joined by ring bonds ----------------
    rw = Chem.RWMol(mol)
    for bond in mol.GetBonds():
        if not bond.IsInRing():
            rw.RemoveBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
    ring_systems = sum(
        1 for frag in Chem.GetMolFrags(rw) if len(frag) >= 3
    )

    # --- hydrogen-bond groups (paper Eq. 9 + Compernolle extension) ---
    counts = {}
    for key in MYRDAL_YALKOWSKY_HBOND_SMARTS:
        counts[key] = len(mol.GetSubstructMatches(patterns[key], uniquify=True))
    hbond_root = (
        counts["hydroxyl"] + counts["carboxyl"]
        + counts["hydroperoxide"] + counts["peracid"]
    )
    mw = None if mol is None else molecular_weight(mol) * 1000.0  # g/mol
    tau = myrdal_yalkowsky_tau(sp3, sp2, ring_systems)
    hbn = myrdal_yalkowsky_hbn(hbond_root, counts["primary_amine"], mw)
    descriptors = dict(
        sp3=sp3, sp2=sp2, ring_systems=ring_systems, tau=tau,
        hbond_root=hbond_root, hbn=hbn, mw=mw, **counts
    )
    return descriptors, True, "OK"
