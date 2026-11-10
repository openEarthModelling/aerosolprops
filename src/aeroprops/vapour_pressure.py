# SPDX-License-Identifier: MIT
"""Vapour-pressure estimation methods (object-oriented wrapper layer).

This module owns everything that needs RDKit: molecule construction,
SMARTS matching, and graph algorithms. The coefficient catalogues and
the pure estimators live in :mod:`aeroprops.methods.simpol1` and
:mod:`aeroprops.methods.myrdal_yalkowsky`; the graph primitives live
in :mod:`aeroprops.topology`.

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

import warnings

from .methods.simpol1 import (
    SIMPOL1_GROUPS,
    SIMPOL1_SMARTS,
    simpol1_bk,
    simpol1_log10_vapour_pressure,
    simpol1_vapourisation_enthalpy,
)
from .methods.myrdal_yalkowsky import (
    myrdal_yalkowsky_hbn,
    myrdal_yalkowsky_log10_vapour_pressure,
    myrdal_yalkowsky_tau,
)
from .structure import (
    load_rdkit,
    mol_from_smiles,
    molecular_weight,
    normalize_mol,
)
from .topology import (
    fragments_excluding_bond,
    ring_arc_distance,
    separate_without_atom,
)

__all__ = [
    "SIMPOL1",
    "MyrdalYalkowsky",
    "myrdal_yalkowsky_descriptors",
    "simpol1_fragmentation",
    "EVAPORATION",
    "evaporation_fragmentation",
    "NANNOOLAL_VP_SMARTS",
    "nannoolal_vp_fragmentation",
    "NannoolalVP",
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


def _my_joback_boiling_point(smi):
    """Normal boiling point [K] from ``thermo``'s Joback method, or None.

    thermo is imported lazily (the Joback route is its only consumer in
    aeroprops; the class is re-exported at the package top level in
    thermo 0.6.1). ``Joback.Tb_estimated`` holds a number even when the
    SMARTS fragmentation failed, so the route gates on ``success``.
    """
    try:
        from thermo import Joback
    except ImportError:
        return None
    try:
        joback = Joback(smi)
    except Exception:
        return None
    if not getattr(joback, "success", False):
        return None
    tb = getattr(joback, "Tb_estimated", None)
    try:
        return tb if tb > 0.0 else None
    except TypeError:
        return None


def _resolve_normal_boiling_point(mol, smi, tb_source):
    """Resolve a normal boiling point [K] from the selected route.

    Shared Tb routing of the *input* vapour-pressure methods (Myrdal-
    Yalkowsky Eq. 12, Nannoolal Part 3 Eq. 6): an experimental value
    must be passed in, or the estimate comes from Nannoolal Part 1
    (:class:`aeroprops.boiling_point.NannoolalBP`, the default route)
    or from ``thermo``'s Joback method (aeroprops does not reimplement
    Joback).

    Returns ``(tb, error)``: ``tb`` is ``None`` exactly when ``error``
    names why the route refused (the string the wrapper layer stores in
    its ``status``); on success the error is ``None``.
    """
    if tb_source == "experimental":
        return None, "experimental Tb required: pass tb=..."
    if tb_source == "joback":
        if smi is None:
            loaded = load_rdkit()
            smi = loaded[0].MolToSmiles(mol)
        tb = _my_joback_boiling_point(smi)
        if tb is None:
            return None, "Joback Tb unavailable (thermo refused)"
        return tb, None
    if tb_source == "nannoolal":
        from .boiling_point import NannoolalBP

        nnp = NannoolalBP(mol=mol)
        tb = nnp.boiling_point()
        if tb is None:
            return None, "Nannoolal Tb refused: " + nnp.status
        return tb, None
    return None, "unknown tb_source: %r" % (tb_source,)


class MyrdalYalkowsky:
    """Myrdal-Yalkowsky vapour pressure (Myrdal & Yalkowsky, IECR 36,
    2494, 1997), with the hydrogen-bond extension of Compernolle et al.
    (ACP 10, 6271, 2010).

    An *input* method: paper Eq. 12 needs the normal boiling point
    ``Tb`` and, for the solid branch, the melting point ``Tm``.
    ``Tb`` is routed from Nannoolal Part 1
    (:class:`aeroprops.boiling_point.NannoolalBP`, the default),
    ``thermo``'s built-in Joback method, or experimental data —
    aeroprops does not reimplement Joback. With ``Tm`` unset the
    estimate is the sub-cooled liquid vapour pressure (the paper's own
    convention for its Table 1 liquid rows), annotated in ``status``.

    Parameters
    ----------
    mol : Mol, optional
        Normalisation is applied internally; pass any parsed RDKit
        molecule.
    smi : str, optional
        SMILES string, used when `mol` is not given.
    tb : float, optional
        Experimental normal boiling point [K]; overrides `tb_source`
        (recorded as ``"experimental"``).
    tb_source : str, optional
        Where an absent ``tb`` comes from: ``"nannoolal"`` (default),
        ``"joback"`` or ``"experimental"``.
    tm : float, optional
        Melting point [K]; ``None`` (default) selects the sub-cooled
        liquid branch.
    sigma : float, optional
        External rotational symmetry number (>= 1, default 1); enters
        only the solid branch.

    Attributes
    ----------
    descriptors : dict
        The output of :func:`myrdal_yalkowsky_descriptors` (empty on
        failure).
    tb, tb_source, tm, sigma
        The resolved inputs (``tb`` ``None`` if its route failed).
    success : bool
    status : str
        ``"OK"`` on success, or ``"OK (sub-cooled liquid: Tm not
        supplied)"``; otherwise the failure reason.

    Examples
    --------
    Paper Table 1 anchor - resorcinol, Tb 549 K, T 425 K, predicted
    -1.76 log10 atm:

    >>> my = MyrdalYalkowsky(smi="Oc1cccc(O)c1", tb=549.0)
    >>> my.success, my.tb_source
    (True, 'experimental')
    >>> round(my.log10_vapour_pressure(425.0), 2)
    3.24
    >>> my.status
    'OK (sub-cooled liquid: Tm not supplied)'

    The solid branch (biphenyl, Tm 342.4 K, sigma 2, Table 1 -2.97):

    >>> my = MyrdalYalkowsky(smi="c1ccccc1c1ccccc1", tb=528.0, tm=342.4,
    ...                      sigma=2)
    >>> round(my.log10_vapour_pressure(342.0), 2)
    2.04
    """

    def __init__(self, mol=None, smi=None, tb=None, tb_source="nannoolal",
                 tm=None, sigma=1.0):
        # set on every path (including failures): the attributes of a
        # thermo-style wrapper always exist, never AttributeError
        self.tb_source = "experimental" if tb is not None else tb_source
        self.tm = tm
        self.sigma = sigma
        if mol is None:
            mol = None if smi is None else mol_from_smiles(smi)
        else:
            mol = normalize_mol(mol)
        if mol is None:
            self.mol = None
            self.descriptors, self.success, self.status = (
                {}, False, "Failed to construct mol"
            )
            self.tb = None
            return
        self.mol = mol
        self.descriptors, self.success, self.status = (
            myrdal_yalkowsky_descriptors(mol)
        )
        self.tb = float(tb) if tb is not None and tb > 0.0 else None
        if tb is not None and self.tb is None:
            self.success, self.status = False, "invalid experimental Tb"
            return
        if not self.success:
            return
        if tm is not None and not tm > 0.0:
            self.success, self.status = False, "invalid Tm"
            return
        if not sigma >= 1.0:
            self.success, self.status = False, "invalid sigma (< 1)"
            return
        if self.tb is None:
            self.tb = self._resolve_tb(smi)
            if self.tb is None:
                self.success = False
                return
        if self.tm is None:
            self.status = "OK (sub-cooled liquid: Tm not supplied)"

    def _resolve_tb(self, smi):
        """Fill ``self.tb`` from the selected route; sets the failure
        status and returns ``None`` when the route refuses."""
        tb, error = _resolve_normal_boiling_point(
            self.mol, smi, self.tb_source
        )
        if error is not None:
            self.status = error
        return tb

    def log10_vapour_pressure(self, T=298.15):
        """Base-10 logarithm of the vapour pressure [log10 Pa].

        Emits a warning above the normal boiling point, where Eq. 12
        extrapolates outside its fitted range.

        >>> my = MyrdalYalkowsky(smi="CC(C)CC(C)(C)CC", tb=404.0)
        >>> round(my.log10_vapour_pressure(323.0), 2)  # Table 1: -1.28 atm
        3.72
        """
        if not self.success:
            return None
        if not T > 0.0:
            return None
        if T > self.tb:
            warnings.warn(
                "T = %.1f K exceeds Tb = %.1f K: Eq. 12 extrapolates "
                "beyond the normal boiling point" % (T, self.tb),
                stacklevel=2,
            )
        d = self.descriptors
        return myrdal_yalkowsky_log10_vapour_pressure(
            T, self.tb, d["tau"], d["hbn"], Tm=self.tm, sigma=self.sigma
        )

    def vapour_pressure(self, T=298.15):
        """Vapour pressure [Pa].

        >>> my = MyrdalYalkowsky(smi="Oc1cccc(O)c1", tb=549.0)
        >>> round(my.vapour_pressure(425.0), 1)
        1756.9
        """
        log10p = self.log10_vapour_pressure(T)
        return None if log10p is None else 10.0 ** log10p


# ======================================================================
# EVAPORATION (Compernolle, Ceulemans & Mueller, ACP 11, 9431-9450,
# 2011, doi:10.5194/acp-11-9431-2011)
# ======================================================================
#
# Why not reuse thermo's fragmenter here either: EVAPORATION counts are
# even more aggressively overlapping than SIMPOL.1 -- the carbon count
# (k2) overlaps every functional group, the "12" ring group fires on
# ring oxygens *and* on ring carbons bearing OH/OOH, and k19 counts
# per-hydroxyl neighbourhoods that share atoms. aeroprops therefore
# counts EVAPORATION descriptors with explicit graph traversal (one
# shared atom classification plus two counting rule sets), which also
# lets the two c3 branching conventions live side by side.

from math import comb as _binomial

from .methods.evaporation import (
    evaporation_boiling_temperature,
    evaporation_log10_vapour_pressure,
)

#: Elements the EVAPORATION parameter set covers (paper data set:
#: C/H/O compounds plus alkyl nitrates and PANs; N is only accepted in
#: those two environments in paper mode).
_EVAPORATION_ELEMENTS = frozenset(("C", "H", "O", "N"))


def _evaporation_atom_maps(mol, Chem):
    """Per-atom topology maps shared by both counting modes.

    Returns a dict with, keyed by atom index:

    ``symbol``, ``nH`` (implicit + explicit total), ``sp3`` (RDKit
    SP3 hybridisation), ``in_ring``, ``minring`` (smallest ring size
    3-8 the atom belongs to, 0 if none -- rings of 9+ atoms get no
    ring correction in either mode), ``ncc`` (number of single-bonded
    carbon neighbours; aromatic C-C bonds are not single bonds),
    ``no0`` (number of single-bonded O neighbours carrying no H:
    ether/peroxide/ester/nitrate oxygens, but not OH), ``alkene``
    (carbon taking part in a non-carbonyl C=C double bond).
    """
    maps = {
        "symbol": {}, "nH": {}, "sp3": {}, "in_ring": {}, "minring": {},
        "ncc": {}, "no0": {}, "alkene": {},
    }
    for atom in mol.GetAtoms():
        i = atom.GetIdx()
        maps["symbol"][i] = atom.GetSymbol()
        maps["nH"][i] = atom.GetTotalNumHs()
        maps["sp3"][i] = atom.GetHybridization() == Chem.HybridizationType.SP3
        maps["in_ring"][i] = atom.IsInRing()
        maps["minring"][i] = 0
        for n in range(3, 9):
            if atom.IsInRingSize(n):
                maps["minring"][i] = n
                break
        maps["ncc"][i] = 0
        maps["no0"][i] = 0
        maps["alkene"][i] = False
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtom(), bond.GetEndAtom()
        ia, ib = a.GetIdx(), b.GetIdx()
        btype = bond.GetBondType()
        single = btype == Chem.BondType.SINGLE
        if a.GetSymbol() == "C" and b.GetSymbol() == "C":
            if single:
                maps["ncc"][ia] += 1
                maps["ncc"][ib] += 1
            elif btype == Chem.BondType.DOUBLE:
                maps["alkene"][ia] = True
                maps["alkene"][ib] = True
        if single:
            for x, y in ((a, b), (b, a)):
                if x.GetSymbol() == "C" and y.GetSymbol() == "O" \
                        and y.GetTotalNumHs() == 0:
                    maps["no0"][x.GetIdx()] += 1
    return maps


def _evaporation_groups(mol, Chem, maps):
    """Functional-group classification shared by both counting modes.

    Returns ``(groups, nbrs, btype)`` -- the group dict of index sets (see the source for the exact
    semantics; every rule was verified against the live UManSysProp
    oracle for the upstream mode and against the paper's worked
    examples for the paper mode):

    ``ketald``   carbonyl C without a single-bond O neighbour
                 (ketones, aldehydes, formaldehyde; Table 4 footnote e)
    ``aldehyde`` / ``ketone`` / ``formaldehyde``  ketald split by H count
    ``acyl_type``  carbonyl C -> 'acid' | 'ester' | 'peracid' |
                 'pan' | 'other' (carbonyl C carrying a single-bond O)
    ``ester_bridge``  ester acyl C -> its bridge O
    ``ester_o_h0``  ester bridge O whose acyl C carries no H (upstream
                 counts only these; formate bridge O is invisible
                 upstream but a full ester in the paper mode)
    ``ether_o``  O without H, two carbon neighbours, none of them a
                 carbonyl C (chain and ring ethers, acetals)
    ``oh``       hydroxyl O -> the carbon it sits on (OH on sp2 C
                 included; the paper mode filters enols out)
    ``ooh_c_side``  carbon-side O of hydroperoxides (C-O-O-H where the
                 carbon is not an acyl carbon; the acyl case is the
                 peracid)
    ``nitrate_n`` alkyl-nitrate nitrogens (three O neighbours, charge
                 form or N=O bond, and an O continuing into a non-acyl
                 carbon); ``pan_n`` the peroxynitrate nitrogens of PANs
                 (same N environment, but the O side continues into the
                 acyl group -- these fire k7, not k4)
    ``ring_oo``  O-O single bonds inside a ring with a carbon flank on
                 both sides (one frozenset per pair)
    ``peroxide_bridges``  number of acyclic C-O-O-C bridges
    ``alkene_bonds``  frozenset(a, b) per non-carbonyl C=C double bond
    """
    single, double = Chem.BondType.SINGLE, Chem.BondType.DOUBLE
    sym, nh = maps["symbol"], maps["nH"]
    nbrs = {a.GetIdx(): [n.GetIdx() for n in a.GetNeighbors()]
            for a in mol.GetAtoms()}
    btype = {}
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        btype[(i, j)] = btype[(j, i)] = bond.GetBondType()

    g = {
        "ketald": set(), "aldehyde": set(), "ketone": set(),
        "formaldehyde": set(), "acyl_type": {}, "ester_bridge": {},
        "ether_o": set(), "oh": {}, "ooh_c_side": set(),
        "nitrate_n": set(), "pan_n": set(), "ring_oo": [],
        "peroxide_bridges": 0, "alkene_bonds": [],
    }
    # Pass 1: candidate nitrate nitrogens and H-carrying / ether
    # oxygens (the carbonyl typing below needs the candidates for PAN
    # detection; alkyl nitrate vs PAN is decided after the carbonyl
    # set is known).
    for i, sy in sym.items():
        if sy == "N":
            o_nbrs = [j for j in nbrs[i] if sym[j] == "O"]
            charged = mol.GetAtomWithIdx(i).GetFormalCharge() >= 1
            if len(o_nbrs) == 3 and (charged or any(
                    btype[(i, j)] == double for j in o_nbrs)):
                g["nitrate_n"].add(i)  # candidates until the post-pass
        elif sy == "O" and nh[i] >= 1:
            if len(nbrs[i]) == 1 and sym[nbrs[i][0]] == "C":
                g["oh"][i] = nbrs[i][0]
            elif len(nbrs[i]) == 1 and sym[nbrs[i][0]] == "O":
                o2 = nbrs[i][0]
                if any(sym[j] == "C" for j in nbrs[o2]):
                    # C-side oxygen of C-O-O-H; the acyl case (peracid)
                    # is excluded in pass 2 via acyl_carbons below.
                    g["ooh_c_side"].add(o2)
        elif sy == "O" and nh[i] == 0:
            if len(nbrs[i]) == 2 and all(sym[j] == "C" for j in nbrs[i]):
                g["ether_o"].add(i)  # acyl check applied after pass 2

    # Pass 2: carbonyl carbons, typed by their single-bond oxygen.
    acyl_carbons = set()
    for i, sy in sym.items():
        if sy != "C" or not any(
            sym[j] == "O" and btype[(i, j)] == double for j in nbrs[i]
        ):
            continue
        acyl_carbons.add(i)
        o_single = [j for j in nbrs[i]
                    if sym[j] == "O" and btype[(i, j)] == single]
        if not o_single:
            g["ketald"].add(i)
            if nh[i] == 0:
                g["ketone"].add(i)
            elif nh[i] == 1:
                g["aldehyde"].add(i)
            else:
                g["formaldehyde"].add(i)
            continue
        bridge = o_single[0]
        others = [j for j in nbrs[bridge] if j != i]
        if nh[bridge] >= 1:
            g["acyl_type"][i] = "acid"
        elif any(sym[j] == "C" for j in others):
            g["acyl_type"][i] = "ester"
            g["ester_bridge"][i] = bridge
        elif any(sym[j] == "O" for j in others):
            o2 = next(j for j in others if sym[j] == "O")
            if nh[o2] >= 1:
                g["acyl_type"][i] = "peracid"
            elif any(j in g["nitrate_n"] for j in nbrs[o2]):
                g["acyl_type"][i] = "pan"
            else:
                g["acyl_type"][i] = "other"  # peroxy ester: uncounted
        else:
            g["acyl_type"][i] = "other"

    # Post-pass fixes that need the completed carbonyl set. An O-H on
    # an acyl carbon is the acid's own hydroxyl: claimed by k9 (and by
    # k8 neither upstream -- 801-804 never match it, adipic raw has no
    # 80x -- nor in the paper mode).
    g["oh"] = {o: c for o, c in g["oh"].items() if c not in acyl_carbons}
    g["ether_o"] = {i for i in g["ether_o"] if i not in acyl_carbons
                    and not any(j in acyl_carbons for j in nbrs[i])}
    g["ooh_c_side"] = {i for i in g["ooh_c_side"]
                       if i not in acyl_carbons
                       and not any(j in acyl_carbons for j in nbrs[i])}
    for (i, j), bt in btype.items():
        if i < j and sym[i] == "C" and sym[j] == "C" and bt == double:
            g["alkene_bonds"].append(frozenset((i, j)))
        if i < j and sym[i] == "O" and sym[j] == "O" and bt == single:
            if not any(sym[k] == "C" for k in nbrs[i]):
                continue
            if not any(sym[k] == "C" for k in nbrs[j]):
                continue
            # The upstream peroxide-bridge count includes ring O-O
            # bonds (1,2-dioxane fires 203 *and* 1202).
            g["peroxide_bridges"] += 1
            if mol.GetBondBetweenAtoms(i, j).IsInRing():
                g["ring_oo"].append(frozenset((i, j)))
    # Alkyl nitrate vs PAN: the N counts as a nitrate nitrogen only
    # when one of its oxygens continues into a *non-acyl* carbon
    # (R-O-NO2). A peroxynitrate N (RCO-O-O-NO2) belongs to the PAN
    # group instead and must not also fire k4.
    pan_n = set()
    for c, t in g["acyl_type"].items():
        if t == "pan":
            bridge = next(o for o in nbrs[c] if sym[o] == "O"
                          and btype[(c, o)] == single)
            o2 = next(j for j in nbrs[bridge]
                      if j != c and sym[j] == "O")
            for n in nbrs[o2]:
                if sym[n] == "N":
                    pan_n.add(n)
    g["nitrate_n"] = {
        n for n in g["nitrate_n"]
        if any(sym[j] == "C" and j not in g["acyl_type"]
               and j not in g["ketald"]
               for o in nbrs[n] if sym[o] == "O"
               for j in nbrs[o] if j != n)}
    g["pan_n"] = pan_n
    g["ester_o_h0"] = {o for c, o in g["ester_bridge"].items() if nh[c] == 0}
    return g, nbrs, btype


def _evaporation_counts_upstream(mol, Chem, maps, g, nbrs, btype):
    """UManSysProp counting semantics: submode counts `m` (verified
    against the live oracle, 26 Example_output goldens + 49 controls).

    The upstream nitrate/PAN SMARTS never match under RDKit charge
    normalisation, so c4/c7 and the PAN partner buckets stay zero by
    construction here (the paper mode counts both).
    """
    sym, nh, sp3, minring = (maps["symbol"], maps["nH"],
                             maps["sp3"], maps["minring"])
    single = Chem.BondType.SINGLE
    ester_acyl = {c for c, t in g["acyl_type"].items() if t == "ester"}
    m = {}

    def inc(k, v=1):
        m[k] = m.get(k, 0) + v

    # c2 ingredients
    m["2"] = sum(1 for s in sym.values() if s == "C")
    m["201"] = len(g["ether_o"])
    m["202"] = len(g["ester_o_h0"])
    m["203"] = g["peroxide_bridges"]
    # branching (subset-match model) and ring correction
    for i, s in sym.items():
        if s == "C" and sp3[i]:
            ncc, o0 = maps["ncc"][i], maps["no0"][i]
            if o0 == 0:
                inc("301", _binomial(ncc, 4))
                inc("302", _binomial(ncc, 3))
            elif ncc == 3:
                inc("303")
            elif ncc == 2 and o0 == 1:
                inc("304")
            elif ncc >= 1 and o0 >= 2:
                inc("306")
        if 3 <= minring[i] <= 8:
            inc(str(304 + minring[i]))
    # carbonyl groups
    for i in g["ketone"]:
        inc("502" if maps["in_ring"][i] else "501")
    for i in g["aldehyde"]:
        inc("503")
    m["7001"] = len(g["formaldehyde"])
    for o in g["ester_o_h0"]:
        inc("602" if maps["in_ring"][o] else "601")
    m["9"] = sum(1 for t in g["acyl_type"].values() if t == "acid")
    m["10"] = len(g["ooh_c_side"])
    m["11"] = sum(1 for t in g["acyl_type"].values() if t == "peracid")
    # hydroxyls: primary/secondary/tertiary by heavy-neighbour degree
    # of the carbinol carbon (2-methoxy-2-propanol fires 803 with only
    # two direct C neighbours -- the methoxy oxygen counts), enol OH
    # on a C=C carbon only (phenol OH matches nothing)
    for o, c in g["oh"].items():
        if sp3[c]:
            key = {2: "801", 3: "802", 4: "803"}.get(len(nbrs[c]))
            if key:
                inc(key)
        elif maps["alkene"][c]:
            inc("804")
    # ring group 12
    m["1201"] = sum(1 for i in g["ether_o"] if maps["in_ring"][i])
    m["1202"] = len(g["ring_oo"])
    m["1203"] = sum(1 for c in g["oh"].values()
                    if maps["in_ring"][c] and sp3[c])
    m["1205"] = 0
    for o in g["ooh_c_side"]:
        for j in nbrs[o]:
            if sym[j] == "C" and maps["in_ring"][j]:
                inc("1205")
    # k13: C=C conjugated to a ketone/aldehyde carbon (one per C=C)
    for pair in g["alkene_bonds"]:
        for i in pair:
            if any(sym[j] == "C" and btype[(i, j)] == single
                   and j in g["ketald"] for j in nbrs[i]):
                inc("13")
                break
    # k15: allylic hydroxyls (capped at 1 in the assembly)
    for o, c in g["oh"].items():
        if sp3[c] and any(sym[j] == "C" and maps["alkene"][j]
                          for j in nbrs[c] if btype[(c, j)] == single):
            inc("15")
    # k16 alpha interactions between carbonyls
    ald_ald_pairs = set()
    ketald_ester_pairs = set()
    for (i, j), bt in btype.items():
        if i < j and bt == single and sym[i] == "C" and sym[j] == "C":
            pair = (i, j)
            both_ketald = i in g["ketald"] and j in g["ketald"]
            if both_ketald and (i in g["aldehyde"] or j in g["aldehyde"]):
                ald_ald_pairs.add(pair)  # 1601, deduped per bond
            if (i in g["ketald"] and j in ester_acyl) or \
                    (j in g["ketald"] and i in ester_acyl):
                ketald_ester_pairs.add(pair)  # 1603
    m["1601"] = len(ald_ald_pairs)
    for k in g["ketone"]:
        m["1602"] = m.get("1602", 0) + sum(
            1 for j in nbrs[k]
            if sym[j] == "C" and btype[(k, j)] == single
            and j in g["ketald"])
    m["1603"] = len(ketald_ester_pairs)
    # k17 beta interactions (strictly 1,3-dicarbonyl: 1,4 pairs fire
    # nothing, verified on succinaldehyde and 4-oxopentanal)
    ald_beta_pairs, ketald_ester_beta_pairs = set(), set()
    for x in g["ketald"]:
        for j in nbrs[x]:
            if sym[j] != "C" or btype[(x, j)] != single:
                continue
            for beta in nbrs[j]:
                if beta == x or sym[beta] != "C" \
                        or btype[(j, beta)] != single:
                    continue
                if beta in g["ketald"]:
                    if x in g["aldehyde"]:
                        ald_beta_pairs.add(frozenset((x, beta)))  # 1701
                    if x in g["ketone"]:
                        inc("1702")  # per-ketone, not deduped
                if beta in ester_acyl and x in g["ketald"]:
                    ketald_ester_beta_pairs.add(frozenset((x, beta)))
    m["1701"] = len(ald_beta_pairs)
    m["1703"] = len(ketald_ester_beta_pairs)
    # k18: functional oxygens on the alpha carbon of ketones/aldehydes;
    # the alpha carbon itself must not be a carbonyl carbon (methyl
    # pyruvate: the bridge O on the ester carbonyl fires nothing)
    for x in g["ketald"]:
        for j in nbrs[x]:
            if sym[j] != "C" or btype[(x, j)] != single \
                    or j in g["ketald"] or j in g["acyl_type"]:
                continue
            for o in nbrs[j]:
                if sym[o] != "O" or btype[(j, o)] != single:
                    continue
                if nh[o] == 0:
                    inc("1801")  # any H-free O, whatever it bonds to
                elif nh[o] >= 1 and g["oh"].get(o) == j:
                    inc("1802")
    # k19: per-hydroxyl neighbourhoods. Vicinal OH pairs are recorded
    # once per pair (the assembly doubles them, reproducing the
    # per-OH two-sided count: ethylene glycol raw 1902 = 1).
    vicinal_pairs = set()
    for o, c in g["oh"].items():
        for a in nbrs[c]:
            if sym[a] != "C" or btype[(c, a)] != single:
                continue
            if sp3[a]:
                for p in nbrs[a]:
                    if sym[p] != "O" or btype[(a, p)] != single:
                        continue
                    if nh[p] >= 1 and g["oh"].get(p) == a:
                        vicinal_pairs.add(frozenset((o, p)))
                    elif nh[p] == 0:
                        # any H-free O on the adjacent carbon counts,
                        # even O-N(O)O drawn neutral (goldens 13/21)
                        inc("1901")
            t = g["acyl_type"].get(a)
            if a in g["aldehyde"]:
                inc("1903")
            elif a in g["ketone"]:
                inc("1904")
            elif t == "acid":
                inc("1905")
            elif t == "ester":
                inc("1906")
            elif t == "peracid":
                inc("1907")
    m["1902"] = len(vicinal_pairs)
    # k20: carboxylic acids next to another carbonyl
    for a, t in g["acyl_type"].items():
        if t != "acid":
            continue
        for j in nbrs[a]:
            if sym[j] != "C" or btype[(a, j)] != single:
                continue
            if j in g["aldehyde"]:
                inc("2001")
            elif j in g["ketone"]:
                inc("2002")
            elif g["acyl_type"].get(j) == "ester":
                inc("2003")
    return m


def _evaporation_counts_paper(mol, Chem, maps, g, nbrs, btype):
    """Paper-text counting semantics (Table 4 with footnotes).

    Differences to the upstream mode are documented in the method
    docstring of :func:`evaporation_fragmentation`; every rule is
    anchored on a reproduced paper worked example (D1-D10).
    """
    sym, nh, sp3 = maps["symbol"], maps["nH"], maps["sp3"]
    single = Chem.BondType.SINGLE
    ester_acyl = {c for c, t in g["acyl_type"].items() if t == "ester"}
    pan_acyl = {c for c, t in g["acyl_type"].items() if t == "pan"}
    # c2: carbons plus "in-chain" oxygens (Eq. 26): O whose removal
    # breaks the carbon skeleton into two carbon-bearing pieces, or
    # ring O with a carbon neighbour. Acid OH, peroxy and nitrate
    # oxygens fail the test; a peroxide/ester bridge counts twice.
    in_chain = 0
    for i, s in sym.items():
        if s != "O" or nh[i] != 0:
            continue
        frags = separate_without_atom(mol, i)
        if sum(1 for f in frags if any(sym[k] == "C" for k in f)) >= 2 \
                or (maps["in_ring"][i]
                    and any(sym[j] == "C" for j in nbrs[i])):
            in_chain += 1
    counts = {"2": sum(1 for s in sym.values() if s == "C") + in_chain}
    # c3: branching (Eq. 26: single C-C bonds per carbon above two)
    # minus the ring number (fractional, minimum ring size)
    branch = sum(max(0, maps["ncc"][i] - 2)
                 for i, s in sym.items() if s == "C")
    ring_no = sum(1.0 / maps["minring"][i]
                  for i in sym if 3 <= maps["minring"][i] <= 8)
    counts["3"] = float(branch) - ring_no
    counts["4"] = len(g["nitrate_n"])
    counts["5"] = len(g["ketald"])
    counts["6"] = len(ester_acyl)      # formates included (unlike upstream)
    counts["7"] = len(pan_acyl)
    counts["8"] = sum(1 for c in g["oh"].values() if sp3[c])
    counts["9"] = sum(1 for t in g["acyl_type"].values() if t == "acid")
    counts["10"] = len(g["ooh_c_side"])
    counts["11"] = sum(1 for t in g["acyl_type"].values() if t == "peracid")
    # typed k12: ring oxygens (per atom), ring O-O pairs (per pair),
    # ring ketone/lactone carbonyl carbons, ring C-OH / C-OOH carbons
    oo_members = {i for pair in g["ring_oo"] for i in pair}
    counts["12lin"] = sum(
        1 for i, s in sym.items()
        if s == "O" and nh[i] == 0 and maps["in_ring"][i]
        and i not in oo_members
        and any(sym[j] == "C" for j in nbrs[i])) + len(g["ring_oo"])
    counts["12CL"] = sum(1 for i in g["ketone"] | ester_acyl
                         if maps["in_ring"][i])
    counts["12HB"] = sum(1 for c in g["oh"].values()
                         if maps["in_ring"][c] and sp3[c]) \
        + sum(1 for o in g["ooh_c_side"] for j in nbrs[o]
              if sym[j] == "C" and maps["in_ring"][j])
    # k13 / k14 / k15
    counts["13"] = 0
    for pair in g["alkene_bonds"]:
        for i in pair:
            if any(sym[j] == "C" and btype[(i, j)] == single
                   and j in g["ketald"] for j in nbrs[i]):
                counts["13"] += 1
                break
    counts["14"] = sum(max(0, maps["ncc"][c] - 1)
                       for c in g["oh"].values() if sp3[c])
    counts["15"] = min(1, sum(
        1 for c in g["oh"].values()
        if sp3[c] and any(sym[j] == "C" and maps["alkene"][j]
                          for j in nbrs[c]
                          if btype[(c, j)] == single)))
    cl_partners = g["ketald"] | ester_acyl | pan_acyl
    # k16/k17: per-carbonyl alpha/beta carbonyl partners; footnote f
    # drops a beta partner whose middle carbon carries a functional
    # group (any non-carbon heavy neighbour on the sp3 middle carbon)
    counts["16"] = counts["17"] = 0
    for x in g["ketald"]:
        for j in nbrs[x]:
            if sym[j] == "C" and btype[(x, j)] == single and j in cl_partners:
                counts["16"] += 1
            if sym[j] != "C" or btype[(x, j)] != single:
                continue
            clean_middle = sp3[j] and not any(
                sym[k] != "C" for k in nbrs[j])
            if not clean_middle:
                continue
            for beta in nbrs[j]:
                if beta != x and sym[beta] == "C" \
                        and btype[(j, beta)] == single and beta in cl_partners:
                    counts["17"] += 1
    # k18: functional oxygens (OH, OOH, nitrate) on the alpha carbon;
    # ethers are not functional groups here (decision O5)
    counts["18"] = 0
    for x in g["ketald"]:
        for j in nbrs[x]:
            if sym[j] != "C" or btype[(x, j)] != single:
                continue
            for o in nbrs[j]:
                if sym[o] != "O" or btype[(j, o)] != single:
                    continue
                if g["oh"].get(o) == j or o in g["ooh_c_side"] \
                        or any(k in g["nitrate_n"] for k in nbrs[o]
                               if k != j):
                    counts["18"] += 1
    # k19: per-hydroxyl neighbourhoods; the partner oxygen may be any
    # single-bonded O on an adjacent sp3 carbon (footnote c counts
    # ethers as functional) plus any adjacent carbonyl carbon
    counts["19"] = 0
    for o, c in g["oh"].items():
        if not sp3[c]:
            continue
        for a in nbrs[c]:
            if sym[a] != "C" or btype[(c, a)] != single:
                continue
            if sp3[a]:
                for p in nbrs[a]:
                    if sym[p] != "O" or btype[(a, p)] != single:
                        continue
                    counts["19"] += 1  # every single-bonded O (fn. c)
            if a in g["ketald"] or a in g["acyl_type"]:
                counts["19"] += 1
    # k20: carboxylic acid alpha carbonyls
    counts["20"] = 0
    for a, t in g["acyl_type"].items():
        if t != "acid":
            continue
        for j in nbrs[a]:
            if sym[j] == "C" and btype[(a, j)] == single \
                    and j in cl_partners:
                counts["20"] += 1
    return counts


def _evaporation_assemble_upstream(m):
    """Submode -> descriptor counts, the exact upstream arithmetic
    (verified bit-for-bit against the live oracle on 75 molecules)."""
    c = {"1": 1}
    c["2"] = (m.get("2", 0) + m.get("201", 0) + m.get("202", 0)
              + 2 * m.get("203", 0) + m.get("7002", 0))
    c["3"] = (2 * (m.get("301", 0) + m.get("303", 0) + m.get("305", 0))
              + (m.get("302", 0) + m.get("304", 0) + m.get("306", 0))
              - (m.get("307", 0) / 3 + m.get("308", 0) / 4
                 + m.get("309", 0) / 5 + m.get("310", 0) / 6
                 + m.get("311", 0) / 7 + m.get("312", 0) / 8))
    c["4"] = m.get("4", 0) + m.get("7002", 0)
    c["5"] = (m.get("501", 0) + m.get("502", 0) + m.get("503", 0)
              + m.get("7001", 0))
    c["6"] = m.get("601", 0) + m.get("602", 0)
    c["7"] = m.get("7", 0)
    c["8"] = (m.get("801", 0) + m.get("802", 0) + m.get("803", 0)
              + m.get("804", 0))
    c["9"] = m.get("9", 0)
    c["10"] = m.get("10", 0)
    c["11"] = m.get("11", 0)
    c["12"] = (m.get("502", 0) + m.get("602", 0) + m.get("1201", 0)
               + m.get("1202", 0) + m.get("1203", 0) + m.get("1204", 0)
               + m.get("1205", 0))
    c["12lin"] = m.get("1201", 0) + m.get("1202", 0) + m.get("1204", 0)
    c["12CL"] = m.get("1202", 0) + m.get("502", 0)
    c["12HB"] = m.get("1203", 0) + m.get("1205", 0)
    c["13"] = m.get("13", 0)
    c["14"] = m.get("802", 0) + 2 * m.get("803", 0)
    c["15"] = min(1, m.get("15", 0))
    c["16"] = (m.get("1601", 0) + m.get("1602", 0) + m.get("1603", 0)
               + m.get("1604", 0))
    c["17"] = (m.get("1701", 0) + m.get("1702", 0) + m.get("1703", 0)
               + m.get("1704", 0))
    c["18"] = m.get("1801", 0) + m.get("1802", 0)
    c["19"] = (m.get("1901", 0) + 2 * m.get("1902", 0) + m.get("1903", 0)
               + m.get("1904", 0) + m.get("1905", 0) + m.get("1906", 0)
               + m.get("1907", 0) + m.get("1908", 0))
    c["20"] = (m.get("2001", 0) + m.get("2002", 0) + m.get("2003", 0)
               + m.get("2004", 0))
    return c


def evaporation_fragmentation(mol, upstream_compat=False):
    """Count the 20 EVAPORATION descriptors of a normalised molecule.

    Two counting conventions live behind one interface:

    * ``upstream_compat=False`` (default) -- the paper-text semantics.
      c2 counts in-chain oxygens (Eq. 26), c3 is the continuous
      branching minus the fractional ring number, c4/c7 count nitrates
      and PANs (charge form), c6 includes formates, k12 is typed
      (12lin/12CL/12HB; a lactone is 2x, ring peroxides 12lin only),
      k16/k17 count per carbonyl with footnote f, k18 takes only
      functional oxygens (OH/OOH/nitrate, not ethers, decision O5),
      k19 takes every single-bonded O on the adjacent sp3 carbon
      (footnote c: ethers included). The returned dict leaves the
      plain "12" key at zero; k12 enters through its typed split.
    * ``upstream_compat=True`` -- the UManSysProp semantics, verified
      against the live oracle on the 26 Example_output goldens and 49
      control molecules: c3 through subset-matching branching patterns
      (the binomial match-count model), formate bridge oxygens
      invisible, nitrate/PAN SMARTS dead (c4 = c7 = 0 on RDKit charge
      normalisation), no footnote f, beta interactions strictly
      1,3-dicarbonyl. The dict additionally carries the total "12".

    Parameters
    ----------
    mol : Mol
        Normalised RDKit molecule (charge-form nitrates/PANs, as
        produced by :func:`aeroprops.structure.mol_from_smiles`).
    upstream_compat : bool, optional
        Select the counting convention (see above).

    Returns
    -------
    counts : dict
        Non-zero descriptor counts ("1" is implied; both modes key
        "2".."20" and "12lin"/"12CL"/"12HB"). Empty on failure.
    success : bool
    status : str

    Examples
    --------
    >>> counts, ok, status = evaporation_fragmentation(
    ...     mol_from_smiles("CC(C)(C)C"))
    >>> ok, status, counts
    (True, 'OK', {'2': 5, '3': 2.0})
    >>> evaporation_fragmentation(mol_from_smiles("c1ccccc1"))[1:]
    (False, 'aromatic system outside the EVAPORATION data set')
    """
    loaded = load_rdkit()
    if loaded is None:
        return {}, False, "RDKit unavailable"
    Chem = loaded[0]
    if mol is None:
        return {}, False, "Failed to construct mol"
    unsupported = sorted(
        {a.GetSymbol() for a in mol.GetAtoms()} - _EVAPORATION_ELEMENTS)
    if unsupported:
        return {}, False, "unsupported element(s): " + ", ".join(unsupported)
    if not upstream_compat and any(
            a.GetIsAromatic() for a in mol.GetAtoms()):
        return {}, False, "aromatic system outside the EVAPORATION data set"
    maps = _evaporation_atom_maps(mol, Chem)
    g, nbrs, btype = _evaporation_groups(mol, Chem, maps)
    if not upstream_compat:
        stray_n = {i for i, s in maps["symbol"].items() if s == "N"} \
            - g["nitrate_n"] - g["pan_n"]
        if stray_n:
            return {}, False, "nitrogen outside nitrate/PAN environments"
        counts = _evaporation_counts_paper(
            mol, Chem, maps, g, nbrs, btype)
    else:
        m = _evaporation_counts_upstream(mol, Chem, maps, g, nbrs, btype)
        counts = _evaporation_assemble_upstream(m)
    keep = {k: v for k, v in counts.items() if v}
    keep.pop("1", None)  # c1 is implied; the estimator enforces it
    return dict(sorted(keep.items())), True, "OK"


class EVAPORATION:
    """EVAPORATION sub-cooled liquid vapour pressure (Compernolle
    et al., ACP 11, 9431-9450, 2011).

    Wraps :func:`evaporation_fragmentation` and the pure estimators of
    :mod:`aeroprops.methods.evaporation` in the same idiom as
    :class:`SIMPOL1`: fragmentation runs once in the constructor and
    every estimation method returns ``None`` if it failed.

    Parameters
    ----------
    mol : Mol, optional
        Normalisation (including charge-form nitrate/PAN conversion)
        is applied internally.
    smi : str, optional
        SMILES string, used when `mol` is not given.
    upstream_compat : bool, optional
        Use the UManSysProp counting convention and assembly
        arithmetic (see :func:`evaporation_fragmentation`); default is
        the paper-text semantics.

    Attributes
    ----------
    counts : dict
        Non-zero descriptor counts (empty on failure).
    success : bool
    status : str

    Examples
    --------
    >>> eva = EVAPORATION(smi="CC(C)(C(=O)O)C(C(=O)O)CC(=O)O")  # MBTCA
    >>> eva.success, eva.counts
    (True, {'2': 8, '3': 3.0, '9': 3})
    >>> round(eva.log10_vapour_pressure(298.15), 4)  # log10 Pa
    -6.5127

    The upstream convention (bit-for-bit the UManSysProp result,
    here differing through the branching descriptor):

    >>> up = EVAPORATION(smi="CC(C)(C)C", upstream_compat=True)
    >>> round(up.log10_vapour_pressure(298.15), 4)
    6.0153

    The boiling point solves Eq. 13 at 1 atm (Eq. 14 basis; no Tb
    enters the method as an input):

    >>> round(EVAPORATION(smi="CC(C)(C)C").boiling_point(), 1)
    283.6

    Unsupported input fails softly:

    >>> bad = EVAPORATION(smi="CCCl")
    >>> bad.success, bad.counts, bad.vapour_pressure() is None
    (False, {}, True)
    """

    def __init__(self, mol=None, smi=None, upstream_compat=False):
        if mol is None:
            mol = None if smi is None else mol_from_smiles(smi)
        else:
            mol = normalize_mol(mol)
        if mol is None:
            self.counts, self.success, self.status = \
                {}, False, "Failed to construct mol"
            return
        self._upstream_compat = upstream_compat
        self.counts, self.success, self.status = evaporation_fragmentation(
            mol, upstream_compat=upstream_compat)

    def log10_vapour_pressure(self, T=298.15):
        """Base-10 logarithm of the vapour pressure [log10 Pa].

        >>> round(EVAPORATION(smi="O=C(O)CCCCC(=O)O").log10_vapour_pressure(), 4)
        -3.467
        """
        if not self.success:
            return None
        return evaporation_log10_vapour_pressure(
            self.counts, T, self._upstream_compat)

    def vapour_pressure(self, T=298.15):
        """Vapour pressure [Pa]."""
        log10p = self.log10_vapour_pressure(T)
        return None if log10p is None else 10.0 ** log10p

    def boiling_point(self):
        """Boiling temperature [K] (p0 = 1 atm, Eq. 14 basis).

        Returns ``None`` on failure or when no positive solution of
        Eq. 13 exists.

        >>> round(EVAPORATION(smi="CCO").boiling_point(), 1)
        350.6
        """
        if not self.success:
            return None
        return evaporation_boiling_temperature(
            self.counts, self._upstream_compat)


# ---------------------------------------------------------------------------
# Nannoolal Part 3 vapour pressure: fragmentation layer
# ---------------------------------------------------------------------------
# Same fragmenter engine as the Part 1 boiling-point module
# (aeroprops.boiling_point): the priority-claim partition of
# fragmentation.PriorityFragmenter with a Part 3 catalogue. The SMARTS
# patterns are the ones authored clean-room from the Part 1/Part 3 group
# definitions (both papers share the group system); the Part 3 deltas
# are called out inline and are exactly three:
#
# * 118 is a *main* group here (1,2-diketone, PR 1, valueless - "no
#   estimation possible"): Part 1 used ID 118 for the conjugation
#   correction, which Part 3 renumbers to 134;
# * aromatic bridge carbons split off the fused family: 18 keeps the
#   three-aromatic-bond (fused) atoms, new 214 claims the biphenyl-style
#   bridge carbon (two aromatic bonds + one single bond to an aromatic
#   carbon);
# * priorities are Part 3's own (drawn from methods.nannoolal_vp
#   automatically below).

from .fragmentation import PriorityFragmenter
from .methods.nannoolal_vp import (
    NANNOOLAL_VP_GROUPS,
    nannoolal_vp_boiling_point,
    nannoolal_vp_dB,
    nannoolal_vp_log10_vapour_pressure,
)

#: SMARTS catalogue: group ID -> (pattern, claim); see the Part 1
#: catalogue (aeroprops.boiling_point.NANNOOLAL_BP_SMARTS) for the
#: pattern-by-pattern rationale. Groups the paper defines without a dB
#: value (106 imidazole, 108 thiocyanate, 118 1,2-diketone) are kept:
#: they are counted, then refused at estimation.
NANNOOLAL_VP_SMARTS = {
    # --- Carbon skeleton, chain (paper Table 1) ----------------------
    1: ("[CH3;!R;!$([CH3]~[a]);!$([CH3]~[#7,#8,#9,#17])]", (0,)),
    2: ("[CH3;!R;$([CH3]~[#7,#8,#9,#17])]", (0,)),
    3: ("[CH3;!R;$([CH3]~[a])]", (0,)),
    4: ("[CH2;!R]", (0,)),
    5: ("[CH;!R]", (0,)),
    6: ("[C;H0;!R]", (0,)),
    7: ("[C;!R;$([C]~[#7,#8,#9,#17])]", (0,)),
    8: ("[C;!R;$([C]~[c])]", (0,)),
    9: ("[CH2;R]", (0,)),
    10: ("[CH;R]", (0,)),
    11: ("[C;H0;R]", (0,)),
    12: ("[$([C;R]-!@[#7,#8]),$([C;R]~[#9,#17])]", (0,)),
    13: ("[$([C;R]-@[#7,#8])]", (0,)),
    14: ("[C;R;$([C]~[c])]", (0,)),
    # Aromatic carbons; the Part 3 split: 18 requires three *aromatic*
    # bonds (fused, naphthalene), 214 the bridge carbon of biphenyl /
    # m-terphenyl (two aromatic bonds + single bond to an aromatic C).
    15: ("[cH]", (0,)),
    16: ("[cH0;!$([cH0]~[#7,#8,#9,#17])]", (0,)),
    17: ("[cH0;$([cH0]~[#7,#8,#9,#17])]", (0,)),
    18: ("[cH0](:[c])(:[c]):[c]", (0,)),
    214: ("[cH0](:[c])(:[c])-[c]", (0,)),
    # --- Halogens (family patterns verified against Tables 15c/15h) ---
    19: ("[F][C]", (0,)),
    20: ("[F;$([F][C]=[C])]", (0,)),
    21: ("[F;$([F][C;H0]~[F])]", (0,)),
    102: ("[F;$([F][C;H0]~[Cl])]", (0,)),
    22: ("[F;$([F][C;H1]([F,Cl])[!F,!Cl])]", (0,)),
    23: ("[F;$([F][C;H1]([F,Cl])[F,Cl])]", (0,)),
    25: ("[Cl][C]", (0,)),
    26: ("[Cl;$([Cl][C]~[F,Cl])]", (0,)),
    27: ("[Cl;$([Cl][C]([F,Cl])[F,Cl])]", (0,)),
    28: ("[Cl][c]", (0,)),
    29: ("[Cl;$([Cl][C]=[C])]", (0,)),
    30: ("[Br][C]", (0,)),
    31: ("[Br][c]", (0,)),
    32: ("[I][#6]", (0,)),
    # --- Oxygen ------------------------------------------------------
    # 36 = OH short chain (<C5), applied by post-processing below; the
    # pattern claims any remaining alcohol O (methanol).
    33: ("[OX2H;$([OX2H][C;H0])]", (0,)),
    34: ("[OX2H;$([OX2H][C;H1]([#6])[#6])]", (0,)),
    35: ("[OX2H;$([OX2H][C;H2][#6])]", (0,)),
    36: ("[OX2H;$([OX2H][C])]", (0,)),
    37: ("[OX2H;$([OX2H][c])]", (0,)),
    38: ("[OX2;$([OX2]([#6])[#6])]", (0,)),
    39: ("[OX2;r3]([C;r3])[C;r3]", (0, 1, 2)),
    44: ("[CX3](=[OX1])([OX2H1])[#6]", (0, 1, 2)),
    45: ("[CX3](=[OX1])([OX2;!R][#6])[#6]", (0, 1, 2)),
    47: ("[CX3](=[OX1])([OX2;R;$([OX2]~[#6;!$([#6]=[OX1])])])[#6]",
         (0, 1, 2)),
    46: ("[CX3;H1](=[OX1])([OX2][#6])", (0, 1, 2)),
    48: ("[NX3;H0]([#6])([#6])[CX3](=[OX1])", (0, 3, 4)),
    49: ("[NX3;H1]([#6])[CX3](=[OX1])", (0, 2, 3)),
    50: ("[NX3;H2][CX3](=[OX1])", (0, 1, 2)),
    99: ("[CX3](=[OX1])([NX3])[OX2]", (0, 1, 2, 3)),
    100: ("[CX3](=[OX1])([NX3])[NX3]", (0, 1, 2, 3)),
    51: ("[CX3](=[OX1])([C])[C]", (0, 1)),
    92: ("[CX3](=[OX1])([c])[#6]", (0, 1)),
    52: ("[CX3;H1](=[OX1])([C])", (0, 1)),
    90: ("[CX3;H1](=[OX1])([c])", (0, 1)),
    94: ("[#6][OX2][OX2][#6]", (1, 2)),
    75: ("[N;X2]([OX2H1])=[#6]", (0, 1)),
    77: ("[CX3](=[OX1])([Cl])[#6]", (0, 1)),
    76: ("[CX3](=[OX1])([OX2;!R][CX3](=[OX1]))", (0, 1, 2, 3, 4)),
    96: ("[CX3](=[OX1])([OX2;R][CX3](=[OX1]))", (0, 1, 2, 3, 4)),
    79: ("[CX3](=[OX1])([OX2;!R][#6])([OX2][#6])", (0, 1, 2, 4)),
    103: ("[CX3](=[OX1])([OX2;R])([OX2;R])", (0, 1, 2, 3)),
    72: ("[#6][OX2][N+](=[OX1])([OX1-])", (1, 2, 3, 4)),
    74: ("[#6][OX2][N;X2](=[OX1])", (1, 2, 3)),
    80: ("[#6][N;X2]=[CX2]=[OX1]", (1, 2, 3)),
    65: ("[o](:[c]):[c]", (0,)),
    115: ("[o;r5]~[n;r5]", (0, 1)),
    # Part 3 main group 118: the 1,2-diketone (O=C-C=O, 2,3-butanedione
    # example; PR 1 claims before every carbonyl family). Valueless:
    # counted, then refused ("no estimation possible").
    118: ("[CX3](=[OX1])-[CX3](=[OX1])", (0, 1, 2, 3)),
    # Compernolle (2010) Table 3 extensions (engine-local order below:
    # only 302-before-301 is load-bearing).
    301: ("[OX2H1][OX2][#6]", (0, 1)),
    302: ("[OX2H1][OX2][CX3](=[OX1])", (0, 1, 2, 3)),
    303: ("[#6][CX3](=[OX1])([OX2][OX2][N+](=[OX1])([OX1-]))",
          (1, 2, 3, 4, 5, 6, 7)),
    # --- Nitrogen ----------------------------------------------------
    # 97 = cyclic secondary amine (Part 3's own revision of the Table 1
    # text; morpholine/pyrrolidine examples). 106 (imidazole) and 108
    # (thiocyanate) are valueless: counted, then refused.
    40: ("[NX3;H2][#6]", (0,)),
    41: ("[NX3;H2][c]", (0,)),
    42: ("[#7X3;H1;!R]([#6])[#6]", (0,)),
    97: ("[#7X3;H1;R]([#6;R])[#6;R]", (0,)),
    43: ("[NX3;H0]([#6])([#6])[#6]", (0,)),
    101: ("[N+;X4;H0]([#6])([#6])([#6])[#6]", (0,)),
    66: ("[n;r5;H0]", (0,)),
    67: ("[n;r6;H0]", (0,)),
    57: ("[#6][CX2]#[NX1]", (1, 2)),
    108: ("[#16][CX2]#[NX1]", (1, 2)),
    111: ("[#7][CX2]#[NX1]", (1, 2)),
    91: ("[NX2]=[#6]", (0,)),
    68: ("[C][N+](=[OX1])([OX1-])", (1, 2, 3)),
    69: ("[c][N+](=[OX1])([OX1-])", (1, 2, 3)),
    106: ("[n;r5]1[c;r5][n;r5][c;r5][c;r5]1", (0, 1, 2, 3, 4)),
    # --- Sulfur ------------------------------------------------------
    53: ("[SX2H1][#6]", (0,)),
    54: ("[SX2]([#6])[#6]", (0,)),
    55: ("[#6][SX2][SX2][#6]", (1, 2)),
    56: ("[s]", (0,)),
    82: ("[SX4](=[OX1])(=[OX1])([#6])[#6]", (0, 1, 2)),
    104: ("[SX4](=[OX1])(=[OX1])([OX2][#6])([OX2][#6])", (0, 1, 2, 3, 5)),
    105: ("[SX4](=[OX1])(=[OX1])([NX3])[#6]", (0, 1, 2)),
    107: ("[S+;X3]([OX1-])([#6])[#6]", (0, 1)),
    # --- Unsaturated carbon (bond-anchored) --------------------------
    87: ("[C]=[C]=[C]", (0, 1, 2)),
    88: ("[C;R]=[C;R][C;R]=[C;R]", (0, 1, 2, 3)),
    89: ("[C;!R]=[C;!R][C;!R]=[C;!R]", (0, 1, 2, 3)),
    95: ("[C;X2]#[C;X2][C;X2]#[C;X2]", (0, 1, 2, 3)),
    64: ("[CH1]#[C]", (0, 1)),
    63: ("[C]#[C]", (0, 1)),
    61: ("[CH2]=[C]", (0, 1)),
    60: ("[$([C]~[#7,#8,#9,#17]),$([C]=[C]~[#7,#8,#9,#17]);!R]"
         "=[C;!R]", (0, 1)),
    59: ("[$([C]~[c]),$([C]=[C]~[c]);!R]=[C;!R]", (0, 1)),
    58: ("[C;!R;X2,X3]=[C;!R;X2,X3]", (0, 1)),
    62: ("[C;R]=[C;R]", (0, 1)),
}

#: Engine-local claiming order for the extension groups (no priority
#: numbers exist for them); identical to the Part 1 layer's.
NANNOOLAL_VP_EXTENSION_PRIORITIES = {303: 1, 302: 2, 301: 3}

#: Elements the aeroprops catalogue parameterises (same scope decision
#: as the Part 1 layer: the Si/Ge/Sn/P/As/Se/B/Al families 70-93, 113,
#: 116, 117 and the Si-halogen correction 217 stay untranscribed).
NANNOOLAL_VP_ELEMENTS = frozenset("C N O S F Cl Br I".split())

_NVP_ENGINE = None


def _nannoolal_vp_engine():
    """Build (once) the priority fragmenter for the Part 3 catalogue."""
    global _NVP_ENGINE
    if _NVP_ENGINE is None:
        catalogue = {}
        for gid, (smarts, claim) in NANNOOLAL_VP_SMARTS.items():
            priority = NANNOOLAL_VP_GROUPS[gid][2]
            if priority is None:
                priority = NANNOOLAL_VP_EXTENSION_PRIORITIES[gid]
            catalogue[gid] = {
                "smarts": smarts, "priority": priority, "claim": claim,
            }
        _NVP_ENGINE = PriorityFragmenter(catalogue)
    return _NVP_ENGINE


#: Carbonyl groups whose acyl neighbours decide corrections 134/119/120.
#: Part 3 renumbers Part 1's conjugation correction 118 -> 134; the
#: 1,2-diketone itself is main group 118 (refused, never corrected).
_NVP_HO_CONJUGATION_GIDS = (44, 45, 46, 47, 51, 52, 90, 92)
_NVP_HALO_CONJUGATION_GIDS = (44, 45, 46, 47, 48, 49, 50, 51, 52, 77, 90, 92)


def _nvp_carbonyl_carbon(mol, claim):
    """Index of the C=O carbon inside a carbonyl group's claim tuple."""
    for idx in claim:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetSymbol() != "C":
            continue
        for bond in atom.GetBonds():
            other = bond.GetOtherAtom(atom)
            if other.GetSymbol() == "O" and bond.GetBondTypeAsDouble() == 2.0:
                return idx
    return None


def _nvp_halogen_count(atom):
    """Number of directly bonded fluorine/chlorine atoms."""
    return sum(
        1 for nb in atom.GetNeighbors() if nb.GetSymbol() in ("F", "Cl")
    )


def _nvp_carbon_neighbours(atom):
    return [nb for nb in atom.GetNeighbors() if nb.GetSymbol() == "C"]


def _nvp_is_sp3(atom):
    """True when every bond of `atom` is a single non-aromatic bond."""
    return not atom.GetIsAromatic() and all(
        bond.GetBondTypeAsDouble() == 1.0 for bond in atom.GetBonds()
    )


def _nvp_carbonyl_corrections(mol, assignments, corrections):
    """Corrections 134/119/120: conjugated and halo-substituted carbonyls.

    134 - a carbonyl (acid, ester, lactone, ketone, aldehyde) whose acyl
           carbon neighbour is sp2/aromatic (Part 3's renumbering of the
           Part 1 correction 118; paper Table 2, benzaldehyde example).
    119 - a carbonyl with one acyl neighbour carrying two or more F/Cl.
    120 - a carbonyl with two or more such neighbours (mutually
          exclusive with 119, paper Table 15c).
    """
    for gid in _NVP_HO_CONJUGATION_GIDS:
        for claim in assignments.get(gid, ()):
            cx = _nvp_carbonyl_carbon(mol, claim)
            if cx is None:
                continue
            for nb in mol.GetAtomWithIdx(cx).GetNeighbors():
                if nb.GetIdx() in claim and nb.GetSymbol() == "O":
                    continue
                if nb.GetIsAromatic() or any(
                    b.GetBondTypeAsDouble() == 2.0 for b in nb.GetBonds()
                ):
                    corrections[134] = corrections.get(134, 0) + 1
                    break
    for gid in _NVP_HALO_CONJUGATION_GIDS:
        for claim in assignments.get(gid, ()):
            cx = _nvp_carbonyl_carbon(mol, claim)
            if cx is None:
                continue
            halo_carbons = sum(
                1
                for nb in mol.GetAtomWithIdx(cx).GetNeighbors()
                if nb.GetSymbol() == "C" and _nvp_halogen_count(nb) >= 2
            )
            if halo_carbons >= 2:
                corrections[120] = corrections.get(120, 0) + 1
            elif halo_carbons == 1:
                corrections[119] = corrections.get(119, 0) + 1


def _nvp_halogen_carbon_corrections(mol, corrections):
    """Corrections 121/122: heavily halogenated carbon atoms.

    121 - a carbon with three or more F/Cl neighbours (per carbon).
    122 - a secondary carbon with exactly two F/Cl neighbours.
    """
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != "C":
            continue
        halo = _nvp_halogen_count(atom)
        if halo >= 3:
            corrections[121] = corrections.get(121, 0) + 1
        elif halo == 2 and len(_nvp_carbon_neighbours(atom)) == 2:
            corrections[122] = corrections.get(122, 0) + 1


def _nvp_steric_corrections(mol, corrections):
    """Corrections 130-133: crowded quaternary carbons and C-C bonds.

    130 - a carbon with four carbon neighbours at least one of which is
          aromatic or double-bonded (tert-butylbenzene crowding).
    131/132/133 - per single C(sp3)-C(sp3) bond, by the number of
          further carbon neighbours of each endpoint: 2-2, 3-2, 3-3.
          Ring bonds qualify (ethyl bornyl ether's 132 sits on the
          bornane bridge bond).
    """
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != "C":
            continue
        nbs = _nvp_carbon_neighbours(atom)
        if len(nbs) != 4:
            continue
        if any(
            nb.GetIsAromatic()
            or any(b.GetBondTypeAsDouble() == 2.0 for b in nb.GetBonds())
            for nb in nbs
        ):
            corrections[130] = corrections.get(130, 0) + 1
    seen_bonds = set()
    for bond in mol.GetBonds():
        a1, a2 = bond.GetBeginAtom(), bond.GetEndAtom()
        key = tuple(sorted((a1.GetIdx(), a2.GetIdx())))
        if key in seen_bonds:
            continue
        seen_bonds.add(key)
        if (
            a1.GetSymbol() != "C"
            or a2.GetSymbol() != "C"
            or bond.GetBondTypeAsDouble() != 1.0
            or not (_nvp_is_sp3(a1) and _nvp_is_sp3(a2))
        ):
            continue
        d1 = sum(
            1
            for nb in a1.GetNeighbors()
            if nb.GetSymbol() == "C" and nb.GetIdx() != a2.GetIdx()
        )
        d2 = sum(
            1
            for nb in a2.GetNeighbors()
            if nb.GetSymbol() == "C" and nb.GetIdx() != a1.GetIdx()
        )
        pair = tuple(sorted((d1, d2)))
        gid = {(2, 2): 131, (2, 3): 132, (3, 3): 133}.get(pair)
        if gid is not None:
            corrections[gid] = corrections.get(gid, 0) + 1


def _nvp_ring_corrections(mol, corrections):
    """Corrections 125-129: small rings and aromatic substitution.

    125/126 - one per SSSR ring of size 3-4 / 5 that is not fully
              aromatic.
    127/128/129 - ortho/meta/para: once per molecule, by the closest
              pair of substituted, non-fused carbons of any aromatic
              ring. Part 3 prints no value for the para row (129): it
              stays counted here so the estimator can refuse.
    """
    loaded = load_rdkit()
    if loaded is None:
        return
    rings = [list(ring) for ring in loaded[0].GetSymmSSSR(mol)]
    for ring in rings:
        size = len(ring)
        if size not in (3, 4, 5):
            continue
        if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
            continue
        gid = 125 if size <= 4 else 126
        corrections[gid] = corrections.get(gid, 0) + 1
    ring_count = {}
    for ring in rings:
        for idx in ring:
            ring_count[idx] = ring_count.get(idx, 0) + 1
    closest = None
    for ring in rings:
        if not all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
            continue
        members = set(ring)
        substituted = [
            idx
            for idx in ring
            if ring_count[idx] == 1
            and any(
                nb.GetIdx() not in members
                for nb in mol.GetAtomWithIdx(idx).GetNeighbors()
            )
        ]
        for pos, i1 in enumerate(substituted):
            for i2 in substituted[pos + 1:]:
                d = ring_arc_distance(mol, i1, i2)
                if d in (1, 2, 3) and (closest is None or d < closest):
                    closest = d
    if closest is not None:
        corrections[{1: 127, 2: 128, 3: 129}[closest]] = 1


def _nvp_hydrogen_count_corrections(mol, corrections):
    """Corrections 123/124: molecules with zero/one hydrogen atom."""
    n_h = sum(atom.GetTotalNumHs() for atom in mol.GetAtoms())
    if n_h == 0:
        corrections[123] = 1
    elif n_h == 1:
        corrections[124] = 1


def nannoolal_vp_fragmentation(mol):
    """Fragment `mol` into Nannoolal Part 3 groups (paper Tables 1-2).

    Parameters
    ----------
    mol : Mol
        Normalised RDKit molecule.

    Returns
    -------
    counts : dict
        Group ID -> count, main groups plus second-order corrections
        119-134 plus extension groups 301-303.
    n : int
        Number of non-hydrogen atoms (paper Eq. (7) normalisation: every
        non-H atom counts).
    success : bool
    status : str
        ``"OK"``, or the reason the molecule cannot be parameterised
        (unsupported element, unclaimed atom, failed parse).

    Examples
    --------
    Paper Table 15a (alpha-pinene; both carbonyl-free, so the OH rule
    is idle and the ring/steric corrections carry the slope):

    >>> counts, n, ok, status = nannoolal_vp_fragmentation(
    ...     mol_from_smiles("CC1=CCC2CC1C2(C)C"))
    >>> (ok, status, n, counts)
    (True, 'OK', 10, {1: 3, 9: 2, 10: 2, 11: 1, 62: 1, 132: 2, 125: 1})

    Paper Table 15b (1,2-ethanediol; two carbons < 5, so both OH go to
    the short-chain group 36):

    >>> counts, n, ok, status = nannoolal_vp_fragmentation(
    ...     mol_from_smiles("OCCO"))
    >>> (ok, status, n, counts)
    (True, 'OK', 4, {7: 2, 36: 2})
    """
    if mol is None:
        return {}, 0, False, "Failed to construct mol"
    unsupported = sorted(
        {
            atom.GetSymbol()
            for atom in mol.GetAtoms()
            if atom.GetSymbol() not in NANNOOLAL_VP_ELEMENTS
        }
    )
    if unsupported:
        return {}, mol.GetNumAtoms(), False, (
            "unsupported element(s) %s: outside the C/N/O/S/F/Cl/Br/I "
            "families aeroprops parameterises" % ", ".join(unsupported)
        )
    counts, assignments, success, status = _nannoolal_vp_engine().fragment(mol)
    n = mol.GetNumAtoms()
    counts = dict(sorted(counts.items()))
    if not success:
        return counts, n, False, status

    # Group 36: every alcohol OH of a molecule with fewer than five
    # carbon atoms reclassifies to the short-chain group (paper Table 1;
    # same rule and same claimed atoms as the Part 1 layer).
    n_carbon = sum(
        1 for atom in mol.GetAtoms() if atom.GetSymbol() == "C"
    )
    if n_carbon < 5:
        n_oh = sum(counts.pop(gid, 0) for gid in (33, 34, 35, 36))
        if n_oh:
            counts[36] = n_oh

    corrections = {}
    _nvp_carbonyl_corrections(mol, assignments, corrections)
    _nvp_halogen_carbon_corrections(mol, corrections)
    _nvp_steric_corrections(mol, corrections)
    _nvp_ring_corrections(mol, corrections)
    _nvp_hydrogen_count_corrections(mol, corrections)
    counts.update(corrections)
    return counts, n, True, "OK"


class NannoolalVP:
    """Vapour pressure of pure components (Nannoolal Part 3, 2008).

    Wraps :func:`nannoolal_vp_fragmentation` and the pure estimators of
    :mod:`aeroprops.methods.nannoolal_vp` in the thermo ``Joback``
    idiom (one-time fragmentation in the constructor, ``counts`` /
    ``success`` / ``status`` attributes, estimators return ``None`` on
    failure - never zero). Like :class:`MyrdalYalkowsky` this is an
    *input* method: paper Eq. (6) needs the normal boiling point
    ``Tb``, routed from Nannoolal Part 1 (the default), ``thermo``'s
    Joback method, or experimental data.

    Parameters
    ----------
    mol : Mol, optional
        Normalisation is applied internally; pass any parsed RDKit
        molecule.
    smi : str, optional
        SMILES string, used when `mol` is not given.
    tb : float, optional
        Experimental normal boiling point [K]; overrides `tb_source`
        (recorded as ``"experimental"``).
    tb_source : str, optional
        Where an absent ``tb`` comes from: ``"nannoolal"`` (default),
        ``"joback"`` or ``"experimental"``.

    Attributes
    ----------
    counts : dict
        Non-zero group counts (empty on parse/element failure; partial
        counts when an atom stays unclaimed).
    n : int
        Number of non-hydrogen atoms (0 on parse failure).
    db, db_status
        The Eq. (7) slope and its status, computed once in the
        constructor (``db`` ``None`` when a valueless group or a
        missing interaction parameter refuses the estimate).
    gi_status : str
        Status of the group-interaction term: ``"OK"``, the no-GI flag
        for Compernolle extension groups, or the refusal reason.
    tb, tb_source
        The resolved inputs.
    success : bool
    status : str

    Examples
    --------
    Paper Table 15a (alpha-pinene, experimental Tb 429.0 K; the table
    prints 31.03 kPa at 388.15 K):

    >>> vp = NannoolalVP(smi="CC1=CCC2CC1C2(C)C", tb=429.0)
    >>> (vp.success, vp.db_status, round(vp.db, 7))
    (True, 'OK', 0.1078596)
    >>> round(10.0 ** vp.log10_vapour_pressure(388.15) / 1000.0, 2)
    31.03

    Glycine carries a COOH and an NH2: Table 6 row 218 prints "Do not
    estimate" for the pair, so the estimate is refused - never zeroed:

    >>> vp = NannoolalVP(smi="NCC(=O)O", tb=501.5)
    >>> vp.db is None
    True
    >>> vp.db_status
    'missing interaction parameter (C, M)'
    """

    def __init__(self, mol=None, smi=None, tb=None, tb_source="nannoolal"):
        # set on every path (including failures): the attributes of a
        # thermo-style wrapper always exist, never AttributeError
        self.tb_source = "experimental" if tb is not None else tb_source
        self.counts, self.n = {}, 0
        self.success, self.status = False, "Failed to construct mol"
        self.db, self.db_status, self.gi_status = None, None, None
        self.tb = None
        if mol is None:
            mol = None if smi is None else mol_from_smiles(smi)
        else:
            mol = normalize_mol(mol)
        if mol is None:
            self.mol = None
            return
        self.mol = mol
        self.counts, self.n, self.success, self.status = (
            nannoolal_vp_fragmentation(mol)
        )
        # the Eq. (7) slope needs no Tb: compute it once, whatever the
        # Tb route later does
        if self.success:
            self.db, self.db_status = nannoolal_vp_dB(self.counts, self.n)
            self.gi_status = self.db_status
        self.tb = float(tb) if tb is not None and tb > 0.0 else None
        if tb is not None and self.tb is None:
            self.success, self.status = False, "invalid experimental Tb"
            return
        if not self.success:
            return
        if self.tb is None:
            self.tb, error = _resolve_normal_boiling_point(
                mol, smi, self.tb_source
            )
            if error is not None:
                self.status = error
            if self.tb is None:
                self.success = False
                return

    def log10_vapour_pressure(self, T=298.15):
        """Base-10 logarithm of the vapour pressure [log10 Pa].

        Returns ``None`` on any failure; refuses (with a status) when
        the slope was refused. Emits a warning above the normal boiling
        point, where Eq. 6 extrapolates outside its fitted range.
        """
        if not self.success or self.db is None:
            if self.success and self.db is None:
                self.status = "refused: " + self.db_status
            return None
        if not T > 0.0:
            return None
        if T > self.tb:
            warnings.warn(
                "T = %.1f K exceeds Tb = %.1f K: Eq. 6 extrapolates "
                "beyond the normal boiling point" % (T, self.tb),
                stacklevel=2,
            )
        return nannoolal_vp_log10_vapour_pressure(
            self.counts, self.n, T, self.tb
        )

    def vapour_pressure(self, T=298.15):
        """Vapour pressure [Pa], or ``None`` on any failure."""
        log10p = self.log10_vapour_pressure(T)
        return None if log10p is None else 10.0 ** log10p

    def boiling_point_from_pressure(self, T, p):
        """Normal boiling point [K] from one measured point (Appendix A).

        Inverts Eq. (6) at pressure `p` [Pa] measured at temperature
        `T` [K] (paper Table 15f: 520.85 K from 0.4 kPa at 374.65 K):

        >>> vp = NannoolalVP(smi="CCCOC(=O)CCC(=O)OCCC", tb=520.85)
        >>> round(vp.boiling_point_from_pressure(374.65, 0.4e3), 2)
        520.85
        """
        if not self.success or self.db is None:
            return None
        if not T > 0.0 or not p > 0.0:
            return None
        return nannoolal_vp_boiling_point(self.counts, self.n, T, p)
