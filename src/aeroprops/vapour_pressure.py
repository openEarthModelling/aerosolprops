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
from .structure import load_rdkit, mol_from_smiles, normalize_mol
from .topology import fragments_excluding_bond, separate_without_atom

__all__ = ["SIMPOL1", "simpol1_fragmentation",
           "EVAPORATION", "evaporation_fragmentation"]

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
