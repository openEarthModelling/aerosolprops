# SPDX-License-Identifier: MIT
"""Normal boiling point estimation (object-oriented wrapper layer).

Nannoolal Part 1 (Nannoolal, Rarey, Ramjugernath & Cordes, Fluid Phase
Equilibria 226, 45-63, 2004) via the priority-claim fragmenter of
:mod:`aeroprops.fragmentation` and the pure coefficient catalogue and
estimator of :mod:`aeroprops.methods.nannoolal_bp`.

This module owns everything RDKit: the SMARTS catalogue (authored
clean-room from the paper's group *definitions*; no third-party pattern
files consulted), the claimed-atom sets of every group, and the
second-order corrections 118-133 that no single SMARTS pattern can
express (ring sizes, ortho/meta/para pairs, steric bond counts,
molecular hydrogen count). Functional groups claim every atom of the
functional pattern except their flanking carbon skeleton (paper Table
17d: the ester group covers exactly the C(=O)O atoms); the skeleton
groups claim single carbons; overlaps are adjudicated by the paper's
priority numbers (Table 1, footnote a: lower number wins).

The aeroprops catalogue parameterises the C/N/O/S/F/Cl/Br/I subset of
the paper's groups. Molecules containing other elements the paper
covers (Si, Ge, Sn, P, As, Se, B, Al - groups 70-93, 113, 116, 117)
fail cleanly with a status naming the missing family: atmospheric SOA
chemistry is outside their scope and transcribing them without test
compounds would be guesswork.

"""

from .fragmentation import PriorityFragmenter
from .methods.nannoolal_bp import (
    NANNOOLAL_BP_GROUPS,
    nannoolal_bp_boiling_point,
    nannoolal_bp_gi,
)
from .structure import load_rdkit, mol_from_smiles, normalize_mol
from .topology import ring_arc_distance

__all__ = ["NannoolalBP", "nannoolal_fragmentation"]

#: SMARTS catalogue: group ID -> (pattern, claim). The pattern is
#: authored from the paper's group definition; ``claim`` lists the
#: pattern-atom indices the group instance owns (its anchor first is not
#: required - the engine only needs the set). Every non-hydrogen atom
#: must end up claimed by exactly one instance.
NANNOOLAL_BP_SMARTS = {
    # --- Carbon skeleton, chain (paper Table 1) ----------------------
    # CH3 by neighbour type; (e) = N/O/F/Cl, (a) = aromatic atom of any
    # element. 2 < 3 < 1 by priority; a methyl on an aromatic N or O is
    # both (e) and (a) and goes to the lower number, 2.
    1: ("[CH3;!R;!$([CH3]~[a]);!$([CH3]~[#7,#8,#9,#17])]", (0,)),
    2: ("[CH3;!R;$([CH3]~[#7,#8,#9,#17])]", (0,)),
    3: ("[CH3;!R;$([CH3]~[a])]", (0,)),
    # Chain sp3 carbons of any hydrogen count; the qualified variants 7
    # (electronegative neighbour) and 8 (aromatic carbon neighbour) have
    # lower priority numbers and claim first.
    4: ("[CH2;!R]", (0,)),
    5: ("[CH;!R]", (0,)),
    6: ("[C;H0;!R]", (0,)),
    7: ("[C;!R;$([C]~[#7,#8,#9,#17])]", (0,)),
    8: ("[C;!R;$([C]~[c])]", (0,)),
    # Ring sp3 carbons; 12/13/14 (heteroatom or aromatic neighbours)
    # have lower priority numbers than the plain 9/10/11.
    9: ("[CH2;R]", (0,)),
    10: ("[CH;R]", (0,)),
    11: ("[C;H0;R]", (0,)),
    # 12: ring C with an external (non-ring-bond) N/O, or any F/Cl.
    12: ("[$([C;R]-!@[#7,#8]),$([C;R]~[#9,#17])]", (0,)),
    # 13: ring C with the N/O inside the same ring (ring bond).
    13: ("[$([C;R]-@[#7,#8])]", (0,)),
    14: ("[C;R;$([C]~[c])]", (0,)),
    # Aromatic carbons. 18 (fused, three aromatic neighbours) beats 16
    # and 17 by priority; 17/16 split on electronegative neighbours.
    15: ("[cH]", (0,)),
    16: ("[cH0;!$([cH0]~[#7,#8,#9,#17])]", (0,)),
    17: ("[cH0;$([cH0]~[#7,#8,#9,#17])]", (0,)),
    18: ("[cH0](:[c])(:[c])~[c]", (0,)),
    # --- Halogens ----------------------------------------------------
    # Fluorine: the specific families (20-23, 102) have lower priority
    # numbers than the generic 19, so the patterns can stay positive -
    # whatever no specific family claims falls to 19 (paper: "F
    # connected to C" as the remainder class). The neighbour-hydrogen
    # conditions follow Table 1: H2 carbon "a" variants take exactly one
    # further halogen, H1/H0 "b" variants take at least two/one.
    19: ("[F][C]", (0,)),
    20: ("[F;$([F][C]=[C])]", (0,)),
    21: ("[F;$([F][C;H0]~[F])]", (0,)),
    102: ("[F;$([F][C;H0]~[Cl])]", (0,)),
    22: ("[F;$([F][C;H1]([F,Cl])[!F,!Cl])]", (0,)),
    23: ("[F;$([F][C;H1]([F,Cl])[F,Cl])]", (0,)),
    # Chlorine: 27 (two or more co-halogens) then 26 (one) then 25
    # (none), descending priority; aromatic and vinyl handled first.
    25: ("[Cl][C]", (0,)),
    26: ("[Cl;$([Cl][C]~[F,Cl])]", (0,)),
    27: ("[Cl;$([Cl][C]([F,Cl])[F,Cl])]", (0,)),
    28: ("[Cl][c]", (0,)),
    29: ("[Cl;$([Cl][C]=[C])]", (0,)),
    30: ("[Br][C]", (0,)),
    31: ("[Br][c]", (0,)),
    32: ("[I][#6]", (0,)),
    # --- Oxygen ------------------------------------------------------
    # Alcohols 33/34/35 (carbon neighbour count; paper Table 1); 36 is
    # algorithmic (molecules with fewer than five C: post-processed from
    # 33/34/35); 37 is the phenol O. 44-52 carbonyl families claim the
    # C, the Oxo O and the O/N single-bonded heteroatom, but never the
    # flanking carbon skeleton (paper Table 17d atom lists).
    33: ("[OX2H;$([OX2H][C;H0])]", (0,)),
    34: ("[OX2H;$([OX2H][C;H1]([#6])[#6])]", (0,)),
    35: ("[OX2H;$([OX2H][C;H2][#6])]", (0,)),
    # 36 claims any remaining alcohol O (e.g. methanol, whose carbon has
    # no C neighbour at all); the <C5 reclassification below moves every
    # 33/34/35/36 instance into 36 for short-chain molecules (paper
    # Table 1, note ruling: the short-chain gate precedes the neighbour
    # count - di-isopropanolamine, Table 17b, keeps 34 at six C).
    36: ("[OX2H;$([OX2H][C])]", (0,)),
    37: ("[OX2H;$([OX2H][c])]", (0,)),
    # 38: plain ether O (two C neighbours, any ring/aromatic status);
    # furan-type O is aromatic ([o], group 65) and never matches [OX2].
    38: ("[OX2;$([OX2]([#6])[#6])]", (0,)),
    # 39: epoxide - the 3-ring O plus both its ring carbons (priority
    # beats the plain ether 38).
    39: ("[OX2;r3]([C;r3])[C;r3]", (0, 1, 2)),
    44: ("[CX3](=[OX1])([OX2H1])[#6]", (0, 1, 2)),
    # 45/47 split acyclic/cyclic on the ester single-bond O; 47 demands
    # a non-carbonyl C neighbour so cyclic carbonates stay whole for 103.
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
    # Compernolle (2010) Table 3 extensions; the paper has no priority
    # numbers for them, so aeroprops assigns engine-local ordering (a
    # separate constant below) - 302 before 301 because the peroxy-acid
    # OH-O also matches the hydroperoxide pattern. Nitro/PAN nitrogen is
    # written in the charge-separated form aeroprops' normaliser emits.
    301: ("[OX2H1][OX2][#6]", (0, 1)),
    302: ("[OX2H1][OX2][CX3](=[OX1])", (0, 1, 2, 3)),
    303: ("[#6][CX3](=[OX1])([OX2][OX2][N+](=[OX1])([OX1-]))",
          (1, 2, 3, 4, 5, 6, 7)),
    # --- Nitrogen ----------------------------------------------------
    # Amines split on hydrogen count and aromatic neighbours; 97 keeps
    # the Part 1 text definition (at least one aromatic C neighbour)
    # over the examples table - see module docs. 66/67 are ring-size
    # classes of aromatic N; 106 claims whole imidazole rings first.
    40: ("[NX3;H2][#6]", (0,)),
    41: ("[NX3;H2][c]", (0,)),
    42: ("[NX3;H1]([#6])[#6]", (0,)),
    97: ("[NX3;H1]([#6])[c]", (0,)),
    43: ("[NX3;H0]([#6])([#6])[#6]", (0,)),
    101: ("[N+;X4;H0]([#6])([#6])([#6])[#6]", (0,)),
    66: ("[n;r5]", (0,)),
    67: ("[n;r6]", (0,)),
    57: ("[#6][CX2]#[NX1]", (1, 2)),
    108: ("[#16][CX2]#[NX1]", (1, 2)),
    111: ("[#7][CX2]#[NX1]", (1, 2)),
    91: ("[NX2]=[#6]", (0,)),
    68: ("[#6][N+](=[OX1])([OX1-])", (1, 2, 3)),
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
    # Sulfoxide: RDKit's normaliser renders S=O sulfoxides charge-
    # separated (CS(C)=O -> C[S+](C)[O-]); sulfones keep their double
    # bonds (82 above), so 107 is written against the charge form only.
    107: ("[S+;X3]([OX1-])([#6])[#6]", (0, 1)),
    # --- Unsaturated carbon (bond-anchored; claims every sp2/sp C of
    # the unit).  The C=C families are written so the substitution
    # condition sits on the sp2 carbon itself (59/60: the OR recursion
    # checks the atom's own substituent or its double-bond partner's),
    # which resolves the paper's cinnamic-alcohol example in favour of
    # 59, not 60.
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

#: Engine-local claiming order for the extension groups 301-303 (the
#: paper and Compernolle (2010) give no priority numbers for them).
#: Values are chosen below every competing native group; only the
#: 302-before-301 relation is load-bearing (a peroxy acid otherwise
#: matches the hydroperoxide pattern first).
NANNOOLAL_BP_EXTENSION_PRIORITIES = {303: 1, 302: 2, 301: 3}

#: Elements the aeroprops catalogue parameterises (paper groups outside
#: this set - Si, Ge, Sn, P, As, Se, B, Al families 70-93, 113, 116,
#: 117 - are deliberately not transcribed; see module docstring).
NANNOOLAL_BP_ELEMENTS = frozenset("C N O S F Cl Br I".split())

#: Paper group 109 ("CO connected to N", example "methyl thioacetate"
#: which contains no nitrogen) has no implementable definition and is
#: excluded from the catalogue; real amides are covered by 48-50, 99,
#: 100.  Unclaimed C(=O)-N fragments of exotic acyl species therefore
#: fail the atom balance cleanly.
_ENGINE = None


def _nannoolal_engine():
    """Build (once) the priority fragmenter for the Nannoolal catalogue."""
    global _ENGINE
    if _ENGINE is None:
        catalogue = {}
        for gid, (smarts, claim) in NANNOOLAL_BP_SMARTS.items():
            priority = NANNOOLAL_BP_GROUPS[gid][2]
            if priority is None:
                priority = NANNOOLAL_BP_EXTENSION_PRIORITIES[gid]
            catalogue[gid] = {
                "smarts": smarts, "priority": priority, "claim": claim,
            }
        _ENGINE = PriorityFragmenter(catalogue)
    return _ENGINE


#: Carbonyl groups whose acyl neighbours decide corrections 118-120.
_HO_CONJUGATION_GIDS = (44, 45, 46, 47, 51, 52, 90, 92)
_HALO_CONJUGATION_GIDS = (44, 45, 46, 47, 48, 49, 50, 51, 52, 77, 90, 92)


def _carbonyl_carbon(mol, claim):
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


def _halogen_count(atom):
    """Number of directly bonded fluorine/chlorine atoms."""
    return sum(
        1 for nb in atom.GetNeighbors() if nb.GetSymbol() in ("F", "Cl")
    )


def _carbon_neighbours(atom):
    return [nb for nb in atom.GetNeighbors() if nb.GetSymbol() == "C"]


def _is_sp3(atom):
    """True when every bond of `atom` is a single non-aromatic bond."""
    return not atom.GetIsAromatic() and all(
        bond.GetBondTypeAsDouble() == 1.0 for bond in atom.GetBonds()
    )


def _carbonyl_corrections(mol, assignments, corrections):
    """Corrections 118-120: conjugated and halo-substituted carbonyls.

    118 - a carbonyl (acid, ester, lactone, ketone, aldehyde) whose acyl
           carbon neighbour is sp2/aromatic (paper Table 11); aryl
           ketones and benzaldehyde always carry it by definition.
    119 - a carbonyl with one acyl neighbour carrying two or more F/Cl.
    120 - a carbonyl with two or more such neighbours (per paper Table
           17c, 119 and 120 are mutually exclusive).
    """
    for gid in _HO_CONJUGATION_GIDS:
        for claim in assignments.get(gid, ()):
            cx = _carbonyl_carbon(mol, claim)
            if cx is None:
                continue
            for nb in mol.GetAtomWithIdx(cx).GetNeighbors():
                if nb.GetIdx() in claim and nb.GetSymbol() == "O":
                    continue
                if nb.GetIsAromatic() or any(
                    b.GetBondTypeAsDouble() == 2.0 for b in nb.GetBonds()
                ):
                    corrections[118] = corrections.get(118, 0) + 1
                    break
    for gid in _HALO_CONJUGATION_GIDS:
        for claim in assignments.get(gid, ()):
            cx = _carbonyl_carbon(mol, claim)
            if cx is None:
                continue
            halo_carbons = sum(
                1
                for nb in mol.GetAtomWithIdx(cx).GetNeighbors()
                if nb.GetSymbol() == "C" and _halogen_count(nb) >= 2
            )
            if halo_carbons >= 2:
                corrections[120] = corrections.get(120, 0) + 1
            elif halo_carbons == 1:
                corrections[119] = corrections.get(119, 0) + 1


def _halogen_carbon_corrections(mol, corrections):
    """Corrections 121/122: heavily halogenated carbon atoms.

    121 - a carbon with three or more F/Cl neighbours (per carbon).
    122 - a carbon with exactly two F/Cl and two carbon neighbours.
    """
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != "C":
            continue
        halo = _halogen_count(atom)
        if halo >= 3:
            corrections[121] = corrections.get(121, 0) + 1
        elif halo == 2 and len(_carbon_neighbours(atom)) == 2:
            corrections[122] = corrections.get(122, 0) + 1


def _steric_corrections(mol, corrections):
    """Corrections 130-133: crowded quaternary carbons and C-C bonds.

    130 - a carbon with four carbon neighbours at least one of which is
          aromatic or double-bonded (tert-butylbenzene crowding; the
          bond from an sp3 carbon to a ring atom is single, so the
          condition sits on the neighbour, not the bond).
    131/132/133 - per *acyclic* single C(sp3)-C(sp3) bond, by the number
          of further carbon neighbours of each endpoint: 2-2, 3-2, 3-3
          (paper Table 2 names C2C-CC2 / C3C-CC2 / C3C-CC3). Ring bonds
          are excluded: the bicyclohexyl junction bond is acyclic even
          though its carbons are ring atoms, while cyclohexane ring
          bonds must not fire (their crowding is parameterised by the
          ring corrections 125/126).
    """
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != "C":
            continue
        nbs = _carbon_neighbours(atom)
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
        if key in seen_bonds or bond.IsInRing():
            continue
        seen_bonds.add(key)
        if (
            a1.GetSymbol() != "C"
            or a2.GetSymbol() != "C"
            or bond.GetBondTypeAsDouble() != 1.0
            or not (_is_sp3(a1) and _is_sp3(a2))
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


def _ring_corrections(mol, corrections):
    """Corrections 125-129: small rings and aromatic substitution patterns.

    125/126 - one per SSSR ring of size 3-4 / 5 that is not fully
              aromatic (cyclopropane/butane/pentane skeleton strain;
              furan-style aromatic rings are excluded).
    127/128/129 - ortho/meta/para: once per molecule, by the closest
              pair of substituted, non-fused carbons of any aromatic
              ring (ring-arc distance 1/2/3). Fused ring atoms carry
              no substituent (naphthalene gets no o/m/p correction);
              1,2,4-trimethylbenzene counts as ortho only.
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


def _hydrogen_count_corrections(mol, corrections):
    """Corrections 123/124: molecules with zero/one hydrogen atom."""
    n_h = sum(atom.GetTotalNumHs() for atom in mol.GetAtoms())
    if n_h == 0:
        corrections[123] = 1
    elif n_h == 1:
        corrections[124] = 1


