# SPDX-License-Identifier: MIT
"""Verification gate 5 for Nannoolal Part 1: live UManSysProp comparison.

Runs the *unmodified* UManSysProp ``boiling_points.nannoolal``
implementation live as the behavioural oracle (UMS's data and SMARTS are
never copied into aerosolprops) over the 26 Example_output.txt SMILES plus
a control battery, and compares group counts and boiling temperatures
with :class:`aerosolprops.boiling_point.NannoolalBP`.

Every Example_output molecule is either bit-for-bit identical to UMS
(EXACT, no extension groups involved) or carries an entry in
``WL_26`` with the live-recorded UMS values, the aerosolprops values and
the paper-based reason for the divergence. UManSysProp deviations from
the paper found this way (paper evidence in the reason strings):

* the N(O)O fragment has no group in the paper set - UMS silently
  drops those atoms (no atom-balance check) and still returns a Tb;
  aerosolprops refuses ("unclaimed atom");
* UMS applies the secondary-alcohol group 34 / fires the aromatic-5-ring
  correction 126 where the paper's OH-family and "non-aromatic ring"
  definitions say otherwise;
* extension groups 301-303: UMS stuffs them into GI classes A/C/F and
  computes interactions; aerosolprops skips GI entirely (no parameters in
  the paper) - identical counts, different Tb;
* several primary groups are missing from the UMS catalogue entirely
  (aromatic N 67 with pyridine as the paper example, bridgehead 18,
  aromatic S 56, amide 48 with DMF as the paper example, urea 100,
  nitrite 74, peroxide 94, carbamate 99, the S family 54/82/107, the
  perfluoro family 21, the methyl of acetonitrile).

Skipped entirely when the UManSysProp source tree is not present
(``AEROSOLPROPS_UMS_REPO``, see :mod:`ums_oracle`).
"""

import pytest

import ums_oracle
from aerosolprops.boiling_point import NannoolalBP

pytestmark = pytest.mark.skipif(
    not ums_oracle.available(), reason="UManSysProp source tree not available"
)


def _ums_counts(ums, smi):
    compound = ums_oracle.readstring(ums, smi)
    raw = dict(ums.groups.nannoolal_primary(compound))
    raw.update(ums.groups.nannoolal_secondary(compound))
    return {k: v for k, v in raw.items() if v}


def _ums_tb(ums, smi):
    return ums.boiling_points.nannoolal(ums_oracle.readstring(ums, smi))


def _ours(smi):
    return NannoolalBP(smi=smi)


# --- the 26 Example_output.txt molecules: exact set --------------------
# No extension groups and no UMS deviation: counts and Tb bit-for-bit.
EXACT_26 = [
    "C1=CC=C(C(=C1)[N+](=O)[O-])O",
    "C1=CC(=CC(=C1)O)[N+](=O)[O-]",
    "C1=CC(=CC=C1[N+](=O)[O-])O",
    "C1=CC=C(C(=C1)C(=O)O)[N+](=O)[O-]",
    "C1=CC(=CC=C1C(=O)O)[N+](=O)[O-]",
    "COO",
    "CC1(C2CCC(O)C1C2)C",
    "C1(C(C(CCO)C1)(C)C)C(O)C",
    "C(O)CO",
]


def test_example_output_exact_bit_for_bit(ums):
    for smi in EXACT_26:
        nnp = _ours(smi)
        assert nnp.success, smi
        assert nnp.counts == _ums_counts(ums, smi), smi
        assert nnp.boiling_point() == _ums_tb(ums, smi), smi


# --- the 26 Example_output.txt molecules: documented deviations --------
# (name, SMILES, UMS counts, UMS Tb, aerosolprops counts or None when the
#  wrapper refuses, aerosolprops Tb or None, reason)
REASON_NOON = (
    "the N(O)O fragment has no group in the paper's set; UMS silently"
    " drops those atoms and still returns a Tb (no atom-balance check),"
    " aerosolprops refuses with 'unclaimed atom'"
)
REASON_OH = (
    "OH family: aerosolprops applies the <5-carbon reassignment to group"
    " 36 (paper Table 1 OH family; the short-chain gate precedes the"
    " neighbour-count classes) and falls back to 36 for gem-diol"
    " oxygens, UMS keeps 33/34 in both situations"
)
REASON_GI = (
    "extension 301 present: aerosolprops skips the group-interaction term"
    " entirely (flagged no-GI; the paper has no GI parameters for"
    " extension groups), UMS stuffs 301 into GI class A and computes"
    " interactions - identical counts, different Tb"
)

WL_26 = [
    # -- UMS drops the N(O)O atoms, aerosolprops refuses ------------------
    (
        "terpene OON(O)O 1",
        "C(C(O)OON(O)O)C1C(C(C(O)C)C1)(C)C",
        {1: 3, 4: 1, 7: 2, 9: 1, 10: 2, 11: 1, 34: 2, 125: 1, 132: 2},
        459.5010232038448, None, None, REASON_NOON,
    ),
    (
        "hydroxyiminoacetaldehydehydrate",
        "C(O)(OON(O)O)CO",
        {7: 2, 34: 1, 36: 1},
        333.20056234333606, None, None, REASON_NOON,
    ),
    (
        "terpene ON(O)O",
        "CC1(C2CCC(CO)(C1C2)ON(O)O)C",
        {1: 2, 7: 1, 9: 3, 10: 2, 11: 1, 12: 1, 35: 1, 125: 1,
         131: 1, 132: 2},
        446.7638910523607, None, None, REASON_NOON,
    ),
    (
        "glycol hydroxylamine adduct",
        "OCCC(O)OON(O)O",
        {4: 1, 7: 2, 34: 1, 36: 1},
        358.3685313338272, None, None, REASON_NOON,
    ),
    (
        "isopropyl hydroxylamine adduct",
        "CC(C(O)OON(O)O)C",
        {1: 2, 5: 1, 7: 1, 34: 1},
        285.0523843468044, None, None, REASON_NOON,
    ),
    (
        "polyol hydroxylamine adduct",
        "CC(O)C(O)CC(O)OON(O)O",
        {1: 1, 4: 1, 7: 3, 34: 3},
        430.2674050180728, None, None, REASON_NOON,
    ),
    (
        "terpene ON(O)O 2",
        "C12C(C(CC(C1(C)O)ON(O)O)C2)(C)C",
        {1: 3, 9: 2, 10: 2, 11: 1, 12: 2, 33: 1, 125: 1, 131: 1,
         132: 2},
        424.2244372058693, None, None, REASON_NOON,
    ),
    (
        "tert hydroperoxy hydroxylamine adduct",
        "CC(C)(C(O)OON(O)O)O",
        {1: 2, 7: 2, 33: 1, 34: 1},
        341.676592674215, None, None, REASON_NOON,
    ),
    (
        "terpene ON(O)O 3",
        "C12C(C(CC(C1(C)OO)ON(O)O)C2)(C)C",
        {1: 3, 9: 2, 10: 2, 11: 1, 12: 2, 125: 1, 131: 1, 132: 2,
         301: 1},
        466.9448047440741, None, None, REASON_NOON,
    ),
    # -- OH family (UMS 33/34 vs aerosolprops 36) -------------------------
    (
        "alpha-hydroxy hydroperoxide",
        "CC(O)OO",
        {1: 1, 7: 1, 34: 1, 301: 1},
        461.85651016792997,
        {1: 1, 7: 1, 36: 1, 301: 1},
        457.7263031668237,
        REASON_OH + " (2 carbons: whole-molecule <C5 reassignment)",
    ),
    (
        "gem-diol on a menthol skeleton",
        "C(C(O)O)C1C(C(C(C)O)C1)(C)C",
        {1: 3, 4: 1, 7: 2, 9: 1, 10: 2, 11: 1, 34: 3, 125: 1, 132: 2},
        574.3724468670971,
        {1: 3, 4: 1, 7: 2, 9: 1, 10: 2, 11: 1, 34: 1, 36: 2, 125: 1,
         132: 2},
        601.9392737001315,
        REASON_OH + " (gem-diol oxygens have one carbon neighbour)",
    ),
    (
        "tetraol hydroperoxide",
        "C(C(C(CO)OO)O)O",
        {7: 4, 34: 1, 36: 2, 301: 1},
        646.8506021054702,
        {7: 4, 36: 3, 301: 1},
        641.4837368270543,
        REASON_OH + " (4 carbons: whole-molecule <C5 reassignment)",
    ),
    # -- identical counts, GI handled differently for extensions -------
    (
        "terpene hydroperoxide 1",
        "C12C(C(CC(C1(C)OO)O)C2)(C)C",
        {1: 3, 9: 2, 10: 2, 11: 1, 12: 2, 34: 1, 125: 1, 131: 1,
         132: 2, 301: 1},
        570.9706450948225,
        {1: 3, 9: 2, 10: 2, 11: 1, 12: 2, 34: 1, 125: 1, 131: 1,
         132: 2, 301: 1},
        564.6462058946723,
        REASON_GI,
    ),
    (
        "terpene hydroperoxide 2",
        "OOC1(CCC2C(C1C2)(C)C)CO",
        {1: 2, 7: 1, 9: 3, 10: 2, 11: 1, 12: 1, 35: 1, 125: 1,
         131: 1, 132: 2, 301: 1},
        589.529486556046,
        {1: 2, 7: 1, 9: 3, 10: 2, 11: 1, 12: 1, 35: 1, 125: 1,
         131: 1, 132: 2, 301: 1},
        583.2050473558958,
        REASON_GI,
    ),
    (
        "hydroperoxy menthol skeleton",
        "C(OO)C1C(C(C(O)C)C1)(C)C",
        {1: 3, 7: 2, 9: 1, 10: 2, 11: 1, 34: 1, 125: 1, 132: 2, 301: 1},
        555.6364088640469,
        {1: 3, 7: 2, 9: 1, 10: 2, 11: 1, 34: 1, 125: 1, 132: 2, 301: 1},
        548.5059133673617,
        REASON_GI,
    ),
    (
        "tertiary hydroperoxy alcohol",
        "OCC(CC)(C)OO",
        {1: 2, 4: 1, 7: 2, 35: 1, 301: 1},
        515.0039341241109,
        {1: 2, 4: 1, 7: 2, 35: 1, 301: 1},
        502.01851343597167,
        REASON_GI,
    ),
    (
        "dihydroperoxy terpene",
        "C(O)C1C(C(C(O)C)(C1)OO)(C)C",
        {1: 3, 7: 2, 9: 1, 10: 1, 11: 1, 12: 1, 34: 1, 35: 1, 125: 1,
         132: 2, 301: 1},
        606.4258253223428,
        {1: 3, 7: 2, 9: 1, 10: 1, 11: 1, 12: 1, 34: 1, 35: 1, 125: 1,
         132: 2, 301: 1},
        596.9391665221175,
        REASON_GI,
    ),
]


@pytest.mark.parametrize(
    "name,smi,ums_counts,ums_tb,ours_counts,ours_tb,reason",
    WL_26,
    ids=[e[0] for e in WL_26],
)
def test_example_output_whitelist(ums, name, smi, ums_counts, ums_tb,
                                  ours_counts, ours_tb, reason):
    """Documented deviations: the live oracle must reproduce the
    recorded UMS values (drift guard) and aerosolprops must reproduce the
    recorded paper-faithful side."""
    assert _ums_counts(ums, smi) == ums_counts, (name, "UMS drifted")
    assert _ums_tb(ums, smi) == ums_tb, (name, "UMS drifted")
    nnp = _ours(smi)
    if ours_counts is None:
        assert not nnp.success and "unclaimed atom" in nnp.status, name
        assert nnp.boiling_point() is None
    else:
        assert nnp.success, (name, nnp.status)
        assert nnp.counts == ours_counts, name
        assert nnp.boiling_point() == ours_tb, name
    assert reason  # every deviation carries a justification


def test_example_output_fully_partitioned(ums):
    """No unexplained behaviour: each of the 26 Example_output SMILES
    is either in the exact set (and identical) or whitelisted."""
    smiles = ums_oracle.example_output_smiles()
    assert len(smiles) == 26
    exact = set(EXACT_26)
    listed = {e[1] for e in WL_26}
    assert len(exact) + len(listed) == 26
    assert not exact & listed
    assert set(smiles) == exact | listed
    for smi in exact:
        nnp = _ours(smi)
        assert nnp.success and nnp.counts == _ums_counts(ums, smi)
        assert nnp.boiling_point() == _ums_tb(ums, smi), smi


# --- control battery outside Example_output ----------------------------
# Molecules whose UMS behaviour was verified identical live (alcohols,
# acyclic amines incl. the 42/97 split, esters, amides, halocarbons).
CONTROL_EXACT = [
    "CO", "CCO", "CCCO", "CCCCO", "OCCO", "Oc1ccccc1", "CCOc1ccccc1",
    "CCOC(C)=O", "O=C1OCCC1", "COC(=O)OC", "CC(=O)C", "O=Cc1ccccc1",
    "CC(N)=O", "CNC(C)=O", "CC(=O)N(C)C",
    "C[N+](=O)[O-]", "O=[N+]([O-])c1ccccc1", "CCO[N+](=O)[O-]",
    "OC(=O)c1ccccc1", "ClCCl", "ClC(Cl)Cl", "Clc1ccccc1",
    "CCNCC", "C1CCNCC1", "C1CCNC1", "C1COCCN1",
    "c1ccccc1Nc1ccccc1", "CNc1ccccc1",
]


def test_control_battery_exact(ums):
    for smi in CONTROL_EXACT:
        nnp = _ours(smi)
        assert nnp.success, smi
        assert nnp.counts == _ums_counts(ums, smi), smi
        assert nnp.boiling_point() == _ums_tb(ums, smi), smi


#: UMS catalogue gaps / deviations on ordinary molecules (live-recorded
#: 2026-11-10). aerosolprops follows the paper in every case; the paper
#: evidence is named per entry.
WL_CONTROL = [
    # (name, SMILES, UMS counts, UMS Tb, ours counts, ours Tb, reason)
    (
        "pyridine",
        "n1ccccc1",
        {15: 5}, 322.5635872766465,
        {15: 5, 67: 1}, 384.6564860636765,
        "UMS misses group 67 (=N(a)(r6)); the paper's own example "
        "molecules for 67 are pyridine and nicotine",
    ),
    (
        "quinoline",
        "n1cccc2ccccc12",
        {15: 7, 17: 1}, 404.172414742249,
        {15: 7, 17: 1, 18: 1, 67: 1}, 512.2970559993364,
        "UMS misses 18 and 67; the paper's examples for 18 "
        "((a)=C(a)2(a)) are naphthalene and quinoline",
    ),
    (
        "naphthalene",
        "c1ccc2ccccc2c1",
        {15: 8}, 386.0738772618168,
        {15: 8, 18: 2}, 504.01465977254634,
        "UMS misses the fused-ring bridgehead group 18 (paper example: "
        "naphthalene)",
    ),
    (
        "thiophene",
        "c1ccsc1",
        {15: 4}, 309.68141237996883,
        {15: 4, 56: 1}, 356.82342190362533,
        "UMS misses aromatic sulfur 56 (paper examples: thiazole, "
        "thiophene)",
    ),
    (
        "thiazole",
        "c1nccs1",
        {15: 3}, 253.34593428497664,
        {15: 3, 56: 1, 66: 1}, 373.993700848729,
        "UMS misses 56 and 66 (=N(a)(r5), paper example: thiazole)",
    ),
    (
        "pyrrole",
        "c1cc[nH]c1",
        {15: 4, 126: 1}, 281.50410862140336,
        {15: 4, 97: 1}, 334.29142683306327,
        "UMS drops the pyrrole NH entirely and fires 126 on the "
        "aromatic five-ring (paper: 126 is a five-membered NON-aromatic "
        "ring); aerosolprops assigns the cyclic secondary amine 97",
    ),
    (
        "acetonitrile",
        "CC#N",
        {57: 1}, 276.2484186801902,
        {1: 1, 57: 1}, 323.5572788185698,
        "UMS drops the methyl; paper group 1 (CH3 (ne)) is a methyl "
        "not bonded to N/O/F/Cl, which the acetonitrile methyl is",
    ),
    (
        "DMF",
        "CN(C)C=O",
        {2: 2}, 194.51005559235523,
        {2: 2, 48: 1}, 424.7538116817152,
        "UMS misses the disubstituted amide 48; the paper's own example "
        "for CON (48) is N,N-dimethylformamide",
    ),
    (
        "urea",
        "NC(=O)N",
        {}, 84.3395,
        {100: 1}, 334.4922830751697,
        "UMS matches nothing (returns the bare base constant); the "
        "carbamide group 100 exists in the paper set",
    ),
    (
        "methyl nitrite",
        "CON=O",
        {2: 1}, 144.6216610184428,
        {2: 1, 74: 1}, 262.93568726778057,
        "UMS misses the nitrite ester group 74",
    ),
    (
        "methyl ethyl sulfide",
        "CSCC",
        {1: 1}, 131.64836013837962,
        {1: 2, 4: 1, 54: 1}, 337.03212419764844,
        "UMS misses the entire sulfur side (54 sulfide, the chain CH2); "
        "also S is not in the (ne) heteroatom list, so both methyls are "
        "group 1 per the paper",
    ),
    (
        "DMSO",
        "CS(C)=O",
        {}, 84.3395,
        {1: 2, 107: 1}, 499.4266988403252,
        "UMS misses the sulfoxide group 107",
    ),
    (
        "dimethyl sulfone",
        "CS(C)(=O)=O",
        {}, 84.3395,
        {1: 2, 82: 1}, 503.1320474252615,
        "UMS misses the sulfone group 82",
    ),
    (
        "perfluoroethane",
        "FC(F)(F)C(F)(F)F",
        {123: 1}, 53.64780647051278,
        {7: 2, 21: 6, 121: 2, 123: 1}, 198.2018610423643,
        "UMS misses the CF3/CF2 family 21 (and 7); the paper worked "
        "example 17c is perfluoro-2-propanone with exactly these "
        "groups",
    ),
    (
        "di-tert-butyl peroxide",
        "CC(C)(C)OOCC(C)(C)",
        {1: 5, 5: 1, 7: 2}, 350.52731737764765,
        {1: 5, 5: 1, 7: 2, 94: 1}, 394.306672735672,
        "UMS misses the peroxide bond group 94",
    ),
    (
        "methyl phenylcarbamate",
        "COC(=O)Nc1ccccc1",
        {2: 1, 15: 5, 17: 1}, 356.25059239874594,
        {2: 1, 15: 5, 17: 1, 99: 1}, 491.95210034873617,
        "UMS misses the carbamate group 99",
    ),
    (
        "glycerol",
        "OCC(O)CO",
        {7: 3, 34: 1, 36: 2}, 552.5838373684409,
        {7: 3, 36: 3}, 572.3905744171899,
        REASON_OH + " (3 carbons: whole-molecule <C5 reassignment)",
    ),
]


@pytest.mark.parametrize(
    "name,smi,ums_counts,ums_tb,ours_counts,ours_tb,reason",
    WL_CONTROL,
    ids=[e[0] for e in WL_CONTROL],
)
def test_control_whitelist(ums, name, smi, ums_counts, ums_tb, ours_counts,
                           ours_tb, reason):
    assert _ums_counts(ums, smi) == ums_counts, (name, "UMS drifted")
    assert _ums_tb(ums, smi) == ums_tb, (name, "UMS drifted")
    nnp = _ours(smi)
    assert nnp.success, (name, nnp.status)
    assert nnp.counts == ours_counts, name
    assert nnp.boiling_point() == ours_tb, name
    assert reason


def test_control_ulp_class(ums):
    """Identical counts with Tb differing only in the last ulp (float
    summation order: UMS sums primary, secondary and GI contributions
    in three separate passes, aerosolprops accumulates in group-ID order -
    both are deterministic, they just round differently)."""
    for smi in ("Nc1ccccc1", "N#Cc1ccccc1", "CC(=O)c1ccccc1"):
        nnp = _ours(smi)
        assert nnp.success, smi
        assert nnp.counts == _ums_counts(ums, smi), smi
        assert abs(nnp.boiling_point() - _ums_tb(ums, smi)) < 1e-9, smi

