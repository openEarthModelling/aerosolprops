# SPDX-License-Identifier: MIT
"""Verification gate 5 for Nannoolal Part 3: live UManSysProp comparison.

Runs the *unmodified* UManSysProp ``vapour_pressures.nannoolal``
implementation live as the behavioural oracle (UMS's data and SMARTS are
never copied into aerosolprops). Protocol = the Example_output.txt
walkthrough: log10(P/atm) at T = 298.15 with Tb taken from UMS's own
``boiling_points.nannoolal`` - the same inputs UMS's example script feeds
the method - so the comparison isolates the Part 3 dB layer.

Every Example_output molecule is either identical to UMS (bit-for-bit or
last-ulp accumulation order), or carries an entry in ``WL_26`` with the
live-recorded UMS and aerosolprops values and the paper-based reason.
UManSysProp deviations from the paper found this way:

* three coefficient misprints in UMS's Part 3 data files: group 7 as
  0.0378487 (paper Table 4 prints 37.8485e-3), group 48 as 0.2692470
  (Table 4 prints 269.2471e-3) and extension 303 as 0.519 (Compernolle
  2010 Table 3 prints 519.03e-3) - aerosolprops follows the printed tables;
* the para-pair group 129 sits in UMS's secondary file with value
  0.0000000; Table 5 has no row for 129, so aerosolprops refuses;
* Table 6 row 218 (C, M) is "Do not estimate" in the paper; UMS's
  interaction matrix zero-fills the pair and estimates amino acids;
* UMS lacks group 214 (aromatic bridge carbon; biphenyl) and the
  1,2-diketone 118 (UMS fragments biacetyl as two plain ketones);
* UMS keeps the Part 1 OH classes (33/34/35) where the Part 3 Table 1
  family reassigns short-chain (<5 C) and gem-diol OHs to 36;
* extension groups: UMS stuffs 301/302/303 into GI classes A/C/F and
  computes interactions; the paper has no GI parameters for extensions,
  so aerosolprops skips the GI term entirely (no-GI branch);
* the N(O)O fragment has no group in the paper set - UMS silently drops
  those atoms (no atom-balance check), aerosolprops refuses;
* UMS's SMARTS catalogue misses the sulfur family and CF3/CF2 group 21
  (Part 1 findings, inherited by the Part 3 data files).

Where UMS and the paper agree, GI terms agree bit-for-bit even for
6-member same-class cases (hexol) and multiple classes (adipic acid in
the exact set). Skipped entirely when the UManSysProp source tree is not
present (``AEROSOLPROPS_UMS_REPO``, see :mod:`ums_oracle`).
"""

import os

import pytest

import ums_oracle
from aerosolprops.methods.nannoolal_vp import LOG10_ATM_TO_PA
from aerosolprops.vapour_pressure import NannoolalVP

pytestmark = [
    pytest.mark.skipif(
        not ums_oracle.available(), reason="UManSysProp source tree not available"
    ),
    # Two whitelist molecules (sulfur/hexafluoro controls) get a Tb below
    # T from UMS's own method, so the wrapper's Eq. 6 extrapolation
    # warning fires; that warning is gate-7 territory, not this file's.
    pytest.mark.filterwarnings("ignore:T = .* extrapolates"),
]

T = 298.15


def _golden():
    """The VP_Nannoolal_BP_Nannoolal column of Example_output.txt."""
    path = os.path.join(
        ums_oracle.USP_REPO, "umansysprop", "examples", "Example_output.txt"
    )
    out = {}
    with open(path) as handle:
        for line in handle.read().splitlines()[1:]:
            if line.strip():
                smiles, value = line.split("\t")[:2]
                out[smiles] = float(value)
    return out


def _ums_parts(ums, smi):
    """(Tb from UMS's own Part 1 method, log10(P/atm) from UMS Part 3)."""
    compound = ums_oracle.readstring(ums, smi)
    tb = ums.boiling_points.nannoolal(compound)
    return tb, ums.vapour_pressures.nannoolal(compound, T, tb)


def _ums_counts(ums, smi):
    compound = ums_oracle.readstring(ums, smi)
    raw = dict(ums.groups.nannoolal_primary(compound))
    raw.update(ums.groups.nannoolal_secondary(compound))
    return {k: v for k, v in raw.items() if v}


def _ours_vp(w):
    """log10(P/atm); None when the wrapper refuses."""
    value = w.log10_vapour_pressure(T)
    return None if value is None else value - LOG10_ATM_TO_PA


def _ours(smi, tb):
    return NannoolalVP(smi=smi, tb=tb, tb_source="experimental")


# --- the 26 Example_output.txt molecules: UMS-side drift lock -----------
def test_ums_reproduces_example_output_column(ums):
    """The live oracle must reproduce the shipped VP_Nannoolal_BP_
    Nannoolal column bit-for-bit (this also pins the protocol: T =
    298.15, Tb = UMS's own nannoolal boiling point)."""
    golden = _golden()
    smiles = ums_oracle.example_output_smiles()
    assert len(smiles) == 26 and set(smiles) == set(golden)
    for smi in smiles:
        _, vp = _ums_parts(ums, smi)
        assert vp == golden[smi], smi


# --- the 26 Example_output.txt molecules: exact set ---------------------
# Counts, dB and the final expression identical to the last bit.
EXACT_26 = [
    "C1=CC(=CC(=C1)O)[N+](=O)[O-]",
    "C1=CC=C(C(=C1)C(=O)O)[N+](=O)[O-]",
    "CC1(C2CCC(O)C1C2)C",
]


def test_example_output_exact_bit_for_bit(ums):
    golden = _golden()
    for smi in EXACT_26:
        tb, ums_vp = _ums_parts(ums, smi)
        assert _ours_vp(_ours(smi, tb)) == ums_vp == golden[smi], smi


def test_example_output_exact_default_route(ums):
    """End-to-end with aerosolprops' own default route (our Part 1 Tb,
    which the Part 1 oracle suite already pins bit-for-bit to UMS's for
    these molecules): the golden column is reproduced without feeding
    UMS anything."""
    golden = _golden()
    for smi in EXACT_26:
        w = NannoolalVP(smi=smi)
        assert w.success and w.tb is not None, (smi, w.status)
        assert _ours_vp(w) == golden[smi], smi


# --- the 26 Example_output.txt molecules: last-ulp class ----------------
# Identical model, different float accumulation order (UMS sums primary,
# secondary and interaction contributions in three separate passes,
# aerosolprops accumulates dB in group-ID order) - a deterministic last-ulp
# rounding difference only.
ULP_26 = [
    "C1=CC=C(C(=C1)[N+](=O)[O-])O",
    "COO",
]


def test_example_output_ulp_class(ums):
    for smi in ULP_26:
        tb, ums_vp = _ums_parts(ums, smi)
        diff = _ours_vp(_ours(smi, tb)) - ums_vp
        assert 0 < abs(diff) < 1e-9, (smi, diff)


# --- the 26 Example_output.txt molecules: documented deviations --------
# (name, SMILES, UMS log10(P/atm), aerosolprops log10(P/atm) or None when
#  refused, refusal-status substring or None, reason)
REASON_129 = (
    "UMS's secondary file carries the para-pair group 129 with value"
    " 0.0000000; Table 5 has no row for 129, so aerosolprops counts it and"
    " then refuses (valueless)"
)
REASON_NOON = (
    "the N(O)O fragment has no group in the paper's set; UMS silently"
    " drops those atoms and still returns a pressure (no atom-balance"
    " check), aerosolprops refuses with 'unclaimed atom'"
)
REASON_OH = (
    "OH family: aerosolprops applies the Part 3 Table 1 short-chain (<5 C)"
    " reassignment to group 36 and the gem-diol fallback, UMS keeps the"
    " Part 1 classes 33/34/35"
)
REASON_GI = (
    "extension 301 present: aerosolprops skips the group-interaction term"
    " entirely (no GI parameters for extension groups), UMS stuffs 301"
    " into GI class A and computes interactions - identical counts,"
    " different dB"
)
REASON_T7 = (
    "counts identical; UMS's Part 3 data file misprints group 7 as"
    " 0.0378487 where Table 4 prints 37.8485e-3, aerosolprops follows the"
    " printed table"
)

WL_26 = [
    (
        "para nitrophenol",
        "C1=CC(=CC=C1[N+](=O)[O-])O",
        -4.250442700439299, None,
        "refused: no dB parameter for group 129", REASON_129,
    ),
    (
        "para nitrobenzoic acid",
        "C1=CC(=CC=C1C(=O)O)[N+](=O)[O-]",
        -7.58580417980072, None,
        "refused: no dB parameter for group 129", REASON_129,
    ),
    (
        "terpene OON(O)O 1",
        "C(C(O)OON(O)O)C1C(C(C(O)C)C1)(C)C",
        -3.8009840780271804, None, "unclaimed atom", REASON_NOON,
    ),
    (
        "hydroxyiminoacetaldehydehydrate",
        "C(O)(OON(O)O)CO",
        -0.7366175484461129, None, "unclaimed atom", REASON_NOON,
    ),
    (
        "terpene ON(O)O",
        "CC1(C2CCC(CO)(C1C2)ON(O)O)C",
        -2.9725199214037374, None, "unclaimed atom", REASON_NOON,
    ),
    (
        "glycol hydroxylamine adduct",
        "OCCC(O)OON(O)O",
        -1.2972155995339851, None, "unclaimed atom", REASON_NOON,
    ),
    (
        "isopropyl hydroxylamine adduct",
        "CC(C(O)OON(O)O)C",
        0.23917168525698887, None, "unclaimed atom", REASON_NOON,
    ),
    (
        "polyol hydroxylamine adduct",
        "CC(O)C(O)CC(O)OON(O)O",
        -3.3803427189506103, None, "unclaimed atom", REASON_NOON,
    ),
    (
        "terpene ON(O)O 2",
        "C12C(C(CC(C1(C)O)ON(O)O)C2)(C)C",
        -2.4903877630281936, None, "unclaimed atom", REASON_NOON,
    ),
    (
        "tert hydroperoxy hydroxylamine adduct",
        "CC(C)(C(O)OON(O)O)O",
        -0.9206870670619538, None, "unclaimed atom", REASON_NOON,
    ),
    (
        "terpene ON(O)O 3",
        "C12C(C(CC(C1(C)OO)ON(O)O)C2)(C)C",
        -3.527517239544599, None, "unclaimed atom", REASON_NOON,
    ),
    # -- OH family and/or GI handling of extension 301 ------------------
    (
        "alpha-hydroxy hydroperoxide",
        "CC(O)OO",
        -3.6759962673833866, -3.827343029971476, None,
        REASON_OH + " (2 carbons: whole-molecule <C5 reassignment); also "
        + REASON_GI,
    ),
    (
        "terpene hydroperoxide 1",
        "C12C(C(CC(C1(C)OO)O)C2)(C)C",
        -6.836867533054316, -6.940793306759183, None, REASON_GI,
    ),
    (
        "terpene hydroperoxide 2",
        "OOC1(CCC2C(C1C2)(C)C)CO",
        -7.332580320435168, -7.444722641758563, None,
        REASON_GI + "; plus the group-7 data misprint",
    ),
    (
        "hydroperoxy menthol skeleton",
        "C(OO)C1C(C(C(O)C)C1)(C)C",
        -6.438228011366474, -6.543595171675318, None,
        REASON_GI + "; plus the group-7 data misprint",
    ),
    (
        "tertiary hydroperoxy alcohol",
        "OCC(CC)(C)OO",
        -5.130551016473907, -5.260768865961851, None,
        REASON_GI + "; plus the group-7 data misprint",
    ),
    (
        "dihydroperoxy terpene",
        "C(O)C1C(C(C(O)C)(C1)OO)(C)C",
        -8.837166501626118, -9.01682459926845, None,
        REASON_GI + "; plus the group-7 data misprint",
    ),
    (
        "tetraol hydroperoxide",
        "C(C(C(CO)OO)O)O",
        -11.21071216429776, -11.607440300081954, None,
        REASON_OH + " (4 carbons: whole-molecule <C5 reassignment); also "
        + REASON_GI,
    ),
    # -- UMS data misprint only (group 7) --------------------------------
    (
        "hydroxy terpene alcohol",
        "C1(C(C(CCO)C1)(C)C)C(O)C",
        -5.917706528534178, -5.91770610472788, None,
        REASON_T7 + " (two group-7 occurrences)",
    ),
    (
        "ethylene glycol",
        "C(O)CO",
        -4.012097717735033, -4.012097411014363, None,
        REASON_T7 + " (two group-7 occurrences)",
    ),
    (
        "gem-diol on a menthol skeleton",
        "C(C(O)O)C1C(C(C(C)O)C1)(C)C",
        -7.767815979797785, -7.76210588301405, None,
        REASON_OH + " (gem-diol oxygens have one carbon neighbour)",
    ),
]


@pytest.mark.parametrize(
    "name,smi,ums_vp,ours_vp,status,reason",
    WL_26,
    ids=[e[0] for e in WL_26],
)
def test_example_output_whitelist(ums, name, smi, ums_vp, ours_vp, status,
                                  reason):
    """Documented deviations: the live oracle must reproduce the recorded
    UMS value (drift guard) and aerosolprops must reproduce the recorded
    paper-faithful side - value or refusal."""
    tb, live_ums = _ums_parts(ums, smi)
    assert live_ums == ums_vp, (name, "UMS drifted")
    w = _ours(smi, tb)
    live_ours = _ours_vp(w)
    if ours_vp is None:
        assert live_ours is None, (name, "refusal drifted")
        assert status in w.status, (name, w.status)
    else:
        assert live_ours == ours_vp, (name, "aerosolprops drifted")
        assert status is None
    assert reason  # every deviation carries a justification


def test_example_output_fully_partitioned():
    """No unexplained behaviour: each of the 26 Example_output SMILES is
    exact, last-ulp, or whitelisted - and only once."""
    smiles = ums_oracle.example_output_smiles()
    assert len(smiles) == 26
    assert len(EXACT_26) + len(ULP_26) + len(WL_26) == 26
    parts = [set(EXACT_26), set(ULP_26), {e[1] for e in WL_26}]
    assert sum(len(p) for p in parts) == 26
    assert not parts[0] & parts[1] and not parts[0] & parts[2]
    assert not parts[1] & parts[2]
    assert set(smiles) == parts[0] | parts[1] | parts[2]


# --- control battery outside Example_output ----------------------------
# Molecules whose UMS behaviour is identical live: bit-for-bit when the
# accumulation orders coincide, last-ulp otherwise (same lists as the
# Part 1 oracle controls, re-run through the Part 3 dB layer; adipic
# acid additionally pins a two-member same-class GI term and the cyclic
# ether/ester set pins cross-class zero interactions).
CONTROL_EXACT = [
    "O=C1OCCC1", "O=Cc1ccccc1", "CC(N)=O", "CNC(C)=O",
    "O=[N+]([O-])c1ccccc1", "OC(=O)c1ccccc1", "C1CCNC1", "C1COCCN1",
    "c1ccccc1Nc1ccccc1", "CNc1ccccc1", "CCC(=O)O", "OC(=O)CC(=O)O",
    "Cc1ccccc1",
]

CONTROL_ULP = [
    "CO", "Oc1ccccc1", "COC(=O)OC", "CC(=O)C", "C[N+](=O)[O-]",
    "Clc1ccccc1", "C1CCNCC1", "C1CCCCC1", "CCCCCCCCCC", "c1ccc(O)cc1",
    "c1ccc(Cl)cc1", "CC(=O)OO",
]


def test_control_battery_exact(ums):
    for smi in CONTROL_EXACT:
        tb, ums_vp = _ums_parts(ums, smi)
        assert _ours_vp(_ours(smi, tb)) == ums_vp, smi


def test_control_ulp_class(ums):
    for smi in CONTROL_ULP:
        tb, ums_vp = _ums_parts(ums, smi)
        diff = _ours_vp(_ours(smi, tb)) - ums_vp
        assert abs(diff) < 1e-9, (smi, diff)


# --- control battery: documented deviations ----------------------------
#: (name, SMILES, UMS log10(P/atm), aerosolprops log10(P/atm) or None when
#:  refused, refusal-status substring or None, reason)
WL_CONTROL = [
    (
        "biphenyl",
        "c1ccc(-c2ccccc2)cc1",
        -2.2829193816792883, -2.3762746946777424, None,
        "UMS lacks the aromatic bridge-carbon group 214 (Table 1 row 214,"
        " paper example biphenyl): it counts only 15x10 and drops both"
        " bridge carbons, aerosolprops counts {15: 10, 214: 2}",
    ),
    (
        "biacetyl",
        "CC(=O)C(C)=O",
        -1.4242427781778901, None,
        "refused: no dB parameter for group 118",
        "UMS lacks the 1,2-diketone group 118 (Table 1, PR 1, paper"
        " example 2,3-butanedione) and fragments it as two plain"
        " ketones; aerosolprops counts it and refuses - the paper prints no"
        " value for it",
    ),
    (
        "p-xylene",
        "Cc1ccc(C)cc1",
        -2.1924455054407987, None,
        "refused: no dB parameter for group 129", REASON_129,
    ),
    (
        "asparagine-like",
        "NC(=O)CC(N)C(=O)O",
        -10.87597707842977, None,
        "refused: missing interaction parameter (C, M)",
        "Table 6 row 218 (C, M) is printed 'Do not estimate'; UMS's"
        " interaction matrix zero-fills the pair and still returns a"
        " pressure, aerosolprops refuses",
    ),
    (
        "beta-alanine",
        "O=C(O)CCN",
        -3.849733825698069, None,
        "refused: missing interaction parameter (C, M)",
        "Table 6 row 218 (C, M) is printed 'Do not estimate'; UMS's"
        " interaction matrix zero-fills the pair and still returns a"
        " pressure, aerosolprops refuses",
    ),
    (
        "methyl ethyl sulfide",
        "CSCC",
        2.3279115095153444, 2.4433484446747062, None,
        "UMS's SMARTS catalogue misses the sulfur side (sulfide 54 and"
        " the chain CH2, Part 1 finding inherited by the Part 3 data);"
        " aerosolprops counts {1: 2, 4: 1, 54: 1}",
    ),
    (
        "heptane-2,4-dithiol",
        "CC(S)CC(S)CC",
        0.6867341050137132, 0.7540352353801216, None,
        "UMS's SMARTS catalogue misses the thiol family (Part 1 finding);"
        " aerosolprops counts {1: 2, 4: 2, 5: 2, 53: 2}",
    ),
    (
        "hexafluoroacetone",
        "FC(F)(F)C(=O)C(F)(F)F",
        2.06993716220592, 2.1536861964134095, None,
        "UMS misses the CF3/CF2 family 21 (Part 1 finding); this exact"
        " molecule is the paper's worked example 15c, whose group set"
        " the Part 3 gates verify end-to-end",
    ),
    (
        "isopropanol",
        "CC(C)O",
        -0.5234719130660291, -0.5232139712716641, None,
        REASON_OH + " (3 carbons: whole-molecule <C5 reassignment)",
    ),
    (
        "glycerol",
        "OCC(O)CO",
        -6.695824346639082, -6.69322532916196, None,
        REASON_OH + " (3 carbons: whole-molecule <C5 reassignment)",
    ),
]


@pytest.mark.parametrize(
    "name,smi,ums_vp,ours_vp,status,reason",
    WL_CONTROL,
    ids=[e[0] for e in WL_CONTROL],
)
def test_control_whitelist(ums, name, smi, ums_vp, ours_vp, status, reason):
    tb, live_ums = _ums_parts(ums, smi)
    assert live_ums == ums_vp, (name, "UMS drifted")
    w = _ours(smi, tb)
    live_ours = _ours_vp(w)
    if ours_vp is None:
        assert live_ours is None, (name, "refusal drifted")
        assert status in w.status, (name, w.status)
    else:
        assert live_ours == ours_vp, (name, "aerosolprops drifted")
        assert status is None
    assert reason


# --- control battery: the group-7/48 data misprints ---------------------
#: Identical group counts on both sides; the only divergence is UMS's
#: Part 3 data file misprinting group 7 (0.0378487 vs Table 4's
#: 37.8485e-3) and group 48 (0.2692470 vs Table 4's 269.2471e-3). The
#: recorded pairs lock both sides against drift; the delta must equal
#: the misprint scaled by the Eq. 6 slope (2e-7 x number of group-7/48
#: occurrences x slope), i.e. small but nonzero.
TYPO_7_48 = [
    # (SMILES, UMS log10(P/atm), aerosolprops log10(P/atm))
    ("CCO", -0.6445294154175775, -0.6445293881784675),
    ("CCCO", -1.2636720283151528, -1.2636719755195909),
    ("CCCCO", -1.8444251078148124, -1.8444250316255841),
    ("OCCO", -4.012097717735033, -4.012097411014363),
    ("CCOc1ccccc1", -2.518401106983425, -2.51840098957376),
    ("CCOC(C)=O", -0.6931975705686472, -0.6931975386633482),
    ("CC(=O)N(C)C", -2.470877802915269, -2.4708778591828606),
    ("CCO[N+](=O)[O-]", -1.12987627035553, -1.1298762210141744),
    ("ClCCl", -0.16819675693783726, -0.16819674873033552),
    ("ClC(Cl)Cl", -0.20580358583837635, -0.20580357554859852),
    ("CCNCC", -0.4809739342232625, -0.4809738886465107),
    ("CCCCCCCCO", -3.7313025629192937, -3.731302413815132),
    ("CCCCCCCCCCO", -4.7046845139710225, -4.704684329988518),
    ("CC(O)C(C)C", -1.4459821568251, -1.4459820966645562),
    ("CCCCOCCCCO", -4.064879736175726, -4.064879240452171),
    ("CCOCCOCC", -1.6721955240699753, -1.6721952169843002),
    ("CCOC(=O)CC", -1.210099867148592, -1.2100998121440627),
    ("CCCC(N)C", -1.2791526989910602, -1.2791526401907998),
    ("CC(C)N", -0.19356717481397498, -0.19356716568673615),
    ("NCCCN", -2.2451395278861668, -2.245139348047472),
    ("CC(C)COO", -2.3515726971781064, -2.3515726016589134),
    ("O=C(O)CCO", -5.403815387482973, -5.403815202668671),
    ("CC(N)C(=O)N(C)C", -3.956995632231873, -3.9569955476999783),
    ("OCC(O)C(O)C(O)C(O)CO", -15.823403971691643, -15.823401685791584),
]


@pytest.mark.parametrize("smi,ums_vp,ours_vp", TYPO_7_48)
def test_control_data_misprint_class(ums, smi, ums_vp, ours_vp):
    """Counts identical, values locked on both sides, and a small
    nonzero delta - the signature of the UMS group-7/48 data misprints
    (hexol carries six group-7 occurrences and shows the largest
    delta, ~2.3e-6)."""
    tb, live_ums = _ums_parts(ums, smi)
    assert live_ums == ums_vp, (smi, "UMS drifted")
    w = _ours(smi, tb)
    assert w.success, (smi, w.status)
    assert _ours_vp(w) == ours_vp, (smi, "aerosolprops drifted")
    assert w.counts == _ums_counts(ums, smi), smi
    assert ours_vp != ums_vp and abs(ours_vp - ums_vp) < 1e-5, smi


def test_control_extension_303_misprint(ums):
    """PAN: identical counts {1: 1, 303: 1}; the delta is exactly the
    UMS extension-303 misprint (0.519 vs Compernolle 2010 Table 3's
    519.03e-3) scaled by the Eq. 6 slope at Tb = 350 K (dTb chosen so
    the slope, and hence the visible delta, is sizeable)."""
    tb = 350.0
    compound = ums_oracle.readstring(ums, "CC(=O)OON(=O)=O")
    ums_vp = ums.vapour_pressures.nannoolal(compound, T, tb)
    assert ums_vp == -0.9084860452240575
    w = _ours("CC(=O)OON(=O)=O", tb)
    assert w.success and w.counts == _ums_counts(ums, "CC(=O)OON(=O)=O")
    assert _ours_vp(w) == -0.9084921596108497
    slope = (T / tb - 1.0) / (T / tb - 0.125)
    assert abs((_ours_vp(w) - ums_vp) - 3e-5 * slope) < 1e-12
