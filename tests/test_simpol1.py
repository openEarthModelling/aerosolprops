# SPDX-License-Identifier: MIT
"""Verification gates 1-4 for the SIMPOL.1 implementation.

Gate 1: coefficient completeness (31 groups x 4 coefficients, table-
    driven against an independent transcription of paper Table 5 via
    the method notes section 4) plus the two derivable columns of paper
    Table 6 (b_k at 293.15 K, and the d(dh_vap,k)/dT < 0 flag).
Gate 2: the paper's section 3.3 whole-molecule gold standards at
    293.15 K (cyclohexene, adipic acid) and the COOH transform.
Gate 3: the full group-accounting suite (paper Table 5 footnotes and
    section 3.3 examples; 33 cases verified in the notes plus 16 edge
    cases probed during implementation).
Gate 4: numerical consistency of the enthalpy estimator with the
    pressure estimator (Clausius-Clapeyron identity by finite
    differences; doctest-level b_k samples live in the modules).
"""

import math

import pytest

from aerosolprops.methods.simpol1 import (
    LOG10_ATM_TO_PA,
    SIMPOL1_GAS_CONSTANT,
    SIMPOL1_GROUPS,
    SIMPOL1_SMARTS,
    SIMPOL1_ALGORITHMIC_GROUPS,
    simpol1_bk,
    simpol1_log10_vapour_pressure,
    simpol1_vapourisation_enthalpy,
)
from aerosolprops.structure import mol_from_smiles
from aerosolprops.vapour_pressure import SIMPOL1, simpol1_fragmentation


# ---------------------------------------------------------------------------
# Gate 1: coefficient catalogue
# ---------------------------------------------------------------------------

#: Independent transcription of paper Table 5 (B_k1, B_k2, B_k3, B_k4),
#: taken from the method notes section 4 table (itself transcribed from
#: the paper text layer, p. 2783). The source module holds its own
#: transcription; equality of the two is the check.
PAPER_TABLE5 = {
    0: (-4.26938e02, 2.89223e-01, 4.42057e-03, 2.92846e-01),
    1: (-4.11248e02, 8.96919e-01, -2.48607e-03, 1.40312e-01),
    2: (-1.46442e02, 1.54528e00, 1.71021e-03, -2.78291e-01),
    3: (3.50262e01, -9.20839e-01, 2.24399e-03, -9.36300e-02),
    4: (-8.72770e01, 1.78059e00, -3.07187e-03, -1.04341e-01),
    5: (5.73335e00, 1.69764e-02, -6.28957e-04, 7.55434e-03),
    6: (-2.61268e02, -7.63282e-01, -1.68213e-03, 2.89038e-01),
    7: (-7.25373e02, 8.26326e-01, 2.50957e-03, -2.32304e-01),
    8: (-7.29501e02, 9.86017e-01, -2.92664e-03, 1.78077e-01),
    9: (-1.37456e01, 5.23486e-01, 5.50298e-04, -2.76950e-01),
    10: (-7.98796e02, -1.09436e00, 5.24132e-03, -2.28040e-01),
    11: (-3.93345e02, -9.51778e-01, -2.19071e-03, 3.05843e-01),
    12: (-1.44334e02, -1.85617e00, -2.37491e-05, 2.88290e-01),
    13: (4.05265e01, -2.43780e00, 3.60133e-03, 9.86422e-02),
    14: (-7.07406e01, -1.06674e00, 3.73104e-03, -1.44003e-01),
    15: (-7.83648e02, -1.03439e00, -1.07148e-03, 3.15535e-01),
    16: (-5.63872e02, -7.18416e-01, 2.63016e-03, -4.99470e-02),
    17: (-4.53961e02, -3.26105e-01, -1.39780e-04, -3.93916e-02),
    18: (3.71375e01, -2.66753e00, 1.01483e-03, 2.14233e-01),
    19: (-5.03710e02, 1.04092e00, -4.12746e-03, 1.82790e-01),
    20: (-3.59763e01, -4.08458e-01, 1.67264e-03, -9.98919e-02),
    21: (-6.09432e02, 1.50436e00, -9.09024e-04, -1.35495e-01),
    22: (-1.02367e02, -7.16253e-01, -2.90670e-04, -5.88556e-01),
    23: (-1.93802e03, 6.48262e-01, 1.73245e-03, 3.47940e-02),
    24: (-5.26919e00, 3.06435e-01, 3.25397e-03, -6.81506e-01),
    25: (-2.84042e02, -6.25424e-01, -8.22474e-04, -8.80549e-02),
    26: (1.50093e02, 2.39875e-02, -3.37969e-03, 1.52789e-02),
    27: (-2.03387e01, -5.48718e00, 8.39075e-03, 1.07884e-01),
    28: (-8.38064e02, -1.09600e00, -4.24385e-04, 2.81812e-01),
    29: (-5.27934e01, -4.63689e-01, -5.11647e-03, 3.84965e-01),
    30: (-1.61520e03, 9.01669e-01, 1.44536e-03, 2.66889e-01),
}

#: Paper Table 6 (p. 2784): b_k(293.15 K) to three significant figures.
PAPER_TABLE6_BK = {
    0: 1.79, 1: -0.438, 2: -0.0338, 3: -0.675, 4: -0.0104, 5: -0.105,
    6: -0.506, 7: -2.23, 8: -1.35, 9: -0.935, 10: -3.58, 11: -1.20,
    12: -0.718, 13: -0.683, 14: -1.03, 15: -2.23, 16: -2.15, 17: -2.14,
    18: -1.03, 19: -0.849, 20: -0.608, 21: -1.61, 22: -4.49, 23: -5.26,
    24: -2.63, 25: -2.34, 26: -0.368, 27: -2.48, 28: -2.48, 29: 0.0432,
    30: -2.67,
}

#: Paper Table 6, "d(dh_vap,k)/dT < 0 at T = 293.15 K?" column. By paper
#: Eq. 11 the per-group derivative is 2.303*R*(2*B_k3*T + B_k4), so this
#: column is fully derivable from the Table 5 coefficients - a strong
#: independent check on the B_k3/B_k4 columns.
PAPER_TABLE6_DHDH_NEGATIVE = {
    0: False, 1: True, 2: False, 3: False, 4: True, 5: True, 6: True,
    7: False, 8: True, 9: False, 10: False, 11: True, 12: False,
    13: False, 14: False, 15: True, 16: False, 17: True, 18: False,
    19: True, 20: False, 21: True, 22: True, 23: False, 24: False,
    25: True, 26: True, 27: False, 28: False, 29: True, 30: False,
}


def test_gate1_catalogue_structure():
    assert sorted(SIMPOL1_GROUPS) == list(range(31))
    for record in SIMPOL1_GROUPS.values():
        assert len(record) == 5
        assert isinstance(record[0], str)
        assert all(isinstance(v, float) for v in record[1:])


@pytest.mark.parametrize("k", list(range(31)))
def test_gate1_coefficients_match_paper_table5(k):
    assert SIMPOL1_GROUPS[k][1:] == PAPER_TABLE5[k]


@pytest.mark.parametrize("k", list(range(31)))
def test_gate1_bk_matches_paper_table6(k):
    # Table 6 gives three significant figures; compare at that precision.
    assert "{:.3g}".format(simpol1_bk(k, 293.15)) == "{:.3g}".format(
        PAPER_TABLE6_BK[k]
    )


@pytest.mark.parametrize("k", list(range(31)))
def test_gate1_enthalpy_slope_flag_matches_paper_table6(k):
    _, _b1, _b2, b3, b4 = SIMPOL1_GROUPS[k]
    slope_negative = (2.0 * b3 * 293.15 + b4) < 0.0
    assert slope_negative == PAPER_TABLE6_DHDH_NEGATIVE[k]


def test_gate1_smarts_catalogue_markers():
    # Groups 1-30 have catalogue entries: 27 SMARTS patterns plus the 3
    # algorithmic placeholders (2, 3, 4); group 0 is the fitted constant.
    assert sorted(SIMPOL1_SMARTS) == list(range(1, 31))
    algorithmic = sorted(k for k, v in SIMPOL1_SMARTS.items() if v is None)
    assert algorithmic == sorted(SIMPOL1_ALGORITHMIC_GROUPS) == [2, 3, 4]
    assert all(isinstance(v, str) for k, v in SIMPOL1_SMARTS.items() if v)


# ---------------------------------------------------------------------------
# Gate 2: paper section 3.3 whole-molecule gold standards (T = 293.15 K)
# ---------------------------------------------------------------------------

def test_gate2_cyclohexene_gold_standard():
    # Paper section 3.3: calculated -0.94 (log10 atm), measured -1.08;
    # method-notes re-computation: -0.9490 (Pa basis +4.0567).
    log10p = simpol1_log10_vapour_pressure({1: 6, 4: 1, 5: 1}, 293.15)
    assert log10p == pytest.approx(4.0567, abs=5e-4)
    atm = log10p - LOG10_ATM_TO_PA
    assert atm == pytest.approx(-0.9490, abs=5e-4)
    assert abs(atm - (-0.94)) < 0.01  # paper's displayed precision
    # Measured -1.08: within the paper's overall fit sigma of 0.34.
    assert abs(atm - (-1.08)) < 0.34


def test_gate2_adipic_acid_gold_standard():
    # Paper section 3.3: calculated -7.99 (log10 atm), measured -8.49;
    # method-notes re-computation: -7.9900 (Pa basis -2.9843).
    log10p = simpol1_log10_vapour_pressure({1: 6, 10: 2}, 293.15)
    assert log10p == pytest.approx(-2.9843, abs=5e-4)
    atm = log10p - LOG10_ATM_TO_PA
    assert atm == pytest.approx(-7.9900, abs=5e-4)
    assert abs(atm - (-7.99)) < 0.01  # paper's displayed precision
    # Measured -8.49: 0.50 above the fit; the paper itself flags growing
    # errors at the low-volatility end (its Figs. 10-12), so 0.6 bounds
    # this documented example rather than the global sigma of 0.34.
    assert abs(atm - (-8.49)) < 0.6


def test_gate2_cooh_transform():
    # Paper section 3.3: adding 2 COOH to cyclohexene changes log10 p by
    # 2*b10 = -7.16 at 293.15 K.
    assert 2.0 * simpol1_bk(10, 293.15) == pytest.approx(-7.16, abs=0.005)


# ---------------------------------------------------------------------------
# Gate 3: group accounting (nu vectors)
# ---------------------------------------------------------------------------

#: (name, SMILES, expected non-zero counts without the zeroth group).
#: Cases 1-33 are the paper Table 5 footnote / section 3.3 examples
#: verified in the method notes section 5.3; the remainder are edge
#: cases exercised during implementation (small amides, lactams, ring
#: heterocycles, bridged rings).
ACCOUNTING_CASES = [
    ("cyclohexene", "C1=CCCCC1", {1: 6, 4: 1, 5: 1}),
    ("adipic acid", "OC(=O)CCCCC(=O)O", {1: 6, 10: 2}),
    ("methyl acetate", "COC(C)=O", {1: 3, 11: 1}),
    ("methyl nitroacetate", "COC(=O)C[N+](=O)[O-]", {1: 3, 16: 1, 30: 1}),
    ("methyl 2-nitropropionate", "COC(=O)C(C)[N+](=O)[O-]", {1: 4, 16: 1, 30: 1}),
    ("THF", "C1CCOC1", {1: 4, 4: 1, 13: 1}),
    ("anisole", "COc1ccccc1", {1: 7, 3: 1, 14: 1}),
    ("acetamide", "CC(N)=O", {1: 2, 2: 2, 22: 1}),
    ("n-propylbutamide", "CCCC(NCCC)=O", {1: 7, 2: 4, 23: 1}),
    ("o-cresol", "Cc1ccccc1O", {1: 7, 3: 1, 17: 1}),
    ("m-cresol", "Cc1cccc(O)c1", {1: 7, 3: 1, 17: 1}),
    ("2-nitrophenol", "Oc1ccccc1[N+](=O)[O-]", {1: 6, 3: 1, 16: 1, 29: 1}),
    ("3-nitrophenol", "Oc1cccc([N+](=O)[O-])c1", {1: 6, 3: 1, 16: 1, 29: 1}),
    ("4-nitrophenol", "Oc1ccc([N+](=O)[O-])cc1", {1: 6, 3: 1, 16: 1, 29: 1}),
    ("methylamine", "CN", {1: 1, 18: 1}),
    ("dimethylamine", "CNC", {1: 2, 19: 1}),
    ("trimethylamine", "CN(C)C", {1: 3, 20: 1}),
    ("aniline", "Nc1ccccc1", {1: 6, 3: 1, 21: 1}),
    ("DMF", "CN(C)C=O", {1: 3, 2: 1, 24: 1}),
    ("cyclohex-2-enone", "O=C1C=CCCC1", {1: 6, 4: 1, 5: 1, 6: 1, 9: 1}),
    ("peroxyacetyl nitrate (PAN)", "CC(=O)OON(=O)=O", {1: 2, 25: 1}),
    ("tert-butyl hydroperoxide", "CC(C)(C)OO", {1: 4, 27: 1}),
    ("peracetic acid", "CC(=O)OO", {1: 2, 28: 1}),
    ("di-tert-butyl peroxide", "CC(C)(C)OOC(C)(C)C", {1: 8, 26: 1}),
    ("ethyl nitrate", "CCO[N+](=O)[O-]", {1: 2, 15: 1}),
    ("benzaldehyde", "O=Cc1ccccc1", {1: 7, 3: 1, 8: 1}),
    ("methyl benzoate", "COC(=O)c1ccccc1", {1: 8, 3: 1, 11: 1}),
    ("2-nitropropane", "CC(C)[N+](=O)[O-]", {1: 3, 16: 1}),
    ("nitrobenzene", "O=[N+]([O-])c1ccccc1", {1: 6, 3: 1, 16: 1}),
    ("benzoic acid", "OC(=O)c1ccccc1", {1: 7, 3: 1, 10: 1}),
    ("diphenyl ether", "c1ccc(Oc2ccccc2)cc1", {1: 12, 3: 2, 14: 1}),
    ("ethylenediamine", "NCCN", {1: 2, 18: 2}),
    ("levoglucosan", "C1C2C(C(C(C(O1)O2)O)O)O", {1: 6, 4: 2, 7: 3, 13: 2}),
    ("naphthalene", "c1ccc2ccccc2c1", {1: 10, 3: 2}),
    ("anthracene", "c1ccc2cc3ccccc3cc2c1", {1: 14, 3: 3}),
    ("biphenyl", "c1ccc(-c2ccccc2)cc1", {1: 12, 3: 2}),
    ("indane", "c1ccc2c(c1)CCC2", {1: 9, 3: 1, 4: 1}),
    ("decalin", "C1CCC2CCCCC2C1", {1: 10, 4: 2}),
    ("phenol", "Oc1ccccc1", {1: 6, 3: 1, 17: 1}),
    # --- implementation edge cases (documented conservative choices) ---
    # Urea: two primary amide matches on one carbonyl; both acid-side
    # cuts return the whole (still connected) molecule.
    ("urea", "NC(N)=O", {1: 1, 2: 2, 22: 1}),
    # Lactams (OQ-1, undefined in the paper): the C-N graph cut leaves
    # the ring connected, so nu_2 counts every carbon of the component
    # containing the acyl carbon - the whole ring here (plus the
    # N-substituent's carbons, which remain attached to the ring).
    ("piperidin-2-one", "O=C1CCCCN1", {1: 5, 2: 5, 4: 1, 23: 1}),
    ("1-methyl-piperidin-2-one", "CN1CCCCC1=O", {1: 6, 2: 6, 4: 1, 24: 1}),
    ("formamide", "NC=O", {1: 1, 2: 1, 22: 1}),
    ("N-methylformamide", "CNC=O", {1: 2, 2: 1, 23: 1}),
    ("N-methylacetamide", "CNC(C)=O", {1: 3, 2: 2, 23: 1}),
    # Pyridine N is aromatic (not an "amine, aromatic" instance).
    ("pyridine", "c1ccncc1", {1: 5, 3: 1}),
    # Furan is an all-aromatic SSSR ring (group 3).
    ("furan", "c1ccoc1", {1: 4, 3: 1}),
    ("vinyl alcohol", "C=CO", {1: 2, 5: 1, 7: 1}),
    # Bridged ring: minimal (not symmetrized) SSSR -> cyclomatic count.
    ("fenchyl alcohol", "CC1(C2CCC(O)C1C2)C", {1: 9, 4: 2, 7: 1}),
    # No SIMPOL.1 group beyond the carbon count matches hydroxylamine;
    # the method simply has no parameterisation for it.
    ("hydroxylamine", "NO", {}),
]


@pytest.mark.parametrize(
    "name, smi, expected",
    ACCOUNTING_CASES,
    ids=[case[0] for case in ACCOUNTING_CASES],
)
def test_gate3_group_accounting(name, smi, expected):
    counts, success, status = simpol1_fragmentation(mol_from_smiles(smi))
    assert success, status
    assert counts == {0: 1, **expected}


# ---------------------------------------------------------------------------
# Gate 4: estimator consistency
# ---------------------------------------------------------------------------

CYCLOHEXENE = {1: 6, 4: 1, 5: 1}
ADIPIC = {1: 6, 10: 2}


@pytest.mark.parametrize("counts", [CYCLOHEXENE, ADIPIC])
@pytest.mark.parametrize("T", [273.15, 298.15, 333.15, 393.15])
def test_gate4_enthalpy_clausius_clapeyron_identity(counts, T):
    # Paper Eq. 9: d(log10 p)/d(1/T) = -dHvap/(2.303 R), verified by
    # central finite differences against the pressure estimator.
    dT = 1e-4
    hi = simpol1_log10_vapour_pressure(counts, T + dT)
    lo = simpol1_log10_vapour_pressure(counts, T - dT)
    slope = (hi - lo) / (1.0 / (T + dT) - 1.0 / (T - dT))
    hvap = simpol1_vapourisation_enthalpy(counts, T)
    assert hvap == pytest.approx(-2.303 * SIMPOL1_GAS_CONSTANT * slope, rel=1e-6)


def test_gate4_enthalpy_positive_and_decreasing_cyclohexene():
    # Paper: "for any real compound i in the liquid state, dHvap > 0
    # and d(dHvap)/dT < 0". All groups of cyclohexene (0, 1, 4, 5) are
    # flagged YES in its Table 6, so both properties hold.
    prev = None
    for T in (273.15, 298.15, 333.15, 353.15, 393.15):
        hvap = simpol1_vapourisation_enthalpy(CYCLOHEXENE, T)
        assert hvap > 0.0
        if prev is not None:
            assert hvap < prev
        prev = hvap


def test_gate4_enthalpy_adipic_increases_as_table6_flags():
    # Groups 0 (zeroth) and 10 (carboxylic acid) are flagged "d(dh/dT)
    # < 0? NO" in paper Table 6, and their positive slopes dominate
    # adipic acid (2 COOH): its predicted dHvap rises with T. This is a
    # property of the paper's own published coefficients, not a bug;
    # asserted here so any change in sign handling is noticed.
    prev = None
    for T in (273.15, 298.15, 333.15, 353.15, 393.15):
        hvap = simpol1_vapourisation_enthalpy(ADIPIC, T)
        assert hvap > 0.0
        if prev is not None:
            assert hvap > prev
        prev = hvap


def test_gate4_estimator_input_validation():
    assert simpol1_bk(31, 298.15) is None
    assert simpol1_bk(0, 0.0) is None
    assert simpol1_log10_vapour_pressure({99: 1}) is None
    assert simpol1_vapourisation_enthalpy({99: 1}) is None
    # A count under key 0 is ignored (nu_0 is always exactly 1).
    base = simpol1_log10_vapour_pressure(CYCLOHEXENE, 293.15)
    padded = simpol1_log10_vapour_pressure({0: 7, **CYCLOHEXENE}, 293.15)
    assert padded == pytest.approx(base)


# ---------------------------------------------------------------------------
# Wrapper class
# ---------------------------------------------------------------------------

def test_wrapper_smiles_and_mol_construction_agree():
    from_mol = SIMPOL1(mol=mol_from_smiles("Oc1ccccc1"))
    from_smi = SIMPOL1(smi="Oc1ccccc1")
    assert from_mol.success and from_smi.success
    assert from_mol.counts == from_smi.counts == {0: 1, 1: 6, 3: 1, 17: 1}


def test_wrapper_normalizes_valence_forms():
    # The neutral valence form of nitrobenzene must give the same counts
    # as the charge-separated form the SMARTS catalogue targets.
    assert SIMPOL1(smi="c1ccc([N](=O)=O)cc1").counts == SIMPOL1(
        smi="O=[N+]([O-])c1ccccc1"
    ).counts


def test_wrapper_invalid_smiles_fails_cleanly():
    sim = SIMPOL1(smi="not a smiles")
    assert not sim.success
    assert sim.counts == {}
    assert sim.status == "Failed to construct mol"
    assert sim.log10_vapour_pressure() is None
    assert sim.vapour_pressure() is None
    assert sim.vapourisation_enthalpy() is None


def test_wrapper_no_input_fails_cleanly():
    sim = SIMPOL1()
    assert not sim.success
    assert sim.vapour_pressure() is None


def test_wrapper_unsupported_element():
    # SIMPOL.1 parameterises C, H, N, O only.
    sim = SIMPOL1(smi="CCCl")
    assert not sim.success
    assert "unsupported element" in sim.status
    assert sim.log10_vapour_pressure() is None


def test_wrapper_matches_pure_estimator():
    sim = SIMPOL1(smi="C1=CCCCC1")
    assert sim.log10_vapour_pressure(293.15) == simpol1_log10_vapour_pressure(
        sim.counts, 293.15
    )
    assert sim.vapour_pressure(293.15) == pytest.approx(
        10.0 ** sim.log10_vapour_pressure(293.15)
    )
    assert sim.vapourisation_enthalpy(300.0) == simpol1_vapourisation_enthalpy(
        sim.counts, 300.0
    )
    assert math.isfinite(sim.vapour_pressure(293.15))
