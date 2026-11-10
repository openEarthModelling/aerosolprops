# SPDX-License-Identifier: MIT
"""Verification gates 1-3 for the Myrdal-Yalkowsky pure layer.

Gate 1: coefficient triple guard - the module's constants against an
    independent transcription of the paper text (ie950242l text layer,
    Eqs. 2/7-12), plus the paper-rounding identity 19.1 ~ 2.303 R.
Gate 2: paper Table 1, all 19 rows, through the pure estimator with the
    descriptor conventions locked in the wrapper gate (SP3/SP2/RING
    walk, sqrt-structure HBN, RDKit average MWs); every row reproduces
    the paper's predicted column to its printed 2 decimals.
Gate 3: the same 19 rows against the *observed* column (the paper's
    compiled experimental values, Boublik et al. 1973 and refs
    therein): the average absolute error and its maximum reproduce the
    paper's own summary statistics (0.09 log units / factor 1.23, all
    but two within 0.15).
"""

import math

import pytest

from aerosolprops.methods.myrdal_yalkowsky import (
    LOG10_ATM_TO_PA,
    MY_DCPGL_INTERCEPT,
    MY_DCPGL_TAU,
    MY_DENOMINATOR_CONSTANT,
    MY_NH2_ROOT_COEFFICIENT,
    MY_RING_WEIGHT,
    MY_SB_HBN,
    MY_SB_INTERCEPT,
    MY_SB_TAU,
    MY_SM_INTERCEPT,
    MY_SM_LOG_SIGMA,
    MY_SM_TAU,
    MY_SP2_WEIGHT,
    MY_TAU_OFFSET,
    myrdal_yalkowsky_entropy_of_boiling,
    myrdal_yalkowsky_entropy_of_melting,
    myrdal_yalkowsky_hbn,
    myrdal_yalkowsky_heat_capacity_change_boiling,
    myrdal_yalkowsky_log10_vapour_pressure,
    myrdal_yalkowsky_tau,
)

# ---------------------------------------------------------------------------
# Gate 1: coefficient catalogue (independent transcription, paper text)
# ---------------------------------------------------------------------------

#: Independent transcription of the in-line coefficients of the paper
#: text layer (pages 2495-2497). The source module holds its own
#: transcription; equality of the two is the check.
PAPER_COEFFICIENTS = {
    "sm_intercept": 56.5,        # Eq. 2 and Eq. 10
    "sm_log_sigma": -19.2,       # Eq. 2 and Eq. 10
    "sm_tau": 9.2,               # Eq. 10 (Eq. 2 has 9.2*(N-3))
    "sb_intercept": 86.0,        # Eq. 7 ("86"; Eq. 12 prints 86.0)
    "sb_tau": 0.4,               # Eq. 7
    "sb_hbn": 1421.0,            # Eq. 7
    "dcpgl_intercept": -90.0,    # Eq. 11
    "dcpgl_tau": -2.1,           # Eq. 11
    "denominator": 19.1,         # every denominator of Eq. 12
    "nh2_root": 0.33,            # Eq. 9
    "sp2_weight": 0.5,           # Eq. 8
    "ring_weight": 0.5,          # Eq. 8
    "tau_offset": -1.0,          # Eq. 8
}


def test_coefficients_triple_guard():
    module_view = {
        "sm_intercept": MY_SM_INTERCEPT,
        "sm_log_sigma": MY_SM_LOG_SIGMA,
        "sm_tau": MY_SM_TAU,
        "sb_intercept": MY_SB_INTERCEPT,
        "sb_tau": MY_SB_TAU,
        "sb_hbn": MY_SB_HBN,
        "dcpgl_intercept": MY_DCPGL_INTERCEPT,
        "dcpgl_tau": MY_DCPGL_TAU,
        "denominator": MY_DENOMINATOR_CONSTANT,
        "nh2_root": MY_NH2_ROOT_COEFFICIENT,
        "sp2_weight": MY_SP2_WEIGHT,
        "ring_weight": MY_RING_WEIGHT,
        "tau_offset": MY_TAU_OFFSET,
    }
    assert module_view == PAPER_COEFFICIENTS


def test_denominator_is_the_papers_rounding_of_2303r():
    """19.1 is the authors' rounding of 2.303 R (not 19.147): the
    identity must hold at the paper's precision while the constant
    itself stays what reproduces Table 1."""
    assert abs(MY_DENOMINATOR_CONSTANT - 2.303 * 8.314) < 0.06
    assert MY_DENOMINATOR_CONSTANT == 19.1


def test_appendix_a_intermediates():
    """Paper Appendix A (3,5-diethylphenol, tau=1.5, HBN=0.006 66):
    the tabulated intermediate quantities 96.06 and -93.15 J/(mol K)."""
    assert round(myrdal_yalkowsky_entropy_of_boiling(1.5, 0.00666), 2) == 96.06
    assert (
        round(myrdal_yalkowsky_heat_capacity_change_boiling(1.5), 2) == -93.15
    )


def test_invalid_inputs_return_none():
    for bad in (
        lambda: myrdal_yalkowsky_tau(-1, 0, 0),
        lambda: myrdal_yalkowsky_hbn(-1, 0, 100.0),
        lambda: myrdal_yalkowsky_hbn(1, 0, 0.0),
        lambda: myrdal_yalkowsky_entropy_of_boiling(-1.0, 0.0),
        lambda: myrdal_yalkowsky_entropy_of_melting(0.0, 0.5),
        lambda: myrdal_yalkowsky_heat_capacity_change_boiling(-1.0),
        lambda: myrdal_yalkowsky_log10_vapour_pressure(-5.0, 400.0),
        lambda: myrdal_yalkowsky_log10_vapour_pressure(300.0, 0.0),
        lambda: myrdal_yalkowsky_log10_vapour_pressure(
            300.0, 400.0, Tm=-1.0),
        lambda: myrdal_yalkowsky_log10_vapour_pressure(
            300.0, 400.0, sigma=0.0),
        lambda: myrdal_yalkowsky_log10_vapour_pressure(
            300.0, 400.0, tau=-1.0),
        lambda: myrdal_yalkowsky_log10_vapour_pressure(
            300.0, 400.0, hbn=-1e-6),
    ):
        assert bad() is None


# ---------------------------------------------------------------------------
# Gate 2: paper Table 1 - all 19 rows through the pure estimator
# ---------------------------------------------------------------------------

#: (name, Tb [K], T [K], Tm [K] or None, sigma, sp3, sp2, ring,
#:  n_hbond_root, n_primary_amine, MW [g/mol], paper predicted, paper
#:  observed) - every row of Table 1. Descriptors follow the wrapper
#: conventions (nonterminal = heavy-atom degree >= 2; ring atoms
#: excluded from SP3/SP2; tert-butyl centres excluded; RING = ring
#: systems; first HBN root = alcohols+phenols+acids, Compernolle (2010)
#: hydroperoxides/peracids included); the observed column is the
#: experimental literature the paper compiled.
TABLE1 = [
    ("1-undecene", 466, 379, None, 1, 8, 1, 0, 0, 0, 154.297, -1.20, -1.20),
    ("2-methyl-5-ethylpyridine", 452, 325, None, 1, 1, 0, 1, 0, 0,
     121.183, -2.05, -2.10),
    ("3,5-diethylphenol", 521, 388, None, 1, 2, 0, 1, 1, 0, 150.221,
     -1.96, -1.88),
    ("cis-decalin", 469, 373, None, 1, 0, 0, 1, 0, 0, 138.254, -1.29, -1.26),
    ("beta-naphthol", 556, 415, None, 1, 0, 0, 1, 1, 0, 144.173, -1.93,
     -1.88),
    ("isopropylcyclohexane", 428, 344, None, 1, 1, 0, 1, 0, 0, 126.243,
     -1.22, -1.20),
    ("methyl salicylate", 494, 352, None, 1, 1, 1, 1, 1, 0, 152.149,
     -2.33, -2.25),
    ("n-ethylaniline", 478, 323, None, 1, 2, 0, 1, 0, 0, 121.183, -2.60,
     -2.50),
    ("resorcinol", 549, 425, None, 1, 0, 0, 1, 2, 0, 110.112, -1.76, -1.88),
    ("valeric acid", 446, 345, None, 1, 3, 1, 0, 1, 0, 102.133, -1.73,
     -1.86),
    ("1-nonanol", 486, 365, None, 1, 8, 0, 0, 1, 0, 144.258, -1.96, -2.13),
    ("1,4-dioxane", 374, 303, None, 1, 0, 0, 1, 0, 0, 88.106, -1.17, -1.20),
    ("1-hexanethiol", 426, 354, None, 1, 5, 0, 0, 0, 0, 118.245, -1.03,
     -1.02),
    ("1-tetradecanol", 569, 425, None, 1, 13, 0, 0, 1, 0, 214.393, -2.01,
     -2.16),
    ("1-chlorotetradecane", 570, 415, None, 1, 13, 0, 0, 0, 0, 232.839,
     -2.11, -2.16),
    ("biphenyl", 528, 342, 342.4, 2, 0, 0, 2, 0, 0, 154.212, -2.97, -2.99),
    ("beta-pinene", 439, 300, None, 1, 0, 0, 1, 0, 0, 136.238, -2.48,
     -2.36),
    ("2,4,4-trimethylhexane", 404, 323, None, 1, 4, 0, 0, 0, 0, 128.259,
     -1.28, -1.20),
    ("octacosane", 705, 408, None, 1, 26, 0, 0, 0, 0, 394.772, -5.01, -5.40),
]


def _table1_prediction(row):
    (_, Tb, T, Tm, sigma, sp3, sp2, ring, n_hb, n_nh2, mw, _, _) = row
    tau = myrdal_yalkowsky_tau(sp3, sp2, ring)
    hbn = myrdal_yalkowsky_hbn(n_hb, n_nh2, mw)
    lp = myrdal_yalkowsky_log10_vapour_pressure(
        float(T), float(Tb), tau, hbn, Tm=Tm, sigma=float(sigma)
    )
    assert lp is not None
    return lp - LOG10_ATM_TO_PA


@pytest.mark.parametrize(
    "row", TABLE1, ids=[r[0] for r in TABLE1]
)
def test_table1_predictions(row):
    """Every Table 1 row reproduces the paper's predicted value at its
    printed two decimals."""
    assert round(_table1_prediction(row), 2) == row[11]


# ---------------------------------------------------------------------------
# Gate 3: the observed column - the paper's own summary statistics
# ---------------------------------------------------------------------------

def test_table1_literature_accuracy():
    """Against the observed (experimental literature) column the 19
    predictions reproduce the paper's summary: average absolute error
    0.09 log units (factor error 1.23), every row except octacosane
    within 0.15, maximum 0.39."""
    errors = [abs(_table1_prediction(r) - r[12]) for r in TABLE1]
    assert len(errors) == 19
    assert abs(sum(errors) / len(errors) - 0.09) < 0.005
    assert max(errors) == pytest.approx(0.39, abs=0.005)
    # "except two ... within 0.15 log units": 1-nonanol (0.16) and
    # octacosane (0.39, the outlier the paper discusses).
    assert sum(e <= 0.1501 for e in errors) == 17
    # The paper also reports eq 12 beating the Yalkowsky-Mishra eq 6 on
    # 17 of the 19; eq 6 is not implemented (lineage only).


def test_hbn_root_structure_is_load_bearing():
    """The sqrt structure of Eq. 9 is what Table 1's resorcinol row
    locks: a linear count of 2 hydroxyls would give the *observed*
    -1.88, the sqrt(2) count gives the *predicted* -1.76."""
    row = TABLE1[8]
    assert row[0] == "resorcinol"
    (_, Tb, T, Tm, sigma, sp3, sp2, ring, n_hb, n_nh2, mw, _, _) = row
    linear = (n_hb + MY_NH2_ROOT_COEFFICIENT * n_nh2) / mw
    lp = myrdal_yalkowsky_log10_vapour_pressure(
        float(T), float(Tb), 0.0, linear)
    assert round(lp - LOG10_ATM_TO_PA, 2) == -1.88  # wrong model -> obs
    assert round(_table1_prediction(row), 2) == -1.76  # sqrt -> pred


def test_solid_vs_liquid_branch():
    """Below Tm the melting term is present; dropping Tm removes it
    (the sub-cooled liquid, paper Appendix A ruling)."""
    solid = myrdal_yalkowsky_log10_vapour_pressure(
        300.0, 500.0, 2.0, 0.0, Tm=350.0, sigma=3.0)
    liquid = myrdal_yalkowsky_log10_vapour_pressure(
        300.0, 500.0, 2.0, 0.0)
    melting_term = myrdal_yalkowsky_entropy_of_melting(2.0, 3.0) * (
        350.0 - 300.0) / (MY_DENOMINATOR_CONSTANT * 300.0)
    assert solid is not None and liquid is not None
    assert solid == pytest.approx(liquid - melting_term)
    # exactly at the melting point the compound is liquid (Appendix A)
    at_tm = myrdal_yalkowsky_log10_vapour_pressure(
        350.0, 500.0, 2.0, 0.0, Tm=350.0, sigma=3.0)
    assert at_tm == pytest.approx(
        myrdal_yalkowsky_log10_vapour_pressure(350.0, 500.0, 2.0, 0.0))
