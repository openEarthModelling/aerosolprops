# SPDX-License-Identifier: MIT
"""Verification gate 5 for Myrdal-Yalkowsky: live UManSysProp comparison.

Runs the *unmodified* UManSysProp ``vapour_pressures.myrdal_and_yalkowsky``
live as the behavioural oracle (UMS's SMARTS and code are never copied
into aeroprops) over all 19 paper Table 1 SMILES plus a control battery,
feeding both implementations the same experimental Tb (the method's
input) and comparing log10(P/atm).

Result: every Table 1 row and 14 control molecules agree with UMS
bit-for-bit or at the last-ulp level (float summation order). UMS has
no melting-point input, so the comparison is always on the sub-cooled
liquid branch - the paper's own convention for its Table 1 rows; the
biphenyl row is compared at its Table 1 liquid reading.

The whitelist carries the live-recorded UMS and aeroprops values (both
sides drift-locked) and the paper evidence for every divergence:

* fused aromatics: UMS feeds the SSSR ring count into tau (anthracene/
  phenanthrene tau = 0.5); the paper's RING is *independent ring
  systems* (one fused system, tau = 0);
* tert-butyl centres: UMS counts them (tert-butyl hydroperoxide
  tau = 1); the paper states tert-butyl groups do not contribute to
  flexibility (tau = 0);
* sulfur: UMS omits S atoms from the SP3 walk (diethyl sulfide
  tau = 1); the paper's Eq. 8 explicitly includes S;
* enols and organic nitrates: UMS's wide OH/sp2 tables accept them and
  return a number; aeroprops refuses both (outside the published
  parameterisation, Compernolle 2010 Section 2.6).

Example_output.txt has no Myrdal-Yalkowsky column, so no exact-26
partition is possible; the Table 1 battery plays that role. Skipped
entirely when the UManSysProp source tree is not present
(``AEROPROPS_UMS_REPO``, see :mod:`ums_oracle`).
"""

import pytest

import ums_oracle
from aeroprops.methods.simpol1 import LOG10_ATM_TO_PA
from aeroprops.vapour_pressure import MyrdalYalkowsky

pytestmark = pytest.mark.skipif(
    not ums_oracle.available(), reason="UManSysProp source tree not available"
)


def _ums(ums, smi, T, tb):
    return ums.vapour_pressures.myrdal_and_yalkowsky(
        ums_oracle.readstring(ums, smi), T, tb
    )


def _ours(smi, T, tb):
    """log10(P/atm), or None when aeroprops refuses the molecule.

    The T > Tb extrapolation warning (tested in gate 7) is silenced:
    two controls sit above their boiling points."""
    import warnings

    my = MyrdalYalkowsky(smi=smi, tb=tb)
    if not my.success:
        return None
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        return my.log10_vapour_pressure(T) - LOG10_ATM_TO_PA


# --- all 19 paper Table 1 rows + control battery: agree with UMS -------
# (SMILES, Tb [K], T [K], live-recorded UMS log10 atm)
EXACT = [
    ("CCCCCCCCCC=C", 466.0, 379.0, -1.196438115047803),
    ("CCc1ccc(C)nc1", 452.0, 325.0, -2.053945891688173),
    ("CCc1cc(O)cc(CC)c1", 521.0, 388.0, -1.958240373472998),
    ("C1CCC2CCCCC2C1", 469.0, 373.0, -1.292429096290343),
    ("Oc1ccc2ccccc2c1", 556.0, 415.0, -1.927866419231353),
    ("CC(C)C1CCCCC1", 428.0, 344.0, -1.224567156560805),
    ("COC(=O)c1ccccc1O", 494.0, 352.0, -2.333145747436456),
    ("CCNc1ccccc1", 478.0, 323.0, -2.604545241606901),
    ("Oc1cccc(O)c1", 549.0, 425.0, -1.760973888510557),
    ("CCCCC(=O)O", 446.0, 345.0, -1.726163460295989),
    ("CCCCCCCCCO", 486.0, 365.0, -1.959961397870369),
    ("C1COCCO1", 374.0, 303.0, -1.167217529892387),
    ("CCCCCCS", 426.0, 354.0, -1.026832097069824),
    ("CCCCCCCCCCCCCCO", 569.0, 425.0, -2.011990613163108),
    ("CCCCCCCCCCCCCCCl", 570.0, 415.0, -2.114143235563706),
    ("c1ccccc1c1ccccc1", 528.0, 342.0, -2.965109897784553),
    ("CC1(C2CCC(=C)C1C2)C", 439.0, 300.0, -2.475504827694206),
    ("CC(C)CC(C)(C)CC", 404.0, 323.0, -1.281083447249223),
    ("C" * 28, 705.0, 408.0, -5.009233300527750),
    # controls: alcohols, aromatics, amines, acids, ester, ether,
    # peracid, sulfide, sulfoxide, pyridine
    ("CCO", 351.5, 298.15, -1.162147359548780),
    ("CO", 337.8, 298.15, -0.945875104765052),
    ("c1ccc2ccccc2c1", 491.0, 298.15, -3.609655346547967),
    ("Cc1ccccc1", 383.8, 298.15, -1.457201425264114),
    ("CCN", 289.8, 298.15, 0.139468178145859),
    ("Nc1ccccc1", 457.5, 298.15, -3.048186028109160),
    ("OCC(O)CO", 563.0, 298.15, -6.525834226667049),
    ("OC(=O)CCCCC(=O)O", 610.0, 298.15, -7.250717828725447),
    ("CCOCC", 307.7, 298.15, -0.148041576936148),
    ("CCOC(=O)CC", 350.2, 298.15, -0.863356983942789),
    ("CC(=O)OO", 380.0, 298.15, -1.659834636080126),
    ("CSC", 311.0, 298.15, -0.198313396532082),
    ("CS(C)=O", 462.0, 298.15, -3.000243379089377),
    ("n1ccccc1", 388.4, 298.15, -1.543229590164430),
]

_EXACT_IDS = [
    "undecene", "methylpyridine", "diethylphenol", "decalin", "naphthol",
    "iprcyclohexane", "salicylate", "ethylaniline", "resorcinol",
    "valeric", "nonanol", "dioxane", "hexanethiol", "tetradecanol",
    "chlorotetradecane", "biphenyl", "betapinene", "trimethylhexane",
    "octacosane",
    "ethanol", "methanol", "naphthalene", "toluene", "ethylamine",
    "aniline", "glycerol", "adipic", "ether", "ethylacetate",
    "peracetic", "dimethylsulfide", "dmso", "pyridine",
]


@pytest.mark.parametrize(
    "smi,tb,T,ums_val", EXACT, ids=_EXACT_IDS
)
def test_agrees_with_ums(ums, smi, tb, T, ums_val):
    """Live UMS must reproduce the recorded value (drift guard) and
    aeroprops must match it - bit-for-bit or at the last ulp (the two
    implementations accumulate Eq. 12 in different orders)."""
    assert abs(_ums(ums, smi, T, tb) - ums_val) < 1e-12, (
        smi,
        "UMS drifted",
    )
    ours = _ours(smi, T, tb)
    assert ours is not None, smi
    assert abs(ours - ums_val) < 1e-9, smi


# --- documented divergences ---------------------------------------------
# (name, SMILES, Tb [K], T [K], live-recorded UMS log10 atm, recorded
#  aeroprops log10 atm or None when refused, reason)
WL = [
    (
        "anthracene",
        "c1ccc2cc3ccccc3cc2c1", 615.0, 298.15,
        -6.390881167010023, -6.380965402536875,
        "UMS feeds the SSSR ring count into tau (3 rings -> 0.5); the "
        "paper's RING is independent ring systems: one fused system, "
        "tau = 0",
    ),
    (
        "phenanthrene",
        "c1ccc2c(c1)ccc3ccccc23", 615.0, 298.15,
        -6.390881167010023, -6.380965402536875,
        "same SSSR-vs-ring-systems divergence as anthracene (SSSR 3 "
        "rings, one fused system)",
    ),
    (
        "tert-butyl hydroperoxide",
        "CC(C)(C)OO", 380.0, 298.15,
        -1.622534901654939, -1.613272645899507,
        "UMS counts the tert-butyl centre (tau = 1); the paper states "
        "tert-butyl groups do not contribute to flexibility (tau = 0)",
    ),
    (
        "diethyl sulfide",
        "CCSCC", 308.0, 298.15,
        -0.152020054341123, -0.152770641417282,
        "UMS omits S atoms from the SP3 walk (tau = 1); paper Eq. 8 "
        "explicitly includes S ('NH, N, O, and S', tau = 2)",
    ),
    (
        "prop-1-en-1-ol (enol)",
        "CC=CO", 350.0, 298.15,
        -1.069735868121966, None,
        "UMS's OH family counts the enol OH and its sp2 table accepts "
        "the C=C-C=O pattern, returning a number; aeroprops refuses "
        "enols - outside the published basis (Compernolle 2010 2.6)",
    ),
    (
        "ethyl nitrate",
        "CCO[N+](=O)[O-]", 350.0, 298.15,
        -0.846970307540820, None,
        "UMS's wide tables accept organic nitrates; aeroprops refuses "
        "them - outside the published basis (Compernolle 2010 2.6)",
    ),
]


@pytest.mark.parametrize(
    "name,smi,tb,T,ums_val,ours_val,reason",
    WL,
    ids=[e[0].split(" ")[0] for e in WL],
)
def test_whitelist(ums, name, smi, tb, T, ums_val, ours_val, reason):
    """Documented deviations: the live oracle must reproduce the
    recorded UMS value (drift guard) and aeroprops must reproduce the
    recorded paper-faithful side (or the recorded refusal)."""
    assert abs(_ums(ums, smi, T, tb) - ums_val) < 1e-12, (
        name,
        "UMS drifted",
    )
    ours = _ours(smi, T, tb)
    if ours_val is None:
        assert ours is None, name
        my = MyrdalYalkowsky(smi=smi, tb=tb)
        assert my.status.startswith("unsupported functional group"), name
    else:
        assert abs(ours - ours_val) < 1e-12, name
    assert reason


def test_no_ums_melting_point_input(ums):
    """The oracle API takes no melting point: every comparison above is
    the sub-cooled liquid branch, which is both the paper's Table 1
    convention and aeroprops' Tm default."""
    import inspect

    params = inspect.signature(
        ums.vapour_pressures.myrdal_and_yalkowsky
    ).parameters
    assert list(params) == ["compound", "temperature", "boiling_point"]
