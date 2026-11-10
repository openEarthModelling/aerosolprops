# SPDX-License-Identifier: MIT
"""Verification gate 7: the MyrdalYalkowsky wrapper.

The three Tb routes (Nannoolal Part 1, thermo's Joback, experimental
passthrough), the Tm default (sub-cooled liquid, annotated in status),
the T > Tb warning, and the failure paths (bad SMILES, refused
functional groups, refused elements, invalid Tb/Tm/sigma, route
refusals) - every failure leaves ``success == False`` and both
estimators returning ``None``.
"""

import pytest

from aerosolprops.methods.simpol1 import LOG10_ATM_TO_PA
from aerosolprops.structure import mol_from_smiles
from aerosolprops.vapour_pressure import MyrdalYalkowsky

# ---------------------------------------------------------------------------
# End to end: paper Table 1 rows from SMILES through the wrapper
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "smi,tb,T,tm,sigma,pred",
    [
        ("Oc1cccc(O)c1", 549.0, 425.0, None, 1, -1.76),  # resorcinol
        ("CCc1cc(O)cc(CC)c1", 521.0, 388.0, None, 1, -1.96),
        ("CCCCC(=O)O", 446.0, 345.0, None, 1, -1.73),    # valeric acid
        ("CCc1ccc(C)nc1", 452.0, 325.0, None, 1, -2.05),
        ("c1ccccc1c1ccccc1", 528.0, 342.0, 342.4, 2, -2.97),  # solid branch
        ("C" * 28, 705.0, 408.0, None, 1, -5.01),        # octacosane
    ],
    ids=["resorcinol", "diethylphenol", "valeric", "methylpyridine",
         "biphenyl-solid", "octacosane"],
)
def test_table1_end_to_end(smi, tb, T, tm, sigma, pred):
    """SMILES in, paper's predicted column out (2 dp, log10 atm)."""
    my = MyrdalYalkowsky(smi=smi, tb=tb, tm=tm, sigma=sigma)
    assert my.success, my.status
    lp = my.log10_vapour_pressure(T)
    assert lp is not None
    assert round(lp - LOG10_ATM_TO_PA, 2) == pred


def test_experimental_tb_passthrough():
    my = MyrdalYalkowsky(smi="Oc1cccc(O)c1", tb=549.0, tm=342.0)
    assert my.success
    assert my.tb == 549.0
    assert my.tb_source == "experimental"
    assert my.status == "OK"
    d = my.descriptors
    from aerosolprops.methods.myrdal_yalkowsky import (
        myrdal_yalkowsky_log10_vapour_pressure,
    )

    assert my.log10_vapour_pressure(425.0) == (
        myrdal_yalkowsky_log10_vapour_pressure(
            425.0, 549.0, d["tau"], d["hbn"], Tm=342.0, sigma=1.0
        )
    )


# ---------------------------------------------------------------------------
# Tb route 1: Nannoolal Part 1 (default)
# ---------------------------------------------------------------------------


def test_nannoolal_route_default():
    from aerosolprops.boiling_point import NannoolalBP

    my = MyrdalYalkowsky(smi="CCCCCC")
    assert my.success, my.status
    assert my.tb_source == "nannoolal"
    nnp = NannoolalBP(smi="CCCCCC")
    assert my.tb == nnp.boiling_point()
    # live-recorded drift locks (Nannoolal Part 1 through aerosolprops)
    assert round(my.tb, 2) == 350.03
    assert round(my.log10_vapour_pressure(298.15), 2) == 4.14


def test_nannoolal_route_extension_groups_still_estimate():
    """Compernolle extension groups (hydroperoxide) lose the GI term but
    still produce a Tb - the reduced-confidence branch flows through."""
    my = MyrdalYalkowsky(smi="CC(C)(C)OO")
    assert my.success, my.status
    assert my.tb is not None
    assert round(my.tb, 2) == 402.9


def test_nannoolal_route_refused():
    """Carbon disulfide: an unclaimed sulfur atom refuses Nannoolal,
    which fails the whole wrapper (no silent Tb fallback)."""
    my = MyrdalYalkowsky(smi="S=C=S")
    assert not my.success
    assert my.tb is None
    assert my.status.startswith("Nannoolal Tb refused")
    assert my.log10_vapour_pressure() is None
    assert my.vapour_pressure() is None


# ---------------------------------------------------------------------------
# Tb route 2: thermo's built-in Joback
# ---------------------------------------------------------------------------


def test_joback_route():
    from thermo import Joback

    my = MyrdalYalkowsky(smi="CCCCCC", tb_source="joback")
    assert my.success, my.status
    assert my.tb_source == "joback"
    assert my.tb == Joback("CCCCCC").Tb_estimated
    assert round(my.tb, 2) == 336.88  # drift lock (thermo 0.6.1)


def test_joback_route_from_mol_object():
    """Without `smi` the route derives canonical SMILES from the mol."""
    my = MyrdalYalkowsky(mol=mol_from_smiles("CCCCCC"), tb_source="joback")
    assert my.success, my.status
    assert round(my.tb, 2) == 336.88


def test_joback_route_gates_on_success():
    """thermo's Joback returns a numeric Tb_estimated even when its
    fragmentation failed (nitric acid: 291.08 with success False); the
    route must refuse, not propagate the number."""
    my = MyrdalYalkowsky(smi="O=[N+]([O-])O", tb_source="joback")
    assert not my.success
    assert my.tb is None
    assert my.status == "Joback Tb unavailable (thermo refused)"


# ---------------------------------------------------------------------------
# Tb route 3: experimental without data, and unknown routes
# ---------------------------------------------------------------------------


def test_experimental_route_without_tb_refused():
    my = MyrdalYalkowsky(smi="CCCCCC", tb_source="experimental")
    assert not my.success
    assert my.tb is None
    assert my.status == "experimental Tb required: pass tb=..."


def test_unknown_tb_source_refused():
    my = MyrdalYalkowsky(smi="CCCCCC", tb_source="bogus")
    assert not my.success
    assert my.status == "unknown tb_source: 'bogus'"


# ---------------------------------------------------------------------------
# Tm handling and the extrapolation warning
# ---------------------------------------------------------------------------


def test_tm_default_is_subcooled_liquid():
    my = MyrdalYalkowsky(smi="c1ccccc1c1ccccc1", tb=528.0)
    assert my.success
    assert my.status == "OK (sub-cooled liquid: Tm not supplied)"
    liquid = my.log10_vapour_pressure(342.0)
    # supplying the melting point activates the solid branch (lower P)
    my.tm = 342.4
    my.sigma = 2.0
    solid = my.log10_vapour_pressure(342.0)
    assert solid < liquid
    assert round(solid - LOG10_ATM_TO_PA, 2) == -2.97  # Table 1


def test_tm_below_T_is_liquid_branch():
    """A supplied Tm <= T leaves the liquid branch (Appendix A): sigma
    then has nothing to multiply."""
    my = MyrdalYalkowsky(smi="c1ccccc1c1ccccc1", tb=528.0, tm=300.0,
                         sigma=2.0)
    assert my.success
    assert my.status == "OK"
    bare = MyrdalYalkowsky(smi="c1ccccc1c1ccccc1", tb=528.0)
    assert my.log10_vapour_pressure(342.0) == (
        bare.log10_vapour_pressure(342.0)
    )


def test_warning_above_boiling_point():
    my = MyrdalYalkowsky(smi="CCCCCC")
    with pytest.warns(UserWarning, match="exceeds Tb"):
        assert my.log10_vapour_pressure(400.0) is not None
    # at and below Tb: silent
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        my.log10_vapour_pressure(350.0)
        my.log10_vapour_pressure(298.15)


def test_invalid_temperature_returns_none():
    my = MyrdalYalkowsky(smi="CCCCCC")
    assert my.log10_vapour_pressure(0.0) is None
    assert my.log10_vapour_pressure(-50.0) is None
    assert my.vapour_pressure(-50.0) is None


# ---------------------------------------------------------------------------
# Failure paths: construction, descriptors, inputs
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "smi,status",
    [
        ("this is not a smiles", "Failed to construct mol"),
        ("CCO[N+](=O)[O-]", "unsupported functional group: organic nitrate"),
        ("[Si](C)(C)(C)C", None),  # startswith "unsupported element(s):"
    ],
    ids=["bad-smiles", "nitrate", "silicon"],
)
def test_descriptor_failures(smi, status):
    my = MyrdalYalkowsky(smi=smi, tb=400.0)
    assert not my.success
    assert my.descriptors == {}
    if status is None:
        assert my.status.startswith("unsupported element(s):")
    else:
        assert my.status == status
    # thermo-style contract: attributes exist on failed objects too
    assert my.tb_source == "experimental"  # tb was given
    assert my.tm is None
    assert my.sigma == 1.0
    assert my.log10_vapour_pressure() is None
    assert my.vapour_pressure() is None


@pytest.mark.parametrize(
    "kw,status",
    [
        ({"tb": 0.0}, "invalid experimental Tb"),
        ({"tb": -5.0}, "invalid experimental Tb"),
        ({"tb": 400.0, "tm": 0.0}, "invalid Tm"),
        ({"tb": 400.0, "tm": -1.0}, "invalid Tm"),
        ({"tb": 400.0, "sigma": 0.5}, "invalid sigma (< 1)"),
    ],
    ids=["tb-zero", "tb-negative", "tm-zero", "tm-negative", "sigma-low"],
)
def test_invalid_inputs(kw, status):
    my = MyrdalYalkowsky(smi="CCCCCC", **kw)
    assert not my.success
    assert my.status == status
    # attributes exist (round-tripping the given inputs) on failures
    assert my.tb_source == "experimental"
    assert my.tm == kw.get("tm")
    assert my.sigma == kw.get("sigma", 1.0)
    assert my.log10_vapour_pressure() is None
    assert my.vapour_pressure() is None
