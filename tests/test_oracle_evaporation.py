# SPDX-License-Identifier: MIT
"""Verification gate 5: behavioural comparison against UManSysProp.

Runs the *unmodified* UManSysProp ``evaporation`` implementation live as
the behavioural oracle only (none of its SMARTS or data content is
copied into aeroprops - the aeroprops fragmentation was written
clean-room from the ACP 2011 paper text) and compares upstream-mode
group counts and log10 vapour pressures with
:func:`aeroprops.vapour_pressure.evaporation_fragmentation` and
:func:`aeroprops.methods.evaporation.evaporation_log10_vapour_pressure`.

The divergence whitelist is empty: every probe molecule reproduces the
live oracle's assembled counts exactly and its value (on the oracle's
log10 atm basis) to better than 1e-12; the residual is one ulp of
floating-point association order. Tests that need the UManSysProp tree
are skipped when it is absent (``AEROPROPS_UMS_REPO``); the recorded
drift lock below runs regardless.
"""

import pytest

import ums_oracle
from aeroprops.methods.evaporation import (
    LOG10_ATM_TO_PA,
    evaporation_log10_vapour_pressure,
)
from aeroprops.structure import mol_from_smiles
from aeroprops.vapour_pressure import (
    _evaporation_assemble_upstream,
    evaporation_fragmentation,
)

live = pytest.mark.skipif(
    not ums_oracle.available(), reason="UManSysProp source tree not available"
)


def _ours(smi):
    counts, success, status = evaporation_fragmentation(
        mol_from_smiles(smi), upstream_compat=True)
    assert success, status
    return counts


def _live_assembled(ums, smi):
    """Assembled 22-key counts as produced by the live UMS oracle."""
    compound = ums_oracle.readstring(ums, smi)
    raw = {
        str(k): v
        for k, v in ums.groups.matches(
            ums.groups.data.EVAPORATION_SMARTS, compound
        ).items()
    }
    return _evaporation_assemble_upstream(raw)


#: 51 probe molecules covering every decoded counting rule, with the
#: assembled upstream-mode counts recorded from the live oracle
#: (2026-11-10, UManSysProp public master). Key "1" (carbon count,
#: which the oracle reports but aeroprops drops) removed; zero counts
#: omitted.
CONTROLS = [
    ('1-2-dioxane', 'C1COOCC1', {'12': 1, '12CL': 1, '12lin': 1, '2': 6, '3': -1.0}),
    ('1-4-dioxane', 'C1COCCO1', {'12': 2, '12lin': 2, '2': 6, '3': -1.0}),
    ('2-3-butanediol', 'CC(O)C(O)C', {'14': 2, '19': 2, '2': 4, '8': 2}),
    ('2-3-pentanedione', 'CC(=O)C(=O)CC', {'16': 2, '2': 5, '5': 2}),
    ('2-4-pentanedione', 'CC(=O)CC(=O)C', {'17': 2, '2': 5, '5': 2}),
    ('2-ethoxyethanol', 'CCOCCO', {'19': 1, '2': 5, '8': 1}),
    ('2-hydroxyethyl-acetate', 'CC(=O)OCCO', {'19': 1, '2': 5, '6': 1, '8': 1}),
    ('2-methoxy-2-propanol', 'OC(C)(C)OC', {'14': 2, '2': 5, '3': 1.0, '8': 1}),
    ('3-hydroxy-2-4-pentanedione', 'CC(=O)C(O)C(=O)C',
     {'14': 1, '17': 2, '18': 2, '19': 2, '2': 5, '5': 2, '8': 1}),
    ('3-methylpentane', 'CCC(C)CC', {'2': 6, '3': 1.0}),
    ('acrolein', 'C=CC=O', {'13': 1, '2': 3, '5': 1}),
    ('adipic-acid', 'O=C(O)CCCCC(=O)O', {'2': 6, '9': 2}),
    ('allyl-alcohol', 'C=CCO', {'15': 1, '2': 3, '8': 1}),
    ('alpha-pinene', 'CC1=CCC2CC1C2(C)C', {'2': 10, '3': 6.5}),
    ('butane', 'CCCC', {'2': 4}),
    ('citric-acid', 'OC(=O)CC(O)(CC(=O)O)C(=O)O',
     {'14': 2, '19': 1, '2': 6, '3': 1.0, '8': 1, '9': 3}),
    ('cyclohexane', 'C1CCCCC1', {'2': 6, '3': -1.0}),
    ('cyclohexanone', 'O=C1CCCCC1', {'12': 1, '12CL': 1, '2': 6, '3': -1.0, '5': 1}),
    # sqrt-branch k12 coverage: n_HB > 1 with typed 12HB ring hydroxyls
    # (the A parameter takes k12 only through its typed split there)
    ('cyclohexanetriol', 'OC1CC(O)CCC1O',
     {'12': 3, '12HB': 3, '14': 3, '19': 2, '2': 6, '3': -1.0, '8': 3}),
    ('ethyl-acetate', 'CC(=O)OCC', {'2': 5, '6': 1}),
    ('ethyl-pyruvate', 'CCOC(=O)C(=O)C', {'16': 1, '2': 6, '5': 1, '6': 1}),
    ('ethylene-glycol', 'OCCO', {'19': 2, '2': 2, '8': 2}),
    ('fenchyl-alcohol', 'C1(C)(C)C(O)C2(C)CCC1C2',
     {'12': 1, '12HB': 1, '14': 1, '2': 10, '3': 11.6, '8': 1}),
    ('formaldehyde', 'C=O', {'2': 1, '5': 1}),
    ('gamma-butyrolactone', 'O=C1OCCC1', {'12': 1, '2': 5, '3': -1.0, '6': 1}),
    ('glycolic-acid', 'OCC(=O)O', {'19': 1, '2': 2, '8': 1, '9': 1}),
    ('glyoxal', 'O=CC=O', {'16': 1, '2': 2, '5': 2}),
    ('glyoxylic-acid', 'O=CC(=O)O', {'2': 2, '20': 1, '5': 1, '9': 1}),
    ('hydroxyacetone', 'CC(=O)CO', {'18': 1, '19': 1, '2': 3, '5': 1, '8': 1}),
    ('hydroxyethyl-hydroperoxide', 'OCCOO', {'10': 1, '19': 1, '2': 2, '8': 1}),
    ('hydroxyethyl-nitrate', 'OCCO[N+](=O)[O-]', {'19': 1, '2': 2, '8': 1}),
    # sqrt-branch k12 coverage: n_CL > 1 where the ester acyl carbon is
    # on the ring (upstream submode 602: aggregate "12" only, no typed
    # split - it drops out of the A parameter entirely there)
    ('lactone-ketone', 'O=C1OCC(=O)CC1',
     {'12': 2, '12CL': 1, '18': 1, '2': 6, '3': -1.0, '5': 1, '6': 1}),
    ('mbtca', 'CC(C)(C(=O)O)C(C(=O)O)CC(=O)O', {'2': 8, '3': 7.0, '9': 3}),
    ('methane', 'C', {'2': 1}),
    ('methoxyacetaldehyde', 'COCC=O', {'18': 1, '2': 4, '5': 1}),
    ('methyl-formate', 'COC=O', {'2': 2}),
    ('methyl-nitrate', 'CO[N+](=O)[O-]', {'2': 1}),
    ('methyl-pyruvate', 'COC(=O)C(=O)C', {'16': 1, '2': 5, '5': 1, '6': 1}),
    ('mtbe', 'CCOC(C)(C)C', {'2': 7, '3': 2.0}),
    ('mvk', 'C=CC(=O)C', {'13': 1, '2': 4, '5': 1}),
    ('neopentane', 'C(C)(C)(C)C', {'2': 5, '3': 6.0}),
    ('oxalic-acid', 'O=C(O)C(=O)O', {'2': 2, '9': 2}),
    ('pan', 'CC(=O)OO[N+](=O)[O-]', {'2': 2}),
    ('peracetic-acid', 'CC(=O)OO', {'11': 1, '2': 2}),
    ('pyruvic-acid', 'CC(=O)C(=O)O', {'2': 3, '20': 1, '5': 1, '9': 1}),
    ('succinic-acid', 'O=C(O)CC(=O)O', {'2': 3, '9': 2}),
    ('t-butanol', 'CC(O)(C)C', {'14': 2, '2': 4, '3': 1.0, '8': 1}),
    ('t-butyl-hydroperoxide', 'CC(C)(C)OO', {'10': 1, '2': 4, '3': 2.0}),
    ('tartaric-acid', 'OC(=O)C(O)C(O)C(=O)O', {'14': 2, '19': 4, '2': 4, '8': 2, '9': 2}),
    ('tetrahydrofuran', 'C1CCOC1', {'12': 1, '12lin': 1, '2': 5, '3': -1.0}),
    ('vinyl-alcohol', 'C=CO', {'2': 2, '8': 1}),
]


def _assert_same_counts(recorded, produced, label):
    # key "1" is the plain carbon count the oracle reports but
    # aeroprops drops from its counts dict
    for k in set(recorded) | set(produced):
        if k == "1":
            continue
        assert recorded.get(k, 0) == pytest.approx(
            produced.get(k, 0), abs=1e-12), (label, k)


@live
def test_oracle_example_output_counts_and_log10p(ums):
    """All 26 Example_output.txt SMILES: identical upstream-mode counts
    and identical log10 vapour pressure on the oracle's atm basis."""
    smiles = ums_oracle.example_output_smiles()
    assert len(smiles) == 26
    worst = 0.0
    for smi in smiles:
        want = _live_assembled(ums, smi)
        ours = _ours(smi)
        _assert_same_counts(want, ours, smi)
        live_atm = ums.vapour_pressures.evaporation(
            ums_oracle.readstring(ums, smi), 298.15)
        ours_atm = evaporation_log10_vapour_pressure(
            ours, 298.15, True) - LOG10_ATM_TO_PA
        worst = max(worst, abs(ours_atm - live_atm))
    assert worst < 1e-12


@pytest.mark.parametrize(
    "name, smi, recorded", CONTROLS, ids=[row[0] for row in CONTROLS]
)
def test_control_drift_lock(name, smi, recorded):
    """aeroprops must keep reproducing the recorded oracle counts even
    without the UManSysProp tree (regression lock on every decoded
    counting rule)."""
    _assert_same_counts(recorded, _ours(smi), name)


@live
@pytest.mark.parametrize(
    "name, smi, recorded", CONTROLS, ids=[row[0] for row in CONTROLS]
)
def test_control_matches_live_oracle(ums, name, smi, recorded):
    """Two-sided drift lock: the live oracle must still produce the
    recorded counts (upstream change detection) and aeroprops must
    match the live oracle exactly."""
    live_counts = _live_assembled(ums, smi)
    _assert_same_counts(recorded, live_counts, (name, "oracle drifted"))
    _assert_same_counts(live_counts, _ours(smi), name)

