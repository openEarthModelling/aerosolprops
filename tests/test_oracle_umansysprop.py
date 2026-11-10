# SPDX-License-Identifier: MIT
"""Verification gate 5: behavioural comparison against UManSysProp.

Runs the *unmodified* UManSysProp ``simpol`` implementation live (as the
behavioural oracle only - none of its data or SMARTS content is copied
into aerosolprops) through the RDKit-based ``openbabel.pybel`` /
``pkg_resources`` shim from the UManSysProp RDKit tutorial, and compares
group counts and log10 vapour pressures with
:func:`aerosolprops.vapour_pressure.simpol1_fragmentation`.

Skipped entirely when the UManSysProp source tree is not present (set
``AEROSOLPROPS_UMS_REPO`` to point at a checkout; the default is the local
path used during development).
"""

import pytest

import ums_oracle
from aerosolprops.methods.simpol1 import LOG10_ATM_TO_PA, simpol1_bk
from aerosolprops.structure import mol_from_smiles
from aerosolprops.vapour_pressure import simpol1_fragmentation

USP_REPO = ums_oracle.USP_REPO

pytestmark = pytest.mark.skipif(
    not ums_oracle.available(), reason="UManSysProp source tree not available"
)


def _ours(smi):
    counts, success, status = simpol1_fragmentation(mol_from_smiles(smi))
    assert success, status
    return counts


def _readstring_smi(smi):
    return ums_oracle.readstring(None, smi)


def _ums_log10p_atm(ums_modules, smi, T=298.15):
    return ums_modules.vapour_pressures.simpol(_readstring_smi(smi), T)


def _ours_log10p_atm(counts, T=298.15):
    # aerosolprops returns log10 Pa; subtract the conversion constant to
    # compare on the paper's (and UMS's) log10 atm basis.
    total = simpol1_bk(0, T)
    for k, nu in counts.items():
        if k and nu:
            total += nu * simpol1_bk(k, T)
    return total


def _example_output_smiles():
    return ums_oracle.example_output_smiles()


#: Known UManSysProp deviations from the paper, whitelisted per the
#: method notes section 6 / open questions. Each entry records the UMS
#: group counts and the aerosolprops counts (source of the aerosolprops
#: values: live run on 2026-11-10 against UManSysProp public master).
WHITELIST = [
    # name, SMILES, {k: (ums, ours)}, justification
    (
        "acetamide",
        "CC(N)=O",
        {2: (1, 2)},
        "notes 3.2/6.1: UMS group-2 SMARTS under-counts methylamide "
        "(paper semantics: 1 methyl + 1 carbonyl carbon = 2)",
    ),
    (
        "n-propylbutamide",
        "CCCC(NCCC)=O",
        {2: (5, 4)},
        "notes 3.2: paper example nu_2 = 4 (butyl side only); UMS "
        "counts 5; the raw-notes value of 10 is the non-deduplicated "
        "OpenBabel output, the uniquified count is 5",
    ),
    (
        "DMF",
        "CN(C)C=O",
        {2: (0, 1), 24: (0, 1)},
        "notes 6.3: UMS misses both the tertiary amide group and its "
        "acid-side methyl for N(C)=O written without an explicit H "
        "(paper: DMF is a tertiary amide, nu_2 = 1)",
    ),
    (
        "formamide",
        "NC=O",
        {2: (0, 1), 22: (0, 1)},
        "notes 6.3: UMS finds no primary amide when the carbonyl carbon "
        "bears two hydrogens",
    ),
    (
        "N-methylformamide",
        "CNC=O",
        {2: (0, 1), 23: (0, 1)},
        "notes 6.3: same formyl edge case, secondary amide missed",
    ),
    (
        "N-methylacetamide",
        "CNC(C)=O",
        {2: (1, 2)},
        "notes 3.2: same under-count as acetamide",
    ),
    (
        "piperidin-2-one",
        "O=C1CCCCN1",
        {2: (4, 5)},
        "OQ-1: lactam acid-side carbons undefined in the paper; "
        "aerosolprops counts the whole acyl-side component (5 ring "
        "carbons), UMS's fixed-length pattern captures 4",
    ),
    (
        "1-methyl-piperidin-2-one",
        "CN1CCCCC1=O",
        {2: (4, 6)},
        "OQ-1: as above; the ring stays connected after the C-N cut, so "
        "the N-methyl carbons are included (6)",
    ),
    (
        "methyl benzoate",
        "COC(=O)c1ccccc1",
        {11: (0, 1)},
        "notes 6.2: UMS misses aromatic (benzoate-type) esters entirely",
    ),
    (
        "methyl 2-nitropropionate",
        "COC(=O)C(C)[N+](=O)[O-]",
        {11: (1, 0), 30: (0, 1)},
        "notes 6.8: the nitroester group 30 must replace the plain "
        "ester 11 (alpha-carbon nitro); UMS keeps 11 and misses 30",
    ),
    (
        "2,4-dinitrophenol",
        "Oc1ccc([N+](=O)[O-])cc1[N+](=O)[O-]",
        {29: (2, 1)},
        "OQ-2: group 29 is per-OH (one nitrophenol unit regardless of "
        "the number of nitro ringsmates); UMS double-counts",
    ),
    # --- additional deviations found by the extended amide/ester probe
    # --- (same classes as above, live run 2026-11-10)
    (
        "ethyl benzoate",
        "CCOC(=O)c1ccccc1",
        {11: (0, 1)},
        "notes 6.2 (benzoate class): UMS misses aromatic esters of "
        "benzoic acid, ethyl benzoate behaves like methyl benzoate",
    ),
    (
        "N,N-dimethylacetamide",
        "CC(=O)N(C)C",
        {2: (1, 2)},
        "notes 3.2 (methylamide under-count class): acid side = acetyl "
        "(2 carbons incl. the carbonyl C); UMS counts 1",
    ),
    (
        "succinamide",
        "O=C(N)CCC(N)=O",
        {2: (3, 8)},
        "OQ-1 family: two primary amides on one chain; acid-side is "
        "counted per amide (each sees the whole 4-carbon chain, so "
        "nu_2 = 4 + 4 = 8); UMS counts 3",
    ),
    (
        "propionamide",
        "CCC(=O)N",
        {2: (2, 3)},
        "notes 3.2/6.1 class (acetamide): acid side includes the "
        "carbonyl carbon (propionyl = 3 C); UMS omits it",
    ),
    (
        "2-acetamidoacetamide",
        "NC(=O)CNC(C)=O",
        {2: (2, 6)},
        "per-amide acid-side graph cut: the terminal formyl amide sees "
        "4 carbons, the acetamide end sees 2 (nu_2 = 6); UMS counts 2",
    ),
]


def test_oracle_example_output_counts_and_log10p(ums):
    """All 26 Example_output.txt SMILES: identical counts, identical
    log10 vapour pressure on the atm basis."""
    smiles = _example_output_smiles()
    assert len(smiles) == 26
    worst = 0.0
    for smi in smiles:
        raw = dict(zip([str(k) for k in range(31)], [0] * 31))
        raw.update(ums.groups.simpol(_readstring_smi(smi)))
        ums_counts = {int(k): v for k, v in raw.items() if v}
        ours = _ours(smi)
        assert ours == ums_counts, smi
        delta = abs(
            _ums_log10p_atm(ums, smi) - _ours_log10p_atm(ours)
        )
        worst = max(worst, delta)
    assert worst < 1e-12  # identical counts must give identical values


def test_oracle_amide_probe_set_counts_match(ums):
    """Amide/ester probe molecules that are NOT whitelisted must match
    exactly (guards the whitelist against both stale entries and new
    regressions)."""
    whitelisted = {entry[1] for entry in WHITELIST}
    probe = [
        # non-whitelisted amides/esters and controls:
        "COC(=O)CC(N)=O",  # ester + primary amide, both counted by UMS
        "OC(=O)c1ccccc1",  # benzoic acid (aromatic COOH is found by UMS)
        "CC(=O)NCC",  # N-ethylacetamide: acid side = acetyl (2), agrees
        "CC1(C2CCC(O)C1C2)C",  # fenchyl alcohol (bridged ring control)
    ]
    matched = 0
    for smi in probe:
        assert smi not in whitelisted
        raw = ums.groups.simpol(_readstring_smi(smi))
        ums_counts = {int(k): v for k, v in raw.items() if v}
        ours = _ours(smi)
        assert ours == ums_counts, smi
        matched += 1
    assert matched >= 3


@pytest.mark.parametrize(
    "name, smi, deviations, justification",
    WHITELIST,
    ids=[entry[0] for entry in WHITELIST],
)
def test_oracle_whitelisted_deviations(ums, name, smi, deviations, justification):
    """Documented UMS deviations: the live UMS oracle must reproduce the
    recorded UMS values (guards against oracle drift), and aerosolprops
    must reproduce the paper-semantics values."""
    raw = ums.groups.simpol(_readstring_smi(smi))
    ums_counts = {int(k): v for k, v in raw.items() if v}
    ours = _ours(smi)
    for k, (ums_value, ours_value) in deviations.items():
        assert ums_counts.get(k, 0) == ums_value, (name, k, "UMS drifted")
        assert ours.get(k, 0) == ours_value, (name, k)
    # Every other group must agree (the deviation is local).
    for k in set(ums_counts) | set(ours):
        if k in deviations:
            continue
        assert ums_counts.get(k, 0) == ours.get(k, 0), (name, k)


def test_oracle_whitelist_covers_all_amide_diffs(ums):
    """No unexplained diffs: over the full probe set, any count mismatch
    must belong to a whitelisted (molecule, group) pair."""
    allowed = {(entry[1], k) for entry in WHITELIST for k in entry[2]}
    molecules = {entry[1] for entry in WHITELIST} | {
        "COC(=O)CC(N)=O",
        "OC(=O)c1ccccc1",
        "CC(=O)NCC",
    }
    for smi in molecules:
        raw = ums.groups.simpol(_readstring_smi(smi))
        ums_counts = {int(k): v for k, v in raw.items() if v}
        ours = _ours(smi)
        for k in set(ums_counts) | set(ours):
            if ums_counts.get(k, 0) != ours.get(k, 0):
                assert (smi, k) in allowed, (smi, k)
