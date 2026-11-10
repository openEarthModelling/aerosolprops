# SPDX-License-Identifier: MIT
"""Verification gates for the Nannoolal Part 1 implementation.

Gate 1  coefficient triple-guard: independent re-transcription of paper
        Table 1/3 rows against the committed catalogue, including the
        two documented transcription rulings (IDs 12 and 24).
Gate 2  paper Tables 17a-17d worked examples end-to-end (SMILES ->
        fragmentation -> Tb).
Gate 3  literature normal boiling points (CRC Handbook values).
Gate 4  atom balance over a broad molecule battery: every non-H atom
        claimed exactly once, n consistent.
Gate 6  priority-number adjudication on the real catalogue (fluorine
        and chlorine families, ester vs lactone, epoxide vs ether).
Gate 7  the no-GI extension branch and the two printed-without-value
        interactions (Table 5 rows 153/211) refuse, never zero.

Gate 5 (live UManSysProp behaviour comparison) lives in
``test_oracle_umansysprop.py``.
"""

import pytest

from aeroprops.boiling_point import (
    NANNOOLAL_BP_SMARTS,
    NannoolalBP,
    _nannoolal_engine,
    nannoolal_fragmentation,
)
from aeroprops.methods.nannoolal_bp import (
    GI_STATUS_NO_GI,
    GI_STATUS_OK,
    NANNOOLAL_BP_BASE,
    NANNOOLAL_BP_EXTENSIONS,
    NANNOOLAL_BP_GROUPS,
    NANNOOLAL_BP_INTERACTION_IDS,
    NANNOOLAL_BP_INTERACTIONS,
    nannoolal_bp_gi,
)
from aeroprops.structure import mol_from_smiles


# --- Gate 1: independent transcription spot-check --------------------
# Rows transcribed a second time from the paper text layer (Tables 1
# and 2; statistics cross-checked against Table 3/4 lines), before
# comparing with the committed catalogue.
PAPER_ROWS = [
    # (id, name, contribution [K], priority or None)
    (1, "CH3-(ne)", 177.3066, 104),
    (2, "CH3-(e)", 251.8338, 102),
    (3, "CH3-(a)", 157.9527, 103),
    (4, "C(c)H2", 239.4531, 111),
    (5, "C(c)H", 240.6785, 117),
    (6, "C(c)", 249.5809, 119),
    (7, "C(c)(e)", 266.8769, 107),
    (8, "C(c)(a)", 201.0115, 108),
    (9, "C(r)H2", 239.4957, 112),
    (10, "C(r)H", 222.1163, 116),
    (11, "C(r)", 209.9749, 118),
    (13, "C(r)(e,r)", 291.2291, 110),
    (14, "C(r)(Ca)", 244.3581, 106),
    (15, "=C(a)H", 235.3462, 105),
    (16, "=C(a)(ne)", 315.4128, 115),
    (17, "=C(a)(e)", 348.2779, 113),
    (18, "(a)=C(a)2(a)", 367.9649, 114),
    (58, "C(c)=C(c)", 475.7958, 63),
    (59, "C(c)=C(c)(C(a))", 586.1413, 60),
    (60, "-(e)C(c)=C(c)", 500.2434, 59),
    (61, "H2C(c)=C (1-ene)", 412.6276, 58),
    (62, "C(r)=C(r)", 475.9623, 61),
    (63, "C≡C", 512.2893, 62),
    (64, "HC≡C (1-yne)", 422.2307, 57),
    (87, "C=C=C (cumulated)", 664.0903, 6),
    (88, "C=C-C=C (conjugated, ring)", 957.6388, 7),
    (89, "C=C-C=C (conjugated, chain)", 928.9954, 8),
    (95, "C≡C-C≡C (conjugated)", 1218.1878, 9),
    (118, "C=C-C=O", 40.4205, None),
    (133, "C3C-CC3", 111.8372, None),
]


@pytest.mark.parametrize("gid,name,value,priority", PAPER_ROWS)
def test_catalogue_matches_independent_transcription(gid, name, value,
                                                     priority):
    assert NANNOOLAL_BP_GROUPS[gid] == (name, value, priority)


def test_catalogue_structure_counts():
    # 132 groups = 113 paper main + 16 corrections + 3 extensions.
    assert len(NANNOOLAL_BP_GROUPS) == 132
    assert sorted(NANNOOLAL_BP_EXTENSIONS) == [301, 302, 303]
    # 76 valued interactions of the 78 printed rows.
    assert len(NANNOOLAL_BP_INTERACTIONS) == 76
    assert len(NANNOOLAL_BP_INTERACTION_IDS) == 78
    assert NANNOOLAL_BP_BASE == {"a": 0.6583, "b": 1.6868, "c": 84.3395}


def test_documented_transcription_rulings():
    # Table 1 misprint 492.0707 (duplicate of ID 109) resolved to the
    # Table 3 regression value 250.9584 (method note 04, OQ-1/OQ-15).
    assert NANNOOLAL_BP_GROUPS[12][1] == 250.9584
    # Table 1 omits the minus sign of ID 24; Table 3 prints -19.5575.
    assert NANNOOLAL_BP_GROUPS[24][1] == -19.5575


def test_every_catalogued_smarts_group_has_a_priority():
    # The engine cannot order PR-less groups; only corrections and
    # extensions are PR-less and they must not appear in the catalogue.
    engineless = [
        gid
        for gid in NANNOOLAL_BP_SMARTS
        if NANNOOLAL_BP_GROUPS[gid][2] is None
        and gid not in (301, 302, 303)
    ]
    assert engineless == []
    assert 109 not in NANNOOLAL_BP_SMARTS  # unimplementable, see docs


