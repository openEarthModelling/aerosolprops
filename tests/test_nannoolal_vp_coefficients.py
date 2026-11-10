# SPDX-License-Identifier: MIT
"""Verification gate 1 for Nannoolal Part 3: coefficient triple-guard.

1.  An independent second transcription of paper Table 4 (main dB
    contributions x 1e3 with the Table 1 priority numbers), Table 5
    (second-order) and Table 6 (interactions), re-typed from the paper
    text layer (1-s2.0-S0378381208001611), against the committed
    catalogue.
2.  Paper-derivable column cross-checks: single coefficients read off
    the printed per-group rows of the worked examples Tables 15a-15h -
    a source independent of the statistics tables.
3.  Structure rulings locked against the paper: the Table 4 row census
    (IDs 106/108/109/112/114 carry no dB row; adjudication below), the
    defined-but-valueless set, the 98/110 no-definition rows, the
    questionable-value flags, and the interaction-class membership
    agreement with the Part 1 tables (same letters A-S, Part 3's own
    OH/OH(a) split).

Adjudication recorded here (paper wins over the method notes): the
notes' Table 4 census lists the valueless IDs as {106, 108, 112, 114};
the paper text has no row for 109 either ("CO connected to N", the same
unusable definition and valueless entry as in Part 1), so the valueless
census is {106, 108, 109, 112, 114} and 109 is excluded from the
catalogue entirely.
"""

import pytest

from aeroprops.methods.nannoolal_bp import (
    NANNOOLAL_BP_INTERACTION_CLASSES,
    NANNOOLAL_BP_INTERACTION_IDS,
)
from aeroprops.methods.nannoolal_vp import (
    NANNOOLAL_VP_BASE,
    NANNOOLAL_VP_EXTENSIONS,
    NANNOOLAL_VP_GROUPS,
    NANNOOLAL_VP_INTERACTION_CLASSES,
    NANNOOLAL_VP_INTERACTION_IDS,
    NANNOOLAL_VP_INTERACTIONS,
    NANNOOLAL_VP_NO_PARAMETER,
    NANNOOLAL_VP_UNDEFINED_VALUE_GROUPS,
)

# --- independent transcription: paper Table 4 (value x 1e3, PR) -------
# (id, dB x 1e3, priority). Rows re-typed from the Table 4 text; PRs
# from the Table 1 (ID, PR) column pairs. ID 41's PR cell is blank in
# print; 95 is the recorded inference (Part 1 PR + 1, unused otherwise).
# IDs 98/110 print values but have no Table 1 definition: priority None,
# value parked in NANNOOLAL_VP_UNDEFINED_VALUE_GROUPS, never counted.
PAPER_TABLE_4 = [
    (1, 13.3063, 105), (2, 91.8000, 103), (3, 50.1939, 104),
    (4, 54.6564, 112), (5, 45.7437, 119), (6, -31.7531, 121),
    (7, 37.8485, 108), (8, 96.1386, 109), (9, 22.2573, 113),
    (10, 32.8162, 118), (11, 4.8500, 120), (12, 23.6411, 110),
    (13, 49.8237, 111), (14, -3.6950, 107), (15, 32.7177, 106),
    (16, 69.8796, 117), (17, 41.5534, 114), (18, 43.7191, 115),
    (19, 79.5429, 87), (20, 51.2880, 85), (21, 42.0887, 81),
    (22, 56.9998, 84), (23, 142.1060, 83), (24, 45.9652, 86),
    (25, 93.6679, 71), (26, 67.8082, 70), (27, 55.9304, 68),
    (28, 46.0435, 72), (29, 84.9162, 69), (30, 104.9291, 65),
    (31, -40.1837, 66), (32, 134.3501, 63), (33, 719.3666, 91),
    (34, 758.4218, 90), (35, 700.7226, 88), (36, 756.0824, 92),
    (37, 441.8437, 89), (38, 108.4964, 94), (39, 286.9731, 50),
    (40, 251.9212, 96), (41, 361.7760, 95), (42, 193.7667, 100),
    (43, -102.7252, 101), (44, 1074.1000, 23), (45, 355.7381, 24),
    (46, 350.5184, 26), (47, 292.8046, 25), (48, 269.2471, 49),
    (49, 736.9540, 48), (50, 1216.0700, 27), (51, 255.8480, 4),
    (52, 252.9059, 53), (53, 123.2143, 73), (54, 127.3380, 74),
    (55, 222.2789, 51), (56, 20.1604, 75), (57, 226.1873, 55),
    (58, 86.4601, 62), (59, 224.2238, 59), (60, 134.9382, 58),
    (61, 34.2541, 57), (62, 97.4210, 60), (63, 206.6665, 61),
    (64, 128.0247, 56), (65, 48.8839, 93), (66, 375.0486, 98),
    (67, 126.3340, 97), (68, 375.8217, 20), (69, 238.2066, 21),
    (70, 2.8992, 80), (71, 9.3624, 79), (72, 603.5347, 13),
    (73, 662.0582, 9), (74, 510.9666, 22), (75, 1317.4360, 29),
    (76, 681.3525, 11), (77, 564.1116, 18), (78, 391.3697, 15),
    (79, 318.2350, 14), (80, 435.8446, 28), (81, 218.6012, 19),
    (82, 381.5945, 17), (83, 80.2735, 64), (84, 231.3919, 16),
    (85, 186.9204, 12), (86, 48.5026, 67), (87, 168.3007, 5),
    (88, 108.5277, 6), (89, 213.7165, 7), (90, 183.1130, 52),
    (91, 1178.1950, 102), (92, 158.3258, 54), (93, -47.3420, 78),
    (94, 186.7950, 31), (95, 392.2006, 8), (96, 595.1778, 10),
    (97, 268.7081, 99), (98, 183.3467, None), (99, 612.9546, 2),
    (100, 258.9924, 3),
    (101, -316.4392, 32), (102, 64.0566, 82), (103, 660.2247, 33),
    (104, 554.7941, 34), (105, 420.7591, 35), (107, -237.2124, 37),
    (110, 37.0590, None),
    (111, 319.4879, 41), (113, 118.8412, 43), (115, 305.1341, 45),
    (116, 191.5058, 46), (117, 423.5251, 47), (214, 86.9450, 116),
    (215, -66.5670, 77), (216, -81.1543, 76),
]

# Table 4 rows defined in Table 1 but carrying no dB value (paper
# Table 4 "do not estimate" markers / absent rows). 106, 108 and 118
# have definitions; 129 is defined as a second-order group but absent
# from Table 5 (the paper's Table 5 has no para row); 109 and the
# numbering gaps 112/114 have neither row nor usable definition.
PAPER_VALUELESS_DEFINED = {106, 108, 118, 129}
PAPER_ABSENT_IDS = {109, 112, 114}

# --- independent transcription: paper Table 5 (x 1e3) ------------------
PAPER_TABLE_5 = [
    (119, 34.3545), (120, 2.5030), (121, -83.3326), (122, -64.4854),
    (123, -125.9208), (124, -47.2962), (125, 33.9765), (126, -7.0982),
    (127, -45.0531), (128, -3.2036), (130, -20.6706), (131, -36.3170),
    (132, -1.1994), (133, 123.7433), (134, -15.9694), (217, 36.7574),
]

# --- independent transcription: paper Table 6 (x 1e3) ------------------
# (row id, class pair, value x 1e3). Row 218 is printed with a class
# pair but no value: refused, never zero.
PAPER_TABLE_6 = [
    (135, "A", "A", -561.5153), (136, "A", "M", 1067.6660),
    (137, "A", "N", 42.4825), (140, "A", "D", -799.5332),
    (141, "A", "E", -618.2760), (142, "A", "F", -1797.6930),
    (143, "A", "G", -1181.5990), (145, "A", "P", 1431.2430),
    (146, "A", "I", 1132.0400), (148, "B", "B", -97.6205),
    (151, "B", "D", -751.6676), (152, "B", "F", 548.4352),
    (157, "M", "M", 1085.8320), (158, "M", "N", -206.7811),
    (159, "D", "M", -198.2791), (160, "F", "M", -1676.4770),
    (162, "M", "Q", 1659.0340), (165, "N", "N", -307.1018),
    (166, "D", "N", 65.4421), (170, "L", "L", 240.3037),
    (172, "C", "C", -2601.7090), (175, "C", "G", -787.8563),
    (176, "O", "O", -3929.1300), (178, "D", "D", 144.6074),
    (179, "D", "E", 1118.9580), (180, "D", "F", -225.7802),
    (181, "D", "G", 1981.2980), (182, "D", "H", 362.7540),
    (183, "D", "J", -1425.0170), (185, "D", "P", 743.3353),
    (187, "E", "E", -3748.8180), (189, "F", "F", 920.3138),
    (190, "F", "G", 1594.1640), (192, "F", "P", 108.1305),
    (193, "F", "I", 1590.3210), (194, "G", "G", -1270.0830),
    (201, "H", "H", 946.7309), (204, "H", "I", 705.3049),
    (205, "J", "J", 838.3372), (206, "Q", "Q", -1501.3550),
    (208, "K", "R", 675.0414), (209, "P", "S", 994.4996),
    (210, "I", "R", 135.5896), (212, "S", "S", -29.6785),
    (218, "C", "M", None),
]

# --- paper-derivable anchors: single coefficients printed in the ------
# per-group rows of worked examples Tables 15a-15h (a source independent
# of the statistics tables above). (group id, dB x 1e3 as printed).
WORKED_EXAMPLE_ANCHORS = [
    # 15a alpha-pinene
    (1, 13.3063), (9, 22.2573), (10, 32.8162), (11, 4.8500),
    (62, 97.4210), (125, 33.9765), (132, -1.1994),
    # 15b 1,2-ethanediol
    (7, 37.8485), (36, 756.0824),
    # 15c perfluoro-2-propanone
    (21, 42.0887), (51, 255.8480), (120, 2.5030),
    (121, -83.3326), (123, -125.9208),
    # 15d acrylic acid
    (44, 1074.1000), (61, 34.2541), (134, -15.9694),
    # 15e glycol monoacetate
    (45, 355.7381),
    # 15f dipropyl succinate
    (4, 54.6564),
    # 15g diethanolamine
    (42, 193.7667),
    # 15h R122 (1,2,2-trichloro-1,1-difluoroethane)
    (26, 67.8082), (27, 55.9304), (124, -47.2962),
]

# interaction rows quoted inside the worked examples
WORKED_EXAMPLE_INTERACTIONS = [
    (135, -561.5153),  # 15b/15g (A,A) OH-OH
    (137, 42.4825),  # 15g (A,N) OH-amine
    (142, -1797.6930),  # 15e (A,F) OH-ester
    (189, 920.3138),  # 15f (F,F) ester-ester
]


@pytest.mark.parametrize(
    "gid,value,priority", PAPER_TABLE_4, ids=[str(r[0]) for r in PAPER_TABLE_4]
)
def test_table4_matches_independent_transcription(gid, value, priority):
    if priority is None:
        # valued row, but no Table 1 definition exists: parked, never
        # catalogued for counting
        assert gid not in NANNOOLAL_VP_GROUPS, gid
        assert NANNOOLAL_VP_UNDEFINED_VALUE_GROUPS[gid] == pytest.approx(
            value / 1000.0, abs=1e-9
        ), gid
        return
    name, db, pr = NANNOOLAL_VP_GROUPS[gid]
    assert db == pytest.approx(value / 1000.0, abs=1e-9), gid
    assert pr == priority, gid
    assert name


def test_table4_census():
    """Row counts and the valueless/absent/undefined rulings."""
    assert len(PAPER_TABLE_4) == 115
    assert len([r for r in PAPER_TABLE_4 if r[0] <= 117]) == 112
    transcribed = {r[0] for r in PAPER_TABLE_4 if r[2] is not None}
    valued = {g for g, r in NANNOOLAL_VP_GROUPS.items() if r[1] is not None}
    second_order = {r[0] for r in PAPER_TABLE_5}
    extensions = {301, 302, 303}
    # every valued catalogue entry is one of the three paper sources
    assert valued == transcribed | second_order | extensions
    # Table 4's missing-row census over the 1..117 ID range: 106/108
    # defined-valueless; 109 unusable definition + absent (adjudication
    # in the module docstring); 112/114 numbering gaps. 118 (defined
    # valueless, "do not estimate") sits outside the range and is
    # covered by the NO_PARAMETER loop below.
    missing = set(range(1, 118)) - {r[0] for r in PAPER_TABLE_4}
    assert missing == {106, 108} | PAPER_ABSENT_IDS
    valueless = {g for g, r in NANNOOLAL_VP_GROUPS.items() if r[1] is None}
    assert valueless == NANNOOLAL_VP_NO_PARAMETER
    # the four valueless groups stay catalogued (counted, then refused);
    # 109/112/114/98/110 are absent from the catalogue entirely
    for gid in PAPER_VALUELESS_DEFINED:
        assert gid in NANNOOLAL_VP_GROUPS and NANNOOLAL_VP_GROUPS[gid][1] is None
    # 129 (para pair) is defined but absent from Table 5: counted by the
    # fragmenter, then refused at estimation
    assert 129 not in second_order and 129 in NANNOOLAL_VP_NO_PARAMETER
    for gid in PAPER_ABSENT_IDS | {98, 110}:
        assert gid not in NANNOOLAL_VP_GROUPS, gid
    assert NANNOOLAL_VP_UNDEFINED_VALUE_GROUPS == {98: 0.1833467, 110: 0.0370590}


def test_priorities_unique():
    """Table 1's priority numbers are unique among the main groups."""
    prs = [r[2] for r in NANNOOLAL_VP_GROUPS.values() if r[2] is not None]
    assert len(prs) == len(set(prs))
    # the minimum is the valueless 1,2-diketone (118, PR 1); the maximum
    # is CH C(c) (6, PR 121)
    assert min(prs) == 1 and max(prs) == 121
    # second-order corrections and extensions carry no priority number
    second_order = {r[0] for r in PAPER_TABLE_5} | {301, 302, 303}
    for gid, record in NANNOOLAL_VP_GROUPS.items():
        if gid in second_order:
            assert record[2] is None, gid


@pytest.mark.parametrize(
    "gid,value", PAPER_TABLE_5, ids=[str(r[0]) for r in PAPER_TABLE_5]
)
def test_table5_matches_independent_transcription(gid, value):
    assert len(PAPER_TABLE_5) == 16
    assert NANNOOLAL_VP_GROUPS[gid][1] == pytest.approx(value / 1000.0, abs=1e-9), gid


@pytest.mark.parametrize(
    "row,ci,cj,value", PAPER_TABLE_6, ids=[str(r[0]) for r in PAPER_TABLE_6]
)
def test_table6_matches_independent_transcription(row, ci, cj, value):
    assert len(PAPER_TABLE_6) == 45
    assert NANNOOLAL_VP_INTERACTION_IDS[(ci, cj)] == row
    if value is None:
        # printed without a parameter: refused, never zeroed
        assert (ci, cj) not in NANNOOLAL_VP_INTERACTIONS, row
    else:
        assert NANNOOLAL_VP_INTERACTIONS[(ci, cj)] == pytest.approx(
            value / 1000.0, abs=1e-9
        ), row
    # interactions are symmetric; the catalogue stores each pair once
    assert (cj, ci) not in NANNOOLAL_VP_INTERACTIONS or ci == cj


@pytest.mark.parametrize(
    "gid,value", WORKED_EXAMPLE_ANCHORS,
    ids=[str(r[0]) for r in WORKED_EXAMPLE_ANCHORS]
)
def test_worked_example_group_anchors(gid, value):
    """Third source: the per-group rows of Tables 15a-15h must agree
    with the Table 4/5 statistics transcriptions."""
    assert NANNOOLAL_VP_GROUPS[gid][1] == pytest.approx(value / 1000.0, abs=1e-9), gid


@pytest.mark.parametrize(
    "row,value", WORKED_EXAMPLE_INTERACTIONS,
    ids=[str(r[0]) for r in WORKED_EXAMPLE_INTERACTIONS]
)
def test_worked_example_interaction_anchors(row, value):
    by_row = {v: k for k, v in NANNOOLAL_VP_INTERACTION_IDS.items()}
    assert NANNOOLAL_VP_INTERACTIONS[by_row[row]] == pytest.approx(
        value / 1000.0, abs=1e-9
    ), row


def test_base_constants():
    """Paper Eq. (6) slope constant 4.1012 and the 1/8 shift; Eq. (7)
    intercept -0.176055."""
    assert NANNOOLAL_VP_BASE == {
        "log10_term": 4.1012,
        "trb_shift": 0.125,
        "constant": -0.176055,
    }


def test_questionable_values_transcribed_verbatim():
    """Both values the paper flags with a footnote-a question mark are
    carried verbatim, so a future correction stays visible."""
    assert NANNOOLAL_VP_GROUPS[91][1] == pytest.approx(1.1781950, abs=1e-9)
    assert NANNOOLAL_VP_INTERACTIONS[("E", "E")] == pytest.approx(-3.7488180, abs=1e-9)


def test_interaction_classes_match_part1():
    """Part 3 reuses Part 1's 19 interaction classes (letters A-S) with
    identical memberships - including Part 3's own OH/OH(a) split of
    Part 1's 37-row double listing (P3 Table 3: OH = 34,35,36)."""
    p1 = NANNOOLAL_BP_INTERACTION_CLASSES
    assert set(NANNOOLAL_VP_INTERACTION_CLASSES) == set(p1)
    for letter, (_name, members) in NANNOOLAL_VP_INTERACTION_CLASSES.items():
        assert members == p1[letter][1], letter
        for gid in members:
            assert gid in NANNOOLAL_VP_GROUPS, (letter, gid)


def test_table6_row_ids_match_part1():
    """The 44 valued Table 6 rows share Part 1's row numbering; row 218
    (COOH-NH2, printed "Do not estimate") exists only in Part 3."""
    assert len(NANNOOLAL_VP_INTERACTIONS) == 44
    for pair in NANNOOLAL_VP_INTERACTIONS:
        assert NANNOOLAL_VP_INTERACTION_IDS[pair] == NANNOOLAL_BP_INTERACTION_IDS[
            pair
        ], pair
    assert set(NANNOOLAL_VP_INTERACTION_IDS) - set(NANNOOLAL_VP_INTERACTIONS) == {
        ("C", "M")
    }
    assert NANNOOLAL_VP_INTERACTION_IDS[("C", "M")] == 218


def test_extension_values():
    """Compernolle et al. (2010) Table 3 dB column (x 1e3): 888.42,
    926.02, 519.03. Extension groups take the no-GI branch; their dB
    values ride in the main catalogue."""
    assert {gid: NANNOOLAL_VP_GROUPS[gid][1] for gid in NANNOOLAL_VP_EXTENSIONS} == {
        301: 0.88842,
        302: 0.92602,
        303: 0.51903,
    }
    for gid, (_name, provenance) in NANNOOLAL_VP_EXTENSIONS.items():
        assert "Compernolle" in provenance, gid
