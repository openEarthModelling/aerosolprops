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


# --- Gate 2: paper Tables 17a-17d worked examples ---------------------
WORKED_EXAMPLES = [
    # (name, SMILES, expected counts, n, paper-calculated Tb [K])
    ("17a 3,3,4,4-tetramethylhexane", "CCC(C)(C)C(C)(C)CC",
     {1: 6, 4: 2, 6: 2, 133: 1}, 10, 429.5),
    ("17b di-isopropanolamine", "CC(O)CNCC(O)C",
     {1: 2, 7: 4, 34: 2, 42: 1}, 9, 509.3),
    ("17c perfluoro-2-propanone", "FC(F)(F)C(=O)C(F)(F)F",
     {7: 2, 21: 6, 51: 1, 120: 1, 121: 2, 123: 1}, 10, 246.3),
    # 17d: the arithmetic gives 490.65 -> 490.7; the paper prints
    # 490.6 (display truncation of the same group sum 2655.112).
    ("17d methyl m-toluate", "COC(=O)c1cccc(C)c1",
     {2: 1, 3: 1, 15: 4, 16: 2, 45: 1, 118: 1, 128: 1}, 11, 490.7),
]


@pytest.mark.parametrize(
    "name,smi,expected,n,tb_expected", WORKED_EXAMPLES)
def test_worked_examples_end_to_end(name, smi, expected, n, tb_expected):
    nnp = NannoolalBP(smi=smi)
    assert (nnp.success, nnp.status, nnp.n) == (True, "OK", n)
    assert nnp.counts == expected
    assert nnp.boiling_point() == pytest.approx(tb_expected, abs=0.06)


# --- Gate 3: literature values ----------------------------------------
LITERATURE = [
    # (name, SMILES, experimental Tb [K], source)
    ("furfural", "O=Cc1ccoc1", 434.7, "CRC Handbook, 97th ed."),
    ("anisole", "COc1ccccc1", 426.9, "CRC Handbook, 97th ed."),
    ("morpholine", "C1COCCN1", 401.8, "CRC Handbook, 97th ed."),
    ("acetophenone", "CC(=O)c1ccccc1", 475.6, "CRC Handbook, 97th ed."),
    ("quinoline", "n1cccc2ccccc12", 511.2, "CRC Handbook, 97th ed."),
]


@pytest.mark.parametrize("name,smi,exp,_src", LITERATURE)
def test_literature_boiling_points(name, smi, exp, _src):
    # The paper reports the method as a regression with a few-K spread;
    # these five sit within 1.2 K of experiment (they also pin the
    # o/m/p, ring and conjugation machinery to real numbers).
    assert NannoolalBP(smi=smi).boiling_point() == pytest.approx(
        exp, abs=1.5
    )


# --- Gate 4: atom balance over a broad battery ------------------------
BATTERY = [
    "CC", "CCC", "CCCC", "CCCCCCCCCC", "CC(C)C", "CC(C)(C)C",
    "C1CCCCC1", "C1CCCCCC1", "C1CC1", "C1CCC1",
    "C=CC", "C=CCC", "CC=CC", "C=CC=C", "C=C(C)C=C", "C#CC", "CC#CC",
    "C=C=C", "C1=CCCCC1", "c1ccccc1", "Cc1ccccc1", "CCc1ccccc1",
    "CC(C)c1ccccc1", "CC(C)(C)c1ccccc1", "c1ccc(-c2ccccc2)cc1",
    "c1ccc2ccccc2c1", "Cc1ccccc1C", "Cc1cccc(C)c1", "Cc1ccc(C)cc1",
    "c1ccc2[nH]ccc2c1", "n1ccccc1", "n1cccc2ccccc12", "c1cc[nH]c1",
    "c1ncc[nH]1", "Cn1ccnc1", "c1ccsc1", "c1nccs1", "c1ccoc1",
    "CC(=O)Oc1ccccc1C(=O)O",  # aspirin: ortho pair + ester + acid + ether
    "OCC", "CCO", "CCCO", "CCCCO", "CCCCCO", "CC(C)(C)O", "OCCO",
    "OCC(O)CO", "Oc1ccccc1", "CCOc1ccccc1", "OC(=O)c1ccccc1",
    "CCOC(C)=O", "O=C1OCCC1", "CC(=O)OC(C)=O", "O=C1OC(=O)C=C1",
    "CC(=O)C", "CCC(=O)C", "O=Cc1ccccc1", "CC(=O)c1ccccc1",
    "O=C(c1ccccc1)c1ccccc1", "CCOCC", "C1COCCO1", "C1COC1", "C1C(O1)CCl",
    "CN(C)C=O", "CC(N)=O", "CN1CCCC1=O", "NC(=O)N", "COC(=O)NC",
    "COC(=O)OC", "C1COC(=O)O1", "Nc1ccccc1", "CCN(CC)CC", "C1COCCN1",
    "C[N+](C)(C)C", "CC#N", "N#Cc1ccccc1", "C[N+](=O)[O-]",
    "O=[N+]([O-])c1ccccc1", "CCO[N+](=O)[O-]", "CON=O", "CC(=NO)C",
    "CN=C=O", "CSCC", "CSSCC", "CS(C)=O", "CS(C)(=O)=O", "C1CS(=O)(=O)CC1",
    "ClCCl", "ClC(Cl)Cl", "ClC(Cl)(Cl)Cl", "ClCCCl", "C=CCl",
    "ClCc1ccccc1", "Clc1ccccc1", "FC(F)(F)C(F)(F)F", "FC(F)(F)Cl",
    "FC(F)=CF", "CC(C)(C)OO", "CC(C)(C)OOCC(C)(C)", "CC(=O)OO",
    "CC(=O)OON(=O)=O", "CCCC(=O)OON(=O)=O",
]


def test_battery_fragments_with_full_atom_balance():
    assert len(BATTERY) >= 50
    engine = _nannoolal_engine()
    for smi in BATTERY:
        mol = mol_from_smiles(smi)
        counts, n, ok, status = nannoolal_fragmentation(mol)
        assert ok, "%s: %s" % (smi, status)
        assert counts and n == mol.GetNumAtoms()
        # The engine's exclusive claim set must cover exactly n atoms.
        _c, assignments, _ok, _s = engine.fragment(mol)
        claimed = sum(
            len(claim) for insts in assignments.values() for claim in insts
        )
        assert claimed == n
        assert NannoolalBP(mol=mol).boiling_point() is not None


# --- Gate 6: priority adjudication on the real catalogue --------------


@pytest.mark.parametrize(
    "smi,expected",
    [
        # 22 (CH carbon, exactly one co-halogen) vs 23 (two co-halogens):
        # the CHF2 end of FC(Cl)C(F)F carries one co-fluorine only.
        ("FC(Cl)C(F)F", {7: 2, 22: 3, 26: 1}),
        ("FC(F)Cl", {7: 1, 23: 2, 27: 1, 121: 1, 124: 1}),
        # 21 beats 19 (quaternary carbon with co-fluorine).
        ("FC(F)(F)CCl", {7: 2, 21: 3, 25: 1, 121: 1}),
        # 26/27 split by co-halogen count, never plain 25.
        ("ClCCl", {7: 1, 26: 2}),
        ("ClC(Cl)Cl", {7: 1, 27: 3, 121: 1, 124: 1}),
        ("ClCCCl", {7: 2, 25: 2}),
        # 45 acyclic vs 47 cyclic ester split on the ring O.
        ("CCOC(C)=O", {1: 2, 7: 1, 45: 1}),
        ("O=C1OCCC1", {47: 1, 9: 2, 13: 1, 126: 1}),
        # 39 claims the epoxide triangle before the plain ether 38.
        ("C1C(O1)CCl", {7: 1, 25: 1, 39: 1, 125: 1}),
        # 61 claims a terminal C=C before the substituted families.
        ("C=CCl", {29: 1, 61: 1}),
        ("CC=CC", {1: 2, 58: 1}),
    ],
)
def test_priority_adjudication(smi, expected):
    counts, _n, ok, status = nannoolal_fragmentation(mol_from_smiles(smi))
    assert (ok, status) == (True, "OK")
    assert counts == expected


# --- Gate 7: no-GI branch and refusals --------------------------------


@pytest.mark.parametrize(
    "smi",
    [
        "CC(C)(C)OO",  # 301 hydroperoxide
        "CC(=O)OO",  # 302 peroxy acid
        "CC(=O)OON(=O)=O",  # 303 PAN
        "CCCC(=O)OON(=O)=O",  # 303 PBN
    ],
)
def test_extension_groups_take_the_no_gi_branch(smi):
    nnp = NannoolalBP(smi=smi)
    assert nnp.success
    assert any(gid in nnp.counts for gid in (301, 302, 303))
    # GI is skipped (zero) and flagged, not silently estimated.
    gi, status = nannoolal_bp_gi(nnp.counts, nnp.n)
    assert (gi, status) == (0.0, GI_STATUS_NO_GI)
    assert nnp.boiling_point() is not None
    assert nnp.gi_status == GI_STATUS_NO_GI


@pytest.mark.parametrize(
    "smi,pair",
    [
        # Table 5 row 153: B-G printed without a value.
        ("CC(=O)c1ccc(O)cc1", ("B", "G")),
        # Table 5 row 211: R-R printed without a value.
        ("c1cn[nH]c1", ("R", "R")),
    ],
)
def test_missing_interaction_parameters_refuse(smi, pair):
    nnp = NannoolalBP(smi=smi)
    assert nnp.success
    gi, status = nannoolal_bp_gi(nnp.counts, nnp.n)
    assert gi is None
    assert status == "missing interaction parameter (%s, %s)" % pair
    assert nnp.boiling_point() is None
    assert nnp.status == "refused: missing interaction parameter"


def test_plain_molecules_get_gi_status_ok():
    nnp = NannoolalBP(smi="CC(O)CNCC(O)C")
    assert nnp.boiling_point() is not None
    assert nnp.gi_status == GI_STATUS_OK


def test_out_of_scope_molecules_fail_cleanly():
    # Element gate: the paper covers Si/P but aeroprops does not
    # parameterise those families (documented scope cut).
    nnp = NannoolalBP(smi="C[Si](C)(C)C")
    assert not nnp.success
    assert "unsupported element" in nnp.status
    assert nnp.boiling_point() is None
    # No group parameterises pyridine N-oxide: unclaimed atom, and the
    # wrapper refuses rather than estimating something meaningless.
    nnp = NannoolalBP(smi="c1cc[n+]([O-])cc1")
    assert not nnp.success
    assert "unclaimed atom" in nnp.status
