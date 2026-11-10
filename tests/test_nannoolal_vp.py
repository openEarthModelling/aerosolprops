# SPDX-License-Identifier: MIT
"""Verification gates 2/3/4/6 for the Nannoolal Part 3 fragmentation.

Gate 2 - all eight worked examples of paper Tables 15a-15h reproduce
end-to-end through the real fragmenter: group counts, atom number, dB
(including the 15d/15g misprint adjudications) and vapour pressure /
the 15f Appendix-A boiling-point inversion.
Gate 3 - every estimate sits within a stated tolerance of the
experimental value the paper quotes alongside it.
Gate 4 - a broad battery fragments with a full atom balance.
Gate 6 - priority adjudication cases (18/214 fused/bridge split, 118
1,2-diketone, 134 conjugation, halogen families, ester split,
valueless groups counted then refused) and read-only fragmentation.
"""

import pytest

from aeroprops.methods.nannoolal_vp import (
    LOG10_ATM_TO_PA,
    nannoolal_vp_boiling_point,
    nannoolal_vp_dB,
    nannoolal_vp_log10_vapour_pressure,
)
from aeroprops.structure import mol_from_smiles
from aeroprops.vapour_pressure import (
    _nannoolal_vp_engine,
    nannoolal_vp_fragmentation,
)


def _frag(smi):
    return nannoolal_vp_fragmentation(mol_from_smiles(smi))


def _kPa(smi, T, Tb):
    counts, n, ok, status = _frag(smi)
    assert ok, status
    log10_pa = nannoolal_vp_log10_vapour_pressure(counts, n, T, Tb)
    return 10.0 ** (log10_pa - LOG10_ATM_TO_PA) * 101.325


# --- Gate 2: paper Tables 15a-15h worked examples -----------------------
# (name, SMILES, expected counts, n, dB, T [K], Tb,exp [K], P [kPa])
# 15g pins the formula-derived dB 1.6131919; the printed 1.7493490 is
# the adjudicated misprint (module docstring, note OQ-3).
WORKED = [
    ("15a alpha-pinene", "CC1=CCC2CC1C2(C)C",
     {1: 3, 9: 2, 10: 2, 11: 1, 62: 1, 125: 1, 132: 2}, 10,
     0.1078596, 388.15, 429.00, 31.03),
    ("15b ethanediol", "OCCO",
     {7: 2, 36: 2}, 4,
     1.1310491, 410.65, 470.50, 13.05),
    ("15c hexafluoroacetone", "FC(F)(F)C(=O)C(F)(F)F",
     {7: 2, 21: 6, 51: 1, 120: 1, 121: 2, 123: 1}, 10,
     0.1179392, 210.16, 245.90, 14.63),
    ("15d acrylic acid", "C=CC(=O)O",
     {44: 1, 61: 1, 134: 1}, 5,
     0.9163297, 344.15, 413.60, 6.52),
    ("15e glycol monoacetate", "CC(=O)OCCO",
     {1: 1, 7: 2, 36: 1, 45: 1}, 7,
     0.5111422, 352.65, 458.65, 2.24),
    ("15g diethanolamine", "N(CCO)CCO",
     {7: 4, 36: 2, 42: 1}, 7,
     1.6131919, 401.13, 541.15, 0.404),
    # 15h's printed "number of atoms: 5" contradicts its own formula
    # (CHCl2CClF2 = 7 heavy atoms) and every other table of the series;
    # n plays no role in 15h anyway (no interaction classes, GI = 0).
    ("15h R122", "C(Cl)(Cl)C(F)(F)Cl",
     {7: 2, 21: 2, 26: 2, 27: 1, 121: 1, 124: 1}, 7,
     0.0447374, 297.46, 344.25, 17.5093),
]


@pytest.mark.parametrize(
    "name,smi,expected,n,db,T,Tb,p_kpa",
    WORKED,
    ids=[w[0].split(" ")[0] for w in WORKED],
)
def test_worked_examples_end_to_end(name, smi, expected, n, db, T, Tb, p_kpa):
    counts, got_n, ok, status = _frag(smi)
    assert (ok, status) == (True, "OK"), name
    assert counts == expected, name
    assert got_n == n, name
    got_db, db_status = nannoolal_vp_dB(counts, n)
    assert got_db == pytest.approx(db, abs=5e-8), name
    assert db_status == "OK", name
    assert _kPa(smi, T, Tb) == pytest.approx(p_kpa, abs=5e-3), name


def test_15f_boiling_point_inversion():
    """Table 15f: dipropyl succinate, Tb from P = 0.4 kPa at 374.65 K.

    The printed inversion formula is garbled (~210 K read literally);
    the Appendix-A form reproduces the printed 520.85 K."""
    counts, n, ok, status = _frag("CCCOC(=O)CCC(=O)OCCC")
    assert (ok, counts, n) == (True, {1: 2, 4: 4, 7: 2, 45: 2}, 14)
    db, _ = nannoolal_vp_dB(counts, n)
    assert db == pytest.approx(0.9878298, abs=5e-8)
    tb = nannoolal_vp_boiling_point(counts, n, 374.65, 0.4e3)
    assert tb == pytest.approx(520.85, abs=0.01)


# --- Gate 3: estimates vs the experimental values the paper quotes ------
# (name, SMILES, T, Tb,exp, P_exp kPa, relative tolerance). 15e's own
# deviation (2.24 vs 1.6 kPa) is the paper's, not ours.
LITERATURE = [
    ("15a", "CC1=CCC2CC1C2(C)C", 388.15, 429.00, 31.68, 0.05),
    ("15b", "OCCO", 410.65, 470.50, 12.27, 0.07),
    ("15c", "FC(F)(F)C(=O)C(F)(F)F", 210.16, 245.90, 14.05, 0.05),
    ("15d", "C=CC(=O)O", 344.15, 413.60, 6.67, 0.05),
    ("15e", "CC(=O)OCCO", 352.65, 458.65, 1.60, 0.45),
    ("15g", "N(CCO)CCO", 401.13, 541.15, 0.410, 0.05),
    ("15h", "C(Cl)(Cl)C(F)(F)Cl", 297.46, 344.25, 16.804, 0.05),
]


@pytest.mark.parametrize(
    "name,smi,T,Tb,p_exp,rtol",
    LITERATURE,
    ids=[r[0] for r in LITERATURE],
)
def test_estimates_vs_experimental(name, smi, T, Tb, p_exp, rtol):
    assert _kPa(smi, T, Tb) == pytest.approx(p_exp, rel=rtol), name


def test_15f_vs_experimental_boiling_point():
    # Mustafaev & Ragimov report Tb = 523.95 K; the estimate is 0.6% off.
    counts, n, ok, _ = _frag("CCCOC(=O)CCC(=O)OCCC")
    assert ok
    tb = nannoolal_vp_boiling_point(counts, n, 374.65, 0.4e3)
    assert tb == pytest.approx(523.95, rel=0.01)


# --- Gate 4: atom balance over a broad battery --------------------------
# The Part 1 battery plus the Part 3 deltas (diketone, bridge carbons,
# conjugation). Molecules whose only obstacle is a valueless group
# (118/129/106/108) fragment fine and are refused later - they live in
# BATTERY_VALUELESS.
BATTERY = [
    "CC", "CCC", "CCCC", "CCCCCCCCCC", "CC(C)C", "CC(C)(C)C",
    "C1CCCCC1", "C1CCCCCC1", "C1CC1", "C1CCC1",
    "C=CC", "C=CCC", "CC=CC", "C=CC=C", "C=C(C)C=C", "C#CC", "CC#CC",
    "C=C=C", "C1=CCCCC1", "c1ccccc1", "Cc1ccccc1", "CCc1ccccc1",
    "CC(C)c1ccccc1", "CC(C)(C)c1ccccc1", "c1ccc(-c2ccccc2)cc1",
    "c1ccc2ccccc2c1", "Cc1ccccc1C", "Cc1cccc(C)c1",
    # imidazoles (c1ncc[nH]1, Cn1ccnc1) live in BATTERY_VALUELESS:
    # Part 3 prints no dB for group 106 (Part 1 did). Pyrazole
    # (c1cn[nH]c1) stays deliberately absent: its NH bonds to a ring N
    # and no group covers it (unclaimed atom, refused).
    "c1ccc2[nH]ccc2c1", "n1ccccc1", "n1cccc2ccccc12", "c1cc[nH]c1",
    "c1ccsc1", "c1nccs1", "c1ccoc1",
    "CC(=O)Oc1ccccc1C(=O)O",
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
    "FC(F)(F)C(=O)C(F)(F)F", "O=Cc1ccoc1", "O=Cc1ccccc1",
    "c1cc(-c2ccccc2)cc(-c3ccccc3)c1",
]

BATTERY_VALUELESS = [
    "c1ncc[nH]1",  # 106 imidazole
    "Cn1ccnc1",  # 106 N-methylimidazole
    "CC(=O)C(C)=O",  # 118 1,2-diketone
    "CC(=O)C(=O)C",  # 118 + leftover ketone 51
    "Cc1ccc(C)cc1",  # 129 para pair
    "c1cnc[nH]1",  # 106 imidazole
    "CSC#N",  # 108 thiocyanate
]


def _assert_atom_balance(smi):
    mol = mol_from_smiles(smi)
    counts, n, ok, status = nannoolal_vp_fragmentation(mol)
    assert ok, "%s: %s" % (smi, status)
    assert counts and n == mol.GetNumAtoms()
    _c, assignments, _ok, _s = _nannoolal_vp_engine().fragment(mol)
    claimed = sum(
        len(claim) for insts in assignments.values() for claim in insts
    )
    assert claimed == n, smi
    return counts, n


def test_battery_fragments_with_full_atom_balance():
    assert len(BATTERY) >= 90
    for smi in BATTERY:
        counts, n = _assert_atom_balance(smi)
        db, status = nannoolal_vp_dB(counts, n)
        assert db is not None, "%s: %s" % (smi, status)


def test_valueless_battery_fragments_then_refuses():
    """Valueless groups (no dB row) are counted, then refused with a
    status naming the group - never silently zeroed."""
    expected_gid = {118: None, 129: None, 106: None, 108: None}
    for smi in BATTERY_VALUELESS:
        counts, n = _assert_atom_balance(smi)
        db, status = nannoolal_vp_dB(counts, n)
        assert db is None, smi
        valueless = [g for g in counts if g in expected_gid]
        assert valueless, smi
        assert status.startswith(
            "no dB parameter for group %d " % valueless[0]
        ), smi


# --- Gate 6: priority adjudication on the real catalogue ----------------
@pytest.mark.parametrize(
    "smi,expected",
    [
        # Part 3's 18/214 split: three aromatic bonds = fused (18,
        # naphthalene); two aromatic bonds + bridge bond = 214 (biphenyl,
        # m-terphenyl - the paper's own examples).
        ("c1ccc2ccccc2c1", {15: 8, 18: 2}),
        ("c1ccccc1-c1ccccc1", {15: 10, 214: 2}),
        ("c1cc(-c2ccccc2)cc(-c3ccccc3)c1", {15: 14, 214: 4, 128: 1}),
        # 118 (PR 1) claims the 1,2-diketone before any carbonyl family;
        # a third carbonyl of a triketone falls back to the ketone 51.
        ("CC(=O)C(C)=O", {1: 2, 118: 1}),
        ("CC(=O)C(=O)C(=O)C", {1: 2, 118: 1, 51: 1, 134: 1}),
        # 134 conjugation (Part 3's renumbering of the Part 1 118): the
        # paper's own examples benzaldehyde and furfural.
        ("O=Cc1ccccc1", {15: 5, 16: 1, 90: 1, 134: 1}),
        ("O=Cc1ccoc1", {15: 3, 16: 1, 65: 1, 90: 1, 134: 1}),
        # Halogen families (verified against Tables 15c/15h): 21 (F, H0
        # carbon, co-F) beats 102 (co-Cl) by priority; 26/27 split by
        # co-halogen count; 22/23 split by co-halogen count on CH.
        # CF3Cl: zero hydrogens -> 123 (not 124)
        ("FC(F)(F)Cl", {7: 1, 21: 3, 27: 1, 121: 1, 123: 1}),
        ("FC(Cl)C(F)F", {7: 2, 22: 3, 26: 1}),
        ("ClCCl", {7: 1, 26: 2}),
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
    counts, _n, ok, status = _frag(smi)
    assert (ok, status) == (True, "OK")
    assert counts == expected, smi


def test_fragmentation_is_read_only():
    """The fragmenter never modifies the input molecule."""
    from aeroprops.structure import load_rdkit

    to_smiles = load_rdkit()[0].MolToSmiles
    for smi in BATTERY[:12] + BATTERY_VALUELESS:
        mol = mol_from_smiles(smi)
        before = to_smiles(mol)
        nannoolal_vp_fragmentation(mol)
        assert to_smiles(mol) == before, smi


def test_out_of_scope_molecules_fail_cleanly():
    counts, n, ok, status = _frag("C[Si](C)(C)C")
    assert not ok and "unsupported element" in status
    assert (counts, n) == ({}, 5)
    # pyridine N-oxide: no group parameterises the exocyclic O(-);
    # the partial counts of the claimed atoms come back with the
    # failure (same contract as the Part 1 layer).
    counts, n, ok, status = _frag("c1cc[n+]([O-])cc1")
    assert not ok and "unclaimed atom" in status
    assert counts == {15: 5, 67: 1} and n == 7
