# SPDX-License-Identifier: MIT
"""EVAPORATION verification gates 1, 2, 3, 4, 6 and 7.

Gate 1: coefficient catalogue transcription + g_k cross-check.
Gate 2: the paper's worked examples (D-values), reproduced through the
        full SMILES -> counts -> estimator pipeline.
Gate 3: literature boiling points (independent of the paper's fit).
Gate 4: corpus sanity (element conservation, structural bounds).
Gate 6: counting-semantics units for both conventions.
Gate 7: mode selection and failure paths.

Gate 5 (live UManSysProp comparison) lives in test_oracle_evaporation.py.
"""

from math import log10

import pytest

from aeroprops.methods.evaporation import (
    EVAPORATION_B20_UPSTREAM_TYPO,
    EVAPORATION_COEFFICIENTS,
    EVAPORATION_GK_PAPER,
    EVAPORATION_KAPPA,
    evaporation_boiling_temperature,
    evaporation_diacid_factor,
    evaporation_gk,
    evaporation_log10_vapour_pressure,
)
from aeroprops.structure import mol_from_smiles
from aeroprops.vapour_pressure import EVAPORATION, evaporation_fragmentation

LOG10_ATM = log10(101325.0)


def paper_value(smi, T=298.15):
    """log10(p0/atm) of a SMILES under paper counting."""
    counts, ok, status = evaporation_fragmentation(mol_from_smiles(smi))
    assert ok, status
    return evaporation_log10_vapour_pressure(counts, T) - LOG10_ATM


def paper_counts(smi):
    """Rounded counts (exact dict equality over fractional c3)."""
    counts, ok, status = evaporation_fragmentation(mol_from_smiles(smi))
    assert ok, status
    return {k: round(v, 9) for k, v in counts.items()}


# ----------------------------------------------------------------------
# Gate 1: coefficients
# ----------------------------------------------------------------------
class TestCoefficients:
    def test_spot_transcription(self):
        # A few exact spot checks against paper Table 4 (k, A, B).
        assert EVAPORATION_COEFFICIENTS["1"][:2] == (2.6255, -1986.56)
        assert EVAPORATION_COEFFICIENTS["9"][:2] == (0.98567, -23671.0)
        assert EVAPORATION_COEFFICIENTS["12"][:2] == (0.18704, -2509.37)
        assert EVAPORATION_COEFFICIENTS["20"][:2] == (0.46023, 1817.69)
        assert len(EVAPORATION_COEFFICIENTS) == 20

    def test_b20_upstream_typo(self):
        # The upstream source carries 1817.60 where the paper prints
        # 1817.69; upstream_compat assembly uses the typo.
        assert EVAPORATION_B20_UPSTREAM_TYPO == 1817.60

    def test_gk_cross_check(self):
        # Group-contribution value g_k(T) computed the long way (Eq. 5
        # with the pure descriptors of one monofunctional molecule)
        # vs the paper's Table 2. k5 and k9 are off by 0.006/0.008
        # because the paper's own Table 2 rounding is inconsistent
        # with its Table 4 coefficients; 0.009 covers both.
        for k in range(2, 21):
            assert evaporation_gk(k) == pytest.approx(
                EVAPORATION_GK_PAPER[str(k)], abs=0.009), k

    def test_kappa(self):
        assert EVAPORATION_KAPPA == 1.5

    def test_diacid_factor(self):
        # Eq. 32: the factor needs c9 > 1 AND N_CL + N_HB > 2, so two
        # acids alone stay at 1.
        assert evaporation_diacid_factor({"9": 3}) == 2.6 / 3.0
        assert evaporation_diacid_factor({"9": 2}) == 1.0
        assert evaporation_diacid_factor({"9": 1, "8": 5}) == 1.0
        assert evaporation_diacid_factor({"9": 2, "8": 1}) == 2.6 / 3.0
        assert evaporation_diacid_factor({"9": 2, "8": 1}) == 2.6 / 3.0


# ----------------------------------------------------------------------
# Gate 2: paper worked examples
# ----------------------------------------------------------------------
class TestPaperAnchors:
    def test_d2_neopentane_branching(self):
        # Paper Eq. 26 discussion: t = 2 for neopentane; the branching
        # descriptor alone is worth +2 x 0.199 log units over n-pentane.
        neo = paper_value("CC(C)(C)C")
        pent = paper_value("CCCCC")
        assert paper_counts("CC(C)(C)C") == {"2": 5, "3": 2.0}
        assert neo - pent == pytest.approx(2 * 0.199, abs=2e-3)

    def test_d3_dicarbonyl_interactions(self):
        # alpha/beta interactions are counted per carbonyl (glyoxal
        # k16 = 2, not the pair-deduplicated 1 of the notes draft).
        assert paper_counts("O=CC=O")["16"] == 2
        assert paper_counts("CC(=O)C=O")["16"] == 2
        assert paper_counts("COC(=O)C(=O)C")["16"] == 1  # methyl pyruvate
        assert paper_counts("CC(=O)C(=O)O").get("16", 0) == 0  # acid
        assert paper_counts("CC(=O)CC(=O)C")["17"] == 2  # 2,4-PD

    def test_d5_mbtca(self):
        assert paper_counts("CC(C)(C(=O)O)C(C(=O)O)CC(=O)O") == \
            {"2": 8, "3": 3.0, "9": 3}
        assert paper_value("CC(C)(C(=O)O)C(C(=O)O)CC(=O)O") == \
            pytest.approx(-11.5184, abs=5e-5)

    def test_d6_nitrate_and_pan(self):
        assert paper_value("CO[N+](=O)[O-]") == \
            pytest.approx(-0.6114, abs=5e-5)
        # The notes' D6 PAN row (-2.2247) double-counts the ester
        # group; the correct counts {"2": 2, "7": 1} give -1.3561.
        assert paper_counts("CC(=O)OO[N+](=O)[O-]") == {"2": 2, "7": 1}
        assert paper_value("CC(=O)OO[N+](=O)[O-]") == \
            pytest.approx(-1.3561, abs=5e-5)

    def test_d7_ring_branching(self):
        # Fractional ring number in t: alpha-pinene 4 - 6/4, fenchyl
        # alcohol 5 - 7/5.
        assert paper_counts("CC1=CCC2CC1C2(C)C") == {"2": 10, "3": 2.5}
        assert paper_counts("C1(C)(C)C(O)C2(C)CCC1C2") == \
            {"2": 10, "3": 3.6, "8": 1, "12HB": 1, "14": 1}
        assert paper_value("CC1=CCC2CC1C2(C)C") == \
            pytest.approx(-2.1131, abs=5e-5)

    def test_d8_functional_panel(self):
        cases = [
            ("O=C(O)CCC(=O)O", -7.5026),        # succinic acid
            ("OC(=O)C(O)C(O)C(=O)O", -7.4042),  # tartaric acid, c19 = 4
            ("CC(=O)CO", -2.1556),              # hydroxyacetone
            ("CC(O)C(O)C", -3.3489),            # 2,3-butanediol
            ("CC(=O)OO", -1.4258),              # peracetic acid
            ("CCOCCO", -2.3602),                # 2-ethoxyethanol
            ("COC(C)(C)C", -0.4716),            # MTBE (see note below)
        ]
        for smi, want in cases:
            assert paper_value(smi) == pytest.approx(want, abs=5e-5), smi

    def test_mtbe_note(self):
        # MTBE c2 = 6 (5 C + the in-chain ether O of Eq. 26) and the
        # resulting -0.4716 matches the notes' D8 row after all; an
        # earlier working note had pinned -0.9567 by omitting the O,
        # which contradicts the paper's own 2-ethoxyethanol entry.
        assert paper_counts("COC(C)(C)C") == {"2": 6, "3": 1.0}

    def test_pinonaldehyde_pinonic(self):
        # Paper Table 9 molecules. Pinonaldehyde reproduces exactly.
        assert paper_counts("O=CCC1CC(C)(C)C1C(C)=O") == \
            {"2": 10, "3": 3.0, "5": 2}
        assert paper_value("O=CCC1CC(C)(C)C1C(C)=O") == \
            pytest.approx(-4.5172, abs=5e-5)
        # Pinonic acid (same skeleton, COOH): the chemically correct
        # counts give -6.8198. The -7.3048 seen in working notes
        # corresponds to c2 = 11 (one extra) and is not reproducible.
        assert paper_counts("O=C(O)CC1CC(C)(C)C1C(C)=O") == \
            {"2": 10, "3": 3.0, "5": 1, "9": 1}
        assert paper_value("O=C(O)CC1CC(C)(C)C1C(C)=O") == \
            pytest.approx(-6.8198, abs=5e-5)

    def test_footnote_f_beta_interaction(self):
        # Footnote f: a beta partner does not count when the carbon in
        # between carries a functional group. Upstream has no footnote.
        smi = "CC(=O)C(O)C(=O)C"  # 3-hydroxy-2,4-pentanedione
        assert paper_counts(smi).get("17", 0) == 0
        counts, ok, _ = evaporation_fragmentation(
            mol_from_smiles(smi), upstream_compat=True)
        assert ok and counts["17"] == 2


# ----------------------------------------------------------------------
# Gate 3: independent literature values (boiling points, Eq. 14 basis)
# ----------------------------------------------------------------------
class TestLiterature:
    # Experimental 1 atm boiling points [K] (CRC Handbook values);
    # the paper itself reports 30-50 % vapour-pressure uncertainty for
    # multifunctional species, so 2 % on Tb is a fair independent gate.
    BOILING_POINTS = [
        ("CCCCCC", 341.9),   # n-hexane
        ("CC(=O)C", 329.2),  # acetone
        ("CCO", 351.5),      # ethanol
        ("CC(C)(C)C", 282.7),  # neopentane
        ("CCCCO", 390.0),    # 1-butanol
    ]

    @pytest.mark.parametrize("smi,tb", BOILING_POINTS)
    def test_boiling_point(self, smi, tb):
        est = EVAPORATION(smi=smi).boiling_point()
        assert est == pytest.approx(tb, rel=0.02), smi

    def test_adipic_acid_pressure(self):
        # Sub-cooled liquid log10 p0 of adipic acid at 298 K is about
        # -3.5 Pa in the diacid literature (e.g. Soonsin et al., 2010;
        # Bilde et al., 2003); the paper-mode estimate lands at -3.47.
        assert paper_value("O=C(O)CCCCC(=O)O", T=298.15) + LOG10_ATM == \
            pytest.approx(-3.5, abs=0.1)


# ----------------------------------------------------------------------
# Gate 4: corpus sanity
# ----------------------------------------------------------------------
CORPUS = [
    "C", "C(C)(C)(C)C", "C1(C)(C)C(O)C2(C)CCC1C2", "C1CCCCC1", "C1CCOC1",
    "C1COCCO1", "C1COOCC1", "C=CC(=O)C", "C=CC=O", "C=CCO", "C=CO", "C=O",
    "CC(=O)C(=O)CC", "CC(=O)C(=O)O", "CC(=O)C(O)C(=O)C", "CC(=O)CC(=O)C",
    "CC(=O)CO", "CC(=O)OCC", "CC(=O)OCCO", "CC(=O)OO",
    "CC(=O)OO[N+](=O)[O-]", "CC(C)(C(=O)O)C(C(=O)O)CC(=O)O", "CC(C)(C)OO",
    "CC(O)(C)C", "CC(O)C(O)C", "CC1=CCC2CC1C2(C)C", "CCC(C)CC", "CCCC",
    "CCOC(=O)C(=O)C", "CCOC(C)(C)C", "CCOCCO", "COC(=O)C(=O)C", "COC=O",
    "COCC=O", "CO[N+](=O)[O-]", "O=C(O)C(=O)O", "O=C(O)CC(=O)O",
    "O=C(O)CCCCC(=O)O", "O=C1CCCCC1", "O=C1OCCC1", "O=CC(=O)O", "O=CC=O",
    "OC(=O)C(O)C(O)C(=O)O", "OC(=O)CC(O)(CC(=O)O)C(=O)O", "OC(C)(C)OC",
    "OCC(=O)O", "OCCO", "OCCOO", "OCCO[N+](=O)[O-]",
]


class TestCorpusSanity:
    @pytest.mark.parametrize("smi", CORPUS)
    def test_both_modes_succeed_and_finite(self, smi):
        for upstream in (False, True):
            counts, ok, status = evaporation_fragmentation(
                mol_from_smiles(smi), upstream_compat=upstream)
            assert ok, (smi, status)
            value = evaporation_log10_vapour_pressure(
                counts, 298.15, upstream)
            assert value == value and abs(value) < 100.0, (smi, value)

    @pytest.mark.parametrize("smi", CORPUS)
    def test_descriptor_bounds(self, smi):
        mol = mol_from_smiles(smi)
        n_c = sum(1 for a in mol.GetAtoms() if a.GetSymbol() == "C")
        n_o = sum(1 for a in mol.GetAtoms() if a.GetSymbol() == "O")
        for upstream in (False, True):
            counts, ok, _ = evaporation_fragmentation(
                mol, upstream_compat=upstream)
            assert ok
            # c2 = carbons + counted oxygens, so it bounds below by #C
            # and above by #C + #O
            assert n_c <= counts["2"] <= n_c + n_o, (smi, counts["2"])
            # oxygen-based descriptors cannot exceed #O
            for k in ("5", "6", "7", "8", "9", "10", "11", "13"):
                assert counts.get(k, 0) <= n_o, (smi, k)

    def test_molecule_not_modified(self):
        # fragmentation is read-only on the input molecule
        from aeroprops.structure import load_rdkit
        to_smiles = load_rdkit()[0].MolToSmiles
        for smi in CORPUS[:12]:
            mol = mol_from_smiles(smi)
            before = to_smiles(mol)
            evaporation_fragmentation(mol)
            assert to_smiles(mol) == before, smi


# ----------------------------------------------------------------------
# Gate 6: counting-semantics units (every rule decoded from the live
# oracle in upstream mode, or from the paper anchors in paper mode)
# ----------------------------------------------------------------------
def upstream_counts(smi):
    counts, ok, status = evaporation_fragmentation(
        mol_from_smiles(smi), upstream_compat=True)
    assert ok, status
    return {k: round(v, 9) for k, v in counts.items()}


class TestCountingSemantics:
    def test_in_chain_oxygen_paper(self):
        # Eq. 26: ethers/esters/peroxides break the C skeleton on
        # removal and count; acid OH, peroxy and nitrate O do not.
        assert paper_counts("CCOCC") == {"2": 4 + 1}
        assert paper_counts("CC(C)(C)OOC(C)(C)C")["2"] == 8 + 2
        assert paper_counts("CO[N+](=O)[O-]") == {"2": 1, "4": 1}

    def test_lactone_k12_split(self):
        # Paper: a lactone counts twice (12lin bridge O + 12CL acyl C);
        # upstream's 12CL aggregates only 1202+502, so the ring-ester
        # acyl carbon (602) lands in "12" alone.
        assert paper_counts("O=C1OCCC1") == \
            {"2": 5, "3": -1.0, "6": 1, "12lin": 1, "12CL": 1}
        assert upstream_counts("O=C1OCCC1") == \
            {"2": 5, "3": -1.0, "6": 1, "12": 1}

    def test_ring_peroxide(self):
        # Ring O-O pair: upstream fires 1202 (12lin + 12CL) *and* the
        # peroxide-bridge c2 submode (203); paper takes 12lin only.
        assert upstream_counts("C1COOCC1") == \
            {"2": 6, "3": -1.0, "12": 1, "12lin": 1, "12CL": 1}
        assert paper_counts("C1COOCC1") == \
            {"2": 6, "3": -1.0, "12lin": 1}

    def test_formate(self):
        # The formate bridge O is a full ester in the paper (c6, c2)
        # but invisible upstream (acyl-C-with-H submodes never match).
        assert paper_counts("COC=O") == {"2": 3, "6": 1}
        assert upstream_counts("COC=O") == {"2": 2}

    def test_branching_binomial_model(self):
        # Upstream c3 via subset-matching: neopentane's quaternary
        # centre fires m301 = C(4,4) = 1 and m302 = C(4,3) = 4, so
        # c3 = 2 + 4; DTBP's two peroxide-bearing centres fire m303
        # (w2) each: c3 = 2 x 2.
        assert upstream_counts("CC(C)(C)C")["3"] == 6.0
        assert upstream_counts("CC(C)(C)OOC(C)(C)C") == \
            {"2": 10, "3": 4.0}

    def test_hydroxyl_classification(self):
        # Upstream primary/secondary/tertiary by heavy-neighbour degree
        # of the carbinol carbon (2-methoxy-2-propanol is 803 with
        # only two C neighbours); paper k14 by carbon neighbours.
        assert upstream_counts("OC(C)(C)OC") == \
            {"2": 5, "3": 1.0, "8": 1, "14": 2}
        # paper: the two methyls are the only C neighbours (ncc=2, no
        # c3) and carry no oxygen (no c19); k14 still fires
        assert paper_counts("OC(C)(C)OC") == \
            {"2": 5, "8": 1, "14": 1}
        # acid OH belongs to the acid, not k8, in both modes
        assert upstream_counts("O=C(O)CCCCC(=O)O") == {"2": 6, "9": 2}
        assert paper_counts("O=C(O)CCCCC(=O)O") == {"2": 6, "9": 2}
        # enol OH counts upstream (804), not in the paper mode; an
        # aromatic OH matches nothing in either mode
        assert upstream_counts("C=C(O)") == {"2": 2, "8": 1}
        assert paper_counts("C=C(O)") == {"2": 2}
        assert upstream_counts("c1ccccc1O") == {"2": 6, "3": -1.0}

    def test_vicinal_hydroxyl_pairs(self):
        # 1902 counts pairs (assembly doubles them): ethylene glycol
        # c19 = 2, glycerol has two pairs so c19 = 4.
        assert upstream_counts("OCCO")["19"] == 2
        assert upstream_counts("OCC(O)CO") == \
            {"2": 3, "8": 3, "14": 1, "19": 4}

    def test_dicarbonyl_interaction_submodes(self):
        # alpha: ketald pairs involving an aldehyde deduped (1601),
        # per-ketone partners (1602), ester partners (1603); beta:
        # strictly 1,3-dicarbonyl -- 1,4 pairs (succinaldehyde,
        # 4-oxopentanal) fire nothing.
        assert upstream_counts("CC(=O)C=O")["16"] == 2
        assert upstream_counts("CC(=O)C(=O)CC")["16"] == 2
        assert upstream_counts("CCOC(=O)C=O")["16"] == 1
        assert upstream_counts("O=CCC=O")["17"] == 1   # malonaldehyde
        assert upstream_counts("CCOC(=O)CC(=O)C")["17"] == 1  # 1703
        assert upstream_counts("O=CCCC=O").get("17", 0) == 0
        assert upstream_counts("CC(=O)CCC=O").get("17", 0) == 0
        # pyruvic acid: acid is not an alpha partner (k16 = 0) but the
        # acid counts its ketone neighbour (k20, submode 2002)
        assert upstream_counts("CC(=O)C(=O)O") == \
            {"2": 3, "5": 1, "9": 1, "20": 1}

    def test_nitrate_pan_dead_upstream(self):
        # The upstream nitrate/PAN SMARTS never match under RDKit
        # charge normalisation: c4 = c7 = 0 regardless of input form.
        assert upstream_counts("CO[N+](=O)[O-]") == {"2": 1}
        assert upstream_counts("CC(=O)OO[N+](=O)[O-]") == {"2": 2}

    def test_alpha_oxygen_partners(self):
        # k18 takes the OOH oxygen on the alpha carbon (1801); an
        # oxygen on a carbonyl alpha-carbon (methyl pyruvate) fires
        # nothing; paper k18 counts OH/OOH/nitrate but not ethers.
        assert upstream_counts("CC(=O)COO") == \
            {"2": 3, "5": 1, "10": 1, "18": 1}
        assert upstream_counts("COC(=O)C(=O)C").get("18", 0) == 0
        assert paper_counts("CC(=O)COCC").get("18", 0) == 0  # ether
        assert paper_counts("CC(=O)COO")["18"] == 1


# ----------------------------------------------------------------------
# Gate 7: mode selection and failure paths
# ----------------------------------------------------------------------
class TestModesAndFailures:
    def test_modes_agree_where_semantics_coincide(self):
        # adipic acid has no branching/ring/interaction differences:
        # both conventions yield identical counts and values.
        paper = EVAPORATION(smi="O=C(O)CCCCC(=O)O")
        up = EVAPORATION(smi="O=C(O)CCCCC(=O)O", upstream_compat=True)
        assert paper.counts == {"2": 6, "9": 2}
        assert paper.log10_vapour_pressure() == \
            pytest.approx(up.log10_vapour_pressure())

    def test_unsupported_element(self):
        for smi in ("CCCl", "CS", "CCBr"):
            eva = EVAPORATION(smi=smi)
            assert not eva.success
            assert eva.counts == {}
            assert eva.vapour_pressure() is None
            assert eva.log10_vapour_pressure() is None
            assert eva.boiling_point() is None

    def test_paper_mode_refuses_aromatics(self):
        eva = EVAPORATION(smi="c1ccccc1")
        assert not eva.success
        assert "aromatic" in eva.status
        # upstream mode accepts aromatics (SMARTS semantics there)
        assert EVAPORATION(smi="c1ccccc1O", upstream_compat=True).success

    def test_paper_mode_refuses_non_nitrate_nitrogen(self):
        assert not EVAPORATION(smi="CCN").success
        assert EVAPORATION(smi="CO[N+](=O)[O-]").success
        assert EVAPORATION(smi="CC(=O)OO[N+](=O)[O-]").success

    def test_invalid_smiles(self):
        eva = EVAPORATION(smi="this is not a molecule")
        assert not eva.success
        assert eva.counts == {}

    def test_non_positive_temperature(self):
        eva = EVAPORATION(smi="CC(C)(C)C")
        assert eva.log10_vapour_pressure(0.0) is None
        assert eva.log10_vapour_pressure(-10.0) is None
        assert eva.vapour_pressure(0.0) is None

    def test_boiling_point_no_solution(self):
        # A and B of the same sign have no positive Eq. 14 solution
        # (b14 > 0 flips the sign of the otherwise-negative B).
        assert evaporation_boiling_temperature({"14": 1}) is None

    def test_direct_mol_input(self):
        eva = EVAPORATION(mol=mol_from_smiles("CC(C)(C)C"))
        assert eva.success
        assert eva.counts == {"2": 5, "3": 2.0}
