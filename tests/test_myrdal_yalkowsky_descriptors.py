# SPDX-License-Identifier: MIT
"""Verification gates 4 and 6 for the Myrdal-Yalkowsky descriptor layer.

Gate 6: the counting conventions of
    :func:`aeroprops.vapour_pressure.myrdal_yalkowsky_descriptors`,
    each pinned to its paper evidence - the Eq. 8 word definitions
    (nonterminal, nonring, tert-butyl/triple-bond exclusions, the
    paper's propane/butane/pentane/1-pentene examples), the ring-system
    (connected component) reading of RING against the SSSR reading, the
    sqrt-root HBN group catalogue of Eq. 9 with the Compernolle (2010)
    hydroperoxide/peracid extension, and the refusal scope.
Gate 4: a 50+-molecule battery - element bookkeeping (formula against
    hand-written truth, molecular weight against RDKit's), ring-system
    counts against hand-checked molecular graphs, tau against hand
    counts - plus a cross-method sanity check of the descriptor MW
    against the structure layer.
"""

import pytest

from rdkit.Chem import Descriptors

from aeroprops.structure import element_counts, mol_from_smiles
from aeroprops.vapour_pressure import myrdal_yalkowsky_descriptors


def _desc(smi):
    descriptors, success, status = myrdal_yalkowsky_descriptors(
        mol_from_smiles(smi)
    )
    assert success, status
    return descriptors


# ---------------------------------------------------------------------------
# Gate 6: counting conventions
# ---------------------------------------------------------------------------

@pytest.mark.parametrize(
    "smi, sp3, sp2, ring, tau",
    [
        # the paper's own Eq. 8 examples
        ("CCC", 1, 0, 0, 0.0),          # propane
        ("CCCC", 2, 0, 0, 1.0),         # butane: "2 SP3 atoms - 1"
        ("CCCCC", 3, 0, 0, 2.0),        # pentane
        ("C=CCCC", 2, 1, 0, 1.5),       # 1-pentene (paper: 1.5)
        # terminal atoms are not counted
        ("CO", 0, 0, 0, 0.0),           # methanol carbon: degree 1
        ("CC=O", 0, 1, 0, 0.0),         # acetaldehyde: floored
        # heteroatoms of degree >= 2 count (paper: "NH, N, O, and S")
        ("CCOCC", 3, 0, 0, 2.0),        # diethyl ether: 2 CH2 + O
        ("CSC", 1, 0, 0, 0.0),          # dimethyl sulfide: S only
        ("CCNc1ccccc1", 2, 0, 1, 1.5),  # n-ethylaniline: CH2 + NH
        ("CCCCCCS", 5, 0, 0, 4.0),      # 1-hexanethiol: terminal -SH not
        # triple bonds are not counted
        ("C#CCCCC", 3, 0, 0, 2.0),      # hex-1-yne: sp carbons skipped
        # tert-butyl centres excluded; plain quaternary C kept
        ("CC(C)(C)CC", 1, 0, 0, 0.0),   # 2,2-dimethylbutane: CH2 only
        ("CC(C)(C)C(C)CC", 2, 0, 0, 1.0),   # 2,2,4-trimethylpentane
        ("CC(C)CC(C)(C)CC", 4, 0, 0, 3.0),  # 2,4,4-trimethylhexane
        ("C(C)(C)(C)C", 0, 0, 0, 0.0),  # neopentane
        # ring atoms never count; substituent CH/CH2 do
        ("CC(C)C1CCCCC1", 1, 0, 1, 0.5),    # isopropylcyclohexane
        ("CCC1CCCCC1", 1, 0, 1, 0.5),       # ethylcyclohexane
        ("C1COCCO1", 0, 0, 1, 0.0),         # 1,4-dioxane: ring O excluded
        # sp2 carbons: chain alkenes and every carbonyl carbon
        ("CC=CCC", 1, 2, 0, 1.0),       # pent-2-ene
        ("CC(=O)C", 0, 1, 0, 0.0),      # acetone
        ("CCOC(=O)CC", 3, 1, 0, 2.5),   # ethyl acetate: O + 2 CH2 + C=O
        ("CCCCC(=O)O", 3, 1, 0, 2.5),   # valeric acid (Table 1)
        ("COC(=O)c1ccccc1O", 1, 1, 1, 1.0),  # methyl salicylate (Table 1)
    ],
    ids=[
        "propane", "butane", "pentane", "1-pentene",
        "methanol-terminal", "acetaldehyde-floor",
        "ether-oxygen", "sulfide-sulfur", "ethylaniline-nh", "hexanethiol",
        "hexyne-triple",
        "tbu-22dmb", "isooctane", "244tmh", "neopentane",
        "ipr-cyclohexane", "ethyl-cyclohexane", "dioxane",
        "pent2ene", "acetone", "ethyl-acetate", "valeric", "methyl-salicylate",
    ],
)
def test_eq8_conventions(smi, sp3, sp2, ring, tau):
    d = _desc(smi)
    assert d["sp3"] == sp3
    assert d["sp2"] == sp2
    assert d["ring_systems"] == ring
    assert d["tau"] == tau


@pytest.mark.parametrize(
    "smi, systems, sssr",
    [
        # RING = independent ring systems (connected components of the
        # ring atoms through ring bonds) - NOT the SSSR ring count.
        # Fused and spiro systems are one system; linked rings are two.
        ("c1ccccc1", 1, 1),                    # benzene
        ("c1ccc2ccccc2c1", 1, 2),              # naphthalene (fused)
        ("c1ccc2cc3ccccc3cc2c1", 1, 3),        # anthracene (fused)
        ("c1ccc2c(c1)ccc3ccccc23", 1, 3),      # phenanthrene (fused)
        ("c1ccccc1c1ccccc1", 2, 2),            # biphenyl (linked)
        ("c1ccc(-c2ccccc2-c2ccccc2)cc1", 3, 3),    # p-terphenyl
        ("C1CCC2(C1)CCCCC2", 1, 2),            # spiro[4.5]decane
        ("C1C2CC3CC1CC(C2)C3", 1, 3),          # adamantane
        ("C1CC2CCC1C2", 1, 2),                 # norbornane
        ("C1CCCCC1C1CCCCC1", 2, 2),            # bicyclohexyl
        ("CC1(C2CCC(=C)C1C2)C", 1, 2),         # beta-pinene (Table 1)
    ],
    ids=[
        "benzene", "naphthalene", "anthracene", "phenanthrene", "biphenyl",
        "terphenyl", "spiro", "adamantane", "norbornane", "bicyclohexyl",
        "beta-pinene",
    ],
)
def test_ring_systems_vs_sssr(smi, systems, sssr):
    """The two RING readings agree on simple rings and diverge exactly
    on fused systems (anthracene 1 vs 3 - the reason the paper's
    "independent ring systems" wording is implemented, not UMS's SSSR)."""
    from rdkit import Chem

    mol = mol_from_smiles(smi)
    d = _desc(smi)
    assert d["ring_systems"] == systems
    # the contrast quantity is the plain SSSR ring count (UMS's reading);
    # GetSymmSSSR would inflate adamantane to 4 by re-adding its
    # symmetry-equivalent fourth chair ring
    assert len(Chem.GetSSSR(mol)) == sssr


@pytest.mark.parametrize(
    "smi, expected",
    [
        # Eq. 9 group catalogue. Table 1 anchors first.
        ("Oc1cccc(O)c1", {"hydroxyl": 2}),          # resorcinol
        ("CCc1cc(O)cc(CC)c1", {"hydroxyl": 1}),     # 3,5-diethylphenol
        ("Oc1ccc2ccccc2c1", {"hydroxyl": 1}),       # beta-naphthol (phenol)
        ("CCCCC(=O)O", {"carboxyl": 1, "hydroxyl": 0}),
        ("OC(=O)CCCCC(=O)O", {"carboxyl": 2}),      # adipic acid
        # amide nitrogens are not primary amines
        ("NC(=O)C", {"primary_amine": 0}),
        ("NC(=O)c1ccccc1", {"primary_amine": 0}),   # benzamide
        ("Nc1ccccc1", {"primary_amine": 1}),        # aniline
        ("CCN", {"primary_amine": 1}),              # ethylamine
        ("CCNCC", {"primary_amine": 0}),            # diethylamine (2deg)
        # Compernolle (2010) extension: -OOH and peracids in the root
        ("CC(C)(C)OO", {"hydroperoxide": 1, "peracid": 0, "hydroxyl": 0}),
        ("COO", {"hydroperoxide": 1, "peracid": 0, "hydroxyl": 0}),
        ("CC(=O)OO", {"peracid": 1, "hydroperoxide": 0, "carboxyl": 0}),
        ("OC(=O)CC(=O)OO", {"carboxyl": 1, "peracid": 1}),
        # no double counting of peracid as hydroperoxide or acid; a
        # di-peracid counts two peracid groups
        ("OOC(=O)CCCCC(=O)OO", {"peracid": 2, "hydroperoxide": 0,
                                "carboxyl": 0}),
        ("CC(=O)OOC(C)=O", {"peracid": 0, "hydroperoxide": 0}),  # no OOH H
        # hydroxylamine (O on N) is neither alcohol nor peroxide
        ("CON", {"hydroxyl": 0, "hydroperoxide": 0}),
    ],
    ids=[
        "resorcinol", "diethylphenol", "naphthol", "valeric", "adipic",
        "acetamide", "benzamide", "aniline", "ethylamine", "diethylamine",
        "tbuoooh", "meooh", "peracetic", "peroxyadipic", "diperoxyacid",
        "diacetyl-peroxide", "hydroxylamine",
    ],
)
def test_eq9_hbond_catalogue(smi, expected):
    d = _desc(smi)
    for key, value in expected.items():
        assert d[key] == value, key
    root = d["hydroxyl"] + d["carboxyl"] + d["hydroperoxide"] + d["peracid"]
    assert d["hbond_root"] == root


def test_hbn_uses_all_element_molecular_weight():
    """Unlike UMS (whose composition mass omits S), the MW dividing HBN
    covers every element: 1-hexanethiol HBN is zero here, but a
    sulfuric OH carrier shows the convention through the g/mol scale."""
    d = _desc("CCCCCCS")
    assert d["hbn"] == 0.0
    d2 = _desc("OCCS")  # 2-hydroxyethanethiol: OH on MW 78.136 g/mol
    assert d2["hbn"] == pytest.approx(1.0 / d2["mw"])
    assert d2["mw"] == pytest.approx(78.136, abs=5e-4)


@pytest.mark.parametrize(
    "smi, group",
    [
        ("CCO[N+](=O)[O-]", "organic nitrate"),
        ("CC(C)O[N+](=O)[O-]", "organic nitrate"),
        ("CC(=O)OON(=O)=O", "peroxy acyl nitrate"),
        ("OC=C", "enol"),
        ("CC=CO", "enol"),
    ],
    ids=["ethyl-nitrate", "nitrate-ester", "pan-notation", "vinyl-alcohol",
         "enol-methyl"],
)
def test_refused_functional_groups(smi, group):
    """Nitrates, PAN and enols are outside the published basis; the
    descriptor layer refuses them (SIMPOL.1/EVAPORATION territory)."""
    descriptors, success, status = myrdal_yalkowsky_descriptors(
        mol_from_smiles(smi)
    )
    assert not success
    assert descriptors == {}
    assert status == "unsupported functional group: " + group


@pytest.mark.parametrize("smi", ["[Si](C)(C)(C)C", "OP(=O)(O)O"])
def test_refused_elements(smi):
    descriptors, success, status = myrdal_yalkowsky_descriptors(
        mol_from_smiles(smi)
    )
    assert not success
    assert status.startswith("unsupported element(s):")


def test_failed_mol():
    descriptors, success, status = myrdal_yalkowsky_descriptors(None)
    assert (descriptors, success) == ({}, False)
    assert status == "Failed to construct mol"


# ---------------------------------------------------------------------------
# Gate 4: 50+-molecule descriptor battery
# ---------------------------------------------------------------------------

#: (name, SMILES, formula, ring systems, tau) - hand-checked molecular
#: graphs and paper-convention counts across the method's domain.
BATTERY = [
    ("propane", "CCC", "C3H8", 0, 0.0),
    ("butane", "CCCC", "C4H10", 0, 1.0),
    ("pentane", "CCCCC", "C5H12", 0, 2.0),
    ("hexane", "CCCCCC", "C6H14", 0, 3.0),
    ("octacosane", "C" * 28, "C28H58", 0, 25.0),
    ("neopentane", "C(C)(C)(C)C", "C5H12", 0, 0.0),
    ("2,2-dimethylbutane", "CC(C)(C)CC", "C6H14", 0, 0.0),
    ("2,2,4-trimethylpentane", "CC(C)(C)C(C)CC", "C8H18", 0, 1.0),
    ("2,4,4-trimethylhexane", "CC(C)CC(C)(C)CC", "C9H20", 0, 3.0),
    ("1-pentene", "C=CCCC", "C5H10", 0, 1.5),
    ("pent-2-ene", "CC=CCC", "C5H10", 0, 1.0),
    ("1-undecene", "CCCCCCCCCC=C", "C11H22", 0, 7.5),
    ("styrene", "C=Cc1ccccc1", "C8H8", 1, 0.0),
    ("hex-1-yne", "C#CCCCC", "C6H10", 0, 2.0),
    ("cyclohexane", "C1CCCCC1", "C6H12", 1, 0.0),
    ("methylcyclohexane", "CC1CCCCC1", "C7H14", 1, 0.0),
    ("ethylcyclohexane", "CCC1CCCCC1", "C8H16", 1, 0.5),
    ("isopropylcyclohexane", "CC(C)C1CCCCC1", "C9H18", 1, 0.5),
    ("cis-decalin", "C1CCC2CCCCC2C1", "C10H18", 1, 0.0),
    ("cyclopropane", "C1CC1", "C3H6", 1, 0.0),
    ("norbornane", "C1CC2CCC1C2", "C7H12", 1, 0.0),
    ("adamantane", "C1C2CC3CC1CC(C2)C3", "C10H16", 1, 0.0),
    ("spiro[4.5]decane", "C1CCC2(C1)CCCCC2", "C10H18", 1, 0.0),
    ("bicyclohexyl", "C1CCCCC1C1CCCCC1", "C12H22", 2, 0.0),
    ("beta-pinene", "CC1(C2CCC(=C)C1C2)C", "C10H16", 1, 0.0),
    ("benzene", "c1ccccc1", "C6H6", 1, 0.0),
    ("toluene", "Cc1ccccc1", "C7H8", 1, 0.0),
    ("naphthalene", "c1ccc2ccccc2c1", "C10H8", 1, 0.0),
    ("anthracene", "c1ccc2cc3ccccc3cc2c1", "C14H10", 1, 0.0),
    ("phenanthrene", "c1ccc2c(c1)ccc3ccccc23", "C14H10", 1, 0.0),
    ("biphenyl", "c1ccccc1c1ccccc1", "C12H10", 2, 0.0),
    ("p-terphenyl", "c1ccc(-c2ccccc2-c2ccccc2)cc1", "C18H14", 3, 0.5),
    ("resorcinol", "Oc1cccc(O)c1", "C6H6O2", 1, 0.0),
    ("3,5-diethylphenol", "CCc1cc(O)cc(CC)c1", "C10H14O", 1, 1.5),
    ("beta-naphthol", "Oc1ccc2ccccc2c1", "C10H8O", 1, 0.0),
    ("methanol", "CO", "CH4O", 0, 0.0),
    ("ethanol", "CCO", "C2H6O", 0, 0.0),
    ("1-nonanol", "CCCCCCCCCO", "C9H20O", 0, 7.0),
    ("1-tetradecanol", "CCCCCCCCCCCCCCO", "C14H30O", 0, 12.0),
    ("glycerol", "OCC(O)CO", "C3H8O3", 0, 2.0),
    ("2-methoxyethanol", "COCCO", "C3H8O2", 0, 2.0),
    ("diethyl ether", "CCOCC", "C4H10O", 0, 2.0),
    ("1,4-dioxane", "C1COCCO1", "C4H8O2", 1, 0.0),
    ("acetone", "CC(=O)C", "C3H6O", 0, 0.0),
    ("ethyl acetate", "CCOC(=O)CC", "C5H10O2", 0, 2.5),
    ("methyl salicylate", "COC(=O)c1ccccc1O", "C8H8O3", 1, 1.0),
    ("valeric acid", "CCCCC(=O)O", "C5H10O2", 0, 2.5),
    ("adipic acid", "OC(=O)CCCCC(=O)O", "C6H10O4", 0, 4.0),
    ("tBu hydroperoxide", "CC(C)(C)OO", "C4H10O2", 0, 0.0),
    ("peracetic acid", "CC(=O)OO", "C2H4O3", 0, 0.5),
    ("methylamine", "CN", "CH5N", 0, 0.0),
    ("ethylamine", "CCN", "C2H7N", 0, 0.0),
    ("aniline", "Nc1ccccc1", "C6H7N", 1, 0.0),
    ("n-ethylaniline", "CCNc1ccccc1", "C8H11N", 1, 1.5),
    ("2-methyl-5-ethylpyridine", "CCc1ccc(C)nc1", "C8H11N", 1, 0.5),
    ("morpholine", "C1COCCN1", "C4H9NO", 1, 0.0),
    ("N-methylpiperidine", "CN1CCCCC1", "C6H13N", 1, 0.0),
    ("formamide", "NC=O", "CH3NO", 0, 0.0),
    ("dimethyl sulfide", "CSC", "C2H6S", 0, 0.0),
    ("1-hexanethiol", "CCCCCCS", "C6H14S", 0, 4.0),
    ("dipropyl disulfide", "CCCSSCCC", "C6H14S2", 0, 5.0),
    ("thiophene", "c1ccsc1", "C4H4S", 1, 0.0),
    ("chlorobenzene", "Clc1ccccc1", "C6H5Cl", 1, 0.0),
    ("1-chlorotetradecane", "CCCCCCCCCCCCCCCl", "C14H29Cl", 0, 12.0),
    ("bromoethane", "CCBr", "C2H5Br", 0, 0.0),
    ("iodomethane", "CI", "CH3I", 0, 0.0),
]


def _formula(mol):
    counts = element_counts(mol)
    parts = []
    if "C" in counts:
        parts.append(("C", counts["C"]))
        if "H" in counts:
            parts.append(("H", counts["H"]))
    for symbol in sorted(k for k in counts if k not in ("C", "H")):
        parts.append((symbol, counts[symbol]))
    return "".join(
        symbol + (str(n) if n > 1 else "") for symbol, n in parts
    )


@pytest.mark.parametrize(
    "name,smi,formula,ring,tau",
    BATTERY,
    ids=[b[0] for b in BATTERY],
)
def test_battery_descriptor_sanity(name, smi, formula, ring, tau):
    """Element bookkeeping (hand-written formulas), ring-system counts
    (hand-checked molecular graphs), the tau walk, and MW wiring
    against RDKit's own molecular weight."""
    mol = mol_from_smiles(smi)
    d = _desc(smi)
    assert _formula(mol) == formula
    assert d["ring_systems"] == ring
    assert d["tau"] == tau
    assert d["mw"] == pytest.approx(Descriptors.MolWt(mol))  # rel 1e-6
    # the descriptor dictionary is closed: HBN is the group counts and
    # MW combined by the pure layer, nothing else enters
    root = d["hydroxyl"] + d["carboxyl"] + d["hydroperoxide"] + d["peracid"]
    from math import sqrt

    assert d["hbn"] == pytest.approx(
        (sqrt(root) + 0.33 * sqrt(d["primary_amine"])) / d["mw"]
    )


def test_battery_size():
    """Gate 4 is a 50+-molecule battery by design."""
    assert len(BATTERY) >= 50
