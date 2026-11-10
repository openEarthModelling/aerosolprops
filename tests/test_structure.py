# SPDX-License-Identifier: MIT
"""Tests for aeroprops.structure (RDKit entry layer)."""

import pytest

from aeroprops.structure import (
    canonical_smiles,
    element_counts,
    load_rdkit,
    mol_from_smiles,
    molecular_weight,
)


def test_load_rdkit_returns_modules():
    assert load_rdkit() is not None
    Chem, Descriptors = load_rdkit()
    assert hasattr(Chem, "MolFromSmiles")
    assert hasattr(Descriptors, "MolWt")


def test_nitro_valence_forms_agree():
    # The neutral valence form and the charge form must map to the same
    # canonical molecule (the form aeroprops SMARTS catalogues target).
    assert canonical_smiles("CN(=O)=O") == canonical_smiles("C[N+](=O)[O-]")


def test_mol_from_smiles_invalid_returns_none():
    assert mol_from_smiles("not a smiles") is None


def test_element_counts_includes_implicit_hydrogens():
    counts = element_counts(mol_from_smiles("CCO"))
    assert counts == {"C": 2, "H": 6, "O": 1}
    counts = element_counts(mol_from_smiles("c1ccccc1"))
    assert counts == {"C": 6, "H": 6}


def test_molecular_weight_ethanol_si():
    # 46.07 g/mol expressed as kg/mol
    assert molecular_weight(mol_from_smiles("CCO")) == pytest.approx(0.046068, abs=1e-5)


def test_none_molecule_returns_none():
    assert element_counts(None) is None
    assert molecular_weight(None) is None
