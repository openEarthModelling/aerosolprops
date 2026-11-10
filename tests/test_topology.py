# SPDX-License-Identifier: MIT
"""Tests for aerosolprops.topology (graph primitives over RDKit built-ins)."""

from aerosolprops.structure import mol_from_smiles
from aerosolprops.topology import fragments_excluding_bond, ring_arc_distance, separate_without_atom


def test_fragments_excluding_bond_propane():
    propane = mol_from_smiles("CCC")
    assert fragments_excluding_bond(propane, 0, 1) == ((0,), (1, 2))
    # deleting a terminal bond of a chain leaves the rest connected
    assert fragments_excluding_bond(propane, 1, 2) == ((0, 1), (2,))


def test_separate_without_atom_ether_vs_ethanol():
    ether = mol_from_smiles("CCOCC")
    # original atom indices; bridging oxygen is index 2
    assert separate_without_atom(ether, 2) == ((0, 1), (3, 4))
    # terminal oxygen removal leaves a single fragment
    assert separate_without_atom(mol_from_smiles("CCO"), 2) == ((0, 1),)


def test_separate_without_atom_indices_are_original():
    # atoms above the deleted index keep their original numbering
    butane = mol_from_smiles("CCCC")
    assert separate_without_atom(butane, 1) == ((0,), (2, 3))
    assert separate_without_atom(butane, 2) == ((0, 1), (3,))


def test_ring_arc_distance_benzene():
    benzene = mol_from_smiles("c1ccccc1")
    assert [ring_arc_distance(benzene, 0, j) for j in range(6)] == [0, 1, 2, 3, 2, 1]


def test_ring_arc_distance_unshared_returns_none():
    ethylbenzene = mol_from_smiles("CCc1ccccc1")
    assert ring_arc_distance(ethylbenzene, 0, 4) is None
    # two ring atoms always share the benzene ring
    assert ring_arc_distance(ethylbenzene, 3, 6) == 3  # para


def test_none_molecule_returns_none():
    assert fragments_excluding_bond(None, 0, 1) is None
    assert separate_without_atom(None, 0) is None
    assert ring_arc_distance(None, 0, 1) is None
