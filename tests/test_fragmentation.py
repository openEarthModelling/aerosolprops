# SPDX-License-Identifier: MIT
"""Unit tests for the generic priority-claim fragmentation engine.

The engine implements the Nannoolal Table 1 footnote a semantics ("if
two different groups can be used, the one with the lower priority
number must be used") as group-level claiming, not as an a priori
mutually-exclusive split of the SMARTS catalogue.  The tests use small
synthetic catalogues so every adjudication rule is exercised in
isolation before the Nannoolal catalogue rides on top.
"""

import pytest

from aerosolprops.fragmentation import PriorityFragmenter
from aerosolprops.structure import mol_from_smiles


def test_lower_priority_number_wins_contested_atom():
    # Two groups both describe the methyl of methanol; the one with the
    # lower priority number must claim it (paper footnote a).
    engine = PriorityFragmenter({
        "plain": {"smarts": "[CH3]", "priority": 2, "claim": (0,)},
        "on-oxygen": {"smarts": "[CH3][OX2]", "priority": 1, "claim": (0,)},
    })
    for smi, expected in (("CC", {"plain": 2}), ("COC", {"on-oxygen": 2})):
        counts, _assign, ok, status = engine.fragment(mol_from_smiles(smi))
        # COC does not close the atom balance (the O is unparameterised)
        # but the carbon claims are still adjudicated correctly.
        if ok:
            assert (status, counts) == ("OK", expected)
        else:
            assert counts == expected and "unclaimed atom" in status


def test_multi_atom_claim_is_all_or_nothing():
    # A two-atom group loses its whole instance when a lower priority
    # number takes one of its atoms - the surviving atom is NOT left
    # half-claimed (functional groups count once or not at all).
    engine = PriorityFragmenter({
        # claims both atoms of a terminal C-C bond
        "pair": {"smarts": "[CH3][CH2]", "priority": 5, "claim": (0, 1)},
        # claims only the CH2
        "ch2": {"smarts": "[CH2]", "priority": 1, "claim": (0,)},
        "ch3": {"smarts": "[CH3]", "priority": 2, "claim": (0,)},
    })
    counts, _assign, ok, _status = engine.fragment(mol_from_smiles("CCC"))
    assert ok
    # both CH2 and the CH3s are taken by the single-atom groups, so the
    # pair group never fires
    assert counts == {"ch2": 1, "ch3": 2}


def test_unclaimed_atom_fails_with_status():
    engine = PriorityFragmenter({
        "me": {"smarts": "[CH3]", "priority": 1, "claim": (0,)},
    })
    counts, _assign, ok, status = engine.fragment(mol_from_smiles("CCO"))
    assert not ok
    assert counts == {"me": 1}  # the CH2 and the OH oxygen stay free
    assert status.startswith("unclaimed atom")
    assert "(" in status  # the element symbol of the free atom is named


def test_symmetric_matches_deduplicate_on_claim():
    # Both direction matches of the ether pattern yield the same claim
    # set; the group counts once.
    engine = PriorityFragmenter({
        "ether-o": {"smarts": "[#6][OX2][#6]", "priority": 1,
                    "claim": (1,)},
        "me": {"smarts": "[CH3]", "priority": 2, "claim": (0,)},
    })
    counts, assignments, ok, _status = engine.fragment(
        mol_from_smiles("COC"))
    assert ok
    assert counts == {"ether-o": 1, "me": 2}
    assert len(assignments["ether-o"]) == 1


def test_equal_priority_is_deterministic_by_group_id():
    engine = PriorityFragmenter({
        "b": {"smarts": "[CH3]", "priority": 1, "claim": (0,)},
        "a": {"smarts": "[CH3]", "priority": 1, "claim": (0,)},
    })
    counts, _assign, ok, _status = engine.fragment(mol_from_smiles("CC"))
    assert ok
    assert counts == {"a": 2}


def test_invalid_smarts_raises_at_fragmentation_time():
    engine = PriorityFragmenter({
        "bad": {"smarts": "not-a-smarts", "priority": 1, "claim": (0,)},
    })
    with pytest.raises(ValueError):
        engine.fragment(mol_from_smiles("CC"))


def test_none_molecule_fails_cleanly():
    engine = PriorityFragmenter({
        "me": {"smarts": "[CH3]", "priority": 1, "claim": (0,)},
    })
    counts, _assign, ok, status = engine.fragment(None)
    assert (counts, ok) == ({}, False)
    assert status == "Failed to construct mol"
