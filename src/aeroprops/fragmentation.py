# SPDX-License-Identifier: MIT
"""Priority-claim SMARTS fragmentation engine.

A small, method-agnostic fragmenter implementing the priority semantics
used by the Nannoolal paper series (Part 1, Table 1, footnote a):

    "Priority numbers define the importance of a group. If two different
    groups can be used, the one with the lower priority number must be
    used."

Groups are SMARTS patterns that *claim* a set of molecular atoms (a
subset of the pattern's atoms - functional groups typically leave their
flanking neighbours to the carbon-skeleton groups). Every non-hydrogen
atom of the molecule must end up claimed by exactly one group instance.
Overlaps are adjudicated strictly by ascending priority number: an
instance counts only if it wins every atom it claims, i.e. no group with
a lower priority number has already claimed any of them.

Why not ``thermo.group_contribution.smarts_fragment_priority``? Its
priority order (highest value first) is compatible in spirit - and the
same greedy "accept a match only if disjoint from already-claimed atoms"
rule is used here - but its claim model takes the *whole* pattern match
as the claim, cannot express anchor-subset claims (e.g. an ester group
that claims C, O, O but not the flanking carbons), has no place for
non-atom-claiming second-order corrections (Nannoolal 118-133), for
groups without a priority number, or for refusing molecules with atoms
no group parameterises. aeroprops therefore ships this engine, keeping
the thermo idiom of one-time fragmentation with ``counts``/``success``/
``status`` results.

The engine itself knows nothing about any estimation method; method
modules provide a catalogue (see
:func:`aeroprops.boiling_point.nannoolal_fragmentation` for the first
user) and post-processing hooks.

Examples
--------
>>> engine = PriorityFragmenter({
...     "me": {"smarts": "[CH3]", "priority": 2, "claim": (0,)},
...     "me-on-o": {"smarts": "[CH3][OX2]", "priority": 1, "claim": (0,)},
...     "o": {"smarts": "[OX2]([#6])[#6]", "priority": 3, "claim": (0,)},
... })
>>> counts, assignments, ok, status = engine.fragment(
...     mol_from_smiles("COC"))
>>> ok, status, counts
(True, 'OK', {'me-on-o': 2, 'o': 1})
"""

from .structure import load_rdkit, mol_from_smiles  # noqa: F401 (doctest)

__all__ = ["PriorityFragmenter"]


class PriorityFragmenter(object):
    """SMARTS fragmenter with paper-style lowest-priority-wins claims.

    Parameters
    ----------
    catalogue : dict
        Group definitions keyed by an arbitrary hashable group ID.
        Each value is a dict with the keys:

        * ``smarts`` - str (compiled once per process); no ``None``
          placeholders here: purely algorithmic groups are computed by
          the method module, not the engine.
        * ``priority`` - int or float. Lower wins conflicts. Groups
          without a priority number (Nannoolal corrections/extensions)
          must not enter the catalogue.
        * ``claim`` - tuple of pattern-atom indices that the group
          instance claims when it matches. Instances are deduplicated
          on their claim tuple (symmetric matches of one functional
          group collapse to a single count).

    Attributes
    ----------
    catalogue : dict
    """

    def __init__(self, catalogue):
        self.catalogue = catalogue
        self._patterns = None

    def _compiled(self):
        """Compile the catalogue SMARTS once per fragmenter instance."""
        if self._patterns is None:
            loaded = load_rdkit()
            if loaded is None:
                return None
            patterns = {}
            for gid, spec in self.catalogue.items():
                pattern = loaded[0].MolFromSmarts(spec["smarts"])
                if pattern is None:
                    raise ValueError(
                        "SMARTS for group %r does not parse: %r"
                        % (gid, spec["smarts"])
                    )
                patterns[gid] = pattern
            self._patterns = patterns
        return self._patterns

    def fragment(self, mol):
        """Fragment `mol` into exclusive priority-claimed groups.

        Parameters
        ----------
        mol : Mol
            Normalised RDKit molecule (as produced by
            :func:`aeroprops.structure.mol_from_smiles`).

        Returns
        -------
        counts : dict
            Group ID -> number of winning instances.
        assignments : dict
            Group ID -> list of claim-atom tuples (original atom
            indices), for post-processing hooks and debugging.
        success : bool
        status : str
            ``"OK"`` when every non-hydrogen atom is claimed exactly
            once; otherwise the first reason found, e.g. ``"unclaimed
            atom 4 (O)"``.

        Examples
        --------
        Priority adjudication (paper footnote a): with two applicable
        groups, the one with the lower priority number claims the atom;

        >>> engine = PriorityFragmenter({
        ...     "plain-me": {"smarts": "[CH3]", "priority": 2, "claim": (0,)},
        ...     "me-on-o": {"smarts": "[CH3][OX2]", "priority": 1, "claim": (0,)},
        ... })
        >>> counts, assignments, ok, status = engine.fragment(
        ...     mol_from_smiles("CC"))
        >>> (ok, status, counts)
        (True, 'OK', {'plain-me': 2})
        """
        patterns = self._compiled()
        if mol is None or patterns is None:
            return {}, {}, False, "Failed to construct mol"

        # Candidate instances per group: deduplicated on the claim tuple.
        instances = {}
        for gid, spec in self.catalogue.items():
            claims = set()
            for match in mol.GetSubstructMatches(patterns[gid], uniquify=True):
                claims.add(tuple(sorted(match[i] for i in spec["claim"])))
            if claims:
                instances[gid] = claims

        # Lowest priority number claims first; an instance survives only
        # if every atom it claims is still free (footnote a semantics).
        order = sorted(
            instances, key=lambda gid: (self.catalogue[gid]["priority"], gid)
        )
        counts, assignments, taken = {}, {}, set()
        for gid in order:
            for claim in sorted(instances[gid]):
                if any(a in taken for a in claim):
                    continue
                counts[gid] = counts.get(gid, 0) + 1
                assignments.setdefault(gid, []).append(claim)
                taken.update(claim)

        if len(taken) != mol.GetNumAtoms():
            free = sorted(set(range(mol.GetNumAtoms())) - taken)
            return counts, assignments, False, "unclaimed atom %d (%s)" % (
                free[0],
                mol.GetAtomWithIdx(free[0]).GetSymbol(),
            )
        return counts, assignments, True, "OK"
