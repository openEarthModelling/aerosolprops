# EVAPORATION vapour pressure

**Status:** implemented (`aeroprops.methods.evaporation` for the pure
estimators and coefficient catalogue, `aeroprops.vapour_pressure.EVAPORATION`
for the object-oriented wrapper and the graph fragmenter).

## Reference

Compernolle, S., Ceulemans, K., & Müller, J.-F. (2011). EVAPORATION: a new
vapour pressure estimation method for organic molecules including
non-additivity and intramolecular interactions. *Atmospheric Chemistry and
Physics*, 11, 9431–9450.
[doi:10.5194/acp-11-9431-2011](https://doi.org/10.5194/acp-11-9431-2011)

## Scope

- Structure-only method: no boiling point input (the boiling point is an
  *output* of the method's Eq. 14, not a required input).
- log10(p°/Pa) = A + B/T^1.5 (paper Eq. 13) with A, B linear in 27 count
  classes: k = 1..20 plus the ring-functionality split `12lin`/`12CL`/`12HB`
  (paper Table 4 coefficients, intercepts, and the empirical diacid factor).
- Elements covered: C, H, O, and N only as alkyl nitrate or PAN; the paper
  mode refuses aromatics, other elements, and other nitrogen forms.
- Fitted against mostly 298 K data over zero-, mono- and polyfunctional
  SOA-relevant molecules (paper §1, §5).

## Quick start

```python
>>> from aeroprops.vapour_pressure import EVAPORATION
>>> eva = EVAPORATION(smi="CC(C)(C(=O)O)C(C(=O)O)CC(=O)O")   # MBTCA
>>> eva.counts
{'2': 8, '3': 3.0, '9': 3}
>>> round(eva.log10_vapour_pressure(298.15), 4)   # log10 Pa
-6.5127
>>> round(eva.boiling_point(), 1)                 # Eq. 14, K
687.8
>>> eva = EVAPORATION(smi="this is not a molecule")
>>> eva.success, eva.counts, eva.log10_vapour_pressure()
(False, {}, None)
```

Units are pinned as follows: `log10_vapour_pressure` returns log10(P/Pa)
and `vapour_pressure` returns Pa; the paper (and UManSysProp) work on an
atm basis and the conversion constant log10(101325) is applied inside the
pure estimator (`LOG10_ATM_TO_PA`). The pure layer accepts hand-built
count dictionaries (`evaporation_parameters`, `evaporation_boiling_temperature`,
`evaporation_log10_vapour_pressure` in `aeroprops.methods.evaporation`).

Failure modes return `None` (and set `success = False` on the wrapper):
unparseable SMILES, elements outside C/H/N/O, aromatic input (paper mode),
nitrogen that is not an alkyl nitrate or PAN, non-positive temperature,
and boiling points with no positive Eq. 14 solution.

## Counting modes

`evaporation_fragmentation(mol, upstream_compat=False)` implements the
paper's counting definitions (default) and, with `upstream_compat=True`,
the counting behaviour of the UManSysProp implementation. The differences,
all decoded from the paper text and from live UManSysProp behaviour (never
from its SMARTS/data files):

- k2 "in-chain oxygen": paper adds every single-bond O whose removal
  splits the carbon skeleton (ethers, ester bridges, peroxides — the
  latter twice); upstream adds two per ether-type oxygen instead.
- k3 branching: paper t = Σ max(0, n_C-neighbours − 2) minus ring-size
  corrections; upstream counts subsets of neighbour sets.
- k12 on rings: paper counts a ring functionality once per role (a lactone
  fires both `12lin` and `12CL`); upstream aggregates differently and a
  ring ester fires only the aggregate `12`.
- Nitrates/PAN are visible only in the paper mode (the upstream SMARTS
  never match RDKit-normalised charge forms).
- The paper's footnote f (k17 middle-carbon condition) exists only in the
  paper mode.

The default mode is the paper's; `upstream_compat` exists to reproduce
UManSysProp outputs exactly (see the oracle section below).

## Accuracy

Paper-level expectations (log10 units): for hydrocarbons and
monofunctional compounds the fit MAD is 0.085 and the predicted MAD 0.087
(paper §5.1); bifunctional compounds perform somewhat worse, especially
dinitrates, diacids, keto acids and hydroxy nitrates; functionalised
diacids motivated the empirical factor of Eq. 31. aeroprops reproduces
every paper-stated count/value anchor we could test (D-values for
neopentane, glyoxal, MBTCA, methyl nitrate/PAN, α-pinene, fenchyl alcohol,
pinonaldehyde, the Table 5 MTBE entry, footnote f) and predicted boiling
points of hexane, acetone, ethanol, neopentane and 1-butanol within 3% of
experiment.

## Implementation notes

- One graph walk over atom/bond neighbourhoods computes all count classes;
  `aeroprops.fragmentation.PriorityFragmenter` is deliberately *not* reused
  because EVAPORATION counting is overlapping/subset semantics (branching
  neighbour subsets, deduplicated pairs, per-role ring counting), not a
  mutually-exclusive SMARTS partition.
- "In-chain oxygen" is the removal-connectivity test of paper Eq. 26,
  computed with `aeroprops.topology.separate_without_atom`.
- Structures are normalised on entry (`aeroprops.structure.mol_from_smiles`,
  which applies RDKit's Normalizer) so that nitrate/PAN charge forms are
  detected; a neutral-written ON(O)O stays neutral and is still classified
  by neighbour topology.
- Table 4 transcription ruling: b20 = 1817.69 (paper) where the
  UManSysProp source tree carries 1817.60 — a typo there; aeroprops uses
  the paper value in both modes (within fit noise).

## Known limitations / open questions

- Upstream submodes that could not be observed from any probe molecule
  (305, 1604, 1704) are implemented per decoded semantics but never
  exercised by an oracle molecule.
- Aromatics are refused in the paper mode; in `upstream_compat` mode they
  are accepted but only partially validated (phenol verified).
- Two paper-anchor values could not be reproduced from the stated counts
  (pinonic acid −7.3048; the D4 oxaloacetic panel) and are documented, not
  forced, in the tests; the paper's supplement (D10 table) was not
  available to disambiguate.
