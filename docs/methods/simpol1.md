# SIMPOL.1 vapour pressure

**Status:** planned for v0.1 (not yet implemented).

## Reference

Pankow, J. F., & Asher, B. E. (2008). SIMPOL.1: a simple group contribution
method for predicting vapor pressures and enthalpies of vaporization of
multifunctional organic compounds. *Atmospheric Chemistry and Physics*, 8,
2773–2796. [doi:10.5194/acp-8-2773-2008](https://doi.org/10.5194/acp-8-2773-2008)

## Scope

- Structure-only method: no boiling point or melting point input.
- `log10(p°/Pa)` as a linear function of 31 group counts.

## Implementation notes

- Group 2 counts the "acid-side" carbons of amides — a *global* graph
  property (unbounded acyl-chain length) that cannot be expressed as any
  single SMARTS pattern. It is computed with the connected-component
  primitive in `aeroprops.topology.fragments_excluding_bond`.
