# EVAPORATION vapour pressure

**Status:** planned for v0.1 (not yet implemented).

## Reference

Compernolle, S., Ceulemans, K., & Müller, J.-F. (2011). EVAPORATION: a new
vapour pressure estimation method for organic molecules including
non-additivity and intramolecular interactions. *Atmospheric Chemistry and
Physics*, 11, 9431–9450.
[doi:10.5194/acp-11-9431-2011](https://doi.org/10.5194/acp-11-9431-2011)

## Scope

- Structure-only method: no boiling point input (the boiling point is an
  *output* of the method's Eq. 14, not a required input).
- 27 counting-type groups with linear, halogen, and H-bond coefficient
  branches; group counts enter as `log10(p°/Pa) = intercept + Σ g_k·c_k`.

## Implementation notes

- The "in-chain oxygen" definition (an oxygen whose removal splits the
  carbon skeleton, as in ethers, esters, peroxides) is a graph-connectivity
  test; computed with `aeroprops.topology.separate_without_atom`.
- Structures are normalised on entry (`aeroprops.structure.mol_from_smiles`)
  so that nitrate/PAN charge forms match the SMARTS catalogues.
