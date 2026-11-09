# Nannoolal boiling point (Part 1) and vapour pressure (Part 3)

**Status:** planned for v0.1 (not yet implemented).

## References

- Part 1 (boiling point): Nannoolal, Y., Rarey, J., Ramjugernath, D., &
  Cordes, W. (2004). Estimation of pure component properties Part 1.
  Estimation of the normal boiling point of non-electrolyte organic
  compounds via group contributions and group interactions. *Fluid Phase
  Equilibria*, 226, 45–63.
  [doi:10.1016/j.fluid.2004.09.001](https://doi.org/10.1016/j.fluid.2004.09.001)
- Part 3 (vapour pressure): Nannoolal, Y., Rarey, J., & Ramjugernath, D.
  (2008). Estimation of pure component properties Part 3. Estimation of the
  vapor pressure of non-electrolyte organic compounds via group
  contributions and group interactions. *Fluid Phase Equilibria*, 269,
  117–133. [doi:10.1016/j.fluid.2008.04.020](https://doi.org/10.1016/j.fluid.2008.04.020)
- SOA extension groups: Compernolle, S., et al. (2010). *Atmospheric
  Chemistry and Physics*, 10, 6271–6282.
  [doi:10.5194/acp-10-6271-2010](https://doi.org/10.5194/acp-10-6271-2010)

## Scope

- Part 1: structure-only normal boiling point (second-order groups plus
  group interactions). Used both as a standalone method and as the Tb input
  of Part 3.
- Part 3: vapour pressure from the normal boiling point and group
  contributions (requires Tb; Part 1 can supply it).
- The extension groups of Compernolle et al. (2010) add peroxide (−OO−),
  hydroperoxide (−OOH), peracid (−C(=O)OOH), and PAN (−C(=O)OONO₂)
  coverage; per that paper no group-interaction parameters exist for these
  groups, so they enter without interaction contributions.

## Implementation notes

- Overlapping groups are adjudicated by the papers' priority numbers
  (Table 1 footnote: "If two different groups can be used, the one with the
  lower priority number must be used"), not by pre-splitting SMARTS into
  mutually exclusive patterns.
- Aromatic ortho/meta/para pair groups use ring arc distances
  (`aeroprops.topology.ring_arc_distance`).
- The alcohol boiling-point correction of Part 3 Eq. (9) as printed is
  physically unrealisable (the sigmoid transition sits at ~2 K) and
  contradicts the paper's own reported statistics; it is therefore **not**
  implemented — matching the reference implementation's behaviour.
