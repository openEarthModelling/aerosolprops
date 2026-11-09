# Myrdal–Yalkowsky vapour pressure

**Status:** planned for v0.1 (not yet implemented).

## References

- Myrdal, P. B., & Yalkowsky, S. H. (1997). Estimating pure component vapor
  pressures of complex organic molecules. *Industrial & Engineering
  Chemistry Research*, 36, 2494–2499.
  [doi:10.1021/ie950242l](https://doi.org/10.1021/ie950242l)
- Hydrogen-bonding extension: Compernolle, S., et al. (2010). *Atmospheric
  Chemistry and Physics*, 10, 6271–6282.
  [doi:10.5194/acp-10-6271-2010](https://doi.org/10.5194/acp-10-6271-2010)

## Scope

- Requires the normal boiling point `Tb` and melting point `Tm` as inputs
  (both obtainable from other routes; without `Tm` the method returns the
  subcooled-liquid vapour pressure, matching the reference implementation's
  structural behaviour).
- H-bonding terms via the count of hydrogen-bonding groups per unit mass;
  the extension of Compernolle et al. (2010, §2.6) brings hydroperoxides
  and peracids into the hydroxyl hydrogen-bonding root. Molecules
  containing groups outside the published parameterisation (nitrates, PAN,
  enols) are reported as unsupported and fall through to other methods.

## Implementation notes

- Composition and molecular weight over all elements are computed by
  `aeroprops.structure.element_counts` / `molecular_weight`.
