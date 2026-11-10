# Myrdal–Yalkowsky vapour pressure

**Status:** implemented (`aerosolprops.methods.myrdal_yalkowsky` for the
pure estimators and coefficient catalogue,
`aerosolprops.vapour_pressure.MyrdalYalkowsky` for the object-oriented
wrapper).

## References

- Myrdal, P. B., & Yalkowsky, S. H. (1997). Estimating pure component
  vapor pressures of complex organic molecules. *Industrial &
  Engineering Chemistry Research*, 36, 2494–2499.
  [doi:10.1021/ie950242l](https://doi.org/10.1021/ie950242l)
- Hydrogen-bonding extension: Compernolle, S., et al. (2010).
  *Atmospheric Chemistry and Physics*, 10, 6271–6282.
  [doi:10.5194/acp-10-6271-2010](https://doi.org/10.5194/acp-10-6271-2010)

## Scope

- **Input method**: paper Eq. 12 needs the normal boiling point `Tb`,
  and the melting point `Tm` for the solid branch. `Tb` is routed from
  `NannoolalBP` (default), `thermo`'s built-in Joback method, or
  experimental data; without `Tm` the estimate is the subcooled-liquid
  vapour pressure (the paper's own Table 1 convention), annotated in
  the wrapper `status`.
- Descriptor method: Eq. 8 torsion count τ (SP3/SP2/RING walk) and the
  Eq. 9 hydrogen-bond number HBN on the 2D graph — element counts, ring
  systems and a handful of clean-room SMARTS; no group catalogue.
- Elements covered: C, H, N, O, S, F, Cl, Br, I. The Compernolle (2010)
  §2.6 extension brings hydroperoxides and peracids into the
  hydroxyl hydrogen-bonding root; enols, organic nitrates and PAN are
  outside the published parameterisation and are refused.
- log10(P/Pa) output; the paper fitted on an atm basis and the
  conversion constant log10(101325) is applied inside the pure
  estimator.

## Quick start

```python
>>> from aerosolprops.vapour_pressure import MyrdalYalkowsky
>>> my = MyrdalYalkowsky(smi="Oc1cccc(O)c1", tb=549.0)   # resorcinol
>>> my.status
'OK (sub-cooled liquid: Tm not supplied)'
>>> round(my.log10_vapour_pressure(425.0), 2)   # log10 Pa
3.24
>>> my = MyrdalYalkowsky(smi="CCCCCC")   # Tb from Nannoolal Part 1
>>> round(my.tb, 2)
350.03
>>> my = MyrdalYalkowsky(smi="CCCCCC", tb_source="joback")
>>> round(my.tb, 2)
336.88
```

The solid branch activates when a melting point is supplied and
`T < Tm` (σ, the external rotational symmetry number, enters only
there):

```python
>>> my = MyrdalYalkowsky(smi="c1ccccc1c1ccccc1", tb=528.0, tm=342.4, sigma=2)
>>> round(my.log10_vapour_pressure(342.0), 2)
2.04
```

## Accuracy

Paper Table 1 (19 compounds, experimental `Tb`): aerosolprops reproduces
every predicted value at the paper's printed two decimals, and against
the observed column the average absolute error is 0.09 log units
(factor 1.23) with a maximum of 0.39 — the paper's own summary
statistics. The estimator is expected to be rougher than SIMPOL.1 or
EVAPORATION (three descriptor terms instead of a group catalogue); its
niche is molecules whose groups those methods do not parameterise.

## Implementation notes

- The paper's Eq. 12 denominators print 19.1, the authors' rounding of
  2.303 R (19.147). The literal 19.1 is what reproduces Table 1 and is
  kept, with the identity pinned in tests.
- Descriptor conventions are locked against Table 1: SP3 = nonring,
  nonterminal sp3 atoms among C/N/O/S with tert-butyl centres and
  triple-bond atoms excluded; SP2 = nonring sp2 carbons (carbonyl
  carbons at the half weight); RING = **independent ring systems**
  (connected components of ring atoms through ring bonds — naphthalene
  1, biphenyl 2; deliberately not the SSSR count UManSysProp feeds in).
- HBN = (√(n_OH-family) + 0.33·√(n_NH2)) / MW with MW over **all**
  elements in g/mol (the paper divides by MW, not by a partial
  composition mass).
- The thermo Joback route gates on `Joback.success`:
  `Tb_estimated` holds a number even when the fragmentation failed.
- `Tm` exactly equal to `T` selects the liquid branch (paper Appendix A
  ruling on their biphenyl row).

## Known limitations / open questions

- Table 1's 2-methyl-5-ethylpyridine row does not discriminate τ = 0
  vs τ = 0.5 (both print −2.05); aerosolprops uses the paper-convention
  walk (τ = 0.5).
- Appendix A prints −1.95 for 3,5-diethylphenol (rounded stepwise from
  −1.72 − 0.23) while Table 1 prints −1.96 for the same compound; the
  unrounded value is −1.958.
- The σ (rotational symmetry number) is a caller input; no automatic
  perception is attempted.

The test suite also runs the unmodified UManSysProp
`myrdal_and_yalkowsky` as a behavioural oracle (skipped when its source
tree is absent): all 19 Table 1 molecules and 14 controls agree
bit-for-bit or at the last ulp with the same `Tb` inputs; the known
UMS deviations (SSSR ring counts, tert-butyl centres, sulfur in the
SP3 walk, enols/nitrates accepted) are whitelisted with paper evidence
in `tests/test_oracle_myrdal_yalkowsky.py`.
