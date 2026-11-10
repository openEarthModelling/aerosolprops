# SIMPOL.1 vapour pressure

**Status:** implemented (`aerosolprops.methods.simpol1` for the pure
estimators and coefficient catalogue, `aerosolprops.vapour_pressure.SIMPOL1`
for the object-oriented wrapper).

## Reference

Pankow, J. F., & Asher, W. E. (2008). SIMPOL.1: a simple group contribution
method for predicting vapor pressures and enthalpies of vaporization of
multifunctional organic compounds. *Atmospheric Chemistry and Physics*, 8,
2773–2796. [doi:10.5194/acp-8-2773-2008](https://doi.org/10.5194/acp-8-2773-2008)

## Scope

- Structure-only method: no boiling point or melting point input, no
  descriptors — just group counts from the 2D graph.
- log10(P/Pa) as a linear function of 31 group counts (paper Table 5;
  indices `k = 0..30` follow the paper exactly).
- Vapourisation enthalpy ΔH_vap(T) in J/mol from the same coefficients
  (paper Eqs. 10–11).
- Elements covered: C, H, N, O only; anything else fails cleanly.
- Fitted temperature range 273.15–393.15 K
  (`SIMPOL1_TEMPERATURE_RANGE`); errors grow markedly below 300 K for
  low-volatility compounds (paper Figs. 10–12).

## Quick start

```python
>>> from aerosolprops.vapour_pressure import SIMPOL1
>>> sim = SIMPOL1(smi="C1=CCCCC1")            # cyclohexene
>>> sim.counts
{0: 1, 1: 6, 4: 1, 5: 1}
>>> round(sim.log10_vapour_pressure(293.15), 4)   # log10 Pa
4.0567
>>> round(sim.vapour_pressure(293.15), 1)         # Pa
11395.4
>>> round(sim.vapourisation_enthalpy(333.15) / 1000.0, 1)   # kJ/mol
33.4
```

Units are pinned as follows: `log10_vapour_pressure` returns
log10(P/Pa), `vapour_pressure` returns Pa, and `vapourisation_enthalpy`
returns J/mol. The paper fitted on an atm basis; the conversion constant
log10(101325) is applied inside the pure estimator
(`LOG10_ATM_TO_PA`).

The pure layer can be used with hand-built count dictionaries (the same
interface `thermo` uses for RDKit-optional group-contribution methods):

```python
>>> from aerosolprops.methods.simpol1 import simpol1_log10_vapour_pressure
>>> round(simpol1_log10_vapour_pressure({1: 6, 10: 2}, T=293.15) - 5.005716612413731, 2)  # adipic acid, log10 atm
-7.99
```

Failure modes return `None` (and set `success = False` on the wrapper):
unparseable SMILES, elements outside C/H/N/O, unknown group indices in a
hand-built counts dictionary.

## Accuracy

Paper-level expectations (on the log10 atm basis, T = 293.15 K):

- Cyclohexene: predicted −0.95, measured −1.08 (paper §3.3).
- Adipic acid: predicted −7.99, measured −8.49 (paper §3.3).
- Overall fit: σFIT = 0.34 (456 compounds), ΔHvap σ = 8.9 kJ/mol at
  333.15 K. Individual deviations at the low-volatility end
  (< 10⁻¹⁰ atm) are larger.

## Implementation notes

- The 124 coefficients live in `SIMPOL1_GROUPS` (one in-source literal,
  keyed by the paper's group index). They are cross-verified in tests
  against the independent b_k(293.15 K) column of paper Table 6 and its
  d(Δh_vap)/dT < 0 flags.
- Group 2 ("acid-side carbons" of amides) is a *global* graph property —
  unbounded acyl-chain length, no single SMARTS can express it. It is
  computed by cutting each amide C–N bond and counting the carbons of
  the connected component containing the acyl carbon
  (`aerosolprops.topology.fragments_excluding_bond`).
- Groups 3/4 (aromatic/non-aromatic rings) use minimal SSSR ring
  perception, which reproduces the paper's cyclomatic-number semantics
  (naphthalene 2, anthracene 3, biphenyl 2; fused bridged systems are
  *not* inflated — deliberately `GetSSSR`, not `GetSymmSSSR`).
- Counting semantics per paper Table 5 footnotes: match sets are
  deduplicated on counting anchors (atom or bond); carbonyl-class groups
  (10/11/25/28) claim their oxygens so ether/hydroxyl groups (7, 12/13/14)
  do not double-count them; peroxy-acid O–O bonds are claimed against
  groups 26/27; nitrophenol (29) *replaces* aromatic hydroxyl (17) and
  nitroester (30) *replaces* ester (11).

## Known limitations / open questions

Cases where the paper itself does not define the count; aerosolprops picks
the conservative choice and each is exercised by a test:

- **Lactams (cyclic amides)**: the C–N graph cut leaves the ring
  connected, so ν₂ counts every carbon of the component containing the
  acyl carbon (the whole ring, plus carbons of N-substituents that stay
  attached to it). Piperidin-2-one → ν₂ = 5,
  1-methyl-piperidin-2-one → ν₂ = 6.
- **Diamides / chain amides**: ν₂ is summed per amide, each seeing its
  own acyl component (succinamide → ν₂ = 8).
- **Compounds with no parameterised group** (e.g. hydroxylamine): only
  the zeroth group and carbon count contribute — the estimate is an
  extrapolation with no warning beyond the counts themselves.

The test suite also runs the unmodified UManSysProp `simpol` as a
behavioural oracle (skipped when its source tree is absent): all 26
molecules of its `Example_output.txt` agree on every group count and on
log10 P to 1e-12; the known UManSysProp deviations (amide/benzoate/
nitroester classes) are whitelisted with justifications in
`tests/test_oracle_umansysprop.py`.
