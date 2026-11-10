# Nannoolal boiling point (Part 1) and vapour pressure (Part 3)

**Status:** both parts are implemented. Part 1:
`aeroprops.methods.nannoolal_bp` (pure estimator and coefficient
catalogue) and `aeroprops.boiling_point.NannoolalBP` (object-oriented
wrapper). Part 3: `aeroprops.methods.nannoolal_vp` (pure estimator and
coefficient catalogue) and `aeroprops.vapour_pressure.NannoolalVP`
(wrapper), which takes its boiling-point input from Part 1 by default,
from `thermo`'s Joback method, or from experimental data.

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
  of Part 3. Contributions are in kelvin; `boiling_point()` returns K.
- Part 1 has 132 groups: 113 paper main groups (IDs 1–117), 16
  second-order corrections (118–133), and the 3 SOA extension groups of
  Compernolle et al. (301–303). Extension groups enter without
  interaction contributions (per that paper no GI parameters exist for
  them); a molecule containing one is estimated with GI = 0 and flagged
  (`gi_status == 'no-GI (extension groups present)'`).
- Part 1 has 76 group-interaction parameters between 19 interaction
  classes (paper Table 5); the two pairs printed without a value
  (B–G, R–R) make the estimation *fail* rather than silently
  contributing zero.
- Part 3 (vapour pressure): paper Eq. 6,
  log10(P/1 atm) = (4.1012 + dB)·(T_rb − 1)/(T_rb − 1/8) with
  T_rb = T/Tb, where dB (Eq. 7) is a group sum plus the Table 6 group
  interactions minus 0.176055. Output is log10(P/Pa) (the atm→Pa
  constant is applied inside the estimator); Appendix A inverts Eq. 6
  for the boiling point at a given pressure
  (`boiling_point_from_pressure`).
- Part 3 group set: the Part 1 catalogue with exactly three
  fragmentation deltas — the 1,2-diketone main group 118 (valueless:
  counted, then refused), the conjugation correction renumbered
  118 → 134, and the aromatic carbon split 18 (three aromatic bonds,
  fused) vs 214 (two aromatic bonds plus an aliphatic bridge bond,
  e.g. biphenyl). The para-pair group 129 is likewise counted and then
  refused (Table 5 prints no value for it); Table 6 row 218 (C, M) is
  printed "Do not estimate" and fails rather than zero-filling.
- Elements covered: C, N, O, S, F, Cl, Br, I. The paper also
  parameterises Si/Ge/Sn/P/As families — out of scope for aeroprops and
  rejected by the element gate.

## Quick start

```python
>>> from aeroprops.boiling_point import NannoolalBP
>>> nnp = NannoolalBP(smi="CCO")                  # ethanol
>>> nnp.counts
{1: 1, 7: 1, 36: 1}
>>> nnp.n                                          # non-H atoms, all elements
3
>>> round(nnp.boiling_point(), 4)                  # K
333.0861
>>> nnp.gi_status
'OK'
```

Paper Table 17 worked examples reproduce through the real fragmenter
(methyl *m*-toluate, 17d; the paper prints 490.6, a display truncation
of the same group sum):

```python
>>> from aeroprops.boiling_point import nannoolal_fragmentation
>>> from aeroprops.structure import mol_from_smiles
>>> nannoolal_fragmentation(mol_from_smiles("COC(=O)c1cccc(C)c1"))
({2: 1, 3: 1, 15: 4, 16: 2, 45: 1, 118: 1, 128: 1}, 11, True, 'OK')
>>> from aeroprops.methods.nannoolal_bp import nannoolal_bp_boiling_point
>>> round(nannoolal_bp_boiling_point({1: 6, 4: 2, 6: 2, 133: 1}, 10), 1)
429.5
```

The pure layer takes plain count dictionaries (n = number of non-H
atoms of all elements — the paper counts sulfur too):

```python
>>> from aeroprops.methods.nannoolal_bp import nannoolal_bp_gi
>>> nannoolal_bp_gi({34: 2, 42: 1}, 9)     # Table 17b footnote: GI = 96.193 K
(96.19312222222223, 'OK')
```

Failure modes return `None` (and set `success = False` on the wrapper):
unparseable SMILES, elements outside the parameterised set, an atom no
group claims (out of scope by construction, e.g. the pyrazole NH or a
sulfonate ester), and the two parameter-less interaction pairs.

## Part 3 quick start

```python
>>> from aeroprops.vapour_pressure import NannoolalVP
>>> w = NannoolalVP(smi="CCO")    # ethanol; Tb routed from Part 1
>>> w.counts, round(w.db, 7)
({1: 1, 7: 1, 36: 1}, 0.6311822)
>>> round(w.tb, 4)                # K, the Part 1 estimate
333.0861
>>> round(w.log10_vapour_pressure(298.15), 4)   # log10(Pa)
4.3612
>>> round(w.vapour_pressure(298.15), 1)         # Pa
22971.4
```

Paper worked example 15a (α-pinene at 388.15 K, Tb = 429.0 K)
reproduces through the real fragmenter:

```python
>>> vp = NannoolalVP(smi="CC1=CCC2CC1C2(C)C", tb=429.0)
>>> (vp.success, round(vp.db, 7))
(True, 0.1078596)
>>> round(10.0 ** vp.log10_vapour_pressure(388.15) / 1000.0, 2)  # kPa
31.03
```

Refusals never return zero: glycine carries a COOH and an NH2, and
Table 6 row 218 (C, M) is printed "Do not estimate":

```python
>>> g = NannoolalVP(smi="NCC(=O)O", tb=501.5)
>>> g.log10_vapour_pressure(298.15) is None
True
>>> g.status
'refused: missing interaction parameter (C, M)'
```

The same refusal discipline applies to the valueless groups 106/108/118
(imidazole-type N, thiocyanate carbon, 1,2-diketone) and the para pair
129 — all counted by the fragmenter, then refused by the estimator —
and to `Tb` routing failures (experimental `tb` missing, Joback
unavailable, Part 1 refusal), which propagate into `status` with the
route named.

## Accuracy

- Paper worked examples 17a–17d: 429.5 / 509.3 / 246.3 / 490.7 K
  (each within 0.1 K of the printed value through the real fragmenter).
- Literature checks (CRC Handbook): furfural, anisole, acetophenone,
  quinoline within 1.5 K; morpholine 4.4 K low (the cyclic-amine group
  97; the paper's mean absolute error per group is 4–8 K).
- Part 1 overall (paper §4): absolute average deviation 6.52 K over
  2812 regression compounds; 6.37 K mean absolute deviation on the
  199-compound external validation set.
- Part 3: all eight worked examples of Tables 15a–15h reproduce
  end-to-end through the real fragmenter (dB to the printed seven
  decimals, pressures within the paper's own display rounding); the
  15f boiling-point inversion via the Appendix A form returns 520.85 K
  as printed. The paper's §5 validation reports average relative
  deviations against the Antoine equation per compound family and
  pressure range (e.g. hydrocarbons 4.0% overall vs Antoine's 3.4%,
  with this work clearly better in the very-low-pressure range);
  aeroprops inherits that accuracy by construction of the coefficients.

## Implementation notes

- Overlapping groups are adjudicated by the paper's priority numbers
  (Table 1 footnote a: "If two different groups can be used, the one
  with the lower priority number must be used") through the reusable
  engine `aeroprops.fragmentation.PriorityFragmenter` — group-level
  claiming with multi-atom all-or-nothing claims, not pre-split SMARTS.
- Corrections 118–133 that no SMARTS can express (conjugation,
  o/m/p ring arc distances, ring sizes, steric bond crowding) are
  computed by graph helpers after the engine pass; 131–133 count the
  further carbon neighbours of *both endpoints* of single C(sp3)–C(sp3)
  bonds, ring bonds included (the paper's own example for 132, ethyl
  bornyl ether, reaches 3+2 only on the bornane bridge bond).
- The secondary-amine split follows the Part 3 revision of the group
  table: 42 = acyclic R₂NH, 97 = cyclic (ring) R₂NH. The Part 1 text
  ("at least one aromatic neighbour") contradicts its own example
  molecules (morpholine, pyrrolidine); 66/67 are the H0 "free electron
  pair" aromatic N only, so a pyrrole-type ring NH is 97.
- Aromatic ortho/meta/para pair groups use ring arc distances
  (`aeroprops.topology.ring_arc_distance`), fused atoms excluded,
  counted once per molecule with the minimum distance winning.
- Coefficient transcription rulings (documented in tests): Table 1
  ID 12 misprint 492.0707 (duplicate of ID 109) resolved to the Table 3
  regression value 250.9584; ID 24 minus sign restored (−19.5575).
- The test suite runs the unmodified UManSysProp `nannoolal` as a
  behavioural oracle (skipped when its source tree is absent): 9 of the
  26 `Example_output.txt` molecules are bit-for-bit identical (counts
  and Tb), the remainder plus an 18-molecule control whitelist are
  documented UManSysProp deviations with paper-based reasons in
  `tests/test_oracle_nannoolal.py` — including UMS silently dropping
  the N(O)O fragment and entire group families (aromatic N 67, bridge
  head 18, amide 48, S groups, perfluoro family).
- Part 3 reuses the same `PriorityFragmenter` engine with a Part 3
  catalogue: the Part 1 SMARTS plus the three deltas above, priorities
  taken from the Part 3 tables (extension groups 303/302/301 keep
  priorities 1/2/3 as in Part 1). The OH-family Table 1 short-chain
  rule (fewer than 5 carbons → group 36) is applied per the Part 3
  group definitions before the neighbour-count classes, plus the
  gem-diol fallback.
- Part 3 transcription adjudications (documented in tests): Table 4
  rows 74/75 print 510.9666 (nitrite) and 1317.4360 (oxime); the
  assignment is confirmed by the independent Part 1 fits. Questionable
  values (ID 91, the (E,E) interaction 187) are recorded verbatim.
  Groups 98/110 (printed without values in the Part 1 numbering)
  carry parked values in `NANNOOLAL_VP_UNDEFINED_VALUE_GROUPS` but stay
  excluded.
- The Part 3 oracle suite (`tests/test_oracle_nannoolal_vp.py`) runs
  UMS's `vapour_pressures.nannoolal` live with the Example_output.txt
  protocol (T = 298.15, Tb from UMS's own Part 1): the oracle
  reproduces the shipped golden column bit-for-bit for all 26
  molecules; 3 are bit-for-bit with aeroprops, 2 last-ulp (dB
  accumulation order), and every remaining molecule plus 60 controls is
  either identical or whitelisted with paper evidence — including
  three UMS data-file misprints (group 7 as 0.0378487 vs the printed
  37.8485e-3, group 48 as 0.2692470 vs 269.2471e-3, extension 303 as
  0.519 vs Compernolle's 519.03e-3), the para pair 129 carried as 0.0,
  Table 6 row 218 zero-filled, and the extension groups stuffed into
  GI classes A/C/F where the paper has no parameters for them.

## Known limitations / open questions

- **Group 109 ("N–(C=O)", thioester-side nitrogen) is not
  implemented**: the paper's example molecule for it contains no
  nitrogen, so its SMARTS cannot be written from the paper text; such
  molecules refuse. Group 48 (disubstituted amide) covers ordinary
  amides.
- **Sulfonate esters R–SO₂–O–R′** match no group (82 sulfone needs two
  carbons on S, 104/105 need two oxygens/N) — a gap in the paper's set;
  clean refusal.
- **Anhydrides** are covered by groups 76/96 without the conjugation
  correction 118 (the paper prints 118 only for aldehyde/ketone-style
  C=O connected to sp2 carbon).
- **Mixed ring/exocyclic alkenes** (e.g. isopropylidene cycloalkanes,
  =C(CH₃)₂ on a ring) fall between the C(c)=C(c) and C(r)=C(r)
  families; they currently fragment to the chain family with the ring
  carbon as (c)-substituted.
- **Part 3 Eq. (9)** (alcohol boiling-point correction) as printed is
  physically unrealisable and contradicts the paper's own statistics;
  it is not implemented (locked ruling, recorded in tests).
- **Part 3 valueless groups**: 106 (imidazole-type N), 108 (thiocyanate
  C), 118 (1,2-diketone) and 129 (para pair) are counted and then
  refused — never silently zeroed; ID 109 is excluded entirely (its
  paper definition names a molecule containing no nitrogen, so no
  SMARTS can be written from the paper text, mirroring the Part 1
  ruling).
- **Eq. 6 extrapolation**: for T above the supplied Tb the wrapper
  warns (Eq. 6 is a normal-boiling-point-anchored correlation) but
  still returns the extrapolated value; the Appendix A inversion
  `boiling_point_from_pressure` accepts any sub-critical pressure.
