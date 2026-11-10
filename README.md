# aerosolprops

Structure-derived pure-component properties for atmospheric organic aerosol.

`aerosolprops` is a modern successor to the pure-component property methods of
[UManSysProp](https://github.com/loftytopping/UManSysProp_public)
(Topping et al., GMD 9, 899–911, 2016), rebuilt on
[RDKit](https://www.rdkit.org/) as the structure engine and following the
conventions of the [`thermo`](https://github.com/CalebBell/thermo) package
family: strict SI units, estimators as pure functions of group counts where
possible, and failures that return `None` rather than raise.

**Status: pre-alpha.** The repository currently contains the packaging
skeleton, CI, documentation scaffold, the structure/topology foundation
layers, and the first estimator (SIMPOL.1). The remaining methods below
land one at a time, each with coefficients entered from the published
papers and tests against literature values.

## Methods (v0.1)

| Method | Property | Status | Reference |
|---|---|---|---|
| SIMPOL.1 | Vapour pressure | implemented | Pankow & Asher, [ACP 8, 2773 (2008)](https://doi.org/10.5194/acp-8-2773-2008) |
| EVAPORATION | Vapour pressure | implemented | Compernolle et al., [ACP 11, 9431 (2011)](https://doi.org/10.5194/acp-11-9431-2011) |
| Nannoolal Part 1 (+ SOA extension groups) | Normal boiling point | implemented | Nannoolal et al., [FPE 226, 45 (2004)](https://doi.org/10.1016/j.fluid.2004.09.001); Compernolle et al., [ACP 10, 6271 (2010)](https://doi.org/10.5194/acp-10-6271-2010) |
| Nannoolal Part 3 | Vapour pressure | implemented | Nannoolal et al., [FPE 269, 117 (2008)](https://doi.org/10.1016/j.fluid.2008.04.020) |
| Myrdal–Yalkowsky (+ H-bonding extension) | Vapour pressure | implemented | Myrdal & Yalkowsky, [IECR 36, 2494 (1997)](https://doi.org/10.1021/ie950242l); Compernolle et al., [ACP 10, 6271 (2010)](https://doi.org/10.5194/acp-10-6271-2010) |

Boiling points used as inputs (by Myrdal–Yalkowsky and Nannoolal Part 3) can
come from the self-implemented Nannoolal Part 1 route, from `thermo`'s
built-in Joback method, or from experimental data — `aerosolprops` does **not**
reimplement Joback.

## Installation

```bash
pip install aerosolprops
```

RDKit and `thermo` are required dependencies.

## Design notes

- **Units are strictly SI** (Pa, K, mol); property failures return `None`.
- **RDKit is the structure engine and a required dependency.** It is
  imported lazily, and estimators ultimately consume group-count
  dictionaries — the same interface `thermo` uses for its own
  RDKit-optional group-contribution methods — so pieces remain extractable
  upstream.
- **No molecular descriptors.** Every method needs only 2D graph quantities:
  SMARTS substructure counts, ring perception, composition, and a few
  connectivity queries. No fingerprints, no 3D, no descriptor packages.
- **No home-grown graph algorithms.** The graph primitives (connected
  components, ring distances, neighbour queries) are RDKit built-ins;
  `aerosolprops.topology` only wraps them. What is method-specific is the
  translation of paper semantics (e.g. SIMPOL.1 "acid-side carbons",
  EVAPORATION "in-chain oxygen") onto those primitives.
- **`thermo` compatibility.** aerosolprops hard-depends on `thermo` and follows
  its group-contribution conventions (`{group_id: coefficient record}`
  catalogues, counts-based static estimators, Joback-style one-time
  fragmenting wrappers). Where a method's counting semantics do not fit
  `thermo`'s mutually-exclusive SMARTS partition (SIMPOL.1's groups overlap
  and replace one another additively), aerosolprops ships its own fragmenter
  but keeps the same interface, so pieces can still be contributed upstream
  with minimal rework.

## Licensing and data provenance

`aerosolprops` is MIT-licensed. Coefficients are entered from the published
papers cited above (facts, not copyrightable expression), and SMARTS
patterns are authored from the group definitions in those papers.
UManSysProp is used **only** as a behavioural test oracle during development;
no code or data files are copied from it.

## Development

```bash
uv venv && uv pip install -e ".[test,docs]"
pytest                        # unit tests + doctests
mkdocs serve                  # documentation
```

CI runs the full test matrix (with doctests) on Linux/macOS/Windows for
Python 3.9–3.14, plus a strict docs build.

## License

MIT — see [LICENSE](LICENSE).
