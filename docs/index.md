# aeroprops

Structure-derived pure-component properties for atmospheric organic aerosol.

`aeroprops` is a modern successor to the pure-component property methods of
[UManSysProp](https://github.com/loftytopping/UManSysProp_public) (Topping et
al., GMD 9, 899–911, 2016), rebuilt on [RDKit](https://www.rdkit.org/) as the
structure engine and following the conventions of the
[`thermo`](https://github.com/CalebBell/thermo) package family: strict SI
units (Pa, K, mol), estimators as pure functions of group counts, and
failures that return `None` rather than raise.

## Methods

| Method | Property | Needs | Status |
|---|---|---|---|
| [SIMPOL.1](methods/simpol1.md) | Vapour pressure | structure only | implemented |
| [EVAPORATION](methods/evaporation.md) | Vapour pressure | structure only | implemented |
| [Nannoolal Part 1](methods/nannoolal.md) | Normal boiling point | structure only | implemented |
| [Nannoolal Part 3](methods/nannoolal.md) | Vapour pressure | + boiling point | implemented |
| [Myrdal–Yalkowsky](methods/myrdal_yalkowsky.md) | Vapour pressure | + Tb, Tm | implemented |

Boiling points used as inputs can come from Nannoolal Part 1, from
`thermo`'s built-in Joback method, or from experimental data; aeroprops does
not reimplement Joback.

## Conventions

- Strict SI units; failures return `None`.
- RDKit is the structure engine (required dependency, imported lazily).
  Estimators consume group-count dictionaries — the interface `thermo`
  uses for its own RDKit-optional group-contribution methods.
- Coefficients are entered from the published papers; SMARTS catalogues are
  authored from the papers' group definitions. UManSysProp is used only as a
  behavioural test oracle during development.

## License

MIT — see [LICENSE](license.md).
