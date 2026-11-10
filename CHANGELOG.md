# Changelog

All notable changes to this project are documented in this file.
The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/)
and the project adheres to [semantic versioning](https://semver.org/).

## [Unreleased]

## [0.1.0] - 2026-11-10

### Changed

- Renamed the package from `aeroprops` to `aerosolprops` (import name and
  project name) ahead of the first public release; the oracle-suite
  environment variable is now `AEROSOLPROPS_UMS_REPO`. No release was
  ever published under the old name.

### Added

- Nannoolal Part 3 vapour pressure (Nannoolal, Rarey & Ramjugernath,
  FPE 269, 117-133, 2008; SOA extension groups from Compernolle et
  al., ACP 10, 6271, 2010): the fifth and final estimator method.
  - `aerosolprops.methods.nannoolal_vp` (pure layer): the Eq. 6/7
    estimators (`nannoolal_vp_dB`, `nannoolal_vp_log10_vapour_pressure`,
    `nannoolal_vp_boiling_point` via the Appendix A inversion), the
    Table 4 dB catalogue (115 groups; the three Part 3 fragmentation
    deltas vs Part 1 — the valueless 1,2-diketone 118, conjugation
    renumbered to 134, aromatic carbon split 18/214), Table 5
    second-order values, the Table 6 interaction matrix (45 rows, 19
    classes; row 218 (C, M) "Do not estimate" fails rather than
    zero-fills) and the Compernolle extension values on the no-GI
    branch.
  - `aerosolprops.vapour_pressure` (wrapper layer):
    `nannoolal_vp_fragmentation` — the Part 3 catalogue on the shared
    `PriorityFragmenter` engine with the OH-family short-chain (<5 C)
    and gem-diol reassignments — and the `NannoolalVP` class routing
    `Tb` from Nannoolal Part 1 (the default), `thermo`'s Joback or
    experimental data, with `db`/`db_status` eager in the constructor
    and `T > Tb` extrapolation warnings. Valueless groups (106, 108,
    118, 129) are counted then refused, never zeroed.
  - Verification: coefficient triple-guard (two independent
    transcriptions of Tables 4/5/6 plus census/priority tests), all
    eight worked examples 15a-15h through the real fragmenter (dB to
    seven decimals; the 15g printed dB is a misprint resolved by the
    formula; the 15h atom-count line is a paper typo), literature
    regression rows, counting/atom-balance semantics over ~100
    molecules, read-only fragmentation, three-route/failure gates on
    the wrapper, and a live UManSysProp oracle comparison (all 26
    Example_output goldens reproduce the shipped column; 60 controls;
    evidence-backed whitelist recording three UMS data misprints for
    groups 7/48/303, the 129-as-zero and row-218 zero-fill deviations,
    and the extension-GI stuffing).
- Myrdal-Yalkowsky vapour pressure (Myrdal & Yalkowsky, IECR 36, 2494,
  1997; hydrogen-bonding extension from Compernolle et al., ACP 10,
  6271, 2010): the third estimator method.
  - `aerosolprops.methods.myrdal_yalkowsky` (pure layer): the Eq. 2/7-12
    coefficient literals (the paper's printed 19.1 denominators, the
    authors' rounding of 2.303 R — kept because it is what reproduces
    Table 1), the tau/HBN/entropy descriptors and the Eq. 12
    log10-vapour-pressure estimator with solid/liquid branch selection.
  - `aerosolprops.vapour_pressure` (wrapper layer):
    `myrdal_yalkowsky_descriptors` — the Eq. 8 SP3/SP2 walk (nonring,
    nonterminal, tert-butyl centres and triple bonds excluded,
    heteroatoms included), RING as independent ring systems
    (connected components, not SSSR), the Eq. 9 hydrogen-bond group
    catalogue as clean-room SMARTS with the Compernolle hydroperoxide/
    peracid extension in the first root, all-element MW in g/mol, and
    refusal of enols/nitrates/PAN/out-of-domain elements — plus the
    `MyrdalYalkowsky` class routing `Tb` from Nannoolal Part 1 (the
    default), `thermo`'s Joback (gated on `success`:
    `Tb_estimated` holds a number even on failed fragmentation) or
    experimental data; `Tm` defaults to the sub-cooled liquid branch
    with a status annotation; `T > Tb` warns before extrapolating.
  - Verification: coefficient triple-guard, all 19 Table 1 rows through
    the pure estimator (exact at the printed two decimals) and the
    paper's own summary statistics against the observed column
    (AAE 0.09, max 0.39), a 67-molecule descriptor battery with
    hand-checked ring systems and Hill formulas, convention unit tests
    (including the ring-system-vs-SSSR contrast via `GetSSSR`, not the
    symmetry-inflated `GetSymmSSSR`), the three Tb routes with drift
    locks, failure paths, and a live UManSysProp oracle comparison
    (19/19 Table 1 molecules plus 14 controls matching bit-for-bit or
    at the last ulp; 6 whitelisted UMS deviations with paper evidence:
    SSSR ring counts, tert-butyl centres, sulfur in the SP3 walk,
    enols and nitrates accepted by UMS's wider tables).

- EVAPORATION vapour pressure (Compernolle, Ceulemans & Müller, ACP 11,
  9431, 2011): the fourth estimator method.
  - `aerosolprops.methods.evaporation` (pure layer): the 40+3 coefficients
    of paper Table 4 (27 count classes incl. the 12lin/12CL/12HB ring
    split, intercepts, diacid factor) cross-verified against the paper's
    independent g_k column, the A/B parameter and Eq. 14 boiling-point
    estimators, and `evaporation_log10_vapour_pressure` (log10 Pa basis).
    Transcription rulings recorded in tests: b20 = 1817.69 (paper) vs
    1817.60 in the UManSysProp tree (typo there).
  - `aerosolprops.vapour_pressure` (wrapper layer):
    `evaporation_fragmentation` with two counting conventions — the
    paper text (default: Eq. 26 in-chain oxygen via removal
    connectivity, branching descriptor t, per-role ring functionalities,
    nitrate/PAN classification, footnote f) and an `upstream_compat`
    mode reproducing UManSysProp counting bit-for-bit — plus the
    Joback-style `EVAPORATION` class. `PriorityFragmenter` deliberately
    not reused: the counting is overlapping/subset semantics, not a
    SMARTS partition.
  - Verification: paper D-value anchors (Table 5 and §6 panels), 5
    literature boiling points, 50-molecule atom-accounting and
    element-conservation suite, per-rule counting-semantics units, and a
    live UManSysProp oracle gate — 26/26 Example_output molecules with
    identical upstream-mode counts and log10 P agreement to 1e-12 (empty
    divergence whitelist), plus 51 probe molecules locked twice against
    recorded and live oracle counts.
- Nannoolal Part 1 normal boiling point (Nannoolal et al., FPE 226, 45,
  2004; SOA extension groups from Compernolle et al., ACP 10, 6271,
  2010): the second estimator method.
  - `aerosolprops.methods.nannoolal_bp` (pure layer): the 132-group
    coefficient catalogue (113 main + 16 corrections + 3 extensions;
    the two Table 1 transcription rulings are documented in tests),
    the 76 group-interaction parameters of Table 5 with the
    counting-form GI estimator `nannoolal_bp_gi` (parameter-less pairs
    B-G and R-R fail; extension groups switch to a flagged no-GI
    branch), and `nannoolal_bp_boiling_point`.
  - `aerosolprops.fragmentation`: reusable priority-number fragmenter
    (paper Table 1 footnote a semantics — lower priority number wins,
    group-level claiming, multi-atom all-or-nothing claims, atom
    balance check) for any catalogue of (SMARTS, priority, claim).
  - `aerosolprops.boiling_point`: clean-room SMARTS catalogue plus graph
    helpers for the corrections that SMARTS cannot express
    (conjugation, o/m/p ring arc distances, ring sizes, steric bond
    crowding), and the Joback-style `NannoolalBP` wrapper with
    one-time fragmentation and None-on-failure estimators.
  - Verification: paper Tables 17a-17d end-to-end, coefficient
    triple-guard, 5 literature boiling points, 50+-molecule atom
    balance, priority adjudication and refusal paths, and a live
    UManSysProp behavioural oracle (9 of 26 Example_output molecules
    bit-for-bit; all divergences whitelisted with paper-based reasons).
- SIMPOL.1 vapour pressure (Pankow & Asher, ACP 8, 2773, 2008): the
  first estimator method.
  - `aerosolprops.methods.simpol1` (pure layer): the 124 coefficients of
    paper Table 5 as an in-source catalogue, the SMARTS catalogue with
    algorithmic-group placeholders, and pure estimators
    `simpol1_bk`, `simpol1_log10_vapour_pressure` (log10 Pa basis),
    and `simpol1_vapourisation_enthalpy` (J/mol).
  - `aerosolprops.vapour_pressure` (wrapper layer): `simpol1_fragmentation`
    and the Joback-style `SIMPOL1` class with one-time fragmentation and
    `None`-on-failure estimators. Fragmentation implements the paper's
    overlapping-additive counting semantics (O-claims of carbonyl
    classes, per-OH nitrophenol replacement of phenol, nitroester
    replacement of ester, amide acid-side graph cut, SSSR ring counts).
  - `aerosolprops.structure.normalize_mol`: public normalisation entry
    point for caller-supplied RDKit molecules (extracted from
    `mol_from_smiles`).
- Verification: coefficient completeness against paper Tables 5 and 6
  (b_k column and d(Δh_vap)/dT flags), section 3.3 gold standards
  (cyclohexene, adipic acid), 50-molecule group-accounting suite,
  Clausius-Clapeyron finite-difference identity for ΔHvap, and a live
  UManSysProp oracle comparison (26/26 molecules identical counts,
  log10 P agreement to 1e-12, 16 whitelisted known-UMS deviations).

## [0.1.0.dev0] - 2026-11-09

### Added

- Repository skeleton: packaging (`pyproject.toml`, MIT licence, required
  dependencies `thermo>=0.6.1` and `rdkit`), CI (pytest incl. doctests on
  Linux/macOS/Windows × Python 3.9–3.14, strict docs build), and an mkdocs
  scaffold with one page per planned method.
- `aerosolprops.structure`: the RDKit entry layer — cached optional import,
  normalising SMILES parser (charge-form correction for nitrate/PAN),
  element counts, molecular weight.
- `aerosolprops.topology`: molecular-graph primitives as thin wrappers over
  RDKit built-ins (fragmentation on atom/bond deletion, ring arc distance).
