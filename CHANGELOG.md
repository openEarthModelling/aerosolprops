# Changelog

All notable changes to this project are documented in this file.
The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/)
and the project adheres to [semantic versioning](https://semver.org/).

## [Unreleased]

### Added

- EVAPORATION vapour pressure (Compernolle, Ceulemans & Müller, ACP 11,
  9431, 2011): the third estimator method.
  - `aeroprops.methods.evaporation` (pure layer): the 40+3 coefficients
    of paper Table 4 (27 count classes incl. the 12lin/12CL/12HB ring
    split, intercepts, diacid factor) cross-verified against the paper's
    independent g_k column, the A/B parameter and Eq. 14 boiling-point
    estimators, and `evaporation_log10_vapour_pressure` (log10 Pa basis).
    Transcription rulings recorded in tests: b20 = 1817.69 (paper) vs
    1817.60 in the UManSysProp tree (typo there).
  - `aeroprops.vapour_pressure` (wrapper layer):
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
  - `aeroprops.methods.nannoolal_bp` (pure layer): the 132-group
    coefficient catalogue (113 main + 16 corrections + 3 extensions;
    the two Table 1 transcription rulings are documented in tests),
    the 76 group-interaction parameters of Table 5 with the
    counting-form GI estimator `nannoolal_bp_gi` (parameter-less pairs
    B-G and R-R fail; extension groups switch to a flagged no-GI
    branch), and `nannoolal_bp_boiling_point`.
  - `aeroprops.fragmentation`: reusable priority-number fragmenter
    (paper Table 1 footnote a semantics — lower priority number wins,
    group-level claiming, multi-atom all-or-nothing claims, atom
    balance check) for any catalogue of (SMARTS, priority, claim).
  - `aeroprops.boiling_point`: clean-room SMARTS catalogue plus graph
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
  - `aeroprops.methods.simpol1` (pure layer): the 124 coefficients of
    paper Table 5 as an in-source catalogue, the SMARTS catalogue with
    algorithmic-group placeholders, and pure estimators
    `simpol1_bk`, `simpol1_log10_vapour_pressure` (log10 Pa basis),
    and `simpol1_vapourisation_enthalpy` (J/mol).
  - `aeroprops.vapour_pressure` (wrapper layer): `simpol1_fragmentation`
    and the Joback-style `SIMPOL1` class with one-time fragmentation and
    `None`-on-failure estimators. Fragmentation implements the paper's
    overlapping-additive counting semantics (O-claims of carbonyl
    classes, per-OH nitrophenol replacement of phenol, nitroester
    replacement of ester, amide acid-side graph cut, SSSR ring counts).
  - `aeroprops.structure.normalize_mol`: public normalisation entry
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
- `aeroprops.structure`: the RDKit entry layer — cached optional import,
  normalising SMILES parser (charge-form correction for nitrate/PAN),
  element counts, molecular weight.
- `aeroprops.topology`: molecular-graph primitives as thin wrappers over
  RDKit built-ins (fragmentation on atom/bond deletion, ring arc distance).
