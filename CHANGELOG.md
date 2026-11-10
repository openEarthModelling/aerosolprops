# Changelog

All notable changes to this project are documented in this file.
The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/)
and the project adheres to [semantic versioning](https://semver.org/).

## [Unreleased]

### Added

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
