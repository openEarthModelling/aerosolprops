# Changelog

All notable changes to this project are documented in this file.
The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/)
and the project adheres to [semantic versioning](https://semver.org/).

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
