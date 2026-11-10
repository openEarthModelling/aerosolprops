# SPDX-License-Identifier: MIT
"""Pure estimator layer, one submodule per property method.

Modules in this package contain only coefficient catalogues, SMARTS
catalogues, and estimators written as pure functions of group-count
dictionaries. They import neither RDKit nor ``thermo``: the structure
engine lives in :mod:`aeroprops.structure`, the graph primitives in
:mod:`aeroprops.topology`, and the object-oriented wrappers in the
property modules (e.g. :mod:`aeroprops.vapour_pressure`).
"""
