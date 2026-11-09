# SPDX-License-Identifier: MIT
"""Structure-derived pure-component properties for atmospheric organic
aerosol.

Successor to the pure-component property methods of UManSysProp
(Topping et al., Geosci. Model Dev. 9, 899-911, 2016), rebuilt on RDKit as
the structure engine and following the conventions of the ``thermo``
package family: strict SI units (Pa, K, mol), estimators as pure functions
of group counts where possible, and failures returned as ``None``.

Planned method coverage for v0.1:

- SIMPOL.1 vapour pressure (Pankow & Asher, ACP 8, 2773, 2008)
- EVAPORATION vapour pressure (Compernolle et al., ACP 11, 9431, 2011)
- Nannoolal Part 1 normal boiling point (Nannoolal et al., FPE 226, 45,
  2004), including the SOA extension groups of Compernolle et al. (ACP 10,
  6271, 2010)
- Nannoolal Part 3 vapour pressure (Nannoolal et al., FPE 269, 117, 2008)
- Myrdal-Yalkowsky vapour pressure (Myrdal & Yalkowsky, IECR 36, 2494,
  1997) with the hydrogen-bonding extension of Compernolle et al. (2010)

Boiling-point inputs may instead come from ``thermo``'s built-in Joback
method or from experimental data; aeroprops does not reimplement Joback.

RDKit is an optional dependency: importing this package pulls in neither
RDKit nor ``thermo``'s group-contribution machinery, and every estimator
consumes plain group-count dictionaries on its no-RDKit face.
"""

__version__ = "0.1.0.dev0"

from . import structure, topology

__all__ = ["__version__", "structure", "topology"]
