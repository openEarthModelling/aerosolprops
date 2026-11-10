# SPDX-License-Identifier: MIT
"""EVAPORATION estimation of pure-component vapour pressure.

EVAPORATION (Compernolle, Ceulemans & Müller, Atmos. Chem. Phys. 11,
9431-9450, 2011, doi:10.5194/acp-11-9431-2011) is a group-contribution
method parameterised for multifunctional atmospheric oxidation
products. It needs *only* the molecular structure and the temperature -
no boiling point input (the normal boiling point is an output, paper
Eq. 14). The sub-cooled liquid vapour pressure is

.. math::

    \\log_{10}(p^0/\\mathrm{atm}) = A + B/T^{\\kappa}, \\qquad
    \\kappa = 1.5

(paper Eq. 13) with the A term split into additive ("lin"),
carbonyl-like ("CL") and hydrogen-bonding ("HB") parts divided by
:math:`\\sqrt{N}` (paper Eqs. 27-30, :math:`r = 0.5`), and a
diacid correction scaling all CL/HB counts by
:math:`N_{\\mathrm{eff}}/(N_{\\mathrm{CL}}+N_{\\mathrm{HB}})` with
:math:`N_{\\mathrm{eff}} = 2.6` when the molecule has more than one
acid and more than two CL/HB groups in total (paper Eq. 31, applied
to both A and B).

The paper fitted :math:`\\log_{10}(p^0/\\mathrm{atm})`; this module
returns :math:`\\log_{10}(P/\\mathrm{Pa})` (strict SI; the conversion
constant is added by :func:`evaporation_log10_vapour_pressure`).

The 40 coefficients below are transcribed from Table 4 of the paper
and cross-verified against the g_k column of the same table (paper
Eq. 22 at 298.15 K; 18 of 20 rows reproduce the printed value to
within 0.005 - rows k5 (-1.194 recomputed vs -1.20 printed) and k9
(-3.6123 vs -3.62) sit 0.006/0.008 away, a paper-side rounding gap
consistent with g_k having been computed from the unrounded fitted
parameters). Descriptor counts (the ``counts``
dictionaries consumed here, keys "1".."20" plus the typed ring
splits "12lin"/"12CL"/"12HB") are produced by the graph fragmenter in
:mod:`aeroprops.vapour_pressure`.

Two assembly modes
------------------
``upstream_compat=False`` (default) follows the paper equations
directly: the sqrt/non-additivity form of Eqs. 29-30 is *continuous*
(it is applied to every molecule, :math:`\\sqrt{1} = 1`), the k12
ring counts enter the B parameter typed (12lin unscaled, 12CL/12HB
diacid-scaled) per the Table 4 typing, and b20 = 1817.69 as printed.

``upstream_compat=True`` reproduces the arithmetic of the
UManSysProp implementation bit-for-bit (verified against the live
oracle on 75 molecules): the sqrt branch is entered only when
N_CL > 1 or N_HB > 1 (otherwise a flat sum is used, which differs
from the continuous form only in float rounding), the *entire* k12
count is treated as CL/HB in the B parameter, and b20 = 1817.60
(the transcription typo of the upstream data file). The counting
semantics that feed this mode are likewise replicated; see
:func:`aeroprops.vapour_pressure.evaporation_fragmentation`.

References
----------
.. [1] Compernolle, S., Ceulemans, K., Müller, J.-F. "EVAPORATION: a
   new vapour pressure estimation method for organic molecules
   including non-additivity and intramolecular interactions",
   Atmospheric Chemistry and Physics 11, 9431-9450 (2011).
   doi:10.5194/acp-11-9431-2011
"""

from math import log10, sqrt

__all__ = [
    "LOG10_ATM_TO_PA",
    "EVAPORATION_KAPPA",
    "EVAPORATION_R",
    "EVAPORATION_N_EFF",
    "EVAPORATION_B20_UPSTREAM_TYPO",
    "EVAPORATION_COEFFICIENTS",
    "EVAPORATION_GK_PAPER",
    "evaporation_gk",
    "evaporation_diacid_factor",
    "evaporation_log10_vapour_pressure",
]

#: log10 of the Pascal content of one atmosphere; converts the paper's
#: fitted log10(p/atm) basis to the aeroprops log10(P/Pa) basis.
LOG10_ATM_TO_PA = log10(101325.0)

#: Temperature exponent kappa of paper Eq. 13 (optimised in paper
#: Sect. 4.1; implies dHvap proportional to T^0.5).
EVAPORATION_KAPPA = 1.5

#: Non-additivity exponent r of paper Eqs. 29-30 (A_CL and A_HB are
#: divided by sqrt(N)).
EVAPORATION_R = 0.5

#: Effective functional-group number of paper Eq. 31 for functionalised
#: diacids (N_eff / (N_CL + N_HB) scales every CL/HB count).
EVAPORATION_N_EFF = 2.6

#: The upstream UManSysProp data file mistypes the b coefficient of
#: descriptor k20 as 1817.60 (paper Table 4: 1817.69);
#: ``upstream_compat=True`` reproduces the typo.
EVAPORATION_B20_UPSTREAM_TYPO = 1817.60

#: Coefficients of paper Table 4, keyed by the paper's descriptor index
#: k (as a string, the key used by the counts dictionaries). Each record
#: is ``(a_k, b_k, type, name)``; k1 is the zero-point constant with
#: c1 = 1 for every compound. k12 ("X on ring") is typed by the
#: functional group it modifies, so its count arrives split across the
#: "12lin"/"12CL"/"12HB" keys.
EVAPORATION_COEFFICIENTS = {
    "1": (2.6255, -1986.56, "lin", "zero point"),
    "2": (0.06298, -2821.46, "lin", "#C + #in-chain O"),
    "3": (-0.00293, 1040.69, "lin", "t = branching - ring number"),
    "4": (0.71114, -15841.13, "lin", "nitrate ester"),
    "5": (0.19747, -7163.72, "CL", "ketone/aldehyde carbonyl"),
    "6": (0.32257, -5208.53, "CL", "ester"),
    "7": (0.29030, -15011.33, "CL", "peroxy acyl nitrate (PAN)"),
    "8": (0.95537, -16699.73, "HB", "hydroxyl"),
    "9": (0.98567, -23671.00, "HB", "carboxylic acid"),
    "10": (0.78348, -18583.48, "HB", "hydroperoxide"),
    "11": (0.81498, -18071.50, "HB", "peroxy acid"),
    "12": (0.18704, -2509.37, "lin,CL,HB", "X on ring (typed)"),
    "13": (-0.18596, 14.21, "CL", "C=C conjugated with carbonyl"),
    "14": (-0.28012, 4201.34, "HB", "0/1/2 for prim/sec/tert OH"),
    "15": (-0.34191, 2961.95, "HB", "alkenoic (allylic) alcohol flag"),
    "16": (0.26830, 1602.62, "CL", "per carbonyl: CL C=O at alpha"),
    "17": (0.11716, 939.71, "CL", "per carbonyl: CL C=O at beta"),
    "18": (-0.30373, 3769.87, "CL", "per carbonyl: non-CL non-acid group at alpha"),
    "19": (-0.04143, 800.12, "HB", "per hydroxyl: functional group at alpha"),
    "20": (0.46023, 1817.69, "HB", "per acid: CL C=O at alpha"),
}

#: The g_k column of paper Table 4 (group contribution at 298 K, paper
#: Eq. 22), used as an independent transcription cross-check by
#: :func:`evaporation_gk` and the test suite. Printed precision varies
#: by row; k5 and k9 are 0.006/0.008 away from a_k + b_k/298.15^1.5
#: (the paper's own rounding), every other row is within 0.005.
EVAPORATION_GK_PAPER = {
    "1": 2.24, "2": -0.485, "3": 0.20, "4": -2.37, "5": -1.20,
    "6": -0.69, "7": -2.63, "8": -2.29, "9": -3.62, "10": -2.83,
    "11": -2.70, "12": -0.30, "13": -0.18, "14": 0.54, "15": 0.23,
    "16": 0.58, "17": 0.30, "18": 0.43, "19": 0.11, "20": 0.81,
}

#: Descriptor keys of each typed part (paper Table 4, footnote d): the
#: additive part (k1-k4 plus k12-on-lin-groups), the carbonyl-like part
#: (k5-7, 13, 16-18 plus k12-on-CL-groups) and the hydrogen-bonding
#: part (k8-11, 14, 15, 19, 20 plus k12-on-HB-groups). N_CL = c5+c6+c7
#: and N_HB = c8+c9+c10+c11 (paper Sect. 4.3.1); the second-order k16
#: to k20 corrections share the type of the group they modify.
_LIN_KEYS = ("1", "2", "3", "4", "12lin")
_CL_KEYS = ("5", "6", "7", "13", "16", "17", "18", "12CL")
_HB_KEYS = ("8", "9", "10", "11", "14", "15", "19", "20", "12HB")


def evaporation_gk(k, T=298.15):
    """Group contribution g_k = a_k + b_k/T^kappa at temperature `T`.

    Paper Eq. 22 tabulated at 298 K in the g_k column of Table 4.

    Parameters
    ----------
    k : str or int
        Descriptor index "1".."20".
    T : float, optional
        Temperature [K].

    Returns
    -------
    gk : float or None
        ``None`` for an unknown descriptor or non-positive `T`.

    Examples
    --------
    Rows of Table 4 recomputed at 298.15 K (printed values in
    parentheses). k3, the branching/ring topology term ("each extra
    branch raises log10 p0 by ~0.2"):

    >>> round(evaporation_gk("3"), 2)
    0.2
    >>> round(evaporation_gk("9"), 4)   # strongest HB term; prints -3.62
    -3.6123
    >>> round(evaporation_gk(4), 2)     # int keys are accepted
    -2.37
    >>> evaporation_gk("21") is None
    True
    """
    try:
        a_k, b_k = EVAPORATION_COEFFICIENTS[str(k)][:2]
        return a_k + b_k / T ** EVAPORATION_KAPPA
    except (KeyError, TypeError, ZeroDivisionError):
        return None


