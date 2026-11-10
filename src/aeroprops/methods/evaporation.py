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

from math import isfinite, log10, sqrt

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


def _typed_maps(b20):
    """(a, b) coefficient maps including the typed k12 duplicate rows."""
    a = {k: rec[0] for k, rec in EVAPORATION_COEFFICIENTS.items()}
    b = {k: rec[1] for k, rec in EVAPORATION_COEFFICIENTS.items()}
    for split in ("12lin", "12CL", "12HB"):
        a[split] = a["12"]
        b[split] = b["12"]
    b["20"] = b20
    return a, b


#: Coefficient maps with k12 duplicated into its typed keys; the B map
#: carries the paper's b20 = 1817.69.
_A_PAPER, _B_PAPER = _typed_maps(EVAPORATION_COEFFICIENTS["20"][1])

#: The upstream-compatibility B map reproduces the b20 = 1817.60 typo
#: of the UManSysProp data file (A is unaffected by the typo).
_A_UPSTREAM, _B_UPSTREAM = _typed_maps(EVAPORATION_B20_UPSTREAM_TYPO)

#: Every accepted key of a counts dictionary ("1" is forced to 1 and
#: never read from the input).
_COUNT_KEYS = tuple(
    [str(k) for k in range(1, 21)] + ["12lin", "12CL", "12HB"]
)


def _validated_counts(counts):
    """Normalise a counts mapping to float values, or ``None``.

    Accepts the descriptor keys "1".."20" plus the typed ring splits;
    unknown keys (typos) and non-numeric or non-finite values reject
    the whole input, mirroring the unknown-group discipline of
    :func:`aeroprops.methods.simpol1.simpol1_log10_vapour_pressure`.
    """
    try:
        out = {}
        for key, value in counts.items():
            if key not in _COUNT_KEYS:
                return None
            out[key] = float(value)
    except (AttributeError, TypeError, ValueError):
        return None
    if any(not isfinite(v) for v in out.values()):
        return None
    out["1"] = 1.0
    for key in _COUNT_KEYS:
        out.setdefault(key, 0.0)
    return out


def evaporation_diacid_factor(counts):
    """Diacid correction factor of paper Eq. 31.

    Functionalised diacids - more than one carboxylic acid (c9 > 1)
    *and* more than two CL/HB groups in total - have their CL/HB terms
    scaled by ``N_eff / (N_CL + N_HB)`` with ``N_eff = 2.6`` (paper
    Sect. 4.3.2); every other molecule gets the neutral factor 1.0.

    Parameters
    ----------
    counts : dict
        Descriptor counts (see
        :func:`evaporation_log10_vapour_pressure`). Only the
        functional-group keys "5"-"11" are read.

    Returns
    -------
    factor : float or None
        ``None`` on invalid input.

    Examples
    --------
    >>> round(evaporation_diacid_factor({"8": 1, "9": 3}), 4)  # citric
    0.65
    >>> round(evaporation_diacid_factor({"9": 3}), 6)  # MBTCA: 2.6/3
    0.866667
    >>> evaporation_diacid_factor({"9": 2})  # plain diacid: neutral
    1.0
    >>> evaporation_diacid_factor({"9": 1, "8": 2})  # monoacid: neutral
    1.0
    """
    c = _validated_counts(counts)
    if c is None:
        return None
    n_cl = c["5"] + c["6"] + c["7"]
    n_hb = c["8"] + c["9"] + c["10"] + c["11"]
    if c["9"] > 1 and n_cl + n_hb > 2:
        return EVAPORATION_N_EFF / (n_cl + n_hb)
    return 1.0


def _aggregate(c, coeffs, keys):
    """Sum(c[k] * coeffs[k]) over `keys` in order (float-exact role of
    the upstream ``aggregate_matches`` helper, same left-to-right sum).
    """
    return sum(c[k] * coeffs[k] for k in keys)


def evaporation_log10_vapour_pressure(counts, T=298.15,
                                      upstream_compat=False):
    """Estimate log10 of the sub-cooled liquid vapour pressure [log10 Pa].

    Assembles paper Eq. 13 from the 20 descriptor counts plus the typed
    ring counts (``upstream_compat=False``) or reproduces the
    UManSysProp assembly bit-for-bit (``upstream_compat=True``); see
    the module docstring for the exact differences.

    Parameters
    ----------
    counts : dict
        Descriptor counts keyed "1".."20" (c1 is forced to 1) plus
        "12lin"/"12CL"/"12HB", as produced by
        :func:`aeroprops.vapour_pressure.evaporation_fragmentation`.
        Counts may be fractional (the topology descriptor c3 carries a
        fractional ring number).
    T : float, optional
        Temperature [K].
    upstream_compat : bool, optional
        Select the upstream assembly arithmetic (see module docstring).

    Returns
    -------
    log10P : float or None
        Base-10 logarithm of the vapour pressure [Pa], or ``None`` on
        invalid input (unknown descriptor key, non-numeric count,
        non-positive temperature).

    Examples
    --------
    Neopentane, paper counting (5 carbons, t = 2; the branching
    descriptor alone raises log10 p0 by 2 x 0.199 over n-pentane):

    >>> lp = evaporation_log10_vapour_pressure({"2": 5, "3": 2.0})
    >>> round(lp, 4)
    5.2184
    >>> round(lp - LOG10_ATM_TO_PA, 4)   # paper basis, log10 atm
    0.2127

    MBTCA (3-methyl-1,2,3-butanetricarboxylic acid), paper counting
    (t = 3, N_HB = 3, diacid factor 2.6/3):

    >>> lp = evaporation_log10_vapour_pressure({"2": 8, "3": 3.0, "9": 3})
    >>> round(lp - LOG10_ATM_TO_PA, 4)
    -11.5184

    The same molecule under upstream counting (t = 7 through the
    subset-matching branching patterns) and upstream arithmetic; the
    value equals the live UManSysProp oracle bit-for-bit (paper Table 9
    lists -10.86, between the two c3 conventions):

    >>> lp = evaporation_log10_vapour_pressure(
    ...     {"2": 8, "3": 7.0, "9": 3}, upstream_compat=True)
    >>> round(lp, 4)
    -5.7158
    >>> round(lp - LOG10_ATM_TO_PA, 4)
    -10.7215

    Invalid input returns ``None`` rather than raising:

    >>> evaporation_log10_vapour_pressure({"2": 1}, T=0.0) is None
    True
    >>> evaporation_log10_vapour_pressure({"foo": 1}) is None
    True
    """
    try:
        temperature = float(T)
    except (TypeError, ValueError):
        return None
    if not temperature > 0.0:
        return None
    c = _validated_counts(counts)
    if c is None:
        return None

    n_cl = c["5"] + c["6"] + c["7"]
    n_hb = c["8"] + c["9"] + c["10"] + c["11"]
    if not upstream_compat:
        # Paper Eqs. 27-31 in their continuous form: sqrt(1) = 1, and
        # the B parameter takes k12 typed (12lin additive, 12CL/12HB
        # diacid-scaled, per the Table 4 typing of the group it
        # modifies).
        factor = evaporation_diacid_factor(c)
        a_lin = _aggregate(c, _A_PAPER, _LIN_KEYS)
        a_cl = _aggregate(c, _A_PAPER, _CL_KEYS) / sqrt(n_cl) \
            if n_cl else 0.0
        a_hb = _aggregate(c, _A_PAPER, _HB_KEYS) / sqrt(n_hb) \
            if n_hb else 0.0
        A = a_lin + (a_cl + a_hb) * factor
        b_cl_hb = _aggregate(c, _B_PAPER, tuple(
            str(k) for k in range(5, 21))) \
            + _B_PAPER["12"] * (c["12CL"] + c["12HB"])
        B = _aggregate(c, _B_PAPER, ("1", "2", "3", "4", "12lin")) \
            + b_cl_hb * factor
    elif n_cl > 1 or n_hb > 1:
        # Upstream sqrt branch: identical trigger, grouping and
        # summation order to the UManSysProp source (verified
        # bit-for-bit against the live oracle); the B parameter scales
        # the *total* k12 count and uses the b20 typo.
        factor = evaporation_diacid_factor(c)
        a_lin = _aggregate(c, _A_UPSTREAM, _LIN_KEYS)
        try:
            a_cl = _aggregate(c, _A_UPSTREAM, _CL_KEYS) / sqrt(n_cl)
        except ZeroDivisionError:
            a_cl = 0.0
        try:
            a_hb = _aggregate(c, _A_UPSTREAM, _HB_KEYS) / sqrt(n_hb)
        except ZeroDivisionError:
            a_hb = 0.0
        A = a_lin + (a_cl + a_hb) * factor
        B = _aggregate(c, _B_UPSTREAM, ("1", "2", "3", "4")) + (
            _aggregate(c, _B_UPSTREAM, tuple(str(k) for k in range(5, 21)))
            * factor
        )
    else:
        keys = tuple(str(k) for k in range(1, 21))
        A = _aggregate(c, _A_UPSTREAM, keys)
        B = _aggregate(c, _B_UPSTREAM, keys)
    return A + B / temperature ** EVAPORATION_KAPPA + LOG10_ATM_TO_PA




