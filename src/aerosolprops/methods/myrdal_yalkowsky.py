# SPDX-License-Identifier: MIT
"""Myrdal-Yalkowsky vapour pressure estimation (transition-temperature input method).

Myrdal & Yalkowsky (*Ind. Eng. Chem. Res.* 36, 2494-2499, 1997,
doi:10.1021/ie950242l) estimate the vapour pressure of an organic
non-electrolyte from two transition temperatures and three structure
level descriptors: the effective number of torsional bonds (flexibility)
:math:`\\tau`, the hydrogen bond number HBN, and the external rotational
symmetry number :math:`\\sigma`. No parameter of the equation was fitted
to vapour-pressure data. The final equation (paper Eq. 12) is

.. math::

    \\log_{10}\\frac{P}{\\mathrm{atm}} =
    -\\frac{[56.5 - 19.2\\log_{10}\\sigma + 9.2\\tau]\\,(T_m - T)}{19.1\\,T}
    -\\frac{[86.0 + 0.4\\tau + 1421\\,\\mathrm{HBN}]\\,(T_b - T)}{19.1\\,T}
    +\\frac{-90.0 - 2.1\\tau}{19.1}
    \\left(\\frac{T_b - T}{T} - \\ln\\frac{T_b}{T}\\right)

with the subsidiary entropies/heat capacity of paper Eqs. 7-11 (see
:func:`myrdal_yalkowsky_entropy_of_melting`,
:func:`myrdal_yalkowsky_entropy_of_boiling`,
:func:`myrdal_yalkowsky_heat_capacity_change_boiling`) and the
descriptors

.. math::

    \\tau = \\max\\{0,\\ \\mathrm{SP3} + 0.5\\,\\mathrm{SP2}
    + 0.5\\,\\mathrm{RING} - 1\\},

    \\mathrm{HBN} = \\frac{\\sqrt{\\mathrm{OH} + \\mathrm{COOH}}
    + 0.33\\sqrt{\\mathrm{NH_2}}}{\\mathrm{MW}} .

Conventions locked against the paper (details and evidence in the
method notes):

* the denominators use the paper's rounded ``2.303 R = 19.1`` J/(mol K),
  not 19.147; Table 1 reproduces digit-for-digit only with 19.1;
* :math:`\\log\\sigma` is base 10 (Yalkowsky-school convention); the
  melting term is the only place it appears;
* sub-cooled liquid branch: with ``Tm=None`` or ``T >= Tm`` the melting
  term is dropped entirely (paper Appendix A: "the first term of eq 12
  is eliminated since the compound is liquid") -- the aerosol-
  thermodynamics quantity; a supplied ``Tm`` below ``T`` activates the
  solid branch (all three terms, solid vapour pressure);
* HBN's first root counts alcohols+phenols and carboxylic acids (paper
  definition), extended -- following Compernolle, Ceulemans & Muller
  (ACP 10, 6271-6282, 2010, section 2.6, doi:10.5194/acp-10-6271-2010)
  -- with hydroperoxides and peracids; the 0.33 root counts primary
  amines (amides are not amines here);
* output is :math:`\\log_{10}(P/\\mathrm{Pa})`: the paper's
  :math:`\\log_{10}(P/\\mathrm{atm})` plus :data:`LOG10_ATM_TO_PA`.

The estimators here are pure functions of the descriptors; the SMILES to
descriptor path (SP3/SP2/RING walk, hydrogen-bond group SMARTS, MW) lives
in :func:`aerosolprops.vapour_pressure.myrdal_yalkowsky_descriptors`.

Examples
--------
Paper Appendix A worked example: 3,5-diethylphenol at 388 K with the
paper's inputs Tb = 521 K, tau = 1.5, HBN = 0.006 66:

>>> round(myrdal_yalkowsky_entropy_of_boiling(1.5, 0.00666), 2)
96.06
>>> round(myrdal_yalkowsky_heat_capacity_change_boiling(1.5), 2)
-93.15
>>> lp = myrdal_yalkowsky_log10_vapour_pressure(388.0, 521.0, 1.5, 0.00666)
>>> round(lp - LOG10_ATM_TO_PA, 2)  # paper: -1.72 - 0.23 = -1.95 (2-dp steps)
-1.96
"""

from math import log, log10, sqrt

from .simpol1 import LOG10_ATM_TO_PA

__all__ = [
    "LOG10_ATM_TO_PA",
    "MY_DCPGL_INTERCEPT",
    "MY_DCPGL_TAU",
    "MY_DENOMINATOR_CONSTANT",
    "MY_NH2_ROOT_COEFFICIENT",
    "MY_RING_WEIGHT",
    "MY_SB_HBN",
    "MY_SB_INTERCEPT",
    "MY_SB_TAU",
    "MY_SM_INTERCEPT",
    "MY_SM_LOG_SIGMA",
    "MY_SM_TAU",
    "MY_SP2_WEIGHT",
    "MY_TAU_FLOOR",
    "MY_TAU_OFFSET",
    "myrdal_yalkowsky_entropy_of_boiling",
    "myrdal_yalkowsky_entropy_of_melting",
    "myrdal_yalkowsky_hbn",
    "myrdal_yalkowsky_heat_capacity_change_boiling",
    "myrdal_yalkowsky_log10_vapour_pressure",
    "myrdal_yalkowsky_tau",
]

# --- coefficients (paper equations; no table, all printed in-line) -----

#: Entropy of melting intercept [J/(mol K)] (paper Eq. 2 and Eq. 10;
#: the 56.5 and -19.2 log(sigma) terms are carried over unchanged from
#: Yalkowsky-Mishra into Eq. 10).
MY_SM_INTERCEPT = 56.5

#: Entropy of melting log10(sigma) coefficient [J/(mol K)] (Eq. 2/10).
MY_SM_LOG_SIGMA = -19.2

#: Entropy of melting tau coefficient [J/(mol K)] (paper Eq. 10; the
#: chain-length term 9.2*(N-3) of Eq. 2 became 9.2*tau).
MY_SM_TAU = 9.2

#: Entropy of boiling intercept [J/(mol K)] (paper Eq. 7, printed 86;
#: Eq. 12 writes 86.0).
MY_SB_INTERCEPT = 86.0

#: Entropy of boiling tau coefficient [J/(mol K)] (paper Eq. 7).
MY_SB_TAU = 0.4

#: Entropy of boiling hydrogen-bond coefficient [J/(mol K) Da]; HBN has
#: units of 1/Da because the paper divides by the molecular weight in
#: g/mol (paper Eq. 7 and Eq. 9).
MY_SB_HBN = 1421.0

#: Effective heat capacity change upon boiling, intercept [J/(mol K)]
#: (paper Eq. 11; back-solved from experimental vapour pressures, an
#: *effective* value; within the Chickos et al. (1993a) range).
MY_DCPGL_INTERCEPT = -90.0

#: Effective heat capacity change upon boiling, tau coefficient
#: [J/(mol K)] (paper Eq. 11).
MY_DCPGL_TAU = -2.1

#: The paper's rounded 2.303 R [J/(mol K)] used in every denominator of
#: Eq. 12 (2.303 * 8.314 = 19.15 rounded by the authors to 19.1; the
#: exact value shifts Table 1 predictions by ~0.02%).
MY_DENOMINATOR_CONSTANT = 19.1

#: Weight of the primary-amine root relative to the alcohol/acid root in
#: HBN (paper Eq. 9).
MY_NH2_ROOT_COEFFICIENT = 0.33

#: Weight of an SP2 atom in the torsional-bond count (paper Eq. 8).
MY_SP2_WEIGHT = 0.5

#: Weight of an independent ring system in the torsional-bond count
#: (paper Eq. 8).
MY_RING_WEIGHT = 0.5

#: The "-1" of the torsional-bond count (paper Eq. 8: one bond fewer
#: than the atom count).
MY_TAU_OFFSET = -1.0

#: The count is floored at zero: the paper sets tau to 0 "when equation
#: 5 is -1 or -0.5" (the equation number is the paper's own misprint
#: for Eq. 8; the semantics are max(0, .)).
MY_TAU_FLOOR = 0.0


def myrdal_yalkowsky_tau(sp3, sp2, ring):
    """Effective number of torsional bonds (paper Eq. 8).

    .. math::

        \\tau = \\max\\{0,\\ \\mathrm{SP3} + 0.5\\,\\mathrm{SP2}
        + 0.5\\,\\mathrm{RING} - 1\\}

    SP3/SP2 are the numbers of nonring, nonterminal sp3 / sp2 atoms and
    RING the number of independent ring systems, as counted by
    :func:`aerosolprops.vapour_pressure.myrdal_yalkowsky_descriptors`.

    Parameters
    ----------
    sp3, sp2, ring : int or float
        Descriptor counts (non-negative; the arguments are not rounded,
        half-integer SP2/RING contributions are legitimate).

    Returns
    -------
    tau : float or None
        The floored count, or ``None`` on invalid input (negative
        counts).

    Examples
    --------
    The paper's own examples (sp3, sp2, ring): propane (1, 0, 0) gives
    the "-1 or -0.5" floor, butane (2, 0, 0) has one torsional bond,
    pentane (3, 0, 0) two, 1-pentene (2, 1, 0) counts the internal
    =CH- at half weight:

    >>> myrdal_yalkowsky_tau(1, 0, 0)
    0.0
    >>> myrdal_yalkowsky_tau(2, 0, 0)
    1.0
    >>> myrdal_yalkowsky_tau(3, 0, 0)
    2.0
    >>> myrdal_yalkowsky_tau(2, 1, 0)  # 1-pentene
    1.5
    >>> myrdal_yalkowsky_tau(0, 0, 1)  # single ring system: floored
    0.0
    >>> myrdal_yalkowsky_tau(-1, 0, 0) is None
    True
    """
    try:
        if sp3 < 0 or sp2 < 0 or ring < 0:
            return None
        tau = (
            sp3 + MY_SP2_WEIGHT * sp2 + MY_RING_WEIGHT * ring + MY_TAU_OFFSET
        )
    except TypeError:
        return None
    return float(max(MY_TAU_FLOOR, tau))


def myrdal_yalkowsky_hbn(n_hbond, n_primary_amine, mw):
    """Hydrogen bond number HBN (paper Eq. 9).

    .. math::

        \\mathrm{HBN} = \\frac{\\sqrt{N_{\\mathrm{HB}}}
        + 0.33\\sqrt{N_{\\mathrm{NH_2}}}}{\\mathrm{MW}}

    The square-root structure (not a linear count) is the pixel-level
    adjudication of the paper's Eq. 9 (the ACS scan OCR layer reads the
    two radical signs as "x"); it is what reproduces Table 1's
    resorcinol prediction (HBN from sqrt(2), not 2).

    Parameters
    ----------
    n_hbond : int
        Number of hydrogen-bonding groups in the first root: alcohols
        (incl. phenols) + carboxylic acids (paper definition), plus
        hydroperoxides and peracids (Compernolle et al., 2010,
        section 2.6 extension adopted by aerosolprops).
    n_primary_amine : int
        Number of primary amines (amides excluded).
    mw : float
        Molecular weight [g/mol] (the paper's dimensional convention;
        the 1421 coefficient of Eq. 7 carries the reciprocal Da).

    Returns
    -------
    hbn : float or None
        HBN [1/Da], or ``None`` on invalid input (negative counts or
        non-positive molecular weight).

    Examples
    --------
    Resorcinol (two phenolic OH, paper MW 110.1) has HBN = sqrt(2)/MW:

    >>> round(myrdal_yalkowsky_hbn(2, 0, 110.1), 6)
    0.012845
    >>> myrdal_yalkowsky_hbn(0, 0, 100.0)
    0.0
    >>> myrdal_yalkowsky_hbn(1, 0, 0.0) is None
    True
    """
    try:
        if n_hbond < 0 or n_primary_amine < 0 or mw <= 0.0:
            return None
        root = sqrt(n_hbond) + MY_NH2_ROOT_COEFFICIENT * sqrt(n_primary_amine)
        return root / mw
    except TypeError:
        return None


def myrdal_yalkowsky_entropy_of_boiling(tau, hbn):
    """Entropy of boiling [J/(mol K)] (paper Eq. 7).

    .. math::

        \\Delta S_b = 86 + 0.4\\tau + 1421\\,\\mathrm{HBN}

    The modified Trouton rule of Myrdal, Krzyzaniak & Yalkowsky
    (IECR 35, 1788, 1996); replaces the chain-length form
    83.7 + 0.7(N-5) of the Yalkowsky-Mishra equation.

    Examples
    --------
    >>> myrdal_yalkowsky_entropy_of_boiling(0.0, 0.0)
    86.0
    >>> round(myrdal_yalkowsky_entropy_of_boiling(1.5, 0.00666), 2)
    96.06
    >>> myrdal_yalkowsky_entropy_of_boiling(-1.0, 0.0) is None
    True
    """
    try:
        if tau < 0.0 or hbn < 0.0:
            return None
        return (
            MY_SB_INTERCEPT + MY_SB_TAU * tau + MY_SB_HBN * hbn
        )
    except TypeError:
        return None


def myrdal_yalkowsky_entropy_of_melting(tau, sigma=1.0):
    """Entropy of melting [J/(mol K)] (paper Eq. 10).

    .. math::

        \\Delta S_m = 56.5 - 19.2\\log_{10}\\sigma + 9.2\\tau

    The chain-length term 9.2*(N-3) of Yalkowsky-Mishra (Eq. 2) with the
    flexibility descriptor in its place; the 56.5 and -19.2 log10(sigma)
    terms are carried over unchanged. ``log`` is base 10 (Yalkowsky
    school convention, cf. Dannenfelser et al., 1993); Table 1's only
    sigma != 1 row (biphenyl, sigma = 2) is insensitive to the base at
    its precision, so the base is documented rather than testable
    there.

    Parameters
    ----------
    tau : float
        Effective number of torsional bonds (>= 0).
    sigma : float, optional
        External rotational symmetry number (>= 1; naphthalene 4,
        benzene 12; 1 for unsymmetric molecules - the SOA default).

    Examples
    --------
    >>> myrdal_yalkowsky_entropy_of_melting(0.0, 1.0)
    56.5
    >>> round(myrdal_yalkowsky_entropy_of_melting(0.0, 12.0), 2)  # benzene
    35.78
    >>> myrdal_yalkowsky_entropy_of_melting(0.0, 0.5) is None
    True
    """
    try:
        if tau < 0.0 or sigma < 1.0:
            return None
        return (
            MY_SM_INTERCEPT + MY_SM_LOG_SIGMA * log10(sigma) + MY_SM_TAU * tau
        )
    except (TypeError, ValueError):
        return None


def myrdal_yalkowsky_heat_capacity_change_boiling(tau):
    """Effective heat capacity change upon boiling [J/(mol K)] (Eq. 11).

    .. math::

        \\Delta C_p^{(gl)} = -90.0 - 2.1\\tau

    Back-solved by the authors from experimental vapour pressures
    through Eq. 1 (an *effective* value that absorbs model error; the
    paper discusses it against the Chickos et al. (1993a) means of
    -64.2 J/(mol K) for liquids and -134 J/(mol K) for solids).

    Examples
    --------
    >>> myrdal_yalkowsky_heat_capacity_change_boiling(0.0)
    -90.0
    >>> round(myrdal_yalkowsky_heat_capacity_change_boiling(1.5), 2)
    -93.15
    >>> myrdal_yalkowsky_heat_capacity_change_boiling(-0.5) is None
    True
    """
    try:
        if tau < 0.0:
            return None
        return MY_DCPGL_INTERCEPT + MY_DCPGL_TAU * tau
    except TypeError:
        return None


def myrdal_yalkowsky_log10_vapour_pressure(
    T, Tb, tau=0.0, hbn=0.0, Tm=None, sigma=1.0
):
    """Estimate log10 of the vapour pressure [log10 Pa] (paper Eq. 12).

    Parameters
    ----------
    T : float
        Temperature of interest [K] (> 0). The equation is validated for
    Tb : float
        Normal boiling point [K] (> 0); an input, not an estimate of
        this module (route it from
        :class:`aerosolprops.boiling_point.NannoolalBP`, ``thermo``'s
        Joback method, or experiment through the
        :class:`aerosolprops.vapour_pressure.MyrdalYalkowsky` wrapper).
    tau : float, optional
        Effective number of torsional bonds (>= 0,
        :func:`myrdal_yalkowsky_tau`).
    hbn : float, optional
        Hydrogen bond number (>= 0, :func:`myrdal_yalkowsky_hbn`).
    Tm : float, optional
        Melting point [K]. ``None`` (default) selects the sub-cooled
        liquid branch (melting term dropped, the aerosol-thermodynamics
        quantity); a value with ``T < Tm`` selects the solid branch
        (all three terms, the solid vapour pressure); a value with
        ``T >= Tm`` is the liquid branch of the paper's Appendix A.
    sigma : float, optional
        External rotational symmetry number (>= 1); only enters the
        melting term, i.e. only the solid branch.

    Returns
    -------
    log10P : float or None
        Base-10 logarithm of the vapour pressure [Pa], or ``None`` on
        invalid input.

    Examples
    --------
    Paper Table 1 rows at the paper's own precision (predicted column;
    HBN built through :func:`myrdal_yalkowsky_hbn` with the paper's
    molecular weights). Resorcinol at 425 K (Tb 549 K, tau 0, two
    phenolic OH on MW 110.1):

    >>> hbn = myrdal_yalkowsky_hbn(2, 0, 110.1)
    >>> lp = myrdal_yalkowsky_log10_vapour_pressure(425.0, 549.0, 0.0, hbn)
    >>> round(lp - LOG10_ATM_TO_PA, 2)
    -1.76

    Octacosane at 408 K (Tb 705 K, tau 25):

    >>> lp = myrdal_yalkowsky_log10_vapour_pressure(408.0, 705.0, 25.0, 0.0)
    >>> round(lp - LOG10_ATM_TO_PA, 2)
    -5.01

    The solid branch (biphenyl at 342 K, just below its 342.4 K melting
    point, sigma = 2; Table 1 predicts -2.97):

    >>> lp = myrdal_yalkowsky_log10_vapour_pressure(
    ...     342.0, 528.0, 0.0, 0.0, Tm=342.4, sigma=2.0)
    >>> round(lp - LOG10_ATM_TO_PA, 2)
    -2.97

    Without the melting point the same call is the sub-cooled liquid
    (melting term absent -- here only -0.003 log units below the solid):

    >>> lp = myrdal_yalkowsky_log10_vapour_pressure(342.0, 528.0, 0.0, 0.0)
    >>> round(lp - LOG10_ATM_TO_PA, 2)
    -2.97

    Invalid input returns ``None`` rather than raising:

    >>> myrdal_yalkowsky_log10_vapour_pressure(0.0, 373.0) is None
    True
    """
    try:
        if T <= 0.0 or Tb <= 0.0 or (Tm is not None and Tm <= 0.0):
            return None
        if tau < 0.0 or hbn < 0.0 or sigma < 1.0:
            return None
        denom = MY_DENOMINATOR_CONSTANT
        logp = 0.0
        if Tm is not None and T < Tm:
            sm = myrdal_yalkowsky_entropy_of_melting(tau, sigma)
            if sm is None:
                return None
            logp -= sm * (Tm - T) / (denom * T)
        sb = myrdal_yalkowsky_entropy_of_boiling(tau, hbn)
        dcp = myrdal_yalkowsky_heat_capacity_change_boiling(tau)
        if sb is None or dcp is None:
            return None
        logp -= sb * (Tb - T) / (denom * T)
        logp += (dcp / denom) * ((Tb - T) / T - log(Tb / T))
    except TypeError:
        return None
    return logp + LOG10_ATM_TO_PA
