# SPDX-License-Identifier: MIT
"""SIMPOL.1 group-contribution estimation of pure-component vapour pressure.

SIMPOL.1 (Pankow & Asher, Atmos. Chem. Phys. 8, 2773-2796, 2008,
doi:10.5194/acp-8-2773-2008) estimates the sub-cooled liquid vapour
pressure of a multifunctional organic compound as a linear function of
31 group counts:

.. math::

    \\log_{10}(P/\\mathrm{Pa}) = \\sum_{k=0}^{30} \\nu_k\\, b_k(T)
    + \\log_{10}(101325)

with the temperature-dependent group contributions (paper Eq. 8)

.. math::

    b_k(T) = B_{k,1}/T + B_{k,2} + B_{k,3} T + B_{k,4} \\ln T .

:math:`\\nu_0 \\equiv 1` is a constant "zeroth group". The paper fitted
:math:`\\log_{10}(p/\\mathrm{atm})`; this module returns
:math:`\\log_{10}(P/\\mathrm{Pa})` (strict SI, the conversion constant is
added by :func:`simpol1_log10_vapour_pressure`). No boiling or melting
point is required.

The 124 coefficients below are transcribed from Table 5 of the paper and
were cross-verified against the independent b_k(293.15 K) column of its
Table 6 (all 31 groups reproduce the tabulated 3 significant figures).
Group counts (the ``counts`` dictionaries consumed here) are produced by
the fragmenter in :mod:`aerosolprops.vapour_pressure`; group indices k=0..30
follow the paper's Table 5 exactly.

Examples
--------
Literature values from the worked examples of the paper (section 3.3,
T = 293.15 K). Cyclohexene, 6 carbons / 1 alicyclic ring / 1 C=C:

>>> nu = {1: 6, 4: 1, 5: 1}
>>> lp = simpol1_log10_vapour_pressure(nu, T=293.15)
>>> round(lp, 4)
4.0567
>>> round(lp - LOG10_ATM_TO_PA, 2)  # paper reports -0.94 (log10 atm)
-0.95

Adipic acid, 6 carbons / 2 carboxylic acids (paper: -7.99):

>>> nu = {1: 6, 10: 2}
>>> round(simpol1_log10_vapour_pressure(nu, T=293.15) - LOG10_ATM_TO_PA, 2)
-7.99

Group terms b_k at 293.15 K against Table 6 of the paper:

>>> round(simpol1_bk(0, 293.15), 2)
1.79
>>> round(simpol1_bk(1, 293.15), 3)
-0.438
>>> round(simpol1_bk(10, 293.15), 2)
-3.58
>>> round(simpol1_bk(29, 293.15), 4)
0.0432

References
----------
.. [1] Pankow, J. F., Asher, W. E. "SIMPOL.1: a simple group contribution
   method for predicting vapor pressures and enthalpies of vaporization of
   multifunctional organic compounds." Atmospheric Chemistry and Physics 8,
   2773-2796 (2008). doi:10.5194/acp-8-2773-2008
"""

from math import log, log10

__all__ = [
    "LOG10_ATM_TO_PA",
    "SIMPOL1_GROUPS",
    "SIMPOL1_SMARTS",
    "SIMPOL1_ALGORITHMIC_GROUPS",
    "SIMPOL1_TEMPERATURE_RANGE",
    "SIMPOL1_GAS_CONSTANT",
    "simpol1_bk",
    "simpol1_log10_vapour_pressure",
    "simpol1_vapourisation_enthalpy",
]

#: log10 of the Pascal content of one atmosphere; converts the paper's
#: fitted log10(p/atm) basis to the aerosolprops log10(P/Pa) basis.
LOG10_ATM_TO_PA = log10(101325.0)

#: Gas constant used by the paper (J/mol/K) in the Clausius-Clapeyron
#: derivation of the vaporisation enthalpy (its Eq. 9-11).
SIMPOL1_GAS_CONSTANT = 8.314

#: Temperature grid of the fit (K); outside this range, and especially
#: below 300 K for low-volatility compounds, errors grow markedly
#: (paper Figs. 10-12).
SIMPOL1_TEMPERATURE_RANGE = (273.15, 393.15)

#: Coefficients of paper Table 5, keyed by the paper's group index k.
#: Each record is ``(name, B_k1, B_k2, B_k3, B_k4)``; group 0 is the
#: constant "zeroth group" with nu_0 = 1 for every compound.
SIMPOL1_GROUPS = {
    0: ("zeroth group (constant term)", -4.26938e+02, 2.89223e-01, 4.42057e-03, 2.92846e-01),
    1: ("carbon number", -4.11248e+02, 8.96919e-01, -2.48607e-03, 1.40312e-01),
    2: ("carbon number, acid-side of amide", -1.46442e+02, 1.54528e+00, 1.71021e-03, -2.78291e-01),
    3: ("aromatic ring", 3.50262e+01, -9.20839e-01, 2.24399e-03, -9.36300e-02),
    4: ("non-aromatic ring", -8.72770e+01, 1.78059e+00, -3.07187e-03, -1.04341e-01),
    5: ("C=C (non-aromatic)", 5.73335e+00, 1.69764e-02, -6.28957e-04, 7.55434e-03),
    6: ("C=C-C=O in non-aromatic ring", -2.61268e+02, -7.63282e-01, -1.68213e-03, 2.89038e-01),
    7: ("hydroxyl (alkyl)", -7.25373e+02, 8.26326e-01, 2.50957e-03, -2.32304e-01),
    8: ("aldehyde", -7.29501e+02, 9.86017e-01, -2.92664e-03, 1.78077e-01),
    9: ("ketone", -1.37456e+01, 5.23486e-01, 5.50298e-04, -2.76950e-01),
    10: ("carboxylic acid", -7.98796e+02, -1.09436e+00, 5.24132e-03, -2.28040e-01),
    11: ("ester", -3.93345e+02, -9.51778e-01, -2.19071e-03, 3.05843e-01),
    12: ("ether (open-chain dialkyl)", -1.44334e+02, -1.85617e+00, -2.37491e-05, 2.88290e-01),
    13: ("ether (alicyclic)", 4.05265e+01, -2.43780e+00, 3.60133e-03, 9.86422e-02),
    14: ("ether, aromatic", -7.07406e+01, -1.06674e+00, 3.73104e-03, -1.44003e-01),
    15: ("nitrate", -7.83648e+02, -1.03439e+00, -1.07148e-03, 3.15535e-01),
    16: ("nitro", -5.63872e+02, -7.18416e-01, 2.63016e-03, -4.99470e-02),
    17: ("aromatic hydroxyl (e.g. phenol)", -4.53961e+02, -3.26105e-01, -1.39780e-04, -3.93916e-02),
    18: ("amine, primary", 3.71375e+01, -2.66753e+00, 1.01483e-03, 2.14233e-01),
    19: ("amine, secondary", -5.03710e+02, 1.04092e+00, -4.12746e-03, 1.82790e-01),
    20: ("amine, tertiary", -3.59763e+01, -4.08458e-01, 1.67264e-03, -9.98919e-02),
    21: ("amine, aromatic", -6.09432e+02, 1.50436e+00, -9.09024e-04, -1.35495e-01),
    22: ("amide, primary", -1.02367e+02, -7.16253e-01, -2.90670e-04, -5.88556e-01),
    23: ("amide, secondary", -1.93802e+03, 6.48262e-01, 1.73245e-03, 3.47940e-02),
    24: ("amide, tertiary", -5.26919e+00, 3.06435e-01, 3.25397e-03, -6.81506e-01),
    25: ("carbonylperoxynitrate (PAN-type)", -2.84042e+02, -6.25424e-01, -8.22474e-04, -8.80549e-02),
    26: ("peroxide (C-O-O-C)", 1.50093e+02, 2.39875e-02, -3.37969e-03, 1.52789e-02),
    27: ("hydroperoxide (C-O-O-H)", -2.03387e+01, -5.48718e+00, 8.39075e-03, 1.07884e-01),
    28: ("carbonylperoxyacid (C(=O)-O-O-H)", -8.38064e+02, -1.09600e+00, -4.24385e-04, 2.81812e-01),
    29: ("nitrophenol", -5.27934e+01, -4.63689e-01, -5.11647e-03, 3.84965e-01),
    30: ("nitroester", -1.61520e+03, 9.01669e-01, 1.44536e-03, 2.66889e-01),
}

#: SMARTS catalogue keyed by the paper's group index k. Patterns are
#: authored from the group definitions and footnotes of paper Table 5
#: (clean-room, no third-party pattern files consulted). ``None`` marks
#: the algorithmic groups whose counts cannot be expressed as a single
#: SMARTS pattern (see :data:`SIMPOL1_ALGORITHMIC_GROUPS`).
SIMPOL1_SMARTS = {
    1: "[#6]",
    2: None,
    3: None,
    4: None,
    5: "[C]=[C]",
    6: "[CR]=[CR]-[CR](=[OX1])",
    7: "[OX2;H1][C]",
    8: "[CX3;H1](=[OX1])[#6]",
    9: "[#6][CX3](=[OX1])[#6]",
    10: "[CX3](=[OX1])[OX2;H1]",
    11: "[CX3](=[OX1])[OX2][#6]",
    12: "[#6;!a]-[OX2;!R]-[#6;!a]",
    13: "[OX2;r]([#6])[#6]",
    14: "[a]-[OX2;!R]-[#6]",
    15: "[$([#6][OX2][NX3](=[OX1])(=[OX1])),$([#6][OX2][N+](=[OX1])([OX1-]))]",
    16: "[$([#6][NX3](=[OX1])=[OX1]),$([#6][N+](=[OX1])[O-])]",
    17: "[OX2;H1][c]",
    18: "[CX4][NX3;H2]",
    19: "[CX4][NX3;H1][CX4]",
    20: "[CX4][NX3;H0]([CX4])[CX4]",
    21: "[a][NX3]",
    22: "[CX3](=[OX1])[NX3;H2]",
    23: "[CX3](=[OX1])[NX3;H1]",
    24: "[CX3](=[OX1])[NX3;H0]",
    25: "[$([CX3](=[OX1])[OX2][OX2][NX3](=[OX1])(=[OX1])),$([CX3](=[OX1])[OX2][OX2][N+](=[OX1])([OX1-]))]",
    26: "[#6][OX2][OX2][#6]",
    27: "[#6][OX2][OX2;H1]",
    28: "[CX3](=[OX1])[OX2][OX2;H1]",
    29: "[$([OX2;H1][c][c][NX3](=[OX1])=[OX1]),$([OX2;H1][c][c][c][NX3](=[OX1])=[OX1]),$([OX2;H1][c][c][c][c][NX3](=[OX1])=[OX1]),$([OX2;H1][c][c][N+](=[OX1])[O-]),$([OX2;H1][c][c][c][N+](=[OX1])[O-]),$([OX2;H1][c][c][c][c][N+](=[OX1])[O-])]",
    30: "[$([NX3](=[OX1])(=[OX1])[CX4][CX3](=[OX1])[OX2]),$([N+](=[OX1])([O-])[CX4][CX3](=[OX1])[OX2])]",
}

#: Groups counted algorithmically in the wrapper layer (no SMARTS):
#:   * 2 - "acid-side carbons" of an amide is a global graph-cut property
#:     (unbounded acyl-chain length) that no fixed-length SMARTS pattern
#:     can express; computed with
#:     :func:`aerosolprops.topology.fragments_excluding_bond`.
#:   * 3/4 - aromatic / non-aromatic ring counts via SSSR ring
#:     perception (the paper's footnote examples equal the cyclomatic
#:     number: naphthalene 2, anthracene 3, biphenyl 2).
SIMPOL1_ALGORITHMIC_GROUPS = {
    2: "amide acid-side carbon count (graph cut on the C-N bond)",
    3: "number of aromatic SSSR rings",
    4: "number of non-aromatic SSSR rings",
}


def simpol1_bk(k, T=298.15):
    """Group contribution term b_k(T) of SIMPOL.1 (paper Eq. 8).

    Parameters
    ----------
    k : int
        Group index, 0-30 as in paper Table 5.
    T : float, optional
        Temperature [K].

    Returns
    -------
    bk : float or None
        b_k(T) = B_k1/T + B_k2 + B_k3*T + B_k4*ln(T), or ``None`` if `k`
        is not a valid group index or `T` is not positive.

    Examples
    --------
    Table 6 of the paper lists b_k(293.15 K) to three significant
    figures; e.g. b1 = -0.438 ("each carbon lowers the vapour pressure
    by about a third of a decade") and b10 = -3.58:

    >>> round(simpol1_bk(1, 293.15), 3)
    -0.438
    >>> round(simpol1_bk(10, 293.15), 2)
    -3.58
    >>> simpol1_bk(31, 293.15) is None
    True
    """
    try:
        _, b1, b2, b3, b4 = SIMPOL1_GROUPS[k]
        return b1 / T + b2 + b3 * T + b4 * log(T)
    except (KeyError, ValueError, ZeroDivisionError, TypeError):
        return None


def simpol1_log10_vapour_pressure(counts, T=298.15):
    """Estimate log10 of the sub-cooled liquid vapour pressure [log10 Pa].

    Implements paper Eq. 6 and 8: the zeroth-group constant is always
    added exactly once (a count under key 0, if present, is ignored),
    and the atm-to-Pa conversion constant is included:

    .. math::

        \\log_{10}(P/\\mathrm{Pa}) = b_0(T) + \\sum_{k=1}^{30}
        \\nu_k\\, b_k(T) + \\log_{10}(101325)

    Parameters
    ----------
    counts : dict
        SIMPOL.1 group counts as produced by
        :func:`aerosolprops.vapour_pressure.simpol1_fragmentation`,
        mapping group index (int, 0-30) to count.
    T : float, optional
        Temperature [K]. The method was fitted on 273.15-393.15 K; see
        :data:`SIMPOL1_TEMPERATURE_RANGE`.

    Returns
    -------
    log10P : float or None
        Base-10 logarithm of the vapour pressure [Pa], or ``None`` on
        invalid input (unknown group index, non-positive temperature).

    See Also
    --------
    simpol1_bk : the per-group temperature dependence.

    Examples
    --------
    The paper's section 3.3 worked examples at 293.15 K. Cyclohexene
    (6 carbons, one alicyclic ring, one C=C; paper calculated -0.94 in
    log10 atm, measured -1.08):

    >>> lp = simpol1_log10_vapour_pressure({1: 6, 4: 1, 5: 1}, T=293.15)
    >>> round(lp, 4)
    4.0567
    >>> round(lp - LOG10_ATM_TO_PA, 2)
    -0.95

    Adipic acid (paper calculated -7.99, measured -8.49):

    >>> lp = simpol1_log10_vapour_pressure({1: 6, 10: 2}, T=293.15)
    >>> round(lp, 4)
    -2.9843
    >>> round(lp - LOG10_ATM_TO_PA, 2)
    -7.99

    Invalid input returns ``None`` rather than raising:

    >>> simpol1_log10_vapour_pressure({99: 1}) is None
    True
    """
    try:
        total = simpol1_bk(0, T)
        if total is None:
            return None
        for k, nu in counts.items():
            if k == 0 or not nu:
                continue
            bk = simpol1_bk(k, T)
            if bk is None:
                return None
            total += nu * bk
        return total + LOG10_ATM_TO_PA
    except (AttributeError, TypeError):
        return None


def simpol1_vapourisation_enthalpy(counts, T=298.15):
    """Estimate the enthalpy of vaporisation [J/mol] (paper Eq. 10).

    By Clausius-Clapeyron (paper Eq. 9), d(log10 p)/d(1/T) =
    -dHvap/(2.303 R). Differentiating b_k(T) with respect to 1/T gives
    B_k1 - B_k3*T**2 - B_k4*T, so (paper Eq. 10)

    .. math::

        \\Delta H_{\\mathrm{vap}}(T) = -2.303\\,R \\sum_k \\nu_k
        \\left(B_{k,1} - B_{k,3} T^2 - B_{k,4} T\\right)
        = 2.303\\,R \\sum_k \\nu_k
        \\left(B_{k,3} T^2 + B_{k,4} T - B_{k,1}\\right)

    with R = 8.314 J/mol/K as used by the paper. The atm-to-Pa constant
    is T-independent and drops out of the derivative. The overall sign
    is fixed by the paper twice over: physically "for any real compound
    in the liquid state dHvap > 0, and d(dHvap)/dT < 0", and formally by
    its Eq. 11, d(dHvap)/dT = 2.303*R*sum(nu_k*(2*B_k3*T + B_k4)), which
    only the minus-sign form above reproduces. (The minus sign was
    dropped in one intermediate transcription of the method notes; the
    paper text governs.)

    Parameters
    ----------
    counts : dict
        SIMPOL.1 group counts (see
        :func:`simpol1_log10_vapour_pressure`).
    T : float, optional
        Temperature [K].

    Returns
    -------
    Hvap : float or None
        Enthalpy of vaporisation [J/mol], or ``None`` on invalid input.
        The paper reports a fit standard deviation of 8.9 kJ/mol at
        333.15 K for this property.

    Examples
    --------
    Consistency with the temperature slope of
    :func:`simpol1_log10_vapour_pressure` (Clausius-Clapeyron identity,
    checked by finite differences for cyclohexene):

    >>> counts = {1: 6, 4: 1, 5: 1}  # cyclohexene
    >>> T = 300.0
    >>> h = simpol1_vapourisation_enthalpy(counts, T)
    >>> dT = 1e-4
    >>> lp_hi = simpol1_log10_vapour_pressure(counts, T + dT)
    >>> lp_lo = simpol1_log10_vapour_pressure(counts, T - dT)
    >>> slope = (lp_hi - lp_lo) / (1.0 / (T + dT) - 1.0 / (T - dT))
    >>> round(h / (-2.303 * SIMPOL1_GAS_CONSTANT * slope), 6)
    1.0

    Positive and decreasing in T, as the paper requires:

    >>> h1 = simpol1_vapourisation_enthalpy(counts, 333.15)
    >>> h2 = simpol1_vapourisation_enthalpy(counts, 353.15)
    >>> round(h1 / 1000.0, 1)
    33.4
    >>> h2 < h1 > 0
    True
    """
    try:
        total = 0.0
        for k, nu in counts.items():
            if k == 0 or not nu:
                continue
            _, b1, _, b3, b4 = SIMPOL1_GROUPS[k]
            total += nu * (b3 * T * T + b4 * T - b1)
        # The zeroth group contributes like any other group (nu_0 = 1);
        # a count under key 0 is ignored, as in
        # simpol1_log10_vapour_pressure.
        _, b1, _, b3, b4 = SIMPOL1_GROUPS[0]
        total += b3 * T * T + b4 * T - b1
        return 2.303 * SIMPOL1_GAS_CONSTANT * total
    except (KeyError, TypeError):
        return None
