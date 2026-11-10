# SPDX-License-Identifier: MIT
"""Nannoolal Part 3 group-contribution coefficients for the vapour pressure.

Catalogue of pure coefficients for the method of Nannoolal, Rarey &
Ramjugernath, "Estimation of pure component properties. Part 3. Estimation
of the vapor pressure of non-electrolyte organic compounds via group
contributions and group interactions", Fluid Phase Equilibria 269 (2008)
117-133, doi:10.1016/j.fluid.2008.04.020. This module contains data and
pure estimators only; no RDKit or thermo import.

Model (paper Eq. (6)-(8)); n is the number of non-hydrogen atoms, T_rb the
reduced temperature T/T_b:

.. math::

    \\log_{10}\\frac{P^s}{1\\ atm} = (4.1012 + dB)\\,
        \\frac{T_{rb} - 1}{T_{rb} - 1/8}

.. math::

    dB = \\sum_i N_i C_i + GI - 0.176055

with the group-interaction term (paper Eq. (8), ordered double sum over the
m non-additive groups of the molecule, :math:`C_{i-i} = 0`,
:math:`C_{i-j} = C_{j-i}`):

.. math::

    GI = \\frac{1}{n (m-1)} \\sum_{i=1}^{m} \\sum_{j=1}^{m} C_{i-j} .

All dB contributions are printed x 1e3 in paper Tables 4/5/6 and are
stored here divided by 1000. The estimator returns log10(P/Pa) (strict
SI, the paper's atm basis plus log10(101325)), following the EVAPORATION
module convention; Appendix A also inverts Eq. (6) to a boiling point.

Transcription provenance and rulings (all open questions closed against
the paper text; the method notes lose any conflict):

- The two numbers printed per Table 1 row are (ID, PR) per its footnote
  a - the priority number, not a second ID.
- ID 41 (NH2-(Ca)): the PR cell is blank in print; PR 95 is inferred from
  the series renumbering (Part 1 PR + 1 throughout the amine family; no
  other group of either paper holds 95).
- Table 4 carries 115 value rows. Defined but valueless: 106 (imidazole)
  and 108 ((S)C=N) among the main groups and 129 (para pair) among the
  corrections - molecules counting them are refused, never zeroed.
  ID 109 ("CO connected to N", example "methyl thioacetate" which
  contains no nitrogen - the same unusable definition as in Part 1) also
  has no dB value and is excluded from the catalogue entirely; real
  amides are covered by 48-50, 99, 100. IDs 112 and 114 are numbering
  gaps. IDs 98 and 110 have dB values but no group definition in either
  Part 1 or Part 3; their values are parked in
  :data:`NANNOOLAL_VP_UNDEFINED_VALUE_GROUPS` (unmatchable, never
  counted).
- Table 4 footnote a flags ID 91 as a questionable contribution; the
  value is recorded verbatim.
- The second-order correction for a carbonyl connected to an sp2 carbon
  is ID 134 here (Part 1 numbered it 118; Part 3 reassigns 118 to the
  1,2-diketone main group, which the paper declares "no estimation
  possible" - counted by the fragmenter, refused by the estimator).
- ID 214 is new in Part 3: the aromatic bridge carbon of biphenyl-type
  systems (three aromatic neighbours, but only two aromatic bonds).
- ID 217 (Si-(F,Cl,Br,I) correction) and the Si groups 70/71/93/215/216
  have values but lie outside the C/N/O/S/F/Cl/Br/I element scope the
  aeroprops fragmenter parameterises; the values are transcribed for
  completeness.
- Paper Eq. (9), a low-pressure correction for mono-functional alcohols,
  is NOT implemented: as printed it is a single temperature-dependent
  sigmoid term that is physically unrealisable (it does not vanish in the
  limits and no auxiliary quantity in the paper defines its shape), and
  UManSysProp does not implement it either. Locked ruling from the method
  notes (OQ-2); recorded here so the omission is deliberate.
- Interactions are transcribed from Table 6 (row IDs 135-218): 44 valued
  pairs plus row 218 (COOH-NH2) which is printed "Do not estimate" and
  must refuse, never zero. A pair absent from Table 6 contributes zero
  (the Part 1 convention for unregressed pairs). Table 6 footnotes flag
  row 187 (Epox-epox) as questionable; the value is recorded verbatim.
- The three SOA extension groups (IDs 301-303, aeroprops slot numbering)
  carry dB contributions from Compernolle, Ceulemans & Muller, Atmos.
  Chem. Phys. 10 (2010) 6271-6282, Table 3 (Part 3 dB column); no group
  interactions can be derived for them, so molecules containing them take
  the no-GI branch with reduced confidence.
- Worked-example adjudications (paper Tables 15a-15h): 15d prints
  "dB = 0.2939942 - 1.0923847" (0.2939942 is 15c's sum carried over by
  mistake); the arithmetic of 15d's own rows gives dB = 0.9163297, which
  reproduces the printed pressure (6.52 kPa). 15g prints dB = 1.7493490,
  contradicting its own formula, which evaluates to 1.6131919 (P = 0.404
  kPa against the experimental 0.410 kPa; the printed dB gives 0.35421);
  the formula wins (note OQ-3). 15f's printed inversion formula, read
  literally, yields ~210 K instead of the printed 520.85 K; the Appendix
  A form Tb = T (L - K) / (L/8 - K) with L = log10(P/atm), K = 4.1012 +
  dB reproduces 520.85 K and is what :func:`nannoolal_vp_boiling_point`
  implements.

References
----------
.. [1] Nannoolal, Y., Rarey, J., Ramjugernath, D. "Estimation of pure
   component properties. Part 3." Fluid Phase Equilibria 269, 117-133
   (2008). doi:10.1016/j.fluid.2008.04.020
.. [2] Compernolle, S., Ceulemans, K., Muller, J.-F. "EVAPORATION: a new
   vapour pressure estimation method for organic molecules including
   non-additivity and intramolecular hydrogen bonding." Atmospheric
   Chemistry and Physics 10, 6271-6282 (2010). doi:10.5194/acp-10-6271-2010
"""

from math import log10

__all__ = [
    "NANNOOLAL_VP_BASE",
    "NANNOOLAL_VP_GROUPS",
    "NANNOOLAL_VP_UNDEFINED_VALUE_GROUPS",
    "NANNOOLAL_VP_NO_PARAMETER",
    "NANNOOLAL_VP_INTERACTION_CLASSES",
    "NANNOOLAL_VP_INTERACTIONS",
    "NANNOOLAL_VP_INTERACTION_IDS",
    "NANNOOLAL_VP_EXTENSIONS",
    "nannoolal_vp_group_sum",
    "nannoolal_vp_interaction_counts",
    "nannoolal_vp_gi",
    "nannoolal_vp_dB",
    "nannoolal_vp_log10_vapour_pressure",
    "nannoolal_vp_boiling_point",
]

#: log10(101325): the paper fits log10(P/1 atm); estimators return
#: log10(P/Pa) (strict SI, EVAPORATION module convention).
LOG10_ATM_TO_PA = log10(101325.0)

#: Constants of the Part 3 equations: paper Eq. (6) slope constant 4.1012
#: (mean from several hundred nonpolar components), the 1/8 reduced-
#: temperature shift of the fitted C ~ Tb/8, and the Eq. (7) intercept.
NANNOOLAL_VP_BASE = {
    "log10_term": 4.1012,
    "trb_shift": 0.125,
    "constant": -0.176055,
}

#: Group catalogue, keyed by the paper's group ID.
#: Each record is ``(name, dB contribution, priority number)``. The dB
#: contribution is the paper Table 4 / Table 5 value divided by 1000
#: (both tables print x 1e3); records with no value in the paper carry
#: ``None`` and belong to :data:`NANNOOLAL_VP_NO_PARAMETER`. The priority
#: number (PR, paper Table 1: lower wins when two groups are applicable)
#: is ``None`` for second-order corrections (IDs 119-134, 217) and
#: extension groups (IDs 301-303).
NANNOOLAL_VP_GROUPS = {
    # --- Main groups, hydrocarbon skeleton (paper Table 1 / Table 4) ---
    1: ("CH3-(ne)", 0.0133063, 105),
    2: ("CH3-(e)", 0.0918000, 103),
    3: ("CH3-(a)", 0.0501939, 104),
    4: ("C(c)H2", 0.0546564, 112),
    5: ("C(c)H", 0.0457437, 119),
    6: ("C(c)", -0.0317531, 121),
    7: ("C(c)(e)", 0.0378485, 108),
    8: ("C(c)(a)", 0.0961386, 109),
    9: ("C(r)H2", 0.0222573, 113),
    10: ("C(r)H", 0.0328162, 118),
    11: ("C(r)", 0.0048500, 120),
    12: ("C(r)(e,c)", 0.0236411, 110),
    13: ("C(r)(e,r)", 0.0498237, 111),
    14: ("C(r)(Ca)", -0.0036950, 107),
    15: ("C(a)H aromatic CH", 0.0327177, 106),
    16: ("C(a)(ne)", 0.0698796, 117),
    17: ("C(a)(e)", 0.0415534, 114),
    18: ("(a)=C(a)2(a) fused aromatic C", 0.0437191, 115),
    214: ("C(a)C(a)C2(a) bridge aromatic C", 0.0869450, 116),
    # --- Main groups, halogens (paper Table 1 / Table 4) ---
    19: ("F-(C,Si)", 0.0795429, 87),
    20: ("FC=C", 0.0512880, 85),
    21: ("F-((C,Si)-([F]))-b", 0.0420887, 81),
    22: ("F-((C,Si)-([F,Cl]))-a", 0.0569998, 84),
    23: ("F-((C,Si)-([F,Cl]2))", 0.1421060, 83),
    24: ("F-(C(a))", 0.0459652, 86),
    102: ("F-((C,Si)-(Cl))-b", 0.0640566, 82),
    25: ("Cl-(C,Si)", 0.0936679, 71),
    26: ("Cl-((C,Si)-([F,Cl]))", 0.0678082, 70),
    27: ("Cl-((C,Si)-([F,Cl]2))", 0.0559304, 68),
    28: ("Cl-(C(a))", 0.0460435, 72),
    29: ("CCl=C", 0.0849162, 69),
    30: ("Br-(C/Si(na))", 0.1049291, 65),
    31: ("Br-(C(a))", -0.0401837, 66),
    32: ("I-(C,Si)", 0.1343501, 63),
    # --- Main groups, oxygen (paper Table 1 / Table 4) ---
    33: ("OH tert", 0.7193666, 91),
    34: ("OH ((C,Si)2H(C,Si)) secondary", 0.7584218, 90),
    35: ("OH (>C4) primary", 0.7007226, 88),
    36: ("OH (<C5)", 0.7560824, 92),
    37: ("OH (Ca) phenol", 0.4418437, 89),
    38: ("(C,Si)-O-(C,Si) ether", 0.1084964, 94),
    39: ("(OC2) epoxide", 0.2869731, 50),
    44: ("COOH-(C)", 1.0741000, 23),
    45: ("(C)-COO-(C) ester", 0.3557381, 24),
    46: ("HCOO-(C) formate ester", 0.3505184, 26),
    47: ("C(c)OO- (lactone)", 0.2928046, 25),
    48: ("CON (disubstituted amide)", 0.2692471, 49),
    49: ("CONH (monosubstituted amide)", 0.7369540, 48),
    50: ("CONH2 (amide)", 1.2160700, 27),
    51: ("O=C-(Cna)2 ketone", 0.2558480, 4),
    52: ("CHO-(Cna) aldehyde", 0.2529059, 53),
    90: ("CHO-(Ca) aldehyde", 0.1831130, 52),
    65: ("(C(a))-O(a)-(C(a)) aromatic ether", 0.0488839, 93),
    72: ("NO3 (nitrate ester)", 0.6035347, 13),
    74: ("O=N-O-(C) nitrite", 1.3174360, 22),
    73: ("PO(O)3 (phosphate)", 0.6620582, 9),
    75: ("ONC (oxime)", 0.5109666, 29),
    76: ("C=O-O-C=O (anhydride)", 0.6813525, 11),
    79: ("O=C(O)2 (carbonate)", 0.3182350, 14),
    80: ("OCN (isocyanate)", 0.4358446, 28),
    94: ("O-O (peroxide)", 0.1867950, 31),
    96: ("(C=O-O-C=O)r (cyclic anhydride)", 0.5951778, 10),
    99: ("OCON (carbamate)", 0.6129546, 2),
    100: ("N-(C=O)-N (urea)", 0.2589924, 3),
    103: ("OCOO (C,O carbonate)", 0.6602247, 33),
    # --- Main groups, sulfur / selenium (paper Table 1 / Table 4) ---
    53: ("SH-(C) thiol", 0.1232143, 73),
    54: ("(C)-S-(C) sulfide", 0.1273380, 74),
    55: ("(C)-S-S-(C) disulfide", 0.2222789, 51),
    56: ("S(a) aromatic sulfur", 0.0201604, 75),
    81: ("SCN-(C) thiocyanate", 0.2186012, 19),
    82: ("(C)-SO2-(C) sulfone", 0.3815945, 17),
    104: ("SO4 (sulfate)", 0.5547941, 34),
    105: ("SO2-N (sulfonamide)", 0.4207591, 35),
    107: ("S=O (sulfoxide)", -0.2372124, 37),
    116: ("Se", 0.1915058, 46),
    # --- Main groups, nitrogen (paper Table 1 / Table 4) ---
    40: ("NH2-(C,Si)", 0.2519212, 96),
    41: ("NH2-(Ca)", 0.3617760, 95),  # PR blank in print; 95 inferred
    42: ("(C,Si)-NH-(C,Si) secondary amine", 0.1937667, 100),
    43: ("(C,Si)2-N-(C,Si) tertiary amine", -0.1027252, 101),
    57: ("(C)-C=N cyanide", 0.2261873, 55),
    66: ("N(a)(r5) aromatic N, five-ring", 0.3750486, 98),
    67: ("N(a)(r6) aromatic N, six-ring", 0.1263340, 97),
    68: ("NO2-(C)", 0.3758217, 20),
    69: ("NO2-(C(a))", 0.2382066, 21),
    77: ("COCl (acid chloride)", 0.5641116, 18),
    91: ("(C,Si)=N", 1.1781950, 102),  # Table 4 footnote a: questionable
    92: ("(O=C-(C)2)a ketone, aromatic side", 0.1583258, 54),
    97: ("(C,Si)r-NH-(Cr,Si)r cyclic secondary amine", 0.2687081, 99),
    101: ("(C,Si)2-N-(C,Si)2 quaternary amine", -0.3164392, 32),
    111: ("(N)-C=N (cyanamide)", 0.3194879, 41),
    115: ("ON=(C,Si) (isoxazole)", 0.3051341, 45),
    # --- Main groups defined but valueless (Table 1, no Table 4 row) ---
    106: ("imidazole (CNC=NC=)", None, 36),
    108: ("(S)-C=N", None, 38),
    118: ("O=C-C=O (1,2-diketone)", None, 1),  # "no estimation possible"
    # --- Main groups, unsaturated carbon (paper Table 1 / Table 4) ---
    58: ("C(c)=C(c)", 0.0864601, 62),
    59: ("C(c)=C(c)(C(a))", 0.2242238, 59),
    60: ("-(e)C(c)=C(c)", 0.1349382, 58),
    61: ("H2C=C (1-ene)", 0.0342541, 57),
    62: ("C(r)=C(r)", 0.0974210, 60),
    63: ("C=C", 0.2066665, 61),
    64: ("HC=C (1-yne)", 0.1280247, 56),
    87: ("C=C=C (cumulated)", 0.1683007, 5),
    88: ("C=C-C=C (conjugated, ring)", 0.1085277, 6),
    89: ("C=C-C=C (conjugated, chain)", 0.2137165, 7),
    95: ("C=C-C=C (conjugated triple)", 0.3922006, 8),
    # --- Si / Ge / Sn / P / As / Se / B / Al (values only, out of the
    #     element scope of the aeroprops fragmenter) ---
    70: ("Si", 0.0028992, 80),
    216: ("Si(C,H)0", -0.0811543, 76),
    215: ("Si(C,H)1", -0.0665670, 77),
    93: ("Si(C,H)2", -0.0473420, 78),
    71: ("Si(C,H)3", 0.0093624, 79),
    86: ("(C)2-Ge-(C)2", 0.0485026, 67),
    85: ("GeCl3", 0.1869204, 12),
    83: ("(C)2-Sn-(C)2", 0.0802735, 64),
    84: ("AsCl2", 0.2313919, 16),
    113: ("P (phosphine)", 0.1188412, 43),
    117: ("Al", 0.4235251, 47),
    78: ("B(O-)3 (borate ester)", 0.3913697, 15),
    # --- Second-order corrections (paper Table 2 / Table 5, no PR) ---
    119: ("(C=O)-C([F,Cl]2,3)", 0.0343545, None),
    120: ("(C=O)-(C([F,Cl]2,3))2", 0.0025030, None),
    121: ("C-[F,Cl]3", -0.0833326, None),
    122: ("(C)2C-[F,Cl]2", -0.0644854, None),
    123: ("no hydrogen", -0.1259208, None),
    124: ("one hydrogen", -0.0472962, None),
    125: ("3/4 ring", 0.0339765, None),
    126: ("five-ring", -0.0070982, None),
    127: ("ortho pair (exclusive)", -0.0450531, None),
    128: ("meta pair (exclusive)", -0.0032036, None),
    129: ("para pair (exclusive)", None, None),  # no Table 5 row
    130: ("(C=)(C)C-CC3", -0.0206706, None),
    131: ("C2C-CC2", -0.0363170, None),
    132: ("C3C-CC2", -0.0011994, None),
    133: ("C3C-CC3", 0.1237433, None),
    134: ("C=C-C=O (carbonyl on sp2 carbon)", -0.0159694, None),
    217: ("Si-(F,Cl,Br,I)", 0.0367574, None),
    # --- SOA extension groups (Compernolle 2010 Table 3, no PR, no GI) ---
    301: ("OOH (hydroperoxide)", 0.88842, None),
    302: ("C(=O)OOH (peroxy acid)", 0.92602, None),
    303: ("C(=O)OONO2 (PAN)", 0.51903, None),
}

#: Table 4 rows whose group is defined in neither Part 1 nor Part 3
#: (IDs 98, 110): values transcribed for completeness, but no SMARTS can
#: be authored from a missing definition and the groups are never counted.
NANNOOLAL_VP_UNDEFINED_VALUE_GROUPS = {
    98: 0.1833467,
    110: 0.0370590,
}

#: Groups the paper defines but gives no dB contribution: counting one of
#: them is a refusal to estimate, never a zero contribution. IDs 106/108/
#: 118 have no Table 4 row (118 additionally prints "do not estimate");
#: ID 129 (para pair) is absent from Table 5.
NANNOOLAL_VP_NO_PARAMETER = frozenset((106, 108, 118, 129))

#: Non-additive interaction classes of paper Table 3, keyed by the class
#: letter of the Nannoolal series (Part 3 prints long names: OH, OH(a),
#: COOH, EtherO, Epox, Ester, Ketone, Alde, AO, Teth, Ats, SH, NH2, NH,
#: OCN, CN, Nitro, AN5, AN6 - the letters A-S of Part 1 Table 9 with the
#: identical memberships, Part 3's own OH/OH(a) split included).
NANNOOLAL_VP_INTERACTION_CLASSES = {
    "A": ("OH alkanol (34, 35, 36)", (34, 35, 36)),
    "B": ("OH(a) phenol (37)", (37,)),
    "C": ("COOH carboxylic acid (44)", (44,)),
    "D": ("EtherO ether (38)", (38,)),
    "E": ("Epox (39)", (39,)),
    "F": ("Ester (45, 46, 47)", (45, 46, 47)),
    "G": ("Ketone (51, 92)", (51, 92)),
    "H": ("Alde aldehyde (52, 90)", (52, 90)),
    "I": ("AO aromatic oxygen (65)", (65,)),
    "J": ("Teth sulfide (54)", (54,)),
    "K": ("Ats aromatic sulfur (56)", (56,)),
    "L": ("SH thiol (53)", (53,)),
    "M": ("NH2 primary amine (40, 41)", (40, 41)),
    "N": ("NH secondary amine (42, 97)", (42, 97)),
    "O": ("OCN isocyanate (80)", (80,)),
    "P": ("CN cyanide (57)", (57,)),
    "Q": ("Nitro (69)", (69,)),
    "R": ("AN5 aromatic N in 5-ring (66)", (66,)),
    "S": ("AN6 aromatic N in 6-ring (67)", (67,)),
}

#: Group-interaction dB contributions of paper Table 6, keyed by the
#: sorted ``(class_i, class_j)`` letter pair (interactions are symmetric,
#: ``C_i-j = C_j-i``). Paper values are printed x 1e3 and stored divided
#: by 1000. Only the 44 pairs with a regressed parameter are listed; the
#: pair of row 218 (C-M, COOH-NH2) is printed "Do not estimate" and is
#: represented by its absence here while its row ID remains in
#: :data:`NANNOOLAL_VP_INTERACTION_IDS` - it must be treated as missing,
#: not as zero.
NANNOOLAL_VP_INTERACTIONS = {
    ("A", "A"): -0.5615153,
    ("A", "M"): 1.0676660,
    ("A", "N"): 0.0424825,
    ("A", "D"): -0.7995332,
    ("A", "E"): -0.6182760,
    ("A", "F"): -1.7976930,
    ("A", "G"): -1.1815990,
    ("A", "P"): 1.4312430,
    ("A", "I"): 1.1320400,
    ("B", "B"): -0.0976205,
    ("B", "D"): -0.7516676,
    ("B", "F"): 0.5484352,
    ("M", "M"): 1.0858320,
    ("M", "N"): -0.2067811,
    ("D", "M"): -0.1982791,
    ("F", "M"): -1.6764770,
    ("M", "Q"): 1.6590340,
    ("N", "N"): -0.3071018,
    ("D", "N"): 0.0654421,
    ("L", "L"): 0.2403037,
    ("C", "C"): -2.6017090,
    ("C", "G"): -0.7878563,
    ("O", "O"): -3.9291300,
    ("D", "D"): 0.1446074,
    ("D", "E"): 1.1189580,
    ("D", "F"): -0.2257802,
    ("D", "G"): 1.9812980,
    ("D", "H"): 0.3627540,
    ("D", "J"): -1.4250170,
    ("D", "P"): 0.7433353,
    ("E", "E"): -3.7488180,  # Table 6 footnote a: questionable
    ("F", "F"): 0.9203138,
    ("F", "G"): 1.5941640,
    ("F", "P"): 0.1081305,
    ("F", "I"): 1.5903210,
    ("G", "G"): -1.2700830,
    ("H", "H"): 0.9467309,
    ("H", "I"): 0.7053049,
    ("J", "J"): 0.8383372,
    ("Q", "Q"): -1.5013550,
    ("K", "R"): 0.6750414,
    ("P", "S"): 0.9944996,
    ("I", "R"): 0.1355896,
    ("S", "S"): -0.0296785,
}

#: Paper Table 6 row IDs (135-218) for every tabulated interaction pair,
#: including the row printed without a value (218, COOH-NH2).
NANNOOLAL_VP_INTERACTION_IDS = {
    ("A", "A"): 135,
    ("A", "M"): 136,
    ("A", "N"): 137,
    ("A", "D"): 140,
    ("A", "E"): 141,
    ("A", "F"): 142,
    ("A", "G"): 143,
    ("A", "P"): 145,
    ("A", "I"): 146,
    ("B", "B"): 148,
    ("B", "D"): 151,
    ("B", "F"): 152,
    ("M", "M"): 157,
    ("M", "N"): 158,
    ("D", "M"): 159,
    ("F", "M"): 160,
    ("M", "Q"): 162,
    ("N", "N"): 165,
    ("D", "N"): 166,
    ("L", "L"): 170,
    ("C", "C"): 172,
    ("C", "G"): 175,
    ("O", "O"): 176,
    ("D", "D"): 178,
    ("D", "E"): 179,
    ("D", "F"): 180,
    ("D", "G"): 181,
    ("D", "H"): 182,
    ("D", "J"): 183,
    ("D", "P"): 185,
    ("E", "E"): 187,
    ("F", "F"): 189,
    ("F", "G"): 190,
    ("F", "P"): 192,
    ("F", "I"): 193,
    ("G", "G"): 194,
    ("H", "H"): 201,
    ("H", "I"): 204,
    ("J", "J"): 205,
    ("Q", "Q"): 206,
    ("K", "R"): 208,
    ("P", "S"): 209,
    ("I", "R"): 210,
    ("S", "S"): 212,
    ("C", "M"): 218,  # printed "Do not estimate" (no parameter)
}

#: SOA extension groups (aeroprops slot IDs 301-303). The dB
#: contributions themselves are merged into :data:`NANNOOLAL_VP_GROUPS`;
#: membership in this table marks a group as having no interaction
#: parameters (the no-GI branch of the estimator).
NANNOOLAL_VP_EXTENSIONS = {
    301: (
        "OOH (hydroperoxide)",
        "Compernolle et al. (2010), Atmos. Chem. Phys. 10, 6271-6282, "
        "Table 3, Part 3 dB column: 888.42 (x 1e3); no group-interaction "
        "parameters derivable.",
    ),
    302: (
        "C(=O)OOH (peroxy acid)",
        "Compernolle et al. (2010), Table 3, Part 3 dB column: 926.02 "
        "(x 1e3); no group-interaction parameters derivable.",
    ),
    303: (
        "C(=O)OONO2 (PAN)",
        "Compernolle et al. (2010), Table 3, Part 3 dB column: 519.03 "
        "(x 1e3); no group-interaction parameters derivable.",
    ),
}


def nannoolal_vp_group_sum(counts):
    """Sum of the first-order and second-order dB contributions.

    Parameters
    ----------
    counts : dict
        Group counts keyed by group ID as produced by
        :func:`aeroprops.vapour_pressure.nannoolal_vp_fragmentation`.

    Returns
    -------
    total : float or None
        ``sum(N_i * dB_i)`` over all groups in `counts`, or ``None`` if a
        key is not a valid group ID or a group without a dB value
        (see :data:`NANNOOLAL_VP_NO_PARAMETER`) is present.

    Examples
    --------
    Paper Table 15a (alpha-pinene): the printed total is 0.2839146 -

    >>> round(nannoolal_vp_group_sum(
    ...     {1: 3, 9: 2, 10: 2, 11: 1, 62: 1, 125: 1, 132: 2}), 7)
    0.2839146
    >>> nannoolal_vp_group_sum({999: 1}) is None
    True
    >>> nannoolal_vp_group_sum({118: 1}) is None  # 1,2-diketone
    True
    """
    total = 0.0
    for gid, n in counts.items():
        record = NANNOOLAL_VP_GROUPS.get(gid)
        if record is None or record[1] is None:
            return None
        total += n * record[1]
    return total


def nannoolal_vp_interaction_counts(counts):
    """Count the non-additive interaction groups per Table 3 class.

    Examples
    --------
    1,2-Ethanediol (paper Table 15b): two alkanol OH (class A):

    >>> nannoolal_vp_interaction_counts({7: 2, 36: 2})
    {'A': 2}

    Diethanolamine (paper Table 15g): two OH (A) and one secondary
    amine (N), m = 3:

    >>> nannoolal_vp_interaction_counts({7: 4, 36: 2, 42: 1}) == \
{'A': 2, 'N': 1}
    True
    """
    group_to_class = {}
    for letter, (_desc, gids) in NANNOOLAL_VP_INTERACTION_CLASSES.items():
        for gid in gids:
            group_to_class[gid] = letter
    class_counts = {}
    try:
        for gid, n in counts.items():
            letter = group_to_class.get(gid)
            if letter is not None and n:
                class_counts[letter] = class_counts.get(letter, 0) + n
    except (AttributeError, TypeError):
        return None
    return class_counts


#: Status of the group-interaction term of an estimation.
#:
#: * ``"OK"`` - GI computed per Eq. (8), or zero because fewer than two
#:   interaction groups are present.
#: * ``"no-GI (extension groups present)"`` - the molecule contains a
#:   Compernolle (2010) extension group (301-303); no interaction
#:   parameters can be derived for them, so GI is skipped and the result
#:   carries reduced confidence.
#: * ``"missing interaction parameter (C, M)"``-style - the molecule
#:   needs an interaction the paper prints without a value (Table 6 row
#:   218, COOH-NH2); the estimation must be refused, never zeroed.
GI_STATUS_OK = "OK"
GI_STATUS_NO_GI = "no-GI (extension groups present)"


def nannoolal_vp_gi(counts, n):
    """Group-interaction term GI of paper Eq. (8).

    The ordered double sum over the m interaction groups of the molecule
    is evaluated through its counting equivalent: a same-class pair
    contributes ``n_g*(n_g - 1)*C_gg`` and a cross-class pair contributes
    ``2*n_g1*n_g2*C_g1g2``, normalised by ``n*(m - 1)``.
    ``C_i-i = 0``.  Pairs absent from paper Table 6 contribute zero; the
    pair printed without a value (row 218, COOH-NH2) makes the
    estimation fail.

    Parameters
    ----------
    counts : dict
        Group counts keyed by group ID.
    n : int
        Number of non-hydrogen atoms in the molecule (all elements).

    Returns
    -------
    gi : float or None
        GI, or ``None`` when a required interaction parameter does not
        exist in the paper.
    status : str
        One of the GI statuses described with :data:`GI_STATUS_OK`.

    Examples
    --------
    1,2-Ethanediol (paper Table 15b, n = 4, m = 2): the paper's printed
    correction total is -0.2807577 -

    >>> gi, status = nannoolal_vp_gi({7: 2, 36: 2}, 4)
    >>> status
    'OK'
    >>> round(gi, 7)
    -0.2807577

    Dipropyl succinate (paper Table 15f, n = 14, two esters): GI is the
    printed (1/14) x 1.8406276 = 0.1314734 -

    >>> round(nannoolal_vp_gi({1: 2, 4: 4, 7: 2, 45: 2}, 14)[0], 7)
    0.1314734

    Glycine (carboxylic acid + primary amine) needs row 218, which the
    paper prints "Do not estimate"; a molecule with an extension group
    (hydroperoxide) takes the no-GI branch; an unregressed pair such as
    two five-ring aromatic nitrogens contributes zero (no Table 6 row):

    >>> nannoolal_vp_gi({44: 1, 40: 1}, 5)
    (None, 'missing interaction parameter (C, M)')
    >>> nannoolal_vp_gi({301: 1, 34: 1}, 20)[1]
    'no-GI (extension groups present)'
    >>> nannoolal_vp_gi({15: 3, 66: 2}, 5)  # R-R not in Table 6: zero
    (0.0, 'OK')
    """
    if any(gid in NANNOOLAL_VP_EXTENSIONS for gid in counts):
        return 0.0, GI_STATUS_NO_GI
    class_counts = nannoolal_vp_interaction_counts(counts)
    if class_counts is None or n is None or n <= 0:
        return None, "invalid input"
    m = sum(class_counts.values())
    if m <= 1:
        return 0.0, GI_STATUS_OK
    total = 0.0
    letters = sorted(class_counts)
    for i, c1 in enumerate(letters):
        n1 = class_counts[c1]
        for c2 in letters[i:]:
            n2 = class_counts[c2]
            if c1 == c2 and n1 < 2:
                continue  # a single group pairs with nothing (C_i-i = 0)
            pair = (c1, c2)
            if pair not in NANNOOLAL_VP_INTERACTION_IDS:
                continue  # not in Table 6: zero contribution
            if pair not in NANNOOLAL_VP_INTERACTIONS:
                return None, "missing interaction parameter (%s, %s)" % pair
            c = NANNOOLAL_VP_INTERACTIONS[pair]
            if c1 == c2:
                total += n1 * (n1 - 1) * c
            else:
                total += 2.0 * n1 * n2 * c
    return total / (n * (m - 1)), GI_STATUS_OK


def nannoolal_vp_dB(counts, n):
    """Vapour-pressure curve slope correction dB (paper Eq. (7)).

    .. math::

        dB = \\sum_i N_i C_i + GI - 0.176055

    Parameters
    ----------
    counts : dict
        Group counts keyed by group ID as produced by
        :func:`aeroprops.vapour_pressure.nannoolal_vp_fragmentation`.
    n : int
        Number of non-hydrogen atoms in the molecule.

    Returns
    -------
    db : float or None
        dB, or ``None`` on refusal: a group without a dB value
        (:data:`NANNOOLAL_VP_NO_PARAMETER`), an unknown group ID, or a
        required interaction parameter that the paper does not provide
        (Table 6 row 218).
    status : str
        The GI status on success, otherwise the refusal reason.

    Examples
    --------
    Paper Tables 15a/15b/15h print dB = 0.1078596 / 1.1310491 /
    0.0447374; Table 15g's own formula (see the module docstring for the
    printed-value misprint) gives 1.6131919:

    >>> counts = {1: 3, 9: 2, 10: 2, 11: 1, 62: 1, 125: 1, 132: 2}
    >>> db, status = nannoolal_vp_dB(counts, 10)
    >>> status, round(db, 7)
    ('OK', 0.1078596)
    >>> db, _ = nannoolal_vp_dB({7: 2, 36: 2}, 4)
    >>> round(db, 7)
    1.1310491
    >>> db, _ = nannoolal_vp_dB(
    ...     {7: 2, 21: 2, 26: 2, 27: 1, 121: 1, 124: 1}, 5)
    >>> round(db, 7)
    0.0447374
    >>> db, _ = nannoolal_vp_dB({7: 4, 36: 2, 42: 1}, 7)
    >>> round(db, 7)
    1.6131919
    >>> nannoolal_vp_dB({44: 1, 118: 1}, 8)  # 1,2-diketone: refuse
    (None, 'no dB parameter for group 118 (O=C-C=O (1,2-diketone))')
    """
    missing = [
        gid for gid in sorted(counts)
        if counts[gid] and gid in NANNOOLAL_VP_NO_PARAMETER
    ]
    if missing:
        gid = missing[0]
        return None, "no dB parameter for group %d (%s)" % (
            gid, NANNOOLAL_VP_GROUPS[gid][0])
    total = nannoolal_vp_group_sum(counts)
    if total is None:
        return None, "invalid group counts"
    gi, gi_status = nannoolal_vp_gi(counts, n)
    if gi is None:
        return None, gi_status
    return total + gi + NANNOOLAL_VP_BASE["constant"], gi_status


def nannoolal_vp_log10_vapour_pressure(counts, n, T, Tb):
    """Base-10 logarithm of the vapour pressure [log10 Pa] (paper Eq. (6)).

    .. math::

        \\log_{10}\\frac{P^s}{1\\ atm} = (4.1012 + dB)\\,
            \\frac{T_{rb} - 1}{T_{rb} - 1/8}, \\quad T_{rb} = T / T_b

    Returns ``None`` (estimation refused or invalid input) when dB
    cannot be formed, when ``T`` or ``Tb`` is not positive, or when
    ``T <= Tb/8`` (at and below the fitted discontinuity the equation is
    outside its domain).

    Examples
    --------
    The worked examples of paper Tables 15a-15h (pressures in kPa, as
    printed):

    >>> def kPa(counts, n, T, Tb):
    ...     x = nannoolal_vp_log10_vapour_pressure(counts, n, T, Tb)
    ...     return None if x is None else round(10.0 ** x / 1000.0, 2)
    >>> kPa({1: 3, 9: 2, 10: 2, 11: 1, 62: 1, 125: 1, 132: 2},
    ...     10, 388.15, 429.0)  # Table 15a: alpha-pinene, 31.03
    31.03
    >>> kPa({7: 2, 36: 2}, 4, 410.65, 470.50)  # 15b: ethanediol, 13.05
    13.05
    >>> kPa({7: 2, 21: 6, 51: 1, 120: 1, 121: 2, 123: 1},
    ...     10, 210.16, 245.90)  # 15c: perfluoro-2-propanone, 14.63
    14.63
    >>> kPa({44: 1, 61: 1, 134: 1}, 5, 344.15, 413.60)  # 15d, 6.52
    6.52
    >>> kPa({1: 1, 7: 2, 36: 1, 45: 1}, 7, 352.65, 458.65)  # 15e, 2.24
    2.24
    >>> kPa({7: 4, 36: 2, 42: 1}, 7, 401.13, 541.15)  # 15g: 0.40 kPa
    0.4
    >>> kPa({7: 2, 21: 2, 26: 2, 27: 1, 121: 1, 124: 1},
    ...     5, 297.46, 344.25)  # 15h: R122, 17.51
    17.51
    >>> nannoolal_vp_log10_vapour_pressure({118: 1}, 6, 300.0, 400.0)
    """
    db, _status = nannoolal_vp_dB(counts, n)
    if db is None or not T > 0.0 or not Tb > 0.0:
        return None
    base = NANNOOLAL_VP_BASE
    trb = T / Tb
    if trb <= base["trb_shift"]:
        return None  # at/below the C = Tb/8 discontinuity
    return (
        (base["log10_term"] + db) * (trb - 1.0)
        / (trb - base["trb_shift"]) + LOG10_ATM_TO_PA
    )


def nannoolal_vp_boiling_point(counts, n, T, p):
    """Normal boiling point [K] from one vapour pressure point (Appendix A).

    Inverts paper Eq. (6) at ``(T, p)``:

    .. math::

        T_b = T\\,\\frac{L - K}{L/8 - K}, \\quad
        L = \\log_{10}\\frac{p}{1\\ atm}, \\quad K = 4.1012 + dB .

    The formula exactly as printed in Table 15f does not reproduce the
    table's own result (~210 K instead of 520.85 K); this verified form
    does (see the module docstring).

    Parameters
    ----------
    counts, n :
        As for :func:`nannoolal_vp_dB`.
    T : float
        Temperature of the data point [K].
    p : float
        Vapour pressure of the data point [Pa].

    Returns
    -------
    Tb : float or None
        Normal boiling point [K], or ``None`` on refusal (dB cannot be
        formed, non-positive `T`/`p`) or when no physical root exists
        (the reduced temperature solves at or below 1/8).

    Examples
    --------
    Paper Table 15f (dipropyl succinate, p = 0.4 kPa at 374.65 K; the
    paper prints Tb,calc = 520.85 K against an experimental 523.95 K):

    >>> counts = {1: 2, 4: 4, 7: 2, 45: 2}
    >>> round(nannoolal_vp_boiling_point(counts, 14, 374.65, 400.0), 2)
    520.85
    >>> nannoolal_vp_boiling_point(counts, 14, 0.0, 400.0) is None
    True
    """
    db, _status = nannoolal_vp_dB(counts, n)
    if db is None or not T > 0.0 or not p > 0.0:
        return None
    base = NANNOOLAL_VP_BASE
    L = log10(p) - LOG10_ATM_TO_PA
    K = base["log10_term"] + db
    if L == K:
        return None
    trb = (L / 8.0 - K) / (L - K)
    if not trb > base["trb_shift"]:
        return None
    return T / trb
