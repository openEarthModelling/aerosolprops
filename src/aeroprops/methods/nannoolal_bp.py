# SPDX-License-Identifier: MIT
"""Nannoolal Part 1 group-contribution coefficients for the normal boiling point.

Catalogue of pure coefficients for the method of Nannoolal, Rarey,
Ramjugernath & Cordes, "Estimation of pure component properties. Part 1.
Estimation of the normal boiling point of non-electrolyte organic compounds
via group contributions and group interactions", Fluid Phase Equilibria 226
(2004) 45-63, doi:10.1016/j.fluid.2004.09.001. This module contains data
only; the estimator lives elsewhere in the pure layer.

Model (paper Eq. (1)); n is the number of non-hydrogen atoms:

.. math::

    T_b = \\Big( \\sum_i N_i C_i + \\sum_j N_j C_j + GI \\Big)
          / \\big( n^{a} + b \\big) + c

with the group-interaction term (paper Eq. (3), ordered double sum over the
m non-additive groups of the molecule, :math:`C_{i-i} = 0`,
:math:`C_{i-j} = C_{j-i}`):

.. math::

    GI = \\frac{1}{n (m-1)} \\sum_{i=1}^{m} \\sum_{j=1}^{m} C_{i-j} .

Transcription provenance and rulings (method note 04, 2026-11-09, all open
questions closed):

- Group contributions are transcribed from paper Tables 1/3 (main groups)
  and Tables 2/4 (second-order corrections); the statistics tables 3/4 win
  any conflict with the definition tables 1/2 (note OQ-1).
- ID 12: Table 1 prints 492.0707 K (a typeset error duplicating the value of
  ID 109); Table 3 prints 250.9584 K with full regression statistics, which
  is recorded here (note OQ-1/OQ-15).
- ID 24: Table 1 omits the minus sign (19.5575); Table 3 prints -19.5575,
  which is recorded here.
- ID 78: the paper assigns the same ID, priority and contribution to two
  definitions, Si-(F,Cl) and non-cyclic boric acid ester B(O-)3; both are
  recorded in a single entry.
- ID 93: no definition in Part 1; its value 215.3416 K comes from Table 3
  and its name/priority are backfilled from Part 3 Table 1 (Si(C,H)2,
  PR 78) per the method note.
- IDs 98, 110, 112, 114 are unused gaps in the paper numbering; second-order
  corrections (IDs 118-133) carry no priority number.
- Interactions are transcribed from Table 5 (row IDs 135-212, shared across
  the Tb/dB/Tc/Pc/Vc tables of the paper series). Rows 153 (B-G) and 211
  (R-R) are printed without a value ("no parameter"), which is represented
  by their absence from :data:`NANNOOLAL_BP_INTERACTIONS` while their row
  IDs remain in :data:`NANNOOLAL_BP_INTERACTION_IDS`; this must not be
  silently treated as zero.
- Row 171 prints the pair "L-AF", the only two-letter entry of the table
  (glyph-level verified, note OQ-15). It is recorded verbatim with a
  questionable flag.
- The interaction classes of Table 9 are recorded with the Part 3 split of
  the double-listed phenol OH (group 37 belongs to class B only, note
  OQ-15).
- The three SOA extension groups (IDs 301-303, aeroprops slot numbering)
  carry boiling-point contributions from Compernolle, Ceulemans & Muller,
  Atmos. Chem. Phys. 10 (2010) 6271-6282, doi:10.5194/acp-10-6271-2010
  (open access), Table 3. No group-interaction parameters can be derived
  for them (explicitly stated there); molecules containing them must take
  the no-GI branch with reduced confidence. Their anchor boiling points
  are extrapolated values ("presumably of lower quality", their Table 2).
  Coefficients are taken from the open-access literature only.

References
----------
.. [1] Nannoolal, Y., Rarey, J., Ramjugernath, D., Cordes, W. "Estimation of
   pure component properties. Part 1." Fluid Phase Equilibria 226, 45-63
   (2004). doi:10.1016/j.fluid.2004.09.001
.. [2] Compernolle, S., Ceulemans, K., Muller, J.-F. "EVAPORATION: a new
   vapour pressure estimation method for organic molecules including
   non-additivity and intramolecular hydrogen bonding." Atmospheric
   Chemistry and Physics 10, 6271-6282 (2010). doi:10.5194/acp-10-6271-2010
"""

__all__ = [
    "NANNOOLAL_BP_BASE",
    "NANNOOLAL_BP_GROUPS",
    "NANNOOLAL_BP_INTERACTION_CLASSES",
    "NANNOOLAL_BP_INTERACTIONS",
    "NANNOOLAL_BP_INTERACTION_IDS",
    "NANNOOLAL_BP_EXTENSIONS",
    "nannoolal_bp_group_sum",
    "nannoolal_bp_interaction_counts",
    "nannoolal_bp_gi",
    "nannoolal_bp_boiling_point",
]

#: Constants of the Part 1 main equation (paper Eq. (1)),
#: ``Tb = sum / (n**a + b) + c`` with n = number of non-hydrogen atoms.
NANNOOLAL_BP_BASE = {
    "a": 0.6583,
    "b": 1.6868,
    "c": 84.3395,
}

#: Group catalogue, keyed by the paper's group ID.
#: Each record is ``(name, contribution [K], priority number)``. The priority
#: number (PR, paper Table 1: lower wins when two groups are applicable) is
#: ``None`` for second-order corrections (IDs 118-133) and extension groups
#: (IDs 301-303), for which the paper gives none.
NANNOOLAL_BP_GROUPS = {
    # --- Main groups, hydrocarbon skeleton (paper Table 1 / Table 3) ---
    1: ("CH3-(ne)", 177.3066, 104),
    2: ("CH3-(e)", 251.8338, 102),
    3: ("CH3-(a)", 157.9527, 103),
    4: ("C(c)H2", 239.4531, 111),
    5: ("C(c)H", 240.6785, 117),
    6: ("C(c)", 249.5809, 119),
    7: ("C(c)(e)", 266.8769, 107),
    8: ("C(c)(a)", 201.0115, 108),
    9: ("C(r)H2", 239.4957, 112),
    10: ("C(r)H", 222.1163, 116),
    11: ("C(r)", 209.9749, 118),
    12: ("C(r)(e,c)", 250.9584, 109),  # Table 1 misprints 492.0707 (ID 109's value)
    13: ("C(r)(e,r)", 291.2291, 110),
    14: ("C(r)(Ca)", 244.3581, 106),
    15: ("=C(a)H", 235.3462, 105),
    16: ("=C(a)(ne)", 315.4128, 115),
    17: ("=C(a)(e)", 348.2779, 113),
    18: ("(a)=C(a)2(a)", 367.9649, 114),
    # --- Main groups, halogens (paper Table 1 / Table 3) ---
    19: ("F-(C,Si)", 106.5492, 86),
    20: ("CF=C", 49.2701, 84),
    21: ("F-(C-([F,Cl]))-b", 53.1871, 80),
    22: ("F-(C-([F,Cl]))-a", 78.7578, 83),
    23: ("F-(C-([F,Cl]2))", 103.5672, 82),
    24: ("F-(C(a))", -19.5575, 85),  # Table 1 omits the minus sign
    102: ("F-(C-(Cl))-b", 111.0590, 81),
    25: ("Cl-(C,Si)", 330.9117, 72),
    26: ("Cl-((C,Si)-([F,Cl]))", 287.1863, 71),
    27: ("Cl-((C,Si)-([F,Cl]2))", 267.4170, 69),
    28: ("Cl-(C(a))", 205.7363, 73),
    29: ("CCl=C", 292.5816, 70),
    30: ("Br-(C/Si(na))", 419.4959, 66),
    31: ("Br-(C(a))", 377.6775, 67),
    32: ("I-(C,Si)", 556.3944, 64),
    # --- Main groups, oxygen (paper Table 1 / Table 3) ---
    33: ("OH tert", 349.9409, 90),
    34: ("OH ((C,Si)2H(C,Si)) secondary", 390.2446, 89),
    35: ("OH (C4) primary", 443.8712, 87),
    36: ("OH short chain (<C5)", 488.0819, 91),
    37: ("OH (Ca) phenol", 361.4775, 88),
    38: ("(C,Si)-O-(C,Si) ether", 146.4836, 93),
    39: ("(OC2) epoxide", 820.7118, 50),
    44: ("COOH-(C)", 1080.3139, 24),
    45: ("(C)-COO-(C) ester", 636.2020, 25),
    46: ("HCOO-(C) formate ester", 642.0427, 27),
    47: ("C(c)OO- (lactone)", 1142.6119, 26),
    48: ("CON (disubstituted amide)", 1052.6072, 49),
    49: ("CONH (monosubstituted amide)", 1364.5333, 48),
    50: ("CONH2 (amide)", 1487.4109, 28),
    51: ("O=C-(Cna)2 ketone", 618.9782, 55),
    52: ("CHO-(Cna) aldehyde", 553.8090, 53),
    90: ("CHO-(Ca) aldehyde", 560.1024, 52),
    65: ("(C(a))-O(a)-(C(a)) aromatic ether", 37.1936, 92),
    72: ("NO3 (nitrate ester)", 920.3617, 14),
    74: ("O=N-O-(C) nitrite", 494.2668, 23),
    73: ("PO(O)3 (phosphate)", 1153.1344, 10),
    75: ("ONC (oxime)", 1041.0851, 30),
    76: ("C=O-O-C=O (anhydride)", 1251.2675, 12),
    79: ("O=C(O)2 (carbonate)", 879.7062, 15),
    80: ("OCN (cyanate)", 660.4645, 29),
    94: ("O-O (peroxide)", 273.1755, 32),
    96: ("(C=O-O-C=O)r (cyclic anhydride)", 2082.3288, 11),
    99: ("OCON (carbamate)", 886.7613, 1),
    100: ("N-(C=O)-N (urea)", 1045.0343, 2),
    103: ("OCOO (C,O carbonate)", 1573.3769, 34),
    # --- Main groups, sulfur / selenium (paper Table 1 / Table 3) ---
    53: ("SH-(C) thiol", 434.0811, 74),
    54: ("(C)-S-(C) sulfide", 461.5784, 75),
    55: ("(C)-S-S-(C) disulfide", 864.5074, 51),
    56: ("S(a) aromatic sulfur", 304.3321, 76),
    81: ("SCN-(C) thiocyanate", 1018.4865, 20),
    82: ("(C)-SO2-(C) sulfone", 1559.9840, 18),
    104: ("SO4 (sulfate)", 1483.1289, 35),
    105: ("SO2-N", 1506.8136, 36),
    107: ("S=O (sulfoxide)", 1379.4485, 37),
    116: ("Se", 562.1791, 46),
    # --- Main groups, nitrogen (paper Table 1 / Table 3) ---
    40: ("NH2-(C,Si)", 321.1759, 95),
    41: ("NH2-(Ca)", 441.4388, 94),
    42: ("(C,Si)-NH-(C,Si) secondary amine", 223.0992, 99),
    43: ("(C,Si)2-N-(C,Si) tertiary amine", 126.2952, 100),
    57: ("(C)-C≡N cyanide", 719.2462, 56),
    66: ("=N(a)(r5) aromatic N, five-ring", 453.3397, 97),
    67: ("=N(a)(r6) aromatic N, six-ring", 306.7139, 96),
    68: ("NO2-(C)", 866.5843, 21),
    69: ("NO2-(C(a))", 821.4141, 22),
    77: ("COCl (acid chloride)", 778.9151, 19),
    91: ("(C,Si)=N", 229.2288, 101),
    92: ("(O=C-(C)2)a ketone, aromatic side", 606.1797, 54),
    97: ("(C,Si)a-NH-(Ca,Si)a aromatic sec. amine", 201.3224, 98),
    101: ("(C,Si)2-N-(C,Si)2 quaternary amine", -109.6269, 33),
    106: ("imidazole (=CNC=NC=)", 484.6371, 3),
    108: ("(S)-C≡N", 659.7336, 38),
    109: ("N-(C=O)", 492.0707, 39),
    111: ("(N)-C≡N", 971.0365, 41),
    115: ("ON=(C,Si) (isoxazole)", 612.9506, 45),
    # --- Main groups, unsaturated carbon (paper Table 1 / Table 3) ---
    58: ("C(c)=C(c)", 475.7958, 63),
    59: ("C(c)=C(c)(C(a))", 586.1413, 60),
    60: ("-(e)C(c)=C(c)", 500.2434, 59),
    61: ("H2C(c)=C (1-ene)", 412.6276, 58),
    62: ("C(r)=C(r)", 475.9623, 61),
    63: ("C≡C", 512.2893, 62),
    64: ("HC≡C (1-yne)", 422.2307, 57),
    87: ("C=C=C (cumulated)", 664.0903, 6),
    88: ("C=C-C=C (conjugated, ring)", 957.6388, 7),
    89: ("C=C-C=C (conjugated, chain)", 928.9954, 8),
    95: ("C≡C-C≡C (conjugated)", 1218.1878, 9),
    # --- Main groups, Si / Ge / Sn / P / As / B / Al (Table 1 / Table 3) ---
    70: ("Si", 282.0181, 79),
    71: ("Si-(O)", 207.9312, 77),
    78: ("Si-(F,Cl) / B(O-)3 (borate ester)", 540.0895, 16),  # dual definition, one ID in print
    86: ("(C)2-Ge-(C)2", 347.7717, 68),
    85: ("GeCl3", 1209.2972, 13),
    83: ("(C)2-Sn-(C)2", 510.4223, 65),
    84: ("AsCl2", 1149.9670, 17),
    113: ("P (phosphine)", 428.8911, 43),
    117: ("Al", 761.6006, 47),
    # --- Main group without Part 1 definition (Table 3 value, P3 backfill) ---
    93: ("Si(C,H)2", 215.3416, 78),  # name/PR from Part 3 Table 1 (note 3.1)
    # --- Second-order corrections (paper Table 2 / Table 4, no PR) ---
    118: ("C=C-C=O", 40.4205, None),
    119: ("(C=O)-C([F,Cl]2,3)", -82.2328, None),
    120: ("(C=O)-(C([F,Cl]2,3))2", -247.8893, None),
    121: ("C-[F,Cl]3", -20.3996, None),
    122: ("(C)2C-[F,Cl]2", 15.4720, None),
    123: ("no hydrogen", -172.4201, None),
    124: ("one hydrogen", -99.8035, None),
    125: ("3/4 ring", -62.3740, None),
    126: ("five-ring", -40.0058, None),
    127: ("ortho pair (exclusive)", -27.2705, None),
    128: ("meta pair (exclusive)", -3.5075, None),
    129: ("para pair (exclusive)", 16.1061, None),
    130: ("(C=)(C)C-CC3", 25.8348, None),
    131: ("C2C-CC2", 35.8330, None),
    132: ("C3C-CC2", 51.9098, None),
    133: ("C3C-CC3", 111.8372, None),
    # --- SOA extension groups (Compernolle 2010 Table 3, no PR, no GI) ---
    301: ("OOH (hydroperoxide)", 774.75, None),
    302: ("C(=O)OOH (peroxy acid)", 1110.64, None),
    303: ("C(=O)OONO2 (PAN)", 1467.2, None),
}

#: Non-additive interaction classes of paper Table 9, keyed by the paper's
#: class letter. Values are ``(description, group IDs)``. The paper lists
#: group 37 under both A and B; the Part 3 split is applied here (37 is a
#: phenol, class B only) per method note OQ-15.
NANNOOLAL_BP_INTERACTION_CLASSES = {
    "A": ("alcohol -OH", (34, 35, 36)),
    "B": ("phenol -OH(a)", (37,)),
    "C": ("carboxylic acid -COOH", (44,)),
    "D": ("ether -O-", (38,)),
    "E": ("epoxide (OC2)", (39,)),
    "F": ("ester -COO-", (45, 46, 47)),
    "G": ("ketone -CO-", (51, 92)),
    "H": ("aldehyde -CHO", (52, 90)),
    "I": ("aromatic oxygen -O(a)-", (65,)),
    "J": ("sulfide -S(na)-", (54,)),
    "K": ("aromatic sulfur -S(a)-", (56,)),
    "L": ("thiol -SH", (53,)),
    "M": ("primary amine -NH2", (40, 41)),
    "N": ("secondary amine -NH-", (42, 97)),
    "O": ("isocyanate -OCN", (80,)),
    "P": ("cyanide -CN", (57,)),
    "Q": ("nitro group NO2-(C(a))", (69,)),  # P1 Table 9 prints "Nitrate (69)"
    "R": ("aromatic N in five-ring", (66,)),
    "S": ("aromatic N in six-ring", (67,)),
}

#: Group-interaction contributions [K] of paper Table 5, keyed by the
#: sorted ``(class_i, class_j)`` letter pair (interactions are symmetric,
#: ``C_i-j = C_j-i``). Only the 76 pairs with a regressed parameter are
#: listed; pairs printed without a value (rows 153 B-G and 211 R-R, see
#: :data:`NANNOOLAL_BP_INTERACTION_IDS`) have no parameter and must be
#: treated as missing, not as zero.
NANNOOLAL_BP_INTERACTIONS = {
    ("A", "A"): 291.7985,
    ("A", "C"): 146.7286,
    ("A", "D"): 135.3991,
    ("A", "E"): 226.4980,
    ("A", "F"): 211.6814,
    ("A", "G"): 46.3754,
    ("A", "I"): 435.0923,
    ("A", "J"): -74.0193,
    ("A", "L"): 38.6974,
    ("A", "M"): 314.6126,
    ("A", "N"): 286.9698,
    ("A", "P"): 306.3979,
    ("A", "S"): 1334.6747,
    ("B", "B"): 288.6155,
    ("B", "C"): -1477.9671,
    ("B", "D"): 130.3742,
    ("B", "F"): -1184.9784,
    ("B", "H"): 43.9722,
    ("B", "M"): 797.4327,
    ("B", "Q"): -1048.1236,
    ("B", "S"): -614.3624,
    ("C", "C"): 117.2044,
    ("C", "D"): 612.8821,
    ("C", "F"): -183.2986,
    ("C", "G"): -55.9871,
    ("D", "D"): 91.4997,
    ("D", "E"): 178.7845,
    ("D", "F"): 322.5671,
    ("D", "G"): 15.6980,
    ("D", "H"): 17.0400,
    ("D", "I"): 329.0050,
    ("D", "J"): 394.5505,
    ("D", "P"): 293.5974,
    ("D", "Q"): 963.6518,
    ("E", "E"): 1006.3880,
    ("E", "H"): 163.5475,
    ("F", "F"): 431.0990,
    ("F", "G"): 22.5208,
    ("F", "I"): 707.9404,
    ("F", "P"): 517.0677,
    ("F", "Q"): -205.6165,
    ("G", "G"): -303.9653,
    ("G", "H"): -391.3690,
    ("G", "I"): 176.5481,
    ("G", "K"): 381.0107,
    ("G", "P"): -574.2230,
    ("G", "Q"): -3628.9026,
    ("G", "S"): 124.1943,
    ("H", "H"): 582.1763,
    ("H", "I"): 674.6858,
    ("H", "K"): 397.5750,
    ("H", "Q"): 140.9644,
    ("I", "R"): -888.6123,
    ("J", "J"): -11.9406,
    ("K", "P"): -101.2319,
    ("K", "R"): -348.7400,
    ("L", "L"): 217.6360,
    ("L", "AF"): 501.2778,  # printed "L-AF" (row 171), questionable: note OQ-15
    ("D", "M"): 124.3549,
    ("F", "M"): 182.6291,
    ("I", "M"): 395.4093,
    ("J", "M"): -562.3061,
    ("M", "M"): 174.0258,
    ("M", "N"): 510.3473,
    ("M", "Q"): 663.8009,
    ("M", "S"): 27.2735,
    ("D", "N"): 101.8475,
    ("F", "N"): 317.0200,
    ("G", "N"): -215.3532,
    ("N", "N"): 239.8076,
    ("N", "S"): 758.9855,
    ("O", "O"): -356.5017,
    ("O", "Q"): -263.0807,
    ("P", "S"): -370.9729,
    ("Q", "Q"): 65.1432,
    ("S", "S"): -271.9449,
}

#: Paper Table 5 row IDs (135-212) for every tabulated interaction pair,
#: including the two rows printed without a value (153 and 211). The row
#: IDs are property-independent across the Nannoolal paper series.
NANNOOLAL_BP_INTERACTION_IDS = {
    ("A", "A"): 135,
    ("A", "C"): 139,
    ("A", "D"): 140,
    ("A", "E"): 141,
    ("A", "F"): 142,
    ("A", "G"): 143,
    ("A", "I"): 146,
    ("A", "J"): 144,
    ("A", "L"): 138,
    ("A", "M"): 136,
    ("A", "N"): 137,
    ("A", "P"): 145,
    ("A", "S"): 147,
    ("B", "B"): 148,
    ("B", "C"): 150,
    ("B", "D"): 151,
    ("B", "F"): 152,
    ("B", "G"): 153,  # no parameter printed (do not estimate)
    ("B", "H"): 154,
    ("B", "M"): 149,
    ("B", "Q"): 155,
    ("B", "S"): 156,
    ("C", "C"): 172,
    ("C", "D"): 173,
    ("C", "F"): 174,
    ("C", "G"): 175,
    ("D", "D"): 178,
    ("D", "E"): 179,
    ("D", "F"): 180,
    ("D", "G"): 181,
    ("D", "H"): 182,
    ("D", "I"): 186,
    ("D", "J"): 183,
    ("D", "P"): 185,
    ("D", "Q"): 184,
    ("E", "E"): 187,
    ("E", "H"): 188,
    ("F", "F"): 189,
    ("F", "G"): 190,
    ("F", "I"): 193,
    ("F", "P"): 192,
    ("F", "Q"): 191,
    ("G", "G"): 194,
    ("G", "H"): 195,
    ("G", "I"): 199,
    ("G", "K"): 197,
    ("G", "P"): 198,
    ("G", "Q"): 196,
    ("G", "S"): 200,
    ("H", "H"): 201,
    ("H", "I"): 204,
    ("H", "K"): 203,
    ("H", "Q"): 202,
    ("I", "R"): 210,
    ("J", "J"): 205,
    ("K", "P"): 207,
    ("K", "R"): 208,
    ("L", "L"): 170,
    ("L", "AF"): 171,  # printed "L-AF", questionable: note OQ-15
    ("D", "M"): 159,
    ("F", "M"): 160,
    ("I", "M"): 163,
    ("J", "M"): 161,
    ("M", "M"): 157,
    ("M", "N"): 158,
    ("M", "Q"): 162,
    ("M", "S"): 164,
    ("D", "N"): 166,
    ("F", "N"): 167,
    ("G", "N"): 168,
    ("N", "N"): 165,
    ("N", "S"): 169,
    ("O", "O"): 176,
    ("O", "Q"): 177,
    ("P", "S"): 209,
    ("Q", "Q"): 206,
    ("R", "R"): 211,  # no parameter printed (do not estimate)
    ("S", "S"): 212,
}

#: SOA extension groups (aeroprops slot IDs 301-303), keyed like
#: :data:`NANNOOLAL_BP_GROUPS`. Values are ``(name, source note)``. The
#: boiling-point contributions themselves are merged into
#: :data:`NANNOOLAL_BP_GROUPS`; membership in this table marks a group as
#: having no interaction parameters (the no-GI branch of the estimator).
NANNOOLAL_BP_EXTENSIONS = {
    301: (
        "OOH (hydroperoxide)",
        "Compernolle et al. (2010), Atmos. Chem. Phys. 10, 6271-6282, "
        "Table 3, Part 1 Tb column: 774.75 K, sigma_Tb = 1.3 K; anchor Tb "
        "extrapolated (their Table 2); no group-interaction parameters "
        "derivable.",
    ),
    302: (
        "C(=O)OOH (peroxy acid)",
        "Compernolle et al. (2010), Table 3, Part 1 Tb column: 1110.64 K, "
        "sigma_Tb = 12.3 K; anchor Tb extrapolated; no group-interaction "
        "parameters derivable.",
    ),
    303: (
        "C(=O)OONO2 (PAN)",
        "Compernolle et al. (2010), Table 3, Part 1 Tb column: 1467.2 K, "
        "single-component fit (no sigma given); anchor Tb extrapolated; no "
        "group-interaction parameters derivable.",
    ),
}


def nannoolal_bp_group_sum(counts):
    """Sum of the first-order and second-order group contributions [K].

    Parameters
    ----------
    counts : dict
        Group counts keyed by group ID as produced by
        :func:`aeroprops.boiling_point.nannoolal_fragmentation`.

    Returns
    -------
    total : float or None
        ``sum(N_i * C_i)`` over all groups in `counts`, or ``None`` if a
        key is not a valid group ID.

    Examples
    --------
    Paper Table 17a (3,3,4,4-tetramethylhexane): the printed group sum is
    2153.745 K (6 CH3, 2 chain CH2, 2 quaternary chain C, 1 steric
    correction 133):

    >>> round(nannoolal_bp_group_sum({1: 6, 4: 2, 6: 2, 133: 1}), 3)
    2153.745
    >>> nannoolal_bp_group_sum({999: 1}) is None
    True
    """
    total = 0.0
    try:
        for gid, n in counts.items():
            total += n * NANNOOLAL_BP_GROUPS[gid][1]
    except (KeyError, TypeError):
        return None
    return total


def nannoolal_bp_interaction_counts(counts):
    """Count the non-additive interaction groups per Table 9 class.

    Parameters
    ----------
    counts : dict
        Group counts keyed by group ID.

    Returns
    -------
    class_counts : dict or None
        Mapping of class letter to the number of interaction groups of
        that class, or ``None`` on invalid input.

    Examples
    --------
    Di-isopropanolamine (paper Table 17b): two alcohol OH (class A) and
    one secondary amine (class N), m = 3:

    >>> nannoolal_bp_interaction_counts({1: 2, 7: 4, 34: 2, 42: 1}) == \
{'A': 2, 'N': 1}
    True
    """
    group_to_class = {}
    for letter, (_desc, gids) in NANNOOLAL_BP_INTERACTION_CLASSES.items():
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


