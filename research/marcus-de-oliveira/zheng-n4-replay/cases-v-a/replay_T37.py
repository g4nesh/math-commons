"""Independent exact replay of Case T_37 = {0,5,6,7,9,12} (group V) of D. Zheng's proof of the
Marcus-de Oliveira conjecture for n = 4.

Proof document replayed: MOC4_SUPP/proofs/manin_transfer_20261021/independent_moc4_T36_uniform_20261022.pdf
(6 pages; the document calls the set "T36").  Every object is rebuilt from the paper's
definitions by own code (replay_T35_common.py, replay_T35_common_poly.py).  Nothing from the
supplement is imported or executed.  If MOC4_SUPP is set, the PDF is hashed and four JSON data
files are read for comparison only.

Run:   MOC4_SUPP=/path/to/supplement python replay_T37.py [--no-numeric]
Exit status 0 iff every check passes.  Checks tagged EXACT use integer/rational/polynomial
arithmetic only; checks tagged NUMERIC are floating-point falsification attempts.
"""
from math import comb
import os
import random
import sys
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import replay_T35_common as M            # noqa: E402
import replay_T35_common_poly as PP      # noqa: E402
from replay_T35_common_poly import PR, RF  # noqa: E402

NUMERIC = "--no-numeric" not in sys.argv
T = (0, 5, 6, 7, 9, 12)
DOC = "independent_moc4_T36_uniform_20261022.pdf"
# (2): g = q f2 + p f29 + j f1 + c f45 + e f5 + d f35 + b f0 + a f7
ORDER35 = [2, 29, 1, 45, 5, 35, 0, 7]
NAMES35 = list("qpjcedba")
GEN35 = dict(zip(NAMES35, ORDER35))
# (3): g = p f0 + r f39 + a f5 + b f2 + c f29 + q f7 + d f40 + s f45
ORDER40 = [0, 39, 5, 2, 29, 7, 40, 45]
NAMES40 = list("prabcqds")
GEN40 = dict(zip(NAMES40, ORDER40))

ck = M.Checker("Case T_37 = {0,5,6,7,9,12}: replay of independent_moc4_T36_uniform_20261022.pdf")
M.print_versions(ck)
M.check_document_identity(ck, DOC)

M.check_foundations(ck)
M.check_radial_necessity(ck)
M.check_hessian_necessity(ck)
M.crosscheck_with_supplement_data(ck)

M.check_contact_set(ck, 37, T, "T36")
gens = M.check_generators(ck, T, [0, 1, 2, 5, 7, 29, 35, 39, 40, 45], [35, 40])

# ------------------------------------------------------------------------------------------
ck.section("Section 1 of the document: the three-way reduction, as identities on the contact nullspace")
NB = M.nullspace([M.P[s] for s in T], 14)
NV = [M.values(n) for n in NB]                       # the 24 gaps of each basis vector
ck.check("EXACT", "the contact nullspace {c : (Pc)_T = 0} has dimension 8", len(NB) == 8 and M.rank(NV) == 8)


def lin_check(coeff):
    """coeff: {generator j: function of the 24 gaps}; True iff c = sum_j coeff_j(c) f_j on the nullspace."""
    for v in NV:
        tot = [Fr(0)] * 24
        for jx, fn in coeff.items():
            cj = fn(v)
            tot = [x + cj * y for x, y in zip(tot, M.F_VALUES[jx])]
        if tot != [Fr(x) for x in v]:
            return False
    return True


d35 = lambda v: v[3] - v[1]
d40 = lambda v: v[14] - v[16]
ck.check("EXACT", "delta35 + delta40 = -g17 on the nullspace, with delta35 = g3 - g1, delta40 = g14 - g16",
         all(d35(v) + d40(v) == -v[17] for v in NV))
ck.check("EXACT", "(1): c = g18 f0 - delta40 f1 + g13 f2 + g11 f5 + g14 f7 + g2 f29 - delta35 f39 + g3 f45 for every c in "
         "the nullspace; all eight generators are minor vectors",
         lin_check({0: lambda v: v[18], 1: lambda v: -d40(v), 2: lambda v: v[13], 5: lambda v: v[11],
                    7: lambda v: v[14], 29: lambda v: v[2], 39: lambda v: -d35(v), 45: lambda v: v[3]})
         and all(M.F_IS_MINOR[j] for j in (0, 1, 2, 5, 7, 29, 39, 45)))
ck.check("EXACT", "(2): c = g15 f2 + g4 f29 + g17 f1 + g1 f45 + g11 f5 + delta35 f35 + g18 f0 + g14 f7 for every c in the "
         "nullspace; f35 is the only non-minor generator in this list",
         lin_check({2: lambda v: v[15], 29: lambda v: v[4], 1: lambda v: v[17], 45: lambda v: v[1],
                    5: lambda v: v[11], 35: d35, 0: lambda v: v[18], 7: lambda v: v[14]})
         and [j for j in ORDER35 if not M.F_IS_MINOR[j]] == [35] and M.rank([M.F_VALUES[j] for j in ORDER35]) == 8)
ck.check("EXACT", "(3): c = g23 f0 + g17 f39 + g11 f5 + g13 f2 + g2 f29 + g16 f7 + delta40 f40 + g3 f45 for every c in the "
         "nullspace; f40 is the only non-minor generator in this list",
         lin_check({0: lambda v: v[23], 39: lambda v: v[17], 5: lambda v: v[11], 2: lambda v: v[13],
                    29: lambda v: v[2], 7: lambda v: v[16], 40: d40, 45: lambda v: v[3]})
         and [j for j in ORDER40 if not M.F_IS_MINOR[j]] == [40] and M.rank([M.F_VALUES[j] for j in ORDER40]) == 8)
ck.info("consequence (PROOF-READ): for c in N_T all gaps are >= 0, so delta35 + delta40 <= 0; if both are <= 0, (1) is a "
        "nonnegative combination of minor vectors; if delta35 > 0, (2) has d = delta35 > 0 and nonnegative gaps as the "
        "other coefficients; if delta40 > 0, the same for (3).  No list of extreme rays is needed.")
M.check_linear_identities(ck, [([35, 40], [0, 7, 2, 29, 45]), ([35, 39], [1, 2, 29, 45]), ([40, 1], [0, 7, 39])],
                          "the three 'generating identities' of p. 1")
M.check_minor_leibniz(ck, [j for j in gens if M.F_IS_MINOR[j]])
rng = random.Random(3705)
okneg = True
for _ in range(4):
    O = M.random_rational_SO4(rng)
    Om = [[-row[0]] + list(row[1:]) for row in O]            # first column negated: det = -1
    okneg &= M.det4(Om) == -1
    for jx in gens:
        okneg &= M.phi(M.F_VALUES[jx], Om) == -M.phi(M.F_VALUES[jx], O)
ck.check("EXACT", "reading of the document's 'every O in O(4)': negating one column of O in SO(4) changes the sign of "
         "every raw Leibniz sum (each product contains exactly one entry of that column; confirmed for the ten "
         "generators on 4 random rational matrices), so nonnegativity can only be meant on SO(4), or for Phi/det "
         "on O(4); the case statement and the cited inequalities concern SO(4)", okneg)

# ------------------------------------------------------------------------------------------
ck.section("Section 2: radial necessity and the two guards")
rd35 = PP.radial_data(ORDER35, NAMES35)
n35, G35 = PP.check_radial(ck, rd35, "branch delta35 > 0", 983, "independent_moc4_T36_bad35_full_radial_20261022.json")
rd40 = PP.radial_data(ORDER40, NAMES40)
n40, G40 = PP.check_radial(ck, rd40, "branch delta40 > 0", 1385, "independent_moc4_T36_bad40_full_radial_20261022.json",
                           alias={"s": "sa"})
assert n35 == ["q", "p", "j", "c", "e", "b", "a"] and n40 == ["p", "r", "a", "b", "c", "q", "s"]
c35 = PR.terms(rd35["Cstar"])
exd = lambda R, k: tuple(k if nm == "d" else 0 for nm in R.names)
ck.check("EXACT", "branch (2): all 13 coefficients of C* are negative (each equals -1), among them -d^3: C* <= -d^3 < 0 "
         "for d > 0, so C* + hQ = 0 forces Q != 0", len(c35) == 13 and all(v == -1 for v in c35.values())
         and c35.get(exd(rd35["R"], 3)) == -1)
R9 = rd40["R"]
W9 = {nm: RF(R9, R9[nm]) for nm in NAMES40}
p, r, a, b, c, q, d, s = [W9[nm] for nm in NAMES40]
Mm, Uu = b + c, d + r
T0 = b * (d + q) + c * (d + p) + (d + p) * (d + q)
Cst = RF(R9, rd40["Cstar"])
last = Mm * Uu * T0 - d * d * b * c
lt = PR.terms(last.num)
bd3 = tuple({"b": 1, "d": 3}.get(nm, 0) for nm in R9.names)
cd3 = tuple({"c": 1, "d": 3}.get(nm, 0) for nm in R9.names)
ck.check("EXACT", "branch (3): C* = bc(U - a) + U T0 and M C* = bc(dM + d^2 - aM) + rMbc + (M U T0 - d^2 bc), with "
         "M = b + c, U = d + r, T0 = b(d+q) + c(d+p) + (d+p)(d+q); the last bracket has only nonnegative "
         "coefficients and contains b d^3 and c d^3 (identity (5))",
         Cst == b * c * (Uu - a) + Uu * T0
         and Mm * Cst == b * c * (d * Mm + d * d - a * Mm) + r * Mm * b * c + last
         and all(v >= 0 for v in lt.values()) and lt.get(bd3, 0) >= 1 and lt.get(cd3, 0) >= 1)
ck.check("EXACT", "branch (3) at M = b + c = 0: C* = (d+r)(d+p)(d+q) > 0",
         R9.compose({ex: v for ex, v in PR.terms(rd40["Cstar"]).items()},
                    [R9.zero if nm in ("b", "c") else R9[nm] for nm in R9.names])
         == ((d + r) * (d + p) * (d + q)).num)
ck.info("consequence (PROOF-READ): if M > 0 and a < d + d^2/M then dM + d^2 - aM > 0 and M C* > 0; the energy bound "
        "J <= 1 + 1/(b + c) (below) gives a < d + d^2/M whenever a < J.  So C* > 0, hence Q != 0, on the hard set")

# ------------------------------------------------------------------------------------------
ck.section("Section 3: the branch delta35 > 0.  Inequalities used (CITED) and their roles")
ROLES35 = [
    ("SOC", (GEN35["c"], GEN35["j"], 35), "cj >= d^2: f45, f1"),
    ("four", (GEN35["a"], GEN35["b"], GEN35["e"], 35), "K >= 0, i.e. ab(e - d) >= d^2(a+b+e+d): f7, f0; f5"),
    ("old-H", (GEN35["e"], GEN35["c"], 35), "c >= d + d^2/e: t = e on f5, main = c on f45"),
]
M.check_roles(ck, ROLES35)
M.check_coefficient_forms(ck)
SRC24 = (67, 55, 53, 54, 59, 2)
DST24 = (GEN35["a"], GEN35["b"], GEN35["e"], 35, GEN35["c"], GEN35["j"])
ks = M.find_group_elements(SRC24, DST24)
kx = M.find_group_elements(SRC24, (GEN35["b"], GEN35["a"], GEN35["e"], 35, GEN35["c"], GEN35["j"]))
g24 = M.G_ELEMENTS[ks[0]] if ks else None
rng = random.Random(3701)
cov = bool(ks)
if ks:
    for _ in range(5):
        O = M.random_rational_SO4(rng)
        O2 = M.transfer_to_SO4(g24, O)
        cov &= M.is_SO4(O2) and all(M.phi(M.F_VALUES[x], O) == M.phi(M.F_VALUES[y], O2) for x, y in zip(SRC24, DST24))
ck.check("EXACT", "core action of the document: a group element maps (f67, f55, f53, f54, f59, f2) to "
         "(f7, f0, f5, f35, f45, f1) = (a, b, e, d, c, j) position by position; it is unique, equals "
         "(rho = (1,0,3,2), rho' = (2,0,1,3), no inverse), i.e. the document's 'row (1,2,0,3), column (1,0,3,2), no "
         "transpose' read as row = rho'^(-1), column = rho; Phi_{g.f}(O'') = Phi_f(O) exactly on 5 random rational "
         "O in SO(4); a second element realises the same map with f7 and f0 exchanged",
         ks and len(ks) == 1 and g24 == ((1, 0, 3, 2), (2, 0, 1, 3), 0)
         and M.inverse(g24[1]) == (1, 2, 0, 3) and cov and len(kx) == 1)
ck.info("CITED (audit foundations/so4-atoms-part2a.md, item F1; document root_moc4_T24_uniform p. 2): for the cost "
        "q f0 + p f2 + j f59 + c f6 + e f53 + d f54 + b f55 + a f67 with d > 0, Phi >= 0 on SO(4) if  (5) ab > d^2, e > 0, "
        "j >= d^2/e - d(ab - d^2)/((a+d)(b+d));  (6) ab < d^2, a != b, e > 0, j >= d^2/e + b(d - w)^2/(b - w)^2, w = sqrt(ab);  "
        "double-SOC: ab >= d^2, ej >= d^2.  Transferred by the core action: the central coefficient j of f59 becomes "
        "the coefficient c of f45, and (a, b, e, d) keep their names")
st = M.find_group_elements((7, 0, 2, 29, 1, 45, 5, 35), (0, 7, 29, 2, 1, 45, 5, 35))
ck.check("EXACT", "a group element exchanges f7 <-> f0 and f2 <-> f29 and fixes f1, f45, f5, f35 (the document's "
         "stabilizer: row (2,3,0,1), column (3,2,1,0)): it exchanges (a, b, q, p) with (b, a, p, q) in (2)",
         len(st) == 1 and M.G_ELEMENTS[st[0]] == ((3, 2, 1, 0), (2, 3, 0, 1), 0))
R7 = PR(n35)
V7 = {nm: RF(R7, R7[nm]) for nm in n35}
q, p, j, c, e, b, a = [V7[nm] for nm in n35]
G35p = R7.from_terms(G35)
swap = {"a": "b", "b": "a", "q": "p", "p": "q"}
G35s = R7.compose(G35, [R7[swap.get(nm, nm)] for nm in n35])
ck.check("EXACT", "G35 (at d = 1) is invariant under (a, b, q, p) -> (b, a, p, q)", G35s == G35p)
IC = n35.index("c")
kc = PP.coeffs_in(G35, IC)
k0, k1, k2 = kc[0], kc[1], kc[2]
DEG35 = PP.term_dict_degrees(G35, 7)
ck.check("EXACT", "G35 is quadratic in c: G35 = k0 + k1 c + k2 c^2 (%d + %d + %d terms); degrees %s in (q,p,j,c,e,b,a)"
         % (len(k0), len(k1), len(k2), DEG35), sorted(kc) == [0, 1, 2])
Cbar = (a + 1) * (b + 1) + a * q + b * p + (1 + e) * (p * q + p + q) + e
K0p = R7.from_terms(k0)
R90 = K0p // Cbar.num
Kc = (a * b - 1) * e - (a + 1) * (b + 1)
ck.check("EXACT", "k0 = Cbar R90 with Cbar = (a+1)(b+1) + aq + bp + (1+e)(pq+p+q) + e (all coefficients positive, constant "
         "term 1) and a polynomial R90 with %d terms" % len(R90),
         R90 * Cbar.num == K0p and PP.nonneg(Cbar.num) and len(R90) == 90)
ck.check("EXACT", "K = (ab-1)e - (a+1)(b+1) equals ab(e-1) - (a+b+e+1): K >= 0 is the four-coefficient condition",
         (Kc - (a * b * (e - 1) - (a + b + e + 1))).is_zero())

ck.section("Section 3: the strict constant coefficient (8)")
diff8 = PR.terms(R90 + ((a + 1) * (b + 1) * Kc).num)
RE = PR(["A", "T", "V", "U", "q", "p", "j"])
A_, T_, V_, U_, q_, p_, j_ = RE.gens
re_ = lambda z: RF(RE, z)
DOC8 = {(0, 0): 280, (0, 1): 404, (1, 0): 280, (1, 1): 404}
for exch in (0, 1):
    for above in (0, 1):
        Rr = re_(A_ + 1)
        z = re_(T_ + 1) if above else re_(T_) / re_(T_ + 1)
        aa, bb = z * Rr, z / Rr
        if exch:
            aa, bb = bb, aa
        ee = (aa + 1) * (bb + 1) / (aa * bb - 1) * re_(V_) / re_(V_ + 1) if above else re_(V_)
        ch = {"a": aa, "b": bb, "e": ee, "q": re_(q_), "p": re_(p_), "j": re_(j_), "c": re_(RE.one)}
        okc = PP.chart_denominators_positive(RE, ch, ["T"] if above else [])
        okc &= (aa * bb == z * z)
        if above:        # -K = (a+1)(b+1)/(1+V) > 0
            okc &= (-((aa * bb - 1) * ee - (aa + 1) * (bb + 1)) == (aa + 1) * (bb + 1) / re_(V_ + 1))
        Z, dens = PP.cleared(RE, diff8, n35, ch)
        stz = PP.poly_stats(Z)
        ck.check("EXACT", "(8), chart %s, ab %s 1 (z = %s, e = %s): chart valid; the cleared difference "
                 "R90 - (a+1)(b+1)(-K) has %d terms (document: %d), all coefficients >= 0"
                 % ("b = zR, a = z/R" if exch else "a = zR, b = z/R", "><"[1 - above],
                    "1+T" if above else "T/(1+T)", "Emax V/(1+V)" if above else "V", stz[0], DOC8[(exch, above)]),
                 okc and stz[1] == 0 and stz[0] == DOC8[(exch, above)])
ck.info("consequence: k0 >= Cbar (a+1)(b+1)(-K) > 0 whenever K < 0, on the closure of the four charts (all a, b >= 0)")

ck.section("Section 3: the seven intervals (Bernstein tables)")
RS = PR(["A", "T", "V", "U", "q", "p"])
A_, T_, V_, U_, q_, p_ = RS.gens
rs = lambda z: RF(RS, z)
b0 = rs(T_) / rs(T_ + 1)
b1 = rs(T_ + 1)


def emax(x, y):
    return (x + 1) * (y + 1) / (x * y - 1)


def coherent(x, y, ev):
    return 1 / ev - (x * y - 1) / ((x + 1) * (y + 1))


SEVEN = []
aa, bb, ee = rs(A_), rs(T_) / rs(A_ * (T_ + 1)), rs(V_)
jj = rs(V_ + U_)
SEVEN.append(("I-", aa, bb, ee, jj, 1 / jj, "SOC 1/j", ["A"], (969, 1704, 2364)))
aa, bb = rs(A_), rs(T_ + 1) / rs(A_)
ee = emax(aa, bb) * rs(V_) / rs(V_ + 1)
jj = ee + rs(U_)
SEVEN.append(("I+", aa, bb, ee, jj, 1 / jj, "SOC 1/j", ["A", "T"], (2643, 5745, 10406)))
jj = ee * rs(U_) / rs(U_ + 1)
Em = emax(aa, bb)
# the document's fraction for 1/e - 1/Emax is not reduced: (Emax - e)/(e Emax) over the two reduced fractions
thrII = RF(RS, Em.num * ee.den - ee.num * Em.den, ee.num * Em.num, reduce=False)
SEVEN.append(("II", aa, bb, ee, jj, thrII, "coherent 1/e - 1/Emax", ["A", "T", "V"], (4414, 9616, 16507)))
aa, bb, ee = b0 * rs(V_) / rs(V_ + 1), b0, rs(A_)
jj = rs(A_ * U_) / rs((A_ + 1) * (U_ + 1))
SEVEN.append(("IIIa-", aa, bb, ee, jj, 1 + 1 / ee, "old-H 1 + 1/e", ["A"], (3779, 4535, 4535)))
jj = rs(A_) / rs(A_ + 1) + rs(A_ * A_ * U_) / rs((A_ + 1) * (U_ + 1))
SEVEN.append(("IIIa+", aa, bb, ee, jj, 1 / jj, "SOC 1/j", ["A"], (4284, 7308, 9756)))
aa, bb, ee = rs(V_) / (b1 * rs(V_ + 1)), b1, rs(A_)
jj = rs(A_) * b1 / (rs(A_) + b1) + rs(A_ * A_ * U_) / ((rs(A_) + b1) * rs(U_ + 1))
SEVEN.append(("IIIb+", aa, bb, ee, jj, 1 / jj, "SOC 1/j", ["A"], (4428, 8448, 12396)))
ww = rs(V_) / rs(V_ + 1)
aa = ww * ww / b1
jj = rs(A_) * b1 * rs(U_) / ((rs(A_) + b1) * rs(U_ + 1))
SEVEN.append(("IIIb-", aa, bb, ee, jj, 1 / ee + bb * (1 - ww) ** 2 / ((bb - ww) ** 2), "parallel 1/e + b(1-w)^2/(b-w)^2",
              ["A"], (7028, 13410, 18874)))
DEGK = [0 if nm == "c" else dg for nm, dg in zip(n35, DEG35)]
SEVEN_CH = {}


def pos_rf(R, x, posv):
    return PP.positive_on_domain(R, x.num, posv) and PP.positive_on_domain(R, x.den, posv)


def nn_rf(R, x, posv):
    return PP.nonneg(x.num) and PP.positive_on_domain(R, x.den, posv)


for name, aa, bb, ee, jj, thr, tname, posv, doc in SEVEN:
    ch = {"a": aa, "b": bb, "e": ee, "j": jj, "q": rs(q_), "p": rs(p_), "c": rs(RS.one)}
    SEVEN_CH[name] = (ch, thr, posv)
    okc = PP.chart_denominators_positive(RS, ch, posv)
    # J = N/D with N > 0 on the domain and D >= 0 (D = 0 means J = +infinity, e.g. j = 0 in the SOC rows)
    okc &= PP.nonneg(thr.den) and PP.positive_on_domain(RS, thr.num, posv)
    ab = aa * bb
    if name == "I-":
        okc &= ab == rs(T_) / rs(T_ + 1) and (jj - ee) == rs(U_)
    elif name in ("I+", "II"):
        okc &= ab == rs(T_ + 1) and (1 - ee / emax(aa, bb)) == 1 / rs(V_ + 1)
        if name == "I+":
            okc &= (jj - ee) == rs(U_)
        else:
            okc &= (ee - jj) == ee / rs(U_ + 1) and thr == coherent(aa, bb, ee)
    elif name.startswith("IIIa"):
        okc &= bb == rs(T_) / rs(T_ + 1) and (bb - aa) == bb / rs(V_ + 1)       # a <= b < 1
        cut = ee / (1 + ee)
        if name == "IIIa-":
            okc &= (cut - jj) == cut / rs(U_ + 1)
        else:
            okc &= (jj - cut) == (ee - cut) * rs(U_) / rs(U_ + 1) and (ee - jj) == (ee - cut) / rs(U_ + 1)
    else:
        okc &= bb == rs(T_ + 1)                                               # b >= 1
        if name == "IIIb+":
            okc &= (1 - ab) == 1 / rs(V_ + 1)                                 # ab = V/(1+V) < 1
        else:
            okc &= ab == ww * ww and pos_rf(RS, bb - ww, posv) and pos_rf(RS, 1 - ww, posv)   # w = sqrt(ab) < 1 <= b
        cut = ee * bb / (ee + bb)
        if name == "IIIb+":
            okc &= (jj - cut) == (ee - cut) * rs(U_) / rs(U_ + 1) and (ee - jj) == (ee - cut) / rs(U_ + 1)
        else:
            okc &= (cut - jj) == cut / rs(U_ + 1)
    K0, mid, end, W = PP.bernstein(RS, (k0, k1, k2), n35, ch, thr, DEGK)
    s0, s1, s2 = PP.poly_stats(K0), PP.poly_stats(mid), PP.poly_stats(end)
    ck.check("EXACT", "row %s, threshold J = %s: chart valid (nonnegative coefficients, positive denominators for "
             "%s > 0; the stated ranges of ab, e and j hold identically); the three Bernstein numerators have "
             "%d, %d, %d terms (document: %d, %d, %d), all coefficients >= 0"
             % (name, tname, ", ".join(posv), s0[0], s1[0], s2[0], doc[0], doc[1], doc[2]),
             okc and s0[1] == 0 and s1[1] == 0 and s2[1] == 0 and (s0[0], s1[0], s2[0]) == doc)
thr_red = coherent(SEVEN[2][1], SEVEN[2][2], SEVEN[2][3])
K0, mid, end, W = PP.bernstein(RS, (k0, k1, k2), n35, SEVEN_CH["II"][0], thr_red, DEGK)
ck.check("EXACT", "row II again with the threshold as a reduced fraction N/D = AT/(V(A+1)(A+T+1)): the middle and "
         "endpoint numerators (%d and %d terms) also have only nonnegative coefficients (the document's fraction "
         "carries the extra positive factor (A+1)(A+T+1))" % (len(mid), len(end)),
         PP.nonneg(mid) and PP.nonneg(end) and thr_red == SEVEN[2][5])
ie, ij = n35.index("e"), n35.index("j")
g0 = PP.set_var(PP.set_var(G35, ie, 0), ij, 0)
ck.check("EXACT", "closed layer e = j = 0: G35 = [(1+a+b)c + (a+1)(b+1)](p+a+1)^2 (q+b+1)^2 > 0",
         R7.from_terms(g0) == (((1 + a + b) * c + (a + 1) * (b + 1)) * (p + a + 1) ** 2 * (q + b + 1) ** 2).num)
RW = PR(["b", "w"])
bw, wv = RF(RW, RW["b"]), RF(RW, RW["w"])
par = bw * (1 - wv) ** 2 / ((bw - wv) ** 2)
ck.check("EXACT", "the parallel term b(1-w)^2/(b-w)^2: it equals 1/b at w = 0 (chain threshold 1/e + 1/b), equals 1 at "
         "b = 1 (old-H threshold 1 + 1/e), and 1/b - b(1-w)^2/(b-w)^2 = w(b-1)(2b - w - bw)/(b(b-w)^2) >= 0 for "
         "b >= 1 > w >= 0, so for b >= 1 it never exceeds the chain threshold",
         (1 / bw - par) == wv * (bw - 1) * (2 * bw - wv - bw * wv) / (bw * (bw - wv) ** 2)
         and RF(RW, RW.compose(PR.terms(par.num), [RW["b"], RW.zero]),
                RW.compose(PR.terms(par.den), [RW["b"], RW.zero])) == 1 / bw
         and RW.compose(PR.terms(par.num), [RW.one, RW["w"]]) == RW.compose(PR.terms(par.den), [RW.one, RW["w"]]))

ck.check("EXACT", "remark of p. 4: the chain tuples (45;5,0;35) and (45;5,7;35) (centre f45, leaves f5 and f0, resp. f7) "
         "lie in the cited chain orbit, so with a single zero leaf the limit 1/e + 1/b (resp. 1/e + 1/a) of the "
         "parallel threshold is also the ordinary chain condition",
         M.role_group_element("chain", (45, 5, 0, 35)) is not None
         and M.role_group_element("chain", (45, 5, 7, 35)) is not None)

# ------------------------------------------------------------------------------------------
ck.section("Section 4: the branch delta40 > 0.  Inequalities used (CITED) and their roles")
ROLES40 = [
    ("old-H", (GEN40["b"], GEN40["a"], 40), "old-H atom of leaf 1: t on f2 (b), main on f5 (a)"),
    ("old-H", (GEN40["c"], GEN40["a"], 40), "old-H atom of leaf 2: t on f29 (c), main on f5 (a)"),
    ("chain", (GEN40["a"], GEN40["p"], GEN40["b"], 40), "chain of leaf 1: centre f5 (a), leaves f0 (p), f2 (b)"),
    ("chain", (GEN40["a"], GEN40["c"], GEN40["q"], 40), "chain of leaf 2: centre f5 (a), leaves f29 (c), f7 (q)"),
    ("SOC", (GEN40["a"], GEN40["s"], 40), "as >= d^2: f5, f45"),
]
M.check_roles(ck, ROLES40)
k7 = M.find_group_elements((4, 2, 33, 5, 6, 46, 45), (0, 5, 2, 29, 7, 40, 39))
ksoc = M.find_group_elements((1, 9, 22), (5, 45, 40))
src_ok = all(M.role_group_element(f, t) is not None for f, t in
             (("old-H", (33, 2, 46)), ("old-H", (5, 2, 46)), ("chain", (2, 4, 33, 46)), ("chain", (2, 5, 6, 46))))
ck.check("EXACT", "the document's transposed core action exists and is unique: (f4, f2, f33, f5, f6, f46, f45) -> "
         "(f0, f5, f2, f29, f7, f40, f39) by (rho = (3,2,1,0), rho' = (2,3,0,1), inverse); its source tuples "
         "old-H (33,2;46), (5,2;46) and chain (2;4,33;46), (2;5,6;46) lie in the cited orbits and are mapped to the "
         "four tuples above; the SOC action (1,9,22) -> (5,45,40) contains (rho = (0,2,1,3), rho' = (1,0,2,3), inverse)",
         len(k7) == 1 and M.G_ELEMENTS[k7[0]] == ((3, 2, 1, 0), (2, 3, 0, 1), 1) and src_ok
         and ((0, 2, 1, 3), (1, 0, 2, 3), 1) in [M.G_ELEMENTS[k] for k in ksoc])

ck.section("Section 4: the resource energy (10) and its five charts")
RN = PR(["B", "P", "be", "de"])
Bn, Pn, ben, den_ = [RF(RN, g) for g in RN.gens]
eta = ben - den_
BH, BC = Bn * eta / ben, Bn * den_ / ben
ck.check("EXACT", "allocation of one leaf (B, P) with nonminor shares eta = beta - delta (old-H atom, t = B_H = B eta/beta) "
         "and delta (chain, leaves B_C = B delta/beta and P): the required central coefficients "
         "[eta + eta^2/B_H] + [delta^2/B_C + delta^2/P] equal beta - delta + beta^2/B + delta^2/P, and B_H + B_C = B",
         ((eta + eta * eta / BH) + (den_ * den_ / BC + den_ * den_ / Pn)
          - (ben - den_ + ben * ben / Bn + den_ * den_ / Pn)).is_zero() and (BH + BC - Bn).is_zero())
RN2 = PR(["b", "c", "x", "y"])
b2, c2, x2, y2 = [RF(RN2, g) for g in RN2.gens]
ck.check("EXACT", "x^2/b + y^2/c - (x+y)^2/(b+c) = (cx - by)^2/(bc(b+c)) >= 0, so the energy is >= 1/(b+c); and with "
         "delta = 0, beta = (b, c)/(b+c) the objective is 1 + 1/(b+c), so J <= 1 + 1/(b+c)",
         (x2 * x2 / b2 + y2 * y2 / c2 - (x2 + y2) ** 2 / (b2 + c2)
          - (c2 * x2 - b2 * y2) ** 2 / (b2 * c2 * (b2 + c2))).is_zero()
         and ((b2 / (b2 + c2)) ** 2 / b2 + (c2 / (b2 + c2)) ** 2 / c2 - 1 / (b2 + c2)).is_zero())

RH = PR(["T", "L", "U", "W", "r", "Y"])
T_, L_, U_, W_, r_, Y_ = RH.gens
rh = lambda z: RF(RH, z)
th = rh(T_) / rh(T_ + 1)
th1 = 1 - th
half = RF(RH, RH.one, RH.const(2))
CH40 = {}
bb, cc = th / rh(L_), th1 / rh(L_)
pp, qq = 2 * th * rh(U_) / rh(U_ + 1), 2 * th1 * rh(W_) / rh(W_ + 1)
CH40["HH"] = (pp, bb, cc, qq, 1 + rh(L_) - (pp + qq) * half * half, "HH")
qq = 2 * th1 * rh(W_ + 1)
cc = 2 * th1 * rh(W_ + 1) / rh(L_ * 2 * (W_ + 1) + W_)
CH40["HC"] = (pp, bb, cc, qq, rh(L_) + half + th * half / rh(U_ + 1), "HC")
cc = th1 / rh(L_)
qq = 2 * th1 * rh(U_) / rh(U_ + 1)
pp = 2 * th * rh(W_ + 1)
bb = 2 * th * rh(W_ + 1) / rh(L_ * 2 * (W_ + 1) + W_)
CH40["CH"] = (pp, bb, cc, qq, rh(L_) + half + th1 * half / rh(U_ + 1), "CH")
Vv = U_ + W_
pp, qq = 2 * th * rh(U_ + 1), 2 * th1 * rh(Vv + 1)
bb = th / rh(L_)
cc = 2 * th1 * rh((U_ + 1) * (Vv + 1)) / rh(L_ * 2 * (U_ + 1) * (Vv + 1) + W_)
CH40["CC"] = (pp, bb, cc, qq, rh(L_) + half / rh(U_ + 1), "CC")
qq, pp = 2 * th1 * rh(U_ + 1), 2 * th * rh(Vv + 1)
cc = th1 / rh(L_)
bb = 2 * th * rh((U_ + 1) * (Vv + 1)) / rh(L_ * 2 * (U_ + 1) * (Vv + 1) + W_)
CH40["CC'"] = (pp, bb, cc, qq, rh(L_) + half / rh(U_ + 1), "CC")
POS40 = ["T", "L"]
for cname, (pp, bb, cc, qq, Jt, types) in CH40.items():
    beta = (th, th1)
    leaves = ((bb, pp, types[0]), (cc, qq, types[1]))
    okc = all(nn_rf(RH, x, POS40) for x in (pp, bb, cc, qq)) and pos_rf(RH, bb, POS40) and pos_rf(RH, cc, POS40)
    okc &= pos_rf(RH, Jt, POS40)
    energy = RF(RH, RH.zero)
    lam = []
    for (Bi, Pi, ty), be in zip(leaves, beta):
        if ty == "H":
            de = Pi * half                                  # delta = P/2 < beta
            okc &= pos_rf(RH, 2 * be - Pi, POS40)
            lam.append(2 * be / Bi)
            energy = energy + (be - de) + be * be / Bi + Pi * half * half      # delta^2/P = P/4
        else:
            de = be                                         # delta = beta, P >= 2 beta
            okc &= nn_rf(RH, Pi - 2 * be, POS40) and pos_rf(RH, Pi, POS40)
            lam.append(2 * be / Bi + 2 * be / Pi - 1)
            energy = energy + (be - de) + be * be / Bi + de * de / Pi
    okc &= (energy == Jt) and (lam[0] == lam[1])
    ck.check("EXACT", "chart %s (T, L > 0; U, W >= 0; theta = T/(1+T) = beta_1): p, b, c, q have nonnegative coefficients "
             "and b, c > 0; the point beta = (theta, 1-theta), delta_i = %s is feasible; the objective of (10) there "
             "equals the stated energy J; the stationarity quantities of the two leaves agree"
             % (cname, "p_i/2 < beta_i on an H leaf, beta_i <= p_i/2 on a C leaf"), okc)
ck.info("consequence (PROOF-READ): at a chart point the stated J is the value of the objective at a feasible allocation, "
        "so a f5 + ... is safe for a >= J by the four cited atoms; the reduced objective is strictly convex in beta_1 "
        "with derivative < 0 at 0 and > 0 at 1, so every positive (p, b, c, q) has a unique stationary beta and lies "
        "in one of the five charts (leaf i is H iff p_i < 2 beta_i)")

ck.section("Section 4: projective compensation (tables for L0, B1, B2)")
i_s, i_a = n40.index("s"), n40.index("a")
Ls = PP.coeffs_in(G40, i_s)
L0, L1, L2 = Ls[0], Ls[1], Ls[2]


def shift_a(dd, k):
    return {ex[:i_a] + (ex[i_a] + k,) + ex[i_a + 1:]: cf for ex, cf in dd.items()}


def lincomb(*parts):
    out = {}
    for dd, m in parts:
        for ex, cf in dd.items():
            w = out.get(ex, 0) + m * cf
            if w:
                out[ex] = w
            else:
                out.pop(ex, None)
    return out


B1 = lincomb((L1, 1), (shift_a(L0, 1), 2))
B2 = lincomb((L2, 1), (shift_a(L1, 1), 1), (shift_a(L0, 2), 1))
POLY = {"L0": L0, "B1": B1, "B2": B2}
dega = {nm: PP.term_dict_degrees(dd, 7)[i_a] for nm, dd in POLY.items()}
ck.check("EXACT", "G40 (at d = 1) is quadratic in s: G40 = L0 + s L1 + s^2 L2; L0, B1 = L1 + 2a L0, B2 = L2 + a L1 + a^2 L0 "
         "have degrees 2, 3, 4 in a", sorted(Ls) == [0, 1, 2] and dega == {"L0": 2, "B1": 3, "B2": 4})
RI = PR(["L0", "L1", "L2", "a", "S"])
l0, l1, l2, av, Sv = [RF(RI, g) for g in RI.gens]
sv = Sv / (1 + av * Sv)
ck.check("EXACT", "(13): (1 + aS)^2 [L0 + s L1 + s^2 L2] at s = S/(1 + aS) equals L0 + (L1 + 2a L0) S + (L2 + a L1 + a^2 L0) S^2; "
         "and 1 + aS = 1/(1 - as) for S = s/(1 - as)",
         ((1 + av * Sv) ** 2 * (l0 + sv * l1 + sv * sv * l2)
          - (l0 + (l1 + 2 * av * l0) * Sv + (l2 + av * l1 + av * av * l0) * Sv * Sv)).is_zero()
         and ((1 + av * (Sv / (1 - av * Sv))) - 1 / (1 - av * Sv)).is_zero())
L00 = PP.set_var(L0, i_a, 0)
const_ex = tuple([0] * 7)
ck.check("EXACT", "L0 at a = 0 (a polynomial in p, b, c, q, r with %d terms) has only positive coefficients and constant "
         "term 1: L0(0, r) >= 1" % len(L00), all(v > 0 for v in L00.values()) and L00.get(const_ex) == 1)
DOC40 = {
    "HH": ([6372, 7776, 7757], [9408, 11466, 13524, 13502], [8832, 15960, 18979, 21801, 21863]),
    "HC": ([6084, 7380, 7379], [7728, 9366, 11004, 11003], [5456, 11550, 13581, 15591, 15599]),
    "CH": ([6084, 7380, 7379], [7728, 9366, 11004, 11003], [5456, 11550, 13581, 15591, 15599]),
    "CC": ([6645, 8394, 8259], [7029, 8958, 10887, 10887], [4252, 9048, 11157, 13266, 13266]),
    "CC'": ([6645, 8394, 8259], [7029, 8958, 10887, 10887], [4252, 9048, 11157, 13266, 13266]),
}
ntab = 0
for cname, (pp, bb, cc, qq, Jt, types) in CH40.items():
    nJ, dJ = Jt.num, Jt.den
    ch = {"p": pp, "b": bb, "c": cc, "q": qq, "r": rh(r_), "a": rh(RH.one), "s": rh(RH.one)}
    for k, nm in enumerate(("L0", "B1", "B2")):
        dd = POLY[nm]
        m = dega[nm]
        bi = PP.coeffs_in(dd, i_a)
        dg = PP.term_dict_degrees(dd, 7)
        dg[i_a] = 0
        dg[i_s] = 0
        bt = [PP.cleared(RH, bi.get(i, {}), n40, ch, degs=dg)[0] for i in range(m + 1)]
        counts = []
        allok = True
        c0 = None
        for l in range(m + 1):
            cf = RH.zero
            for i in range(l + 1):
                cf = cf + bt[i] * (nJ ** i) * (dJ ** (m - i)) * comb(m - i, l - i)
            if l == 0:
                c0 = cf
            stc = PP.poly_stats(cf)
            counts.append(stc[0])
            allok &= stc[1] == 0
            ntab += 1
        ck.check("EXACT", "chart %s, %s (degree %d in a), a = J Y/(1+Y): the coefficients of Y^0..Y^%d of "
                 "W d0^%d (1+Y)^%d %s have %s terms (document: %s), all coefficients >= 0; the coefficient of Y^0 is "
                 "W d0^%d times the value at a = 0" % (cname, nm, m, m, m, m, nm, counts, DOC40[cname][k], m),
                 allok and counts == DOC40[cname][k] and c0 == bt[0] * dJ ** m)
ck.check("EXACT", "number of coefficient polynomials checked in the five charts: 5 x (3 + 4 + 5) = 60", ntab == 60)
ck.info("consequence: for positive p, b, c, q, r >= 0 and 0 <= a < J: L0 >= L0(0, r)(1 - a/J)^2 >= (1 - a/J)^2 > 0, B1 >= 0, "
        "B2 >= 0; with (13), G40 >= (1 - as)^2 L0 >= (1 - as)^2 (1 - a/J)^2 > 0 whenever also as < 1: inequality (14)")

# ------------------------------------------------------------------------------------------
ck.section("Exact evaluations on the closed hard sets, boundaries included (falsification attempts)")


def rq(rng, zero=0.0, lo=1, hi=30):
    if rng.random() < zero:
        return Fr(0)
    return Fr(rng.randint(lo, hi), rng.randint(lo, hi)) * rng.choice([1, 1, 1, 5, Fr(1, 5)])


def sq_frac(rng):
    """a random positive rational square (so that sqrt(ab) is rational below)."""
    return Fr(rng.randint(1, 12), rng.randint(1, 12))


def kval35(vals):
    pt = [vals[nm] for nm in n35]
    return [PP.eval_terms(k, pt) for k in (k0, k1, k2)]


rng = random.Random(3702)
ok8 = okrow = okid = True
rows_seen = {}
for _ in range(2200):
    # choose (a, b) with rational w = sqrt(ab) when the parallel bound may be needed
    wq = Fr(rng.randint(0, 30), rng.randint(1, 12))
    bbv = rng.choice([Fr(1), sq_frac(rng), sq_frac(rng) * 4, Fr(1, 3)])
    aav = wq * wq / bbv
    if rng.random() < 0.5:
        aav, bbv = bbv, aav
    if rng.random() < 0.1:
        aav = Fr(0)
    if rng.random() < 0.05:
        bbv = Fr(0)
    ab = aav * bbv
    wv_ = wq if ab != 0 else Fr(0)                    # sqrt(ab), rational by construction
    assert wv_ * wv_ == ab
    if ab > 1:
        Em_ = (aav + 1) * (bbv + 1) / (ab - 1)
        eev = Em_ * Fr(rng.randint(0, 99), 100)
    else:
        eev = rq(rng, 0.15)
    jjv = rng.choice([Fr(0), eev, eev / (1 + eev), rq(rng), eev * Fr(rng.randint(0, 100), 100)])
    qqv, ppv = rq(rng, 0.3), rq(rng, 0.3)
    Kv = (ab - 1) * eev - (aav + 1) * (bbv + 1)
    assert Kv < 0
    vals = dict(q=qqv, p=ppv, j=jjv, c=Fr(0), e=eev, b=bbv, a=aav)
    kv = kval35(vals)
    Cb = (aav + 1) * (bbv + 1) + aav * qqv + bbv * ppv + (1 + eev) * (ppv * qqv + ppv + qqv) + eev
    ok8 &= kv[0] >= Cb * (aav + 1) * (bbv + 1) * (-Kv) > 0
    # the row and its threshold (after the exchange that makes b >= a when ab < 1 and j < e)
    if eev == 0 and jjv == 0:
        cv = rq(rng, 0.2)
        g = kv[0] + kv[1] * cv + kv[2] * cv * cv
        okid &= g == ((1 + aav + bbv) * cv + (aav + 1) * (bbv + 1)) * (ppv + aav + 1) ** 2 * (qqv + bbv + 1) ** 2 > 0
        rows_seen["e=j=0"] = rows_seen.get("e=j=0", 0) + 1
        continue
    if jjv >= eev:
        row, Jv = "I", 1 / jjv
    elif ab > 1:
        row, Jv = "II", 1 / eev - (ab - 1) / ((aav + 1) * (bbv + 1))
    elif ab == 1:
        row, Jv = "ab=1", 1 / eev
    else:
        blo, ahi = max(aav, bbv), min(aav, bbv)       # exchange: b >= a
        if blo < 1:
            row, Jv = ("IIIa-", 1 + 1 / eev) if jjv < eev / (1 + eev) else ("IIIa+", 1 / jjv)
        elif jjv >= eev * blo / (eev + blo):
            row, Jv = "IIIb+", 1 / jjv
        else:
            w_ = wv_
            row, Jv = "IIIb-", 1 / eev + blo * (1 - w_) ** 2 / (blo - w_) ** 2
    rows_seen[row] = rows_seen.get(row, 0) + 1
    okrow &= Jv > 0 and 2 * kv[0] + Jv * kv[1] >= 0 and kv[0] + Jv * kv[1] + Jv * Jv * kv[2] >= 0
    tt = Fr(rng.randint(0, 99), 100)
    okrow &= kv[0] + kv[1] * (tt * Jv) + kv[2] * (tt * Jv) ** 2 >= kv[0] * (1 - tt) ** 2 > 0
ck.check("EXACT", "(8) at 2200 random rational points with K < 0 (seed 3702; a = 0, b = 0, ab = 1, e = 0, j = 0, j = e, "
         "the j-cuts and a = b occur with positive probability): k0 >= Cbar (a+1)(b+1)(-K) > 0 every time", ok8)
ck.check("EXACT", "the same points, each assigned to its row with the row's threshold J (rows I, II, IIIa-, IIIa+, IIIb+, "
         "IIIb-; ab = 1 with j < e uses J = 1/e; for ab < 1, j < e the exchange a <-> b, q <-> p makes b >= a): "
         "2k0 + J k1 >= 0, k0 + J k1 + J^2 k2 >= 0 and G35(tJ) >= k0 (1-t)^2 > 0 at a random t in [0,1)",
         okrow and min(rows_seen.get(k, 0) for k in ("I", "II", "IIIa-", "IIIa+", "IIIb+", "IIIb-", "ab=1")) > 5,
         str(dict(sorted(rows_seen.items()))))
ck.check("EXACT", "the points with e = j = 0 among them: G35 equals the closed-layer product and is > 0 at a random c",
         okid and rows_seen.get("e=j=0", 0) > 5)


def energy_exact(pv, bv, cv, qv, th_):
    """value of the objective of (10) at beta = (th_, 1 - th_) with the optimal delta (exact; requires
    b, c > 0 when the corresponding beta is > 0; conventions 0/0 = 0)."""
    tot = Fr(0)
    for (B_, P_), be in (((bv, pv), th_), ((cv, qv), 1 - th_)):
        if be == 0:
            continue
        de = min(P_ / 2, be)
        tot += be - de + be * be / B_ + (de * de / P_ if de else 0)
    return tot


rng = random.Random(3703)
ok14 = okJ = True
npos = 0
for _ in range(900):
    cname = rng.choice(list(CH40))
    pp, bb, cc, qq, Jt, types = CH40[cname]
    pt = [Fr(rng.randint(1, 30), rng.randint(1, 30)),                        # T
          Fr(rng.randint(1, 30), rng.randint(1, 30)) * rng.choice([1, 1, 4, Fr(1, 4)]),   # L
          rq(rng, 0.2), rq(rng, 0.2), rq(rng, 0.3), Fr(0)]                   # U, W, r, Y (Y unused)
    pv, bv, cv, qv, Jv = [x.value(pt) for x in (pp, bb, cc, qq, Jt)]
    thv = pt[0] / (pt[0] + 1)
    okJ &= energy_exact(pv, bv, cv, qv, thv) == Jv
    for dth in (Fr(1, 50), Fr(-1, 50), Fr(1, 7), Fr(-1, 7)):                 # J is the minimum over beta
        t2 = thv + dth
        if 0 <= t2 <= 1:
            okJ &= energy_exact(pv, bv, cv, qv, t2) >= Jv
    okJ &= Jv <= 1 + 1 / (bv + cv)
    av_ = Jv * Fr(rng.randint(0, 99), 100)
    sv_ = (Fr(rng.randint(0, 99), 100) / av_) if (av_ > 0 and rng.random() < 0.8) else rq(rng, 0.3)
    if av_ * sv_ >= 1:
        continue
    vals = dict(p=pv, r=pt[4], a=av_, b=bv, c=cv, q=qv, s=sv_)
    g = PP.eval_terms(G40, [vals[nm] for nm in n40])
    ok14 &= g >= (1 - av_ * sv_) ** 2 * (1 - av_ / Jv) ** 2 > 0
    npos += 1
ck.check("EXACT", "energy at 900 random rational points of the five charts (seed 3703): the stated J equals the objective "
         "of (10) with the optimal delta at beta = theta, is <= the objective at four other values of beta_1, and "
         "is <= 1 + 1/(b+c)", okJ)
ck.check("EXACT", "(14) at %d of these points with random a in [0, J), s with as < 1 and r >= 0 (p = 0 or q = 0 when "
         "U = 0 or W = 0): G40 >= (1 - as)^2 (1 - a/J)^2 > 0 every time" % npos, ok14 and npos > 600)

rng = random.Random(3704)
okz = True
nz = {"b=0": 0, "c=0": 0, "b=c=0": 0}
for _ in range(600):
    mode = rng.choice(["b=0", "c=0", "b=c=0"])
    pv, qv, rv = rq(rng, 0.3), rq(rng, 0.3), rq(rng, 0.3)
    if mode == "b=c=0":
        bv = cv = Fr(0)
        Jv = None
    else:
        Bx = rq(rng)                     # the positive core coefficient; its partner leaf coefficient Px
        Px = qv if mode == "b=0" else pv
        de = min(Px / 2, Fr(1))
        Jv = 1 - de + 1 / Bx + (de * de / Px if de else 0)        # beta = 1 on the surviving leaf
        bv, cv = (Fr(0), Bx) if mode == "b=0" else (Bx, Fr(0))
    av_ = (Jv * Fr(rng.randint(0, 99), 100)) if Jv is not None else rq(rng, 0.2)
    sv_ = (Fr(rng.randint(0, 99), 100) / av_) if (av_ > 0 and rng.random() < 0.8) else rq(rng, 0.3)
    if av_ * sv_ >= 1:
        continue
    vals = dict(p=pv, r=rv, a=av_, b=bv, c=cv, q=qv, s=sv_)
    g = PP.eval_terms(G40, [vals[nm] for nm in n40])
    bound = (1 - av_ * sv_) ** 2 * ((1 - av_ / Jv) ** 2 if Jv is not None else 1)
    okz &= g >= bound > 0
    nz[mode] += 1
ck.check("EXACT", "zero core coefficients (seed 3704): with b = 0 < c the extended minimum is J = 1 - delta + 1/c + delta^2/q, "
         "delta = min(q/2, 1) (beta = (0,1); p arbitrary), symmetrically for c = 0 < b, and J = infinity for b = c = 0; "
         "G40 >= (1 - as)^2 (1 - a/J)^2 > 0 (with 1/J = 0 when J is infinite) at random a < J, as < 1, r, p, q >= 0",
         okz and min(nz.values()) > 100, str(nz))

M.crosscheck_with_audit_library(ck)

# ------------------------------------------------------------------------------------------
if NUMERIC:
    import numpy as np
    ck.section("NUMERIC falsification attempts (floating point; evidence, not proof)")
    M.numeric_pointwise_attack(ck, ROLES35 + ROLES40, seed=37001)
    for jx in (35, 40):
        ctrl = M.numeric_min_phi(M.F_VALUES[jx], 37100 + jx, starts=6)
        ck.check("NUMERIC", "control: for the bare non-minor vector f%d the same minimiser finds a negative value of "
                 "the Leibniz sum on SO(4) (so it can detect a violation)" % jx, ctrl < -0.05, "minimum found %.5f" % ctrl)

    rngn = np.random.default_rng(37002)

    def vec35(co):
        return M.combo_values({GEN35[k]: float(v) for k, v in co.items()})

    def vec40(co):
        return M.combo_values({GEN40[k]: float(v) for k, v in co.items()})

    worst = [np.inf, np.inf, np.inf]
    for _ in range(40):
        # coherent (5): ab > 1
        aa = rngn.exponential(1.5) + 0.05
        bb = (1 + rngn.exponential(1.0)) / aa
        Em_ = (aa + 1) * (bb + 1) / (aa * bb - 1)
        ee = Em_ * rngn.uniform(0.02, 1.0)
        co = dict(a=aa, b=bb, e=ee, d=1.0, c=1 / ee - 1 / Em_, j=0.0, q=0.0, p=0.0)
        worst[0] = min(worst[0], M.numeric_min_phi(vec35(co), int(rngn.integers(1 << 30)), starts=4))
        # parallel (6): ab < 1, a != b (both orientations)
        w_ = rngn.uniform(0.0, 0.98)
        bb = rngn.choice([rngn.uniform(1.0, 5.0), rngn.uniform(w_ + 0.02, 1.0)])
        aa = w_ * w_ / bb
        ee = rngn.exponential(1.0) + 0.05
        thr_ = 1 / ee + bb * (1 - w_) ** 2 / (bb - w_) ** 2
        if rngn.random() < 0.5:
            aa, bb = bb, aa
        co = dict(a=aa, b=bb, e=ee, d=1.0, c=thr_, j=0.0, q=0.0, p=0.0)
        worst[1] = min(worst[1], M.numeric_min_phi(vec35(co), int(rngn.integers(1 << 30)), starts=4))
        # double-SOC at ab = 1: c = 1/e
        aa = rngn.exponential(1.0) + 0.05
        ee = rngn.exponential(1.0) + 0.05
        co = dict(a=aa, b=1 / aa, e=ee, d=1.0, c=1 / ee, j=0.0, q=0.0, p=0.0)
        worst[2] = min(worst[2], M.numeric_min_phi(vec35(co), int(rngn.integers(1 << 30)), starts=4))
    ck.check("NUMERIC", "the cited coherent bound (5), transferred to (2), as an inequality on SO(4): 40 random coefficient "
             "vectors (a, b, e) with ab > 1 and c exactly on the threshold 1/e - (ab-1)/((a+1)(b+1)), all other "
             "coefficients 0 (seed 37002): minimum of the Leibniz sum (4 BFGS starts) never below -1e-8",
             worst[0] > -1e-8, "smallest minimum found %.3e" % worst[0])
    ck.check("NUMERIC", "the cited parallel bound (6), transferred to (2), both orientations a < b and a > b, with c on the "
             "threshold 1/e + b(1-w)^2/(b-w)^2: 40 random coefficient vectors, minimum never below -1e-8",
             worst[1] > -1e-8, "smallest minimum found %.3e" % worst[1])
    ck.check("NUMERIC", "the closed bound c = 1/e at ab = 1: 40 random coefficient vectors, minimum never below -1e-8",
             worst[2] > -1e-8, "smallest minimum found %.3e" % worst[2])

    def energy_num(pv, bv, cv, qv):
        """numerical minimum of (10) over beta_1 (the objective is convex in beta_1)."""
        from scipy.optimize import minimize_scalar

        def f(t):
            tot = 0.0
            for (B_, P_), be in (((bv, pv), t), ((cv, qv), 1 - t)):
                if be <= 0:
                    continue
                if B_ <= 0:
                    return np.inf
                de = min(P_ / 2, be)
                tot += be - de + be * be / B_ + (de * de / P_ if de > 0 else 0.0)
            return tot
        if bv <= 0 and cv <= 0:
            return np.inf
        if bv <= 0:
            return f(0.0)
        if cv <= 0:
            return f(1.0)
        res = minimize_scalar(f, bounds=(0.0, 1.0), method="bounded", options=dict(xatol=1e-13))
        return min(res.fun, f(0.0), f(1.0))

    worst = np.inf
    for _ in range(50):
        pv, bv, cv, qv = [rngn.exponential(1.0) * rngn.choice([0.2, 1.0, 3.0]) + 0.02 for _ in range(4)]
        Jv = energy_num(pv, bv, cv, qv)
        co = dict(p=pv, b=bv, c=cv, q=qv, a=Jv, d=1.0, r=0.0, s=0.0)
        worst = min(worst, M.numeric_min_phi(vec40(co), int(rngn.integers(1 << 30)), starts=4))
    ck.check("NUMERIC", "the energy threshold of (10) as an inequality on SO(4): 50 random positive (p, b, c, q) with a equal "
             "to the numerically minimised energy J and r = s = 0 (seed 37002): minimum of the Leibniz sum never "
             "below -1e-8", worst > -1e-8, "smallest minimum found %.3e" % worst)

    def make_num(G, names):
        ex = np.array(list(G.keys()), float)
        cf = np.array([float(v) for v in G.values()])

        def val(v):
            mon = np.prod(np.array([v[nm] for nm in names], float)[None, :] ** ex, axis=1)
            return float(cf @ mon), float(np.abs(cf) @ np.abs(mon))
        return val

    G35n, G40n = make_num(G35, n35), make_num(G40, n40)
    rngn = np.random.default_rng(37003)
    minratio = np.inf
    for _ in range(5000):
        pv, bv, cv, qv = [rngn.exponential(1.0) * rngn.choice([0.0, 0.05, 1.0, 1.0, 4.0]) for _ in range(4)]
        Jv = energy_num(pv, bv, cv, qv)
        av_ = (Jv if np.isfinite(Jv) else 5.0) * rngn.uniform(0, 1) ** 0.5 * 0.999999
        sv_ = rngn.uniform(0, 1) ** 0.5 * 0.999999 / av_ if av_ > 0 else rngn.exponential(1.0)
        rv = rngn.exponential(1.0) * rngn.choice([0.0, 1.0])
        g, _ = G40n(dict(p=pv, r=rv, a=av_, b=bv, c=cv, q=qv, s=sv_))
        bound = (1 - av_ * sv_) ** 2 * ((1 - av_ / Jv) ** 2 if np.isfinite(Jv) else 1.0)
        minratio = min(minratio, g / bound)
    ck.check("NUMERIC", "(14) in original coordinates with J minimised numerically (seed 37003): at 5000 random points with "
             "p, b, c, q >= 0 (zeros and near-zeros included), a = tJ, as < 1, r >= 0 the ratio "
             "G40 / [(1 - as)^2 (1 - a/J)^2] is never below 1 - 1e-6", minratio > 1 - 1e-6,
             "smallest ratio %.6f" % minratio)

    # realizable points
    free = []
    for seed in (37011, 37012):
        free += M.numeric_stratum_search(T, seed, 300)
    dd35 = [pt["gaps"][3] - pt["gaps"][1] for pt in free]
    dd40 = [pt["gaps"][14] - pt["gaps"][16] for pt in free]
    ck.info("free stratum search (gaps = 0 on T, gaps >= 0, sum of gaps = 1; 600 starts, seeds 37011, 37012): %d converged "
            "points; delta35 in [%.3g, %.3g], delta40 in [%.3g, %.3g]"
            % (len(free), min(dd35), max(dd35), min(dd40), max(dd40)))
    tol = 1e-6
    for label, order, names, dix, Gn, seed in (("delta35 > 0", ORDER35, NAMES35, 5, G35n, 37020),
                                               ("delta40 > 0", ORDER40, NAMES40, 6, G40n, 37030)):
        hard = M.numeric_hard_points(T, order, dix, seed)
        hard = [(co, pt) for co, pt in hard if min(co) > -1e-8]
        ds = [co[dix] for co, pt in hard]
        n_exact = sum(1 for co, pt in hard if pt["contacts"] == tuple(sorted(T)))
        ck.info("branch %s: search with the hinge equation delta >= target (targets 0.004, 0.01; then 300 perturbed "
                "restarts, seed %d): %d stratum points with delta > 5e-4, delta up to %.4g; contact set exactly T at %d "
                "of them" % (label, seed, len(hard), max(ds) if ds else 0.0, n_exact))
        if not hard:
            ck.info("branch %s: no stratum point found: the radial equation is NOT tested on realizable points" % label)
            continue
        maxG = 0.0
        contrast = []
        unsafe = 0
        which = {}
        for co, pt in hard:
            cod = dict(zip(names, co))
            v = {k: max(cod[k], 0.0) / cod["d"] for k in names if k != "d"}
            g, gabs = Gn(v)
            maxG = max(maxG, abs(g) / gabs)
            g2, gabs2 = Gn({k: x * rngn.uniform(0.9, 1.1) for k, x in v.items()})
            contrast.append(abs(g2) / gabs2)
            if label.startswith("delta35"):
                aa, bb, ee, jj, cc = v["a"], v["b"], v["e"], v["j"], v["c"]
                ab = aa * bb
                if (ab - 1) * ee - (aa + 1) * (bb + 1) >= -tol:
                    key = "four (K >= 0)"
                else:
                    if ee < 1e-9 and jj < 1e-9:
                        key, Jv = "NONE (e = j = 0)", np.inf
                    elif jj >= ee:
                        key, Jv = "row I: SOC", 1 / jj
                    elif ab > 1:
                        key, Jv = "row II: coherent", 1 / ee - (ab - 1) / ((aa + 1) * (bb + 1))
                    else:
                        blo = max(aa, bb)
                        w_ = np.sqrt(ab)
                        if blo < 1:
                            key, Jv = ("row IIIa-: old-H", 1 + 1 / ee) if jj < ee / (1 + ee) else ("row IIIa+: SOC", 1 / jj)
                        elif jj >= ee * blo / (ee + blo):
                            key, Jv = "row IIIb+: SOC", 1 / jj
                        else:
                            key, Jv = "row IIIb-: parallel", 1 / ee + blo * (1 - w_) ** 2 / (blo - w_) ** 2
                    if not cc >= Jv * (1 - tol):
                        key = "NONE (" + key + ")"
                        unsafe += 1
            else:
                if v["a"] * v["s"] >= 1 - tol:
                    key = "SOC as >= 1"
                elif v["a"] >= energy_num(v["p"], v["b"], v["c"], v["q"]) * (1 - tol):
                    key = "energy a >= J"
                else:
                    key = "NONE"
                    unsafe += 1
            which[key] = which.get(key, 0) + 1
        ck.check("NUMERIC", "branch %s: radial equation G = 0 at the %d stratum points: |G| <= 1e-5 times the sum of the "
                 "absolute values of its terms (the points are known to about 1e-10/delta in these coordinates)"
                 % (label, len(hard)), maxG < 1e-5,
                 "largest relative |G| = %.3e; for contrast, with every coefficient multiplied by an "
                 "independent random factor in [0.9, 1.1] the relative |G| has median %.2e" % (maxG, float(np.median(contrast))))
        ck.check("NUMERIC", "branch %s: conclusion of the document at these points (%s)"
                 % (label, "K >= 0, or c >= the threshold J of the row of (a, b, e, j)" if label.startswith("delta35")
                    else "as >= 1 or a >= J"), unsafe == 0, str(which))

ck.finish()
