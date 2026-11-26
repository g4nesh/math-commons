"""Independent exact replay of Case T_36 = {2,3,4,6,11,12} (group V) of D. Zheng's proof of the
Marcus-de Oliveira conjecture for n = 4.

Proof document replayed: MOC4_SUPP/proofs/manin_transfer_20261021/final_referee_T35_uniform_20261022.pdf
(4 pages; the document calls the set "T35").  Every object is rebuilt from the paper's
definitions by own code (replay_T35_common.py, replay_T35_common_poly.py).  Nothing from the
supplement is imported or executed.  If MOC4_SUPP is set, the PDF is hashed and four JSON data
files are read for comparison only.

Run:   MOC4_SUPP=/path/to/supplement python replay_T36.py [--no-numeric]
Exit status 0 iff every check passes.  Checks tagged EXACT use integer/rational/polynomial
arithmetic only; checks tagged NUMERIC are floating-point falsification attempts.
"""
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
T = (2, 3, 4, 6, 11, 12)
DOC = "final_referee_T35_uniform_20261022.pdf"
NAMES = list("qpjcedba")
# (68): g = q f0 + p f1 + j f7 + c f9 + e f53 + d f68 + b f67 + a f10
ORDER68 = [0, 1, 7, 9, 53, 68, 67, 10]
# (69): g = q f1 + p f2 + j f8 + c f9 + e f53 + d f69 + b f67 + a f10
ORDER69 = [1, 2, 8, 9, 53, 69, 67, 10]
GEN68 = dict(zip(NAMES, ORDER68))
GEN69 = dict(zip(NAMES, ORDER69))

ck = M.Checker("Case T_36 = {2,3,4,6,11,12}: replay of final_referee_T35_uniform_20261022.pdf")
M.print_versions(ck)
M.check_document_identity(ck, DOC)

M.check_foundations(ck)
M.check_radial_necessity(ck)
M.check_hessian_necessity(ck)
M.crosscheck_with_supplement_data(ck)

M.check_contact_set(ck, 36, T, "T35")
gens = M.check_generators(ck, T, [0, 1, 2, 7, 8, 9, 10, 53, 67, 68, 69], [68, 69])

# ------------------------------------------------------------------------------------------
ck.section("Cancellation and the two normal forms (document p. 1)")
IDENT = [([68, 69], [1, 9, 10, 67, 67]),
         ([68, 2], [0, 9, 67]), ([68, 8], [7, 10, 67]),
         ([69, 0], [1, 2, 10, 67]), ([69, 7], [1, 8, 9, 67])]
M.check_linear_identities(ck, IDENT, "the five identities of p. 1")
ck.check("EXACT", "the eleven generators span a space of dimension 8 (= dimension of the contact nullspace); each of "
         "the two eight-element lists (68) and (69) is linearly independent, so its coefficients are uniquely "
         "determined linear functions of the gap vector; each list has exactly one non-minor generator",
         M.rank([M.F_VALUES[j] for j in gens]) == 8 and M.rank([M.F_VALUES[j] for j in ORDER68]) == 8
         and M.rank([M.F_VALUES[j] for j in ORDER69]) == 8 and len(M.nullspace([M.P[s] for s in T], 14)) == 8
         and [j for j in ORDER68 if not M.F_IS_MINOR[j]] == [68]
         and [j for j in ORDER69 if not M.F_IS_MINOR[j]] == [69])
M.check_elimination(ck, gens, [68, 69],
                    {(68, 69): [1, 9, 10, 67, 67]},
                    {68: [(2, [0, 9, 67]), (8, [7, 10, 67])],
                     69: [(0, [1, 2, 10, 67]), (7, [1, 8, 9, 67])]},
                    {68: ORDER68, 69: ORDER69}, seed=3601, samples=6000)
M.check_minor_leibniz(ck, [j for j in gens if M.F_IS_MINOR[j]])

# ------------------------------------------------------------------------------------------
ck.section("Radial necessity in the two normal forms (document p. 1)")
rd68 = PP.radial_data(ORDER68, NAMES)
names7, G68 = PP.check_radial(ck, rd68, "branch (68)", 1864, "final_referee_T35_bad68_full_radial_20261022.json")
rd69 = PP.radial_data(ORDER69, NAMES)
names7b, G69 = PP.check_radial(ck, rd69, "branch (69)", 1367, "final_referee_T35_bad69_full_radial_20261022.json")
assert names7 == names7b == ["q", "p", "j", "c", "e", "b", "a"]
R9 = rd68["R"]
exd = lambda k: tuple(k if nm == "d" else 0 for nm in R9.names)
c68, q68 = PR.terms(rd68["Cstar"]), PR.terms(rd68["Q"])
ck.check("EXACT", "branch (68): C* (33 terms) and Q (10 terms) have only positive coefficients and contain d^3 and "
         "d^2 with coefficient 1: C* >= d^3 > 0 and Q >= d^2 > 0 for d > 0",
         len(c68) == 33 and len(q68) == 10 and all(v > 0 for v in c68.values()) and all(v > 0 for v in q68.values())
         and c68.get(exd(3)) == 1 and q68.get(exd(2)) == 1)
c69 = PR.terms(rd69["Cstar"])
neg = {ex: v for ex, v in c69.items() if v < 0}
aej = tuple(1 if nm in ("a", "e", "j") else 0 for nm in rd69["R"].names)
R7 = PR(names7)
V7 = {nm: RF(R7, R7[nm]) for nm in names7}
q, p, j, c, e, b, a = [V7[nm] for nm in names7]
_, c69d1 = PP.drop_h_set_d(rd69, rd69["Cstar"])
rest = R7.from_terms(c69d1) - (1 + a + e + a * e - a * e * j).num
ck.check("EXACT", "branch (69): the only negative monomial of C* (22 terms) is -aej; at d = 1, "
         "C* - [1 + a + e + ae(1 - j)] has only nonnegative coefficients: inequality (1) of the document",
         len(c69) == 22 and neg == {aej: -1} and PP.nonneg(rest), "%d terms in the difference" % len(rest))
ck.check("EXACT", "(a + e)^2 - 4ae = (a - e)^2, so ae <= S^2/4 for S = a + e: if j > 1 and S(j-1) < 1 then "
         "ae(j-1) <= ae/S <= S/4 and C* >= 1 + S - S/4 > 0; if j <= 1 then C* >= 1 + a + e > 0.  Hence C* > 0 on "
         "the hard set of branch (69), and C* + hQ = 0 forces Q != 0", ((a + e) ** 2 - 4 * a * e - (a - e) ** 2).is_zero())
DEG68 = PP.term_dict_degrees(G68, 7)
DEG69 = PP.term_dict_degrees(G69, 7)
ck.check("EXACT", "degrees of G68 and G69 (at d = 1) in (q,p,j,c,e,b,a) are (3,2,3,2,2,3,2) and (2,2,2,2,2,3,2), "
         "as the document states", DEG68 == [3, 2, 3, 2, 2, 3, 2] and DEG69 == [2, 2, 2, 2, 2, 3, 2])

# ------------------------------------------------------------------------------------------
ck.section("Universal inequalities used (CITED: audit PROOF.md 7.6); roles and coefficient conditions (p. 2)")
ROLES68 = [
    ("four", (GEN68["q"], GEN68["j"], GEN68["p"], 68), "(3a): qj(p - d) >= d^2(q+j+p+d), f0, f7; f1"),
    ("four", (GEN68["j"], GEN68["b"], GEN68["e"], 68), "(3b): jb(e - d) >= d^2(j+b+e+d), f7, f67; f53"),
]
ROLES69 = [
    ("SOC", (GEN69["p"], GEN69["j"], 69), "(4a): pj >= d^2, f2, f8"),
    ("old-H", (GEN69["c"], GEN69["p"], 69), "(4b): c(p - d) >= d^2, t = c on f9, main = p on f2"),
    ("old-H", (GEN69["a"], GEN69["j"], 69), "H-sum atom 1: t = a on f10, main on f8"),
    ("old-H", (GEN69["e"], GEN69["j"], 69), "H-sum atom 2: t = e on f53, main on f8"),
    ("chain", (GEN69["j"], GEN69["e"], GEN69["b"], 69), "(4d): j >= d^2/e + d^2/b, centre f8, leaves f53, f67"),
]
ck.check("EXACT", "the role tuples are the document's: four (0,7;1;68), (7,67;53;68); SOC (2,8;69); old-H (9,2;69), "
         "(10,8;69), (53,8;69); chain centre 8 with leaves 53, 67",
         [r[1] for r in ROLES68 + ROLES69] == [(0, 7, 1, 68), (7, 67, 53, 68), (2, 8, 69), (9, 2, 69), (10, 8, 69),
                                                 (53, 8, 69), (8, 53, 67, 69)])
M.check_roles(ck, ROLES68 + ROLES69)
M.check_coefficient_forms(ck)
S_ = a + e
d1, d2 = a / S_, e / S_
ck.check("EXACT", "old-H sum: with the non-minor coefficient split as d1 = a/(a+e), d2 = e/(a+e), the atoms "
         "a f10 + j1 f8 + d1 f69 and e f53 + j2 f8 + d2 f69 need j1 >= d1 + d1^2/a and j2 >= d2 + d2^2/e, and "
         "d1 + d1^2/a + d2 + d2^2/e = 1 + 1/(a+e): the condition is (a+e)(j-1) >= 1, as in (4)",
         (d1 + d1 * d1 / a + d2 + d2 * d2 / e - (1 + 1 / S_)).is_zero() and (d1 + d2 - 1).is_zero())

# ------------------------------------------------------------------------------------------
ck.section("Branch (68): the strict inequality (5) (document p. 2-3)")
h1 = (q + 1) * (j + 1) - p * (q * j - 1) - ((q + j + p + 1) - q * j * (p - 1))
h2 = (b + 1) * (j + 1) - e * (j * b - 1) - ((j + b + e + 1) - j * b * (e - 1))
ck.check("EXACT", "(q+1)(j+1) - p(qj-1) = (q+j+p+1) - qj(p-1) and the same with (b, e): both conditions (3) fail iff "
         "1 - L_p p > 0 and 1 - L_e e > 0", h1.is_zero() and h2.is_zero())
ip, ie, ia, ic, ib, iq = [names7.index(x) for x in "peacbq"]
g0 = PP.set_var(PP.set_var(PP.set_var(PP.set_var(G68, ip, 0), ie, 0), ia, 0), ic, 0)
ck.check("EXACT", "G68 at p = e = a = c = 0 equals (j+1)^3 (b+1)^3 (q+1)^3",
         R7.from_terms(g0) == ((j + 1) ** 3 * (b + 1) ** 3 * (q + 1) ** 3).num)

RC = PR(["J", "X", "Y", "P", "E", "A", "C"])
J_, X_, Y_, P_, E_, A_, C_ = RC.gens
one = RC.one
rf = lambda z: RF(RC, z)


def chart6(qj_above, jb_above):
    ch = {"j": rf(J_), "a": rf(A_), "c": rf(C_)}
    if not qj_above:
        ch["q"] = rf(X_) / rf(J_ * (X_ + 1))
        ch["p"] = rf(P_)
    else:
        ch["q"] = rf(X_ + 1) / rf(J_)
        ch["p"] = rf((X_ + J_ + 1) * (J_ + 1) * P_) / rf(J_ * X_ * (P_ + 1))
    if not jb_above:
        ch["b"] = rf(Y_) / rf(J_ * (Y_ + 1))
        ch["e"] = rf(E_)
    else:
        ch["b"] = rf(Y_ + 1) / rf(J_)
        ch["e"] = rf((Y_ + J_ + 1) * (J_ + 1) * E_) / rf(J_ * Y_ * (E_ + 1))
    return ch


Hq = J_ * J_ * (X_ + 1) - X_
NU = [(5, 3, 0, 1, 0, 1, 2), (7, 3, 0, 1, 0, 2, 1), (7, 3, 1, 1, 0, 2, 1)]
DOC_Z = {(0, 0): (6190, 5893), (0, 1): (8949, 8802), (1, 0): (8705, 8439), (1, 1): (12375, 12197)}
for qa in (0, 1):
    for ja in (0, 1):
        tag = "chamber %d%d (qj %s 1, jb %s 1)" % (qa, ja, "><"[1 - qa], "><"[1 - ja])
        pos = ["J"] + (["X"] if qa else []) + (["Y"] if ja else [])
        ch = chart6(qa, ja)
        qj = ch["q"] * ch["j"]
        jb = ch["b"] * ch["j"]
        okc = PP.chart_denominators_positive(RC, ch, pos)
        okc &= (qj == (rf(X_ + 1) if qa else rf(X_) / rf(X_ + 1)))
        okc &= (jb == (rf(Y_ + 1) if ja else rf(Y_) / rf(Y_ + 1)))
        if qa:
            okc &= (1 - ch["p"] * (qj - 1) / ((ch["q"] + 1) * (ch["j"] + 1)) == 1 / rf(P_ + 1))
        if ja:
            okc &= (1 - ch["e"] * (jb - 1) / ((ch["b"] + 1) * (ch["j"] + 1)) == 1 / rf(E_ + 1))
        ck.check("EXACT", "%s: chart (6) valid (nonnegative coefficients, denominators > 0 for %s > 0); qj = %s, "
                 "jb = %s; on a product-above-one chart 1 - L p = 1/(1+P), resp. 1 - L e = 1/(1+E)"
                 % (tag, ", ".join(pos), "1+X" if qa else "X/(1+X)", "1+Y" if ja else "Y/(1+Y)"), okc)
        Z, dens = PP.cleared(RC, G68, names7, ch)
        W = PP.prod_dens(RC, dens)
        if not qa:
            if not ja:
                H = (J_ * Y_ + J_ + Y_) * Hq - E_ * J_ * X_ * (Y_ + 1)
                Mm = A_ * C_ * P_ * J_ * J_ * (X_ + 1) * (Y_ + 1) * (A_ * C_ + A_ * P_ + C_ * P_)
            else:
                H = (J_ + Y_ + 1) * Hq
                Mm = A_ * C_ * P_ * J_ ** 4 * Y_ * Y_ * (X_ + 1) * (A_ * C_ + A_ * P_ + C_ * P_)
        else:
            H = J_ * J_ - X_ - 1
            if not ja:
                F0 = J_ * J_ * (J_ + 1) * (J_ + X_ + 1) * (J_ * Y_ + J_ + Y_) ** 2
                Y0 = Y_ + 1
            else:
                F0 = J_ ** 4 * Y_ * Y_ * (J_ + 1) * (J_ + X_ + 1) * (J_ + Y_ + 1) ** 2
                Y0 = one
            Mm = F0 * Y0 * P_ * A_ * C_ * (P_ * A_ * (C_ * J_ * X_ + (J_ + 1) ** 2 * (X_ + 1))
                                           + P_ * C_ * (J_ + X_ + 1) ** 2 + A_ * C_ * J_ * X_
                                           + P_ * (J_ + X_ + 1) * J_ * (X_ + 1))
        Rm = Z - Mm * H * H
        negs = sorted((ex, cf) for ex, cf in PR.terms(Rm).items() if cf < 0)
        R2 = Rm * 2
        extra_ok = True
        if (qa, ja) == (0, 0):
            extra_ok = negs == sorted((nu, -1) for nu in NU)
            for nu in NU:
                mon = one
                for g_, ex in zip(RC.gens, (nu[0] - 1,) + nu[1:]):
                    mon = mon * g_ ** ex
                R2 = R2 - mon * (J_ - 1) ** 2
                extra_ok &= nu[3] >= 1 and nu[5] >= 1 and nu[6] >= 1     # contains P, A and C
        else:
            extra_ok = not negs
        sZ, sR = PP.poly_stats(Z), PP.poly_stats(R2)
        vanish = len(PP.restrict_zero(RC, Mm, ["P"])) == 0 and len(PP.restrict_zero(RC, Mm, ["A"])) == 0 \
            and len(PP.restrict_zero(RC, Mm, ["C"])) == 0
        Z0 = PP.restrict_zero(RC, Z, ["P", "E", "A", "C"])
        R0 = PP.restrict_zero(RC, R2, ["P", "E", "A", "C"])
        W0 = PP.restrict_zero(RC, W, ["P", "E", "A", "C"])
        base = (ch["j"] + 1) ** 3 * (ch["b"] + 1) ** 3 * (ch["q"] + 1) ** 3
        ck.check("EXACT", "%s: Z = W G68 has %d terms (document: %d); the multiplier M has nonnegative coefficients "
                 "and the factor P A C; Z - M H^2%s has %d terms (document: %d), all coefficients >= 0"
                 % (tag, sZ[0], DOC_Z[(qa, ja)][0],
                    " has exactly the three negative monomials of (7), each with coefficient -1, and twice "
                    "[Z - M H^2] - sum x^(nu - e_J) (J-1)^2" if (qa, ja) == (0, 0) else "", sR[0], DOC_Z[(qa, ja)][1]),
                 PP.nonneg(Mm) and vanish and extra_ok and sR[1] == 0
                 and (sZ[0], sR[0]) == DOC_Z[(qa, ja)])
        ck.check("EXACT", "%s: the part of the remainder free of P, E, A, C equals that of Z (every subtracted term "
                 "contains P, A and C) and equals W0 (j+1)^3 (b+1)^3 (q+1)^3, with W/W0 = (1+P)^%d (1+E)^%d.  "
                 "Hence G68 >= (j+1)^3 (b+1)^3 (q+1)^3 (1 - L_p p)^2 (1 - L_e e)^2 on this chamber: inequality (5)"
                 % (tag, 2 * qa, 2 * ja),
                 R0 == Z0 * 2 and rf(Z0) == base * rf(W0)
                 and rf(W) / rf(W0) == rf((P_ + 1) ** (2 * qa) * (E_ + 1) ** (2 * ja)))

# ------------------------------------------------------------------------------------------
ck.section("Branch (69): the strict inequality (10) (document p. 3-4)")
B0 = (e + 1) * (p + 1) * ((j + 1) * (1 - p * j) + e * (j * j + j + 1 - j * p))
g0 = PP.set_var(PP.set_var(PP.set_var(PP.set_var(G69, ia, 0), ic, 0), ib, 0), iq, 0)
ck.check("EXACT", "G69 at a = c = b = q = 0 equals B0 = (e+1)(p+1)[(j+1)(1-pj) + e(j^2+j+1-jp)]; B0 > 0 when "
         "pj < 1 (then 1 - pj > 0 and j^2 + j + 1 - jp > j^2 + j >= 0)", R7.from_terms(g0) == B0.num)
la = (1 - (j - 1) / (1 - e * (j - 1)) * a) - (1 - (a + e) * (j - 1)) / (1 - e * (j - 1))
lb = (1 - (j * e - 1) / e * b) - (e + b - j * e * b) / e
ck.check("EXACT", "1 - L_a a = (1 - (a+e)(j-1))/(1 - e(j-1)) for j > 1 and 1 - L_b b = (e + b - jeb)/e for je > 1; "
         "with 1 - L_c c = 1 - c(p-1) for p > 1 these are > 0 exactly when the old-H-sum, chain and old-H "
         "conditions of (4) fail (e(j-1) < 1 follows from (a+e)(j-1) < 1)", la.is_zero() and lb.is_zero())

RD = PR(["A", "C", "B", "Q", "t", "x", "u"])
A2, C2, B2, Q2, t2, x2, u2 = RD.gens
rg = lambda z: RF(RD, z)
hh = rg(t2) / rg(t2 + 1)
ROWS = []
j2 = (hh + rg(x2)) / rg(x2 + 1)
j3 = rg(x2 + 1)
ROWS.append(("both low, v < 1", rg(t2) / rg(t2 + 1), rg(x2) / rg(x2 + 1), rg(u2) / rg(u2 + 1), "lo", "lo", 0))
ROWS.append(("both low, v > 1", hh / j2, j2, rg(u2 + 1), "lo", "lo", 1))
ROWS.append(("j > 1, v < 1", hh / j3, j3, rg(u2) / rg(u2 + 1), "j", "lo", 0))
ROWS.append(("j > 1, v > 1", hh / j3, j3, 1 + rg(u2) / rg(x2 * (u2 + 1)), "j", "lo", 1))
ROWS.append(("p > 1, v < 1", rg(x2 + 1), hh / rg(x2 + 1), rg(u2) / rg(u2 + 1), "lo", "p", 0))
ROWS.append(("p > 1, v > 1", rg(x2 + 1), hh / rg(x2 + 1), rg(u2 + 1), "lo", "p", 1))
DOC69 = [(3226, 3199), (41018, 41018), (14406, 14406), (26220, 26220), (6760, 6760), (22059, 22059)]
CHARTS69 = []
for row, (name, pv, jv, vv, amode, cmode, vgt) in enumerate(ROWS):
    ev = vv / jv
    av = rg(A2) if amode == "lo" else (1 / rg(x2) - ev) * rg(A2) / rg(A2 + 1)
    cv = rg(C2) if cmode == "lo" else rg(C2) / rg(x2 * (C2 + 1))
    bv = rg(B2) if not vgt else ev * rg(B2) / ((vv - 1) * rg(B2 + 1))
    ch = {"p": pv, "j": jv, "e": ev, "a": av, "c": cv, "q": rg(Q2), "b": bv}
    CHARTS69.append((name, ch, amode, cmode, vgt))
    okc = PP.chart_denominators_positive(RD, ch, ["t", "x", "u"])
    one_m_pj = 1 - pv * jv
    okc &= PP.positive_on_domain(RD, one_m_pj.num, ["t", "x", "u"]) and PP.positive_on_domain(RD, one_m_pj.den, ["t", "x", "u"])
    # region of (p, j)
    for val, mode in ((pv, "p"), (jv, "j")):
        dlt = (val - 1) if (amode == "j" and mode == "j") or (cmode == "p" and mode == "p") else (1 - val)
        okc &= PP.positive_on_domain(RD, dlt.num, ["t", "x", "u"]) and PP.positive_on_domain(RD, dlt.den, ["t", "x", "u"])
    dv = (vv - 1) if vgt else (1 - vv)
    okc &= PP.positive_on_domain(RD, dv.num, ["t", "x", "u"]) and PP.positive_on_domain(RD, dv.den, ["t", "x", "u"])
    if amode == "j":      # 1 - L_a a = 1/(1+A), and e(j-1) < 1
        okc &= (1 - (jv - 1) / (1 - ev * (jv - 1)) * av == 1 / rg(A2 + 1))
        em = 1 - ev * (jv - 1)
        okc &= PP.positive_on_domain(RD, em.num, ["t", "x", "u"]) and PP.positive_on_domain(RD, em.den, ["t", "x", "u"])
    if cmode == "p":      # 1 - L_c c = 1/(1+C)
        okc &= (1 - (pv - 1) * cv == 1 / rg(C2 + 1))
    if vgt:               # 1 - L_b b = 1/(1+B)
        okc &= (1 - (vv - 1) / ev * bv == 1 / rg(B2 + 1))
    ck.check("EXACT", "row %d (%s) of table (11): chart valid for t, x, u > 0 and A, C, B, Q >= 0 (nonnegative "
             "coefficients, positive denominators); p, j, v = je lie on the stated sides of 1 and pj < 1; "
             "the bounded coefficients satisfy 1 - L a = 1/(1+A) (j > 1), 1 - L c = 1/(1+C) (p > 1), "
             "1 - L b = 1/(1+B) (v > 1)" % (row + 1, name), okc)
    Z, dens = PP.cleared(RD, G69, names7, ch)
    W = PP.prod_dens(RD, dens)
    sZ = PP.poly_stats(Z)
    Rm = Z
    sq_ok = True
    if row == 0:
        sq = x2 * x2 * (t2 - x2) ** 2 * (B2 + 1) ** 2 * A2 * C2 * Q2 * (A2 * C2 + A2 * Q2 + C2 * Q2)
        Rm = Z - sq
        sq_ok = sZ[1] > 0 and len(PP.restrict_zero(RD, sq, ["A"])) == 0
    sR = PP.poly_stats(Rm)
    Z0 = PP.restrict_zero(RD, Z, ["A", "C", "B", "Q"])
    R0 = PP.restrict_zero(RD, Rm, ["A", "C", "B", "Q"])
    W0 = PP.restrict_zero(RD, W, ["A", "C", "B", "Q"])
    B0c = (ev + 1) * (pv + 1) * ((jv + 1) * (1 - pv * jv) + ev * (jv * jv + jv + 1 - jv * pv))
    expo = (2 * (amode == "j"), 2 * (cmode == "p"), 3 * vgt)
    ck.check("EXACT", "row %d (%s): Z = W G69 has %d terms and the remainder%s has %d terms (document: %d, %d), all "
             "coefficients of the remainder >= 0; its part free of A, C, B, Q equals that of Z and equals W0 B0, "
             "with W/W0 = (1+A)^%d (1+C)^%d (1+B)^%d.  Hence G69 >= B0 (1 - L_a a)^2 (1 - L_c c)^2 (1 - L_b b)^3 "
             "on this chamber: inequality (10)"
             % (row + 1, name, sZ[0], " Z - x^2 (t-x)^2 (B+1)^2 ACQ(AC+AQ+CQ)" if row == 0 else " (= Z)", sR[0],
                DOC69[row][0], DOC69[row][1], expo[0], expo[1], expo[2]),
             sq_ok and sR[1] == 0 and (sZ[0], sR[0]) == DOC69[row] and R0 == Z0 and rg(Z0) == B0c * rg(W0)
             and rg(W) / rg(W0) == rg((A2 + 1) ** expo[0] * (C2 + 1) ** expo[1] * (B2 + 1) ** expo[2]))

# ------------------------------------------------------------------------------------------
ck.section("Exact evaluations on the closed hard sets, boundaries included (falsification attempts)")


def rq(rng, zero=0.0, lo=1, hi=30):
    if rng.random() < zero:
        return Fr(0)
    return Fr(rng.randint(lo, hi), rng.randint(lo, hi)) * rng.choice([1, 1, 1, 5, Fr(1, 5)])


rng = random.Random(3602)
ok5 = True
worst = None
for _ in range(1500):
    jj = rng.choice([Fr(0), Fr(1), rq(rng), rq(rng), rq(rng)])
    qq = rng.choice([Fr(0), (1 / jj if jj else Fr(1)), rq(rng), rq(rng)])
    bb = rng.choice([Fr(0), (1 / jj if jj else Fr(1)), rq(rng), rq(rng)])
    Lp = max((qq * jj - 1) / ((qq + 1) * (jj + 1)), 0)
    Le = max((jj * bb - 1) / ((bb + 1) * (jj + 1)), 0)
    pp = (Fr(rng.randint(0, 99), 100) / Lp) if Lp > 0 else rq(rng, 0.2)
    ee = (Fr(rng.randint(0, 99), 100) / Le) if Le > 0 else rq(rng, 0.2)
    aa, cc = rq(rng, 0.25), rq(rng, 0.25)
    vals = dict(q=qq, p=pp, j=jj, c=cc, e=ee, b=bb, a=aa)
    g = PP.eval_terms(G68, [vals[nm] for nm in names7])
    bound = (jj + 1) ** 3 * (bb + 1) ** 3 * (qq + 1) ** 3 * (1 - Lp * pp) ** 2 * (1 - Le * ee) ** 2
    ok5 &= g >= bound > 0
    r_ = g / bound
    worst = r_ if worst is None or r_ < worst else worst
ck.check("EXACT", "(5) at 1500 random rational points of the closed hard set of branch (68) (seed 3602; j = 0, j = 1, "
         "q = 0, b = 0, qj = 1, jb = 1, p = 0, e = 0, a = 0, c = 0 occur with positive probability): "
         "G68 >= (j+1)^3 (b+1)^3 (q+1)^3 (1 - L_p p)^2 (1 - L_e e)^2 > 0 every time", ok5,
         "smallest ratio G68/bound = %.6f" % float(worst))

rng = random.Random(3603)
ok10 = True
worst = None
cntb = n10 = 0
for _ in range(2500):
    # (p, j) with pj < 1
    mode = rng.randint(0, 4)
    if mode == 0:
        pp, jj = Fr(rng.randint(0, 100), 100), Fr(rng.randint(0, 100), 100)
        if pp * jj >= 1:
            continue
    elif mode in (1, 2):
        jj = 1 + rq(rng)
        pp = Fr(rng.randint(0, 99), 100) / jj
    else:
        pp = 1 + rq(rng)
        jj = Fr(rng.randint(0, 99), 100) / pp
    # e with e (j-1) < 1 when j > 1
    if jj > 1:
        ee = Fr(rng.randint(0, 99), 100) / (jj - 1)
        if rng.random() < 0.15:
            ee = 1 / jj                     # the wall je = 1 (note 1/j < 1/(j-1))
        aa = (1 / (jj - 1) - ee) * Fr(rng.randint(0, 99), 100)
    else:
        ee = rng.choice([Fr(0), (1 / jj if jj else Fr(1)), rq(rng), rq(rng)])
        aa = rq(rng, 0.25)
    cc = Fr(rng.randint(0, 99), 100) / (pp - 1) if pp > 1 else rq(rng, 0.25)
    if jj * ee > 1:
        bb = ee / (jj * ee - 1) * Fr(rng.randint(0, 99), 100)
    else:
        bb = rq(rng, 0.25)
    qq = rq(rng, 0.25)
    La = (jj - 1) / (1 - ee * (jj - 1)) if jj > 1 else 0
    Lc = max(pp - 1, 0)
    Lb = (max(jj * ee - 1, 0) / ee) if ee > 0 else 0
    f1, f2, f3 = 1 - La * aa, 1 - Lc * cc, 1 - Lb * bb
    assert f1 > 0 and f2 > 0 and f3 > 0 and pp * jj < 1
    cntb += (pp == 1 or jj == 1 or jj * ee == 1 or ee == 0 or pp == 0 or jj == 0)
    n10 += 1
    vals = dict(q=qq, p=pp, j=jj, c=cc, e=ee, b=bb, a=aa)
    g = PP.eval_terms(G69, [vals[nm] for nm in names7])
    B0v = (ee + 1) * (pp + 1) * ((jj + 1) * (1 - pp * jj) + ee * (jj * jj + jj + 1 - jj * pp))
    bound = B0v * f1 ** 2 * f2 ** 2 * f3 ** 3
    ok10 &= g >= bound > 0
    r_ = g / bound
    worst = r_ if worst is None or r_ < worst else worst
ck.check("EXACT", "(10) at %d random rational points of the closed hard set of branch (69) (seed 3603; %d of them on "
         "the walls p = 1, j = 1, je = 1 or with p, j or e equal to 0; a, c, b, q vanish with positive "
         "probability): G69 >= B0 (1 - L_a a)^2 (1 - L_c c)^2 (1 - L_b b)^3 > 0 every time" % (n10, cntb),
         ok10 and cntb > 100 and n10 > 2000, "smallest ratio G69/bound = %.6f" % float(worst))

M.crosscheck_with_audit_library(ck)

# ------------------------------------------------------------------------------------------
if NUMERIC:
    import numpy as np
    ck.section("NUMERIC falsification attempts (floating point; evidence, not proof)")
    M.numeric_pointwise_attack(ck, ROLES68 + ROLES69, seed=36001)
    for jx in (68, 69):
        ctrl = M.numeric_min_phi(M.F_VALUES[jx], 36100 + jx, starts=6)
        ck.check("NUMERIC", "control: for the bare non-minor vector f%d the same minimiser finds a negative value of "
                 "the Leibniz sum on SO(4) (so it can detect a violation)" % jx, ctrl < -0.05, "minimum found %.5f" % ctrl)

    rngn = np.random.default_rng(36002)
    worst = np.inf
    for _ in range(50):
        aa, ee = rngn.exponential(1.0), rngn.exponential(1.0)
        co = {GEN69["a"]: aa, GEN69["e"]: ee, GEN69["j"]: 1 + 1 / (aa + ee), 69: 1.0,
              GEN69["q"]: rngn.exponential(0.3), GEN69["b"]: rngn.exponential(0.3)}
        worst = min(worst, M.numeric_min_phi(M.combo_values(co), int(rngn.integers(1 << 30)), starts=4))
    ck.check("NUMERIC", "old-H-sum threshold j = 1 + 1/(a+e) as an inequality on SO(4): 50 random coefficient vectors "
             "on the threshold (seed 36002), minimum of the Leibniz sum (4 BFGS starts each) never below -1e-8",
             worst > -1e-8, "smallest minimum found %.3e" % worst)

    def make_rel(G):
        ex = np.array(list(G.keys()), float)
        cf = np.array([float(v) for v in G.values()])

        def rel(v):
            mon = np.prod(np.array([v[nm] for nm in names7], float)[None, :] ** ex, axis=1)
            return abs(float(cf @ mon)) / float(np.abs(cf) @ np.abs(mon))
        return rel

    G68r, G69r = make_rel(G68), make_rel(G69)
    free = []
    for seed in (36011, 36012):
        free += M.numeric_stratum_search(T, seed, 300)
    d68 = [M.numeric_coefficients(ORDER68, pt["gaps"])[0][5] for pt in free]
    d69 = [M.numeric_coefficients(ORDER69, pt["gaps"])[0][5] for pt in free]
    ck.info("free stratum search (gaps = 0 on T, gaps >= 0, sum of gaps = 1; 600 starts, seeds 36011, 36012): %d "
            "converged points; coefficient of f68 in the basis (68) in [%.3g, %.3g], of f69 in the basis (69) in "
            "[%.3g, %.3g]" % (len(free), min(d68), max(d68), min(d69), max(d69)))
    tol = 1e-6
    for label, order, Gr, seed in (("68", ORDER68, G68r, 36020), ("69", ORDER69, G69r, 36030)):
        hard = M.numeric_hard_points(T, order, 5, seed)
        # keep the points whose representation in this basis is a normal form (all coefficients >= 0)
        hard = [(co, pt) for co, pt in hard if min(co) > -1e-8]
        ds = [co[5] for co, pt in hard]
        n_exact = sum(1 for co, pt in hard if pt["contacts"] == tuple(sorted(T)))
        ck.info("branch (%s): search with the hinge equation d >= target (targets 0.004, 0.01; then 300 perturbed "
                "restarts, seed %d): %d stratum points in normal form (%s) with d > 5e-4, d up to %.4g; contact set "
                "exactly T at %d of them" % (label, seed, len(hard), label, max(ds) if ds else 0.0, n_exact))
        if not hard:
            ck.info("branch (%s): no stratum point with d > 0 found: the radial equation is NOT tested there" % label)
            continue
        maxG = 0.0
        contrast = []
        unsafe = 0
        which = {}
        for co, pt in hard:
            cod = dict(zip(NAMES, co))
            v = {k: max(cod[k], 0.0) / cod["d"] for k in names7}
            maxG = max(maxG, Gr(v))
            contrast.append(Gr({k: x * rngn.uniform(0.9, 1.1) for k, x in v.items()}))
            if label == "68":
                if v["q"] * v["j"] * (v["p"] - 1) >= v["q"] + v["j"] + v["p"] + 1 - tol:
                    key = "four (0,7;1;68)"
                elif v["j"] * v["b"] * (v["e"] - 1) >= v["j"] + v["b"] + v["e"] + 1 - tol:
                    key = "four (7,67;53;68)"
                else:
                    key = "NONE"
                    unsafe += 1
            else:
                if v["p"] * v["j"] >= 1 - tol:
                    key = "SOC"
                elif v["c"] * (v["p"] - 1) >= 1 - tol:
                    key = "old-H"
                elif (v["a"] + v["e"]) * (v["j"] - 1) >= 1 - tol:
                    key = "old-H sum"
                elif v["e"] > 0 and v["b"] > 0 and v["j"] >= (1 / v["e"] + 1 / v["b"]) * (1 - tol):
                    key = "chain"
                else:
                    key = "NONE"
                    unsafe += 1
            which[key] = which.get(key, 0) + 1
        ck.check("NUMERIC", "branch (%s): radial equation G%s = 0 at the %d stratum points with d > 0: |G| <= 1e-5 "
                 "times the sum of the absolute values of its terms (the points are known to about 1e-10/d in "
                 "these coordinates)" % (label, label, len(hard)), maxG < 1e-5,
                 "largest relative |G| = %.3e; for contrast, with every coefficient multiplied by an "
                 "independent random factor in [0.9, 1.1] the relative |G| has median %.2e" % (maxG, float(np.median(contrast))))
        ck.check("NUMERIC", "branch (%s): conclusion of the document at these points: one of the conditions (%s) "
                 "holds" % (label, "3" if label == "68" else "4"), unsafe == 0, str(which))

ck.finish()
