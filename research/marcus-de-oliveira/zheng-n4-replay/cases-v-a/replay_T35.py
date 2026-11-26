"""Independent exact replay of Case T_35 = {2,3,4,6,9,12} (group V) of D. Zheng's proof of the
Marcus-de Oliveira conjecture for n = 4.

Proof document replayed: MOC4_SUPP/proofs/manin_transfer_20261021/final_referee_T34_uniform_20261022.pdf
(4 pages; the document calls the set "T34").  Every object is rebuilt from the paper's
definitions by own code (replay_T35_common.py, replay_T35_common_poly.py).  Nothing from the
supplement is imported or executed.  If MOC4_SUPP is set, the PDF is hashed and three JSON data
files are read for comparison only.

Run:   MOC4_SUPP=/path/to/supplement python replay_T35.py [--no-numeric]
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
T = (2, 3, 4, 6, 9, 12)
DOC = "final_referee_T34_uniform_20261022.pdf"
# normal form (1) of the document: g = q f67 + p f5 + j f53 + c f9 + e f1 + d f68 + b f7 + a f0
ORDER = [67, 5, 53, 9, 1, 68, 7, 0]
NAMES = list("qpjcedba")
GEN = dict(zip(NAMES, ORDER))

ck = M.Checker("Case T_35 = {2,3,4,6,9,12}: replay of final_referee_T34_uniform_20261022.pdf")
M.print_versions(ck)
M.check_document_identity(ck, DOC)

# ------------------------------------------------------------------------------------------
M.check_foundations(ck)
M.check_radial_necessity(ck)
M.check_hessian_necessity(ck)
M.crosscheck_with_supplement_data(ck)

# ------------------------------------------------------------------------------------------
M.check_contact_set(ck, 35, T, "T34")
gens = M.check_generators(ck, T, [0, 1, 2, 5, 7, 9, 53, 67, 68], [68])

ck.section("Normal form (document p. 1, 'Source and full-height necessity')")
M.check_linear_identities(ck, [([68, 2], [0, 9, 67])], "f68 + f2 = f0 + f9 + f67")
ck.check("EXACT", "the nine generators span a space of dimension 8 (= dimension of the contact nullspace), and "
         "the eight generators of the normal form (all but f2) are linearly independent: the coefficients "
         "(q,p,j,c,e,d,b,a) of (1) are uniquely determined linear functions of the gap vector",
         M.rank([M.F_VALUES[j] for j in gens]) == 8 and M.rank([M.F_VALUES[j] for j in ORDER]) == 8
         and len(M.nullspace([M.P[s] for s in T], 14)) == 8)
ck.check("EXACT", "f68 is the only non-minor generator of the normal form; it carries the coefficient d",
         [j for j in ORDER if not M.F_IS_MINOR[j]] == [68] and GEN["d"] == 68)
rng = random.Random(3501)
ok = True
n_minor = n_hard = 0
for _ in range(4000):
    mode = rng.random()
    lam = {j: (Fr(rng.randint(0, 40), rng.randint(1, 9)) if rng.random() < mode else Fr(0)) for j in gens}
    before = M.combo_values(lam)
    s = min(lam[2], lam[68])
    new = dict(lam)
    new[2] -= s
    new[68] -= s
    for j in (0, 9, 67):
        new[j] += s
    ok &= M.combo_values(new) == before and all(v >= 0 for v in new.values())
    if new[68] == 0:
        n_minor += 1
    else:
        ok &= new[2] == 0 and new[68] == lam[68] - lam[2]
        n_hard += 1
ck.check("EXACT", "absorbing min(lambda_2, lambda_68) by the identity keeps the vector and all coefficients >= 0 and "
         "ends either with lambda_68 = 0 (a nonnegative combination of minor vectors) or with lambda_2 = 0 and "
         "d = lambda_68 - lambda_2 > 0 (the normal form (1)); 4000 random nonnegative rational coefficient "
         "vectors, seed 3501", ok and n_minor > 0 and n_hard > 0, "%d minor, %d normal form" % (n_minor, n_hard))
M.check_minor_leibniz(ck, [j for j in gens if M.F_IS_MINOR[j]])

# ------------------------------------------------------------------------------------------
ck.section("Radial necessity in the normal form (document p. 1)")
rd = PP.radial_data(ORDER, NAMES)
R9 = rd["R"]
names7, G1 = PP.check_radial(ck, rd, "T_35", 1475, "final_referee_T34_full_radial_20261022.json")
V9 = {nm: RF(R9, R9[nm]) for nm in NAMES}
q, p, j, c, e, d, b, a = [V9[nm] for nm in NAMES]
Qdoc = c * d + c * e + c * j - c * p + d * d + d * e + d * j + e * j
ck.check("EXACT", "Q = cd + ce + cj - cp + d^2 + de + dj + ej (the document's Q at d = 1)", Qdoc.num == rd["Q"])
cst = PR.terms(rd["Cstar"])
d3 = tuple(3 if nm == "d" else 0 for nm in R9.names)
ck.check("EXACT", "C* has 22 terms, all with coefficient +1, among them d^3: C* >= d^3 > 0 on the closed "
         "orthant with d > 0.  With C(u) = C* + hQ = 0 this forces Q != 0 and h = -C*/Q",
         len(cst) == 22 and all(v == 1 for v in cst.values()) and cst.get(d3) == 1)
assert names7 == ["q", "p", "j", "c", "e", "b", "a"]
R7 = PR(names7)
V7 = {nm: RF(R7, R7[nm]) for nm in names7}
q, p, j, c, e, b, a = [V7[nm] for nm in names7]
G1p = R7.from_terms(G1)
IP = names7.index("p")
kp = PP.coeffs_in(G1, IP)
k0, k1, k2 = kp[0], kp[1], kp[2]
DEG = PP.term_dict_degrees(G1, 7)
ck.check("EXACT", "at d = 1, G has degrees %s in (q,p,j,c,e,b,a); it is quadratic in p: G = k0 + k1 p + k2 p^2 "
         "with %d + %d + %d terms" % (DEG, len(k0), len(k1), len(k2)), sorted(kp) == [0, 1, 2] and DEG[IP] == 2)

# ------------------------------------------------------------------------------------------
ck.section("Universal inequalities used (CITED: audit PROOF.md 7.6); roles and coefficient conditions (p. 1-2)")
ROLES = [
    ("old-H", (GEN["c"], GEN["p"], 68), "condition (3): c(p - d) >= d^2, t = c on f9, main = p on f5"),
    ("chain", (GEN["p"], GEN["b"], GEN["c"], 68), "condition (4): p >= d^2/b + d^2/c, centre f5, leaves f7, f9"),
    ("four", (GEN["b"], GEN["a"], GEN["e"], 68), "condition (5): ab(e - d) >= d^2(a+b+e+d), f7, f0; f1"),
    ("four", (GEN["b"], GEN["q"], GEN["j"], 68), "condition (6): bq(j - d) >= d^2(b+q+j+d), f7, f67; f53"),
]
ck.check("EXACT", "the role tuples are the document's: old-H (9,5,68), chain (5;7,9;68), four (7,0;1;68) and "
         "(7,67;53;68)", [r[1] for r in ROLES] == [(9, 5, 68), (5, 7, 9, 68), (7, 0, 1, 68), (7, 67, 53, 68)])
M.check_roles(ck, ROLES)
M.check_coefficient_forms(ck)
hard_e = (a + 1) * (b + 1) - e * (a * b - 1) - ((a + b + e + 1) - a * b * (e - 1))
hard_j = (q + 1) * (b + 1) - j * (q * b - 1) - ((q + b + j + 1) - q * b * (j - 1))
ck.check("EXACT", "(a+1)(b+1) - e(ab-1) = (a+b+e+1) - ab(e-1), and the same with (q, j): (5) fails iff "
         "L_e e < 1 and (6) fails iff L_j j < 1, with L_e = max((ab-1)/((a+1)(b+1)), 0), L_j likewise",
         hard_e.is_zero() and hard_j.is_zero())

# ------------------------------------------------------------------------------------------
ck.section("Strict constant coefficient (7) and the corner c = 0 (9) (document p. 2)")
G00 = (1 + b) * (c + 1 + b) ** 2 * (a + 1) ** 3 * (q + 1) ** 3
ip, ie, ij, ic, ib = [names7.index(x) for x in "pejcb"]
g00 = PP.set_var(PP.set_var(PP.set_var(G1, ip, 0), ie, 0), ij, 0)
ck.check("EXACT", "G at p = e = j = 0 equals G00 = (1+b)(c+1+b)^2 (a+1)^3 (q+1)^3", R7.from_terms(g00) == G00.num)
k0b0 = PP.set_var(k0, ib, 0)
ck.check("EXACT", "k0 at b = 0 equals (a+1)(q+1)(e+1+a)^2 (c+1)^2 (q+j+1)^2 > 0",
         R7.from_terms(k0b0) == ((a + 1) * (q + 1) * (e + 1 + a) ** 2 * (c + 1) ** 2 * (q + j + 1) ** 2).num)
Gc0 = PP.set_var(G1, ic, 0)
Gcb0 = PP.set_var(Gc0, ib, 0)
ck.check("EXACT", "G at c = b = 0 equals (qp + q + qa + p + pa + 1 + a)(e+1+a)^2 (q+j+1)^2 > 0",
         R7.from_terms(Gcb0) == ((q * p + q + q * a + p + p * a + 1 + a) * (e + 1 + a) ** 2 * (q + j + 1) ** 2).num)

RC = PR(["B", "X", "Y", "C", "E", "J", "P"])
B_, X_, Y_, C_, E_, J_, P_ = RC.gens
one = RC.one
rf = lambda z: RF(RC, z)


def chart8(ab_above, bq_above, bval):
    """Table (8) of the document, with b = bval."""
    ch = {"b": bval, "c": rf(C_), "p": rf(P_)}
    if not ab_above:
        ch["a"] = rf(X_) / (bval * rf(X_ + 1))
        ch["e"] = rf(E_)
    else:
        ch["a"] = rf(X_ + 1) / bval
        ch["e"] = (ch["a"] + 1) * (bval + 1) / rf(X_) * rf(E_) / rf(E_ + 1)
    if not bq_above:
        ch["q"] = rf(Y_) / (bval * rf(Y_ + 1))
        ch["j"] = rf(J_)
    else:
        ch["q"] = rf(Y_ + 1) / bval
        ch["j"] = (ch["q"] + 1) * (bval + 1) / rf(Y_) * rf(J_) / rf(J_ + 1)
    return ch


def chart_ok(ch, ab_above, bq_above):
    """the chart parametrises the stated chamber: products and the hard intervals."""
    ok = True
    ab = ch["a"] * ch["b"]
    bq = ch["q"] * ch["b"]
    ok &= (ab == (rf(X_ + 1) if ab_above else rf(X_) / rf(X_ + 1)))
    ok &= (bq == (rf(Y_ + 1) if bq_above else rf(Y_) / rf(Y_ + 1)))
    if ab_above:      # 1 - L_e e = 1/(1+E)
        ok &= (1 - ch["e"] * (ab - 1) / ((ch["a"] + 1) * (ch["b"] + 1)) == 1 / rf(E_ + 1))
    if bq_above:
        ok &= (1 - ch["j"] * (bq - 1) / ((ch["q"] + 1) * (ch["b"] + 1)) == 1 / rf(J_ + 1))
    return ok


DOC_K0 = {(0, 0): 2227, (0, 1): 3268, (1, 0): 3268, (1, 1): 4945}
DOC_C0 = {(0, 0): 1419, (0, 1): 2110, (1, 0): 2110, (1, 1): 3211}
DOC_B1 = {(0, 0): 3612, (0, 1): 7230, (1, 0): 7230, (1, 1): 12258}
rng = random.Random(3502)
for ab_above in (0, 1):
    for bq_above in (0, 1):
        tag = "chart (8), ab %s 1, bq %s 1" % ("><"[1 - ab_above], "><"[1 - bq_above])
        pos = ["B"] + (["X"] if ab_above else []) + (["Y"] if bq_above else [])
        ch = chart8(ab_above, bq_above, rf(B_))
        ck.check("EXACT", "%s, b = B: every chart value is a quotient of polynomials with nonnegative "
                 "coefficients whose denominator is > 0 for B%s > 0 and the other parameters >= 0; "
                 "ab = %s, bq = %s; on a product-above-one chart 1 - L e = 1/(1+E)"
                 % (tag, "".join(", " + v for v in pos[1:]), "1+X" if ab_above else "X/(1+X)",
                    "1+Y" if bq_above else "Y/(1+Y)"),
                 PP.chart_denominators_positive(RC, ch, pos) and chart_ok(ch, ab_above, bq_above))
        # (7): k0
        Z, dens = PP.cleared(RC, k0, names7, ch)
        W = PP.prod_dens(RC, dens)
        st = PP.poly_stats(Z)
        Z0 = PP.restrict_zero(RC, Z, ["E", "J"])
        W0 = PP.restrict_zero(RC, W, ["E", "J"])
        G00c = (1 + ch["b"]) * (ch["c"] + 1 + ch["b"]) ** 2 * (ch["a"] + 1) ** 3 * (ch["q"] + 1) ** 3
        ratio = rf(W) / rf(W0)
        want = rf((E_ + 1) ** (2 * ab_above) * (J_ + 1) ** (2 * bq_above))
        pt, vals = PP.chart_point(RC, ch, rng, pos)
        direct = PP.eval_terms(k0, [vals[nm] for nm in names7])
        spot = direct * PP.eval_terms(PR.terms(W), pt) == PP.eval_terms(PR.terms(Z), pt)
        ck.check("EXACT", "%s: W k0 (W = product of the positive denominators) has %d terms (document: %d), all "
                 "coefficients >= 0; its part free of E, J equals W0 G00 with W/W0 = (1+E)^%d (1+J)^%d.  Hence "
                 "k0 >= G00 (1 - L_e e)^2 (1 - L_j j)^2 on this chamber: inequality (7)"
                 % (tag, st[0], DOC_K0[(ab_above, bq_above)], 2 * ab_above, 2 * bq_above),
                 st[1] == 0 and st[0] == DOC_K0[(ab_above, bq_above)] and rf(Z0) == G00c * rf(W0)
                 and ratio == want and spot, "spot check of the clearing at a random rational point: %s" % spot)
        # (9): all of G at c = 0, p = P free
        Z, dens = PP.cleared(RC, Gc0, names7, ch)
        W = PP.prod_dens(RC, dens)
        st = PP.poly_stats(Z)
        Z0 = PP.restrict_zero(RC, Z, ["E", "J", "P"])
        W0 = PP.restrict_zero(RC, W, ["E", "J", "P"])
        base9 = (1 + ch["b"]) ** 3 * (ch["a"] + 1) ** 3 * (ch["q"] + 1) ** 3
        ck.check("EXACT", "%s, c = 0, p = P >= 0 free: W G has %d terms (document: %d), all coefficients >= 0; its "
                 "part free of E, J, P equals W0 (1+b)^3 (a+1)^3 (q+1)^3 with W/W0 = (1+E)^%d (1+J)^%d: "
                 "inequality (9)" % (tag, st[0], DOC_C0[(ab_above, bq_above)], 2 * ab_above, 2 * bq_above),
                 st[1] == 0 and st[0] == DOC_C0[(ab_above, bq_above)] and rf(Z0) == base9 * rf(W0)
                 and rf(W) / rf(W0) == want)

# ------------------------------------------------------------------------------------------
ck.section("Bernstein criterion (10) and the region b <= 1 with the old-H threshold P* = 1 + 1/c (p. 2-3)")
RB = PR(["k0", "k1", "k2", "N", "D", "s"])
k0_, k1_, k2_, N_, D_, s_ = RB.gens
lhs = (k0_ * D_ * D_ + k1_ * N_ * s_ * D_ + k2_ * N_ * N_ * s_ * s_) * 2
rhs = (k0_ * D_ * D_ * 2) * (-s_ + 1) ** 2 + (k0_ * D_ * 2 + k1_ * N_) * D_ * s_ * (-s_ + 1) * 2 \
    + (k0_ * D_ * D_ + k1_ * N_ * D_ + k2_ * N_ * N_) * s_ * s_ * 2
ck.check("EXACT", "for p = s N/D: D^2 G(p) = D^2 k0 (1-s)^2 + D (2D k0 + N k1) s(1-s) + (D^2 k0 + N D k1 + N^2 k2) s^2, "
         "so k0 > 0 and nonnegative middle and endpoint numerators give G(p) >= k0 (1 - p/P*)^2 > 0 on "
         "0 <= p < P* = N/D", lhs == rhs)
DEGK = [0 if nm == "p" else dg for nm, dg in zip(names7, DEG)]
for ab_above in (0, 1):
    for bq_above in (0, 1):
        tag = "chart (8), ab %s 1, bq %s 1" % ("><"[1 - ab_above], "><"[1 - bq_above])
        pos = ["B", "C"] + (["X"] if ab_above else []) + (["Y"] if bq_above else [])
        ch = chart8(ab_above, bq_above, rf(B_) / rf(B_ + 1))
        ch["p"] = rf(one)
        thr = 1 + 1 / rf(C_)
        K0, mid, end, W = PP.bernstein(RC, (k0, k1, k2), names7, ch, thr, DEGK)
        s0, s1, s2 = PP.poly_stats(K0), PP.poly_stats(mid), PP.poly_stats(end)
        ck.check("EXACT", "%s, b = B/(1+B) < 1, c = C > 0, P* = (C+1)/C: chart valid; base, middle and endpoint "
                 "numerators have %d, %d, %d terms (document: middle and endpoint %d each), all coefficients >= 0"
                 % (tag, s0[0], s1[0], s2[0], DOC_B1[(ab_above, bq_above)]),
                 PP.chart_denominators_positive(RC, ch, pos) and chart_ok(ch, ab_above, bq_above)
                 and s0[1] == 0 and s1[1] == 0 and s2[1] == 0
                 and s1[0] == s2[0] == DOC_B1[(ab_above, bq_above)])

# ------------------------------------------------------------------------------------------
ck.section("Shared-coefficient allocation for b >= 1 (document p. 3)")
RA = PR(["a", "b", "q", "c", "dl", "d1", "d2"])
a_, b_, q_, c_, dl_, d1_, d2_ = [RF(RA, g) for g in RA.gens]


def Ecal(x, bb, dl):
    return dl * (bb + 1) * (x + dl) / (x * bb - dl)


# derivative of E_a(delta) with respect to delta, by the quotient rule on num/den
num = dl_ * (b_ + 1) * (a_ + dl_)
den = a_ * b_ - dl_
dnum = (b_ + 1) * (a_ + 2 * dl_)
dden = RF(RA, -RA.one)
deriv = (dnum * den - num * dden) / (den * den)
ck.check("EXACT", "(11): d/d delta of E_a(delta) = delta(b+1)(a+delta)/(ab-delta) equals "
         "(b+1)(a^2 b + 2ab delta - delta^2)/(ab-delta)^2; the numerator is a^2 b + delta(2ab - delta) > 0 for "
         "a, b > 0, 0 <= delta < ab, so E_a is strictly increasing there, E_a(0) = 0",
         deriv == (b_ + 1) * (a_ * a_ * b_ + 2 * a_ * b_ * dl_ - dl_ * dl_) / ((a_ * b_ - dl_) ** 2)
         and (a_ * a_ * b_ + 2 * a_ * b_ * dl_ - dl_ * dl_) == a_ * a_ * b_ + dl_ * (2 * a_ * b_ - dl_))
E1 = (b_ + 1) * (a_ + 1) / (a_ * b_ - 1)
ck.check("EXACT", "E_a(1) = (a+1)(b+1)/(ab-1): for ab > 1, E_a maps [0,1) onto the hard interval "
         "[0, (a+1)(b+1)/(ab-1)); for ab <= 1 the numerator at delta = ab is ab(b+1)(a+ab) > 0 while the "
         "denominator tends to 0+, so E_a maps [0, ab) onto [0, infinity)",
         Ecal(a_, b_, RF(RA, RA.one)) == E1
         and (dl_ * (b_ + 1) * (a_ + dl_)).num == (RA["dl"] * (RA["b"] + 1) * (RA["a"] + RA["dl"])))
Ea = Ecal(a_, b_, dl_)
bd13 = a_ * (dl_ * b_) * (Ea - dl_) - dl_ * dl_ * (a_ + dl_ * b_ + Ea + dl_)
ck.check("EXACT", "(13): a (delta b)(E_a(delta) - delta) - delta^2 (a + delta b + E_a(delta) + delta) = 0: the piece "
         "a f0 + (delta b) f7 + E_a(delta) f1 + delta f68 lies on the boundary of the four-coefficient "
         "condition (roles (7,0;1;68)); the same with (q, f67, f53) for roles (7,67;53;68)", bd13.is_zero())
# bookkeeping of the allocation: pieces add up to the normal form
r_ = 1 - d1_ - d2_
piece1 = {0: a_, 7: d1_ * b_, 1: Ecal(a_, b_, d1_), 68: d1_}
piece2 = {67: q_, 7: d2_ * b_, 53: Ecal(q_, b_, d2_), 68: d2_}
chainp = {5: r_ * r_ / c_ + r_ / b_, 7: r_ * b_, 9: c_, 68: r_}
tot = {}
for pc in (piece1, piece2, chainp):
    for k, v in pc.items():
        tot[k] = tot.get(k, 0) + v
want = {0: a_, 67: q_, 7: b_, 68: RF(RA, RA.one), 1: Ecal(a_, b_, d1_), 53: Ecal(q_, b_, d2_),
        9: c_, 5: r_ * r_ / c_ + r_ / b_}
chain_cond = chainp[5] - (r_ * r_ / (r_ * b_) + r_ * r_ / c_)
ck.check("EXACT", "allocation: [a f0 + d1 b f7 + E_a(d1) f1 + d1 f68] + [q f67 + d2 b f7 + E_q(d2) f53 + d2 f68] + "
         "[P* f5 + r b f7 + c f9 + r f68] with r = 1 - d1 - d2 and P* = r^2/c + r/b has total coefficients "
         "(a, q, b, 1, e, j, c, P*) on (f0, f67, f7, f68, f1, f53, f9, f5); the third piece is on the boundary "
         "of the chain condition (centre f5: P* = r^2/(rb) + r^2/c), which is (14)",
         all((tot[k] - want[k]).is_zero() for k in want) and set(tot) == set(want) and chain_cond.is_zero())

RK = PR(["B", "U", "V", "X", "Y", "C"])
Bk, Uk, Vk, Xk, Yk, Ck = RK.gens
rk = lambda z: RF(RK, z)
bK = rk(Bk + 1)
DK = rk(Uk + Vk + 1)
chK = {"b": bK, "c": rk(Ck), "p": rk(RK.one)}
chK["a"] = (rk(Uk) + rk(Xk) * DK) / (bK * DK)
chK["q"] = (rk(Vk) + rk(Yk) * DK) / (bK * DK)
d1K, d2K, rK = rk(Uk) / DK, rk(Vk) / DK, 1 / DK
chK["e"] = Ecal(chK["a"], bK, d1K)
chK["j"] = Ecal(chK["q"], bK, d2K)
PstarK = rK * rK / rk(Ck) + rK / bK
ck.check("EXACT", "chart (15): b = 1+B, D = 1+U+V, delta1 = U/D, delta2 = V/D, r = 1/D, a = (U+XD)/(bD), "
         "q = (V+YD)/(bD), c = C gives ab - delta1 = X, bq - delta2 = Y, e = E_a(delta1) = U(b+1)(XD+U(1+b))/(bD^2 X), "
         "j = E_q(delta2) = V(b+1)(YD+V(1+b))/(bD^2 Y) and P* = r^2/c + r/b = (b+CD)/(CbD^2)  (16); all chart "
         "values have nonnegative coefficients and denominators > 0 for X, Y, C > 0",
         (chK["a"] * bK - d1K) == rk(Xk) and (chK["q"] * bK - d2K) == rk(Yk)
         and chK["e"] == rk(Uk) * (bK + 1) * (rk(Xk) * DK + rk(Uk) * (1 + bK)) / (bK * DK * DK * rk(Xk))
         and chK["j"] == rk(Vk) * (bK + 1) * (rk(Yk) * DK + rk(Vk) * (1 + bK)) / (bK * DK * DK * rk(Yk))
         and PstarK == (bK + rk(Ck) * DK) / (rk(Ck) * bK * DK * DK)
         and PP.chart_denominators_positive(RK, chK, ["X", "Y", "C"])
         and PP.positive_on_domain(RK, PstarK.den, ["C"]) and PP.nonneg(PstarK.num))
K0, mid, end, W = PP.bernstein(RK, (k0, k1, k2), names7, chK, PstarK, DEGK)
s0, s1, s2 = PP.poly_stats(K0), PP.poly_stats(mid), PP.poly_stats(end)
ck.check("EXACT", "allocation chamber b >= 1, r > 0, threshold (16): base, middle and endpoint numerators have "
         "%d, %d, %d terms (document: 117110, 209489, 289848), all coefficients >= 0"
         % (s0[0], s1[0], s2[0]),
         s0[1] == 0 and s1[1] == 0 and s2[1] == 0 and (s0[0], s1[0], s2[0]) == (117110, 209489, 289848))
rng = random.Random(3503)
spot = True
for _ in range(3):
    pt, vals = PP.chart_point(RK, chK, rng, ["X", "Y", "C"])
    kv = [PP.eval_terms(k, [vals[nm] for nm in names7]) for k in (k0, k1, k2)]
    Wv = PP.eval_terms(PR.terms(W), pt)
    Nv, Dv = PP.eval_terms(PR.terms(PstarK.num), pt), PP.eval_terms(PR.terms(PstarK.den), pt)
    spot &= PP.eval_terms(PR.terms(K0), pt) == Wv * kv[0]
    spot &= PP.eval_terms(PR.terms(mid), pt) == Wv * (2 * Dv * kv[0] + Nv * kv[1])
    spot &= PP.eval_terms(PR.terms(end), pt) == Wv * (Dv * Dv * kv[0] + Nv * Dv * kv[1] + Nv * Nv * kv[2])
ck.check("EXACT", "the three cleared numerators equal W times the Bernstein numerators computed directly from "
         "k0, k1, k2 at 3 random rational chart points (guards the substitution code)", spot)

# ------------------------------------------------------------------------------------------
ck.section("Exact evaluations on the closed hard region, boundaries included (falsification attempts)")


def kvals(vals):
    pt = [vals[nm] for nm in names7]
    return [PP.eval_terms(k, pt) for k in (k0, k1, k2)]


def rq(rng, zero=0.0, lo=1, hi=30):
    if rng.random() < zero:
        return Fr(0)
    return Fr(rng.randint(lo, hi), rng.randint(lo, hi)) * rng.choice([1, 1, 1, 5, Fr(1, 5)])


def hard_coeff(rng, x, bb):
    """random e in the hard interval of the pair (x, b)."""
    if x * bb > 1:
        emax = (x + 1) * (bb + 1) / (x * bb - 1)
        return emax * Fr(rng.randint(0, 99), 100)
    return rq(rng, 0.2)


rng = random.Random(3504)
ok7 = ok9 = okb = True
worst7 = None
cnt = 0
for _ in range(1500):
    mode = rng.random()
    bb = rng.choice([Fr(0), Fr(1), rq(rng), rq(rng)]) if mode < 0.5 else rq(rng)
    aa = rng.choice([Fr(0), (1 / bb if bb else Fr(1)), rq(rng), rq(rng)])
    qq = rng.choice([Fr(0), (1 / bb if bb else Fr(1)), rq(rng), rq(rng)])
    cc = rq(rng, 0.15)
    ee, jj = hard_coeff(rng, aa, bb), hard_coeff(rng, qq, bb)
    vals = dict(a=aa, b=bb, q=qq, c=cc, e=ee, j=jj, p=Fr(0))
    Le = max((aa * bb - 1) / ((aa + 1) * (bb + 1)), 0)
    Lj = max((qq * bb - 1) / ((qq + 1) * (bb + 1)), 0)
    assert Le * ee < 1 and Lj * jj < 1
    kv = kvals(vals)
    bound = (1 + bb) * (cc + 1 + bb) ** 2 * (aa + 1) ** 3 * (qq + 1) ** 3 * (1 - Le * ee) ** 2 * (1 - Lj * jj) ** 2
    ok7 &= kv[0] >= bound > 0
    ratio = kv[0] / bound
    worst7 = ratio if worst7 is None or ratio < worst7 else worst7
    # (9) at c = 0 with random p
    pp = rq(rng, 0.2)
    v0 = dict(vals, c=Fr(0), p=pp)
    g = PP.eval_terms(G1, [v0[nm] for nm in names7])
    ok9 &= g >= (1 + bb) ** 3 * (aa + 1) ** 3 * (qq + 1) ** 3 * (1 - Le * ee) ** 2 * (1 - Lj * jj) ** 2 > 0
    # b <= 1, c > 0: Bernstein with P* = 1 + 1/c
    if bb <= 1 and cc > 0:
        cnt += 1
        Ps = 1 + 1 / cc
        okb &= kv[0] > 0 and 2 * kv[0] + Ps * kv[1] >= 0 and kv[0] + Ps * kv[1] + Ps * Ps * kv[2] >= 0
        t = Fr(rng.randint(0, 99), 100)
        okb &= kv[0] + kv[1] * (t * Ps) + kv[2] * (t * Ps) ** 2 > 0
ck.check("EXACT", "(7) at 1500 random rational points of the closed hard set {L_e e < 1, L_j j < 1} (seed 3504; "
         "b = 0, b = 1, a = 0, q = 0, ab = 1, bq = 1, c = 0, e = 0, j = 0 occur with positive probability): "
         "k0 >= G00 (1 - L_e e)^2 (1 - L_j j)^2 > 0 every time", ok7, "smallest ratio k0/bound = %.6f" % float(worst7))
ck.check("EXACT", "(9) at the same points with c = 0 and random p >= 0: G >= (1+b)^3 (a+1)^3 (q+1)^3 "
         "(1 - L_e e)^2 (1 - L_j j)^2 > 0 every time", ok9)
ck.check("EXACT", "b <= 1, c > 0 at those of the points with b <= 1 (%d points, b = 0 and b = 1 included): "
         "k0 > 0, 2k0 + P* k1 >= 0, k0 + P* k1 + P*^2 k2 >= 0 for P* = 1 + 1/c, and G(p) > 0 at a random "
         "p in [0, P*)" % cnt, okb and cnt > 300)

rng = random.Random(3505)
oka = True
cnt = cnt0 = 0
for _ in range(1200):
    bb = 1 + rq(rng, 0.25)
    aa = rq(rng, 0.2)
    qq = rq(rng, 0.2)
    cc = rq(rng, 0.0)
    # allocation parameters; delta = 0 is forced when the side coefficient is 0
    lim1 = min(Fr(1), aa * bb)
    lim2 = min(Fr(1), qq * bb)
    dd1 = lim1 * Fr(rng.randint(0, 99), 100) * rng.choice([0, 1, 1, 1])
    dd2 = lim2 * Fr(rng.randint(0, 99), 100) * rng.choice([0, 1, 1, 1])
    if dd1 + dd2 >= 1:
        continue
    ee = dd1 * (bb + 1) * (aa + dd1) / (aa * bb - dd1) if aa > 0 else rq(rng, 0.3)
    jj = dd2 * (bb + 1) * (qq + dd2) / (qq * bb - dd2) if qq > 0 else rq(rng, 0.3)
    rr = 1 - dd1 - dd2
    Ps = rr * rr / cc + rr / bb
    vals = dict(a=aa, b=bb, q=qq, c=cc, e=ee, j=jj, p=Fr(0))
    kv = kvals(vals)
    cnt += 1
    cnt0 += (aa == 0 or qq == 0)
    oka &= kv[0] > 0 and 2 * kv[0] + Ps * kv[1] >= 0 and kv[0] + Ps * kv[1] + Ps * Ps * kv[2] >= 0
    t = Fr(rng.randint(0, 99), 100)
    oka &= kv[0] + kv[1] * (t * Ps) + kv[2] * (t * Ps) ** 2 > 0
ck.check("EXACT", "b >= 1, c > 0, r = 1 - delta1 - delta2 > 0 at %d random rational points given by (a, b, q, c, "
         "delta1, delta2) with e = E_a(delta1), j = E_q(delta2) (seed 3505; b = 1, delta = 0 and the limit "
         "layers a = 0 or q = 0 with delta = 0 and arbitrary e or j, %d points, are included): k0 > 0, "
         "2k0 + P* k1 >= 0, k0 + P* k1 + P*^2 k2 >= 0 for P* = r^2/c + r/b, and G(p) > 0 at a random p in "
         "[0, P*)" % (cnt, cnt0), oka and cnt > 500 and cnt0 > 100)

# ------------------------------------------------------------------------------------------
M.crosscheck_with_audit_library(ck)

# ------------------------------------------------------------------------------------------
if NUMERIC:
    import numpy as np
    ck.section("NUMERIC falsification attempts (floating point; evidence, not proof)")
    M.numeric_pointwise_attack(ck, ROLES, seed=35001)
    for jx in (68,):
        ctrl = M.numeric_min_phi(M.F_VALUES[jx], 35100 + jx, starts=6)
        ck.check("NUMERIC", "control: for the bare non-minor vector f%d the same minimiser finds a negative value of "
                 "the Leibniz sum on SO(4) (so it can detect a violation)" % jx, ctrl < -0.05, "minimum found %.5f" % ctrl)

    # (i) the coefficient thresholds as inequalities on SO(4): boundary points of each condition
    rngn = np.random.default_rng(35002)

    def vec(co):
        return M.combo_values({GEN[k]: float(v) for k, v in co.items()})

    worst = np.inf
    n_pts = 0
    for _ in range(60):
        bb = 1.0 + rngn.exponential(1.0) * rngn.choice([0.0, 1.0, 1.0])
        aa, qq, cc = rngn.exponential(1.0), rngn.exponential(1.0), rngn.exponential(1.0) + 0.05
        d1 = min(1.0, aa * bb) * rngn.uniform(0, 0.98) * rngn.choice([0.0, 1.0, 1.0])
        d2 = min(1.0, qq * bb) * rngn.uniform(0, 0.98) * rngn.choice([0.0, 1.0, 1.0])
        if d1 + d2 >= 1:
            s = (d1 + d2) * 1.02
            d1, d2 = d1 / s, d2 / s
        ee = d1 * (bb + 1) * (aa + d1) / (aa * bb - d1)
        jj = d2 * (bb + 1) * (qq + d2) / (qq * bb - d2)
        rr = 1 - d1 - d2
        co = dict(a=aa, b=bb, q=qq, c=cc, e=ee, j=jj, d=1.0, p=rr * rr / cc + rr / bb)
        worst = min(worst, M.numeric_min_phi(vec(co), int(rngn.integers(1 << 30)), starts=4))
        n_pts += 1
    ck.check("NUMERIC", "allocation threshold (14) as an inequality on SO(4): at %d random coefficient vectors with "
             "b >= 1 and p exactly on the threshold P* = r^2/c + r/b (e = E_a(delta1), j = E_q(delta2), seed 35002), "
             "the minimum of the Leibniz sum over SO(4) (4 BFGS starts each) is never below -1e-8" % n_pts,
             worst > -1e-8, "smallest minimum found %.3e" % worst)
    worst = np.inf
    for _ in range(40):
        cc = rngn.exponential(1.0) + 0.05
        co = dict(a=rngn.exponential(0.5), b=rngn.uniform(0, 1), q=rngn.exponential(0.5), c=cc, e=rngn.exponential(0.5),
                  j=rngn.exponential(0.5), d=1.0, p=1 + 1 / cc)
        worst = min(worst, M.numeric_min_phi(vec(co), int(rngn.integers(1 << 30)), starts=4))
    ck.check("NUMERIC", "old-H threshold p = 1 + 1/c (condition (3)) as an inequality on SO(4): 40 random coefficient "
             "vectors on the threshold, minimum of the Leibniz sum never below -1e-8", worst > -1e-8,
             "smallest minimum found %.3e" % worst)

    # (ii) the hard region in floating point: G > 0 below the selected threshold
    Gex = np.array(list(G1.keys()), float)
    Gcf = np.array([float(v) for v in G1.values()])

    def Gnum(v):
        x = np.array([v[nm] for nm in names7], float)
        return float(Gcf @ np.prod(x[None, :] ** Gex, axis=1))

    def solve_delta(x, bb, ee):
        """delta in [0, min(1, x b)) with E_x(delta) = ee (positive root of the quadratic)."""
        if x * bb <= 0 or ee <= 0:
            return 0.0
        A2, B2, C2 = (bb + 1), (x * (bb + 1) + ee), -ee * x * bb
        return (-B2 + np.sqrt(B2 * B2 - 4 * A2 * C2)) / (2 * A2)

    rngn = np.random.default_rng(35003)
    minratio = np.inf
    nb = 0
    for _ in range(4000):
        bb = 1.0 + rngn.exponential(2.0) * rngn.choice([0.0, 1.0, 1.0, 1.0])
        aa = rngn.exponential(1.0) * rngn.choice([0.0, 1e-6, 1.0, 1.0, 1.0])
        qq = rngn.exponential(1.0) * rngn.choice([0.0, 1e-6, 1.0, 1.0, 1.0])
        cc = rngn.exponential(1.0) + 1e-3
        ee = rngn.exponential(1.0) * rngn.choice([0.0, 1.0, 1.0])
        jj = rngn.exponential(1.0) * rngn.choice([0.0, 1.0, 1.0])
        if aa * bb > 1:
            ee = (aa + 1) * (bb + 1) / (aa * bb - 1) * rngn.uniform(0, 1)
        if qq * bb > 1:
            jj = (qq + 1) * (bb + 1) / (qq * bb - 1) * rngn.uniform(0, 1)
        d1, d2 = solve_delta(aa, bb, ee), solve_delta(qq, bb, jj)
        rr = 1 - d1 - d2
        if rr <= 1e-9:
            continue
        Ps = rr * rr / cc + rr / bb
        tt = rngn.uniform(0, 1) ** 0.3 * 0.999999
        v = dict(a=aa, b=bb, q=qq, c=cc, e=ee, j=jj, p=tt * Ps)
        g = Gnum(v)
        g0 = Gnum(dict(v, p=0.0))
        minratio = min(minratio, g / (g0 * (1 - tt) ** 2))
        nb += 1
    ck.check("NUMERIC", "b >= 1 in original coordinates (delta1, delta2 obtained by solving (12) in floating point, "
             "seed 35003): at %d random hard points with r > 0 and p = t P*, t in [0,1), G(p) >= k0 (1-t)^2 "
             "(ratio never below 1 - 1e-6); a = 0, q = 0, e = 0, j = 0, b = 1 and near-zero side coefficients "
             "are sampled" % nb, minratio > 1 - 1e-6 and nb > 2000, "smallest ratio G/(k0 (1-t)^2) = %.6f" % minratio)

    # (iii) realizable points of the stratum with d > 0
    Gab = np.abs(Gcf)

    def Grel(v):
        x = np.array([v[nm] for nm in names7], float)
        mon = np.prod(x[None, :] ** Gex, axis=1)
        return abs(float(Gcf @ mon)) / float(Gab @ np.abs(mon))

    free = []
    for seed in (35011, 35012):
        free += M.numeric_stratum_search(T, seed, 300)
    dfree = [M.numeric_coefficients(ORDER, pt["gaps"])[0][NAMES.index("d")] for pt in free]
    ck.info("free stratum search (gaps = 0 on T, gaps >= 0, sum of gaps = 1; 600 starts, seeds 35011, 35012): %d "
            "converged points (max residual 1e-10); their coefficient d lies in [%.3g, %.3g]: the free search "
            "finds no point of the hard branch" % (len(free), min(dfree), max(dfree)))
    hard = M.numeric_hard_points(T, ORDER, NAMES.index("d"), 35020)
    ds = [co[NAMES.index("d")] for co, pt in hard]
    n_exact = sum(1 for co, pt in hard if pt["contacts"] == tuple(sorted(T)))
    ck.info("search with the hinge equation d >= target (targets 0.004, 0.01; then 300 perturbed restarts, seed "
            "35020): %d stratum points with d > 5e-4, d up to %.4g (normalisation: sum of gaps = 1); the contact "
            "set is exactly T at %d of them" % (len(hard), max(ds) if ds else 0.0, n_exact))
    if hard:
        maxG = 0.0
        contrast = []
        unsafe = 0
        which = {}
        mincoef = 0.0
        for co, pt in hard:
            cod = dict(zip(NAMES, co))
            dd = cod["d"]
            mincoef = min(mincoef, min(co) / dd)
            v = {k: max(cod[k], 0.0) / dd for k in names7}
            maxG = max(maxG, Grel(v))
            contrast.append(Grel({k: x * rngn.uniform(0.9, 1.1) for k, x in v.items()}))
            aa, bb, qq, cc, ee, jj, ppv = v["a"], v["b"], v["q"], v["c"], v["e"], v["j"], v["p"]
            tol = 1e-6
            if aa * bb * (ee - 1) >= aa + bb + ee + 1 - tol:
                key = "(5)"
            elif bb * qq * (jj - 1) >= bb + qq + jj + 1 - tol:
                key = "(6)"
            elif bb <= 1 + 1e-9:
                if cc > 0 and cc * (ppv - 1) >= 1 - tol:
                    key = "(3), b <= 1"
                else:
                    key = "NONE (b <= 1)"
                    unsafe += 1
            else:
                d1, d2 = solve_delta(aa, bb, ee), solve_delta(qq, bb, jj)
                rr = 1 - d1 - d2
                if rr <= tol:
                    key = "two four pieces"
                elif cc > 0 and ppv >= (rr * rr / cc + rr / bb) * (1 - tol):
                    key = "(14), b > 1"
                else:
                    key = "NONE (b > 1)"
                    unsafe += 1
            which[key] = which.get(key, 0) + 1
        ck.check("NUMERIC", "normal form at the %d points: all eight coefficients are >= -1e-6 d" % len(hard),
                 mincoef > -1e-6, "smallest coefficient/d = %.2e" % mincoef)
        ck.check("NUMERIC", "radial equation G = 0 at the %d stratum points with d > 0: |G| <= 1e-5 times the sum of "
                 "the absolute values of its 1475 terms, at the normalised coefficients (the points are known to "
                 "about 1e-10/d in these coordinates)" % len(hard), maxG < 1e-5,
                 "largest relative |G| = %.3e; for contrast, with every coefficient multiplied by an "
                 "independent random factor in [0.9, 1.1] the relative |G| has median %.2e" % (maxG, float(np.median(contrast))))
        ck.check("NUMERIC", "conclusion of the document at these points: each satisfies (5) or (6); else for b <= 1 "
                 "the old-H condition (3); else for b > 1 either delta1 + delta2 >= 1 (two four pieces) or the "
                 "allocation threshold (14)", unsafe == 0, str(which))
    else:
        ck.info("no stratum point with d > 0 was found: the radial equation is NOT tested on realizable points")

ck.finish()
