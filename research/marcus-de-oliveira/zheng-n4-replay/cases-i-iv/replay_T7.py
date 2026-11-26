"""Independent replay of Case T_7 = {0,1,2,4,7,8} (Zheng, MOC n = 4; group II).

Document replayed: root_moc4_T6_uniform_20261022.pdf (3 pages; historical label "T6").
Nothing from the supplement is imported or executed.  If MOC4_SUPP is set, the PDF is hashed
(sha256) to pin the document; no other use.  All objects are re-derived in replay_T7_common.py
from the paper's definitions.

Run:  MOC4_SUPP=/path/to/supplement python replay_T7.py     (exit code 0 iff every check passes)
Packages: sympy (exact part), numpy + scipy (NUMERIC part only).

Sections 1-5 are EXACT (integers, rationals, rational functions over Q).  Section 6 is NUMERIC
(falsification attempts; seeds and sample sizes are printed).

Method for the "kernel" identities.  In the chart all parameters except v are real, and
kappa (F_j - F_0) = i mu W_j(v) with W_j a rational function of v with real coefficients.  Hence
g_j = -mu Im W_j(v) and, with rho = mu Im v,
        g_j / rho = -(W_j(V) - W_j(Vb)) / (V - Vb)      at V = v, Vb = conj(v).
The right-hand side is computed in the field Q(h,k,l,s,V,Vb) with V, Vb independent; an identity
there holds for all complex v.  |p(v)|^2 becomes p(V) p(Vb) for a real polynomial p.  Section 5
re-checks this at exact Gaussian-rational points directly from the 24 products.
"""
import math
import os
import random
import sys
from fractions import Fraction

sys.dont_write_bytecode = True           # keep the output directory free of __pycache__
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import replay_T7_common as H            # noqa: E402
from replay_T7_common import GQ         # noqa: E402

ck = H.Checks("replay_T7: Case T_7 = {0,1,2,4,7,8}, document root_moc4_T6_uniform_20261022.pdf")
T = H.PAPER_T[7]
DOC_T = (0, 1, 2, 4, 7, 8)                      # the document's "T6"
DOC_GENS = [0, 1, 2, 3, 4, 5, 6, 11, 21, 22, 23, 26]
A_COEF = {3: -1, 9: 1, 16: 1}                   # A = -g3 + g9 + g16
B_COEF = {3: -1, 6: 1, 9: 1, 12: -1, 13: 1, 18: 1}

# ---------------------------------------------------------------------------------------------
ck.section("1. Objects (EXACT)")
H.document_hash_check(ck, "root_moc4_T6_uniform_20261022.pdf")
H.facet_table_sanity(ck)
H.crosscheck_with_audit_lib(ck)
ck.check("document's contact set 'T6' = paper's T_7 (identical set, same indexing; group element = identity)",
         tuple(sorted(DOC_T)) == tuple(sorted(T)))
ck.check("the reference permutation pi_0 = identity belongs to the contact set", 0 in T)

# ---------------------------------------------------------------------------------------------
ck.section("2. Contact cone and exact minor-cone criterion, doc p.1 (EXACT)")
gens, minor, nonminor, facets, cols = H.minor_cone_report(
    ck, T, DOC_GENS, [22, 26], None, {"A": A_COEF, "B": B_COEF}, 14)
allv = [H.F_VALUES[j] for j in gens]
ident2 = H.functional_values({19: 1, 6: -1, 12: 1, 13: -1, 18: -1}, allv)
ck.check("identity (2): g19 = g6 - g12 + g13 + g18 on all 12 generators", all(x == 0 for x in ident2))
ck.check("f_22 and f_26 are positive at permutation 3 (so g_3 = 0 kills both non-minor coefficients)",
         H.F_VALUES[22][3] > 0 and H.F_VALUES[26][3] > 0,
         "f_22(pi_3) = %d, f_26(pi_3) = %d" % (H.F_VALUES[22][3], H.F_VALUES[26][3]))
ck.check("A(f_22) = 1, A(f_26) = -1, B(f_22) = -1, B(f_26) = 1",
         H.functional_values(A_COEF, [H.F_VALUES[22], H.F_VALUES[26]]) == [1, -1]
         and H.functional_values(B_COEF, [H.F_VALUES[22], H.F_VALUES[26]]) == [-1, 1])

# ---------------------------------------------------------------------------------------------
ck.section("3. The chart is exhaustive: identities for general labels, doc p.2 (EXACT)")
K, g = H.make_field("z,w,v,r,c")
z, w, v, r, c = [g[n] for n in ["z", "w", "v", "r", "c"]]
# A = (inf, 0, 1, z), C = (w, v, r, c); F_pi = -C_pi(1) (1 - C_pi(2)) (z - C_pi(3))
F = H.products([None, K(0), K(1), z], [w, v, r, c], inf_A=0)
C_ = [w, v, r, c]
ck.check("products agree with the document's formula F_pi = -C_pi(1)(1-C_pi(2))(z-C_pi(3))",
         all(F[i] == -C_[p[1]] * (1 - C_[p[2]]) * (z - C_[p[3]]) for i, p in enumerate(H.PERMS)))
ck.check("F0 - F1 = v (r - c)(z - 1)   (nonzero for distinct labels)", F[0] - F[1] == v * (r - c) * (z - 1))
hh = z * (v - 1) / (v * (z - 1))
ck.check("(F2 - F4)/(F0 - F1) = z(v-1)/(v(z-1)) =: h   (so contacts 0,1,2,4 make h real)",
         (F[2] - F[4]) / (F[0] - F[1]) == hh)
ck.check("v - z = (h - 1) v (z - 1)   (so h = 1 iff v = z)", v - z == (hh - 1) * v * (z - 1))
ck.check("F0 - F2 = (z - c)(r - v)   (so kappa = i mu/((z-c)(v-r)) gives kappa(F2 - F0) = i mu)",
         F[0] - F[2] == (z - c) * (r - v))
tt = (c - z) * (v - r) / ((c - r) * (v - z))
ck.check("(F0 - F2)/(F0 - F1) = -(h - 1) t with t = (c-z)(v-r)/((c-r)(v-z))   (so t is real)",
         (F[0] - F[2]) / (F[0] - F[1]) == -(hh - 1) * tt)
ck.check("1 - t = (v - c)(z - r)/((c - r)(v - z))   (so t = 1 iff c = v or z = r)",
         1 - tt == (v - c) * (z - r) / ((c - r) * (v - z)))
den = (z - c) * (v - r)
ss = r * (w - v) / (v - r)
ck.check("contact 8: (F8 - F0)/((z-c)(v-r)) = 1 + s with s = r(w-v)/(v-r)   (so s is real)",
         (F[8] - F[0]) / den == 1 + ss)
ck.check("contact 7: (F7 - F0)/((z-c)(v-r)) = -s (1-c)(z-r)/(r(z-c)) - (F0-F1)/((z-c)(v-r))",
         (F[7] - F[0]) / den == -ss * (1 - c) * (z - r) / (r * (z - c)) - (F[0] - F[1]) / den)
ck.info("=> given s real and nonzero, contact 7 holds iff (1-c)(z-r)/(r(z-c)) is real")

K2, g2 = H.make_field("h,t,v,r")
h, t, v, r = [g2[n] for n in ["h", "t", "v", "r"]]
z = h * v / (1 + (h - 1) * v)
c = (z * (v - r) - t * r * (v - z)) / ((v - r) - t * (v - z))
ck.check("solved chart: z = hv/(1+(h-1)v) gives back h; c-formula gives back t",
         z * (v - 1) / (v * (z - 1)) == h and (c - z) * (v - r) / ((c - r) * (v - z)) == t)
k = 1 + t * (h - 1)
ck.check("(1-c)(z-r)/(r(z-c)) = (k/r - 1/v)/(k-1) - 1 with k = 1 + t(h-1)",
         (1 - c) * (z - r) / (r * (z - c)) == (k / r - 1 / v) / (k - 1) - 1)
ck.check("D = h - k = (h-1)(1-t)   (nonzero since h != 1, t != 1)", h - k == (h - 1) * (1 - t))
ck.info("=> l := k/r - 1/v real, i.e. r = kv/(1+lv);  w = v + s(v-r)/r;  t = (k-1)/(h-1)")

# ---------------------------------------------------------------------------------------------
ck.section("4. Kernel identities (4)-(7), ratios and the interpolation identity, doc p.2-3 (EXACT)")
K3, g3 = H.make_field("h,k,l,s,V,Vb")
h, k, l, s, V, Vb = [g3[n] for n in ["h", "k", "l", "s", "V", "Vb"]]
t = (k - 1) / (h - 1)


def chart(v):
    z = h * v / (1 + (h - 1) * v)
    r = k * v / (1 + l * v)
    w = v + s * (v - r) / r
    c = (z * (v - r) - t * r * (v - z)) / ((v - r) - t * (v - z))
    F = H.products([None, K3(0), K3(1), z], [w, v, r, c], inf_A=0)
    den = (z - c) * (v - r)
    return [(F[j] - F[0]) / den for j in range(24)]      # W_j, with kappa (F_j - F_0) = i mu W_j


W1 = chart(V)
W2 = chart(Vb)
G = [-(a - b) / (V - Vb) for a, b in zip(W1, W2)]        # G_j = g_j / rho, rho = mu Im v
ck.check("the chart satisfies all six contacts: g_j = 0 identically for j in T_7",
         all(G[j] == 0 for j in T))
ck.check("no other gap vanishes identically on the chart",
         all(G[j] != 0 for j in range(24) if j not in T))
D = h - k
m = h - l - 1
n = k - l - 1
alpha = h * k - h * l - k
M1 = alpha * V - D
M2 = alpha * Vb - D
ck.check("(4a) g3 = rho h D n k / ((k-1)|M(v)|^2)", G[3] == h * D * n * k / ((k - 1) * M1 * M2))
ck.check("(4b) g6 = -rho s / (k |v|^2)", G[6] == -s / (k * V * Vb))
ck.check("(4c) g10 = rho h s n / ((k-1) k |v-1|^2)", G[10] == h * s * n / ((k - 1) * k * (V - 1) * (Vb - 1)))
ck.check("(5a) e := g14 = -rho (k + l s)/k", G[14] == -(k + l * s) / k)
ck.check("(5b) g18 = -D/(k-1) e", G[18] == -D / (k - 1) * G[14])
ck.check("(5c) g9 = (1 + s m / D) g3", G[9] == (1 + s * m / D) * G[3])
Aexp = -G[3] + G[9] + G[16]
Bexp = -G[3] + G[6] + G[9] - G[12] + G[13] + G[18]
ck.check("(6) A = g16 + (s m/D) g3", Aexp == G[16] + s * m / D * G[3])
ck.check("(6) A = g23 + e/t", Aexp == G[23] + G[14] / t)
ck.check("(7) B = g19 + (s m/D) g3", Bexp == G[19] + s * m / D * G[3])
ck.check("(7) B = g15 - e/t", Bexp == G[15] - G[14] / t)
ck.check("ratio G/g10 = D k^2 |v-1|^2 / (s |M|^2)",
         G[3] / G[10] == D * k ** 2 * (V - 1) * (Vb - 1) / (s * M1 * M2))
ck.check("ratio g6/e = s / (|v|^2 (k + l s))", G[6] / G[14] == s / (V * Vb * (k + l * s)))
ck.check("ratio g10/g6 = -h (k-l-1) |v|^2 / ((k-1)|v-1|^2)",
         G[10] / G[6] == -h * (k - l - 1) * V * Vb / ((k - 1) * (V - 1) * (Vb - 1)))
ck.check("observation (not in the document): g5 = g3 identically on the chart", G[5] == G[3])

K4, g4 = H.make_field("h,k,s,W")
h4, k4, s4, W4 = [g4[n] for n in ["h", "k", "s", "W"]]
l4 = h4 - 1 + W4                       # m = h - l - 1 = -W
D4 = h4 - k4
al4 = h4 * k4 - h4 * l4 - k4
zeta = s4 / D4 * W4
ck.check("with m = -W: k - l - 1 = -D - W", k4 - l4 - 1 == -D4 - W4)
ck.check("with m = -W: alpha/D = 1 - h(1 + W/D)", al4 / D4 == 1 - h4 * (1 + W4 / D4))
ck.check("interpolation identity: (k + l s)/D = (1 - zeta) k/D - zeta alpha/(D W), zeta = (s/D) W",
         (k4 + l4 * s4) / D4 == (1 - zeta) * k4 / D4 - zeta * al4 / (D4 * W4))

# ---------------------------------------------------------------------------------------------
ck.section("5. Exact evaluation at Gaussian-rational points, directly from the 24 products (EXACT)")


def exact_gaps(hq, kq, lq, sq, vq, mu):
    """All 24 gaps g_j = Re(kappa (F_j - F_0)) for the chart point; exact."""
    tq = (kq - 1) / (hq - 1)
    zq = hq * vq / (1 + (hq - 1) * vq)
    rq = kq * vq / (1 + lq * vq)
    wq = vq + sq * (vq - rq) / rq
    cq = (zq * (vq - rq) - tq * rq * (vq - zq)) / ((vq - rq) - tq * (vq - zq))
    Cq = [wq, vq, rq, cq]
    Fq = [-Cq[p[1]] * (1 - Cq[p[2]]) * (zq - Cq[p[3]]) for p in H.PERMS]
    kap = GQ(0, mu) / ((zq - cq) * (vq - rq))
    labels = [GQ(0), GQ(1), zq] + Cq
    return [(kap * (f - Fq[0])).re for f in Fq], labels


def eval_frac(elem, point):
    """Evaluate an element of Q(h,k,l,s,V,Vb) at GQ values (dict generator-name -> GQ)."""
    names = ["h", "k", "l", "s", "V", "Vb"]

    def ev(poly):
        tot = GQ(0)
        for mon, coef in poly.terms():
            term = GQ(Fraction(int(coef.numerator), int(coef.denominator)))
            for nm, e in zip(names, mon):
                for _ in range(e):
                    term = term * point[nm]
            tot = tot + term
        return tot
    return ev(elem.numer) / ev(elem.denom)


def conj(x):
    return GQ(x.re, -x.im)


exact_points = [
    # (h, k, l, s, Re v, Im v, mu); the first three are supported (all gaps >= 0), the last two are not
    (Fraction(9, 10), Fraction(7, 10), Fraction(-1, 2), Fraction(1, 5), Fraction(9, 10), Fraction(-3, 10), 1),
    (Fraction(3, 2), Fraction(3, 10), Fraction(-7, 5), Fraction(1, 5), Fraction(1, 2), Fraction(-4, 5), 1),
    (Fraction(21, 10), Fraction(9, 10), Fraction(-1, 5), Fraction(7, 10), Fraction(13, 10), Fraction(3, 5), -1),
    (Fraction(-4, 3), Fraction(5, 2), Fraction(2, 7), Fraction(-3, 5), Fraction(1, 3), Fraction(2, 9), 1),
    (Fraction(1, 3), Fraction(-2, 1), Fraction(5, 4), Fraction(7, 3), Fraction(-5, 2), Fraction(-3, 7), -1),
]
ok_all = True
n_supp = 0
ab_pos = True
for (hq, kq, lq, sq, pr, pi_, mu) in exact_points:
    vq = GQ(pr, pi_)
    gq, labels = exact_gaps(hq, kq, lq, sq, vq, mu)
    rho = mu * pi_
    pt = {"h": GQ(hq), "k": GQ(kq), "l": GQ(lq), "s": GQ(sq), "V": vq, "Vb": conj(vq)}
    for j in range(24):
        val = eval_frac(G[j], pt)
        if not (val.im == 0 and val.re * rho == gq[j]):
            ok_all = False
    if not all(gq[j] == 0 for j in T):
        ok_all = False
    if len(set((x.re, x.im) for x in labels)) != 7:
        ok_all = False                              # seven finite labels, all distinct (A_0 = infinity)
    if min(gq) >= 0:
        n_supp += 1
        Aq = -gq[3] + gq[9] + gq[16]
        Bq = -gq[3] + gq[6] + gq[9] - gq[12] + gq[13] + gq[18]
        ab_pos = ab_pos and Aq > 0 and Bq > 0 and hq - lq - 1 >= 0
ck.check("at 5 Gaussian-rational chart points: 24 gaps from the products = rho * (polarized formula), "
         "contacts vanish exactly, labels distinct", ok_all)
ck.check("3 of these points are supported (all 24 gaps >= 0); there A > 0, B > 0 and m >= 0",
         n_supp == 3 and ab_pos)

# ---------------------------------------------------------------------------------------------
ck.section("5b. Observation, not in the document: for eight distinct labels the real-v layer is trivial (EXACT)")
# (i) v real and k != 0: by the identity of Section 3, 1/r = ((k-1)(X+1) + 1/v)/k with X real, so r is
#     real; then z, c, w are real, every W_j is real and all 24 gaps vanish: c(x) = 0.
real_points = [
    (Fraction(3, 2), Fraction(1, 2), Fraction(-1, 3), Fraction(2, 5), Fraction(-4, 7), 1),
    (Fraction(-5, 3), Fraction(7, 2), Fraction(2, 9), Fraction(-3, 4), Fraction(5, 2), -1),
]
ok_real = True
for (hq, kq, lq, sq, vr, mu) in real_points:
    gq, labels = exact_gaps(hq, kq, lq, sq, GQ(vr, 0), mu)
    ok_real = ok_real and all(x == 0 for x in gq) and all(lab.im == 0 for lab in labels) \
        and len(set((x.re, x.im) for x in labels)) == 7
ck.check("at 2 chart points with v real (all labels real, distinct): all 24 gaps are exactly 0", ok_real)
# (ii) v real and k = 0 (t = 1/(1-h)): contact 7 holds for every complex r.  Polarize in r.
K5, g5 = H.make_field("h,s,v,R,Rb")
h5, s5, v5, R5, Rb5 = [g5[n] for n in ["h", "s", "v", "R", "Rb"]]
t5 = 1 / (1 - h5)
ck.check("k = 1 + t(h-1) = 0 for t = 1/(1-h)", 1 + t5 * (h5 - 1) == 0)


def chart_k0(r):
    z = h5 * v5 / (1 + (h5 - 1) * v5)
    w = v5 + s5 * (v5 - r) / r
    c = (z * (v5 - r) - t5 * r * (v5 - z)) / ((v5 - r) - t5 * (v5 - z))
    F = H.products([None, K5(0), K5(1), z], [w, v5, r, c], inf_A=0)
    den = (z - c) * (v5 - r)
    return [(F[j] - F[0]) / den for j in range(24)]


Wk1 = chart_k0(R5)
Wk2 = chart_k0(Rb5)
Gk = [-(a - b) / (R5 - Rb5) for a, b in zip(Wk1, Wk2)]       # g_j / (mu Im r)
ck.check("sub-layer v real, k = 0, r complex: all six contacts hold identically", all(Gk[j] == 0 for j in T))
ck.check("there g6 = -mu Im r s/|r|^2,  g12 = -(v-1) g6,  g18 = mu Im r h s v/|r|^2,  g10 (v-1) = g18",
         Gk[6] == -s5 / (R5 * Rb5) and Gk[12] == -(v5 - 1) * Gk[6]
         and Gk[18] == h5 * s5 * v5 / (R5 * Rb5) and Gk[10] * (v5 - 1) == Gk[18])


def exact_gaps_k0(hq, sq, vr, rq, mu):
    """Gaps on the sub-layer v real, k = 0, directly from the products (exact)."""
    tq = 1 / (1 - hq)
    vq = GQ(vr)
    zq = hq * vq / (1 + (hq - 1) * vq)
    wq = vq + sq * (vq - rq) / rq
    cq = (zq * (vq - rq) - tq * rq * (vq - zq)) / ((vq - rq) - tq * (vq - zq))
    Cq = [wq, vq, rq, cq]
    Fq = [-Cq[p_[1]] * (1 - Cq[p_[2]]) * (zq - Cq[p_[3]]) for p_ in H.PERMS]
    kap = GQ(0, mu) / ((zq - cq) * (vq - rq))
    return [(kap * (f - Fq[0])).re for f in Fq], [GQ(0), GQ(1), zq] + Cq


ok_k0 = True
for (hq, sq, vr, rre, rim, mu) in [
        (Fraction(3, 2), Fraction(-2, 5), Fraction(-4, 7), Fraction(1, 3), Fraction(5, 4), 1),
        (Fraction(-5, 3), Fraction(7, 2), Fraction(5, 2), Fraction(-2, 3), Fraction(-1, 6), -1),
        (Fraction(2, 7), Fraction(-3, 1), Fraction(-1, 9), Fraction(4, 1), Fraction(2, 5), 1)]:
    gq, labels = exact_gaps_k0(hq, sq, vr, GQ(rre, rim), mu)
    rn = rre * rre + rim * rim
    ok_k0 = ok_k0 and all(gq[j] == 0 for j in T) and len(set((u.re, u.im) for u in labels)) == 7 \
        and gq[6] == -mu * rim * sq / rn and gq[12] == -(vr - 1) * gq[6] \
        and gq[18] == mu * rim * hq * sq * vr / rn and gq[10] * (vr - 1) == gq[18] and min(gq) < 0
ck.check("at 3 exact points of this sub-layer (r Gaussian rational, labels distinct): contacts vanish, the four "
         "formulas hold from the products, and some gap is negative", ok_k0)
ck.info("=> if r is nonreal: g6 != 0, so g6 > 0 and g12 >= 0 force v < 1; g18 != 0, so g18 > 0 and then")
ck.info("   g10 = g18/(v-1) < 0.  No supported point.  If r is real all labels are real and c(x) = 0.")
ck.info("=> for eight distinct labels and v real the gap vector is 0: no circle theorem is needed there.")

# ---------------------------------------------------------------------------------------------
ck.section("6. NUMERIC falsification attempts (floating point prefilter, exact confirmation)")
import numpy as np   # noqa: E402

np.seterr(all="ignore")     # degenerate random samples are filtered explicitly below

SEED = 20261125
rng = random.Random(SEED)


def float_gaps(hf, kf, lf, sf, vf, mu):
    tf = (kf - 1) / (hf - 1)
    zf = hf * vf / (1 + (hf - 1) * vf)
    rf = kf * vf / (1 + lf * vf)
    wf = vf + sf * (vf - rf) / rf
    cf = (zf * (vf - rf) - tf * rf * (vf - zf)) / ((vf - rf) - tf * (vf - zf))
    Cf = [wf, vf, rf, cf]
    Ff = [-Cf[p[1]] * (1 - Cf[p[2]]) * (zf - Cf[p[3]]) for p in H.PERMS]
    kap = 1j * mu / ((zf - cf) * (vf - rf))
    return np.array([(kap * (f - Ff[0])).real for f in Ff])


def rnd_real():
    u = rng.random()
    if u < 0.4:
        return rng.choice([-1, 1]) * math.exp(rng.uniform(-3, 3))
    if u < 0.5:
        return rng.choice([0.0, 1.0]) + rng.choice([-1, 1]) * 10 ** rng.uniform(-6, -1)
    return rng.uniform(-3, 3)


def to_q(x):
    return Fraction(x).limit_denominator(10 ** 6)


def exact_tests(hq, kq, lq, sq, pr, pi_, mu, stats):
    """Exact evaluation at a rational point.  Returns False iff a claim of the document fails."""
    try:
        gq, labels = exact_gaps(hq, kq, lq, sq, GQ(pr, pi_), mu)
    except ZeroDivisionError:
        return True
    if not all(gq[j] == 0 for j in T):
        stats["contact_fail"] += 1
        return False
    if min(gq) < 0:
        return True
    stats["supported"] += 1
    A = -gq[3] + gq[9] + gq[16]
    B = -gq[3] + gq[6] + gq[9] - gq[12] + gq[13] + gq[18]
    ok = A >= 0 and B >= 0
    Dq = hq - kq
    mq = hq - lq - 1
    if gq[3] > 0 and gq[14] > 0:
        stats["generic_branch"] += 1
        stats["minA_gen"] = min(stats["minA_gen"], A / max(gq))
        stats["minB_gen"] = min(stats["minB_gen"], B / max(gq))
        ok = ok and sq / Dq > 0 and Dq / (kq + lq * sq) > 0 and Dq / (kq - 1) < 0 \
            and -hq * (kq - lq - 1) / (kq - 1) > 0 and mq >= 0 and gq[6] > 0 and gq[10] > 0
    if gq[3] == 0:
        stats["g3_zero"] += 1
    if gq[14] == 0:
        stats["e_zero"] += 1
        ok = ok and A == gq[23] and B == gq[15]
    # all 14 facet inequalities, i.e. membership in the minor cone
    for key, nrm in facets.items():
        val = sum(a * gq[c] for a, c in zip(nrm, cols))
        if val < 0:
            ok = False
    if not ok:
        stats["violations"] += 1
    stats["minA"] = min(stats["minA"], A / max(gq))
    stats["minB"] = min(stats["minB"], B / max(gq))
    return ok


stats = dict(supported=0, generic_branch=0, g3_zero=0, e_zero=0, violations=0, contact_fail=0,
             minA=Fraction(10 ** 9), minB=Fraction(10 ** 9),
             minA_gen=Fraction(10 ** 9), minB_gen=Fraction(10 ** 9))
N1 = 150000
cand = 0
for it in range(N1):
    hf, kf, lf, sf = rnd_real(), rnd_real(), rnd_real(), rnd_real()
    vf = complex(rng.uniform(-3, 3), rng.choice([-1, 1]) * 10 ** rng.uniform(-4, 0.5))
    mu = rng.choice([-1, 1])
    try:
        gf = float_gaps(hf, kf, lf, sf, vf, mu)
    except ZeroDivisionError:
        continue
    sc = np.max(np.abs(gf))
    if not np.isfinite(sc) or sc == 0:
        continue
    if gf.min() >= -1e-9 * sc:
        cand += 1
        exact_tests(to_q(hf), to_q(kf), to_q(lf), to_q(sf), to_q(vf.real), to_q(vf.imag), mu, stats)
ck.info("N1a: %d unrestricted chart samples (seed %d), %d float-supported candidates re-evaluated exactly"
        % (N1, SEED, cand))
# N1b: targeted sampling in the sign pattern in which supported points were found
N1b = 60000
cand_b = 0
for it in range(N1b):
    kf = rng.uniform(0.001, 0.999)
    hf = kf + math.exp(rng.uniform(-5, 2.5))
    lf = hf - 1 - rng.choice([1, -1]) * math.exp(rng.uniform(-6, 2))
    sf = math.exp(rng.uniform(-5, 2.5))
    vf = complex(rng.uniform(-3, 3), rng.choice([-1, 1]) * 10 ** rng.uniform(-4, 0.5))
    mu = rng.choice([-1, 1])
    try:
        gf = float_gaps(hf, kf, lf, sf, vf, mu)
    except ZeroDivisionError:
        continue
    sc = np.max(np.abs(gf))
    if not np.isfinite(sc) or sc == 0:
        continue
    if gf.min() >= -1e-9 * sc:
        cand_b += 1
        exact_tests(to_q(hf), to_q(kf), to_q(lf), to_q(sf), to_q(vf.real), to_q(vf.imag), mu, stats)
ck.info("N1b: %d targeted samples (0<k<1, h>k, s>0), %d float-supported candidates re-evaluated exactly"
        % (N1b, cand_b))
# N1c: the two degenerate branches of the proof, placed exactly: e = 0 (k + l s = 0) and g3 = 0 (n = 0)
N1c = 20000
for it in range(N1c):
    kq = Fraction(rng.randint(1, 999), 1000)
    hq = kq + Fraction(rng.randint(1, 4000), 1000)
    mu = rng.choice([-1, 1])
    pr = Fraction(rng.randint(-3000, 3000), 1000)
    pi_ = Fraction(rng.choice([-1, 1]) * rng.randint(1, 3000), 1000)
    if it % 2 == 0:                                  # e = 0
        lq = -Fraction(rng.randint(1, 4000), 1000)
        sq = -kq / lq
    else:                                            # n = k - l - 1 = 0, hence g3 = 0
        lq = kq - 1
        sq = Fraction(rng.randint(1, 4000), 1000)
    exact_tests(hq, kq, lq, sq, pr, pi_, mu, stats)
ck.info("N1c: %d exact rational samples on the branches e = 0 and g3 = 0" % N1c)
ck.check("N1 (EXACT at sampled points): contacts vanish; every supported sample has A >= 0, B >= 0, "
         "all 14 facet inequalities, and on the generic branch the sign claims (8) and m >= 0",
         stats["violations"] == 0 and stats["contact_fail"] == 0 and stats["supported"] > 500,
         "supported %d (generic branch %d, g3 = 0: %d, e = 0: %d); min A/max gap = %.3e, min B/max gap = %.3e "
         "overall; %.3e, %.3e on the generic branch"
         % (stats["supported"], stats["generic_branch"], stats["g3_zero"], stats["e_zero"],
            float(stats["minA"]), float(stats["minB"]), float(stats["minA_gen"]), float(stats["minB_gen"])))

# N2: the pure sign implication of p.3 on random real parameters (exact rationals)
N2 = 200000
neg_m = 0
four = 0
bad = 0
for it in range(N2):
    hq = Fraction(rng.randint(-4000, 4000), 1000)
    kq = Fraction(rng.randint(-4000, 4000), 1000)
    lq = Fraction(rng.randint(-4000, 4000), 1000)
    sq = Fraction(rng.randint(-4000, 4000), 1000)
    Dq = hq - kq
    if Dq == 0 or kq == 1 or kq + lq * sq == 0 or sq == 0:
        continue
    mq = hq - lq - 1
    if mq >= 0:
        continue
    neg_m += 1
    if sq / Dq > 0 and Dq / (kq + lq * sq) > 0 and Dq / (kq - 1) < 0 and -hq * (kq - lq - 1) / (kq - 1) > 0:
        four += 1
        if 1 + sq * mq / Dq >= 0:
            bad += 1
ck.check("N2 (EXACT at sampled points): for m < 0 the four inequalities (8) and g9 >= 0 never hold together",
         bad == 0 and four > 200,
         "%d random rational (h,k,l,s); m < 0 at %d; (8) holds at %d of them; with g9 >= 0 as well: %d"
         % (N2, neg_m, four, bad))

# N3: points satisfying the six contact equalities found WITHOUT the chart (least squares on the
# spectra from random starts): the chart parameters must come out real, and the chart formulas
# must reproduce all 24 gaps.  If such a point happens to be supported, A and B must be >= 0.
nprng = np.random.default_rng(SEED)
N3 = 400
found = 0
used = 0
worst_real = 0.0
worst_rebuild = 0.0
n_supp3 = 0
worst_AB = 0.0
for it in range(N3):
    res = H.find_contact_point(T, nprng, tries=3)
    if res is None:
        continue
    found += 1
    a, cc, kap, gf = res
    if not H.labels_well_separated(a, cc):
        continue
    mob = H.moebius_to(a[0], a[1], a[2])
    zc = mob(a[3])
    wc, vc, rc, c4 = [mob(x) for x in cc]
    if max(abs(x) for x in (zc, wc, vc, rc, c4)) > 1e4:
        continue
    used += 1
    hc = zc * (vc - 1) / (vc * (zc - 1))
    tc = (c4 - zc) * (vc - rc) / ((c4 - rc) * (vc - zc))
    kc = 1 + tc * (hc - 1)
    lc = kc / rc - 1 / vc
    sc_ = rc * (wc - vc) / (vc - rc)
    worst_real = max(worst_real, max(abs(x.imag) / (1 + abs(x)) for x in (hc, tc, lc, sc_)))
    # rebuild the gaps from the chart with the real parts of the parameters; the ratio of the two
    # gap vectors is one real constant (the chart normalizes kappa (F_2 - F_0) = i mu)
    g_chart = float_gaps(hc.real, kc.real, lc.real, sc_.real, vc, 1)
    i0 = int(np.argmax(np.abs(gf)))
    scale = gf[i0] / g_chart[i0]
    worst_rebuild = max(worst_rebuild, np.max(np.abs(gf - scale * g_chart)) / np.max(np.abs(gf)))
    if gf.min() > -1e-9 and abs(vc.imag) > 1e-6:
        n_supp3 += 1
        A = -gf[3] + gf[9] + gf[16]
        B = -gf[3] + gf[6] + gf[9] - gf[12] + gf[13] + gf[18]
        worst_AB = max(worst_AB, -min(A, B) / np.max(np.abs(gf)))
ck.check("N3 (NUMERIC): at chart-free solutions of the six contact equalities with well separated labels, "
         "h, t, l, s are real and the chart reproduces all 24 gaps",
         used >= 50 and worst_real < 1e-6 and worst_rebuild < 1e-6 and worst_AB < 1e-7,
         "%d searches, %d converged, %d with separated labels; max |Im(h,t,l,s)|/(1+|.|) = %.1e; "
         "max relative gap mismatch = %.1e; supported among them: %d (worst -min(A,B) = %.1e)"
         % (N3, found, used, worst_real, worst_rebuild, n_supp3, worst_AB))

# N4: the dual criterion against nonnegative least squares on random points of the 12-ray cone
N4 = 3000
mism = 0
nin = 0
Mall = np.array([H.F_VALUES[j] for j in gens], float)
for it in range(N4):
    lam = nprng.exponential(size=len(gens)) * (nprng.random(len(gens)) < 0.6)
    lam[gens.index(22)] *= nprng.choice([0.0, 0.3, 1.0, 3.0])
    lam[gens.index(26)] *= nprng.choice([0.0, 0.3, 1.0, 3.0])
    gv = lam @ Mall
    if gv.max() == 0:
        continue
    gv = gv / gv.max()
    A = -gv[3] + gv[9] + gv[16]
    B = -gv[3] + gv[6] + gv[9] - gv[12] + gv[13] + gv[18]
    rn = H.minor_cone_nnls_residual(gv, minor)
    inside_dual = (A >= -1e-12 and B >= -1e-12)
    inside_nnls = rn < 1e-9
    if abs(A) < 1e-7 or abs(B) < 1e-7:
        continue
    nin += inside_dual
    if inside_dual != inside_nnls:
        mism += 1
ck.check("N4 (NUMERIC): on %d random points of the contact cone, 'A >= 0 and B >= 0' agrees with "
         "a nonnegative-least-squares test of membership in the minor cone" % N4, mism == 0,
         "%d inside, mismatches %d" % (nin, mism))

# N5: chart-free probe of observation O1 (Section 5b): solutions of the six contact equalities with a
# REAL spectrum A (four real a_i) and complex C, found by least squares on (a, C, phase of kappa).
# The gaps are normalized to mean square 1 off T_7, which excludes the all-real solutions (zero gaps).
from scipy.optimize import least_squares   # noqa: E402

T_sorted = sorted(T)
OFF = [j for j in range(24) if j not in T]


def unpack_realA(xv):
    return xv[0:4].astype(complex), xv[4:8] + 1j * xv[8:12], np.exp(1j * xv[12])


def gaps_realA(xv):
    a_, c_, kap_ = unpack_realA(xv)
    Fv = H.float_products(a_, c_)
    return (kap_ * (Fv - Fv[T_sorted[0]])).real


def res_realA(xv):
    gv = gaps_realA(xv)
    return np.array(list(gv[T_sorted[1:]]) + [np.mean(gv[OFF] ** 2) - 1.0])


N5 = 800
usable = n_k0 = n_other = n_supp5 = 0
best5 = -1.0
for it in range(N5):
    try:
        sol = least_squares(res_realA, nprng.normal(size=13), xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=2000)
    except Exception:
        continue
    if np.max(np.abs(res_realA(sol.x))) > 1e-11:
        continue
    a_, c_, kap_ = unpack_realA(sol.x)
    if not H.labels_well_separated(a_, c_):
        continue
    usable += 1
    gv = gaps_realA(sol.x)
    mob = H.moebius_to(a_[0], a_[1], a_[2])
    zc = mob(a_[3])
    wc, vc, rc, c4 = [mob(u) for u in c_]
    hc = zc * (vc - 1) / (vc * (zc - 1))
    tc = (c4 - zc) * (vc - rc) / ((c4 - rc) * (vc - zc))
    kc = 1 + tc * (hc - 1)
    if abs(kc) < 1e-6 and abs(vc.imag) < 1e-6 * (1 + abs(vc)):
        n_k0 += 1
    else:
        n_other += 1
    ratio5 = gv.min() / np.max(np.abs(gv))
    best5 = max(best5, ratio5)
    if ratio5 > -1e-9:
        n_supp5 += 1
ck.check("N5 (NUMERIC): chart-free solutions of the contact equalities with real spectrum A, nonzero gaps and "
         "well separated labels all have v real and k = 0, and none has all gaps >= 0",
         usable >= 5 and n_other == 0 and n_supp5 == 0,
         "%d searches, %d usable, %d with k = 0, %d other, %d supported; largest (min gap)/(max |gap|) = %.2e"
         % (N5, usable, n_k0, n_other, n_supp5, best5))

ck.finish()
