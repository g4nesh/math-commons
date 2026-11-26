"""Independent replay of Case T_19 = {0,1,4,6,7,8} (Zheng, MOC n = 4; group II).

Document replayed: root_moc4_T18_uniform_20261022.pdf (3 pages; historical label "T18").
Nothing from the supplement is imported or executed.  If MOC4_SUPP is set, the PDF is hashed
(sha256) to pin the document; no other use.  All objects are re-derived in replay_T7_common.py.

Run:  MOC4_SUPP=/path/to/supplement python replay_T19.py     (exit code 0 iff every check passes)
Packages: sympy (exact part), numpy + scipy (NUMERIC part only).

Sections 1-6 are EXACT.  Section 7 is NUMERIC (falsification attempts; seeds printed).

Method for the kernel identities.  In the chart kappa F_0 = i mu and R_j = F_j/F_0 is a rational
function of v with real coefficients (x, h, k, l real), so g_j = -mu Im R_j(v) and
        g_j / (mu Im v) = -(R_j(V) - R_j(Vb)) / (V - Vb)     at V = v, Vb = conj(v),
computed in Q(x,h,k,l,V,Vb) with V, Vb independent.  |p(v)|^2 = p(V) p(Vb) for real polynomials p.
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

ck = H.Checks("replay_T19: Case T_19 = {0,1,4,6,7,8}, document root_moc4_T18_uniform_20261022.pdf")
T = H.PAPER_T[19]
DOC_T = (0, 1, 4, 6, 7, 8)                       # the document's "T18"
DOC_GENS = [0, 1, 2, 3, 4, 5, 6, 21, 22, 23, 33, 34]
DOC_COORD = [2, 5, 9, 10, 14, 15, 16, 17, 18, 19, 22, 23]
A_COEF = {13: 1, 2: -1}                          # A = g13 - g2
B_COEF = {2: 1, 3: -1, 9: 1, 12: -1, 13: 1, 18: 1}

# ---------------------------------------------------------------------------------------------
ck.section("1. Objects (EXACT)")
H.document_hash_check(ck, "root_moc4_T18_uniform_20261022.pdf")
H.facet_table_sanity(ck)
H.crosscheck_with_audit_lib(ck)
ck.check("document's contact set 'T18' = paper's T_19 (identical set, same indexing; group element = identity)",
         tuple(sorted(DOC_T)) == tuple(sorted(T)))

# ---------------------------------------------------------------------------------------------
ck.section("2. Contact cone and complete minor dual, doc p.1 (EXACT)")
gens, minor, nonminor, facets, cols = H.minor_cone_report(
    ck, T, DOC_GENS, [22, 34], DOC_COORD, {"A": A_COEF, "B": B_COEF}, 14)

# ---------------------------------------------------------------------------------------------
ck.section("3. Rectangle argument and exhaustiveness of the chart, doc p.1-2 (EXACT)")
K, A, C, F = H.rectangle_checks(ck)
ck.check("F4/F1 = (A1-C3)(A2-C1)/((A1-C1)(A2-C3))",
         F[4] / F[1] == (A[1] - C[3]) * (A[2] - C[1]) / ((A[1] - C[1]) * (A[2] - C[3])))
ck.check("F4 - F1 = (A0-C0)(A3-C2)(A1-A2)(C3-C1)   (nonzero for distinct labels, so h != 1)",
         F[4] - F[1] == (A[0] - C[0]) * (A[3] - C[2]) * (A[1] - A[2]) * (C[3] - C[1]))

Kn, gn = H.make_field("x,A2,A3,v,w")
x, A2, A3, v, w = [gn[n] for n in ["x", "A2", "A3", "v", "w"]]
Fn = H.products([x, Kn(1), A2, A3], [None, Kn(0), v, w], inf_C=0)
ck.check("normalized labels A = (x,1,A2,A3), C = (inf,0,v,w): F6/F0 = x", Fn[6] / Fn[0] == x)
ck.check("F8/F6 = (1 - v)/(A2 - v)", Fn[8] / Fn[6] == (1 - v) / (A2 - v))
ck.check("F4/F1 = (1 - w) A2/(A2 - w)", Fn[4] / Fn[1] == (1 - w) * A2 / (A2 - w))
ck.check("F1/F0 = (A2 - w)(A3 - v)/((A2 - v)(A3 - w))",
         Fn[1] / Fn[0] == (A2 - w) * (A3 - v) / ((A2 - v) * (A3 - w)))

Kc, gc = H.make_field("x,h,k,l,v")
x, h, k, l, v = [gc[n] for n in ["x", "h", "k", "l", "v"]]
a2 = 1 + k * (v - 1)
w_doc = (1 - h) * a2 / (a2 - h)
a3_doc = a2 * (l * (k - 1) * (1 - h) - k * v) / (l * (k - 1) * (a2 - h) - k * a2)
a3_solved = (l * w_doc * (a2 - v) - v * (a2 - w_doc)) / (l * (a2 - v) - (a2 - w_doc))
ck.check("A2 = 1 + k(v-1) is the solution of F8/F6 = -1/(k-1); w = (1-h)a2/(a2-h) the solution of F4/F1 = h",
         (1 - v) / (a2 - v) == -1 / (k - 1) and (1 - w_doc) * a2 / (a2 - w_doc) == h)
ck.check("the document's A3 equals the solution of F1/F0 = l (linear equation in A3)", a3_doc == a3_solved)
Fc = H.products([x, Kc(1), a2, a3_doc], [None, Kc(0), v, w_doc], inf_C=0)
Rc = [f / Fc[0] for f in Fc]
ck.check("(R0,R1,R4,R6,R7,R8) = (1, l, hl, x, lx, -x/(k-1))",
         [Rc[j] for j in (0, 1, 4, 6, 7, 8)] == [Kc(1), l, h * l, x, l * x, -x / (k - 1)])

# ---------------------------------------------------------------------------------------------
ck.section("4. Kernel identities (3) and Poisson identities (5), doc p.2 (EXACT)")
K3, g3 = H.make_field("x,h,k,l,V,Vb")
NAMES = ["x", "h", "k", "l", "V", "Vb"]
x, h, k, l, V, Vb = [g3[n] for n in NAMES]


def chart(v):
    a2 = 1 + k * (v - 1)
    w = (1 - h) * a2 / (a2 - h)
    a3 = a2 * (l * (k - 1) * (1 - h) - k * v) / (l * (k - 1) * (a2 - h) - k * a2)
    Fv = H.products([x, K3(1), a2, a3], [None, K3(0), v, w], inf_C=0)
    return [f / Fv[0] for f in Fv]


R1 = chart(V)
R2 = chart(Vb)
G = [-(a - b) / (V - Vb) for a, b in zip(R1, R2)]         # G_j = g_j / (mu y), y = Im v
ck.check("the chart satisfies all six contacts: Im R_j = 0 identically for j in T_19",
         all(G[j] == 0 for j in T))
ck.check("no other gap vanishes identically on the chart",
         all(G[j] != 0 for j in range(24) if j not in T))
D1 = k * V + (k - 1) * (h - 1)
D2 = k * Vb + (k - 1) * (h - 1)
ck.check("(3) g2 = mu k y/(k-1)", G[2] == k / (k - 1))
ck.check("(3) g5/g2 = h (k-1)^2 (h-1)(l-1)/|D|^2", G[5] / G[2] == h * (k - 1) ** 2 * (h - 1) * (l - 1) / (D1 * D2))
ck.check("(3) g9/g5 = -k x/((h-1)(k-1))", G[9] / G[5] == -k * x / ((h - 1) * (k - 1)))
ck.check("(3) g10/g2 = h l x (k-1)/|a2|^2, a2 = kv - k + 1",
         G[10] / G[2] == h * l * x * (k - 1) / ((k * V - k + 1) * (k * Vb - k + 1)))
ck.check("(3) g14/g2 = (x-1)/(k |v-1|^2)", G[14] / G[2] == (x - 1) / (k * (V - 1) * (Vb - 1)))
ck.check("(3) g18/g14 = -l (h-1)(k-1)/k", G[18] / G[14] == -l * (h - 1) * (k - 1) / k)
e = 1 + h * (k - 1)
P = e + k * (x - 1)
L = l * e + k * (1 - l)
Aexp = G[13] - G[2]
Bexp = G[2] - G[3] + G[9] - G[12] + G[13] + G[18]
ck.check("(5) A/g2 = h(h-1)(k-1)^2(l-1)P/(e|D|^2) + (x-1)L/(k e |v-1|^2)",
         Aexp / G[2] == h * (h - 1) * (k - 1) ** 2 * (l - 1) * P / (e * D1 * D2)
         + (x - 1) * L / (k * e * (V - 1) * (Vb - 1)))
ck.check("(5) B/g2 = -h k (k-1)(l-1)P/(e|D|^2) - (h-1)(k-1)(x-1)L/(k^2 e |v-1|^2)",
         Bexp / G[2] == -h * k * (k - 1) * (l - 1) * P / (e * D1 * D2)
         - (h - 1) * (k - 1) * (x - 1) * L / (k ** 2 * e * (V - 1) * (Vb - 1)))

# ---------------------------------------------------------------------------------------------
ck.section("5. Sign analysis, doc p.2-3 (EXACT: complete enumeration of sign cells; certificates)")
# The five ratios of (3) are positive.  After removing the squared moduli, each ratio is a
# rational function of (x,h,k,l) alone; its sign is constant on each of the 81 cells of
# (x,h,k,l) relative to 0 and 1 and is certified cell by cell.
Kp, gp = H.make_field("a,b,c,d")
pa, pb, pc, pd = [gp[n] for n in ["a", "b", "c", "d"]]
prefactors = [G[5] / G[2] * D1 * D2,                                  # g5/g2  * |D|^2
              G[9] / G[5],                                            # g9/g5
              G[10] / G[2] * (k * V - k + 1) * (k * Vb - k + 1),      # g10/g2 * |a2|^2
              G[14] / G[2] * (V - 1) * (Vb - 1),                      # g14/g2 * |v-1|^2
              G[18] / G[14]]                                          # g18/g14
table = H.cell_signs(prefactors, NAMES, ["x", "h", "k", "l"], Kp, [pd, pa, pb, pc])
ck.check("the five kernel ratios, stripped of squared moduli, are free of v and have a certified sign on "
         "each of the 81 cells", len(table) == 81 and all(0 not in sg for _, sg in table))
good = [lab for lab, sg in table if min(sg) > 0]
regI = {"x": ">1", "h": "(0,1)", "k": ">1", "l": "(0,1)"}
regII = {"x": "(0,1)", "h": "<0", "k": "<0", "l": ">1"}
ck.check("of the 81 sign cells of (x,h,k,l) exactly two are compatible with the five positive ratios: "
         "I (k>1, 0<h,l<1, x>1) and II (h,k<0, l>1, 0<x<1)",
         len(good) == 2 and regI in good and regII in good, str(good))


def region_values(reg):
    if reg == "I":      # k > 1, 0 < h, l < 1, x > 1
        return dict(h=1 / (1 + pa), k=1 + pb, l=1 / (1 + pc), x=1 + pd)
    return dict(h=-pa, k=-pb, l=1 + pc, x=1 / (1 + pd))         # II: h, k < 0, l > 1, 0 < x < 1


for reg in ("I", "II"):
    val = region_values(reg)
    hh, kk, ll, xx = val["h"], val["k"], val["l"], val["x"]
    ee = 1 + hh * (kk - 1)
    PP = ee + kk * (xx - 1)
    LL = ll * ee + kk * (1 - ll)
    ck.check("region %s: e > 1, P > e, L > 0 (positive-coefficient certificate in a,b,c,d > 0)" % reg,
             [H.certified_sign(q) for q in (ee - 1, PP - ee, LL)] == [1, 1, 1])
    terms = [hh * (hh - 1) * (kk - 1) ** 2 * (ll - 1) * PP / ee,
             (xx - 1) * LL / (kk * ee),
             -hh * kk * (kk - 1) * (ll - 1) * PP / ee,
             -(hh - 1) * (kk - 1) * (xx - 1) * LL / (kk ** 2 * ee)]
    ck.check("region %s: the four coefficients of the positive kernels in (5) are positive (certificate)" % reg,
             [H.certified_sign(q) for q in terms] == [1, 1, 1, 1])
    ck.info("region %s: certified sign of k/(k-1) is %+d, so g2 > 0 fixes the sign of mu Im v"
            % (reg, H.certified_sign(kk / (kk - 1))))

# ---------------------------------------------------------------------------------------------
ck.section("5b. Observation, not in the document: the supported eight-distinct nonreal layer is EMPTY (EXACT)")
# g15, g17, g19, g22 share one numerator N, which is a real quadratic in v:
#     N = alpha |v|^2 + beta (v + conj v) + gamma,   alpha, beta, gamma in Q[x,h,k,l].
Nfrac = G[15] / G[2] * k * (V - 1) * (Vb - 1) * D1 * D2
ck.check("N := (g15/g2) k |v-1|^2 |D|^2 is a polynomial", Nfrac.denom == 1)
Npoly = Nfrac.numer
ck.check("g17 = h g15   and   g22 = -g19/(k-1)   identically on the chart",
         G[17] == h * G[15] and G[22] == -G[19] / (k - 1))
ck.check("g19/g2 = -(h-1)(k-1) N/(k^2 |v-1|^2 |D|^2),   g22/g2 = (h-1) N/(k^2 |v-1|^2 |D|^2)",
         G[19] / G[2] == -(h - 1) * (k - 1) * Nfrac / (k ** 2 * (V - 1) * (Vb - 1) * D1 * D2)
         and G[22] / G[2] == (h - 1) * Nfrac / (k ** 2 * (V - 1) * (Vb - 1) * D1 * D2))
ring5 = Npoly.ring
coef = {}
for mon, cf in Npoly.terms():
    key = (mon[4], mon[5])                      # exponents of V, Vb
    rest = ring5({(mon[0], mon[1], mon[2], mon[3], 0, 0): cf})
    coef[key] = coef.get(key, ring5.zero) + rest
ck.check("N has only the monomials V Vb, V, Vb, 1 and is symmetric in V, Vb",
         set(coef) <= {(1, 1), (1, 0), (0, 1), (0, 0)} and coef[(1, 0)] == coef[(0, 1)])
alpha_, beta_, gamma_ = K3(coef[(1, 1)]), K3(coef[(1, 0)]), K3(coef[(0, 0)])
ck.check("alpha gamma - beta^2 = -h k^2 (k-1)(l-1)(x-1) P L   (P, L as in the document)",
         alpha_ * gamma_ - beta_ ** 2 == -h * k ** 2 * (k - 1) * (l - 1) * (x - 1) * P * L)
for reg, sgn_alpha in (("I", 1), ("II", -1)):
    val = region_values(reg)

    def sub(elem, val=val):
        """Substitute the region parametrization into an element of Q(x,h,k,l) (no V, Vb)."""
        return H.substitute(elem, NAMES, val, Kp)
    ck.check("region %s: alpha is %s and alpha gamma - beta^2 > 0 (certificates), hence N = "
             "alpha |v + beta/alpha|^2 + (alpha gamma - beta^2)/alpha has the sign of alpha for EVERY complex v"
             % (reg, "positive" if sgn_alpha > 0 else "negative"),
             H.certified_sign(sub(alpha_)) == sgn_alpha
             and H.certified_sign(sub(alpha_ * gamma_ - beta_ ** 2)) == 1)
    if reg == "I":
        ck.check("region I: (h-1)/k^2 < 0, so g22/g2 < 0: g22 < 0 whenever g2 > 0",
                 H.certified_sign(sub((h - 1) / k ** 2)) == -1)
    else:
        ck.check("region II: h/k > 0 and N < 0, so g17/g2 = h N/(k |v-1|^2 |D|^2) < 0: g17 < 0 whenever g2 > 0",
                 H.certified_sign(sub(h / k)) == 1
                 and G[17] / G[2] == h * Nfrac / (k * (V - 1) * (Vb - 1) * D1 * D2))
ck.info("=> with Section 5: no configuration with eight distinct labels and nonreal v has all gaps >= 0;")
ck.info("   every realizable point of the contact cone of T_19 has a repeated or shared label or v real.")

# ---------------------------------------------------------------------------------------------
ck.section("6. Exact evaluation at Gaussian-rational points, directly from the 24 products (EXACT)")


def labels_q(xq, hq, kq, lq, vq):
    a2q = 1 + kq * (vq - 1)
    wq = (1 - hq) * a2q / (a2q - hq)
    a3q = a2q * (lq * (kq - 1) * (1 - hq) - kq * vq) / (lq * (kq - 1) * (a2q - hq) - kq * a2q)
    return [GQ(xq), GQ(1), a2q, a3q], [None, GQ(0), vq, wq]


exact_points = [
    # (x, h, k, l, Re v, Im v, mu): two points in region I, two in region II, two outside
    (Fraction(3, 2), Fraction(1, 2), Fraction(7, 3), Fraction(2, 5), Fraction(1, 4), Fraction(3, 5), 1),
    (Fraction(5, 1), Fraction(1, 5), Fraction(9, 4), Fraction(1, 7), Fraction(-2, 1), Fraction(1, 3), 1),
    (Fraction(1, 2), Fraction(-3, 1), Fraction(-2, 1), Fraction(5, 2), Fraction(1, 3), Fraction(4, 7), 1),
    (Fraction(1, 5), Fraction(-3, 2), Fraction(-1, 4), Fraction(7, 2), Fraction(-1, 1), Fraction(2, 3), 1),
    (Fraction(-4, 3), Fraction(5, 2), Fraction(2, 7), Fraction(-3, 5), Fraction(1, 3), Fraction(2, 9), 1),
    (Fraction(1, 3), Fraction(-2, 1), Fraction(5, 4), Fraction(7, 3), Fraction(-5, 2), Fraction(-3, 7), -1),
]
ok_all = True
region_ok = True
for (xq, hq, kq, lq, pr, pi_, mu) in exact_points:
    vq = GQ(pr, pi_)
    Aq, Cq = labels_q(xq, hq, kq, lq, vq)
    gq = H.exact_ratio_gaps(Aq, Cq, mu)
    pt = {"x": GQ(xq), "h": GQ(hq), "k": GQ(kq), "l": GQ(lq), "V": vq, "Vb": H.gq_conj(vq)}
    for j in range(24):
        val = H.eval_frac(G[j], NAMES, pt)
        if not (val.im == 0 and val.re * mu * pi_ == gq[j]):
            ok_all = False
    if not all(gq[j] == 0 for j in T):
        ok_all = False
    labs = Aq + Cq[1:]
    if len(set((u.re, u.im) for u in labs)) != 7:
        ok_all = False
    if (kq > 1 and 0 < hq < 1) or (hq < 0 and kq < 0):
        Aval = gq[13] - gq[2]
        Bval = gq[2] - gq[3] + gq[9] - gq[12] + gq[13] + gq[18]
        bad_gap = gq[22] if kq > 1 else gq[17]
        region_ok = region_ok and Aval / gq[2] > 0 and Bval / gq[2] > 0 and bad_gap / gq[2] < 0
ck.check("at 6 Gaussian-rational chart points: 24 gaps from the products = mu y * (polarized formula), "
         "contacts vanish exactly, labels distinct", ok_all)
ck.check("at the 4 points in regions I, II: A/g2 > 0, B/g2 > 0 and g22/g2 < 0 (I), g17/g2 < 0 (II)", region_ok)

# ---------------------------------------------------------------------------------------------
ck.section("6b. Observation, not in the document: for eight distinct labels the real-v layer is trivial (EXACT)")
# x, h, k, l are real for every configuration with eight distinct labels (Section 3).  The chart
# expresses A2, A3, w as rational functions of (x,h,k,l,v) with rational coefficients, so v real
# makes all eight labels real; then every ratio R_j is real and g_j = -mu Im R_j = 0 for all j.
real_points = [
    (Fraction(3, 2), Fraction(1, 2), Fraction(7, 3), Fraction(2, 5), Fraction(-4, 7), 1),
    (Fraction(1, 2), Fraction(-3, 1), Fraction(-2, 1), Fraction(5, 2), Fraction(7, 3), -1),
]
ok_real = True
for (xq, hq, kq, lq, vr, mu) in real_points:
    Aq, Cq = labels_q(xq, hq, kq, lq, GQ(vr, 0))
    gq = H.exact_ratio_gaps(Aq, Cq, mu)
    labs = Aq + Cq[1:]
    ok_real = ok_real and all(u == 0 for u in gq) and all(u.im == 0 for u in labs) \
        and len(set((u.re, u.im) for u in labs)) == 7
ck.check("at 2 chart points with v real (all labels real, distinct): all 24 gaps are exactly 0", ok_real)
ck.info("=> for eight distinct labels and v real the gap vector is 0: no circle theorem is needed there.")

# ---------------------------------------------------------------------------------------------
ck.section("7. NUMERIC falsification attempts (floating point prefilter, exact confirmation)")
import numpy as np   # noqa: E402

np.seterr(all="ignore")     # degenerate random samples are filtered explicitly below

SEED = 20261125
rng = random.Random(SEED)


def labels_f(xf, hf, kf, lf, vf):
    a2f = 1 + kf * (vf - 1)
    wf = (1 - hf) * a2f / (a2f - hf)
    a3f = a2f * (lf * (kf - 1) * (1 - hf) - kf * vf) / (lf * (kf - 1) * (a2f - hf) - kf * a2f)
    return [xf, 1.0, a2f, a3f], [None, 0.0, vf, wf]


def rnd_real():
    u = rng.random()
    if u < 0.4:
        return rng.choice([-1, 1]) * math.exp(rng.uniform(-3, 3))
    if u < 0.5:
        return rng.choice([0.0, 1.0]) + rng.choice([-1, 1]) * 10 ** rng.uniform(-6, -1)
    return rng.uniform(-3, 3)


def rnd_unit():
    """Random point of (0,1), with mass near both ends."""
    u = rng.random()
    if u < 0.2:
        return 10 ** rng.uniform(-6, -1)
    if u < 0.4:
        return 1 - 10 ** rng.uniform(-6, -1)
    return rng.uniform(0.001, 0.999)


def rnd_above_one():
    u = rng.random()
    if u < 0.3:
        return 1 + 10 ** rng.uniform(-6, -1)
    return 1 + math.exp(rng.uniform(-3, 4))


def to_q(val):
    return Fraction(val).limit_denominator(10 ** 7)


def in_region(xq, hq, kq, lq):
    if kq > 1 and 0 < hq < 1 and 0 < lq < 1 and xq > 1:
        return "I"
    if hq < 0 and kq < 0 and lq > 1 and 0 < xq < 1:
        return "II"
    return None


stats = dict(supported=0, I=0, II=0, violations=0, contact_fail=0, minA=None, minB=None)


def exact_tests(xq, hq, kq, lq, pr, pi_, mu):
    """Exact evaluation at a rational chart point."""
    try:
        Aq, Cq = labels_q(xq, hq, kq, lq, GQ(pr, pi_))
        gq = H.exact_ratio_gaps(Aq, Cq, mu)
    except ZeroDivisionError:
        return
    if not all(gq[j] == 0 for j in T):
        stats["contact_fail"] += 1
        return
    reg = in_region(xq, hq, kq, lq)
    if min(gq) >= 0:
        stats["supported"] += 1          # Section 5b says this never happens
        if reg is None:
            stats["violations"] += 1
    if reg is None:
        return
    # the document's intermediate claims on its stated regions (they do not need support):
    stats[reg] += 1
    g2 = gq[2]
    Aval = gq[13] - gq[2]
    Bval = gq[2] - gq[3] + gq[9] - gq[12] + gq[13] + gq[18]
    ok = Aval / g2 > 0 and Bval / g2 > 0
    ok = ok and all(gq[j] / g2 > 0 for j in (5, 9, 10, 14, 18))
    # Section 5b: the gap that is negative
    ok = ok and (gq[22] / g2 < 0 if reg == "I" else gq[17] / g2 < 0)
    if not ok:
        stats["violations"] += 1
    ra, rb = Aval / g2, Bval / g2
    stats["minA"] = ra if stats["minA"] is None else min(stats["minA"], ra)
    stats["minB"] = rb if stats["minB"] is None else min(stats["minB"], rb)


def run_batch(sampler, n, prefilter):
    cand = 0
    for _ in range(n):
        xf, hf, kf, lf = sampler()
        vf = complex(rng.uniform(-3, 3), rng.choice([-1, 1]) * 10 ** rng.uniform(-4, 0.5))
        mu = rng.choice([-1, 1])
        if prefilter:
            try:
                Af, Cf = labels_f(xf, hf, kf, lf, vf)
                gf = H.float_ratio_gaps(Af, Cf, mu)
            except ZeroDivisionError:
                continue
            sc = np.max(np.abs(gf))
            if not np.isfinite(sc) or sc == 0 or gf.min() < -1e-9 * sc:
                continue
        cand += 1
        exact_tests(to_q(xf), to_q(hf), to_q(kf), to_q(lf), to_q(vf.real), to_q(vf.imag), mu)
    return cand


N1 = 100000
c1 = run_batch(lambda: (rnd_real(), rnd_real(), rnd_real(), rnd_real()), N1, True)
ck.info("N1a: %d unrestricted chart samples (seed %d): %d float-supported candidates re-evaluated exactly, "
        "%d of them supported in exact arithmetic" % (N1, SEED, c1, stats["supported"]))
N1b = 4000
c2 = run_batch(lambda: (rnd_above_one(), rnd_unit(), rnd_above_one(), rnd_unit()), N1b, False)
c3 = run_batch(lambda: (rnd_unit(), -math.exp(rng.uniform(-6, 3)), -math.exp(rng.uniform(-6, 3)),
                        rnd_above_one()), N1b, False)
ck.info("N1b: %d + %d exact samples in regions I and II (with mass near the boundary values)" % (c2, c3))
ck.check("N1 (EXACT at sampled points): contacts vanish; on regions I and II: A/g2 > 0, B/g2 > 0, the five "
         "kernel ratios are positive, and g22/g2 < 0 (I), g17/g2 < 0 (II); no supported sample anywhere",
         stats["violations"] == 0 and stats["contact_fail"] == 0 and stats["supported"] == 0
         and stats["I"] > 1000 and stats["II"] > 1000,
         "region I %d, region II %d samples; supported %d; min A/g2 = %.3e, min B/g2 = %.3e"
         % (stats["I"], stats["II"], stats["supported"], float(stats["minA"]), float(stats["minB"])))

# N2: chart-free solutions of the six contact equalities
nprng = np.random.default_rng(SEED)
N2 = 300
found = used = 0
worst_H = worst_real = worst_rebuild = 0.0
for it in range(N2):
    res = H.find_contact_point(T, nprng, tries=3)
    if res is None:
        continue
    found += 1
    a, cc, kap, gf = res
    if not H.labels_well_separated(a, cc):
        continue
    Fv = H.float_products(a, cc)
    # Moebius map: C0 -> infinity, C1 -> 0, A1 -> 1
    mob = H.moebius_to(cc[0], cc[1], a[1])
    xn, _, a2n, a3n = [mob(u) for u in a]
    vn, wn = mob(cc[2]), mob(cc[3])
    if max(abs(u) for u in (xn, a2n, a3n, vn, wn)) > 1e4 or abs(vn - 1) < 1e-3:
        continue
    used += 1
    worst_H = max(worst_H, abs((kap * Fv[0]).real) / abs(Fv[0]))
    kn = (a2n - 1) / (vn - 1)
    hn = (1 - wn) * a2n / (a2n - wn)
    ln = Fv[1] / Fv[0]
    worst_real = max(worst_real, max(abs(u.imag) / (1 + abs(u)) for u in (xn, kn, hn, ln)))
    mu_n = (kap * Fv[0]).imag
    Af, Cf = labels_f(xn.real, hn.real, kn.real, ln.real, vn)
    g_chart = H.float_ratio_gaps(Af, Cf, mu_n)
    worst_rebuild = max(worst_rebuild, np.max(np.abs(gf - g_chart)) / np.max(np.abs(gf)))
ck.check("N2 (NUMERIC): at chart-free solutions of the six contact equalities with well separated labels: "
         "Re(kappa F0) = 0, the parameters x, h, k, l are real, and the chart reproduces all 24 gaps",
         used >= 30 and worst_H < 1e-7 and worst_real < 1e-6 and worst_rebuild < 1e-6,
         "%d searches, %d converged, %d usable; max |Re(kappa F0)|/|F0| = %.1e; max |Im(x,h,k,l)|/(1+|.|) = %.1e; "
         "max relative gap mismatch = %.1e" % (N2, found, used, worst_H, worst_real, worst_rebuild))

# N3: the dual criterion against nonnegative least squares
N3 = 3000
mism = nin = 0
Mall = np.array([H.F_VALUES[j] for j in gens], float)
for it in range(N3):
    lam = nprng.exponential(size=len(gens)) * (nprng.random(len(gens)) < 0.6)
    for j in nonminor:
        lam[gens.index(j)] *= nprng.choice([0.0, 0.3, 1.0, 3.0])
    gv = lam @ Mall
    if gv.max() == 0:
        continue
    gv = gv / gv.max()
    Aval = sum(cf * gv[j] for j, cf in A_COEF.items())
    Bval = sum(cf * gv[j] for j, cf in B_COEF.items())
    if abs(Aval) < 1e-7 or abs(Bval) < 1e-7:
        continue
    inside_dual = (Aval >= 0 and Bval >= 0)
    nin += inside_dual
    if inside_dual != (H.minor_cone_nnls_residual(gv, minor) < 1e-9):
        mism += 1
ck.check("N3 (NUMERIC): on %d random points of the contact cone, 'A >= 0 and B >= 0' agrees with "
         "a nonnegative-least-squares test of membership in the minor cone" % N3, mism == 0,
         "%d inside, mismatches %d" % (nin, mism))

# N4: direct search for a supported point in regions I and II (Section 5b says there is none):
# maximize the smallest off-contact gap, normalized by the largest |gap|, by Nelder-Mead.
from scipy.optimize import minimize   # noqa: E402

OFF = [j for j in range(24) if j not in T]


def neg_min_gap(par, reg, mu):
    ea, eb, ec, ed = np.exp(par[:4])
    if reg == "I":
        hf, kf, lf, xf = 1 / (1 + ea), 1 + eb, 1 / (1 + ec), 1 + ed
    else:
        hf, kf, lf, xf = -ea, -eb, 1 + ec, 1 / (1 + ed)
    try:
        Af, Cf = labels_f(xf, hf, kf, lf, complex(par[4], par[5]))
        gf = H.float_ratio_gaps(Af, Cf, mu)
    except ZeroDivisionError:
        return 1.0
    sc = np.max(np.abs(gf))
    if not np.isfinite(sc) or sc == 0:
        return 1.0
    return -float(np.min(gf[OFF])) / sc


N4 = 12
best = {}
for reg in ("I", "II"):
    b = -1.0
    for it in range(N4):
        p0 = nprng.normal(size=6) * 1.5
        for mu in (1, -1):
            r_ = minimize(neg_min_gap, p0, args=(reg, mu), method="Nelder-Mead",
                          options=dict(maxiter=3000, xatol=1e-10, fatol=1e-13))
            b = max(b, -r_.fun)
    best[reg] = b
ck.check("N4 (NUMERIC): Nelder-Mead maximization of the smallest normalized off-contact gap over regions I, II "
         "(%d starts each, both signs of mu) finds no supported interior point (value stays <= 1e-5)" % N4,
         best["I"] <= 1e-5 and best["II"] <= 1e-5,
         "best value: region I %.2e, region II %.2e" % (best["I"], best["II"]))

ck.finish()
