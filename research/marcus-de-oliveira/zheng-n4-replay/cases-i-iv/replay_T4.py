"""Independent replay of the proof of Case T_4 = {0,1,4,5,7,8}.

Document replayed: supplement, proofs/manin_transfer_20261021/
    root_moc4_T4_uniform_20261022.pdf  (5 pages; read, nothing executed),
and, for the four-coefficient inequality it imports,
    final_referee_T5_four_budget_capacity_20261022.pdf (3 pages; read).

Run:  MOC4_SUPP=/path/to/supplement python replay_T4.py
(MOC4_SUPP is used only to read one JSON data file for the facet-numbering comparison.)
Exit status 0 iff every check passes.

Exact parts: sympy rational-function fields / polynomial rings, Fractions.
Section NUMERIC: numpy floating point, falsification attempt only.
"""
import random
import sys

sys.dont_write_bytecode = True

import replay_T3_common as C
from sympy.polys.domains import QQ, ZZ
from sympy.polys.fields import field
from sympy.polys.rings import ring

T4 = (0, 1, 4, 5, 7, 8)
ck = C.Checker("replay_T4: Case T_4 = {0,1,4,5,7,8} (root_moc4_T4_uniform_20261022.pdf)")

# ---------------------------------------------------------------------------------------
# A. objects
# ---------------------------------------------------------------------------------------
C.crosscheck_numbering(ck)

ck.section("A. contact set, generators, the table of minors (document p. 1)")
ck.check("document's T_4 = {0,1,4,5,7,8} is literally the paper's T_4 (identity of G)", set(T4) == {0, 1, 4, 5, 7, 8})
gens = C.generators(T4)
ck.check("generators vanishing on T_4: {0,1,2,3,4,5,6,11,33,34}", gens == [0, 1, 2, 3, 4, 5, 6, 11, 33, 34], gens)
ck.check("only f34 is not a minor vector", [j for j in gens if not C.IS_MINOR[j]] == [34])
C.verify_contact_cone(ck, T4)
doc_table = {0: ((0,), (3,)), 1: ((0, 1), (2, 3)), 2: ((0, 2), (2, 3)), 3: ((0, 3), (2, 3)), 4: ((0,), (2,)),
             5: ((0, 1), (1, 3)), 6: ((3,), (0,)), 11: ((0, 2), (1, 2)), 33: ((0, 1), (0, 2))}
ck.check("the document's table f_j = m_{I,J} (0-based I, J) agrees with Appendix D for j = 0..6, 11, 33",
         all(C.minor_spec(j) == IJ for j, IJ in doc_table.items()))
V = C.FVAL
ck.check("f34 = f0 - f2 + f6 + f33 (values at all 24 permutations)",
         all(V[34][s] == V[0][s] - V[2][s] + V[6][s] + V[33][s] for s in range(24)))
ck.check("f0 - f1 - f2 - f3 + f4 = 0", all(V[0][s] - V[1][s] - V[2][s] - V[3][s] + V[4][s] == 0 for s in range(24)))
ck.check("span of the generators has dimension 8", C.rank([V[j] for j in gens]) == 8)


def on_span(label, lhs, rhs):
    """Check a linear identity lhs(g) = rhs(g) between 24-vectors or numbers for every
    generator g (hence on the span)."""
    ck.check(label, all(lhs(V[k]) == rhs(V[k]) for k in gens))


d_of = lambda g: g[2] - g[13]

ck.section("A. the face decomposition (2) and the case d <= 0 (document p. 2)")
on_span("(2): g = g19 f0 + g16 f1 + g14 f3 + g10 f5 + g15 f6 + g6 f11 + g13 f33 + d f34, d = g2 - g13, on the span",
        lambda g: tuple(g),
        lambda g: tuple(g[19] * V[0][s] + g[16] * V[1][s] + g[14] * V[3][s] + g[10] * V[5][s] + g[15] * V[6][s]
                        + g[6] * V[11][s] + g[13] * V[33][s] + d_of(g) * V[34][s] for s in range(24)))
on_span("g22 = g16 + g19 + d, g18 = g14 + g19 + d, g9 = g15 + d on the span",
        lambda g: (g[22], g[18], g[9]),
        lambda g: (g[16] + g[19] + d_of(g), g[14] + g[19] + d_of(g), g[15] + d_of(g)))
on_span("d <= 0 replacement at tau = 0: g = (g19 + d) f0 + g16 f1 + (-d) f2 + g14 f3 + g10 f5 + g9 f6 + g6 f11 + g2 f33 on the span",
        lambda g: tuple(g),
        lambda g: tuple((g[19] + d_of(g)) * V[0][s] + g[16] * V[1][s] - d_of(g) * V[2][s] + g[14] * V[3][s]
                        + g[10] * V[5][s] + g[9] * V[6][s] + g[6] * V[11][s] + g[2] * V[33][s] for s in range(24)))
ck.info("PROOF-READ: with A = -d >= 0, b = g19, tau = max(0, A - b), adding tau (f0 - f1 - f2 - f3 + f4) = 0 gives coefficients")
ck.info("  (b-A+tau, g16-tau, A-tau, g14-tau, tau) on (f0..f4): all >= 0 because g22 = g16+b-A >= 0, g18 = g14+b-A >= 0;")
ck.info("  the others are g10, g9, g6, g2 >= 0.  So d <= 0 gives a minor decomposition.")

# ---------------------------------------------------------------------------------------
# B. the universal inequality used
# ---------------------------------------------------------------------------------------
ck.section("B. four-coefficient member used: a f6 + b f0 + c f5 + d f34 (document p. 2, eq. (1))")
rng = random.Random(20261125)
g = C.orbit_check(ck, "four-coefficient", (0, 6, 10, 54), (6, 0, 5, 34), expected_orbit=288)
C.covariance_check(ck, "four (0,6;10;54) -> (6,0;5;34)", g, [(0, 6), (6, 0), (10, 5), (54, 34)], rng)
gdoc = ((0, 1, 2, 3), (0, 2, 1, 3), 1)
gm = C.g_index_map(gdoc)
ck.check("the document's transport 'transposition and the swap (1 2)', i.e. sigma -> sigma^{-1} o (1 2), maps "
         "(f0, f6, f10, f54) to (f6, f0, f5, f34)",
         tuple(C.g_facet(gm, j) for j in (0, 6, 10, 54)) == (6, 0, 5, 34))
C.covariance_check(ck, "the document's group element", gdoc, [(0, 6), (6, 0), (10, 5), (54, 34)], rng)
C.replay_four_coefficient_note(ck)

ck.section("B. the 'arrowhead pencil' of p. 2 (alternative derivation named by the document)")
Rz, a, b, c, d = ring("a,b,c,d", ZZ)
cap = a * b * (c - d) - d * d * (a + b + c + d)
ck.check("d (a+d) c (b+d) - d^2 [c(b+d) + (a+d)(b+d) + (a+d)c] = d [ab(c-d) - d^2(a+b+c+d)]  (Schur complement formula)",
         d * (a + d) * c * (b + d) - d * d * (c * (b + d) + (a + d) * (b + d) + (a + d) * c) == d * cap)
ck.info("The PDF gives no identity linking this 4x4 pencil to the Leibniz sums (it refers to its checker); the")
ck.info("inequality (1) is therefore taken from the four-budget note (replayed above) plus the G-transport.")

# ---------------------------------------------------------------------------------------
# C. the chart
# ---------------------------------------------------------------------------------------
ck.section("C. normalisation and contact ratios (document pp. 2-3)")
C.infinity_convention_check(ck, "C", 2)


def products(A, Cc):
    """F_sigma = - prod_{sigma(i) != 2} (A_i - C_sigma(i))  (the document's sign)."""
    one = A[1] * 0 + 1
    F = []
    for s in C.PERMS:
        p = -one
        for i in range(4):
            if Cc[s[i]] is None:
                continue
            p = p * (A[i] - Cc[s[i]])
        F.append(p)
    return F


K0, z, w, v, x, t = field("z,w,v,x,t", QQ)
F = products([z, K0(0), K0(1), x], [w, v, None, t * v])
r_def = z * (w - v) * (1 - t * v) / (v * (z - w))
u_def = -(x - t * v) * (z * (1 + v - w) - v) / (v * (z - w))
ck.check("contact ratios R_sigma = (F_sigma - F_0)/(v(z-w)) on T_4 = (0,1,4,5,7,8): 0, 1-x, t-x, x(t-1), r-x+1, u "
         "(identities in complex z, w, v, x, t; r, u as defined on p. 3)",
         [(F[s] - F[0]) / (v * (z - w)) for s in T4] == [K0(0), 1 - x, t - x, x * (t - 1), r_def - x + 1, u_def])
ck.info("PROOF-READ: kappa = i m/(v(z-w)), m in C*: contacts 1,4,5 <=> m(1-x), m(t-x), m x(t-1) real; then m(t-1) real,")
ck.info("  x real (t != 1), m real (x != 1), t real; contacts 7, 8 <=> r, u real.")
ck.check("chart-free identity: (F_2 - F_0)/(v(z-w)) = t - x/v, hence g2 = Re(i mu (t - x/v)) = -mu x Im(v)/|v|^2 for ALL z, w",
         (F[2] - F[0]) / (v * (z - w)) == t - x / v)
ck.info("FINDING (a): for real v this gives g2 = 0, so d = g2 - g13 = -g13 <= 0: real-v configurations fall under the case")
ck.info("  d <= 0 already settled by minors on p. 2.  The document's appeal to the three-plus-three circle theorem (p. 3) is")
ck.info("  therefore not needed: d > 0 forces Im v != 0.")
lam = (z - v) / (z - w)
ck.check("with lam = (z-v)/(z-w):  r v/(1 - t v) = z (lam - 1)  and  r v/(1 - t v) - u v/(x - t v) = lam  (identities)",
         r_def * v / (1 - t * v) == z * (lam - 1) and r_def * v / (1 - t * v) - u_def * v / (x - t * v) == lam)
ck.info("FINDING (b): for real v, x, t and real r, u these give lam real and z (lam - 1) real; lam = 1 is the collision w = v,")
ck.info("  so z is real, then w = z - (z - v)/lam is real, kappa is imaginary and ALL gaps vanish.  With eight distinct")
ck.info("  points the real-v case is thus trivial, exactly as in the documents for T_3, T_5, T_6.")
den_z = lambda r, u, v, x, t: r * t * v ** 2 - r * v * x + t ** 2 * v ** 2 - t * u * v ** 2 - t * v * x - t * v + u * v + x
den_w = lambda r, u, v, x, t: -r * t * v + r * x + t * u * v - u
ck.check("exhaustiveness: z * den_z = r v (t v - x) and w * den_w = (x - t v)(r - t v + 1) identically for r, u defined from z, w",
         z * den_z(r_def, u_def, v, x, t) == r_def * v * (t * v - x)
         and w * den_w(r_def, u_def, v, x, t) == (x - t * v) * (r_def - t * v + 1))
ck.check("in terms of the points: den_z = -(v-w)(x-tv)(tv-1)/(z-w) and den_w = -(z-v)(x-tv)(tv-1)/(v(z-w)): nonzero for distinct points",
         den_z(r_def, u_def, v, x, t) == -(v - w) * (x - t * v) * (t * v - 1) / (z - w)
         and den_w(r_def, u_def, v, x, t) == -(z - v) * (x - t * v) * (t * v - 1) / (v * (z - w)))
ck.info("PROOF-READ: den_z = 0 would force r v (tv - x) = 0 (collision); den_w = 0 would force r + 1 = t v, impossible for")
ck.info("  real r, t != 0 and nonreal v.  (For real v the chart is not used.)")

K, x, t, u, r, v, vb = field("x,t,u,r,v,vb", QQ)


def chart(x, t, u, r, v):
    z = r * v * (t * v - x) / den_z(r, u, v, x, t)
    w = (-t * v + x) * (r - t * v + 1) / den_w(r, u, v, x, t)
    return z, w


def R_of(x, t, u, r, v):
    z, w = chart(x, t, u, r, v)
    F = products([z, v * 0, v * 0 + 1, x], [w, v, None, t * v])
    return [(f - F[0]) / (v * (z - w)) for f in F]


R1, R2 = R_of(x, t, u, r, v), R_of(x, t, u, r, vb)
gam = [-(p - q) / (v - vb) for p, q in zip(R1, R2)]          # g_sigma = rho * gam_sigma
zc, wc = chart(x, t, u, r, v)
ck.check("the chart inverts the definitions of r and u",
         zc * (wc - v) * (1 - t * v) / (v * (zc - wc)) == r
         and -(x - t * v) * (zc * (1 + v - wc) - v) / (v * (zc - wc)) == u)
ck.check("on the chart R_sigma on T_4 is 0, 1-x, t-x, x(t-1), r-x+1, u and the gaps vanish identically exactly on T_4",
         [R1[s] for s in T4] == [K(0), 1 - x, t - x, x * (t - 1), r - x + 1, u]
         and [s for s in range(24) if gam[s] == 0] == list(T4))
ck.check("the chart's gap vector satisfies (2) in all 24 entries (so it lies in the span)",
         [gam[19] * V[0][s] + gam[16] * V[1][s] + gam[14] * V[3][s] + gam[10] * V[5][s] + gam[15] * V[6][s]
          + gam[6] * V[11][s] + gam[13] * V[33][s] + d_of(gam) * V[34][s] for s in range(24)] == gam)

ck.section("C. basic identities (document p. 3); g_sigma = rho * gamma_sigma, rho = mu Im v")
N = v * vb
ck.check("g2 = -rho x/|v|^2", gam[2] == -x / N)
ck.check("g6 = -rho r t (x-1)/|tv-1|^2", gam[6] == -r * t * (x - 1) / ((t * v - 1) * (t * vb - 1)))
ck.check("c = g10 = rho t u x/|tv-x|^2", gam[10] == t * u * x / ((t * v - x) * (t * vb - x)))
ck.check("g9 = (x-1) c", gam[9] == (x - 1) * gam[10])
ck.check("a = g15 = g9 - d", gam[15] == gam[9] - d_of(gam))

ck.section("C. the two source identities (3), (4) and E = (h-y)(H+hy) + hyL (document pp. 3-4)")
G = [-q for q in gam]                                   # gaps at rho = -1
a_, b_, c_, d_, e_, k_ = G[15], G[19], G[10], d_of(G), G[14], G[16]
y, h = x - 1, 1 - t
H = h + u
E = H * (1 - t * x) - h * r * y
J = H * (1 - t * x) + h * y
L = y - r + u
ck.check("(3a) k = h c - t (b + d)", k_ == h * c_ - t * (b_ + d_))
ck.check("(3b) E e = t y H b - J d", E * e_ == t * y * H * b_ - J * d_)
ck.check("(4)  J = E + h y (r+1) and J - x h H = -u y", J == E + h * y * (r + 1) and J - x * h * H == -u * y)
ck.check("E = (h - y)(H + h y) + h y L", E == (h - y) * (H + h * y) + h * y * L)
ck.check("(h - y)(1 + y) = 1 - t x - y^2  (used for E >= H(h-y)(1+y) when r <= H y/h)", (h - y) * (1 + y) == 1 - t * x - y * y)
ck.check("y - H y/h + u = -u (y - h)/h  (used for L > 0)", y - H * y / h + u == -u * (y - h) / h)

ck.section("C. the Hermitian forms (document p. 4)")
Qf = lambda s: -t * (r + t - u) * s ** 2 + (r * t + r * x + t * x + t - u) * s - x * (r + 1)
Kq = t * y * b_ - x * h * d_
ck.check("Q(v) v = (z - v) den_z: Q(v) != 0 on the chart (it vanishes only for the collision A_0 = C_1)",
         Qf(v) * v == (zc - v) * den_z(r, u, v, x, t))


def hermitian(expr, name):
    """expr must be a polynomial of bidegree <= (1,1) in (v, vb), symmetric: returns
    (m00, m01, m11) with expr = m00 + m01 (v + vb) + m11 v vb."""
    ok = expr.denom == 1
    num = expr.numer
    co = {}
    for mon, cf in num.terms():
        key = (mon[4], mon[5])
        rest = tuple(mon[:4]) + (0, 0)
        co[key] = co.get(key, 0) + K.ring({rest: cf})
    ok = ok and set(co) <= {(0, 0), (1, 0), (0, 1), (1, 1)}
    z0 = K.ring(0)
    m00, m10, m01, m11 = co.get((0, 0), z0), co.get((1, 0), z0), co.get((0, 1), z0), co.get((1, 1), z0)
    ok = ok and m10 == m01
    ck.check("%s is a polynomial of bidegree <= (1,1) in (v, conj v), symmetric, with coefficients in Q[x,t,u,r]" % name, ok)
    return K(m00), K(m01), K(m11)


m00, m01, m11 = hermitian(d_ * Qf(v) * Qf(vb), "D = d |Q(v)|^2")
ck.check("(6) det M_D = -r^3 t^2 u x y L E", m00 * m11 - m01 * m01 == -r ** 3 * t ** 2 * u * x * y * L * E)
ck.check("(6) D(x/t, x/t) = r u x y (h r - H y)", m00 + 2 * m01 * (x / t) + m11 * (x / t) ** 2 == r * u * x * y * (h * r - H * y))
n00, n01, n11 = hermitian(Kq * Qf(v) * Qf(vb), "P = K |Q(v)|^2, K = t y b - x h d")
ck.check("(7) det M_P = -r^2 t^2 u^2 x^2 h y L E", n00 * n11 - n01 * n01 == -r ** 2 * t ** 2 * u ** 2 * x ** 2 * h * y * L * E)
ck.check("(7) P(x/t, x/t) = r u^2 x^2 y (y - h)", n00 + 2 * n01 * (x / t) + n11 * (x / t) ** 2 == r * u ** 2 * x ** 2 * y * (y - h))
ck.info("PROOF-READ: a real symmetric 2x2 matrix with det >= 0 and a negative (positive) value on a real vector is negative")
ck.info("  semidefinite (positive definite if det > 0); its Hermitian form at (1, v) is m00 + m01(v + conj v) + m11 |v|^2.")

ck.section("C. the closing identities (document p. 4)")
Rc, xx, tt, bb, cc, dd = ring("x,t,b,c,d", ZZ)
yy, hh = xx - 1, 1 - tt
aa = yy * cc - dd
ck.check("with a = y c - d: ab(c-d) - d^2(a+b+c+d) = c {y b c - x d (b+d)}",
         aa * bb * (cc - dd) - dd * dd * (aa + bb + cc + dd) == cc * (yy * bb * cc - xx * dd * (bb + dd)))
ck.check("h {y b c - x d (b+d)} = (b+d)(t y b - x h d) + y b (h c - t(b+d))   [= (b+d) K + y b k]",
         hh * (yy * bb * cc - xx * dd * (bb + dd))
         == (bb + dd) * (tt * yy * bb - xx * hh * dd) + yy * bb * (hh * cc - tt * (bb + dd)))
ck.info("=> capacity = (c/h) [(b+d) K + y b k] > 0 for c, h, d > 0, b, k >= 0, y > 0, K > 0; a, b, c, d >= 0: (1) applies to")
ck.info("   a f6 + b f0 + c f5 + d f34 with a = g15, b = g19, c = g10, d = g2 - g13, the four terms of (2).")

# ---------------------------------------------------------------------------------------
# C'. EXACT instances
# ---------------------------------------------------------------------------------------
ck.section("C'. EXACT instances built from the paper's definitions (finite Gaussian-rational spectra)")
from fractions import Fraction as Fr

GQ = C.GQ


def exact_T4(label, xq, tq, uq, rq, vq, want_E):
    zq = rq * vq * (tq * vq - xq) / den_z(rq, uq, vq, xq, tq)
    wq = (-tq * vq + xq) * (rq - tq * vq + 1) / den_w(rq, uq, vq, xq, tq)
    a_s, b_s, alpha, gp = C.exact_paper_gaps([zq, GQ(0), GQ(1), GQ(xq)], [wq, vq, None, tq * vq], 0, 1)
    C.exact_instance_report(ck, label, gp, T4)
    dq = gp[2] - gp[13]
    aq, bq, cq = gp[15], gp[19], gp[10]
    cf = {0: gp[19], 1: gp[16], 3: gp[14], 5: gp[10], 6: gp[15], 11: gp[6], 33: gp[13], 34: dq}
    capq = aq * bq * (cq - dq) - dq * dq * (aq + bq + cq + dq)
    ck.check("%s: d = %s > 0 (the non-minor term is present); (2) reproduces all 24 gaps with coefficients >= 0; "
             "capacity ab(c-d) - d^2(a+b+c+d) = %s > 0" % (label, C.fl(dq), C.fl(capq)),
             dq > 0 and all(val >= 0 for val in cf.values())
             and [sum(cf[j] * V[j][s] for j in cf) for s in range(24)] == gp and capq > 0)
    yq, hq = xq - 1, 1 - tq
    Hq = hq + uq
    Eq = Hq * (1 - tq * xq) - hq * rq * yq
    Lq = yq - rq + uq
    rho = -gp[2] * vq.norm2() / xq
    ck.check("%s: E = %s (%s); rho := -g2 |v|^2/x < 0; g6, g10 formulas, k = hc - t(b+d), E e = t y H b - J d hold; "
             "K = t y b - x h d > 0%s" % (label, Eq, want_E, "; L > 0 and y > h" if Eq < 0 else ""),
             {"E<0": Eq < 0, "E>0": Eq > 0, "E=0": Eq == 0}[want_E] and rho < 0
             and gp[6] == -rho * rq * tq * (xq - 1) / (tq * vq - 1).norm2()
             and gp[10] == rho * tq * uq * xq / (tq * vq - xq).norm2()
             and gp[16] == hq * cq - tq * (bq + dq)
             and Eq * gp[14] == tq * yq * Hq * bq - (Hq * (1 - tq * xq) + hq * yq) * dq
             and tq * yq * bq - xq * hq * dq > 0 and (Eq >= 0 or (Lq > 0 and yq > hq)))


exact_T4("instance (x,t,u,r,v) = (2,1/2,-1/4,1/4,7/2+i/4)", Fr(2), Fr(1, 2), Fr(-1, 4), Fr(1, 4), GQ(Fr(7, 2), Fr(1, 4)), "E<0")
exact_T4("instance (x,t,u,r,v) = (2,1/4,-3/8,1/8,21/4+3i/4)", Fr(2), Fr(1, 4), Fr(-3, 8), Fr(1, 8), GQ(Fr(21, 4), Fr(3, 4)), "E>0")
exact_T4("boundary instance (x,t,u,r,v) = (2,1/4,-3/8,1/4,21/4+i/2)", Fr(2), Fr(1, 4), Fr(-3, 8), Fr(1, 4), GQ(Fr(21, 4), Fr(1, 2)), "E=0")

# ---------------------------------------------------------------------------------------
# D. NUMERIC falsification attempts
# ---------------------------------------------------------------------------------------
ck.section("D. NUMERIC falsification attempts on the intermediate claims (numpy, seeds printed)")
import numpy as np


def sample(xs, ts, us, rs, vs, rngn):
    """Gaps from FINITE spectra: random Moebius image of (A, C) = ((z,0,1,x), (w,v,inf,tv)),
    a_i = m(A_i), b_j = -m(C_j), F_sigma = prod (a_i + b_sigma(i)), direction with Re kappa(F_1 - F_0) = 0."""
    dz = den_z(rs, us, vs, xs, ts)
    dw = den_w(rs, us, vs, xs, ts)
    zs = rs * vs * (ts * vs - xs) / dz
    ws = (-ts * vs + xs) * (rs - ts * vs + 1) / dw
    n = len(xs)
    al, be, ga, de = [rngn.normal(size=n) + 1j * rngn.normal(size=n) for _ in range(4)]
    mob = lambda s: (al * s + be) / (ga * s + de)
    a = np.stack([mob(zs), mob(0 * zs), mob(0 * zs + 1), mob(xs + 0j)])
    b = -np.stack([mob(ws), mob(vs), al / ga, mob(ts * vs)])
    Fs = np.ones((24, n), dtype=complex)
    for k, s in enumerate(C.PERMS):
        for i in range(4):
            Fs[k] *= a[i] + b[s[i]]
    d1 = Fs[1] - Fs[0]
    kap = 1j * np.conj(d1) / np.abs(d1)
    g = (kap * (Fs - Fs[0])).real
    return g / np.max(np.abs(g), axis=0)


def run_family(name, seed, n, gen):
    rngn = np.random.default_rng(seed)
    xs, ts, us, rs, vs = gen(rngn, n)
    g = sample(xs, ts, us, rs, vs, rngn)
    tol = 1e-9
    cerr = np.max(np.abs(g[list(T4)]), axis=0)
    good = np.isfinite(cerr) & (cerr < 1e-6)
    stats = dict(n=int(good.sum()), supp=0, dpos=0, Eneg=0)
    worst = dict(K=np.inf, cap=np.inf, L=np.inf, yh=np.inf)
    viol = []
    for sign in (1.0, -1.0):
        gg = sign * g
        idx = np.nonzero(good & (np.min(gg, axis=0) >= -tol))[0]
        stats["supp"] += len(idx)
        for i in idx:
            hv = gg[:, i]
            xv, tv, uv, rv = xs[i], ts[i], us[i], rs[i]
            dv = hv[2] - hv[13]
            if dv <= 1e-9:
                continue
            stats["dpos"] += 1
            av, bv, cv = hv[15], hv[19], hv[10]
            if not (xv > 1 and 0 < tv < 1 and rv > 0 and uv < 0):
                viol.append(("sign claims x>1, 0<t<1, r>0, u<0 fail", i, xv, tv, rv, uv))
                continue
            yv, hh_ = xv - 1, 1 - tv
            Hv = hh_ + uv
            Ev = Hv * (1 - tv * xv) - hh_ * rv * yv
            Kv = tv * yv * bv - xv * hh_ * dv
            capv = av * bv * (cv - dv) - dv * dv * (av + bv + cv + dv)
            worst["K"] = min(worst["K"], Kv)
            worst["cap"] = min(worst["cap"], capv)
            if Kv <= 0:
                viol.append(("K <= 0", i))
            if capv < 0:
                viol.append(("capacity < 0", i))
            if Ev < 0:
                stats["Eneg"] += 1
                Lv = yv - rv + uv
                worst["L"] = min(worst["L"], Lv)
                worst["yh"] = min(worst["yh"], yv - hh_)
                if Lv <= 0 or yv <= hh_:
                    viol.append(("E<0 but L<=0 or y<=h", i))
    ck.check("NUMERIC %s: seed %d, %d samples on the contact locus, %d supporting, %d with d > 0 (%d with E < 0); "
             "no violation of the intermediate claims" % (name, seed, stats["n"], stats["supp"], stats["dpos"], stats["Eneg"]),
             not viol, viol[:3])
    ck.info("   min K %.3e | min capacity (normalised gaps) %.3e | min L (E<0) %.3e | min y-h (E<0) %.3e"
            % (worst["K"], worst["cap"], worst["L"], worst["yh"]))
    return stats


def gen_broad(rngn, n):
    sg = lambda: rngn.choice([-1.0, 1.0], size=n)
    mag = lambda: np.exp(rngn.normal(scale=1.2, size=n))
    return (sg() * mag(), sg() * mag(), sg() * mag(), sg() * mag(), rngn.normal(size=n) + 1j * rngn.normal(size=n))


def gen_target(rngn, n):
    """The region forced by d > 0: x > 1, 0 < t < 1, r > 0, u < 0."""
    xs = 1 + np.exp(rngn.normal(scale=1.2, size=n))
    ts = rngn.uniform(0, 1, size=n)
    rs = np.exp(rngn.normal(scale=1.2, size=n))
    us = -np.exp(rngn.normal(scale=1.2, size=n))
    return (xs, ts, us, rs, rngn.normal(size=n) + 1j * rngn.normal(size=n))


def gen_dpos(rngn, n, mode=0):
    """Where supporting sources with d > 0 actually occur (found by a first broad search): x > 1, 0 < t < 1,
    small r > 0, small u < 0, Re v between x and x/t, small |Im v|.  mode 1: r chosen so that E = 0 exactly;
    mode 2: x = 2 - t (y = h exactly); mode 3: r = y + u (L = 0 exactly)."""
    xs = 1 + np.exp(rngn.normal(scale=1.0, size=n))
    ts = rngn.uniform(0, 1, size=n)
    if mode == 2:
        xs = 2 - ts
    rs = np.exp(rngn.normal(loc=-1.5, scale=1.2, size=n))
    us = -np.exp(rngn.normal(loc=-1.2, scale=1.2, size=n))
    if mode == 1:
        rs = (1 - ts + us) * (1 - ts * xs) / ((1 - ts) * (xs - 1))
    if mode == 3:
        rs = (xs - 1) + us
    vre = xs + (xs / ts - xs) * rngn.uniform(0, 1, size=n)
    vs = vre + 1j * rngn.choice([-1.0, 1.0], size=n) * np.exp(rngn.normal(loc=-2, scale=1.5, size=n))
    return (xs, ts, us, rs, vs)


s1 = run_family("broad chart sample", 41001, 200000, gen_broad)
s2 = run_family("region x>1, 0<t<1, r>0, u<0", 41002, 200000, gen_target)
s3 = run_family("d>0 region (Re v in (x, x/t), small r, |u|, |Im v|)", 41003, 400000, gen_dpos)
s4 = run_family("d>0 region with E = 0 exactly", 41004, 200000, lambda g_, n: gen_dpos(g_, n, 1))
s5 = run_family("d>0 region with y = h exactly", 41005, 200000, lambda g_, n: gen_dpos(g_, n, 2))
s6 = run_family("d>0 region with L = 0 exactly", 41006, 200000, lambda g_, n: gen_dpos(g_, n, 3))
ck.check("NUMERIC: tests not vacuous (many supporting samples with d > 0, with E < 0 and with E >= 0)",
         s3["dpos"] > 1000 and s3["Eneg"] > 1000 and s3["dpos"] - s3["Eneg"] > 1000,
         (s3["dpos"], s3["Eneg"]))

ck.finish()
