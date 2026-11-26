"""Independent replay of the proof of Case T_6 = {0,1,2,4,6,8}.

Document replayed: supplement, proofs/manin_transfer_20261021/
    independent_moc4_T2_uniform_20261022.pdf  (5 pages; read, nothing executed).
The document is named "T2" in the supplement's older numbering; its contact set is the
paper's T_6.

Run:  MOC4_SUPP=/path/to/supplement python replay_T6.py
(MOC4_SUPP is used only to read one JSON data file for the facet-numbering comparison.)
Exit status 0 iff every check passes.

Exact parts: sympy rational-function fields / polynomial rings, Fractions.
Section NUMERIC: numpy floating point, falsification attempt only.
"""
import random
import sys

sys.dont_write_bytecode = True

import replay_T3_common as C
from sympy.polys.domains import QQ
from sympy.polys.fields import field

T6 = (0, 1, 2, 4, 6, 8)
ck = C.Checker("replay_T6: Case T_6 = {0,1,2,4,6,8} (independent_moc4_T2_uniform_20261022.pdf)")
V = C.FVAL

# ---------------------------------------------------------------------------------------
# A. objects
# ---------------------------------------------------------------------------------------
C.crosscheck_numbering(ck)

ck.section("A. contact set, generators, finite linear algebra (document pp. 1-2)")
ck.check("document's T = {0,1,2,4,6,8} is literally the paper's T_6 (identity of G); it is NOT the paper's T_2 = {0,1,2,3,6,8}",
         set(T6) == {0, 1, 2, 4, 6, 8})
orb2 = any(C.g_set(gm, {0, 1, 2, 3, 6, 8}) == frozenset(T6) for gm in C.G_MAPS)
ck.check("(the paper's T_2 is in a different G-orbit, so the file name refers to the older numbering)", not orb2)
gens = C.generators(T6)
ck.check("generators vanishing on T: {0,1,2,3,4,5,6,9,10,21,22,23,25}", gens == [0, 1, 2, 3, 4, 5, 6, 9, 10, 21, 22, 23, 25], gens)
ck.check("only f22, f25 are not minor vectors", [j for j in gens if not C.IS_MINOR[j]] == [22, 25])
C.verify_contact_cone(ck, T6)
COORD = (3, 7, 9, 10, 12, 13, 16, 18)
ck.check("span has dimension 8 and q = (g3,g7,g9,g10,g12,g13,g16,g18) are coordinates on it",
         C.rank([V[j] for j in gens]) == 8 and C.rank([[V[j][s] for s in COORD] for j in gens]) == 8)
REL = {5: (1, 0, 0, 0, 0, 0, 0, 0), 11: (0, -1, 1, 1, 0, 0, 0, 0), 14: (0, 0, 0, 0, 1, 0, 0, 0),
       15: (-1, -1, 1, 0, 0, 1, 0, 0), 17: (-1, -1, 1, 0, 0, 0, 1, 0), 19: (0, -1, 0, 0, -1, 1, 0, 1),
       20: (0, -1, 0, 1, 0, 0, 0, 1), 21: (-1, -2, 1, 1, -1, 1, 0, 1), 22: (0, -1, 0, 0, -1, 0, 1, 1),
       23: (-1, -1, 1, 0, -1, 0, 1, 1)}
qof = lambda g: [g[s] for s in COORD]
ck.check("the ten linear formulas for g5, g11, g14, g15, g17, g19, g20, g21, g22, g23 in terms of q hold on the span",
         all(V[k][idx] == C.dot(cf, qof(V[k])) for k in gens for idx, cf in REL.items()))
A_co = (-1, -2, 1, 1, -1, 0, 1, 1)
B_co = (-1, -1, 1, 0, -1, 1, 0, 1)
A_of = lambda g: C.dot(A_co, qof(g))
B_of = lambda g: C.dot(B_co, qof(g))
minors = [j for j in gens if C.IS_MINOR[j]]
ck.check("the eleven minor generators have rank 8", C.rank([V[j] for j in minors]) == 8)
facets = C.cone_facets([qof(V[j]) for j in minors])
ck.check("cone of the eleven minors has exactly 18 facets", len(facets) == 18, len(facets))
ck.check("two of the facet functionals are A and B of (2)", C.primitive(A_co) in facets and C.primitive(B_co) in facets)
ck.check("each of the other sixteen is >= 0 on all thirteen generators",
         all(C.dot(l, qof(V[j])) >= 0 for l in facets if l not in (C.primitive(A_co), C.primitive(B_co)) for j in gens))
ck.info("=> for g in N_T: g in cone(minors) <=> A(g) >= 0 and B(g) >= 0.")
ck.check("on the two non-minor generators: A(f22) = 1, A(f25) = -1, B(f22) = -1, B(f25) = 1",
         (A_of(V[22]), A_of(V[25]), B_of(V[22]), B_of(V[25])) == (1, -1, -1, 1))

# ---------------------------------------------------------------------------------------
# B. the universal inequality used
# ---------------------------------------------------------------------------------------
ck.section("B. SOC member a f1 + b f9 + c f22, ab >= c^2 (document p. 1, eq. (1)): the representative itself")
rng = random.Random(20261125)
C.orbit_check(ck, "SOC", (1, 9, 22), (1, 9, 22), expected_orbit=144)
SWAP = ((0, 1, 2, 3), (0, 2, 1, 3), 0)              # sigma -> sigma o (1 2): relabel A_1 <-> A_2
swm = C.g_index_map(SWAP)
ck.check("the label swap A_1 <-> A_2 (sigma -> sigma o (1 2)) preserves T and fixes E02 = f4, E31 = f21, f3, f9, f23; "
         "it exchanges f1<->f2, f5<->f10, f22<->f25",
         C.g_set(swm, T6) == frozenset(T6)
         and [C.g_facet(swm, j) for j in (4, 21, 3, 9, 23, 1, 2, 5, 10, 22, 25)] == [4, 21, 3, 9, 23, 2, 1, 10, 5, 25, 22])
C.covariance_check(ck, "swap: (f1,f9;f22) <-> (f2,f9;f25)", SWAP, [(1, 2), (9, 9), (22, 25)], rng)
C.replay_SOC_representative(ck)

# ---------------------------------------------------------------------------------------
# C. collision strata
# ---------------------------------------------------------------------------------------
ck.section("C. closed collision strata (document p. 2)")
blocks = {(i, j): [k for k, s in enumerate(C.PERMS) if s[i] == j] for i in range(4) for j in range(4)}
used = {(i, C.PERMS[k][i]) for k in T6 for i in range(4)}
ck.check("the pairs (i, j) not used by a contact are exactly (0,2), (0,3), (3,0), (3,1)",
         sorted(set(blocks) - used) == [(0, 2), (0, 3), (3, 0), (3, 1)])
ck.check("for every block {sigma(i) = j}, all facet vectors vanishing on the block are minors (all 16 blocks)",
         all(all(C.IS_MINOR[j] for j in C.generators(b)) for b in blocks.values()))
C.verify_contact_cone(ck, blocks[(0, 0)])
ck.check("G is transitive on the 16 blocks (so the block cone of every (i,j) is generated by its minor facet vectors)",
         {C.g_set(gm, blocks[(0, 0)]) for gm in C.G_MAPS} == {frozenset(b) for b in blocks.values()})
ck.info("PROOF-READ: A_i = C_j with (i,j) used by a contact sigma_0 gives F_sigma0 = 0, so Re kappa F_0 = 0 and the whole block is")
ck.info("  contact; g then lies in the cone of that block, which is generated by minors.")
span_basis = C.rref([V[j] for j in gens])[0][:8]


def subspace_equal_on(block):
    """Basis of {g in span : g constant on the block}."""
    cons = [[bvec[block[0]] - bvec[k] for bvec in span_basis] for k in block[1:]]
    out = []
    for co in C.nullspace(cons):
        out.append([sum(co[m] * span_basis[m][s] for m in range(8)) for s in range(24)])
    return out


DEC = {(0, 2): {0: 19, 4: 12, 5: 10, 10: 7, 23: 3}, (0, 3): {0: 18, 4: 12, 9: 7, 23: 3},
       (3, 0): {3: 12, 6: 9, 9: 7, 21: 3}, (3, 1): {3: 12, 5: 10, 6: 15, 10: 7, 21: 3}}
for (i, j), dec in DEC.items():
    sub = subspace_equal_on(blocks[(i, j)])
    ok = all(all(C.IS_MINOR[f] for f in dec) for _ in [0])
    for g in sub:
        ok = ok and [sum(g[gi] * V[f][s] for f, gi in dec.items()) for s in range(24)] == list(g)
    ck.check("A_%d = C_%d: g = %s on {g in span: g constant on the block} (dim %d); all f are minors"
             % (i, j, " + ".join("g%d f%d" % (gi, f) for f, gi in dec.items()), len(sub)), ok and len(sub) > 0)

# ---------------------------------------------------------------------------------------
# D. the chart
# ---------------------------------------------------------------------------------------
ck.section("D. the four circles and the chart (document pp. 2-3)")
C.infinity_convention_check(ck, "A", 0)


def products(A, Cc):
    """F_sigma = prod_{i != 0} (A_i - C_sigma(i))  (A_0 at infinity, factor dropped)."""
    one = Cc[1] * 0 + 1
    F = []
    for s in C.PERMS:
        p = one
        for i in range(4):
            if A[i] is None:
                continue
            p = p * (A[i] - Cc[s[i]])
        F.append(p)
    return F


K0, c0, v, c2, c, z = field("c0,v,c2,c,z", QQ)
F = products([None, K0(0), K0(1), z], [c0, v, c2, c])
D = lambda i, j: F[i] - F[j]
ck.check("D06 = (z-c)(v-C0)(C2-1), D28 = C2 (z-c)(v-C0), D02 = -(z-c)(v-C2), D68 = -(z-c)(C0-C2), "
         "D24 = -z (c-C2)(v-1), D01 = -v (c-C2)(z-1), D14 = -(v-c)(z-C2): all nonzero for distinct points",
         D(0, 6) == (z - c) * (v - c0) * (c2 - 1) and D(2, 8) == c2 * (z - c) * (v - c0)
         and D(0, 2) == -(z - c) * (v - c2) and D(6, 8) == -(z - c) * (c0 - c2)
         and D(2, 4) == -z * (c - c2) * (v - 1) and D(0, 1) == -v * (c - c2) * (z - 1)
         and D(1, 4) == -(v - c) * (z - c2))
h_def = z * (v - 1) / (v * (z - 1))
t_def = (c - z) * (v - c2) / ((c - c2) * (v - z))
s_def = (c0 - c2) / (v - c2)
ck.check("D06/D28 = (C2-1)/C2, D02/D68 = 1/s, D24/D01 = h, D14/D02 = 1 - 1/t  (s, h, t the document's cross-ratios)",
         D(0, 6) / D(2, 8) == (c2 - 1) / c2 and D(0, 2) / D(6, 8) == 1 / s_def
         and D(2, 4) / D(0, 1) == h_def and D(1, 4) / D(0, 2) == 1 - 1 / t_def)
ck.info("PROOF-READ: the eight indices are contacts, so kappa D_ij is imaginary and the four ratios are real: C2 = r real,")
ck.info("  s, h, t real.  These are the four circles (A0,A1,A2,C2), (A0,C0,C1,C2), (A1,A2,A3,C1), (A3,C1,C2,C3).")
ck.check("chart (3) is exhaustive: z (1 + (h-1) v) = h v and c ((v-r) - t (v-z)) = z (v-r) - t r (v-z) identically "
         "for h, t defined from the points; 1 + (h-1)v = (v-1)/(z-1), (v-r) - t(v-z) = (v-r)(z-r)/(c-r)",
         z * (1 + (h_def - 1) * v) == h_def * v
         and c * ((v - c2) - t_def * (v - z)) == z * (v - c2) - t_def * c2 * (v - z)
         and 1 + (h_def - 1) * v == (v - 1) / (z - 1)
         and (v - c2) - t_def * (v - z) == (v - c2) * (z - c2) / (c - c2))

K, r, s, h, tau, v, vb, eta = field("r,s,h,tau,v,vb,eta", QQ)


def chart(r, s, h, t, v):
    z = h * v / (1 + (h - 1) * v)
    c = (z * (v - r) - t * r * (v - z)) / ((v - r) - t * (v - z))
    return z, c


def gam_of(r, s, h, tau, v, vb):
    """gamma_sigma with g_sigma = rho * gamma_sigma; kappa = i mu/((z-c)(v-r)), rho = mu Im v, t = 1/tau."""
    t = 1 / tau
    out = []
    for vv in (v, vb):
        z, c = chart(r, s, h, t, vv)
        F = products([None, vv * 0, vv * 0 + 1, z], [r + s * (vv - r), vv, r + 0 * vv, c])
        out.append([(f - F[0]) / ((z - c) * (vv - r)) for f in F])
    return [-(a - b) / (v - vb) for a, b in zip(*out)], out[0]


gam, R1 = gam_of(r, s, h, tau, v, vb)
zc, cc = chart(r, s, h, 1 / tau, v)
ck.check("inverse chart formulas: h = z(v-1)/(v(z-1)), t = (c-z)(v-r)/((c-r)(v-z))",
         zc * (v - 1) / (v * (zc - 1)) == h and (cc - zc) * (v - r) / ((cc - r) * (v - zc)) == 1 / tau)
ck.check("contact ratios R_sigma = (F_sigma - F_0)/((z-c)(v-r)) on T: 0, -tau/(h-1), 1, 1 - h tau/(h-1), (1-r)(1-s), 1 - r(1-s) (all real)",
         [R1[k] for k in T6] == [K(0), -tau / (h - 1), K(1), 1 - h * tau / (h - 1), (1 - r) * (1 - s), 1 - r * (1 - s)])
ck.check("on the chart the gaps vanish identically exactly on T", [k for k in range(24) if gam[k] == 0] == list(T6))
ck.check("the chart's gap vector satisfies the ten linear formulas (so it lies in the span)",
         all(gam[idx] == C.dot(cf, qof(gam)) for idx, cf in REL.items()))


def quantities(r, h, v, vb):
    N = v * vb
    O1 = (v - 1) * (vb - 1)
    a = h * r - h - r
    Dd = (a * v + r) * (a * vb + r)
    Kk = h * r * (r - 1) / Dd
    J = (a * N + r * (v + vb) - r) / Dd
    L = r / ((h - 1) * N)
    M = h * (r - 1) / ((h - 1) * O1)
    return N, O1, a, Dd, Kk, J, L, M


N, O1, a, Dd, Kk, J, L, M = quantities(r, h, v, vb)
al, be = 1 + J, 1 + h * J
E = s * (tau - 1)                                  # E / rho

ck.section("D. gap formulas (document p. 3); g_sigma = rho * gamma_sigma")
ck.check("g3 = rho K (tau-1), g7 = rho L tau (s-1), g10 = rho M tau (s-1), g12 = -rho s",
         gam[3] == Kk * (tau - 1) and gam[7] == L * tau * (s - 1) and gam[10] == M * tau * (s - 1) and gam[12] == -s)
ck.check("g9 = g7 + s g3, g11 = g10 + s g3, g18 = g7 + E, g20 = g10 + E",
         gam[9] == gam[7] + s * gam[3] and gam[11] == gam[10] + s * gam[3] and gam[18] == gam[7] + E and gam[20] == gam[10] + E)
ck.check("g15 = g12 + E J, g17 = g12 + E h J", gam[15] == gam[12] + E * J and gam[17] == gam[12] + E * h * J)
ck.check("g19 = g3 + E (alpha - K), g22 = g3 + E (beta - K)", gam[19] == gam[3] + E * (al - Kk) and gam[22] == gam[3] + E * (be - Kk))
ck.check("g21 = g10 + E alpha, g23 = g7 + E beta", gam[21] == gam[10] + E * al and gam[23] == gam[7] + E * be)
Av, Bv = A_of(gam), B_of(gam)
ck.check("A = g10 + E beta, B = g7 + E alpha", Av == gam[10] + E * be and Bv == gam[7] + E * al)
ck.check("(4): K + a J = 1 and M - L = J D/((h-1) N O1)", Kk + a * J == 1 and M - L == J * Dd / ((h - 1) * N * O1))
ck.check("L M = H/((h-1)^2 N O1) with H = h r (r-1): L M > 0 <=> H > 0 <=> K > 0", L * M == h * r * (r - 1) / ((h - 1) ** 2 * N * O1))
ck.check("A + B = g21 + g23 and A B - g21 g23 = -E (g7/L) J^2 D/(N O1)",
         Av + Bv == gam[21] + gam[23] and Av * Bv - gam[21] * gam[23] == -E * (gam[7] / L) * J ** 2 * Dd / (N * O1))

ck.section("D. Case III (document p. 3): p = -tau, q = 1/(s-1), scale q/(-rho)")
p3, q3 = -tau, 1 / (s - 1)
G3 = [-x * q3 for x in gam]
ck.check("g19 = (p+1)((q+1) alpha - K), g22 = (p+1)((q+1) beta - K)",
         G3[19] == (p3 + 1) * ((q3 + 1) * al - Kk) and G3[22] == (p3 + 1) * ((q3 + 1) * be - Kk))
ck.check("in these units E = (p+1)(q+1) > 0 and A = g10 + E beta, B = g7 + E alpha", -E * q3 == (p3 + 1) * (q3 + 1))

ck.section("D. Case I (document p. 4): p = 1/(tau-1), q = -1/s, l = -L, m = -M, scale p q/rho")
p1, q1 = 1 / (tau - 1), -1 / s
l, m = -L, -M
Rr = (1 + p1) * (1 + q1)
b1, c1 = l * Rr, m * Rr
G1 = [x * p1 * q1 for x in gam]
ck.check("g3 = K q, g12 = p, g7 = b, g10 = c, g9 = b - K, g11 = c - K, g18 = b - 1, g20 = c - 1",
         [G1[k] for k in (3, 12, 7, 10, 9, 11, 18, 20)] == [Kk * q1, p1, b1, c1, b1 - Kk, c1 - Kk, b1 - 1, c1 - 1])
ck.check("g15 = p - J, g17 = p - h J, g19 = K q - (alpha - K), g22 = K q - (beta - K), g21 = c - alpha, g23 = b - beta, A = c - beta, B = b - alpha",
         [G1[k] for k in (15, 17, 19, 22, 21, 23)] == [p1 - J, p1 - h * J, Kk * q1 - (al - Kk), Kk * q1 - (be - Kk), c1 - al, b1 - be]
         and A_of(G1) == c1 - be and B_of(G1) == b1 - al)
E02, E31 = V[4], V[21]
D0vec = [G1[k] - p1 * E02[k] - Kk * q1 * E31[k] - Rr * (l * V[10][k] + m * V[5][k]) for k in range(24)]


def free_of_s_tau(e):
    return all(mon[1] == 0 and mon[3] == 0 for mon in e.numer.monoms()) and all(mon[1] == 0 and mon[3] == 0 for mon in e.denom.monoms())


ck.check("positive-direction identity: g - p E02 - K q E31 - R (l f10 + m f5) does not depend on (p, q) (i.e. on s, tau)",
         all(free_of_s_tau(e) for e in D0vec))
oneb = (1 + (h - 1) * v) * (1 + (h - 1) * vb)
ck.check("K + L alpha = r^2 (r-1) |1+(h-1)v|^2/((h-1) N D) and K + M beta = h r (r-1)^2 |1+(h-1)v|^2/((h-1) O1 D)",
         Kk + L * al == r ** 2 * (r - 1) * oneb / ((h - 1) * N * Dd)
         and Kk + M * be == h * r * (r - 1) ** 2 * oneb / ((h - 1) * O1 * Dd))
# label swap
N2, O12, a2, Dd2, Kk2, J2, L2, M2 = quantities(1 - r, 1 / h, 1 - v, 1 - vb)
ck.check("label swap A_1 <-> A_2: (r, h, v) -> (1-r, 1/h, 1-v) sends (J, L, M, K) to (h J, M, L, K)",
         (J2, L2, M2, Kk2) == (h * J, M, L, Kk))
gam_sw, _ = gam_of(1 - r, s, 1 / h, tau, 1 - v, 1 - vb)
tauperm = (0, 2, 1, 3)
ck.check("label swap: the chart point (1-r, s, 1/h, t, 1-v) with the SAME rho has gaps g'_sigma = g_{sigma o (1 2)} (all 24)",
         all(gam_sw[k] == gam[C.PIDX[C.compose(C.PERMS[k], tauperm)]] for k in range(24)))
T1 = 1 + (2 - h) * J
ck.check("(5): K - l T1 = (l/m)(K - m alpha)", Kk - l * T1 == (l / m) * (Kk - m * al))
ck.check("l T1 - m alpha = (m - l)(K/m - alpha)", l * T1 - m * al == (m - l) * (Kk / m - al))
ck.check("f1 + f3 + f9 - f10 - f22 + f23 = 0 (the eta-freedom in (6) and (8))",
         all(V[1][k] + V[3][k] + V[9][k] - V[10][k] - V[22][k] + V[23][k] == 0 for k in range(24)))
rem = [G1[k] - (p1 - J) * E02[k] - (Kk * q1 - (al - Kk)) * E31[k] for k in range(24)]
rhs6 = [(eta - h * J) * V[1][k] + eta * V[3][k] + (c1 - al) * V[5][k] + (eta + 1) * V[9][k] + (b1 - eta - 1) * V[10][k]
        + (J - eta) * V[22][k] + (eta + 1 - Kk) * V[23][k] for k in range(24)]
ck.check("(6): g - (p-J) E02 - [Kq - (alpha-K)] E31 = (eta - hJ) f1 + eta f3 + (c-alpha) f5 + (eta+1) f9 + (b-eta-1) f10 "
         "+ (J-eta) f22 + (eta+1-K) f23, all 24 entries, eta free", rem == rhs6)
ck.check("(eta - hJ)(eta + 1) - (J - eta)^2 = T1 eta - J (J + h)", (eta - h * J) * (eta + 1) - (J - eta) ** 2 == T1 * eta - J * (J + h))
ck.check("alpha^2/T1 - 1 = J (J + h)/T1 and K - 1 = -a J with -a = r + h (1 - r)",
         al ** 2 / T1 - 1 == J * (J + h) / T1 and Kk - 1 == -a * J and -a == r + h * (1 - r))

ck.section("D. Case II (document p. 5): p = tau/(1-tau), q = s/(1-s), scale (1+p)(1+q)/(-rho)")
p2, q2 = tau / (1 - tau), s / (1 - s)
G2 = [-x * (1 + p2) * (1 + q2) for x in gam]
ck.check("g3 = K (1+q), g12 = q (1+p), g7 = L p, g10 = M p, g15 = q (p + alpha), g17 = q (p + beta), g19 = K + q alpha, g22 = K + q beta",
         [G2[k] for k in (3, 12, 7, 10, 15, 17, 19, 22)]
         == [Kk * (1 + q2), q2 * (1 + p2), L * p2, M * p2, q2 * (p2 + al), q2 * (p2 + be), Kk + q2 * al, Kk + q2 * be])
ck.check("g21 = c + q alpha, g23 = b + q beta, A = c + q beta, B = b + q alpha (b = L p, c = M p)",
         G2[21] == M * p2 + q2 * al and G2[23] == L * p2 + q2 * be and A_of(G2) == M * p2 + q2 * be and B_of(G2) == L * p2 + q2 * al)
xx = -J
P = h * xx - 1
T2 = 2 * P + al
ck.check("alpha = 1 - x, beta = -P, T2 = 2P + alpha = (2h-1) x - 1, P + alpha = (h-1) x  (x = -J, P = h x - 1)",
         al == 1 - xx and be == -P and T2 == (2 * h - 1) * xx - 1 and P + al == (h - 1) * xx)
ck.check("(7): K - M T2 = (M/L)(K - L P)", Kk - M * T2 == (M / L) * (Kk - L * P))
ck.check("M T2 - L P = (L - M)(K/L - P)", M * T2 - L * P == (L - M) * (Kk / L - P))
rem2 = [G2[k] - q2 * (p2 - P) * E02[k] - (Kk - q2 * P) * E31[k] for k in range(24)]
rem2s = [x / q2 for x in C.g_values(swm, rem2)]
Cs, Bs = L * p2 / q2, M * p2 / q2
rhs8 = [(eta + al) * V[1][k] + (eta + 1) * V[3][k] + (Cs - P) * V[5][k] + eta * V[9][k] + (Bs - eta) * V[10][k]
        + (P - eta) * V[22][k] + (eta + Kk) * V[23][k] for k in range(24)]
ck.check("(8): [g - q(p-P) E02 - (K - qP) E31] with A_1, A_2 swapped, divided by q, = (eta+alpha) f1 + (eta+1) f3 + (C*-P) f5 "
         "+ eta f9 + (B*-eta) f10 + (P-eta) f22 + (eta+K) f23, all 24 entries, eta free", rem2s == rhs8)
et = P * P / T2
ck.check("eta = P^2/T2: eta + alpha = (P+alpha)^2/T2, P - eta = P (P+alpha)/T2, eta (eta+alpha) = (P - eta)^2",
         et + al == (P + al) ** 2 / T2 and P - et == P * (P + al) / T2 and et * (et + al) == (P - et) ** 2)

# ---------------------------------------------------------------------------------------
# D'. EXACT instances
# ---------------------------------------------------------------------------------------
ck.section("D'. EXACT instances built from the paper's definitions (finite Gaussian-rational spectra)")
from fractions import Fraction as Fr

GQ = C.GQ
SWL = [C.PIDX[C.compose(C.PERMS[k], tauperm)] for k in range(24)]


def quant_exact(rq, hq, vq):
    Nq, Oq = vq.norm2(), (vq - 1).norm2()
    aq = hq * rq - hq - rq
    Dq = (aq * vq + rq).norm2()
    return (hq * rq * (rq - 1) / Dq, (aq * Nq + 2 * rq * vq.re - rq) / Dq, rq / ((hq - 1) * Nq), hq * (rq - 1) / ((hq - 1) * Oq))


def exact_T6(label, rq, sq, hq, tq, vq, case):
    zq = hq * vq / (1 + (hq - 1) * vq)
    cq = (zq * (vq - rq) - tq * rq * (vq - zq)) / ((vq - rq) - tq * (vq - zq))
    a_s, b_s, alpha, gp = C.exact_paper_gaps([None, GQ(0), GQ(1), zq], [rq + sq * (vq - rq), vq, GQ(rq), cq], 0, 1)
    C.exact_instance_report(ck, label, gp, T6)
    Aq, Bq = A_of(gp), B_of(gp)
    tauq = 1 / tq
    rho = -gp[12] / sq
    Kq, Jq, Lq, Mq = quant_exact(rq, hq, vq)
    ck.check("%s: A = %s, B = %s, one of them < 0 (not in the minor cone); g3, g7, g10 formulas hold with rho := -g12/s"
             % (label, C.fl(Aq), C.fl(Bq)), (Aq < 0 or Bq < 0) and gp[3] == rho * Kq * (tauq - 1)
             and gp[7] == rho * Lq * tauq * (sq - 1) and gp[10] == rho * Mq * tauq * (sq - 1))
    gg = list(gp)
    swapped = False
    if case == "I":
        ok = sq < 0 and tauq > 1 and Lq < 0 and rho > 0 and 0 < rq < 1 and hq < 0
        if Jq < 0:
            rq, hq, vq = 1 - rq, 1 / hq, 1 - vq
            gg = [gp[k] for k in SWL]
            Kq, Jq, Lq, Mq = quant_exact(rq, hq, vq)
            swapped = True
        pq, qq = 1 / (tauq - 1), -1 / sq
        G = [x * pq * qq / rho for x in gg]
        Rq = (1 + pq) * (1 + qq)
        bq, cq_ = -Lq * Rq, -Mq * Rq
        alq = 1 + Jq
        T1q = 1 + (2 - hq) * Jq
        etq = max(Fr(0), Kq - 1, Jq * (Jq + hq) / T1q)
        cf = {1: etq - hq * Jq, 3: etq, 5: cq_ - alq, 9: etq + 1, 10: bq - etq - 1, 22: Jq - etq, 23: etq + 1 - Kq,
              4: pq - Jq, 21: Kq * qq - (alq - Kq)}
        target = G
    else:
        ok = 0 < sq < 1 and 0 < tauq < 1 and Lq > 0 and rho < 0
        if rq < 0:
            rq, hq, vq = 1 - rq, 1 / hq, 1 - vq
            gg = [gp[k] for k in SWL]
            Kq, Jq, Lq, Mq = quant_exact(rq, hq, vq)
            swapped = True
        ok = ok and rq > 1 and hq > 1
        pq, qq = tauq / (1 - tauq), sq / (1 - sq)
        G = [x * (1 + pq) * (1 + qq) / (-rho) for x in gg]
        alq, beq = 1 + Jq, 1 + hq * Jq
        ok = ok and beq < 0
        Pq = -hq * Jq - 1
        T2q = 2 * Pq + alq
        etq = Pq * Pq / T2q
        ok = ok and pq - Pq >= 0 and Kq - qq * Pq >= 0
        rem_ = [G[k] - qq * (pq - Pq) * V[4][k] - (Kq - qq * Pq) * V[21][k] for k in range(24)]
        target = [x / qq for x in C.g_values(swm, rem_)]
        cf = {1: etq + alq, 3: etq + 1, 5: Lq * pq / qq - Pq, 9: etq, 10: Mq * pq / qq - etq, 22: Pq - etq, 23: etq + Kq}
    ck.check("%s: case %s%s; the document's decomposition reproduces all 24 entries exactly, all coefficients >= 0, "
             "SOC margin (coef f1)(coef f9) - (coef f22)^2 = %s >= 0"
             % (label, case, " (after the label swap)" if swapped else "", C.fl(cf[1] * cf[9] - cf[22] ** 2)),
             ok and [sum(cf[j] * V[j][k] for j in cf) for k in range(24)] == target
             and all(val >= 0 for val in cf.values()) and cf[1] * cf[9] >= cf[22] ** 2)


exact_T6("instance (r,s,h,t,v) = (3/4,-1/4,-1,1/2,1+i)", Fr(3, 4), Fr(-1, 4), Fr(-1), Fr(1, 2), GQ(1, 1), "I")
exact_T6("instance (r,s,h,t,v) = (1/2,-1/4,-1,3/4,1+i)", Fr(1, 2), Fr(-1, 4), Fr(-1), Fr(3, 4), GQ(1, 1), "I")
exact_T6("instance (r,s,h,t,v) = (2,1/4,2,2,-1/2-2i)", Fr(2), Fr(1, 4), Fr(2), Fr(2), GQ(Fr(-1, 2), -2), "II")
exact_T6("instance (r,s,h,t,v) = (5/2,1/2,2,7/2,-1/2+3i/2)", Fr(5, 2), Fr(1, 2), Fr(2), Fr(7, 2), GQ(Fr(-1, 2), Fr(3, 2)), "II")

# ---------------------------------------------------------------------------------------
# E. NUMERIC falsification attempts
# ---------------------------------------------------------------------------------------
ck.section("E. NUMERIC falsification attempts (numpy, seeds printed): build the document's certificate for every supporting sample")
import numpy as np

SWIDX = np.array([C.PIDX[C.compose(C.PERMS[k], tauperm)] for k in range(24)])
FV = np.array(V, dtype=float)


def sample(rs, ss, hs, ts, vs, rngn):
    """Gaps from FINITE spectra: random Moebius image of A = (inf, 0, 1, z), C = (r + s(v-r), v, r, c)."""
    zs = hs * vs / (1 + (hs - 1) * vs)
    cs = (zs * (vs - rs) - ts * rs * (vs - zs)) / ((vs - rs) - ts * (vs - zs))
    n = len(rs)
    al_, be_, ga_, de_ = [rngn.normal(size=n) + 1j * rngn.normal(size=n) for _ in range(4)]
    mob = lambda q_: (al_ * q_ + be_) / (ga_ * q_ + de_)
    a_ = np.stack([al_ / ga_, mob(0 * zs), mob(0 * zs + 1), mob(zs)])
    b_ = -np.stack([mob(rs + ss * (vs - rs)), mob(vs), mob(rs + 0j), mob(cs)])
    Fs = np.ones((24, n), dtype=complex)
    for k, sg in enumerate(C.PERMS):
        for i in range(4):
            Fs[k] *= a_[i] + b_[sg[i]]
    d1 = Fs[1] - Fs[0]
    kap = 1j * np.conj(d1) / np.abs(d1)
    g = (kap * (Fs - Fs[0])).real
    return g / np.max(np.abs(g), axis=0)


def quant_num(r_, h_, v_):
    N_ = abs(v_) ** 2
    O_ = abs(v_ - 1) ** 2
    a_ = h_ * r_ - h_ - r_
    D_ = abs(a_ * v_ + r_) ** 2
    return (h_ * r_ * (r_ - 1) / D_, (a_ * N_ + 2 * r_ * v_.real - r_) / D_, r_ / ((h_ - 1) * N_), h_ * (r_ - 1) / ((h_ - 1) * O_))


def certify(g, r_, s_, h_, t_, v_):
    """Returns (kind, None) if the document's certificate is confirmed numerically, else (kind, reason)."""
    tol = 1e-8
    q_ = g[list(COORD)]
    A_ = float(np.dot(A_co, q_))
    B_ = float(np.dot(B_co, q_))
    tau_ = 1 / t_
    rho = -g[12] / s_                               # g12 = -rho s
    K_, J_, L_, M_ = quant_num(r_, h_, v_)
    E_ = rho * s_ * (tau_ - 1)
    if E_ * L_ <= 0:
        return ("EL<=0", None if (A_ >= -tol and B_ >= -tol) else "EL <= 0 but A or B < 0")
    if s_ > 1 and tau_ < 0 and L_ > 0 and rho < 0:
        return ("III", None if (A_ > -tol and B_ > -tol) else "case III but A or B < 0")
    if s_ < 0 and tau_ > 1 and L_ < 0 and rho > 0:
        outer = (r_ > 1 and 0 < h_ < 1) or (r_ < 0 and h_ > 1)
        if outer:
            return ("I outer", None if (A_ >= -tol and B_ >= -tol) else "case I outer but A or B < 0")
        if not (0 < r_ < 1 and h_ < 0):
            return ("I", "case I with (r,h) in none of the three regions")
        if abs(J_) < 1e-12:
            return ("I inner J=0", None if (A_ >= -tol and B_ >= -tol) else "J = 0 but A or B < 0")
        gg = g
        if J_ < 0:                                   # label swap
            r_, h_, v_ = 1 - r_, 1 / h_, 1 - v_
            gg = g[SWIDX]
            K_, J_, L_, M_ = quant_num(r_, h_, v_)
        p_, qq_ = 1 / (tau_ - 1), -1 / s_
        G = gg * (p_ * qq_ / rho)
        l_, m_ = -L_, -M_
        R_ = (1 + p_) * (1 + qq_)
        b_, c_ = l_ * R_, m_ * R_
        al_ = 1 + J_
        T1_ = 1 + (2 - h_) * J_
        eta_ = max(0.0, K_ - 1, J_ * (J_ + h_) / T1_)
        coef = {1: eta_ - h_ * J_, 3: eta_, 5: c_ - al_, 9: eta_ + 1, 10: b_ - eta_ - 1, 22: J_ - eta_, 23: eta_ + 1 - K_,
                4: p_ - J_, 21: K_ * qq_ - (al_ - K_)}
        rec = sum(cf * FV[j] for j, cf in coef.items())
        scale = max(1.0, np.max(np.abs(G)))
        if np.max(np.abs(rec - G)) > 1e-6 * scale:
            return ("I inner", "decomposition (6) does not reproduce the gap vector")
        if min(coef.values()) < -1e-7 * scale:
            return ("I inner", "negative coefficient in (6): %r" % (min(coef.items(), key=lambda kv: kv[1]),))
        if coef[1] * coef[9] - coef[22] ** 2 < -1e-7 * scale ** 2:
            return ("I inner", "SOC condition fails in (6)")
        if not (m_ > l_ > 0 and T1_ > 0 and b_ >= al_ ** 2 / T1_ - 1e-7 * scale):
            return ("I inner", "intermediate claim m > l > 0, T1 > 0, b >= alpha^2/T1 fails")
        return ("I inner (SOC needed)" if (A_ < 0 or B_ < 0) else "I inner", None)
    if 0 < s_ < 1 and 0 < tau_ < 1 and L_ > 0 and rho < 0:
        gg = g
        if not ((r_ > 1 and h_ > 1) or (r_ < 0 and 0 < h_ < 1)):
            return ("II", "case II with (r,h) in neither region")
        if r_ < 0:
            r_, h_, v_ = 1 - r_, 1 / h_, 1 - v_
            gg = g[SWIDX]
            K_, J_, L_, M_ = quant_num(r_, h_, v_)
        p_, qq_ = tau_ / (1 - tau_), s_ / (1 - s_)
        G = gg * ((1 + p_) * (1 + qq_) / (-rho))
        be_ = 1 + h_ * J_
        al_ = 1 + J_
        Ag = float(np.dot(A_co, G[list(COORD)]))
        Bg = float(np.dot(B_co, G[list(COORD)]))
        if be_ >= 0:
            return ("II beta>=0", None if (Ag >= -tol and Bg >= -tol) else "case II, beta >= 0 but A or B < 0")
        x_ = -J_
        P_ = h_ * x_ - 1
        T2_ = 2 * P_ + al_
        eta_ = P_ * P_ / T2_
        Cs_, Bs_ = L_ * p_ / qq_, M_ * p_ / qq_
        coef = {1: eta_ + al_, 3: eta_ + 1, 5: Cs_ - P_, 9: eta_, 10: Bs_ - eta_, 22: P_ - eta_, 23: eta_ + K_}
        rem_ = (G - qq_ * (p_ - P_) * FV[4] - (K_ - qq_ * P_) * FV[21])[SWIDX] / qq_
        rec = sum(cf * FV[j] for j, cf in coef.items())
        scale = max(1.0, np.max(np.abs(rem_)))
        if np.max(np.abs(rec - rem_)) > 1e-6 * scale:
            return ("II beta<0", "decomposition (8) does not reproduce the remainder")
        if min(coef.values()) < -1e-7 * scale or qq_ * (p_ - P_) < -1e-7 or K_ - qq_ * P_ < -1e-7:
            return ("II beta<0", "negative coefficient in (8) or in the singleton costs")
        if coef[1] * coef[9] - coef[22] ** 2 < -1e-7 * scale ** 2:
            return ("II beta<0", "SOC condition fails in (8)")
        if not (L_ > M_ > 0 and T2_ > 0):
            return ("II beta<0", "intermediate claim L > M > 0, T2 > 0 fails")
        return ("II beta<0 (SOC needed)" if (Ag < 0 or Bg < 0) else "II beta<0", None)
    return ("unclassified", "EL > 0 but none of the cases I, II, III")


def run_family(name, seed, n, gen):
    rngn = np.random.default_rng(seed)
    rs, ss, hs, ts, vs = gen(rngn, n)
    g = sample(rs, ss, hs, ts, vs, rngn)
    cerr = np.max(np.abs(g[list(T6)]), axis=0)
    good = np.isfinite(cerr) & (cerr < 1e-6)
    counts, viol, nsupp = {}, [], 0
    for sign in (1.0, -1.0):
        gg = sign * g
        idx = np.nonzero(good & (np.min(gg, axis=0) >= -1e-9))[0]
        nsupp += len(idx)
        for i in idx:
            kind, why = certify(gg[:, i], rs[i], ss[i], hs[i], ts[i], vs[i])
            counts[kind] = counts.get(kind, 0) + 1
            if why is not None:
                viol.append((kind, why, int(i)))
    ck.check("NUMERIC %s: seed %d, %d samples on the contact locus, %d supporting; certificate confirmed for all"
             % (name, seed, int(good.sum()), nsupp), not viol, viol[:3])
    ck.info("   " + ", ".join("%s: %d" % kv for kv in sorted(counts.items())))
    return counts


def gen_broad(rngn, n):
    sg = lambda: rngn.choice([-1.0, 1.0], size=n)
    mag = lambda: np.exp(rngn.normal(scale=1.2, size=n))
    return (sg() * mag(), sg() * mag(), sg() * mag(), sg() * mag(), 2 * rngn.normal(size=n) + 1j * rngn.normal(size=n))


def gen_I(rngn, n):
    mag = lambda: np.exp(rngn.normal(scale=1.2, size=n))
    return (rngn.uniform(0, 1, size=n), -mag(), -mag(), rngn.uniform(0, 1, size=n),
            2 * rngn.normal(size=n) + 1j * rngn.normal(size=n))


def gen_II(rngn, n):
    mag = lambda: np.exp(rngn.normal(scale=1.2, size=n))
    rs, hs = 1 + mag(), 1 + mag()
    half = rngn.integers(0, 2, size=n).astype(bool)             # the unswapped region r < 0, 0 < h < 1
    rs = np.where(half, 1 - rs, rs)
    hs = np.where(half, 1 / hs, hs)
    return (rs, rngn.uniform(0, 1, size=n), hs, 1 + mag(), 2 * rngn.normal(size=n) + 1j * rngn.normal(size=n))


c1_ = run_family("broad chart sample", 61001, 300000, gen_broad)
c2_ = run_family("case I region 0<r<1, h<0, s<0, 0<t<1", 61002, 300000, gen_I)
c3_ = run_family("case II region (both (r,h) regions), 0<s<1, t>1", 61003, 300000, gen_II)
ck.check("NUMERIC: tests not vacuous (supporting samples that need the SOC inequality, in case I and in case II)",
         c2_.get("I inner (SOC needed)", 0) > 100 and c3_.get("II beta<0 (SOC needed)", 0) > 100,
         (c2_.get("I inner (SOC needed)", 0), c3_.get("II beta<0 (SOC needed)", 0)))

ck.finish()
