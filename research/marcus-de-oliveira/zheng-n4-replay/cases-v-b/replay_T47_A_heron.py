"""Independent replay, T_47 part A, documents 4 and 5: the two necessary inequalities H6 >= 0 and
H8 >= 0 on the cone of the paper's T_47 = {1,3,5,6,7,8}.

Proof documents (read as PDF, never executed):
  final_referee_T46_negative_j_Heron_uniform_20261023.pdf   (3 pages; Theorem 1: H6 >= 0)
  final_referee_T46_two_Heron_necessity_20261023.pdf        (2 pages; Theorem 1: H6 >= 0 and H8 >= 0)

Run:   MOC4_SUPP=/path/to/supplement  python replay_T47_A_heron.py [--no-numeric]

Exit code 0 iff every check passes.  Every polynomial (C0, Q, G, the principal minor, W, R8, T8)
is computed here from the paper's definitions; MOC4_SUPP is used only for sha256 comparison of
the PDFs and for comparing two stored polynomials (G and W) with the ones computed here.  Without
MOC4_SUPP these comparisons are skipped and everything else runs unchanged.
"""
import itertools
import json
import os
import random
import sys
import time
from fractions import Fraction as Fr

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import replay_T47_A_common as K          # noqa: E402

NUMERIC = "--no-numeric" not in sys.argv
ck = K.Checker("Replay of T_47 part A, documents 4 and 5: the Heron-type necessary inequalities H6 >= 0, H8 >= 0")
K.print_versions(ck)
import sympy as sp                        # noqa: E402

T = K.PAPER_T[47]
CANON = [27, 1, 59, 5, 2, 28, 6, 4]
CNAMES = ["q", "p", "j", "c", "e", "d", "b", "a"]

ck.section("Provenance")
have_supp = K.check_file_identity(ck, ["final_referee_T46_negative_j_Heron_uniform_20261023.pdf",
                                       "final_referee_T46_two_Heron_necessity_20261023.pdf",
                                       "root_moc4_T46_full_radial_20261022.json",
                                       "final_referee_T46_j_squareclass_20261022.json"])
K.check_foundations(ck)
ck.check("EXACT", "both notes use the tuple {1,3,5,6,7,8} = T_47 and the source q f27 + p f1 + j f59 + c f5 + e f2 + d f28 "
         "+ b f6 + a f4: the eight generators and names of the top-level note (all vanish on T; only f28 is non-minor)",
         all(K.F_VALUES[jx][s] == 0 for jx in CANON for s in T) and [jx for jx in CANON if not K.F_IS_MINOR[jx]] == [28])

# ----------------------------------------------------------------------------------------
ck.section("C0, Q, the radial polynomial G and its discriminant in j  [document 5, p. 1]")
t0 = time.time()
rad = K.Radial(CANON, CNAMES)
R = rad.R
q, p, j, c, e, d, b, a, h = R.gens
D = a + d
M = (c + d) * (a + d + e)
L = M + b * D
U = c + e
C0 = -(q + d) * L - j * M
Q = j * p - j * (c + e + d) - (c + d) * (e + d)
ck.check("EXACT", "C(g + h kappa) = C0 + h Q with C0 = -(q+d) L - j M and Q = j p - j(c+e+d) - (c+d)(e+d) "
         "(D = a+d, M = (c+d)(a+d+e), L = M + bD); all coefficients of C0 are negative and -d^3 occurs, so C0 < 0 for "
         "nonnegative coefficients with d > 0, hence Q != 0 at a point of X",
         rad.F == C0 + h * Q and all(co < 0 for _, co in R.terms(C0))
         and dict(R.terms(C0)).get((0, 0, 0, 0, 0, 3, 0, 0, 0)) == -1)
okr, G, S = rad.radial_factor(-C0, Q, "j")
N = rad.cleared(-C0, Q)
ck.check("EXACT", "Q^4 det Hess C(g - (C0/Q) kappa) = -4 Q^2 G^2 as a polynomial identity (%d terms on the left); G has 1556 "
         "terms, is homogeneous of degree 9 and has degree 2 in j" % len(N),
         okr and N == Q * Q * G * G * (-4) and len(G) == 1556 and K.total_degrees(R, G) == [9]
         and R.degrees(G)[2] == 2, "%.1f s" % (time.time() - t0))
A = K.coeff_in(R, G, "j", 2)
Bb = K.coeff_in(R, G, "j", 1)
Cc = K.coeff_in(R, G, "j", 0)
H6 = (L * (p + q) - b * D * U) ** 2 - p * q * M * L * 4 + c * e * b * D * M * 4
Xi = q * (c + e + d + b) - c * e
W0 = (q + d) * L * p + b * ((q - d) * L - a * e * (c + d))
H8 = W0 * W0 + d * D * b * b * L * Xi * 4
disc = Bb * Bb - A * Cc * 4
Gj = R.derivative(G, 2)
ck.check("EXACT", "(1): with G = A j^2 + B j + C, B^2 - 4AC = (q+d)^2 H6 H8 where "
         "H6 = [L(p+q) - bD(c+e)]^2 - 4pqML + 4cebDM (150 terms, degree 6) and "
         "H8 = W0^2 + 4dDb^2 L Xi, W0 = (q+d)Lp + b{(q-d)L - ae(c+d)}, Xi = q(c+e+d+b) - ce (279 terms, degree 8); "
         "(5): G_j^2 - 4AG = B^2 - 4AC",
         disc == (q + d) ** 2 * H6 * H8 and len(H6) == 150 and len(H8) == 279
         and K.total_degrees(R, H6) == [6] and K.total_degrees(R, H8) == [8] and Gj * Gj - A * G * 4 == disc)
ck.check("EXACT", "(4): gcd(A, B, C) = 1 and gcd(G, G_j) = 1 (flint multivariate gcd)",
         A.gcd(Bb).gcd(Cc) == R.one and G.gcd(Gj) == R.one)
ck.check("EXACT", "smooth wall of H8 (3): dH8/dp = 2 W0 (q+d) L and dH8/dq = 2 W0 dW0/dq + 4dDb^2 L (c+e+d+b) "
         "(L does not involve p, q)",
         R.derivative(H8, 1) == W0 * (q + d) * L * 2
         and R.derivative(H8, 0) == W0 * R.derivative(W0, 0) * 2 + d * D * b * b * L * (c + e + d + b) * 4
         and R.derivative(L, 0) == R.zero and R.derivative(L, 1) == R.zero)

# ----------------------------------------------------------------------------------------
ck.section("The principal minor deleting the coordinates 0, 1 and the norm identities  [document 5, p. 2]")
Hm = [[rad.H[i][k] for k in range(2, 14)] for i in range(2, 14)]
pm, _ = K.laplace_det(Hm, R.zero)
P10, rem = divmod(pm, j * j * (-1))
rng = random.Random(4711)
import flint
okm = True
for _ in range(5):
    x = [rng.randint(-9, 9) for _ in range(9)]
    Mi = flint.fmpz_mat(12, 12, [int(R.evaluate(Hm[i][k], x)) for i in range(12) for k in range(12)])
    okm &= int(Mi.det()) == R.evaluate(pm, x)
ck.check("EXACT", "the 12 x 12 principal minor of Hess C(g + h kappa) with rows and columns 0, 1 deleted does not involve "
         "the height h, has %d terms (the supplement reports 1411) and equals -j^2 P10 with P10 of degree 2 in j and total "
         "degree 10; equal to integer determinants at 5 random points" % len(pm),
         rem == R.zero and R.degrees(pm)[8] == 0 and len(pm) == 1411 and R.degrees(P10)[2] == 2
         and K.total_degrees(R, P10) == [10] and okm)
al = K.coeff_in(R, P10, "j", 2)
be = K.coeff_in(R, P10, "j", 1)
ga = K.coeff_in(R, P10, "j", 0)
u1 = A * be - al * Bb
u0 = A * ga - al * Cc
W2 = u0 * u0 - u1 * (Bb * ga - be * Cc)
ck.check("EXACT", "A P10 - alpha G = u1 j + u0 with u1 = A beta - alpha B, u0 = A gamma - alpha C; and "
         "u0^2 - u1 (B gamma - beta C) (%d terms) is the resultant of G and P10 in j" % len(W2),
         A * P10 - al * G == u1 * j + u0
         and W2 == (A * ga - al * Cc) ** 2 - (A * be - al * Bb) * (Bb * ga - be * Cc))
t1 = time.time()
W = K.fsqrt(R, W2)
okW = W is not None and W * W == W2 and R.degrees(W)[2] == 0 and W != R.zero
if okW:
    Sp = A * u0 * 2 - Bb * u1 + A * W * 2
    if divmod(Sp, H8)[1] != R.zero:
        W = W * (-1)                      # the document's sign of W is the one for which H8 divides S+
ck.check("EXACT", "this resultant is a perfect square W^2 with W of %d terms, homogeneous of degree 17 and independent of j "
         "(own exact square root, verified by squaring)" % (len(W) if W is not None else 0),
         okW and K.total_degrees(R, W) == [17], "%.1f s" % (time.time() - t1))
recs = {}
for sg, name, target in ((1, "-H8 R8^2", H8 * (-1)), (-1, "-(q+d)^2 H6 T8^2", (q + d) ** 2 * H6 * (-1))):
    Sg = A * u0 * 2 - Bb * u1 + A * W * (2 * sg)
    quo, rem = divmod(Sg, target)
    root = K.fsqrt(R, quo) if rem == R.zero else None
    lhs = P10 * Sg - (u1 * j + u0 + W * sg) ** 2
    cq, cr = divmod(lhs, G)
    recs[sg] = (len(root) if root is not None else 0, len(cq))
    ck.check("EXACT", "S%s := 2A u0 - B u1 %s 2AW = %s with a nonzero polynomial root of %d terms that is independent of j; "
             "and P10 S%s - (u1 j + u0 %s W)^2 is divisible by G (complete polynomial quotient, %d terms)"
             % ("+" if sg == 1 else "-", "+" if sg == 1 else "-", name, recs[sg][0], "+" if sg == 1 else "-",
                "+" if sg == 1 else "-", recs[sg][1]),
             rem == R.zero and root is not None and root != R.zero and R.degrees(root)[2] == 0 and root * root == quo
             and cr == R.zero, "%.1f s" % (time.time() - t1))
ck.check("EXACT", "sizes agree with the supplement's report: R8 has 621 terms, T8 has 332 terms, the two quotients have "
         "93806 and 108074 terms", recs.get(1) == (621, 93806) and recs.get(-1) == (332, 108074))
if have_supp:
    with open(os.path.join(K.supp_dir(), "root_moc4_T46_full_radial_20261022.json")) as fh:
        data = json.load(fh)
    with open(os.path.join(K.supp_dir(), "final_referee_T46_j_squareclass_20261022.json")) as fh:
        sdat = json.load(fh)

    def load(terms):
        return R.from_terms([(tuple(ex)[:8] + (0,) * (9 - len(tuple(ex)[:8])), int(co)) for ex, co in terms])
    G_st = [load(f["terms"]) for f in data["radial_factors"] if f["degree"] == 9][0]
    ck.check("EXACT", "supplement data (JSON, comparison only): the stored degree-nine factor equals G and the stored W "
             "equals W, both with the same sign", G_st == G and load(sdat["W_terms"]) == W)
sx, g1, g2, g3, q0, q1, q2 = sp.symbols("x g1 g2 g3 q0 q1 q2")
Gx = g1 * sx + g2 * sx ** 2 + g3 * sx ** 3                  # G along a curve through a zero of G
Qx = q0 + q1 * sx + q2 * sx ** 2
ser = -4 * Gx ** 2 / Qx ** 2
ck.check("EXACT", "along any curve of coefficients through a point with G = 0, Q != 0: -4 G^2/Q^2 vanishes to first order "
         "there and its second derivative is -8 G_x^2/Q^2 (Taylor coefficients)",
         sp.simplify(ser.subs(sx, 0)) == 0 and sp.simplify(sp.diff(ser, sx).subs(sx, 0)) == 0
         and sp.simplify(sp.diff(ser, sx, 2).subs(sx, 0) + 8 * g1 ** 2 / q0 ** 2) == 0)

# ----------------------------------------------------------------------------------------
ck.section("Document 4: the second completion of H6 and the exceptional locus  [p. 1]")
T0 = L * p + (b * D - M) * q - b * D * (c + e)
ck.check("EXACT", "(1) of document 4: H6 = T0^2 + 4bDM (q-c)(q-e), T0 = Lp + (bD - M) q - bD(c+e); dH6/dp = 2 L T0; "
         "dH6/dq = 2 T0 (bD - M) + 4bDM (2q - c - e)",
         H6 == T0 * T0 + b * D * M * (q - c) * (q - e) * 4 and R.derivative(H6, 1) == L * T0 * 2
         and R.derivative(H6, 0) == T0 * (b * D - M) * 2 + b * D * M * (q * 2 - c - e) * 4)
st, sd_, sb_, sa_ = sp.symbols("t d b a")
sy = dict(zip(CNAMES, sp.symbols("q p j c e d b a")))


def to_sympy(f):
    syms = [sy[n] for n in CNAMES]
    return sp.Add(*[co * sp.Mul(*[s ** e_ for s, e_ in zip(syms, ex[:8]) if e_]) for ex, co in R.terms(f)])


exc = {sy["p"]: st, sy["q"]: st, sy["c"]: st, sy["e"]: st, sy["d"]: sd_, sy["b"]: sb_, sy["a"]: sa_}
DH = sa_ + sd_
MH = (st + sd_) * (sa_ + sd_ + st)
LH = MH + sb_ * DH
WH = (st + sd_) * LH * st + sb_ * ((st - sd_) * LH - sa_ * st * (st + sd_))
ck.check("EXACT", "on q = c = e: T0 = L (p - q); on the locus p = q = c = e = t: H8 = W0^2 + 4dDb^2 L t (t + d + b), and at "
         "t = 0, W0 = -b d (a+d)(b+d) = -b d L",
         sp.expand(to_sympy(T0).subs({sy["c"]: sy["q"], sy["e"]: sy["q"]})
                   - (to_sympy(L).subs({sy["c"]: sy["q"], sy["e"]: sy["q"]}) * (sy["p"] - sy["q"]))) == 0
         and sp.expand(to_sympy(H8).subs(exc, simultaneous=True)
                       - (WH ** 2 + 4 * sd_ * DH * sb_ ** 2 * LH * st * (st + sd_ + sb_))) == 0
         and sp.expand(WH.subs(st, 0) + sb_ * sd_ * DH * (sd_ + sb_)) == 0
         and sp.expand(LH.subs(st, 0) - DH * (sd_ + sb_)) == 0)
ck.info("PROOF-READ (document 4): for b, D, M > 0 the zero set of H6 is smooth except possibly on p = q = c = e, where "
        "H8 > 0; at a point of X with H8 < 0 one gets H6 = 0, G_j = 0, A != 0 and grad G != 0, and the mechanism of "
        "document 5 applies on the nearby part G_j != 0 of the radial hypersurface, where H6 < 0 and H8 < 0.")

# ----------------------------------------------------------------------------------------
ck.section("The geometric inputs cited by both notes (paper, Appendix C)")
ck.info("CITED 1 (paper (C.2) from Howard-Millson-Snowden-Vakil, and Corollary C.5): every real point u with C(u) = 0 and "
        "det Hess C(u) = 0 lies on the cone Sigma over the secant variety of X, and rank Hess C <= 12 on Sigma.")
ck.info("CITED 2 (paper, Proposition C.6; PROOF.md 7.4): for x in X, Hess C(Re x) has at most six positive and six negative "
        "eigenvalues; hence every 12 x 12 principal minor of it is >= 0, and one is > 0 if the rank is 12.")
rng = random.Random(4712)


def rq():
    return K.GQ(Fr(rng.randint(-9, 9), rng.randint(1, 5)), Fr(rng.randint(-9, 9), rng.randint(1, 5)))


oki = True
detail = []
for _ in range(2):
    av = [rq() for _ in range(4)]
    bv = [rq() for _ in range(4)]
    x = [K.GQ.of(cc) for cc in K.matching_vector(av, bv, rq())]
    u = [cc.re for cc in x]
    Hx = K.hess_at(u)
    rk, inert = K.rank(Hx), K.exact_inertia(Hx)
    mins = [K.det_fraction([[Hx[i][k] for k in range(14) if k not in dl] for i in range(14) if i not in dl])
            for dl in itertools.combinations(range(14), 2)]
    oki &= rk == 12 and inert == (6, 6, 2) and min(mins) >= 0 and max(mins) > 0
    detail.append("rank %d, inertia %s, %d of 91 minors > 0, %d < 0" % (rk, inert, sum(1 for m in mins if m > 0),
                                                                       sum(1 for m in mins if m < 0)))
ck.check("EXACT", "spot check of CITED 2 at two random Gaussian-rational points x of X (seed 4712): Hess C(Re x) has rank "
         "12, inertia (6,6,2) (exact congruence diagonalisation) and all 91 principal 12 x 12 minors are >= 0",
         oki, "; ".join(detail))
p1 = [(rng.randint(-5, 5), rng.randint(1, 4)) for _ in range(8)]
p2 = [(rng.randint(-5, 5), rng.randint(1, 4)) for _ in range(8)]
r1, r2 = K.tangent_space_rows(p1), K.tangent_space_rows(p2)
ck.check("EXACT", "at two random rational points of X the tangent spaces have dimension 6 and their sum has dimension 12: "
         "the cone over the secant variety has dimension exactly 12 (Terracini; it is at most 6 + 6)",
         K.rank(r1) == 6 and K.rank(r2) == 6 and K.rank(r1 + r2) == 12)
ck.info("PROOF-READ (own remark): CITED 1 can be replaced, where the notes use it, by Corollary C.5 on Sigma and dim Sigma = 12: "
        "at a point u0 of Re X on the slice span(f's, kappa) at which grad G != 0, the germ of Sigma on the slice has "
        "dimension >= 12 + 9 - 14 = 7 and lies in the smooth 7-dimensional germ {C0 + hQ = 0, G = 0}, so the two germs are "
        "equal and rank Hess C <= 12 at all radial points near u0.")
ck.info("PROOF-READ (linear algebra used by the notes): a real symmetric matrix of rank 12 has a nonzero principal 12-minor, "
        "and all its nonzero principal 12-minors have one sign, (-1)^(number of negative eigenvalues); a matrix curve whose "
        "value has corank >= 3 has a determinant vanishing to order >= 3.")

# ----------------------------------------------------------------------------------------
ck.section("An elementary proof of both theorems (not in the documents)")
lhs = (c + e + d + b) ** 2 * H6
rhs = (c + e + d + b) ** 2 * T0 * T0 + b * D * M * (Xi - c * (c + d + b)) * (Xi - e * (e + d + b)) * 4
ck.check("EXACT", "(q - c)(c+e+d+b) = Xi - c(c+d+b) and (q - e)(c+e+d+b) = Xi - e(e+d+b), hence "
         "(c+e+d+b)^2 H6 = (c+e+d+b)^2 T0^2 + 4bDM (Xi - c(c+d+b)) (Xi - e(e+d+b))",
         (q - c) * (c + e + d + b) == Xi - c * (c + d + b) and (q - e) * (c + e + d + b) == Xi - e * (e + d + b)
         and lhs == rhs)
ck.info("PROOF-READ: let all coefficients be >= 0 and b, d > 0 (so D, M, L, c+e+d+b > 0).  If Xi <= 0, both brackets are "
        "<= 0 and H6 >= 0; if moreover Xi < 0 both are < 0 and H6 > 0.  If Xi > 0 then H8 >= 4dDb^2 L Xi > 0.  Hence "
        "H8 < 0 forces Xi < 0 and H6 > 0, and H6 < 0 forces Xi > 0 and H8 > 0: whenever one of H6, H8 is negative, "
        "H6 H8 < 0.  At a point of X with d > 0 we have G = 0, so by (5) (q+d)^2 H6 H8 = G_j^2 >= 0.  Therefore "
        "H6 >= 0 and H8 >= 0.  This uses only (1), (5), the two square completions and the radial identity; it does not "
        "use the principal minor, the norm identities, j > 0, or any result of Appendix C.  In this orthant "
        "{H6 >= 0, H8 >= 0} is exactly the set where the quadratic G in j has nonnegative discriminant.")


def heron_values(qv, pv, cv_, ev, dv, bv, av):
    """H6 and H8 from the compact formulas (equal to the polynomials by the identities above)."""
    Dv = av + dv
    Mv = (cv_ + dv) * (av + dv + ev)
    Lv = Mv + bv * Dv
    v6 = (Lv * (pv + qv) - bv * Dv * (cv_ + ev)) ** 2 - 4 * pv * qv * Mv * Lv + 4 * cv_ * ev * bv * Dv * Mv
    w0 = (qv + dv) * Lv * pv + bv * ((qv - dv) * Lv - av * ev * (cv_ + dv))
    v8 = w0 * w0 + 4 * dv * Dv * bv * bv * Lv * (qv * (cv_ + ev + dv + bv) - cv_ * ev)
    return v6, v8


rng = random.Random(4713)
okf = True
for _ in range(200):
    pt = [Fr(rng.randint(-9, 9), rng.randint(1, 4)) for _ in range(8)]
    okf &= heron_values(pt[0], pt[1], pt[3], pt[4], pt[5], pt[6], pt[7]) \
        == (R.evaluate(H6, pt + [Fr(0)]), R.evaluate(H8, pt + [Fr(0)]))
cnt = {"H6<0": 0, "H8<0": 0, "both<0": 0, "H6<0,H8<=0": 0, "H8<0,H6<=0": 0}
NPT = 300000
for _ in range(NPT):
    pt = [rng.randint(0, 40) for _ in range(7)]              # integer points suffice: H6, H8 are homogeneous
    if rng.random() < 0.3:
        pt[rng.randrange(7)] = 0
    pt[4] = rng.randint(1, 40)                               # d > 0
    pt[5] = rng.randint(1, 40)                               # b > 0
    v6, v8 = heron_values(*pt)
    cnt["H6<0"] += v6 < 0
    cnt["H8<0"] += v8 < 0
    cnt["both<0"] += v6 < 0 and v8 < 0
    cnt["H6<0,H8<=0"] += v6 < 0 and v8 <= 0
    cnt["H8<0,H6<=0"] += v8 < 0 and v6 <= 0
ck.check("EXACT", "exact evaluation at %d random integer points (q,p,c,e,d,b,a) of the closed orthant with b, d > 0 (seed "
         "4713; zero coefficients occur; the evaluator agrees with the polynomials at 200 rational points): H6 < 0 at %d "
         "points and H8 < 0 at %d points, never both, never one negative with the other <= 0"
         % (NPT, cnt["H6<0"], cnt["H8<0"]),
         okf and cnt["both<0"] == 0 and cnt["H6<0,H8<=0"] == 0 and cnt["H8<0,H6<=0"] == 0
         and cnt["H6<0"] > 1000 and cnt["H8<0"] > 1000)

# ----------------------------------------------------------------------------------------
if NUMERIC:
    import numpy as np
    ck.section("NUMERIC falsification attempts (floating point; evidence only)")
    Gt, C0t, Qt = R.terms(G), R.terms(C0), R.terms(Q)
    H6t, H8t, P10t = R.terms(H6), R.terms(H8), R.terms(P10)
    table = K.hess_linear_table()
    FC = np.array([[float(K.F_COORD[jx][k]) for k in range(14)] for jx in CANON])
    Gs = {k: [(ex[:2] + (0,) + ex[3:], co) for ex, co in Gt if ex[2] == k] for k in (0, 1, 2)}

    def hess_num(u):
        Hn = np.zeros((14, 14))
        for i in range(14):
            for k in range(14):
                for cf, idx in table[i][k]:
                    Hn[i, k] += cf * u[idx]
        return Hn

    # (i) CITED 1 on this slice: random real radial points, coefficients of either sign
    rng = np.random.default_rng(47500)
    rows = []
    while len(rows) < 400:
        x = rng.normal(size=8)
        Av, Bv, Cv = (K.float_eval(Gs[k], x) for k in (2, 1, 0))
        dsc = Bv * Bv - 4 * Av * Cv
        if dsc < 0 or abs(Av) < 1e-6:
            continue
        for sg in (1, -1):
            y = x.copy()
            y[2] = (-Bv + sg * np.sqrt(dsc)) / (2 * Av)
            qv = K.float_eval(Qt, y)
            if abs(qv) < 1e-3:
                continue
            u = y @ FC
            u[13] += -K.float_eval(C0t, y) / qv
            sv = np.linalg.svd(hess_num(u), compute_uv=False)
            rows.append((sv[11] / sv[0], sv[12] / sv[0], sv[13] / sv[0]))
    rows = np.array(rows)
    ck.check("NUMERIC", "CITED 1 on this slice: at %d random real points with C0 + hQ = 0 and G = 0 (coefficients of either "
             "sign, seed 47500) the two smallest singular values of Hess C are <= %.1e and <= %.1e times the largest, the "
             "twelfth is >= %.1e times it: rank exactly 12 at every radial point found"
             % (len(rows), rows[:, 1].max(), rows[:, 2].max(), rows[:, 0].min()),
             rows[:, 1].max() < 1e-8 and rows[:, 0].min() > 1e-7)
    # (ii) the two theorems on stratum points with d > 0
    found = []
    for dmin, seed in ((0.002, 47200), (0.01, 47200)):
        found += K.numeric_stratum_search(T, seed, 300, extra=lambda g, dmin=dmin: [10 * min(0.0, (g[10] - g[20]) - dmin)])
    res = []
    for pt in found:
        g = pt["gaps"]
        co = np.array([g[18], g[22], g[0], g[20], g[19], g[10] - g[20], g[9], g[14]])
        x = co / co[5]
        if x[6] < 1e-6 or x[2] < 1e-6:
            continue
        res.append((K.float_eval(H6t, x) / K.float_eval(H6t, x, True), K.float_eval(H8t, x) / K.float_eval(H8t, x, True),
                    abs(K.float_eval(Gt, x)) / K.float_eval(Gt, x, True)))
    res = np.array(res)
    ck.check("NUMERIC", "points of X with contact set T, d > 0 and b, j > 0 from the steered search (seed 47200, 300 starts "
             "for d >= 0.002 and for d >= 0.01): %d points; relative |G| <= %.1e; H6 >= %.2e and H8 >= %.2e (relative to "
             "the sum of the absolute monomials): both are >= 0 up to rounding at every point"
             % (len(res), res[:, 2].max(), res[:, 0].min(), res[:, 1].min()),
             len(res) >= 3 and res[:, 0].min() > -1e-9 and res[:, 1].min() > -1e-9)
    # (iii) the algebraic mechanism of document 5, tested where both factors can be negative (b < 0 or d < 0)
    Xs = np.exp(rng.normal(size=(400000, 8)) * 1.2)
    cntn = {}
    for neg, nm in ((0, "q"), (1, "p"), (3, "c"), (4, "e"), (7, "a")):
        Y = Xs.copy()
        Y[:, neg] *= -1
        v6, v8 = heron_values(Y[:, 0], Y[:, 1], Y[:, 3], Y[:, 4], Y[:, 5], Y[:, 6], Y[:, 7])
        cntn[nm] = int(((v6 < 0) & (v8 < 0)).sum())
    ck.check("NUMERIC", "400000 random points with b, d > 0 and exactly one of q, p, c, e, a negative (seed 47500): the number "
             "with H6 < 0 and H8 < 0 is %s: the case 'both factors negative' of document 5 does not occur near the "
             "domain of its theorem" % cntn, all(vv == 0 for vv in cntn.values()))
    hits = []
    for neg in (5, 6):                                       # d < 0 or b < 0
        Y = Xs[:60000].copy()
        Y[:, neg] *= -1
        v6, v8 = heron_values(Y[:, 0], Y[:, 1], Y[:, 3], Y[:, 4], Y[:, 5], Y[:, 6], Y[:, 7])
        for x in Y[(v6 < 0) & (v8 < 0)][:40]:
            Av, Bv, Cv = (K.float_eval(Gs[k], x) for k in (2, 1, 0))
            dsc = Bv * Bv - 4 * Av * Cv
            if dsc <= 0 or abs(Av) < 1e-9:
                continue
            for sg in (1, -1):
                y = x.copy()
                y[2] = (-Bv + sg * np.sqrt(dsc)) / (2 * Av)
                qv = K.float_eval(Qt, y)
                if abs(qv) < 1e-6 * K.float_eval(Qt, y, True) or abs(y[2]) < 1e-6:
                    continue
                u = y @ FC
                u[13] += -K.float_eval(C0t, y) / qv
                ev = np.linalg.eigvalsh(hess_num(u))
                tol = 1e-8 * np.abs(ev).max()
                hits.append((K.float_eval(P10t, y) / K.float_eval(P10t, y, True), int((ev > tol).sum()),
                             int((ev < -tol).sum())))
    hits_a = np.array(hits) if hits else np.zeros((0, 3))
    ck.check("NUMERIC", "the mechanism itself, at %d real radial points with H6 < 0 and H8 < 0 (d < 0 or b < 0, outside the "
             "domain of the theorem): P10 > 0 at all of them (smallest relative value %.1e), so the principal minor "
             "-j^2 P10 is negative there, and the inertia is never (6,6): the pairs (n+, n-) found are %s"
             % (len(hits), hits_a[:, 0].min() if len(hits) else float("nan"),
                sorted(set((int(r_[1]), int(r_[2])) for r_ in hits))),
             len(hits) >= 20 and hits_a[:, 0].min() > 0 and all((int(r_[1]), int(r_[2])) != (6, 6) for r_ in hits))

ck.finish()
