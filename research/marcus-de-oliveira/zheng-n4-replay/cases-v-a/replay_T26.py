"""Independent replay of case T_26 = {2,3,5,7,8,12} (group V) of D. Zheng's proof of the
Marcus-de Oliveira conjecture for n = 4.

Proof document (read as PDF, never executed): final_referee_T25_uniform_20261022.pdf, 6 pages
(the document's "T25" is the paper's T_26; its "T24" is the paper's T_25).

Run:   MOC4_SUPP=/path/to/supplement  python replay_T26.py [--no-numeric]

Exit code 0 iff every check passes.  EXACT checks use integers, rationals and exact polynomial
arithmetic (sympy, python-flint); lines marked NUMERIC are floating-point falsification attempts.
MOC4_SUPP is needed for: the sha256 comparison of the files that were read, the comparison of
the polynomial G computed here with the stored one, and the sixteen cubic features and the
8 x 8 change of basis of the sum-of-squares certificate behind Lemma 1 and the bounds (6), (7)
(they exist only in JSON data of the supplement).  Without MOC4_SUPP the script FAILS, because
those Gram constructions are then not replayed.
"""
import os
import random
import sys
import time
from fractions import Fraction as Fr

sys.dont_write_bytecode = True          # write nothing, not even bytecode caches
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import replay_T25_common as K          # noqa: E402

NUMERIC = "--no-numeric" not in sys.argv
NPTS = 300 if "--quick" in sys.argv else 4000          # --quick: development only
ck = K.Checker("Replay of case T_26 = {2,3,5,7,8,12}  (document: final_referee_T25_uniform_20261022.pdf)")
K.print_versions(ck)

T = K.PAPER_T[26]
GENS = [0, 1, 2, 5, 6, 39, 67, 70]
GORDER = (0, 1, 2, 5, 6, 39, 67, 70)                  # Gram order of the document: (a, e, p, j, b, c, q, d)
GNAMES = ["a", "e", "p", "j", "b", "c", "q", "d"]
RORDER = (67, 2, 5, 39, 1, 70, 6, 0)                  # order of the radial ring: (q, p, j, c, e, d, b, a)
RNAMES = ["q", "p", "j", "c", "e", "d", "b", "a"]

# ----------------------------------------------------------------------------------------
ck.section("Provenance")
have_supp = K.check_file_identity(ck, ["final_referee_T25_uniform_20261022.pdf",
                                       "final_referee_T25_full_radial_20261022.json",
                                       "final_referee_T25_budget_PSD8_20261022.json",
                                       "final_referee_T25_cubic_reduction_20261022.json"])
ck.check("EXACT", "MOC4_SUPP is set and the four files are the ones that were read (needed for the Gram constructions)",
         have_supp)
ck.info("the document says that a JSON file records the binomial-square certificates of p. 5; no such file "
        "(final_referee_T25_uniform_exact_20261022.json) is in the public supplement; the certificates are found here")

K.check_foundations(ck)
K.check_group_covariance(ck)
K.crosscheck_with_supplement_data(ck)
K.crosscheck_with_audit_library(ck, [25, 26])

# ----------------------------------------------------------------------------------------
K.check_contact_set(ck, 26, (2, 3, 5, 7, 8, 12), "T25", 25)
K.check_generators(ck, T, GENS, [70])
K.check_minor_leibniz(ck, [j for j in GENS if K.F_IS_MINOR[j]])
ck.check("EXACT", "the eight generators have rank 8: the representation (1) g = a f0 + e f1 + p f2 + j f5 + b f6 + c f39 "
         "+ q f67 + d f70 with nonnegative coefficients is unique; d = 0 is a minor decomposition",
         K.rank([K.F_VALUES[j] for j in GENS]) == 8)
stab = [g for g, gm in zip(K.G_ELEMENTS, K.G_MAPS) if K.g_set(gm, T) == tuple(sorted(T))]
ck.info("stabiliser of T_26 in G: %d element(s) (no symmetry is used by the document)" % len(stab))

K.check_radial_necessity(ck)

# ----------------------------------------------------------------------------------------
ck.section("The cubic on g + h kappa, equation (2)  [p. 1]")
t0 = time.time()
rad = K.Radial(RORDER, RNAMES)
R = rad.R
q, p, j, c, e, d, b, a, h = R.gens
Cst = (d + q) * (a + d) * (b + d) + c * d * (a + d + e) + d * e * (d + q)
Q = c * (p - d - e) - d * (d + e) - p * j
ck.check("EXACT", "C(g + h kappa) = -C_* + h Q with C_* = (d+q)(a+d)(b+d) + cd(a+d+e) + de(d+q) and "
         "Q = c(p-d-e) - d(d+e) - pj  (the document's formulas are these at d = 1; C is the paper's cubic (A.2))",
         rad.F == Q * h - Cst)
ck.check("EXACT", "C_* - d^3 has only nonnegative coefficients: C_* >= d^3 > 0, so C(Re x) = 0 forces Q != 0, h = C_*/Q",
         K.coeff_stats(Cst - d ** 3)[1] >= 0)

# ----------------------------------------------------------------------------------------
ck.section("The full 14 x 14 Hessian determinant, equation (3)  [p. 1]")
Ddet = rad.det()
ck.info("det Hess C(g + h kappa): %d terms, %d minors memoised, degrees in (q,p,j,c,e,d,b,a,h) = %s, %.1f s"
        % (len(Ddet), rad.n_minors, R.degrees(Ddet), time.time() - t0))
rng = random.Random(2612)
ok = True
for _ in range(3):
    pt = [Fr(rng.randint(-7, 9), rng.randint(1, 4)) for _ in range(9)]
    u = [R.evaluate(x, pt) for x in rad.u]
    ok &= R.evaluate(Ddet, pt) == K.det_fraction(K.hess_at(u))
ck.check("EXACT", "the determinant polynomial (Laplace expansion, own code) equals the determinant of the rational matrix "
         "Hess C(u) computed by Gaussian elimination at 3 random rational points (seed 2612)", ok)
okf, G, Ssq = rad.radial_factor(Cst, Q, "p")
ck.check("EXACT", "Q^4 det Hess C(g + (C_*/Q) kappa) is divisible by -4 Q^2 and the quotient is the square of "
         "G := gcd(quotient, d quotient / dp): identity (3), Q^4 det = -4 Q^2 G^2 (degree of det in h: %d)" % rad.h_degree, okf)
D9 = (0, 0, 0, 0, 0, 9, 0, 0, 0)
if dict(R.terms(G)).get(D9, 0) < 0:
    G = G * (-1)
tG = R.terms(G)
ck.check("EXACT", "G is homogeneous of degree 9 in the eight coefficients, free of h, of degree 2 in p, with 1201 integer "
         "monomials; the sign is fixed by a positive coefficient of d^9",
         all(sum(ex) == 9 and ex[8] == 0 and ex[1] <= 2 for ex, _ in tG) and len(tG) == 1201 and dict(tG).get(D9, 0) > 0,
         "%d terms, coefficient of d^9 = %s" % (len(tG), dict(tG).get(D9)))
if K.supp_dir() is not None:
    import json
    with open(os.path.join(K.supp_dir(), "final_referee_T25_full_radial_20261022.json")) as fh:
        radj = json.load(fh)
    gt = [f for f in radj["radial_factors"] if f["degree"] == 9][0]["terms"]
    ck.check("EXACT", "supplement data: the stored degree-9 polynomial (variables q,p,j,c,e,d,b,a,z) equals the G computed "
             "here, with the same sign", R.from_terms([(tuple(ex), int(co)) for ex, co in gt]) == G
             and radj["variables"][:8] == RNAMES)
ck.info("consequence: a point of X with coefficients (1), d > 0, has G = 0 (necessary condition only)")

# ----------------------------------------------------------------------------------------
ck.section("Cited capacities (4), (5): SOC, old-H, two four-coefficient inequalities  [p. 1]")
K.orbit_check(ck, "SOC, pj >= 1: (p: f2, j: f5; d: f70)", (1, 9, 22), (2, 5, 70), sym_positions=(0, 1))
K.orbit_check(ck, "old-H, c(p-1) >= 1: (t = c: f39, c-role = p: f2; d: f70)", (12, 2, 70), (39, 2, 70))
K.orbit_check(ck, "four-coefficient, ab(e-1) >= a+b+e+1: (a: f0, b: f6; c-role = e: f1; d: f70)", (0, 6, 10, 54),
              (0, 6, 1, 70), sym_positions=(0, 1))
K.orbit_check(ck, "four-coefficient, bq(c-1) >= b+q+c+1: (b: f6, q: f67; c: f39; d: f70)", (0, 6, 10, 54),
              (6, 67, 39, 70), sym_positions=(0, 1))
ck.info("CITED statements (PROOF.md 7.6): a f1 + b f9 + d f22 has Phi >= 0 on SO(4) if ab >= d^2; t f12 + c f2 + d f70 if "
        "t(c-d) >= d^2; a f0 + b f6 + c f10 + d f54 if ab(c-d) >= d^2(a+b+c+d); each on its whole G-orbit")
k1 = K.find_group_elements((55, 67, 53, 59, 54), (6, 67, 39, 2, 70))
k2 = K.find_group_elements((55, 67, 53, 59, 54), (67, 6, 39, 2, 70))
ck.check("EXACT", "transfer of the T_25 document's family (b: f55, a: f67; e: f53; j: f59; d: f54) to this face: group "
         "elements map it to (b: f6, q: f67; c: f39; p: f2; d: f70) and also to (q: f67, b: f6; ...) with b, q exchanged",
         len(k1) > 0 and len(k2) > 0,
         "e.g. g = %s and g = %s" % (K.G_ELEMENTS[k1[0]], K.G_ELEMENTS[k2[0]]))
src = (0, 2, 5, 6, 39, 67, 70)
hit = [g for g, gm in zip(K.G_ELEMENTS, K.G_MAPS)
       if [K.g_facet(gm, jx) for jx in src] == [6, 59, 2, 55, 53, 67, 54]]
ck.check("EXACT", "the layer e = 0: g = (rho, rho', transpose) = ((0,1,3,2), (2,1,0,3), 1) maps f0, f2, f5, f6, f39, f67, f70 "
         "to f6, f59, f2, f55, f53, f67, f54, i.e. into the canonical generators of the T_25 document with q = 0 "
         "(coefficients (a,p,j,b,c,q,d) -> (c,j,p,b,e,a,d) there); all seven images vanish on T_25",
         ((0, 1, 3, 2), (2, 1, 0, 3), 1) in hit
         and all(K.F_VALUES[jx][s] == 0 for jx in (6, 59, 2, 55, 53, 67, 54) for s in K.PAPER_T[25]))

# ----------------------------------------------------------------------------------------
ck.section("Gram constructions: the bounds (6), (7), the double-SOC endpoint and Lemma 1 (9), (10)  [pp. 2-3]")
if K.supp_dir() is not None:
    import sympy as sp
    gf = K.GramFace(K.supp_dir(), "final_referee_T25")
    ck.check("EXACT", "certificate data (JSON): source order (0,1,2,5,6,39,67,70); the four signed coordinate maps "
             "permute p and q separately", gf.order == list(GORDER) and gf.signed_perms_ok)
    ck.check("EXACT", "quaternion model M(p,q) = matrix of x -> p x q: M^t M = |p|^2 |q|^2 I and det M = (|p|^2 |q|^2)^2 as "
             "polynomial identities (so M(p,q) is in SO(4) for unit p, q)", gf.model_ok)
    key = {k: i for i, k in enumerate(gf.keys)}
    ok = len(gf.relations) == 10
    for m in range(8):
        v = [Fr(0)] * 36
        v[key[(m, m)]] = Fr(1)
        ok &= gf.weight_rows[m] == v
    ck.check("EXACT", "own linear algebra on (*) (see GramFace; %d coefficient equations, rank %d): the solutions are "
             "exactly: 10 independent linear relations on the Gram matrix, and weights = its diagonal "
             "(a, e, p, j, b, c, q, d)" % (gf.n_equations, gf.rank), ok)
    ck.check("EXACT", "the ten relations stored in the supplement span the same space as the ten derived here",
             K.rank(gf.stored_R) == 10 and K.rank(gf.stored_R + gf.relations) == 10)
    ck.info("relations derived here (vectors (A,E,P,J,B,C,Q,D) at positions 0..7):")
    for r in gf.relations:
        ck.info("   " + " ".join("%+d*G%d%d" % (int(cf), gf.keys[i][0], gf.keys[i][1]) for i, cf in enumerate(r) if cf) + " = 0")

    sa, sb, sc, sq, sp_, sj, sx, ss, s2, sal, sbe, st, sz, sw = sp.symbols("a b c q p j x s s2 alpha beta t z w", positive=True)
    # --- (6): B.Q = 1, C.D = C.Q = x, C.B = -x, P.C = 1 - x, |D| = 1
    E2 = (sb + 1) * (sq + 1) / (sb * sq - 1)
    cmin = (1 - sx) ** 2 / sp_ + E2 * sx ** 2
    K.sym_gram_check(ck, gf, "(6), bq > 1, general x, c = (1-x)^2/p + E2 x^2 + s",
                     {(2, 5): 1 - sx, (5, 7): sx, (5, 6): sx, (4, 5): -sx, (4, 6): 1},
                     [0, 0, sp_, 0, sb, cmin + ss, sq, 1], [0, 0, sp_, 0, sb, cmin + ss, sq, 1],
                     [2, 7, 4, 6, 5], [sp_, 1, sb, (sb * sq - 1) / sb, ss],
                     {sb: Fr(3), sq: Fr(2), sp_: Fr(1, 2), sx: Fr(1, 3), ss: Fr(1, 5)},
                     note="  [E2 = (b+1)(q+1)/(bq-1)]")
    ck.check("EXACT", "(6): (1-x)^2/p + E2 x^2 = E2/(1+pE2) + ((1+pE2)/p)(x - 1/(1+pE2))^2: least |C|^2 = 1/(p + 1/E2); "
             "so c >= 1/(p + 1/E2), i.e. p >= 1/c - 1/E2, is safe",
             sp.cancel(cmin - E2 / (1 + sp_ * E2) - (1 + sp_ * E2) / sp_ * (sx - 1 / (1 + sp_ * E2)) ** 2) == 0)
    K.sym_gram_check(ck, gf, "(6) at p = 0 (x = 1): the four-coefficient inequality c >= E2 re-proved on this face",
                     {(5, 7): 1, (5, 6): 1, (4, 5): -1, (4, 6): 1},
                     [0, 0, 0, 0, sb, E2 + ss, sq, 1], [0, 0, 0, 0, sb, E2 + ss, sq, 1],
                     [7, 4, 6, 5], [1, sb, (sb * sq - 1) / sb, ss], {sb: Fr(5, 2), sq: Fr(2), ss: Fr(1, 7)})
    # --- (7): b = beta^2, q = alpha^2, w = alpha beta, B = beta v, Q = alpha v, P = Pv v + C/c + ...
    Pv = (1 - sal * sbe) / (sbe - sal)
    K.sym_gram_check(ck, gf, "(7), b = beta^2, q = alpha^2, w = alpha beta, p = 1/c + (1-w)^2/(beta-alpha)^2 + s",
                     {(2, 5): 1, (2, 4): Pv * sbe, (2, 6): Pv * sal, (4, 6): sal * sbe},
                     [0, 0, 1 / sc + Pv ** 2 + ss, 0, sbe ** 2, sc, sal ** 2, 1],
                     [0, 0, 1 / sc + Pv ** 2 + ss, 0, sbe ** 2, sc, sal ** 2, 1],
                     [4, 5, 2, 7], [sbe ** 2, sc, ss, 1],
                     {sal: Fr(1, 3), sbe: Fr(3, 2), sc: Fr(2), ss: Fr(1, 4)},
                     note="  [vector Q is alpha/beta times vector B: no pivot]")
    ck.check("EXACT", "(7): (1-w)^2/(sqrt b - sqrt q)^2 = b(1-w)^2/(b-w)^2 = q(1-w)^2/(q-w)^2 for w = sqrt(bq): the bound is "
             "symmetric in b, q (the document's 'symmetric version')",
             sp.cancel(Pv ** 2 - sbe ** 2 * (1 - sal * sbe) ** 2 / (sbe ** 2 - sal * sbe) ** 2) == 0
             and sp.cancel(Pv ** 2 - sal ** 2 * (1 - sal * sbe) ** 2 / (sal ** 2 - sal * sbe) ** 2) == 0)
    K.sym_gram_check(ck, gf, "double-SOC endpoint: blocks (P,C) with P.C = 1 and (B,Q) with B.Q = 1, D alone",
                     {(2, 5): 1, (4, 6): 1}, [0, 0, sp_, 0, sb, sc, sq, 1], [0, 0, sp_, 0, sb, sc, sq, 1],
                     [2, 5, 7, 4, 6], [sp_, sc - 1 / sp_, 1, sb, sq - 1 / sb],
                     {sb: Fr(2), sq: Fr(1, 2), sp_: Fr(3), sc: Fr(1, 3)}, note="  [PSD iff pc >= 1 and bq >= 1]")
    K.sym_gram_check(ck, gf, "SOC pj >= 1 re-proved on this face: block (P,J) with P.J = -1, D alone",
                     {(2, 3): -1}, [0, 0, sp_, sj, 0, 0, 0, 1], [0, 0, sp_, sj, 0, 0, 0, 1],
                     [2, 3, 7], [sp_, sj - 1 / sp_, 1], {sp_: Fr(2), sj: Fr(1, 2)})

    # --- Lemma 1.  Base vectors (A, B, Q, D, C) with Gram [[a,t,0],[t,b,z],[0,z,q]] + [1] + [c];
    #     E = al A + be B + t D,  P = ga (B - (t/a) A - Q) + (r/c) C,  J = 0.
    def lemma1_gram(tt, zz, qq_, gamma):
        base = sp.Matrix([[sa, tt, 0, 0, 0], [tt, sb, zz, 0, 0], [0, zz, qq_, 0, 0], [0, 0, 0, 1, 0], [0, 0, 0, 0, sc]])
        den = sa * sb - tt ** 2
        al, be = tt * (sb + tt) / den, -tt * (sa + tt) / den
        rr = 1 - tt
        rows = sp.Matrix([[1, 0, 0, 0, 0],                               # A
                          [al, be, 0, tt, 0],                            # E
                          [-gamma * tt / sa, gamma, -gamma, 0, rr / sc],   # P
                          [0, 0, 0, 0, 0],                               # J
                          [0, 1, 0, 0, 0],                               # B
                          [0, 0, 0, 0, 1],                               # C
                          [0, 0, 1, 0, 0],                               # Q
                          [0, 0, 0, 1, 0]])                              # D
        return (rows * base * rows.T).applyfunc(sp.cancel)

    def run_lemma1(label, Gm, weights, order, pivots, inst, note=""):
        ent = {(i, k): Gm[i, k] for i in range(8) for k in range(i + 1, 8) if Gm[i, k] != 0}
        return K.sym_gram_check(ck, gf, label, ent, [Gm[i, i] for i in range(8)], weights, order, pivots, inst, note)

    hb = sb - st ** 2 / sa
    rr = 1 - st
    e_t = st ** 2 + st ** 2 * (sa + sb + 2 * st) / (sa * sb - st ** 2)
    gam = (rr - sz) / (hb + sq - 2 * sz)
    Gm = lemma1_gram(st, sz, sq, gam)
    run_lemma1("Lemma 1, general z: Gram(A,B,Q) = [[a,t,0],[t,b,z],[0,z,q]], E and P as in the document",
               Gm, [sa, e_t, rr ** 2 / sc + (rr - sz) ** 2 / (hb + sq - 2 * sz), 0, sb, sc, sq, 1],
               [0, 4, 6, 7, 5], [sa, hb, sq - sz ** 2 / hb, 1, sc],
               {sa: Fr(2), sb: Fr(3), sc: Fr(1, 2), sq: Fr(1), st: Fr(1, 3), sz: Fr(1, 4)},
               note="  [E, P are combinations of A, B, Q, D, C: no pivots; h_b = b - t^2/a; PSD iff q h_b >= z^2]")
    ck.check("EXACT", "Lemma 1: E.A = t, E.B = -t, E.D = t, |E|^2 = e(t); |V|^2 = h_b + q - 2z for V = B - (t/a)A - Q; "
             "det Gram(A,B,Q) = a (q h_b - z^2)",
             sp.cancel(Gm[0, 1] - st) == 0 and sp.cancel(Gm[1, 4] + st) == 0 and sp.cancel(Gm[1, 7] - st) == 0
             and sp.cancel(Gm[1, 1] - e_t) == 0
             and sp.cancel(sp.Matrix([[sa, st, 0], [st, sb, sz], [0, sz, sq]]).det() - sa * (sq * hb - sz ** 2)) == 0)
    Gm9 = lemma1_gram(st, rr, sq, 0)
    run_lemma1("Lemma 1 (9): z = r = 1 - t, P = (r/c) C",
               Gm9, [sa, e_t, rr ** 2 / sc, 0, sb, sc, sq, 1],
               [0, 4, 6, 7, 5], [sa, hb, sq - rr ** 2 / hb, 1, sc],
               {sa: Fr(2), sb: Fr(3), sc: Fr(1, 2), sq: Fr(1), st: Fr(1, 3)},
               note="  [PSD iff q h_b >= r^2; then e >= e(t), p >= r^2/c is safe]")
    qpar = sw ** 2 * rr ** 2 / hb
    Gm10 = lemma1_gram(st, sw * rr, qpar, (rr - sw * rr) / (hb + qpar - 2 * sw * rr))
    p10 = rr ** 2 / sc + hb * rr ** 2 * (1 - sw) ** 2 / (hb - sw * rr) ** 2
    run_lemma1("Lemma 1 (10): z = w r, q = w^2 r^2 / h_b (vector Q parallel to the part of B orthogonal to A)",
               Gm10, [sa, e_t, p10, 0, sb, sc, qpar, 1],
               [0, 4, 6, 7, 5], [sa, hb, 0, 1, sc],
               {sa: Fr(2), sb: Fr(3), sc: Fr(1, 2), st: Fr(1, 3), sw: Fr(1, 2)},
               note="  [pivot of Q is 0; needs h_b != w r]")
    ck.check("EXACT", "(10): with z = wr, q = w^2 r^2/h_b: h_b + q - 2z = (h_b - wr)^2/h_b and "
             "r^2/c + (r-z)^2/(h_b+q-2z) = r^2/c + h_b r^2 (1-w)^2/(h_b - wr)^2",
             sp.cancel(hb + qpar - 2 * sw * rr - (hb - sw * rr) ** 2 / hb) == 0
             and sp.cancel(rr ** 2 / sc + (rr - sw * rr) ** 2 / (hb + qpar - 2 * sw * rr) - p10) == 0)
    sup = [2, 5, 7]                                    # Gram matrices supported on (P, C, D), G77 = 1
    cols = [(i, k) for i in sup for k in sup if i <= k and (i, k) != (7, 7)]
    Arows = [[r[key[ck_]] for ck_ in cols] for r in gf.relations]
    brhs = [-r[key[(7, 7)]] for r in gf.relations]
    ck.check("EXACT", "the ten relations have no solution with a Gram matrix supported on the vectors (P, C, D) and "
             "G77 = 1: the old-H capacity c(p-1) >= 1 alone is NOT certified by this Gram representation; it is a "
             "genuine citation", K.rank(Arows) < K.rank([row + [bb] for row, bb in zip(Arows, brhs)]))
    E1x = (sa + 1) * (sb + 1) / (sa * sb - 1)
    K.sym_gram_check(ck, gf, "first inequality of (5) re-proved on this face: E.A = 1, E.B = -1, E.D = 1, A.B = 1 (t = 1), "
                     "e = E1 + s",
                     {(0, 1): 1, (1, 4): -1, (1, 7): 1, (0, 4): 1},
                     [sa, E1x + ss, 0, 0, sb, 0, 0, 1], [sa, E1x + ss, 0, 0, sb, 0, 0, 1],
                     [0, 4, 7, 1], [sa, (sa * sb - 1) / sa, 1, ss], {sa: Fr(5, 2), sb: Fr(2), ss: Fr(1, 7)},
                     note="  [E1 = (a+1)(b+1)/(ab-1)]")
    sV, sZ = sp.symbols("V Z", positive=True)
    ck.check("EXACT", "missing region II: with c = E2 V/(1+V) and t = Z/(1+Z): P0 = 1/c - 1/E2 = 1/(E2 V), "
             "P1 = (1-t)^2/c = (1+V)/((1+Z)^2 E2 V), and P0 - P1 = (Z(Z+2) - V)/((1+Z)^2 E2 V): P1 <= P0 iff V <= Z(Z+2)",
             sp.cancel(1 / (E2 * sV / (1 + sV)) - 1 / E2 - 1 / (E2 * sV)) == 0
             and sp.cancel((1 - sZ / (1 + sZ)) ** 2 / (E2 * sV / (1 + sV)) - (1 + sV) / ((1 + sZ) ** 2 * E2 * sV)) == 0
             and sp.cancel(1 / (E2 * sV) - (1 + sV) / ((1 + sZ) ** 2 * E2 * sV)
                           - (sZ * (sZ + 2) - sV) / ((1 + sZ) ** 2 * E2 * sV)) == 0)
    ck.info("PROOF-READ: pivots >= 0 under a, b, c > 0, 0 <= t < 1, ab > t^2 and q h_b >= z^2; extra norm of E or P is an "
            "extra multiple of f1 or f2 (diagonal entries of E, P occur in no relation); j f5 and any minor vector may be "
            "added; so (9), (10) are safe for e >= e(t) and p >= the displayed bound")

    # --- the function e(t) and the inequalities of p. 3
    de_dt = sp.diff(e_t, st)
    claimed = 2 * st * ((sa * sb - st ** 2) ** 2 + sa * sb * (sa + sb + 3 * st) - st ** 3) / (sa * sb - st ** 2) ** 2
    E1s = (sa + 1) * (sb + 1) / (sa * sb - 1)
    ck.check("EXACT", "e'(t) = 2t{(ab-t^2)^2 + ab(a+b+3t) - t^3}/(ab-t^2)^2; e(0) = 0; e(1) = E1 = (a+1)(b+1)/(ab-1)",
             sp.cancel(de_dt - claimed) == 0 and e_t.subs(st, 0) == 0 and sp.cancel(e_t.subs(st, 1) - E1s) == 0)
    ck.info("PROOF-READ: for ab > 1 and 0 < t < 1: ab(a+b+3t) > 3t > t^3, so e' > 0; e is a continuous increasing bijection "
            "[0,1) -> [0,E1): every e in [0,E1) is e(t) for exactly one t in [0,1)")
    ck.check("EXACT", "q h_b - r^2 = (bq - 1)(1 - t^2/(ab)) + [(1 - t^2/(ab)) - (1-t^2)] + 2t(1-t)  and  "
             "(1 - t^2/(ab)) - (1 - t^2) = t^2 (ab-1)/(ab): so bq >= 1, ab > 1, 0 <= t < 1 give q h_b >= r^2 "
             "(strictly if bq > 1 or t > 0)",
             sp.cancel(sq * hb - rr ** 2 - ((sb * sq - 1) * (1 - st ** 2 / (sa * sb))
                                             + st ** 2 * (sa * sb - 1) / (sa * sb) + 2 * st * (1 - st))) == 0)
    ck.check("EXACT", "h_b - r = (b-1)(1-t^2) + b t^2 (ab-1)/(ab) + t(1-t): so b >= 1, ab > 1, 0 <= t < 1 give h_b >= r, "
             "with equality only for b = 1, t = 0; and h_b - b r^2 = b t^2(ab-1)/(ab) + 2bt(1-t) >= 0 (so r^2/h_b <= 1/b)",
             sp.cancel(hb - rr - ((sb - 1) * (1 - st ** 2) + sb * st ** 2 * (sa * sb - 1) / (sa * sb) + st * (1 - st))) == 0
             and sp.cancel(hb - sb * rr ** 2 - (sb * st ** 2 * (sa * sb - 1) / (sa * sb) + 2 * sb * st * (1 - st))) == 0)
else:
    ck.check("EXACT", "Gram constructions replayed", False, "MOC4_SUPP not available")

# ----------------------------------------------------------------------------------------
ck.section("The constant coefficient k0 = C_* G6, its positivity, and the corner j = c = 0  [p. 3]")
k0h = R.from_terms([(ex, co) for ex, co in tG if ex[1] == 0])
G6, rem = divmod(k0h, Cst)
dg6 = R.degrees(G6)
ck.check("EXACT", "k0 = G|_{p=0} = C_* G6 (exact division, d general); G6 has degree 6, %d terms, and degree <= 1 in each "
         "of e, c, j" % len(G6), rem == R.zero and all(sum(ex) == 6 for ex, _ in R.terms(G6))
         and dg6[4] <= 1 and dg6[3] <= 1 and dg6[2] <= 1, "degrees in (q,p,j,c,e,d,b,a,h): %s" % dg6)
# d = 1, variables (a, b, q, e, c, j)
R6 = K.Ring(["a", "b", "q", "e", "c", "j"])
g6 = R6.from_terms([((ex[7], ex[6], ex[0], ex[4], ex[3], ex[2]), co) for ex, co in R.terms(G6)])
a6, b6, q6, e6, c6, j6 = R6.gens
ck.check("EXACT", "G6|_{e=c=j=0} = (a+1)^2 (b+1)^2 (q+1)^2  (d = 1)",
         R6.from_terms([(ex, co) for ex, co in R6.terms(g6) if ex[3] == 0 and ex[4] == 0 and ex[5] == 0])
         == (a6 + 1) ** 2 * (b6 + 1) ** 2 * (q6 + 1) ** 2)
ck.check("EXACT", "G6|_{b=0} = (e+1+a)(q+c+1)(qj + q + qa + j + ja + 1 + a)  (d = 1)",
         R6.from_terms([(ex, co) for ex, co in R6.terms(g6) if ex[1] == 0])
         == (e6 + 1 + a6) * (q6 + c6 + 1) * (q6 * j6 + q6 + q6 * a6 + j6 + j6 * a6 + 1 + a6))
RO6 = K.Ring(["B", "X", "Y", "E", "C", "J"])
B6, X6, Y6, E6, C6, J6 = RO6.gens
one6 = RO6.one
expect6 = {(0, 0): 269, (0, 1): 337, (1, 0): 337, (1, 1): 456}
for hi1 in (0, 1):
    for hi2 in (0, 1):
        if hi1:
            an, ad = one6 + X6, B6
            en, ed = (one6 + X6 + B6) * (B6 + one6) * E6, B6 * X6 * (one6 + E6)
        else:
            an, ad = X6, B6 * (one6 + X6)
            en, ed = E6, one6
        if hi2:
            qn, qd = one6 + Y6, B6
            cn, cd = (one6 + Y6 + B6) * (B6 + one6) * C6, B6 * Y6 * (one6 + C6)
        else:
            qn, qd = Y6, B6 * (one6 + Y6)
            cn, cd = C6, one6
        nums6, dens6 = [an, B6, qn, en, cn, J6], [ad, one6, qd, ed, cd, one6]
        F6 = K.subst_clear(R6.terms(g6), nums6, dens6, RO6)
        st6 = K.coeff_stats(F6)
        # the part free of E, C, J and the base term (a+1)^2 (b+1)^2 (q+1)^2 after the same clearing
        base = K.subst_clear(R6.terms((a6 + 1) ** 2 * (b6 + 1) ** 2 * (q6 + 1) ** 2), nums6, dens6, RO6, R6.degrees(g6))
        F6const = RO6.from_terms([(ex, co) for ex, co in RO6.terms(F6) if ex[3] == 0 and ex[4] == 0 and ex[5] == 0])
        base0 = RO6.from_terms([(ex, co) for ex, co in RO6.terms(base) if ex[3] == 0 and ex[4] == 0 and ex[5] == 0])
        lowE = min(ex[3] for ex, _ in RO6.terms(base))
        lowC = min(ex[4] for ex, _ in RO6.terms(base))
        basefac = (one6 + E6) ** (1 if hi1 else 0) * (one6 + C6) ** (1 if hi2 else 0)
        ck.check("EXACT", "G6 table (ab %s 1, bq %s 1), chart (14): after clearing positive denominators all coefficients "
                 "are >= 0 (%d terms, as in the document); the part free of E, C, J equals the cleared base term "
                 "(a+1)^2(b+1)^2(q+1)^2 divided by %s: G6 >= (a+1)^2(b+1)^2(q+1)^2 (1 - Le e)(1 - Lc c) > 0 on this chart"
                 % (">" if hi1 else "<", ">" if hi2 else "<", st6[0],
                    "(1+E)(1+C)" if hi1 and hi2 else "(1+E)" if hi1 else "(1+C)" if hi2 else "1"),
                 st6[1] >= 0 and st6[0] == expect6[(hi1, hi2)] and F6const * basefac == base and lowE == 0 and lowC == 0,
                 "min coefficient %d" % st6[1])
ck.info("PROOF-READ: 1 - Le e = 1/(1+E) for e = E1 E/(1+E) (ab > 1) and Le = 0 for ab <= 1; the bound is continuous in "
        "(a,b,q,e,c,j), so it extends to ab = 1 and bq = 1; with the b = 0 factorisation: k0 > 0 on the closed orthant "
        "under strict failure of (5)")
# corner (15)
RK = K.Ring(["a", "b", "q", "e", "p"])
ak, bk, qk, ek, pk = RK.gens
corner = RK.from_terms([((ex[7], ex[6], ex[0], ex[4], ex[1]), co) for ex, co in tG if ex[2] == 0 and ex[3] == 0])
Lk = (ak + 1) * (bk + 1)
dk = Lk - ek * (ak * bk - 1)
ck.check("EXACT", "identity (15): G|_{j=c=0, d=1} = (a+1)(e+L)(q+1)^2 { (e+L)[pq + (p+q)(b+1)] + (b+1) delta }, "
         "L = (a+1)(b+1), delta = L - e(ab-1); it is > 0 when delta > 0 (strict failure of the first inequality of (5)): "
         "no point of X with d > 0 has j = c = 0",
         corner == (ak + 1) * (ek + Lk) * (qk + 1) ** 2 * ((ek + Lk) * (pk * qk + (pk + qk) * (bk + 1)) + (bk + 1) * dk))

# ----------------------------------------------------------------------------------------
ck.section("Bernstein criterion; the forty baseline tables  [p. 4]")
ck.check("EXACT", "for G = k0 + k1 x + k2 x^2 and x = t N/D: D^2 G = (1-t)^2 D^2 k0 + t(1-t) D (2D k0 + N k1) "
         "+ t^2 (D^2 k0 + N D k1 + N^2 k2); so D > 0, k0 > 0 and both numerators >= 0 give G > 0 on [0, N/D)",
         K.bernstein_identity_check())
RC = K.Ring(["q", "b", "c", "j", "a", "e"])           # core order (x, y, z, s, v, r) = (q, b, c, j, a, e), d = 1
cq, cb, cc, cj, ca, ce = RC.gens
kk = [RC.from_terms([((ex[0], ex[6], ex[3], ex[2], ex[7], ex[4]), co) for ex, co in tG if ex[1] == n]) for n in range(3)]
DEGS = [max(RC.degrees(x_)[i] for x_ in kk) for i in range(6)]
ck.info("d = 1: G = k2 p^2 + k1 p + k0 with k0, k1, k2 in Z[q,b,c,j,a,e] of %s terms" % [len(x) for x in kk])
RO = K.Ring(["A", "T", "V", "U", "R", "E"])
A_, T_, V_, U_, Rr, Ee = RO.gens
one = RO.one
oneC = RC.one
LN = (A_ + one) * (A_ + one + T_) * V_                # z = L V/(1+V) for x = A, y = (1+T)/A
LD = A_ * T_ * (one + V_)
Fw = (one + T_) * (one + V_) - V_


def base_chart(name):
    """(nums, dens) of (x, y, z, s) and the threshold P* = N/D of table (17)."""
    if name == "I-":
        return [A_, T_, V_, V_ + U_], [one, A_ * (one + T_), one, one], ("core", lambda x, y, z, s: (oneC, s))
    if name == "I+":
        return [A_, one + T_, LN, LN + U_ * LD], [one, A_, LD, LD], ("core", lambda x, y, z, s: (oneC, s))
    if name == "II":
        return [A_, one + T_, LN, LN * U_], [one, A_, LD, LD * (one + U_)], \
            ("core", lambda x, y, z, s: ((x + 1) * (y + 1) - z * (x * y - 1), z * (x + 1) * (y + 1)))
    if name == "IIIa-":
        return [T_ * V_, T_, A_, A_ * U_], [(one + T_) * (one + V_), one + T_, one, (one + A_) * (one + U_)], \
            ("core", lambda x, y, z, s: (z + 1, z))
    if name == "IIIa+":
        return [T_ * V_, T_, A_, A_ * (one + U_) + A_ * A_ * U_], \
            [(one + T_) * (one + V_), one + T_, one, (one + A_) * (one + U_)], ("core", lambda x, y, z, s: (oneC, s))
    if name == "IIIb+":
        return [V_, one + T_, A_, A_ * (one + T_) * (one + U_) + A_ * A_ * U_], \
            [(one + T_) * (one + V_), one, one, (A_ + one + T_) * (one + U_)], ("core", lambda x, y, z, s: (oneC, s))
    if name == "IIIb-":
        return [V_ * V_, one + T_, A_, A_ * (one + T_) * U_], \
            [(one + V_) ** 2 * (one + T_), one, one, (A_ + one + T_) * (one + U_)], \
            ("out", (Fw * Fw + A_ * (one + T_), A_ * Fw * Fw))
    raise ValueError(name)


def base_point(name, swap, high, A, Tt, V, U, Rv, Ev):
    """Table (17) and the split of ab as plain rational arithmetic: returns (q, b, c, j, a, e), P*, region test."""
    if name in ("I-", "I+", "II"):
        x = A
        y = Tt / (A * (1 + Tt)) if name == "I-" else (1 + Tt) / A
        if name == "I-":
            z, s = V, V + U
            Pst, reg = 1 / s, (s >= z and x * y < 1)
        else:
            L_ = (x + 1) * (y + 1) / (x * y - 1)
            z = L_ * V / (1 + V)
            if name == "I+":
                s = z + U
                Pst, reg = 1 / s, (s >= z and x * y > 1 and z < L_)
            else:
                s = z * U / (1 + U)
                Pst, reg = 1 / z - 1 / L_, (z >= s and x * y > 1 and z < L_)
    elif name in ("IIIa-", "IIIa+"):
        y = Tt / (1 + Tt)
        x = y * V / (1 + V)
        z = A
        if name == "IIIa-":
            s = A * U / ((1 + A) * (1 + U))
            Pst, reg = 1 + 1 / z, (x <= y <= 1 and s <= z / (1 + z))
        else:
            s = A / (1 + A) + A * A * U / ((1 + A) * (1 + U))
            Pst, reg = 1 / s, (x <= y <= 1 and z / (1 + z) <= s <= z)
    else:
        y = 1 + Tt
        z = A
        if name == "IIIb+":
            x = V / (y * (1 + V))
            s = A * y / (A + y) + A * A * U / ((A + y) * (1 + U))
            Pst, reg = 1 / s, (y >= 1 and x * y < 1 and z * y / (z + y) <= s <= z)
        else:
            w_ = V / (1 + V)
            x = w_ * w_ / y
            s = A * y * U / ((A + y) * (1 + U))
            Pst, reg = 1 / z + y * (1 - w_) ** 2 / (y - w_) ** 2, (y >= 1 and x * y < 1 and s <= z * y / (z + y))
    q_, b_ = (x, y) if swap == 0 else (y, x)
    if high:
        a_ = (1 + Rv) / b_
        E1_ = (a_ + 1) * (b_ + 1) / (a_ * b_ - 1)
        e_ = E1_ * Ev / (1 + Ev)
        reg = reg and a_ * b_ > 1 and e_ < E1_
    else:
        a_ = Rv / (b_ * (1 + Rv))
        e_ = Ev
        reg = reg and a_ * b_ < 1
    return (q_, b_, z, s, a_, e_), Pst, reg


def validate(label, Fm, Fe, ring_out, nums, dens, ND, kind, point_fn, nparams, seed, dm, de):
    """Exact evaluation at random positive rational parameters: the plain chart map lands in its region, N/D is the
    stated threshold, and the table polynomials are the Bernstein numerators times positive denominator powers."""
    rg = random.Random(seed)
    ok_ = True
    for _ in range(3):
        prm = [Fr(rg.randint(1, 9), rg.randint(1, 5)) for _ in range(nparams)]
        pt, Pst, reg = point_fn(prm)
        kv = [RC.evaluate(x_, list(pt)) for x_ in kk]
        dv = [ring_out.evaluate(x_, prm) for x_ in dens]
        if kind == "core":
            Nv, Dv = RC.evaluate(ND[0], list(pt)), RC.evaluate(ND[1], list(pt))
        else:
            Nv, Dv = ring_out.evaluate(ND[0], prm), ring_out.evaluate(ND[1], prm)
        scm = sce = Fr(1)
        for i in range(6):
            scm *= dv[i] ** dm[i]
            sce *= dv[i] ** de[i]
        ok_ &= reg and all(x_ > 0 for x_ in dv) and Dv > 0 and Nv > 0 and Nv / Dv == Pst and kv[0] > 0
        ok_ &= ring_out.evaluate(Fm, prm) == (2 * Dv * kv[0] + Nv * kv[1]) * scm
        ok_ &= ring_out.evaluate(Fe, prm) == (Dv * Dv * kv[0] + Nv * Dv * kv[1] + Nv * Nv * kv[2]) * sce
    return ok_


SEL = [(n, 0, hi) for n in ("I-", "I+") for hi in (0, 1)] + [("II", 0, 0)] \
    + [(n, sw, hi) for n in ("IIIa-", "IIIa+", "IIIb+") for sw in (0, 1) for hi in (0, 1)] \
    + [("IIIb-", 0, 0), ("IIIb-", 1, 0), ("IIIb-", 1, 1)]


def baseline_tables(name, swap, high):
    n4, d4, thr = base_chart(name)
    if swap == 0:
        nq, dq, nb, db = n4[0], d4[0], n4[1], d4[1]
    else:
        nq, dq, nb, db = n4[1], d4[1], n4[0], d4[0]
    if high:                                           # ab > 1: a = (1+R)/b, e = E1 E/(1+E)
        na, da = (one + Rr) * db, nb
        ne, de_ = ((one + Rr) * db + nb) * (nb + db) * Ee, nb * db * Rr * (one + Ee)
    else:                                              # ab < 1: a = R/(b(1+R)), e = E
        na, da = Rr * db, nb * (one + Rr)
        ne, de_ = Ee, one
    nums = [nq, nb, n4[2], n4[3], na, ne]
    dens = [dq, db, d4[2], d4[3], da, de_]
    if thr[0] == "core":
        x, y = (cq, cb) if swap == 0 else (cb, cq)
        ND = thr[1](x, y, cc, cj)
        mid, end = K.bernstein_polys(RC, kk, ND[0], ND[1])
        Fm = K.subst_clear(RC.terms(mid), nums, dens, RO)
        Fe = K.subst_clear(RC.terms(end), nums, dens, RO)
        dm, de2 = RC.degrees(mid), RC.degrees(end)
    else:
        ND = thr[1]
        sk = [K.subst_clear(RC.terms(x_), nums, dens, RO, DEGS) for x_ in kk]
        Fm, Fe = K.bernstein_polys(RO, sk, ND[0], ND[1])
        dm = de2 = DEGS
    return Fm, Fe, nums, dens, ND, thr[0], dm, de2


n_tab = 0
all_ok = True
for idx, (name, swap, high) in enumerate(SEL):
    t0 = time.time()
    Fm, Fe, nums, dens, ND, kind, dm, de2 = baseline_tables(name, swap, high)
    sm, se_ = K.coeff_stats(Fm), K.coeff_stats(Fe)
    okv = validate(name, Fm, Fe, RO, nums, dens, ND, kind,
                   lambda prm, name=name, swap=swap, high=high: base_point(name, swap, high, *prm), 6, 2620 + idx, dm, de2)
    n_tab += 2
    all_ok &= ck.check("EXACT", "baseline chart %-5s, (x,y) = %s, ab %s 1: middle and endpoint numerators have only "
                       "nonnegative coefficients (%d / %d terms); the encoding agrees with the plain chart map, region and "
                       "threshold at 3 random rational parameter points"
                       % (name, "(q,b)" if swap == 0 else "(b,q)", ">" if high else "<", sm[0], se_[0]),
                       sm[1] >= 0 and se_[1] >= 0 and okv,
                       "min coefficients %d / %d, %.1f s" % (sm[1], se_[1], time.time() - t0))
ck.check("EXACT", "all forty baseline tables (twenty charts) are nonnegative", all_ok and n_tab == 40)
neg = {}
for name, swap, high in (("II", 0, 1), ("IIIb-", 0, 1)):
    Fm, Fe = baseline_tables(name, swap, high)[:2]
    neg[name] = (K.coeff_stats(Fm)[2], K.coeff_stats(Fe)[2])
ck.info("for the record: with the baseline thresholds the two charts the document leaves out have negative coefficients: "
        "II with ab > 1: %d / %d; IIIb- with ab > 1 and larger core variable b: %d / %d (middle / endpoint)"
        % (neg["II"] + neg["IIIb-"]))

# ----------------------------------------------------------------------------------------
ck.section("Missing region II (ab > 1, bq > 1, c >= j): charts (18) and the binomial squares (19)  [p. 5]")


def e_of_t(a_, b_, t_):
    return t_ * t_ + t_ * t_ * (a_ + b_ + 2 * t_) / (a_ * b_ - t_ * t_)


def missII_point(mode, prm):
    B, X, Y, Z, Wv, U = prm
    b_, a_, q_ = B, (1 + X) / B, (1 + Y) / B
    t_ = Z / (1 + Z)
    e_ = e_of_t(a_, b_, t_)
    E2_ = (q_ + 1) * (b_ + 1) / (b_ * q_ - 1)
    V = Z * (Z + 2) * Wv / (1 + Wv) if mode == 0 else Z * (Z + 2) + Wv
    c_ = E2_ * V / (1 + V)
    j_ = c_ * U / (1 + U)
    P0, P1 = 1 / c_ - 1 / E2_, (1 - t_) ** 2 / c_
    reg = a_ * b_ > 1 and b_ * q_ > 1 and c_ >= j_ and c_ < E2_ and e_ < (a_ + 1) * (b_ + 1) / (a_ * b_ - 1) \
        and ((P1 <= P0) if mode == 0 else (P0 <= P1))
    return (q_, b_, c_, j_, a_, e_), (P1 if mode == 0 else P0), reg


R2 = K.Ring(["B", "X", "Y", "Z", "W", "U"])
B2, X2, Y2, Z2, W2, U2 = R2.gens
o2 = R2.one
Dt2 = (o2 + X2) * (o2 + Z2) ** 2 - Z2 * Z2
Et2 = ((o2 + X2) * B2 + (o2 + X2) + B2 * B2) * (o2 + Z2) ** 2 + B2 * (Z2 * 2 + Z2 * Z2)
eN2, eD2 = Z2 * Z2 * Et2, B2 * (o2 + Z2) ** 2 * Dt2
EXPECT2 = {0: (102873, 163889, 2), 1: (115308, 184864, 6)}


def binomial_squares(poly, ring):
    """For twice the polynomial: every negative coefficient -k at exponent m is paired with two monomials m + delta,
    m - delta (delta in {-1,0,1}^n, nonzero) whose remaining coefficients are >= k/2; k/2 is subtracted from each.
    (k/2)(X^(m+delta) + X^(m-delta)) - k X^m = (k/2) X^(m-|delta|) (X^(delta+) - X^(delta-))^2 >= 0 for X >= 0.
    Returns (certificate list, True iff the residual has only nonnegative coefficients)."""
    import itertools
    dd = {ex: 2 * co for ex, co in ring.terms(poly)}
    cert = []
    n = ring.n
    deltas = sorted((d_ for d_ in itertools.product((-1, 0, 1), repeat=n) if any(d_)),
                    key=lambda d_: (sum(abs(x_) for x_ in d_), d_))
    for ex in sorted(m for m, co in dd.items() if co < 0):
        k = -dd[ex]
        found = None
        for d_ in deltas:
            lo = tuple(x_ - y_ for x_, y_ in zip(ex, d_))
            hi = tuple(x_ + y_ for x_, y_ in zip(ex, d_))
            if min(lo) < 0 or min(hi) < 0:
                continue
            if k % 2 == 0 and dd.get(lo, 0) >= k // 2 and dd.get(hi, 0) >= k // 2:
                found = (d_, lo, hi)
                break
        if found is None:
            return cert, False
        d_, lo, hi = found
        dd[lo] -= k // 2
        dd[hi] -= k // 2
        dd[ex] = 0
        cert.append((ex, d_, k))
    return cert, min(dd.values()) >= 0


for mode in (0, 1):
    t0 = time.time()
    if mode == 0:
        vn, vd = Z2 * (Z2 + 2) * W2, o2 + W2            # V = Z(Z+2) W/(1+W), threshold P1 = (1-t)^2/c
    else:
        vn, vd = Z2 * (Z2 + 2) + W2, o2                 # V = Z(Z+2) + W,     threshold P0 = 1/c - 1/E2
    cn = (o2 + Y2 + B2) * (B2 + o2) * vn                # c = E2 V/(1+V), E2 = (1+Y+B)(B+1)/(BY)
    cd = B2 * Y2 * (vd + vn)
    nums = [o2 + Y2, B2, cn, cn * U2, o2 + X2, eN2]     # (q, b, c, j, a, e)
    dens = [B2, o2, cd, cd * (o2 + U2), B2, eD2]
    ND = (cd, cn * (o2 + Z2) ** 2) if mode == 0 else (B2 * Y2 * vd, (o2 + Y2 + B2) * (B2 + o2) * vn)
    sk = [K.subst_clear(RC.terms(x_), nums, dens, R2, DEGS) for x_ in kk]
    Fm, Fe = K.bernstein_polys(R2, sk, ND[0], ND[1])
    sm, se_ = K.coeff_stats(Fm), K.coeff_stats(Fe)
    okv = validate("missing II", Fm, Fe, R2, nums, dens, ND, "out", lambda prm, mode=mode: missII_point(mode, prm),
                   6, 2650 + mode, DEGS, DEGS)
    ck.check("EXACT", "missing II, %s: the middle numerator has only nonnegative coefficients; term counts %d / %d and %d "
             "negative endpoint coefficients, as in the document; encoding agrees with the plain chart map at 3 points"
             % ("V = Z(Z+2)W/(1+W), threshold P1 = (1-t)^2/c" if mode == 0 else "V = Z(Z+2) + W, threshold P0 = 1/c - 1/E2",
                sm[0], se_[0], se_[2]),
             sm[1] >= 0 and (sm[0], se_[0], se_[2]) == EXPECT2[mode] and okv,
             "min coefficients %d / %d, %.1f s" % (sm[1], se_[1], time.time() - t0))
    cert, okc = binomial_squares(Fe, R2)
    ck.check("EXACT", "missing II, endpoint numerator of this chart: twice the polynomial minus %d binomial squares "
             "(k/2) X^(m-|delta|) (X^(delta+) - X^(delta-))^2 has only nonnegative coefficients (certificate found here; "
             "k in %s): the endpoint numerator is >= 0 on the closed orthant"
             % (len(cert), sorted(set(k for _, _, k in cert))), okc and len(cert) == EXPECT2[mode][2])
    for ex, d_, k in cert:
        ck.info("   centre exponent (B,X,Y,Z,W,U) = %s, coefficient of twice the polynomial -%d, delta = %s" % (ex, k, d_))

# ----------------------------------------------------------------------------------------
ck.section("Missing region IIIb- (b >= 1, ab > 1, bq < 1, j <= cb/(c+b)): charts (20), (21)  [p. 5]")


def missIII_point(mode, prm):
    B, X, Z, Wv, C, U = prm
    b_ = 1 + B
    a_ = (1 + X) / b_
    t_ = Z / (1 + Z)
    e_ = e_of_t(a_, b_, t_)
    c_ = C
    j_ = C * b_ * U / ((C + b_) * (1 + U))
    hb_ = b_ - t_ * t_ / a_
    r_ = 1 - t_
    Dv = (1 + X) * (1 + Z) ** 2 - Z * Z
    if mode == 0:
        q_ = (1 + X + Dv * Wv) / (b_ * Dv * (1 + Wv))
        Pst = r_ * r_ / c_
        reg = q_ * hb_ >= r_ * r_
    else:
        w_ = Wv / (1 + Wv)
        q_ = w_ * w_ * r_ * r_ / hb_
        Pst = r_ * r_ / c_ + hb_ * r_ * r_ * (1 - w_) ** 2 / (hb_ - w_ * r_) ** 2
        reg = q_ * hb_ <= r_ * r_ and hb_ > w_ * r_ and q_ == Wv * Wv * (1 + X) / ((1 + Wv) ** 2 * b_ * Dv)
    reg = reg and b_ >= 1 and a_ * b_ > 1 and b_ * q_ < 1 and j_ <= c_ * b_ / (c_ + b_) \
        and e_ < (a_ + 1) * (b_ + 1) / (a_ * b_ - 1) and hb_ == b_ * Dv / ((1 + X) * (1 + Z) ** 2) and r_ == 1 / (1 + Z)
    return (q_, b_, c_, j_, a_, e_), Pst, reg


R3 = K.Ring(["B", "X", "Z", "W", "C", "U"])
B3, X3, Z3, W3, C3, U3 = R3.gens
o3 = R3.one
by = o3 + B3
ax = o3 + X3
Dt3 = ax * (o3 + Z3) ** 2 - Z3 * Z3
Et3 = (ax * by + ax + by * by) * (o3 + Z3) ** 2 + by * (Z3 * 2 + Z3 * Z3)
eN3, eD3 = Z3 * Z3 * Et3, by * (o3 + Z3) ** 2 * Dt3
jN3, jD3 = C3 * by * U3, (C3 + by) * (o3 + U3)
EXPECT3 = {0: (131291, 147915), 1: (507835, 1066354)}
for mode in (0, 1):
    t0 = time.time()
    if mode == 0:
        qn, qd = ax + Dt3 * W3, (o3 + W3) * by * Dt3                     # (20): q h_b >= r^2
        ND = (o3, C3 * (o3 + Z3) ** 2)                                   # threshold (9): r^2/c
    else:
        qn, qd = W3 * W3 * ax, (o3 + W3) ** 2 * by * Dt3                 # (21): q = w^2 r^2/h_b, w = W/(1+W)
        Ft = by * Dt3 * (o3 + W3) - W3 * ax * (o3 + Z3)
        ND = (Ft * Ft + C3 * (o3 + Z3) ** 2 * by * Dt3 * ax, C3 * (o3 + Z3) ** 2 * Ft * Ft)      # threshold (10)
    nums = [qn, by, C3, jN3, ax, eN3]
    dens = [qd, o3, o3, jD3, by, eD3]
    sk = [K.subst_clear(RC.terms(x_), nums, dens, R3, DEGS) for x_ in kk]
    Fm, Fe = K.bernstein_polys(R3, sk, ND[0], ND[1])
    sm, se_ = K.coeff_stats(Fm), K.coeff_stats(Fe)
    okv = validate("missing IIIb-", Fm, Fe, R3, nums, dens, ND, "out", lambda prm, mode=mode: missIII_point(mode, prm),
                   6, 2660 + mode, DEGS, DEGS)
    ck.check("EXACT", "missing IIIb-, %s: middle and endpoint numerators have only nonnegative coefficients; term counts "
             "%d / %d as in the document; encoding agrees with the plain chart map at 3 points"
             % ("chart (20), threshold (9)" if mode == 0 else "chart (21), threshold (10)", sm[0], se_[0]),
             sm[1] >= 0 and se_[1] >= 0 and (sm[0], se_[0]) == EXPECT3[mode] and okv,
             "min coefficients %d / %d, %.1f s" % (sm[1], se_[1], time.time() - t0))
    if mode == 1:
        # the key factor b D (1+W) - W (1+X)(1+Z) of the denominator is positive on the chart
        Ftx = K.subst_clear(R3.terms(Ft), [B3, X3, Z3, W3, C3, U3], [o3] * 6, R3)
        ck.check("EXACT", "chart (21): the key factor b D (1+W) - W(1+X)(1+Z), b = 1 + B, D = (1+X)(1+Z)^2 - Z^2, expands to a "
                 "polynomial in B, X, Z, W with nonnegative coefficients and constant term 1: it is > 0",
                 K.coeff_stats(Ftx)[1] >= 0 and dict(R3.terms(Ftx)).get((0, 0, 0, 0, 0, 0), 0) == 1)

# ----------------------------------------------------------------------------------------
ck.section("Cover of the coefficient domain, boundaries included: exact randomised test in the ORIGINAL coefficients  [p. 6]")
from math import isqrt       # noqa: E402


def exact_sqrt(x):
    n_, d_ = x.numerator, x.denominator
    rn, rd = isqrt(n_), isqrt(d_)
    return Fr(rn, rd) if rn * rn == n_ and rd * rd == d_ else None


def rnd(rg, zero_prob=0.0):
    if rg.random() < zero_prob:
        return Fr(0)
    return Fr(rg.randint(1, 40), rg.randint(1, 40)) * Fr(10) ** rg.randint(-1, 1)


def domain_sample(rg):
    """A random rational point (q, b, c, j, a, e, t) of the domain {d = 1, both inequalities of (5) fail strictly,
    (c, j) != (0, 0)}, built so that every square root the thresholds need is rational.  t is the parameter with
    e = e(t) when ab > 1 (None otherwise)."""
    kind = rg.randint(0, 9)
    b_ = rnd(rg, 0.1)
    if kind == 0:
        b_ = Fr(1)
    # q: choose w = sqrt(bq) rational
    if b_ > 0:
        if kind == 1:
            q_ = 1 / b_                                # bq = 1
        elif rg.random() < 0.55:
            w_ = Fr(rg.randint(0, 99), 100)            # bq = w^2 < 1
            q_ = w_ * w_ / b_
        else:
            q_ = rnd(rg)
    else:
        q_ = rnd(rg, 0.2)
    # a and e
    t_ = None
    if b_ > 0 and rg.random() < 0.5:
        a_ = (1 + rnd(rg, 0.15)) / b_                  # ab >= 1
    else:
        a_ = rnd(rg, 0.2)
        if b_ > 0 and a_ * b_ > 1:
            a_ = 1 / (b_ * (1 + a_))
    if a_ * b_ > 1:
        t_ = Fr(rg.randint(0, 99), 100)
        e_ = e_of_t(a_, b_, t_)
    else:
        e_ = rnd(rg, 0.2)
    # c
    c_ = rnd(rg, 0.12)
    if b_ * q_ > 1:
        E2_ = (q_ + 1) * (b_ + 1) / (b_ * q_ - 1)
        c_ = E2_ * Fr(rg.randint(0, 99), 100)
    # j relative to c
    jm = rg.randint(0, 6)
    y_ = max(b_, q_)
    if jm == 0:
        j_ = c_
    elif jm == 1:
        j_ = Fr(0)
    elif jm == 2:
        j_ = c_ / (1 + c_)
    elif jm == 3:
        j_ = c_ * y_ / (c_ + y_) if c_ + y_ > 0 else Fr(0)
    elif jm == 4:
        j_ = c_ * Fr(rg.randint(0, 100), 100)
    else:
        j_ = rnd(rg)
    if c_ == 0 and j_ == 0:
        j_ = rnd(rg)
    return q_, b_, c_, j_, a_, e_, t_


def doc_threshold(q_, b_, c_, j_, a_, e_, t_):
    """Region of the document's cover and its threshold P* = N/D in the original coefficients (None if a needed square
    root is irrational: then the point is skipped)."""
    if j_ >= c_:
        return "I", Fr(1), j_
    bq = b_ * q_
    if bq >= 1:
        N0, D0 = (q_ + 1) * (b_ + 1) - c_ * (bq - 1), c_ * (q_ + 1) * (b_ + 1)
        if a_ * b_ <= 1:
            return "II, ab <= 1 (P0)", N0, D0
        r_ = 1 - t_
        if r_ * r_ * D0 <= N0 * c_:                    # P1 = r^2/c <= P0 = N0/D0
            return "II, ab > 1 (P1)", r_ * r_, c_
        return "II, ab > 1 (P0)", N0, D0
    y_ = max(b_, q_)
    if y_ <= 1:
        if j_ <= c_ / (1 + c_):
            return "IIIa-", c_ + 1, c_
        return "IIIa+", Fr(1), j_
    if j_ >= c_ * y_ / (c_ + y_):
        return "IIIb+", Fr(1), j_
    if a_ * b_ <= 1 or q_ > b_:
        w_ = exact_sqrt(bq)
        if w_ is None:
            return None
        return "IIIb- baseline", (y_ - w_) ** 2 + c_ * y_ * (1 - w_) ** 2, c_ * (y_ - w_) ** 2
    r_ = 1 - t_
    hb_ = b_ - t_ * t_ / a_
    if q_ * hb_ >= r_ * r_:
        return "IIIb- new (9)", r_ * r_, c_
    w_ = exact_sqrt(q_ * hb_ / (r_ * r_))
    if w_ is None:
        return None
    return "IIIb- new (10)", r_ * r_ * (hb_ - w_ * r_) ** 2 + c_ * hb_ * r_ * r_ * (1 - w_) ** 2, c_ * (hb_ - w_ * r_) ** 2


rng = random.Random(2614)
count = {}
ok = True
n_g = n_skip = 0
it = 0
while sum(count.values()) < NPTS:
    it += 1
    q_, b_, c_, j_, a_, e_, t_ = domain_sample(rng)
    if it % 4 == 0:                                    # force the new IIIb- charts: b >= 1, ab > 1, bq < 1
        b_ = 1 + rnd(rng, 0.2)
        a_ = (1 + rnd(rng)) / b_
        t_ = Fr(rng.randint(0, 99), 100)
        e_ = e_of_t(a_, b_, t_)
        hb_ = b_ - t_ * t_ / a_
        r_ = 1 - t_
        if rng.random() < 0.5:
            w_ = Fr(rng.randint(0, 100), 100)
            q_ = w_ * w_ * r_ * r_ / hb_
        else:
            lam = Fr(rng.randint(0, 99), 100)
            q_ = (1 - lam) * r_ * r_ / hb_ + lam / b_
        c_ = rnd(rng)
        j_ = c_ * b_ / (c_ + b_) * Fr(rng.randint(0, 100), 100)
    assert (a_ * b_ <= 1 or e_ * (a_ * b_ - 1) < (a_ + 1) * (b_ + 1)) and (b_ * q_ <= 1 or c_ * (b_ * q_ - 1) < (b_ + 1) * (q_ + 1))
    assert c_ > 0 or j_ > 0
    th = doc_threshold(q_, b_, c_, j_, a_, e_, t_)
    if th is None:
        n_skip += 1
        continue
    reg, N_, D_ = th
    pt = [q_, b_, c_, j_, a_, e_]
    kv = [RC.evaluate(x_, pt) for x_ in kk]
    good = D_ > 0 and N_ > 0 and kv[0] > 0 and 2 * D_ * kv[0] + N_ * kv[1] >= 0 \
        and D_ * D_ * kv[0] + N_ * D_ * kv[1] + N_ * N_ * kv[2] >= 0
    for _ in range(3):
        pv = Fr(rng.randint(0, 999), 1000) * N_ / D_
        good &= kv[0] + kv[1] * pv + kv[2] * pv * pv > 0
        n_g += 1
    ok &= good
    count[reg] = count.get(reg, 0) + 1
ck.check("EXACT", "%d random rational points of {d = 1, both inequalities of (5) fail strictly, (c,j) != (0,0)} (seed 2614; "
         "with zeros among the coefficients and the layers b = 1, bq = 1, ab = 1, j = c, j = c/(1+c), j = cy/(c+y), "
         "q h_b = r^2): in the region of the document's cover, with its threshold P* = N/D written in the original "
         "coefficients, D > 0, N > 0, k0 > 0, both Bernstein numerators >= 0, and G(p) > 0 at 3 random p in [0, P*) each "
         "(%d values of G; %d points skipped for an irrational square root)" % (sum(count.values()), n_g, n_skip), ok,
         "points per region: %s" % dict(sorted(count.items())))
ck.info("PROOF-READ: in every region the two numerators are continuous functions of the original coefficients and of t "
        "(ab > 1); they are >= 0 on the image of the open chart, hence on its closure; the closures cover the domain; a "
        "zero of G has p >= P*, where one of (4), (6), (7), (9), (10) or the double-SOC endpoint applies")

# ----------------------------------------------------------------------------------------
if NUMERIC:
    import numpy as np
    ck.section("NUMERIC falsification attempts (floating point; evidence, not proof)")
    nrng = np.random.default_rng(26001)

    def cost_vals(co):
        w = [co.get(n, 0.0) for n in GNAMES]
        return np.array(K.combo_values(list(zip(GORDER, w))), float)

    def lu(lo=-1.5, hi=1.5):
        return float(10 ** nrng.uniform(lo, hi))

    def maybe(v, prob=0.6):
        return v if nrng.random() < prob else 0.0

    tests = {"(4) SOC pj = 1": [], "(4) old-H c(p-1) = 1": [], "(5) ab(e-1) = a+b+e+1": [], "(5) bq(c-1) = b+q+c+1": [],
             "(6)": [], "(7)": [], "double-SOC bq = 1, pc = 1": [], "(9)": [], "(10)": []}
    NS = 40
    for _ in range(NS):
        co = dict(d=1.0, p=lu(), a=maybe(lu()), b=maybe(lu()), e=maybe(lu()), c=maybe(lu()), q=maybe(lu()))
        co["j"] = 1.0 / co["p"]
        tests["(4) SOC pj = 1"].append(co)
        co = dict(d=1.0, c=lu(), a=maybe(lu()), b=maybe(lu()), e=maybe(lu()), j=maybe(lu()), q=maybe(lu()))
        co["p"] = 1.0 + 1.0 / co["c"]
        tests["(4) old-H c(p-1) = 1"].append(co)
        a_ = lu()
        b_ = (1 + lu()) / a_
        co = dict(d=1.0, a=a_, b=b_, e=(a_ + 1) * (b_ + 1) / (a_ * b_ - 1), p=maybe(lu()), j=maybe(lu()), c=maybe(lu()), q=maybe(lu()))
        tests["(5) ab(e-1) = a+b+e+1"].append(co)
        q_ = lu()
        b_ = (1 + lu()) / q_
        E2_ = (q_ + 1) * (b_ + 1) / (q_ * b_ - 1)
        co = dict(d=1.0, q=q_, b=b_, c=E2_, p=maybe(lu()), j=maybe(lu()), a=maybe(lu()), e=maybe(lu()))
        tests["(5) bq(c-1) = b+q+c+1"].append(co)
        c_ = E2_ * nrng.uniform(0.02, 1.0)
        co = dict(d=1.0, q=q_, b=b_, c=c_, p=max(0.0, 1 / c_ - 1 / E2_), j=maybe(lu()), a=maybe(lu()), e=maybe(lu()))
        tests["(6)"].append(co)
        b_ = lu()
        w_ = float(nrng.uniform(0, 1))
        if abs(b_ - w_) < 0.05:
            b_ += 0.2
        q_ = w_ * w_ / b_
        if nrng.random() < 0.5:
            b_, q_ = q_, b_
        c_ = lu()
        co = dict(d=1.0, q=q_, b=b_, c=c_, p=1 / c_ + (1 - w_) ** 2 / (np.sqrt(b_) - np.sqrt(q_)) ** 2,
                  j=maybe(lu()), a=maybe(lu()), e=maybe(lu()))
        tests["(7)"].append(co)
        b_ = lu()
        c_ = lu()
        co = dict(d=1.0, q=1 / b_, b=b_, c=c_, p=1 / c_, j=maybe(lu()), a=maybe(lu()), e=maybe(lu()))
        tests["double-SOC bq = 1, pc = 1"].append(co)
        # Lemma 1
        t_ = float(nrng.uniform(0, 1))
        a_ = lu()
        b_ = t_ * t_ / a_ + lu()
        hb_ = b_ - t_ * t_ / a_
        r_ = 1 - t_
        e_ = t_ * t_ + t_ * t_ * (a_ + b_ + 2 * t_) / (a_ * b_ - t_ * t_)
        c_ = lu()
        co = dict(d=1.0, a=a_, b=b_, c=c_, e=e_, q=r_ * r_ / hb_ * (1 + maybe(lu())), p=r_ * r_ / c_, j=maybe(lu()))
        tests["(9)"].append(co)
        w_ = float(nrng.uniform(0, 1))
        if abs(hb_ - w_ * r_) < 0.05 * hb_:
            w_ *= 0.5
        co = dict(d=1.0, a=a_, b=b_, c=c_, e=e_, q=w_ * w_ * r_ * r_ / hb_,
                  p=r_ * r_ / c_ + hb_ * r_ * r_ * (1 - w_) ** 2 / (hb_ - w_ * r_) ** 2, j=maybe(lu()))
        tests["(10)"].append(co)
    for nm, lst in tests.items():
        worst = 0.0
        for co in lst:
            mval = K.numeric_min_phi(cost_vals(co), nrng, restarts=8)
            worst = min(worst, mval / float(np.sum(np.abs(list(co.values())))))
        ck.check("NUMERIC", "capacity %s on its boundary: %d random coefficient vectors (seed 26001, coefficients log-uniform "
                 "in [10^-1.5, 10^1.5], unused coefficients 0 with probability 0.4), 8 BFGS restarts each over SO(4): "
                 "smallest Phi / (sum of coefficients) >= -1e-8" % (nm, NS), worst >= -1e-8, "smallest value %.3e" % worst)
    neg = 0
    for co in tests["(9)"][:20] + tests["(10)"][:20]:
        co2 = dict(co)
        co2["p"] = 0.0
        co2["j"] = 0.0
        if K.numeric_min_phi(cost_vals(co2), nrng, restarts=8) < -1e-6:
            neg += 1
    ck.info("control of the minimiser: with p and j replaced by 0, Phi takes a negative value in %d of the first 20 + 20 "
            "cases of (9), (10) (the capacities are sufficient conditions; sharpness is neither claimed nor tested)" % neg)

    # stratum points: real points of X with contact set containing T_26
    exps = np.array([ex[:8] for ex, _ in tG], float)
    cofs = np.array([co for _, co in tG], float)

    def G_float(wv):
        mons = np.prod(np.power(np.maximum(wv, 0.0)[None, :], exps), axis=1)
        return float(cofs @ mons), float(np.abs(cofs) @ mons)

    Bm = np.array([K.F_VALUES[jx] for jx in RORDER], float).T
    pinv = np.linalg.pinv(Bm)
    drow = pinv[5]
    pts = []
    for sdx, eps in ((26101, 0.02), (26102, 0.005), (26103, 0.1)):
        pts += K.numeric_stratum_search(T, sdx, 250, eps=eps)
    n_free = len(pts)
    n_forced = 0
    for kx, target in enumerate((0.003, 0.01, 0.03, 0.06, 0.1, 0.2)):
        got = K.numeric_stratum_search(T, 26200 + kx, 80, extra=lambda g, target=target: [float(drow @ g) - target])
        n_forced += len(got)
        pts += got
    n_hard = n_minor_pts = n_bad = 0
    r_cub = r_G = r_det = 0.0
    hard = []
    for pnt in pts:
        gaps = pnt["gaps"]
        wv = pinv @ gaps                                # coefficients (q, p, j, c, e, d, b, a): unique
        if np.max(np.abs(Bm @ wv - gaps)) > 1e-8 or np.min(wv) < -1e-8:
            n_bad += 1
            continue
        wv = np.maximum(wv, 0.0)
        if wv[5] < 1e-7:
            n_minor_pts += 1
            continue
        n_hard += 1
        hval = pnt["H"]
        qv, pv, jv, cvv, ev, dv, bv, av = wv
        Cs = (dv + qv) * (av + dv) * (bv + dv) + cvv * dv * (av + dv + ev) + dv * ev * (dv + qv)
        Qv = cvv * (pv - dv - ev) - dv * (dv + ev) - pv * jv
        r_cub = max(r_cub, abs(-Cs + hval * Qv) / (abs(Cs) + abs(hval * Qv)))
        gval, gscale = G_float(wv)
        r_G = max(r_G, abs(gval) / gscale)
        u = np.array([float(x_) for x_ in K.combo_coord(list(zip(RORDER, [float(x_) for x_ in wv])))])
        u[13] += hval
        Hm = np.array([[sum(cf * u[k] for cf, k in cell) for cell in row] for row in K.hess_linear_table()], float)
        sv = np.linalg.svd(Hm, compute_uv=False)
        r_det = max(r_det, sv[-1] / sv[0])
        hard.append(wv)
    ck.check("NUMERIC", "every converged stratum point is, within 1e-8, a nonnegative combination of the eight generators "
             "(Lemma 6.5 for T_26)", n_bad == 0, "%d exceptions" % n_bad)
    ck.info("stratum search: 750 free starts (seeds 26101-26103): %d converged points of X with gaps >= 0 vanishing on T_26; "
            "480 starts with the extra equation d = target (targets 0.003 ... 0.2; seeds 26200-26205): %d converged.  "
            "%d points with d < 1e-7 (minor decompositions), %d points with d > 1e-7 (d up to %.3g; normalisation: sum of "
            "gaps = 1)" % (n_free, n_forced, n_minor_pts, n_hard, max([float(w[5]) for w in hard] or [0.0])))
    if n_hard:
        ck.check("NUMERIC", "at the %d stratum points with d > 0: |-C_* + hQ| / (|C_*| + |hQ|) <= 1e-6, "
                 "|G| / (sum of |terms of G|) <= 1e-6 and smallest/largest singular value of Hess C(Re x) <= 1e-6 "
                 "(equations (2), (3) on actual points of X)" % n_hard,
                 r_cub <= 1e-6 and r_G <= 1e-6 and r_det <= 1e-6, "%.2e, %.2e, %.2e" % (r_cub, r_G, r_det))
        cap = {"SOC": 0, "old-H": 0, "four-1": 0, "four-2": 0, "(6)": 0, "(7)": 0, "double-SOC": 0, "(9)/(10)": 0, "none": 0}
        tol = 1e-7
        from scipy.optimize import brentq
        for wv in hard:
            qv, pv, jv, cvv, ev, dv, bv, av = wv / wv[5]
            hit = []
            if pv * jv >= 1 - tol:
                hit.append("SOC")
            if cvv * (pv - 1) >= 1 - tol:
                hit.append("old-H")
            if av * bv * (ev - 1) >= av + bv + ev + 1 - tol:
                hit.append("four-1")
            if bv * qv * (cvv - 1) >= bv + qv + cvv + 1 - tol:
                hit.append("four-2")
            if bv * qv > 1 + tol and cvv > 0 and pv >= 1 / cvv - (bv * qv - 1) / ((bv + 1) * (qv + 1)) - tol:
                hit.append("(6)")
            wq = np.sqrt(bv * qv)
            if bv * qv < 1 - tol and abs(bv - qv) > tol and cvv > 0 and \
                    pv >= 1 / cvv + (1 - wq) ** 2 / (np.sqrt(bv) - np.sqrt(qv)) ** 2 - tol:
                hit.append("(7)")
            if abs(bv * qv - 1) <= tol and pv * cvv >= 1 - tol:
                hit.append("double-SOC")
            # Lemma 1 with the largest admissible t: e >= e(t), ab > t^2
            if av > 0 and bv > 0 and cvv > 0:
                tmax = min(1.0, np.sqrt(av * bv)) * (1 - 1e-12)
                ef = lambda t_: t_ * t_ + t_ * t_ * (av + bv + 2 * t_) / (av * bv - t_ * t_) - ev
                tt = tmax if ef(tmax) <= 0 else brentq(ef, 0.0, tmax)
                hb_ = bv - tt * tt / av
                r_ = 1 - tt
                if qv * hb_ >= r_ * r_ - tol:
                    if pv >= r_ * r_ / cvv - tol:
                        hit.append("(9)/(10)")
                else:
                    w_ = np.sqrt(qv * hb_) / r_
                    if hb_ > w_ * r_ and pv >= r_ * r_ / cvv + hb_ * r_ * r_ * (1 - w_) ** 2 / (hb_ - w_ * r_) ** 2 - tol:
                        hit.append("(9)/(10)")
            if not hit:
                cap["none"] += 1
            for hh in hit:
                cap[hh] += 1
        ck.check("NUMERIC", "every stratum point with d > 0 satisfies at least one of the capacities (4)-(7), (9), (10) "
                 "(tolerance 1e-7 after scaling to d = 1)", cap["none"] == 0, "capacities met: %s" % cap)
    else:
        ck.info("no stratum point with d > 0 was found: equations (2), (3) are NOT tested on actual points by this search")

ck.finish()
