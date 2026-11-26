"""Independent replay of case T_25 = {3,4,7,8,10,12} (group V) of D. Zheng's proof of the
Marcus-de Oliveira conjecture for n = 4.

Proof document (read as PDF, never executed): root_moc4_T24_uniform_20261022.pdf, 4 pages
(the document's "T24" is the paper's T_25).

Run:   MOC4_SUPP=/path/to/supplement  python replay_T25.py [--no-numeric]

Exit code 0 iff every check passes.  EXACT checks use integers, rationals and exact polynomial
arithmetic (sympy, python-flint); lines marked NUMERIC are floating-point falsification attempts.
MOC4_SUPP is needed for: the sha256 comparison of the files that were read, the comparison of
the polynomial G computed here with the stored one, and the sixteen cubic features and the
8 x 8 change of basis of the sum-of-squares certificate behind Lemma 1 of the document (they
exist only in JSON data of the supplement).  Without MOC4_SUPP those parts are skipped and the
script FAILS, because Lemma 1 is then not replayed.
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
NPTS = 300 if "--quick" in sys.argv else 6000        # --quick: development only
ck = K.Checker("Replay of case T_25 = {3,4,7,8,10,12}  (document: root_moc4_T24_uniform_20261022.pdf)")
K.print_versions(ck)

T = K.PAPER_T[25]
GENS = [0, 1, 2, 6, 53, 54, 55, 59, 67, 70]
CANON = (0, 2, 59, 6, 53, 54, 55, 67)                 # coefficients (q, p, j, c, e, d, b, a)
CNAMES = ["q", "p", "j", "c", "e", "d", "b", "a"]

# ----------------------------------------------------------------------------------------
ck.section("Provenance")
have_supp = K.check_file_identity(ck, ["root_moc4_T24_uniform_20261022.pdf",
                                       "root_moc4_T24_full_radial_20261022.json",
                                       "root_moc4_T24_budget_PSD8_20261022.json",
                                       "root_moc4_T24_cubic_reduction_20261022.json",
                                       "root_moc4_T24_uniform_exact_20261022.json"])
ck.check("EXACT", "MOC4_SUPP is set and the five files are the ones that were read (needed for Lemma 1)", have_supp)

K.check_foundations(ck)
K.check_group_covariance(ck)
K.crosscheck_with_supplement_data(ck)
K.crosscheck_with_audit_library(ck, [25, 26])

# ----------------------------------------------------------------------------------------
K.check_contact_set(ck, 25, (3, 4, 7, 8, 10, 12), "T24", 26)
K.check_generators(ck, T, GENS, [54, 70])
K.check_minor_leibniz(ck, [j for j in GENS if K.F_IS_MINOR[j]])

# ----------------------------------------------------------------------------------------
ck.section("Reduction to the canonical representation (1)  [document p. 1]")
cv = K.combo_values
ck.check("EXACT", "f54 + f70 = f0 + f6 + f55 + f67 and f54 + f1 = f0 + f6 + f53 (all 24 values)",
         cv([(54, 1), (70, 1)]) == cv([(0, 1), (6, 1), (55, 1), (67, 1)])
         and cv([(54, 1), (1, 1)]) == cv([(0, 1), (6, 1), (53, 1)]))
GSYM = ((1, 3, 0, 2), (1, 3, 0, 2), 1)
gm_sym = K.g_index_map(GSYM)
IMG = {j: K.g_facet(gm_sym, j) for j in GENS}
ck.check("EXACT", "the element g = (rho, rho', transpose) = ((1,3,0,2), (1,3,0,2), 1) of G stabilises T_25, exchanges "
         "f54 and f70, fixes f2 and f59 and permutes the ten generators",
         K.g_set(gm_sym, T) == tuple(sorted(T)) and IMG[54] == 70 and IMG[70] == 54 and IMG[2] == 2
         and IMG[59] == 59 and sorted(IMG.values()) == GENS, "generator map %s" % IMG)
stab = [(g, gm) for g, gm in zip(K.G_ELEMENTS, K.G_MAPS) if K.g_set(gm, T) == tuple(sorted(T))]
GSWAP = ((3, 2, 1, 0), (3, 2, 1, 0), 0)
gm_swap = K.g_index_map(GSWAP)
IMG2 = {j: K.g_facet(gm_swap, j) for j in GENS}
ck.check("EXACT", "the stabiliser of T_25 in G has 4 elements; g' = ((3,2,1,0), (3,2,1,0), 0) is in it, fixes "
         "f1, f2, f53, f54, f59, f70 and exchanges f0 <-> f6, f55 <-> f67  (the symmetry (a,b,q,c) -> (b,a,c,q) of (1))",
         len(stab) == 4 and K.g_set(gm_swap, T) == tuple(sorted(T))
         and IMG2 == {0: 6, 6: 0, 55: 67, 67: 55, 1: 1, 2: 2, 53: 53, 54: 54, 59: 59, 70: 70})
ck.check("EXACT", "the eight value vectors of f_j, j in (0,2,59,6,53,54,55,67), have rank 8 "
         "(the canonical coefficients are unique)", K.rank([K.F_VALUES[j] for j in CANON]) == 8)


def reduce_to_canonical(lam):
    """The document's elimination on a dict {generator: coefficient >= 0}.  Returns
    (dict, used_symmetry): afterwards either f54, f70 are absent (minor decomposition) or the
    support lies in CANON with d > 0."""
    lam = {j: Fr(v) for j, v in lam.items() if v != 0}
    m = min(lam.get(54, 0), lam.get(70, 0))
    if m:
        lam[54] -= m
        lam[70] -= m
        for j in (0, 6, 55, 67):
            lam[j] = lam.get(j, 0) + m
    lam = {j: v for j, v in lam.items() if v != 0}
    used = False
    if lam.get(70, 0) > 0:                   # only f70 survives: pass to the mirror point g.x
        lam = {IMG[j]: v for j, v in lam.items()}
        used = True
    m = min(lam.get(54, 0), lam.get(1, 0))
    if m:
        lam[54] -= m
        lam[1] -= m
        for j in (0, 6, 53):
            lam[j] = lam.get(j, 0) + m
    return {j: v for j, v in lam.items() if v != 0}, used


rng = random.Random(2511)
ok = True
n_minor = n_canon = n_sym = 0
for _ in range(4000):
    mode = rng.random()
    lam = {j: Fr(rng.randint(0, 40), rng.randint(1, 6)) for j in GENS if rng.random() < mode}
    before = cv(lam)
    out, used = reduce_to_canonical(lam)
    target = list(K.g_values(gm_sym, before)) if used else before
    ok &= cv(out) == target and all(v > 0 for v in out.values())
    n_sym += used
    if out.get(54, 0) == 0:
        ok &= 70 not in out
        n_minor += 1
    else:
        ok &= set(out) <= set(CANON)
        n_canon += 1
ck.check("EXACT", "randomised exact run of the elimination on 4000 nonnegative rational coefficient vectors (seed 2511): "
         "the vector is unchanged (or replaced by its g-image when only f70 survives), coefficients stay >= 0, and the "
         "result is a minor decomposition (%d times) or the canonical form (1) with d > 0 (%d times; mirror used %d times)"
         % (n_minor, n_canon, n_sym), ok and n_minor > 0 and n_canon > 0 and n_sym > 0)

# ----------------------------------------------------------------------------------------
K.check_radial_necessity(ck)

# ----------------------------------------------------------------------------------------
ck.section("The cubic on g + h kappa, equation (2)  [p. 1]")
t0 = time.time()
rad = K.Radial(CANON, CNAMES)
R = rad.R
q, p, j, c, e, d, b, a, h = R.gens
I_ = c + d + q
U_ = (a + d) * (b + d)
V_ = U_ + d * e
S_ = a + b + d + e
alpha = I_ * U_ + c * q * (a + b + d)
Cst = alpha + e * (c + d) * (d + q)
Kc = a * b * (e - d) - d * d * (a + b + e + d)
Q = e * j - d * (d + e) - p * j
ck.check("EXACT", "C(g + h kappa) = -C_* + h Q with the document's C_* and Q = ej - d(d+e) - pj "
         "(g = q f0 + p f2 + j f59 + c f6 + e f53 + d f54 + b f55 + a f67; C is the paper's cubic (A.2))",
         rad.F == Q * h - Cst)
ck.check("EXACT", "C_* - d^3 and alpha - d^3 have only nonnegative coefficients: C_* >= alpha >= d^3 > 0 on the closed "
         "orthant with d > 0", K.coeff_stats(Cst - d ** 3)[1] >= 0 and K.coeff_stats(alpha - d ** 3)[1] >= 0
         and K.coeff_stats(Cst - alpha)[1] >= 0)
ck.info("so C(Re x) = 0 forces h Q = C_* > 0: Q != 0 and h = C_*/Q")

# ----------------------------------------------------------------------------------------
ck.section("The full 14 x 14 Hessian determinant, equation (3)  [pp. 1-2]")
Ddet = rad.det()
ck.info("det Hess C(g + h kappa): %d terms, %d minors memoised, degrees in (q,p,j,c,e,d,b,a,h) = %s, %.1f s"
        % (len(Ddet), rad.n_minors, R.degrees(Ddet), time.time() - t0))
rng = random.Random(2512)
ok = True
for _ in range(3):
    pt = [Fr(rng.randint(-7, 9), rng.randint(1, 4)) for _ in range(9)]
    u = [R.evaluate(x, pt) for x in rad.u]
    ok &= R.evaluate(Ddet, pt) == K.det_fraction(K.hess_at(u))
ck.check("EXACT", "the determinant polynomial (Laplace expansion, own code) equals the determinant of the rational matrix "
         "Hess C(u) computed by Gaussian elimination at 3 random rational points of (q,...,a,h) (seed 2512)", ok)
if "--bareiss" in sys.argv:                           # optional second algorithm (about 3 minutes)
    tb = time.time()
    ck.check("EXACT", "optional: the fraction-free Bareiss determinant equals the Laplace determinant term for term",
             K.bareiss_det(rad.H, R.zero, R.one) == Ddet, "%.0f s" % (time.time() - tb))
okf, G, Ssq = rad.radial_factor(Cst, Q, "j")
ck.check("EXACT", "Q^4 det Hess C(g + (C_*/Q) kappa) is divisible by -4 Q^2 and the quotient is the square of "
         "G := gcd(quotient, d quotient / dj): identity (3), Q^4 det = -4 Q^2 G^2 (degree of det in h: %d)" % rad.h_degree, okf)
if okf:
    d9 = [co for ex, co in R.terms(G) if ex == (0, 0, 0, 0, 0, 9, 0, 0, 0)]
    if d9 and d9[0] < 0:
        G = G * (-1)
    tG = R.terms(G)
    ck.check("EXACT", "G is homogeneous of degree 9 in the eight coefficients, free of h, of degree 2 in j, with 2095 "
             "integer monomials; the sign is fixed by a positive coefficient of d^9",
             all(sum(ex) == 9 and ex[8] == 0 and ex[2] <= 2 for ex, _ in tG) and len(tG) == 2095
             and dict(tG).get((0, 0, 0, 0, 0, 9, 0, 0, 0), 0) > 0,
             "%d terms, coefficient of d^9 = %s" % (len(tG), dict(tG).get((0, 0, 0, 0, 0, 9, 0, 0, 0))))
    sd = K.supp_dir()
    if sd is not None:
        import json
        with open(os.path.join(sd, "root_moc4_T24_full_radial_20261022.json")) as fh:
            radj = json.load(fh)
        gt = [f for f in radj["radial_factors"] if f["degree"] == 9][0]["terms"]
        Gj = R.from_terms([(tuple(ex), int(co)) for ex, co in gt])
        ck.check("EXACT", "supplement data: the stored degree-9 polynomial (variables q,p,j,c,e,d,b,a,z) equals the G "
                 "computed here, with the same sign", Gj == G and radj["variables"][:8] == CNAMES)
    ck.info("consequence: a point of X with canonical coefficients (1), d > 0, has G = 0 (necessary condition only)")

# ----------------------------------------------------------------------------------------
ck.section("Cited capacities (4): SOC, old-H, four-coefficient  [p. 2]")
g_soc = K.orbit_check(ck, "SOC, pj >= d^2: (p: f2, j: f59; d: f54)", (1, 9, 22), (2, 59, 54), sym_positions=(0, 1))
g_oh = K.orbit_check(ck, "old-H, e(j-d) >= d^2: (t = e: f53, c = j: f59; d: f54)", (12, 2, 70), (53, 59, 54))
g_four = K.orbit_check(ck, "four-coefficient, K >= 0: (a: f67, b: f55; c = e: f53; d: f54)", (0, 6, 10, 54),
                       (67, 55, 53, 54), sym_positions=(0, 1))
ck.info("CITED statements (PROOF.md 7.6): a f1 + b f9 + d f22 has Phi >= 0 on SO(4) if ab >= d^2; t f12 + c f2 + d f70 if "
        "t(c-d) >= d^2; a f0 + b f6 + c f10 + d f54 if ab(c-d) >= d^2(a+b+c+d); each on its whole G-orbit (6.6, checked above)")

# ----------------------------------------------------------------------------------------
ck.section("Lemma 1: the two coherent Gram capacities (5), (6) and the double-SOC endpoint  [p. 2]")
gf = None
if K.supp_dir() is not None:
    import sympy as sp
    gf = K.GramFace(K.supp_dir(), "root_moc4_T24")
    ck.check("EXACT", "certificate data (JSON): source order (0,2,59,6,53,54,55,67); the four signed coordinate maps "
             "permute p and q separately (so they preserve |p|^2, |q|^2)", gf.order == list(CANON) and gf.signed_perms_ok)
    ck.check("EXACT", "quaternion model M(p,q) = matrix of x -> p x q: M^t M = |p|^2 |q|^2 I and det M = (|p|^2 |q|^2)^2 as "
             "polynomial identities (so M(p,q) is in SO(4) for unit p, q)", gf.model_ok)
    wr = gf.weight_rows
    key = {k: i for i, k in enumerate(gf.keys)}

    def unit(i, k):
        v = [Fr(0)] * 36
        v[key[(i, k)]] = Fr(1)
        return v
    rel_space = gf.relations
    w_expected = {1: unit(1, 1), 2: unit(2, 2), 3: unit(3, 3), 4: unit(4, 4), 6: unit(6, 6), 7: unit(7, 7)}
    w0 = [x - y for x, y in zip(wr[0], unit(0, 0))]
    w5 = [x - y + z for x, y, z in zip(wr[5], unit(5, 5), unit(0, 0))]
    ck.check("EXACT", "own linear algebra on (*) (see GramFace; %d coefficient equations, rank %d): the solutions "
             "(Gram matrix, weights) are exactly: 10 independent linear relations on the Gram matrix, and weights "
             "w = (G00, G11, G22, G33, G44, G55 - G00, G66, G77) modulo the relations, i.e. Gram diagonal "
             "(q, p, j, c, e, q + d, b, a)" % (gf.n_equations, gf.rank),
             len(rel_space) == 10 and all(wr[m] == w_expected[m] for m in w_expected)
             and K.rank(rel_space + [w0]) == 10 and K.rank(rel_space + [w5]) == 10)
    ck.check("EXACT", "the ten relations stored in the supplement span the same space as the ten derived here",
             K.rank(gf.stored_R) == 10 and K.rank(gf.stored_R + rel_space) == 10)
    ck.info("relations derived here (Gik = Gram entries, vectors in the order (f0,f2,f59,f6,f53,f54,f55,f67)):")
    for r in rel_space:
        ck.info("   " + " ".join("%+d*G%d%d" % (int(cf), gf.keys[i][0], gf.keys[i][1]) for i, cf in enumerate(r) if cf) + " = 0")
    ck.info("consequence (PROOF-READ + CITED surjectivity of (p,q) -> M(p,q) onto SO(4)): if a symmetric PSD matrix "
            "satisfies the ten relations, then Phi of the cost with these weights is >= 0 on SO(4)")

    sa, sb, se, sj, sx, ss, sal, sbe, sp_ = sp.symbols("a b e j x s alpha beta p", positive=True)
    L = (sa + 1) * (sb + 1) / (sa * sb - 1)
    # (5), general x: J.E = x-1, E.D = E.B = x, E.A = -x, B.A = 1, |D| = 1
    emin_x = (1 - sx) ** 2 / sj + L * sx ** 2
    K.sym_gram_check(ck, gf, "(5), d = 1, general x, e = (1-x)^2/j + L x^2 + s",
                     {(2, 4): sx - 1, (4, 5): sx, (4, 6): sx, (4, 7): -sx, (6, 7): 1},
                     [0, 0, sj, 0, emin_x + ss, 1, sb, sa], [0, 0, sj, 0, emin_x + ss, 1, sb, sa],
                     [2, 5, 6, 7, 4], [sj, 1, sb, (sa * sb - 1) / sb, ss],
                     {sa: Fr(3), sb: Fr(2), sj: Fr(1, 2), sx: Fr(1, 3), ss: Fr(1, 5)},
                     note="  [L = (a+1)(b+1)/(ab-1)]")
    xopt = 1 / (1 + sj * L)
    ck.check("EXACT", "(5): (1-x)^2/j + L x^2 = L/(1+jL) + ((1+jL)/j)(x - 1/(1+jL))^2, so the least |E|^2 is "
             "L/(1+jL) = 1/(j + 1/L), at x = 1/(1+jL); and e >= 1/(j+1/L) <=> j >= 1/e - 1/L (e, L > 0); in terms of d: "
             "1/L = d(ab-d^2)/((a+d)(b+d)) after scaling by d",
             sp.cancel(emin_x - L / (1 + sj * L) - (1 + sj * L) / sj * (sx - xopt) ** 2) == 0
             and sp.cancel(1 / (sj + 1 / L) - L / (1 + sj * L)) == 0)
    K.sym_gram_check(ck, gf, "(5) at j = 0 (J = 0, x = 1): the four-coefficient inequality e >= L re-proved on this face",
                     {(4, 5): 1, (4, 6): 1, (4, 7): -1, (6, 7): 1},
                     [0, 0, 0, 0, L + ss, 1, sb, sa], [0, 0, 0, 0, L + ss, 1, sb, sa],
                     [5, 6, 7, 4], [1, sb, (sa * sb - 1) / sb, ss],
                     {sa: Fr(5, 2), sb: Fr(2), ss: Fr(1, 7)})
    # (6): a = alpha^2, b = beta^2, w = alpha beta, B = beta v, A = alpha v, J = Jv v - E/e + ...
    Jv = (1 - sal * sbe) / (sbe - sal)
    K.sym_gram_check(ck, gf, "(6), d = 1, a = alpha^2, b = beta^2, w = alpha beta, j = 1/e + (1-w)^2/(beta-alpha)^2 + s",
                     {(2, 4): -1, (2, 6): Jv * sbe, (2, 7): Jv * sal, (6, 7): sal * sbe},
                     [0, 0, 1 / se + Jv ** 2 + ss, 0, se, 1, sbe ** 2, sal ** 2],
                     [0, 0, 1 / se + Jv ** 2 + ss, 0, se, 1, sbe ** 2, sal ** 2],
                     [6, 4, 2, 5], [sbe ** 2, se, ss, 1],
                     {sal: Fr(1, 3), sbe: Fr(3, 2), se: Fr(2), ss: Fr(1, 4)},
                     note="  [vector A is alpha/beta times vector B: no pivot]")
    sw = sp.Symbol("w", positive=True)
    ck.check("EXACT", "(6): (1-w)^2/(sqrt b - sqrt a)^2 = b(1-w)^2/(b-w)^2 for w = sqrt(ab) (with a = w^2/b), and the "
             "expression is symmetric in a, b: the capacity (6) holds for a != b in either order",
             sp.cancel((1 - sal * sbe) ** 2 / (sbe - sal) ** 2
                       - sbe ** 2 * (1 - sal * sbe) ** 2 / (sbe ** 2 - sal * sbe) ** 2) == 0
             and sp.cancel(sal ** 2 * (1 - sal * sbe) ** 2 / (sal ** 2 - sal * sbe) ** 2
                           - sbe ** 2 * (1 - sal * sbe) ** 2 / (sbe ** 2 - sal * sbe) ** 2) == 0)
    K.sym_gram_check(ck, gf, "double-SOC endpoint: blocks (J,E) with J.E = -1 and (B,A) with B.A = 1, D alone",
                     {(2, 4): -1, (6, 7): 1},
                     [0, 0, sj, 0, se, 1, sb, sa], [0, 0, sj, 0, se, 1, sb, sa],
                     [2, 4, 5, 6, 7], [sj, se - 1 / sj, 1, sb, sa - 1 / sb],
                     {sa: Fr(2), sb: Fr(1, 2), sj: Fr(3), se: Fr(1, 3)},
                     note="  [PSD iff ej >= 1 and ab >= 1]")
    K.sym_gram_check(ck, gf, "SOC pj >= 1 re-proved on this face: block (P,J) with P.J = 1, D alone",
                     {(1, 2): 1}, [0, sp_, sj, 0, 0, 1, 0, 0], [0, sp_, sj, 0, 0, 1, 0, 0],
                     [1, 2, 5], [sp_, sj - 1 / sp_, 1], {sp_: Fr(2), sj: Fr(1, 2)})
    # old-H has no certificate of this form when only e, j, d are present: Gram matrices supported on (J, E, D)
    sup = [2, 4, 5]
    cols = [(i, k) for i in sup for k in sup if i <= k and (i, k) != (5, 5)]
    Arows = [[r[key[ck_]] for ck_ in cols] for r in rel_space]
    brhs = [-r[key[(5, 5)]] for r in rel_space]                 # G55 = 1 (d = 1, q = 0)
    ck.check("EXACT", "the ten relations have no solution with a Gram matrix supported on the vectors (J, E, D) and "
             "G55 = 1: the old-H capacity e(j-d) >= d^2 alone (q = p = c = a = b = 0) is NOT certified by this Gram "
             "representation; it is a genuine citation",
             K.rank(Arows) < K.rank([row + [bb] for row, bb in zip(Arows, brhs)]))
    ck.info("PROOF-READ: all pivots are >= 0 under the stated conditions (j > 0, ab > 1, s >= 0; beta > alpha >= 0, e > 0); "
            "limits of PSD matrices are PSD; the weights are homogeneous, so d = 1 loses nothing; adding q f0, p f2, "
            "c f6 (minor vectors, Phi >= 0) back keeps Phi >= 0")
else:
    ck.check("EXACT", "Lemma 1 replayed", False, "MOC4_SUPP not available")

# ----------------------------------------------------------------------------------------
ck.section("The constant coefficient k0, identity (7), and the layer p = e = 0  [p. 3]")
tG = R.terms(G)
k0h = R.from_terms([(ex, co) for ex, co in tG if ex[2] == 0])
beta = Kc * (-1) + (p + c + q) * V_ + c * q * S_
Wc = p * (b + d)
Wq = p * (a + d)
Hh = V_ * (c * (c + d) * Wc + q * (q + d) * Wq) \
    + c * q * (c * q * p * S_ + c * (S_ * Wc + p * V_) + q * (S_ * Wq + p * V_) + p * (a + b + d * 2) * V_ + a * b * (d + e) * p)
ck.check("EXACT", "k0 = G|_{j=0} = C_* { alpha beta + c q d^2 p S + H } with the document's beta, H (d kept general)",
         k0h == Cst * (alpha * beta + c * q * d * d * p * S_ + Hh))
ck.check("EXACT", "H and beta + K have only nonnegative coefficients (so beta >= -K > 0 when K < 0, and k0 >= C_* alpha beta "
         ">= d^6 (-K) > 0 on {K < 0, d > 0}, all other coefficients >= 0)",
         K.coeff_stats(Hh)[1] >= 0 and K.coeff_stats(beta + Kc)[1] >= 0)
Gpe = R.from_terms([(ex, co) for ex, co in tG if ex[1] == 0 and ex[4] == 0])
st = K.coeff_stats(Gpe)
ck.check("EXACT", "G|_{p=e=0} has only nonnegative coefficients (%d terms) and contains d^9: G > 0 there for every j >= 0, "
         "so no point of X with d > 0 has p = e = 0" % st[0],
         st[1] >= 0 and dict(R.terms(Gpe)).get((0, 0, 0, 0, 0, 9, 0, 0, 0), 0) > 0)
swap = [c, p, j, q, e, d, a, b, h]            # (q,p,j,c,e,d,b,a,h) -> q<->c, b<->a
ck.check("EXACT", "G is invariant under (a,b,q,c) -> (b,a,c,q) (the action of g' above); also K, C_*, Q are",
         K.subst_clear(tG, swap, [R.one] * 9, R) == G
         and K.subst_clear(R.terms(Kc), swap, [R.one] * 9, R) == Kc
         and K.subst_clear(R.terms(Cst), swap, [R.one] * 9, R) == Cst
         and K.subst_clear(R.terms(Q), swap, [R.one] * 9, R) == Q)

# ----------------------------------------------------------------------------------------
ck.section("Bernstein criterion and the seven charts  [pp. 3-4]")
ck.check("EXACT", "for G = k0 + k1 x + k2 x^2 and x = t N/D: D^2 G = (1-t)^2 D^2 k0 + t(1-t) D (2D k0 + N k1) "
         "+ t^2 (D^2 k0 + N D k1 + N^2 k2); so D > 0, k0 > 0 and both numerators >= 0 give G > 0 on [0, N/D)",
         K.bernstein_identity_check())
RC = K.Ring(["a", "b", "e", "p", "q", "c"])                 # d = 1
ca, cb, ce, cp, cq, cc = RC.gens
kk = [RC.from_terms([((ex[7], ex[6], ex[4], ex[1], ex[0], ex[3]), co) for ex, co in tG if ex[2] == n]) for n in range(3)]
ck.info("d = 1: G = k2 j^2 + k1 j + k0 with k0, k1, k2 in Z[a,b,e,p,q,c] of %s terms" % [len(x) for x in kk])
RO = K.Ring(["A", "T", "V", "U", "Q", "C"])
A_, T_, V__, U__, Qo, Co = RO.gens
one = RO.one
oneC = RC.one
En = (A_ + one) * (A_ + one + T_) * V__          # e = L V/(1+V) with ab = 1 + T: L = (A+1)(A+1+T)/(A T)
Ed = A_ * T_ * (one + V__)
Fw = (one + T_) * (one + V__) - V__              # (b - w)(1 + V) for b = 1 + T, w = V/(1+V)
LcN = (ca + 1) * (cb + 1)                        # L = LcN / (ab - 1)

# chart: (nums, dens) for (a, b, e, p, q, c);  threshold J* = N/D as polynomials of the ORIGINAL coefficients
# (or of the chart parameters for IIIb-);  expected (middle, endpoint) term counts from p. 4.
CHARTS = [
    ("I-", [A_, T_, V__, V__ + U__, Qo, Co], [one, A_ * (one + T_), one, one, one, one], (oneC, cp), None, (3218, 4620)),
    ("I+", [A_, one + T_, En, En + U__ * Ed, Qo, Co], [one, A_, Ed, Ed, one, one], (oneC, cp), None, (8961, 16478)),
    ("II", [A_, one + T_, En, En * U__, Qo, Co], [one, A_, Ed, Ed * (one + U__), one, one],
     (LcN - ce * (ca * cb - 1), ce * LcN), None, (18810, 31918)),
    ("IIIa-", [T_ * V__, T_, A_, A_ * U__, Qo, Co], [(one + T_) * (one + V__), one + T_, one, (one + A_) * (one + U__), one, one],
     (ce + 1, ce), None, (7979, 7979)),
    ("IIIa+", [T_ * V__, T_, A_, A_ * (one + U__) + A_ * A_ * U__, Qo, Co],
     [(one + T_) * (one + V__), one + T_, one, (one + A_) * (one + U__), one, one], (oneC, cp), None, (7980, 11532)),
    ("IIIb+", [V__, one + T_, A_, A_ * (one + T_) * (one + U__) + A_ * A_ * U__, Qo, Co],
     [(one + T_) * (one + V__), one, one, (A_ + one + T_) * (one + U__), one, one], (oneC, cp), None, (8852, 14016)),
    ("IIIb-", [V__ * V__, one + T_, A_, A_ * (one + T_) * U__, Qo, Co],
     [(one + V__) ** 2 * (one + T_), one, one, (A_ + one + T_) * (one + U__), one, one],
     None, (Fw * Fw + A_ * (one + T_), A_ * Fw * Fw), (23897, 33654)),
]


def chart_point(name, A, Tt, V, U):
    """The chart maps of the table on p. 3, as plain rational arithmetic (independent of the
    polynomial encoding above).  Returns (a, b, e, p, J*) and the region test."""
    if name == "I-":
        a_, b_, e_ = A, Tt / (A * (1 + Tt)), V
        p_ = V + U
        return a_, b_, e_, p_, 1 / p_, (p_ >= e_ and a_ * b_ < 1)
    if name in ("I+", "II"):
        a_, b_ = A, (1 + Tt) / A
        L_ = (a_ + 1) * (b_ + 1) / (a_ * b_ - 1)
        e_ = L_ * V / (1 + V)
        if name == "I+":
            p_ = e_ + U
            return a_, b_, e_, p_, 1 / p_, (p_ >= e_ and a_ * b_ > 1 and e_ < L_)
        p_ = e_ * U / (1 + U)
        return a_, b_, e_, p_, 1 / e_ - 1 / L_, (e_ >= p_ and a_ * b_ > 1 and e_ < L_)
    if name in ("IIIa-", "IIIa+"):
        b_ = Tt / (1 + Tt)
        a_ = b_ * V / (1 + V)
        e_ = A
        if name == "IIIa-":
            p_ = A * U / ((1 + A) * (1 + U))
            return a_, b_, e_, p_, 1 + 1 / e_, (a_ <= b_ <= 1 and p_ <= e_ / (1 + e_))
        p_ = A / (1 + A) + A * A * U / ((1 + A) * (1 + U))
        return a_, b_, e_, p_, 1 / p_, (a_ <= b_ <= 1 and e_ / (1 + e_) <= p_ <= e_)
    b_ = 1 + Tt
    e_ = A
    if name == "IIIb+":
        a_ = V / (b_ * (1 + V))
        p_ = A * b_ / (A + b_) + A * A * U / ((A + b_) * (1 + U))
        return a_, b_, e_, p_, 1 / p_, (b_ >= 1 and a_ * b_ < 1 and e_ * b_ / (e_ + b_) <= p_ <= e_)
    w_ = V / (1 + V)
    a_ = w_ * w_ / b_
    p_ = A * b_ * U / ((A + b_) * (1 + U))
    return a_, b_, e_, p_, 1 / e_ + b_ * (1 - w_) ** 2 / (b_ - w_) ** 2, (b_ >= 1 and a_ * b_ < 1 and p_ <= e_ * b_ / (e_ + b_))


rng = random.Random(2513)
TABLES = {}
for name, nums, dens, ND_orig, ND_out, expect in CHARTS:
    t0 = time.time()
    if ND_orig is not None:
        mid, end = K.bernstein_polys(RC, kk, ND_orig[0], ND_orig[1])
        Fm = K.subst_clear(RC.terms(mid), nums, dens, RO)
        Fe = K.subst_clear(RC.terms(end), nums, dens, RO)
        dm, de = RC.degrees(mid), RC.degrees(end)
    else:
        degs = [max(RC.degrees(x)[i] for x in kk) for i in range(6)]
        sk = [K.subst_clear(RC.terms(x), nums, dens, RO, degs) for x in kk]
        Fm, Fe = K.bernstein_polys(RO, sk, ND_out[0], ND_out[1])
        dm = de = degs
    sm, se_ = K.coeff_stats(Fm), K.coeff_stats(Fe)
    TABLES[name] = (Fm, Fe)
    ck.check("EXACT", "chart %s: the middle and endpoint Bernstein numerators, after the substitution and clearing of "
             "positive denominators, have only nonnegative integer coefficients; term counts %d / %d as in the document"
             % (name, sm[0], se_[0]), sm[1] >= 0 and se_[1] >= 0 and (sm[0], se_[0]) == expect,
             "min coefficients %d / %d, %.1f s, sha256 %s.. / %s.." % (sm[1], se_[1], time.time() - t0,
                                                                     K.poly_sha(Fm)[:12], K.poly_sha(Fe)[:12]))
    # the encoding against the plain chart map: exact evaluation at random positive rational parameters
    ok = True
    for _ in range(4):
        prm = [Fr(rng.randint(1, 9), rng.randint(1, 5)) for _ in range(6)]
        a_, b_, e_, p_, Jst, inreg = chart_point(name, prm[0], prm[1], prm[2], prm[3])
        pt = [a_, b_, e_, p_, prm[4], prm[5]]
        kv = [RC.evaluate(x, pt) for x in kk]
        Nn, Dd = Jst.numerator, Jst.denominator
        mid_v = 2 * Dd * kv[0] + Nn * kv[1]
        end_v = Dd * Dd * kv[0] + Nn * Dd * kv[1] + Nn * Nn * kv[2]
        fm_v, fe_v = RO.evaluate(Fm, prm), RO.evaluate(Fe, prm)
        dv = [RO.evaluate(x, prm) for x in dens]
        ok &= inreg and all(x > 0 for x in dv) and Jst > 0 and kv[0] > 0
        # same sign and same zero set: the table polynomial is a positive multiple of the numerator
        ok &= fm_v > 0 and fe_v > 0 and mid_v > 0 and end_v > 0
        if ND_orig is not None:
            Nv, Dv = RC.evaluate(ND_orig[0], pt), RC.evaluate(ND_orig[1], pt)
            ok &= Dv > 0 and Nv / Dv == Jst
            scale_m = Fr(1)
            scale_e = Fr(1)
            for i in range(6):
                scale_m *= dv[i] ** dm[i]
                scale_e *= dv[i] ** de[i]
            ok &= fm_v == (2 * Dv * kv[0] + Nv * kv[1]) * scale_m
            ok &= fe_v == (Dv * Dv * kv[0] + Nv * Dv * kv[1] + Nv * Nv * kv[2]) * scale_e
        else:
            Nv, Dv = RO.evaluate(ND_out[0], prm), RO.evaluate(ND_out[1], prm)
            sc = Fr(1)
            for i in range(6):
                sc *= dv[i] ** dm[i]
            ok &= Dv > 0 and Nv / Dv == Jst
            ok &= fm_v == (2 * Dv * kv[0] + Nv * kv[1]) * sc and fe_v == (Dv * Dv * kv[0] + Nv * Dv * kv[1] + Nv * Nv * kv[2]) * sc
    ck.check("EXACT", "chart %s: at 4 random positive rational parameter points (seed 2513) the plain chart map of p. 3 "
             "lands in its region, N/D is the stated threshold J*, and the table polynomials equal the Bernstein "
             "numerators times a positive product of denominator powers" % name, ok)

if K.supp_dir() is not None:
    import json
    with open(os.path.join(K.supp_dir(), "root_moc4_T24_uniform_exact_20261022.json")) as fh:
        stored = json.load(fh)["tables"]
    mine = []
    for name, *_ in CHARTS:
        mine += [K.poly_sha(TABLES[name][0]), K.poly_sha(TABLES[name][1])]
    ck.check("EXACT", "supplement data (comparison only): the sha256 of each of my fourteen table polynomials (sorted "
             "'exponents:coefficient' lines, variables A,T,V,U,Q,C) equals the hash recorded in the supplement's result "
             "file, in the order I-, I+, II, IIIa-, IIIa+, IIIb+, IIIb-: the polynomials are identical to the author's",
             len(stored) == 14 and [t["sha256"] for t in stored] == mine
             and [t["terms"] for t in stored] == [len(TABLES[n][i]) for n, *_ in CHARTS for i in (0, 1)])

# ----------------------------------------------------------------------------------------
ck.section("Cover of the coefficient domain, boundaries included: exact randomised test in the ORIGINAL coefficients")


def L_of(a_, b_):
    return (a_ + 1) * (b_ + 1) / (a_ * b_ - 1)


def rnd(rng, zero_prob=0.0):
    if rng.random() < zero_prob:
        return Fr(0)
    return Fr(rng.randint(1, 60), rng.randint(1, 60)) * Fr(10) ** rng.randint(-2, 2)


def region_sample(rng):
    """A random rational point (a,b,e,p,q,c) of the domain {d = 1, K < 0, (p,e) != (0,0)} with many
    boundary cases, the region of the document's cover it is assigned to (after a <-> b if needed), and the
    threshold J* = N/D of that region written in the original coefficients."""
    mode = rng.randint(0, 9)
    q_, c_ = rnd(rng, 0.35), rnd(rng, 0.35)
    a_, b_ = rnd(rng, 0.15), rnd(rng, 0.15)
    if mode == 0 and a_ > 0:
        b_ = 1 / a_                                   # ab = 1 exactly
    if mode == 1:
        b_ = Fr(1)                                    # b = 1 exactly
    if mode == 2:
        b_ = a_                                       # a = b
    e_ = rnd(rng, 0.15)
    if a_ * b_ > 1:
        e_ = L_of(a_, b_) * Fr(rng.randint(0, 99), 100)       # K < 0  <=>  e < L
    pm = rng.randint(0, 5)
    if pm == 0:
        p_ = e_                                       # p = e
    elif pm == 1:
        p_ = Fr(0)
    elif pm == 2:
        p_ = e_ / (1 + e_)                            # the cut of IIIa
    elif pm == 3:
        m_ = max(a_, b_)
        p_ = e_ * m_ / (e_ + m_) if e_ + m_ > 0 else Fr(0)     # the cut of IIIb
    else:
        p_ = rnd(rng)
    if p_ == 0 and e_ == 0:
        e_ = rnd(rng)
        if a_ * b_ > 1:
            e_ = L_of(a_, b_) * Fr(rng.randint(1, 99), 100)
    return a_, b_, e_, p_, q_, c_


def threshold(a_, b_, e_, p_):
    """(region, N, D, extra) of the document's cover for a point with K < 0, (p,e) != (0,0); extra is the rational
    w = sqrt(ab) when region IIIb- needs it (None if ab is not a rational square: then IIIb- is not testable exactly
    and the point is skipped)."""
    if p_ >= e_:
        return "I", Fr(1), p_
    if a_ * b_ >= 1:
        return "II", (a_ + 1) * (b_ + 1) - e_ * (a_ * b_ - 1), e_ * (a_ + 1) * (b_ + 1)
    m_ = max(a_, b_)
    if m_ <= 1:
        if p_ <= e_ / (1 + e_):
            return "IIIa-", e_ + 1, e_
        return "IIIa+", Fr(1), p_
    if p_ >= e_ * m_ / (e_ + m_):
        return "IIIb+", Fr(1), p_
    return "IIIb-", None, None


def exact_sqrt(x):
    from math import isqrt
    n_, d_ = x.numerator, x.denominator
    rn, rd = isqrt(n_), isqrt(d_)
    return Fr(rn, rd) if rn * rn == n_ and rd * rd == d_ else None


rng = random.Random(2514)
count = {}
ok = True
worst = None
n_g = 0
for it in range(NPTS):
    a_, b_, e_, p_, q_, c_ = region_sample(rng)
    if it % 5 == 0:                                   # force IIIb- with a rational w = sqrt(ab)
        m_ = 1 + rnd(rng, 0.2)
        w_ = Fr(rng.randint(0, 99), 100)
        lo = w_ * w_ / m_
        a_, b_ = (lo, m_) if rng.random() < 0.5 else (m_, lo)
        e_ = rnd(rng)
        p_ = e_ * m_ / (e_ + m_) * Fr(rng.randint(0, 100), 100)
    Kv = a_ * b_ * (e_ - 1) - (a_ + b_ + e_ + 1)
    assert Kv < 0 and (p_ > 0 or e_ > 0)
    reg, N_, D_ = threshold(a_, b_, e_, p_)
    if reg == "IIIb-":
        w_ = exact_sqrt(a_ * b_)
        if w_ is None:
            continue
        m_ = max(a_, b_)
        N_ = (m_ - w_) ** 2 + e_ * m_ * (1 - w_) ** 2
        D_ = e_ * (m_ - w_) ** 2
    pt = [a_, b_, e_, p_, q_, c_]
    kv = [RC.evaluate(x, pt) for x in kk]
    mid_v = 2 * D_ * kv[0] + N_ * kv[1]
    end_v = D_ * D_ * kv[0] + N_ * D_ * kv[1] + N_ * N_ * kv[2]
    good = D_ > 0 and N_ > 0 and kv[0] > 0 and mid_v >= 0 and end_v >= 0
    for _ in range(3):                                # G > 0 strictly below the threshold
        tt = Fr(rng.randint(0, 999), 1000)
        jv = tt * N_ / D_
        gv = kv[0] + kv[1] * jv + kv[2] * jv * jv
        good &= gv > 0
        n_g += 1
    ok &= good
    count[reg] = count.get(reg, 0) + 1
ck.check("EXACT", "%d random rational points of {d = 1, K < 0, (p,e) != (0,0)} (seed 2514; with zeros among a, b, q, c, e, p, "
         "and the layers ab = 1, b = 1, a = b, p = e, p = e/(1+e), p = e m/(e+m), both orders of a, b): in the region of "
         "the document's cover, with its threshold J* = N/D written in the original coefficients, D > 0, N > 0, k0 > 0, "
         "both Bernstein numerators >= 0, and G(j) > 0 at 3 random j in [0, J*) each (%d values of G)" % (NPTS, n_g), ok,
         "points per region: %s" % dict(sorted(count.items())))
ck.info("PROOF-READ: each region's numerators are polynomials (IIIb-: continuous functions) of the original coefficients; "
        "they are >= 0 on the image of the open chart, hence on its closure; the closures cover the domain; at G = 0 the "
        "point has j >= J*, where (4), (5), (6) or the double-SOC endpoint applies")

# ----------------------------------------------------------------------------------------
if NUMERIC:
    import numpy as np
    ck.section("NUMERIC falsification attempts (floating point; evidence, not proof)")
    nrng = np.random.default_rng(25001)

    def cost_vals(co):
        """24 values of q f0 + p f2 + j f59 + c f6 + e f53 + d f54 + b f55 + a f67 for co = dict."""
        w = [co.get(n, 0.0) for n in CNAMES]
        return np.array(K.combo_values(list(zip(CANON, w))), float)

    def lu(lo=-1.5, hi=1.5):
        return float(10 ** nrng.uniform(lo, hi))

    def extras():
        return {n: (lu() if nrng.random() < 0.5 else 0.0) for n in ("q", "p", "c")}

    tests = {"(4) SOC pj = d^2": [], "(4) old-H e(j-d) = d^2": [], "(4) four K = 0": [], "(5)": [], "(6)": [],
             "double-SOC ab = 1, ej = 1": []}
    for _ in range(40):
        co = extras()
        co.update(d=1.0, p=lu(), a=lu() * (nrng.random() < 0.7), b=lu() * (nrng.random() < 0.7), e=lu() * (nrng.random() < 0.7))
        co["j"] = 1.0 / co["p"]
        tests["(4) SOC pj = d^2"].append(co)
        co = extras()
        co.update(d=1.0, e=lu(), a=lu() * (nrng.random() < 0.7), b=lu() * (nrng.random() < 0.7))
        co["j"] = 1.0 + 1.0 / co["e"]
        tests["(4) old-H e(j-d) = d^2"].append(co)
        co = extras()
        a_ = lu()
        b_ = (1 + lu()) / a_
        co.update(d=1.0, a=a_, b=b_, e=(a_ + 1) * (b_ + 1) / (a_ * b_ - 1), j=lu() * (nrng.random() < 0.5))
        tests["(4) four K = 0"].append(co)
        co = extras()
        a_ = lu()
        b_ = (1 + lu()) / a_
        L_ = (a_ + 1) * (b_ + 1) / (a_ * b_ - 1)
        e_ = L_ * nrng.uniform(0.02, 1.0)
        co.update(d=1.0, a=a_, b=b_, e=e_, j=max(0.0, 1 / e_ - 1 / L_))
        tests["(5)"].append(co)
        co = extras()
        b_ = lu()
        w_ = float(nrng.uniform(0, 1))
        if abs(b_ - w_) < 0.05:
            b_ += 0.2
        a_ = w_ * w_ / b_
        e_ = lu()
        jj = 1 / e_ + b_ * (1 - w_) ** 2 / (b_ - w_) ** 2
        co.update(d=1.0, e=e_, j=jj)
        if nrng.random() < 0.5:
            co.update(a=a_, b=b_)
        else:
            co.update(a=b_, b=a_)
        tests["(6)"].append(co)
        co = extras()
        a_ = lu()
        e_ = lu()
        co.update(d=1.0, a=a_, b=1 / a_, e=e_, j=1 / e_)
        tests["double-SOC ab = 1, ej = 1"].append(co)
    for nm, lst in tests.items():
        worst = 0.0
        for co in lst:
            v = cost_vals(co)
            mval = K.numeric_min_phi(v, nrng, restarts=8)
            worst = min(worst, mval / float(np.sum(np.abs(list(co.values())))))
        ck.check("NUMERIC", "capacity %s on its boundary: 40 random coefficient vectors (seed 25001, coefficients "
                 "log-uniform in [10^-1.5, 10^1.5], q, p, c and some others set to 0 at random), 8 BFGS restarts each over "
                 "SO(4): smallest Phi / (sum of coefficients) >= -1e-8" % nm, worst >= -1e-8, "smallest value %.3e" % worst)
    # violation control: just below the capacity (5) and (6) the minimum should be allowed to go negative
    neg = neg2 = 0
    for co in tests["(5)"][:20] + tests["(6)"][:20]:
        co2 = dict(co)
        co2["q"] = co2["p"] = co2["c"] = 0.0
        co2["j"] = 0.0
        if K.numeric_min_phi(cost_vals(co2), nrng, restarts=8) < -1e-6:
            neg += 1
        co2["j"] = co["j"] * 0.5
        if K.numeric_min_phi(cost_vals(co2), nrng, restarts=8) < -1e-6:
            neg2 += 1
    ck.info("control of the minimiser: with q = p = c = 0 and j replaced by 0, Phi takes a negative value in %d of the first "
            "20 + 20 cases of (5), (6); with j halved, in %d (the capacities are sufficient conditions; no sharpness is "
            "claimed by the document or tested here)" % (neg, neg2))

    # stratum points: real points of X with contact set containing T_25
    exps = np.array([ex[:8] for ex, _ in tG], float)
    cofs = np.array([co for _, co in tG], float)

    def G_float(wv):
        mons = np.prod(np.power(np.maximum(wv, 0.0)[None, :], exps), axis=1)
        return float(cofs @ mons), float(np.abs(cofs) @ mons)

    pts = []
    for sdx, eps in ((25101, 0.02), (25102, 0.005), (25103, 0.1)):
        pts += K.numeric_stratum_search(T, sdx, 250, eps=eps)
    n_free = len(pts)
    # steer into the part of the stratum with d > 0: d is a linear function of the gap vector (its f54-coordinate in the
    # basis CANON; it is <= 0 on minor decompositions and on mirror points), so "d = target" is one more equation
    Bc = np.array([K.F_VALUES[jx] for jx in CANON], float).T
    drow = np.linalg.pinv(Bc)[5]
    mrow = np.linalg.pinv(np.array([K.F_VALUES[IMG[jx]] for jx in CANON], float).T)[5]
    n_forced = 0
    for kx, target in enumerate((0.005, 0.02, 0.05, 0.1, 0.2)):
        for row in (drow, mrow):
            got = K.numeric_stratum_search(T, 25200 + kx, 60, extra=lambda g, row=row, target=target: [float(row @ g) - target])
            n_forced += len(got)
            pts += got
    MIN8 = [jx for jx in GENS if K.F_IS_MINOR[jx]]
    MIRROR = [IMG[jx] for jx in CANON]                      # the canonical form of the mirror point g.x
    assert K.rank([K.F_VALUES[jx] for jx in MIN8]) == 8 and K.rank([K.F_VALUES[jx] for jx in MIRROR]) == 8
    bases = {"minor": MIN8, "canon": list(CANON), "mirror": MIRROR}
    Bm = {nm: np.array([K.F_VALUES[jx] for jx in bs], float).T for nm, bs in bases.items()}
    n_hard = n_minor_pts = n_unclass = 0
    r_cub = r_G = r_det = 0.0
    dmax = 0.0
    hard = []
    for pnt in pts:
        gaps = pnt["gaps"]
        sol = {}
        for nm in bases:
            w, res, _, _ = np.linalg.lstsq(Bm[nm], gaps, rcond=None)
            if np.max(np.abs(Bm[nm] @ w - gaps)) < 1e-8 and np.min(w) > -1e-8:
                sol[nm] = np.maximum(w, 0.0)
        if "minor" in sol:
            n_minor_pts += 1
            continue
        if "canon" in sol:
            wv = sol["canon"]
        elif "mirror" in sol:
            wv = sol["mirror"]                              # coefficients of g.x in the order of CANON
        else:
            n_unclass += 1
            continue
        dd = float(wv[5])
        n_hard += 1
        dmax = max(dmax, dd)
        hval = pnt["H"]
        qv, pv, jv, cvv, ev, dv, bv, av = wv
        Cs = (cvv + dv + qv) * (av + dv) * (bv + dv) + cvv * qv * (av + bv + dv) + ev * (cvv + dv) * (dv + qv)
        Qv = ev * jv - dv * (dv + ev) - pv * jv
        r_cub = max(r_cub, abs(-Cs + hval * Qv) / (abs(Cs) + abs(hval * Qv)))
        gval, gscale = G_float(wv)
        r_G = max(r_G, abs(gval) / gscale)
        u = np.array([float(x) for x in K.combo_coord(list(zip(CANON, [float(x) for x in wv])))])
        u[13] += hval
        Hm = np.array([[sum(cf * u[k] for cf, k in cell) for cell in row] for row in K.hess_linear_table()], float)
        sv = np.linalg.svd(Hm, compute_uv=False)
        r_det = max(r_det, sv[-1] / sv[0])
        hard.append((wv, hval))
    ck.check("NUMERIC", "every converged stratum point is, within 1e-8, a nonnegative combination of the 8 minor generators, "
             "or of the canonical generators (1), or of their mirror images (the trichotomy of the reduction)",
             n_unclass == 0, "%d unclassified" % n_unclass)
    ck.info("stratum search: 750 free starts (seeds 25101-25103): %d converged points of X with gaps >= 0 vanishing on T_25; "
            "600 starts with the extra equation d = target (targets 0.005 ... 0.2 for (1) and for its mirror; seeds "
            "25200-25204): %d converged.  Classification: %d minor decompositions, %d points of the form (1) or its mirror "
            "with d > 0 (d between %.3g and %.3g; normalisation: sum of gaps = 1)"
            % (n_free, n_forced, n_minor_pts, n_hard, min([float(w[5]) for w, _ in hard] or [0.0]), dmax))
    if n_hard:
        ck.check("NUMERIC", "at the %d canonical stratum points with d > 0: |-C_* + hQ| / (|C_*| + |hQ|) <= 1e-6, "
                 "|G| / (sum of |terms of G|) <= 1e-6 and smallest/largest singular value of Hess C(Re x) <= 1e-6 "
                 "(equations (2), (3) on actual points of X)" % n_hard,
                 r_cub <= 1e-6 and r_G <= 1e-6 and r_det <= 1e-6, "%.2e, %.2e, %.2e" % (r_cub, r_G, r_det))
        cap = {"SOC": 0, "old-H": 0, "four": 0, "(5)": 0, "(6)": 0, "double-SOC": 0, "none": 0}
        tol = 1e-7
        for wv, hval in hard:
            qv, pv, jv, cvv, ev, dv, bv, av = wv / wv[5]
            Kv = av * bv * (ev - 1) - (av + bv + ev + 1)
            hit = []
            if pv * jv >= 1 - tol:
                hit.append("SOC")
            if ev * (jv - 1) >= 1 - tol:
                hit.append("old-H")
            if Kv >= -tol:
                hit.append("four")
            if av * bv > 1 + tol and ev > 0 and jv >= 1 / ev - (av * bv - 1) / ((av + 1) * (bv + 1)) - tol:
                hit.append("(5)")
            wq = np.sqrt(av * bv)
            if av * bv < 1 - tol and abs(av - bv) > tol and ev > 0 and \
                    jv >= 1 / ev + (1 - wq) ** 2 / (np.sqrt(bv) - np.sqrt(av)) ** 2 - tol:
                hit.append("(6)")
            if abs(av * bv - 1) <= tol and ev * jv >= 1 - tol:
                hit.append("double-SOC")
            if not hit:
                cap["none"] += 1
            for hh in hit:
                cap[hh] += 1
        ck.check("NUMERIC", "every canonical stratum point with d > 0 satisfies at least one of the capacities (4), (5), (6) "
                 "(tolerance 1e-7 after scaling to d = 1)", cap["none"] == 0, "capacities met: %s" % cap)
    else:
        ck.info("no stratum point with d > 0 was found: equations (2), (3) are NOT tested on actual points by this search")

ck.finish()
