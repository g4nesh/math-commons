"""Independent exact replay of Case T_46 = {3,4,5,9,10,12} (group V) of D. Zheng's proof of the
Marcus-de Oliveira conjecture for n = 4.

Proof document replayed:
    MOC4_SUPP/proofs/manin_transfer_20261021/independent_moc4_T45_uniform_20261022.pdf
(6 pages; the document calls the set "T45").  Every polynomial is rebuilt here from the
definitions (the matrix P, the cubic C, the vectors f_j of the paper); nothing from the
supplement is executed or imported.  If MOC4_SUPP is set, three JSON files are read as data:
  * root_moc4_T45_full_radial_20261022.json      comparison only (stored C, Q, G_9);
  * root_moc4_T45_budget_PSD8_20261022.json and root_moc4_T45_cubic_reduction_20261022.json
    the data behind the four small Gram certificates of the document's p. 2 (a change of basis
    and sixteen cubic polynomials in quaternion coordinates).  The identity that makes them a
    certificate is re-derived and verified exactly here against Leibniz sums rebuilt from the
    paper.  Without MOC4_SUPP these four checks are skipped; the SOC and chain inequalities
    are then cited (PROOF.md 7.6) with the exact orbit checks of Section D.

Run:   MOC4_SUPP=/path/to/supplement  python replay_T46.py
Exit status 0 iff every check passes.  Checks tagged EXACT use integer/rational/symbolic
arithmetic only and prove what they state; checks tagged NUMERIC are falsification attempts
(floating point, or exact evaluation at random sample points) and are evidence only.
"""
import itertools
import math
import os
import random
import sys
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import replay_T44_common as M                                   # noqa: E402
from replay_T44_common import PR, RF, F_VALUES                   # noqa: E402

DIR = "proofs/manin_transfer_20261021/"
DOC = DIR + "independent_moc4_T45_uniform_20261022.pdf"
RADIAL_JSON = DIR + "root_moc4_T45_full_radial_20261022.json"
PSD_JSON = DIR + "root_moc4_T45_budget_PSD8_20261022.json"
CUB_JSON = DIR + "root_moc4_T45_cubic_reduction_20261022.json"
DOC_SHA256 = {
    DOC: "a1d4f501dc228423d268431fe06eb6086b63b25d44fef7136212016dfc8292f5",
    RADIAL_JSON: "2d75da111d7ff0771cc1bdd5c73ade1f8490e6185d5362b2a8c6fd13e838a90d",
    PSD_JSON: "4f8f3a5100da49159cd338338438b634d492768b410cf1034f2d59b6bc23e36a",
    CUB_JSON: "12189b86a0f359e313a6efc990b960909a8c6818875897492d3dc379b4f67dbd",
}

T = (3, 4, 5, 9, 10, 12)
ORDER = [0, 1, 2, 12, 59, 67, 70, 72]   # g = p f0 + r f1 + a f2 + b f12 + c f59 + q f67 + d f70 + s f72
NAMES = "p,r,a,b,c,q,d,s,h"
EXCL = (18, 16, 13, 8, 2, 14, 11, 7)    # (p,r,a,b,c,q,d,s) = (g18, g16, g13, g8, g2, g14, g11, g7)
SRC = ["p", "r", "s", "q", "a", "b", "c"]   # variable order of G_9 at d = 1 used for the tables
CH = ["M", "T", "U", "V", "A", "P", "Q"]    # chart variables

# the document's term counts, p. 4 (signs of x - 1, y - 1, z - 1 with x = ab, y = cr, z = cs)
DOC_COUNTS = {
    (-1, -1, -1): 7229, (-1, -1, 0): 2719, (-1, -1, 1): 13644, (-1, 0, -1): 2770, (-1, 0, 0): 1258,
    (-1, 0, 1): 20955, (-1, 1, -1): 15211, (-1, 1, 0): 21216, (-1, 1, 1): 289359,
    (0, -1, -1): 3052, (0, -1, 0): 1139, (0, -1, 1): 5596, (0, 0, -1): 1174, (0, 0, 0): 487,
    (0, 0, 1): 8205, (0, 1, -1): 6756, (0, 1, 0): 8458, (0, 1, 1): 116320,
    (1, -1, -1): 20860, (1, -1, 0): 7651, (1, -1, 1): 36870, (1, 0, -1): 8170, (1, 0, 0): 3052,
    (1, 0, 1): 54718, (1, 1, 0): 36336, (1, 1, 1): 499377,
}
# p. 5: term counts of the Bernstein coefficients and the monomial in (T,U,V) named for the last one
DOC_BERN = {
    ("PB", 0): ([0, 3736, 10783, 14484, 16612, 17692], (6, 6, 2), 10),
    ("PB", 1): ([792, 6108, 9995, 12397, 13756], (6, 5, 2), 36),
    ("PR", 0): ([4140, 7956, 11076, 13139, 13756], (3, 2, 5), 36),
    ("PR", 1): ([5400, 8588, 11282, 12732, 12912], (3, 3, 4), 12),
}


def pmul(v, w):
    return (v[0] * w[0], v[1] * w[1])


def pdiv(v, w):
    return (v[0] * w[1], v[1] * w[0])


def base_chart(Rc, xs, ys, zs, pbound=None):
    """Substitutions (11)-(13) of the document for the sign triple (xs, ys, zs), as pairs
    (numerator, denominator) over Z[M,T,U,V,A,P,Q] in the order SRC = (p, r, s, q, a, b, c).
    The representation of each fraction (which common factors are kept) is that of the
    document's verifier, learnt by reading it as text, so that term counts can be compared."""
    Mv, T_, U, V, A, P, Q = Rc.gens
    one = Rc.one
    xv = (U, U + 1) if xs < 0 else (one, one) if xs == 0 else (U + 1, one)
    yv = (V, V + 1) if ys < 0 else (one, one) if ys == 0 else (V + 1, one)
    zv = (A, A + 1) if zs < 0 else (one, one) if zs == 0 else (A + 1, one)
    cv = (Mv, one)
    av = pdiv((T_, T_ + 1), cv)                    # a = h / c,  h = ac = T/(1+T)
    bv = pdiv(xv, av)                              # b = x / a
    rv = pdiv(yv, cv)                              # r = y / c
    sv = pdiv(zv, cv)                              # s = z / c
    if pbound == "PR" or (pbound is None and xs != 1 and ys == 1):
        pv, pkind = pmul(rv, (P, V * (P + 1))), "PR"      # p = (r/(y-1)) P/(1+P)
    elif xs == 1:
        pv, pkind = pmul(bv, (P, U * (P + 1))), "PB"      # p = (b/(x-1)) P/(1+P)
    else:
        pv, pkind = (P, one), "free"
    qv = pmul(sv, (Q, A * (Q + 1))) if zs == 1 else (Q, one)
    return dict(vals=[pv, rv, sv, qv, av, bv, cv], x=xv, y=yv, z=zv, pkind=pkind)


def stage_H(Rc, out, yv, zv):
    """Substitution (14), done as in the document's proof of Lemma 1: divide by the lowest power
    M^m0 of M, then put M = nn/dd with nn = (y+z) M, dd = (y+z-1)(1+M) (numerators over the
    common denominator of y and z) and multiply by dd^(degree in M).  Returns
    (polynomial, m0, degree, nn, dd)."""
    Mv = Rc["M"]
    d = PR.terms(out)
    m0 = min(ex[0] for ex in d)
    parts = {}
    for ex, cf in d.items():
        parts.setdefault(ex[0] - m0, {})[(0,) + ex[1:]] = cf
    dM = max(parts)
    yn = yv[0] * zv[1] + zv[0] * yv[1]             # numerator of y + z over y_d z_d
    yd = yn - yv[1] * zv[1]                        # numerator of y + z - 1
    nn, dd = yn * Mv, yd * (Mv + 1)
    tot = Rc.zero
    for i in sorted(parts):
        tot = tot + Rc.from_terms(parts[i]) * nn ** i * dd ** (dM - i)
    return tot, m0, dM, nn, dd


def strict_monomial(Rc, pol, forbidden):
    """The smallest (total degree, exponents) monomial with a positive coefficient that
    involves none of the variables in `forbidden`; None if there is none."""
    fb = [Rc.idx[nm] for nm in forbidden]
    best = None
    for i in range(len(pol)):
        if pol.coefficient(i) > 0:
            ex = tuple(int(e) for e in pol.exponent_tuple(i))
            if all(ex[v] == 0 for v in fb):
                key = (sum(ex), ex)
                if best is None or key < best[0]:
                    best = (key, int(pol.coefficient(i)))
    return None if best is None else (best[0][1], best[1])


def hard(p, r, a, b, c, q, s):
    """The region (10) of the document at d = 1 (all five capacities (2)-(6) fail), for
    a, c, r > 0 and p, b, q, s >= 0."""
    return (a > 0 and c > 0 and r > 0 and min(p, b, q, s) >= 0
            and a * c < 1 and c * p * r < p + r and c * (r + s) < r + s + 1
            and (a * p * b < p + b or (p == 0 and b == 0))
            and (c * q * s < q + s or (q == 0 and s == 0)))


def main():
    import sympy as sp
    import numpy
    import scipy
    import flint
    ck = M.Checker("Case T_46 = {3,4,5,9,10,12} (document 'T45'): independent exact replay")
    ck.info("python %s, sympy %s, python-flint %s (integer polynomials of the chart tables and of the Gram "
            "certificates; integer determinants at test points), numpy %s (int64 in det_linear_forms; floats in "
            "NUMERIC checks), scipy %s (NUMERIC only)"
            % (sys.version.split()[0], sp.__version__, flint.__version__, numpy.__version__, scipy.__version__))
    root = os.environ.get("MOC4_SUPP")

    # ------------------------------------------------------------------------------------
    ck.section("A. Provenance and foundations")
    for rel in (DOC, RADIAL_JSON, PSD_JSON, CUB_JSON):
        M.check_document(ck, rel, DOC_SHA256[rel])
    M.check_foundations(ck)
    M.check_radial_necessity(ck)
    M.crosscheck_supplement_tables(ck)
    M.crosscheck_audit_library(ck, [46, 44, 7, 15, 34])

    # ------------------------------------------------------------------------------------
    ck.section("B. Contact set (document p. 1)")
    ck.check("EXACT", "the document's set {3,4,5,9,10,12} is the paper's T_46 (Appendix E) as a set of "
             "lexicographic indices: group element = identity; the document's number is the paper's minus one",
             tuple(sorted(T)) == M.PAPER_T_V[46])
    rep = M.orbit_rep(T)
    stab = [t for t in M.G_MAPS if set(t[s] for s in T) == set(T)]
    ck.check("EXACT", "orbit representative of T_46 under G is (0,1,2,7,8,17); T_44, T_34, T_7, T_15 lie in other "
             "orbits", rep == (0, 1, 2, 7, 8, 17)
             and all(M.orbit_rep(M.PAPER_T_V[k]) != rep for k in (44, 34, 7, 15)),
             "T_46 -> %s; stabiliser of T_46 in G has order %d" % (rep, len(stab)))

    # ------------------------------------------------------------------------------------
    ck.section("C. The contact cone and the coordinates (1) (document Section 1)")
    gens = M.generators(T)
    ck.check("EXACT", "the f_j of Appendix D vanishing on T are exactly f0, f1, f2, f12, f59, f67, f70, f72",
             gens == sorted(ORDER), "computed %s" % gens)
    ck.check("EXACT", "f70 = h_{4231,+} is the only non-minor generator; the seven others are minor vectors",
             [j for j in gens if not M.F_IS_MINOR[j]] == [70],
             ", ".join("f%d = %s" % (j, M.f_name(j)) for j in ORDER))
    null = M.nullspace([M.P[s] for s in T], 14)
    ck.check("EXACT", "the six contact equations have rank 6 on V_R; their solution space has dimension 8",
             M.rank([M.P[s] for s in T]) == 6 and len(null) == 8)
    ident = [[F_VALUES[j][s] for j in ORDER] for s in EXCL]
    ck.check("EXACT", "the values of the eight generators at the permutation indices (18,16,13,8,2,14,11,7) form the "
             "identity matrix: the eight generators are a basis of the solution space and "
             "(p,r,a,b,c,q,d,s) = (g18, g16, g13, g8, g2, g14, g11, g7), all gaps",
             ident == [[int(i == k) for k in range(8)] for i in range(8)]
             and M.rank([M.F_COORD[j] for j in ORDER]) == 8)
    d_, rays = M.extreme_rays_of_contact_cone(T)
    ck.check("EXACT", "direct enumeration of the extreme rays of N_T = {Pc >= 0, (Pc)_T = 0} (dimension %d; no use of "
             "the 82-facet list) gives exactly these eight vectors: N_T is the simplicial cone on them" % d_,
             rays == sorted(F_VALUES[j] for j in ORDER), "%d rays" % len(rays))
    M.check_minor_leibniz(ck, [j for j in ORDER if j != 70])
    ck.info("PROOF-READ: so c(x) = p f0 + r f1 + a f2 + b f12 + c f59 + q f67 + d f70 + s f72 with all eight "
            "coefficients >= 0; if d = 0 it is a nonnegative combination of minor vectors and Phi >= 0.")

    # ------------------------------------------------------------------------------------
    ck.section("D. The capacities (2)-(6) (document Section 1, p. 2)")
    rng = random.Random(20261127)
    M.check_transport(ck, "SOC", (2, 59, 70), rng)
    M.check_transport(ck, "chain", (2, 0, 12, 70), rng)
    M.check_transport(ck, "chain", (59, 0, 1, 70), rng)
    M.check_transport(ck, "chain", (59, 67, 72, 70), rng)
    M.check_transport(ck, "old-H", (1, 59, 70), rng)
    M.check_transport(ck, "old-H", (72, 59, 70), rng)
    # the document's two stated actions for old-H (profile action sigma -> col o sigma o row)
    fidx = {v: i for i, v in enumerate(F_VALUES)}
    bases = []
    for target, row, col in (((1, 59, 70), (0, 3, 1, 2), (3, 1, 0, 2)), ((72, 59, 70), (1, 2, 0, 3), (0, 2, 3, 1))):
        act = [M.PIDX[tuple(col[s[row[i]]] for i in range(4))] for s in M.PERMS]
        base = []
        for j in target:
            w = [None] * 24
            for i in range(24):
                w[act[i]] = F_VALUES[j][i]              # f_target(sigma) = f_base(col o sigma o row)
            base.append(fidx.get(tuple(w)))
        bases.append(tuple(base))
    ck.check("EXACT", "the document's two stated old-H actions (rows (0,3,1,2), (1,2,0,3); columns (3,1,0,2), (0,2,3,1); "
             "f_image(sigma) = f_base(col o sigma o row)) send ONE base triple (m_A, m_B, H) to (f1, f59, f70) and to "
             "(f72, f59, f70), and that base triple lies in the orbit of the representative (f12, f2; f70) of "
             "PROOF.md 7.6 with the same roles",
             bases[0] == bases[1] and None not in bases[0] and bool(M.find_transport((12, 2, 70), bases[0])),
             "base triple (m_A, m_B, H) = (f%s, f%s, f%s)" % bases[0])
    p, r, a, b, c, q, d, s, t = sp.symbols("p r a b c q d s t", positive=True)
    d1, d2 = d * r / (r + s), d * s / (r + s)
    ck.check("EXACT", "(6) from two old-H atoms: split d = d1 + d2 in the ratio r : s; r f1 + c1 f59 + d1 f70 and "
             "s f72 + c2 f59 + d2 f70 satisfy t(c - d) >= d^2 with equality for c1 = d1 + d1^2/r, c2 = d2 + d2^2/s, "
             "and c1 + c2 = d + d^2/(r+s); for s = 0 the first atom alone gives c >= d + d^2/r",
             sp.cancel(d1 + d1 ** 2 / r + d2 + d2 ** 2 / s - (d + d ** 2 / (r + s))) == 0
             and sp.cancel(r * ((d1 + d1 ** 2 / r) - d1) - d1 ** 2) == 0
             and sp.cancel(s * ((d2 + d2 ** 2 / s) - d2) - d2 ** 2) == 0)
    gram_ok = None
    if root and os.path.isfile(os.path.join(root, PSD_JSON)) and os.path.isfile(os.path.join(root, CUB_JSON)):
        gf = M.GramFace(os.path.join(root, DIR), "root_moc4_T45")
        ck.check("EXACT", "certificate data (JSON): source order (0,1,2,12,59,67,70,72); the four signed coordinate "
                 "maps permute p and q separately (so they preserve |p|^2, |q|^2); quaternion model M(p,q) = matrix "
                 "of x -> p x q: M^t M = |p|^2 |q|^2 I and det M = (|p|^2 |q|^2)^2 as polynomial identities",
                 gf.order == ORDER and gf.signed_perms_ok and gf.model_ok)
        key = {kk: i for i, kk in enumerate(gf.keys)}

        def unit(i):
            v = [Fraction(0)] * 36
            v[key[(i, i)]] = Fraction(1)
            return v
        ck.check("EXACT", "own linear algebra on (*) (see GramFace; %d coefficient equations, rank %d): the pairs "
                 "(Gram matrix, weights) with (|p|^2+|q|^2) Sigma(Gram) = 4 Phi_w(M(p,q)) mod (|p|^2 - |q|^2) are "
                 "exactly: ten independent linear relations on the Gram matrix, and weights = its diagonal "
                 "(G00,...,G77) in the order (p,r,a,b,c,q,d,s); the ten relations stored in the supplement span "
                 "the same space" % (gf.n_equations, gf.rank),
                 len(gf.relations) == 10 and all(gf.weight_rows[m] == unit(m) for m in range(8))
                 and M.rank(gf.stored_R) == 10 and M.rank(gf.stored_R + gf.relations) == 10)
        ck.info("relations derived here (Gik = Gram entries, vectors in the order (f0,f1,f2,f12,f59,f67,f70,f72)):")
        for rel in gf.relations:
            ck.info("   " + " ".join("%+d*G%d%d" % (int(cf), gf.keys[i][0], gf.keys[i][1])
                                     for i, cf in enumerate(rel) if cf) + " = 0")
        Z = sp.Integer(0)
        Fr = Fraction
        gram_ok = True
        gram_ok &= M.sym_gram_check(
            ck, gf, "(2) SOC: Gram on (a, c, d) with G24 = d, c = d^2/a + t", {(2, 4): d},
            [Z, Z, a, Z, d * d / a + t, Z, d, Z], [Z, Z, a, Z, d * d / a + t, Z, d, Z],
            [2, 4, 6], [a, t, d], {a: Fr(3, 2), d: Fr(2), t: Fr(1, 3)}, note="  [PSD iff ac >= d^2]")
        gram_ok &= M.sym_gram_check(
            ck, gf, "(3) PB chain: Gram on (p, a, b, d) with G02 = G23 = -d, a = d^2/p + d^2/b + t",
            {(0, 2): -d, (2, 3): -d},
            [p, Z, d * d / p + d * d / b + t, b, Z, Z, d, Z], [p, Z, d * d / p + d * d / b + t, b, Z, Z, d, Z],
            [0, 3, 2, 6], [p, b, t, d], {p: Fr(3, 2), b: Fr(1, 2), d: Fr(2), t: Fr(1, 3)},
            note="  [PSD iff a >= d^2/p + d^2/b]")
        gram_ok &= M.sym_gram_check(
            ck, gf, "(4) PR chain: Gram on (p, r, c, d) with G04 = d, G14 = -d, c = d^2/p + d^2/r + t",
            {(0, 4): d, (1, 4): -d},
            [p, r, Z, Z, d * d / p + d * d / r + t, Z, d, Z], [p, r, Z, Z, d * d / p + d * d / r + t, Z, d, Z],
            [0, 1, 4, 6], [p, r, t, d], {p: Fr(3, 2), r: Fr(1, 2), d: Fr(2), t: Fr(1, 3)},
            note="  [PSD iff c >= d^2/p + d^2/r]")
        gram_ok &= M.sym_gram_check(
            ck, gf, "(5) QS chain: Gram on (c, q, d, s) with G45 = d, G47 = -d, c = d^2/q + d^2/s + t",
            {(4, 5): d, (4, 7): -d},
            [Z, Z, Z, Z, d * d / q + d * d / s + t, q, d, s], [Z, Z, Z, Z, d * d / q + d * d / s + t, q, d, s],
            [5, 7, 4, 6], [q, s, t, d], {q: Fr(3, 2), s: Fr(1, 2), d: Fr(2), t: Fr(1, 3)},
            note="  [PSD iff c >= d^2/q + d^2/s]")
        ck.info("PROOF-READ + CITED (every element of SO(4) is x -> p x q for unit quaternions p, q): on unit "
                "quaternions (*) reads Sigma(Gram) = 2 Phi_w(M(p,q)); the pivots are >= 0 for t >= 0 and positive "
                "p, r, a, b, q, s, d, so each Gram matrix is positive semidefinite and Phi_w >= 0 on SO(4); the "
                "weights are homogeneous in (coefficients, d); adding the remaining minor terms keeps Phi >= 0.  "
                "So (2)-(5) hold on this face without citing PROOF.md 7.6.")
    else:
        ck.info("MOC4_SUPP not set (or the two certificate files not found): the four Gram certificates of p. 2 are "
                "not replayed; (2)-(5) are then cited from PROOF.md 7.6 with the orbit checks above")
    ck.info("CITED: old-H (PROOF.md 7.6): t f12 + c f2 + d f70 has Phi >= 0 on SO(4) if t(c - d) >= d^2, on the whole "
            "orbit; used for (f1, f59; f70) and (f72, f59; f70).  SOC and chain (PROOF.md 7.6) as stated in "
            "replay_T44_common.UNIVERSAL%s." % ("" if gram_ok is None else "; here also certified by the Gram matrices"))

    # ------------------------------------------------------------------------------------
    ck.section("E. The three closed zero faces a = 0, c = 0, r = 0 (document p. 2)")
    okz = True
    det = []
    for nm, extra, nu, g, image in (("a", 13, 34, ((1, 3, 2, 0), (1, 3, 2, 0), 1), (3, 5, 9, 10, 12, 13)),
                                    ("c", 2, 7, ((0, 1, 2, 3), (0, 3, 2, 1), 0), (2, 3, 4, 5, 9, 10)),
                                    ("r", 16, 15, ((3, 2, 0, 1), (3, 0, 1, 2), 0), (3, 4, 5, 9, 10, 16))):
        gm = M.g_index_map(g)
        img = tuple(sorted(gm[s_] for s_ in M.PAPER_T_V[nu]))
        okz &= EXCL[NAMES.split(",").index(nm)] == extra            # the coefficient is the gap at `extra`
        okz &= img == image and set(img) <= set(T) | {extra} and M.generic_covariance(g)
        n_el = sum(1 for g2, gm2 in M.G_EXPLICIT if set(gm2[s_] for s_ in M.PAPER_T_V[nu]) <= set(T) | {extra})
        det.append("%s = g%d = 0: g = %s maps T_%d = %s to %s (%d group elements do)" % (nm, extra, g, nu,
                                                                                         M.PAPER_T_V[nu], img, n_el))
    ck.check("EXACT", "a = g13, c = g2, r = g16; the document's three actions (sigma -> row o sigma o col, with "
             "sigma^-1 in the transposed case) map the paper's T_34, T_7, T_15 (the document's 'T33', 'T6', 'T14') "
             "into T u {13}, T u {2}, T u {16} with the stated images; covariance of Phi under each element "
             "(monomial identity for a generic matrix)", okz, "; ".join(det))
    ck.info("PROOF-READ: if a = 0 then c(x) vanishes on g.T_34, so g^-1.x is a point of X with gap vector in "
            "N_{T_34} and Phi_{c(x)}(O) = Phi_{g^-1.c(x)}(O') >= 0 by Case T_34 (PROOF.md 6.6); likewise c = 0 by "
            "Case T_7 and r = 0 by Case T_15.  CITED cases: T_34, T_7, T_15.")

    # ------------------------------------------------------------------------------------
    ck.section("F. Radial necessity in the coordinates (document Section 2, eqs. (7)-(9))")
    rd = M.radial_data(ck, ORDER, NAMES)
    R9, C0, Q, G = rd["R"], rd["C0"], rd["Q"], rd["G"]
    R8 = PR(NAMES.split(","))
    vp, vr, va, vb, vc, vq, vd, vs, vh = R8.gens
    C0f = R8.from_terms(M.packed_terms(R9, C0, R9.names))
    Qf = R8.from_terms(M.packed_terms(R9, Q, R9.names))
    T0 = (vb + vd) * (vd + vr) + vd * vp
    Rs_ = vr + vs
    ck.check("EXACT", "(7): C(g + h kappa) = C0 + hQ with C0 = c r s - (d+q+s)[(b+d)(d+r) + dp] and "
             "Q = (b-c)a + c(r+s) - (b+d)(d+r+s) (the document's cubic is the paper's C, same sign)",
             C0f == vc * vr * vs - (vd + vq + vs) * T0
             and Qf == (vb - vc) * va + vc * Rs_ - (vb + vd) * (vd + vr + vs))
    Ppos = vd ** 3 * Rs_ + vd * vd * (vr * vr + vr * vs + vs * vs) + Rs_ * (vd + vs) * (vb * (vd + vr) + vd * vp) \
        + Rs_ * vq * T0
    st = M.fl_stats(Ppos)
    ck.check("EXACT", "(8): (r+s) C0 = (c(r+s) - d(r+s) - d^2) r s - Ppos, where Ppos = d^3 R + d^2 (r^2 + rs + s^2) + "
             "R(d+s)[b(d+r) + dp] + R q T0 (R = r + s) has only positive coefficients and contains d^3 r",
             Rs_ * C0f == (vc * Rs_ - vd * Rs_ - vd * vd) * vr * vs - Ppos and st[1] == 0 and st[2] >= 1
             and PR.terms(Ppos).get((0, 1, 0, 0, 0, 0, 3, 0, 0)) == 1, "%d terms" % st[0])
    ck.info("PROOF-READ: for d, r > 0, Ppos >= d^3 r > 0; if (6) fails then c(r+s) - d(r+s) - d^2 < 0, so "
            "(r+s) C0 < 0 and C0 < 0; then C0 + hQ = 0 forces Q != 0 and h = -C0/Q.")
    ck.check("EXACT", "G has 1158 terms (the document's count for G_9), is homogeneous of degree 9, and its d^9 "
             "coefficient is +1", len(G) == 1158 and all(sum(R9.exps(k_)) == 9 for k_ in G)
             and G.get(R9.key((0, 0, 0, 0, 0, 0, 9, 0, 0))) == 1)
    M.compare_radial_json(ck, rd, RADIAL_JSON, alias={"sa": "s"})
    ck.check("EXACT", "observation (not used by the proof, which has r, d > 0 at this point): G_9 vanishes identically "
             "on {r = d = 0}", not R9.subs_zero(G, ["r", "d"]),
             "terms left at d = 0: %d; at r = 0: %d" % (len(R9.subs_zero(G, ["d"])), len(R9.subs_zero(G, ["r"]))))
    G1 = R9.subs_one(G, "d")
    Gt = M.packed_terms(R9, G1, SRC)
    degs = [max(ex[v] for ex in Gt) for v in range(7)]
    ck.check("EXACT", "at d = 1, G_9 is a polynomial in (p,r,s,q,a,b,c) with 1158 terms, constant term 1 and degrees "
             "(3,2,2,3,2,2,2)", len(Gt) == 1158 and Gt.get((0,) * 7) == 1 and degs == [3, 2, 2, 3, 2, 2, 2])

    # ------------------------------------------------------------------------------------
    ck.section("G. The hard region (10) and the product charts (11)-(14) (document Section 3)")
    x_, y_, z_ = sp.symbols("x y z", positive=True)
    ck.check("EXACT", "with x = ab, y = cr, z = cs (a, c, r > 0, d = 1): a p b < p + b <=> p (x-1) < b; "
             "c p r < p + r <=> p (y-1) < r; c q s < q + s <=> q (z-1) < s; and c(r+s) < r + s + 1 <=> "
             "c (y+z-1) < y + z (these are (13) and the bound before (14))",
             sp.expand((p + b - a * p * b) - (b - p * (a * b - 1))) == 0
             and sp.expand((p + r - c * p * r) - (r - p * (c * r - 1))) == 0
             and sp.expand((q + s - c * q * s) - (s - q * (c * s - 1))) == 0
             and sp.cancel(((y_ + z_) / c + 1 - (y_ + z_)) * c - ((y_ + z_) - c * (y_ + z_ - 1))) == 0)
    Rc = PR(CH)
    cM, cT, cU, cV, cA, cP, cQ = Rc.gens
    rf = lambda pair: RF(Rc, pair[0], pair[1])
    rM, rT, rU, rV, rA, rP, rQ = [RF(Rc, g_) for g_ in Rc.gens]
    okch = True
    okden = True
    for xs, ys, zs in itertools.product([-1, 0, 1], repeat=3):
        for pb in ((None,) if (xs, ys, zs) != (1, 1, -1) else ("PB", "PR")):
            ch = base_chart(Rc, xs, ys, zs, pb)
            posvars = ["M", "T", "V"] + (["U"] if xs == 1 else []) + (["A"] if zs == 1 else [])
            okden &= all(M.fl_nonneg(v[0]) and M.fl_positive_on(Rc, v[1], posvars) for v in ch["vals"])
            if ys >= 0 and zs >= 0:
                yn = ch["y"][0] * ch["z"][1] + ch["z"][0] * ch["y"][1]
                okden &= M.fl_positive_on(Rc, yn * cM, posvars) \
                    and M.fl_positive_on(Rc, (yn - ch["y"][1] * ch["z"][1]) * (cM + 1), posvars)
            vp_, vr_, vs_, vq_, va_, vb_, vc_ = [rf(v) for v in ch["vals"]]
            xv, yv, zv = rf(ch["x"]), rf(ch["y"]), rf(ch["z"])
            okch &= va_ * vc_ == rT / (1 + rT) and va_ * vb_ == xv and vc_ * vr_ == yv and vc_ * vs_ == zv
            okch &= xv == {-1: rU / (1 + rU), 0: RF(Rc, 1), 1: 1 + rU}[xs]
            okch &= yv == {-1: rV / (1 + rV), 0: RF(Rc, 1), 1: 1 + rV}[ys]
            okch &= zv == {-1: rA / (1 + rA), 0: RF(Rc, 1), 1: 1 + rA}[zs]
            want = "PR" if pb == "PR" else "PB" if (pb == "PB" or xs == 1) else "PR" if ys == 1 else "free"
            okch &= ch["pkind"] == want
            if want == "PB":
                okch &= vp_ == vb_ / (xv - 1) * rP / (1 + rP)
            elif want == "PR":
                okch &= vp_ == vr_ / (yv - 1) * rP / (1 + rP)
            else:
                okch &= vp_ == rP
            okch &= (vq_ == vs_ / (zv - 1) * rQ / (1 + rQ)) if zs == 1 else (vq_ == rQ)
    ck.check("EXACT", "all 27 sign triples (28 charts: two for (+,+,-)): c = M, ac = T/(1+T), ab = x, cr = y, cs = z "
             "with x, y, z as in (11); p = (b/(x-1)) P/(1+P) if x > 1 (PB bound), else p = (r/(y-1)) P/(1+P) if "
             "y > 1 (PR bound), else p = P; q = (s/(z-1)) Q/(1+Q) if z > 1, else q = Q (identities between rational "
             "functions)", okch)
    ck.check("EXACT", "in all 28 charts every numerator has nonnegative coefficients and every denominator, and both "
             "polynomials (y+z) M and (y+z-1)(1+M) of (14) where it is used, has nonnegative coefficients and a "
             "monomial in the variables that are > 0 on the domain (M, T, V; U if x > 1; A if z > 1): all chart "
             "values are >= 0 and all clearing multipliers are > 0 on the domain", okden)
    ck.info("PROOF-READ: the charts are ONTO the hard region.  Given a hard point: T = ac/(1-ac) > 0; U = x/(1-x) >= 0 "
            "or x - 1 > 0; V = y/(1-y) > 0 or y - 1 > 0 (y > 0 because c, r > 0); A likewise from z >= 0; "
            "P = p/(cap - p) >= 0 with the cap of (13), or P = p; Q likewise; M = c > 0, resp. after (14) "
            "M = c (y+z-1)/((y+z) - c (y+z-1)) > 0.  Where a bound of (10) is not used by a chart (p for x > 1 and "
            "y > 1 in the direct domains, c unless y, z >= 1) the chart covers a larger region.")

    # ------------------------------------------------------------------------------------
    ck.section("H. Lemma 1: the 26 direct product domains (document pp. 4)")
    samples = []          # (label, polynomial, chart, staging) for the bookkeeping test
    nstrict = 0
    for xs, ys, zs in itertools.product([-1, 0, 1], repeat=3):
        if (xs, ys, zs) == (1, 1, -1):
            continue
        ch = base_chart(Rc, xs, ys, zs)
        out = Rc.compose(Gt, [v[0] for v in ch["vals"]], [v[1] for v in ch["vals"]], degs)
        stg = None
        if ys >= 0 and zs >= 0:
            out, m0, dM, nn, dd = stage_H(Rc, out, ch["y"], ch["z"])
            stg = (m0, dM, nn, dd)
        st = M.fl_stats(out)
        forb = ["P", "Q"] + (["U"] if xs < 0 else []) + (["A"] if zs < 0 else [])
        mono = strict_monomial(Rc, out, forb)
        nstrict += mono is not None
        sg = "".join("-0+"[v + 1] for v in (xs, ys, zs))
        ck.check("EXACT", "domain (%s) [signs of x-1, y-1, z-1]%s: the cleared polynomial has only nonnegative "
                 "coefficients, %d terms (document: %d), and a monomial with positive coefficient free of %s"
                 % (sg, " with (14)" if stg else "", st[0], DOC_COUNTS[(xs, ys, zs)], ", ".join(forb)),
                 st[1] == 0 and st[0] == DOC_COUNTS[(xs, ys, zs)] and mono is not None,
                 "smallest coefficient %d; strict monomial %s" % (
                     st[2], None if mono is None else "%d * %s" % (mono[1], "*".join(
                         "%s^%d" % (nm, e) for nm, e in zip(CH, mono[0]) if e) or "1")))
        samples.append(((xs, ys, zs, None), out, ch, stg))
    ck.check("EXACT", "26 direct domains in all, each with a strict monomial", len(samples) == 26 and nstrict == 26)
    ck.info("PROOF-READ: in each domain the variables that may vanish are P, Q (p, q = 0), U when x < 1 (b = 0) and A "
            "when z < 1 (s = 0); M, T, V are > 0, and U resp. A are > 0 when x > 1 resp. z > 1.  A polynomial with "
            "nonnegative coefficients and a positive monomial in the strictly positive variables is > 0 there; the "
            "clearing multipliers are > 0 (Section G).  So G_9 > 0 on the hard region for the 26 sign triples.")

    # ------------------------------------------------------------------------------------
    ck.section("I. The remaining domain x > 1, y > 1, z < 1: the square interval (document Section 3.1)")
    Nn = cT * cU * (cV + 1)
    Dd = (cT + 1) * (cU + 1) * cV
    chB = base_chart(Rc, 1, 1, -1, "PB")
    vb_, vr_ = rf(chB["vals"][5]), rf(chB["vals"][1])
    ck.check("EXACT", "(15): r/V - b/U = (N - D M^2)/(M T U V) with N = T U (1+V), D = (1+T)(1+U) V; so the PB bound "
             "b/U is the smaller one iff M^2 <= L = N/D", (vr_ / rV - vb_ / rU) == RF(Rc, Nn - Dd * cM * cM, cM * cT * cU * cV))
    tt, lam = sp.symbols("t lambda")
    okb = True
    for n in (4, 5):
        aa = sp.symbols("a0:%d" % (n + 1))
        bern = sum(sum(aa[i] * sp.Rational(math.comb(k_, i), math.comb(n, i)) * lam ** i for i in range(k_ + 1))
                   * math.comb(n, k_) * tt ** k_ * (1 - tt) ** (n - k_) for k_ in range(n + 1))
        okb &= sp.expand(bern - sum(aa[i] * (lam * tt) ** i for i in range(n + 1))) == 0
    ck.check("EXACT", "(16): for A(Z) = sum a_i Z^i of degree n = 4, 5 and Z = t lambda: A(Z) = sum_k b_k C(n,k) t^k "
             "(1-t)^(n-k) with b_k = sum_{i<=k} a_i C(k,i)/C(n,i) lambda^i (symbolic identity); with lambda = N*/D* "
             "the document's B_k is m D*^n b_k, m = lcm of the binomial coefficients", okb)
    for bound in ("PB", "PR"):
        ch = base_chart(Rc, 1, 1, -1, bound)
        out = Rc.compose(Gt, [v[0] for v in ch["vals"]], [v[1] for v in ch["vals"]], degs)
        samples.append(((1, 1, -1, bound), out, ch, None))
        dct = PR.terms(out)
        lo, hi = min(ex[0] for ex in dct), max(ex[0] for ex in dct)
        nneg = M.fl_stats(out)[1]
        ck.check("EXACT", "P_%s(M): the polynomial obtained from (11)-(13) with the %s bound on p (without (14)) has "
                 "degree %d in M (document: %d), lowest power M^%d, and it is NOT coefficientwise nonnegative "
                 "(%d negative coefficients): this domain needs the interval argument"
                 % (bound, bound, hi, 10 if bound == "PB" else 12, lo, nneg),
                 hi == (10 if bound == "PB" else 12) and nneg > 0, "%d terms" % len(dct))
        if bound == "PR":
            dct = {(hi - ex[0],) + ex[1:]: cf for ex, cf in dct.items()}      # w^12 P_PR(1/w), w = 1/M
            nn, dd = Dd, Nn
        else:
            nn, dd = Nn, Dd
        for parity in (0, 1):
            terms = {ex: cf for ex, cf in dct.items() if ex[0] % 2 == parity}
            n = max(ex[0] // 2 for ex in terms)
            coeff = [Rc.from_terms({(0,) + ex[1:]: cf for ex, cf in terms.items() if ex[0] // 2 == i})
                     for i in range(n + 1)]
            mult = 1
            for i in range(n + 1):
                mult = mult * math.comb(n, i) // math.gcd(mult, math.comb(n, i))
            cnts, neg = [], 0
            last = None
            for k_ in range(n + 1):
                pol = Rc.zero
                for i in range(k_ + 1):
                    pol = pol + coeff[i] * nn ** i * dd ** (n - i) * (mult * math.comb(k_, i) // math.comb(n, i))
                st = M.fl_stats(pol)
                cnts.append(st[0])
                neg += st[1]
                last = pol
            want_cnt, mono, mcoef = DOC_BERN[(bound, parity)]
            got = PR.terms(last).get((0,) + mono + (0, 0, 0))
            ck.check("EXACT", "%s bound, %s part in %s (degree %d in the squared variable, interval [0, %s], m = %d): "
                     "all %d Bernstein coefficients (16) have only nonnegative coefficients; term counts %s = the "
                     "document's; the last one contains %d T^%d U^%d V^%d, a positive monomial in T, U, V alone, as "
                     "the document says"
                     % (bound, "even" if parity == 0 else "odd", "M" if bound == "PB" else "w = 1/M", n,
                        "N/D" if bound == "PB" else "D/N", mult, n + 1, cnts, mcoef, mono[0], mono[1], mono[2]),
                     neg == 0 and cnts == want_cnt and got == mcoef)
    ck.info("PROOF-READ: for 0 < M^2 <= L put t = M^2/L in (0,1]: D^5 E(M^2) = sum_k (B_k/m) C(5,k) t^k (1-t)^(5-k) >= "
            "(B_5/m) t^5 > 0, since all B_k >= 0 and B_5 > 0 (T, U, V > 0); likewise O(M^2) > 0; so "
            "P_PB(M) = E(M^2) + M O(M^2) > 0.  For M^2 >= L the same with w = 1/M, w^2 in (0, D/N].  Both include "
            "M^2 = L.  With the positive clearing multipliers this gives G_9 > 0 on the hard region for (+,+,-).")

    # ------------------------------------------------------------------------------------
    bookkeeping_test(ck, Rc, Gt, degs, samples)
    hard_region_test(ck, Gt)
    numeric_part(ck, rd, Gt)
    ck.finish()


def bookkeeping_test(ck, Rc, Gt, degs, samples):
    """Exact evaluation at random chart points: (polynomial of the table) = G_9(original
    coefficients given by the chart) * (the explicit positive multiplier), and the original
    coefficients lie in the region where the bounds used by the chart hold."""
    ck.section("J. Exact evaluation tests (falsification attempts with exact arithmetic)")
    rng = random.Random(4601)
    bad = npts = 0
    for (xs, ys, zs, bound), out, ch, stg in samples:
        dct = PR.terms(out)
        for _ in range(2):
            th = {}
            for nm in CH:
                zero_ok = nm in ("P", "Q") or (nm == "U" and xs < 0) or (nm == "A" and zs < 0)
                th[nm] = Fraction(0) if (zero_ok and rng.random() < 0.25) else Fraction(rng.randint(1, 6), rng.randint(1, 3))
            pt = [th[nm] for nm in CH]
            val = M.eval_terms_common([dct], pt)[0]
            # the chart, computed directly from the document's formulas with Fractions
            h_ = th["T"] / (1 + th["T"])
            xv = {-1: th["U"] / (1 + th["U"]), 0: Fraction(1), 1: 1 + th["U"]}[xs]
            yv = {-1: th["V"] / (1 + th["V"]), 0: Fraction(1), 1: 1 + th["V"]}[ys]
            zv = {-1: th["A"] / (1 + th["A"]), 0: Fraction(1), 1: 1 + th["A"]}[zs]
            Mv = th["M"]
            mult = Fraction(1)
            if stg is not None:
                m0, dM, nn, dd = stg
                nnv, ddv = Rc.evaluate(nn, pt), Rc.evaluate(dd, pt)
                Mv = Fraction(nnv) / Fraction(ddv)             # = (y+z) M / ((y+z-1)(1+M))
                if Mv != (yv + zv) * th["M"] / ((yv + zv - 1) * (1 + th["M"])):
                    bad += 1
                mult = Fraction(ddv) ** dM / Mv ** m0
            c_ = Mv
            a_ = h_ / c_
            b_ = xv * c_ / h_
            r_ = yv / c_
            s_ = zv / c_
            pk = ch["pkind"]
            p_ = b_ / th["U"] * th["P"] / (1 + th["P"]) if pk == "PB" else \
                r_ / th["V"] * th["P"] / (1 + th["P"]) if pk == "PR" else th["P"]
            q_ = s_ / th["A"] * th["Q"] / (1 + th["Q"]) if zs == 1 else th["Q"]
            orig = [p_, r_, s_, q_, a_, b_, c_]
            # multiplier of the first clearing, at the (possibly substituted) value of M
            pt2 = [Mv] + pt[1:]
            for v, dg in zip(ch["vals"], degs):
                mult *= Fraction(Rc.evaluate(v[1], pt2)) ** dg
            g = M.eval_terms_common([Gt], orig)[0]
            npts += 1
            if not (mult > 0 and val == g * mult and g > 0):
                bad += 1
            # region: the bounds that this chart imposes do hold
            okr = a_ * c_ < 1 and a_ * b_ == xv and c_ * r_ == yv and c_ * s_ == zv
            if pk == "PB":
                okr &= a_ * p_ * b_ < p_ + b_
            if pk == "PR":
                okr &= c_ * p_ * r_ < p_ + r_
            if zs == 1:
                okr &= c_ * q_ * s_ < q_ + s_
            if stg is not None:
                okr &= c_ * (r_ + s_) < r_ + s_ + 1
            if not okr:
                bad += 1
    ck.check("NUMERIC", "bookkeeping of the 28 cleared polynomials, by exact evaluation at 2 random rational chart "
             "points each (seed 4601; P, Q, and U resp. A in the low domains set to 0 at random): the table "
             "polynomial equals G_9(p,r,s,q,a,b,c) at the original coefficients given by the document's formulas "
             "times the explicit positive multiplier, G_9 > 0 there, and the bounds used by the chart hold",
             bad == 0 and npts == 56, "%d points" % npts)


def hard_region_test(ck, Gt):
    """Random rational points of the hard region (10) in the ORIGINAL coefficients, independent
    of all chart polynomials: G_9 > 0."""
    rng = random.Random(4602)

    def pos(lo=1, hi=40, den=12):
        return Fraction(rng.randint(lo, hi), rng.randint(1, den))

    def unit_open(allow_zero):
        if allow_zero and rng.random() < 0.3:
            return Fraction(0)
        k_ = rng.randint(2, 30)
        return Fraction(rng.randint(1, k_ - 1), k_)

    cnt = {}
    bad = n = 0
    gmin = None
    tries = 0
    while n < 5400 and tries < 100000:
        tries += 1
        xs, ys, zs = rng.choice([-1, 0, 1]), rng.choice([-1, 0, 1]), rng.choice([-1, 0, 1])
        x = unit_open(True) if xs < 0 else Fraction(1) if xs == 0 else 1 + pos(1, 30, 10)
        y = unit_open(False) if ys < 0 else Fraction(1) if ys == 0 else 1 + pos(1, 30, 10)
        z = unit_open(True) if zs < 0 else Fraction(1) if zs == 0 else 1 + pos(1, 30, 10)
        hh = unit_open(False)
        if y + z > 1:
            cap = (y + z) / (y + z - 1)
            c = cap * (unit_open(False) if rng.random() < 0.9 else Fraction(999, 1000))
        else:
            c = pos()
        a, b, r, s = hh / c, x * c / hh, y / c, z / c
        caps = []
        if x > 1:
            caps.append(b / (x - 1))
        if y > 1:
            caps.append(r / (y - 1))
        if caps:
            p = min(caps) * (unit_open(True) if rng.random() < 0.9 else Fraction(999, 1000))
        else:
            p = pos() if rng.random() < 0.7 else Fraction(0)
        if z > 1:
            q = s / (z - 1) * (unit_open(True) if rng.random() < 0.9 else Fraction(999, 1000))
        else:
            q = pos() if rng.random() < 0.7 else Fraction(0)
        if not hard(p, r, a, b, c, q, s):
            bad += 1                                   # the sampler is built to produce hard points only
            continue
        g = M.eval_terms_common([Gt], [p, r, s, q, a, b, c])[0]
        n += 1
        cnt[(xs, ys, zs)] = cnt.get((xs, ys, zs), 0) + 1
        if g <= 0:
            bad += 1
        gmin = g if gmin is None or g < gmin else gmin
    ck.check("NUMERIC", "G_9 > 0 on the hard region (10), by exact evaluation at %d random rational hard points in the "
             "ORIGINAL coefficients (seed 4602; all 27 sign triples of (x-1, y-1, z-1) incl. the exact layers x, y, "
             "z = 1; b, s, p, q = 0 occur; p, q and c up to 0.999 of their bounds): no failure; every sampled point "
             "satisfies (10)" % n, bad == 0 and n == 5400 and len(cnt) == 27 and min(cnt.values()) >= 100,
             "points per sign triple between %d and %d; smallest value of G_9 found: %.3e"
             % (min(cnt.values()), max(cnt.values()), float(gmin)))


def numeric_part(ck, rd, Gt):
    """Floating-point attempts to falsify the intermediate claims."""
    import numpy as np
    ck.section("K. NUMERIC falsification attempts (floating point; evidence only)")
    fv = lambda idx: np.array(F_VALUES[idx], float)
    rng = random.Random(4603)
    lu = lambda: math.exp(rng.gauss(0.0, 1.0))
    opt = lambda: lu() * (rng.random() < 0.5)
    worst = {}

    def cost(p=0.0, r=0.0, a=0.0, b=0.0, c=0.0, q=0.0, s=0.0):
        return p * fv(0) + r * fv(1) + a * fv(2) + b * fv(12) + c * fv(59) + q * fv(67) + fv(70) + s * fv(72)

    def run(name, vals, seed):
        worst[name] = min(worst.get(name, 1e9), M.numeric_min_phi(vals, seed, 4000, 6))

    for trial in range(12):
        pv, rv, av, bv, qv, sv = lu(), lu(), lu(), lu(), lu(), lu()
        run("(2) SOC ac = 1", cost(a=av, c=1 / av, p=opt(), r=opt(), b=opt(), q=opt(), s=opt()), 46000 + trial)
        run("(3) PB chain a = 1/p + 1/b", cost(p=pv, b=bv, a=1 / pv + 1 / bv, r=opt(), q=opt(), s=opt()), 46100 + trial)
        run("(4) PR chain c = 1/p + 1/r", cost(p=pv, r=rv, c=1 / pv + 1 / rv, b=opt(), q=opt(), s=opt()), 46200 + trial)
        run("(5) QS chain c = 1/q + 1/s", cost(q=qv, s=sv, c=1 / qv + 1 / sv, b=opt(), p=opt(), r=opt()), 46300 + trial)
        s2 = sv * (trial % 3 != 0)
        run("(6) old-H c = 1 + 1/(r+s)", cost(r=rv, s=s2, c=1 + 1 / (rv + s2), p=opt(), b=opt(), q=opt()), 46400 + trial)
    for nm, wv in worst.items():
        ck.check("NUMERIC", "capacity %s at its threshold on this face (12 draws, coefficients exp(N(0,1)), the unused "
                 "coefficients present at random; 4000 random rotations + 6 BFGS descents each): min Phi found "
                 ">= -1e-9" % nm, wv > -1e-9, "min = %.3e" % wv)
    ck.info("controls (the search does detect negative values): min Phi found for f70 alone %.4f; for "
            "f2 + 0.5 f59 + f70 (ac = 0.5 < 1) %.4f; for f1 + 1.5 f59 + f70 (c = 1.5 < 1 + 1/r = 2) %.4f"
            % (M.numeric_min_phi(fv(70), 98, 4000, 6), M.numeric_min_phi(cost(a=1.0, c=0.5), 99, 4000, 6),
               M.numeric_min_phi(cost(r=1.0, c=1.5), 97, 4000, 6)))

    # genuine points of X: the necessary equations
    R9, G = rd["R"], rd["G"]
    eq = M.numeric_equality_points(T, 46005, 80)
    fit, cub, hes, qmin, gmax = M.numeric_radial_test(rd, ORDER, eq)
    ck.check("NUMERIC", "end-to-end test of the necessary equations on %d generic points x of X with "
             "Re(F_sigma - F_pi*) = 0 on T (seed 46005; eight labels apart by > 0.05 in normalised distance; the "
             "other gaps have either sign, so the coordinates have mixed signs): relative |C0 + H Q| < 1e-9 and "
             "relative |G_9| < 1e-9" % len(eq),
             len(eq) >= 40 and fit < 1e-9 and cub < 1e-9 and gmax < 1e-9,
             "fit %.1e, cubic %.1e, |QG| %.1e, |G| %.1e, smallest relative |Q| %.1e" % (fit, cub, hes, gmax, qmin))
    # the stratum: gaps >= 0, zero on T
    absG = {kk_: abs(vv) for kk_, vv in G.items()}
    pts = []
    for i, eps in enumerate((0.02, 0.1, 0.3, 0.005)):
        pts += M.numeric_stratum_search(T, 46100 + i, 300, eps=eps)
    dpos = cover = face = 0
    dmax = 0.0
    rel = [0.0]
    seps, ncont = [], []
    which = {"(2) SOC": 0, "(3) PB": 0, "(4) PR": 0, "(5) QS": 0, "(6) old-H": 0, "zero face a, c or r": 0}
    for pnt in pts:
        g = pnt["gaps"]
        co = dict(zip(R9.names[:8], [max(float(g[s_]), 0.0) for s_ in EXCL]))
        if co["d"] <= 1e-4:
            continue
        dpos += 1
        dmax = max(dmax, co["d"])
        seps.append(pnt["sep"])
        ncont.append(sum(1 for x_ in g if abs(x_) < 1e-7))
        co["h"] = 0.0
        rel.append(abs(R9.evaluate(G, co)) / R9.evaluate(absG, co))
        p, r, a, b, c, q, s = [co[nm] / co["d"] for nm in "prabcqs"]
        tol = 1e-7
        hit = []
        if min(a, c, r) < 1e-9:
            hit.append("zero face a, c or r")
            face += 1
        if a * c >= 1 - tol:
            hit.append("(2) SOC")
        if p > 0 and b > 0 and a >= 1 / p + 1 / b - tol:
            hit.append("(3) PB")
        if p > 0 and r > 0 and c >= 1 / p + 1 / r - tol:
            hit.append("(4) PR")
        if q > 0 and s > 0 and c >= 1 / q + 1 / s - tol:
            hit.append("(5) QS")
        if r + s > 0 and c >= 1 + 1 / (r + s) - tol:
            hit.append("(6) old-H")
        for h_ in hit:
            which[h_] += 1
        cover += bool(hit)
    ck.info("stratum (gaps >= 0, zero on T, sum 1): %d points from two-stage LM searches (seeds 46100..46103, 300 "
            "starts each, repulsion 0.02, 0.1, 0.3, 0.005; residual < 1e-10); %d have d > 1e-4 (largest d = %.2e); "
            "among these the smallest normalised distance of two labels is %.2e and the numbers of contacts "
            "(|gap| < 1e-7) are %s"
            % (len(pts), dpos, dmax, min(seps or [float("nan")]), sorted(set(ncont))))
    ck.check("NUMERIC", "at every stratum point found with d > 1e-4 (%d points): relative |G_9| < 1e-7 and a capacity "
             "(2)-(6) holds or the point lies on a zero face%s" % (dpos, "" if dpos else "  [VOID: no such point was found]"),
             cover == dpos and max(rel) < 1e-7,
             "max relative |G_9| %.1e; number of points at which each alternative holds: %s" % (max(rel), which))


if __name__ == "__main__":
    main()
