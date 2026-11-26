"""Independent exact replay of the proof of Case T_27 = {2,3,4,7,8,12} (group V) of
D. Zheng, "The Marcus-de Oliveira conjecture for n = 4".

Document replayed (read as a PDF, nothing executed):
    MOC4_SUPP/proofs/manin_transfer_20261021/final_referee_T26_uniform_20261022.pdf   (3 pages)
The document calls the set "T26" (its numbering is one less than the paper's).  It is a direct
proof; the two notes final_referee_T26_boundary10_uniform and final_referee_T26_boundary16_uniform
in the same folder treat seven-contact faces of the same cone and are NOT used by it.

Run:   MOC4_SUPP=/path/to/supplement  python replay_T27.py  [--no-numeric]
MOC4_SUPP is used only (i) to compare sha256 hashes of the files that were read and (ii) to
compare JSON data (P, the cubic, the 82 facet vectors, the stored radial factors) with the
objects rebuilt here.  Nothing from the supplement is executed and no step of the replay uses
its data; without MOC4_SUPP the comparisons are skipped.
Exit status 0 iff every check passes.
Packages: python-flint (exact polynomials), sympy (two copied replays of cited inequalities),
numpy + scipy (NUMERIC part only).  Report: T27.md in this directory.
"""
import os
import random
import sys
from fractions import Fraction

sys.dont_write_bytecode = True          # nothing is written to disk by this script
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import replay_T27_common as K   # noqa: E402

NU = 27
DOC = "proofs/manin_transfer_20261021/final_referee_T26_uniform_20261022.pdf"
RADIAL_JSON = "proofs/manin_transfer_20261021/final_referee_T26_full_radial_20261022.json"
T = (2, 3, 4, 7, 8, 12)                       # the document's "T26"

# ---- statements of the document, restated as data -------------------------------------------
DOC_SOURCES = [0, 1, 2, 5, 6, 53, 54, 55, 67, 70]                    # p. 1
DOC_NONMINOR = [54, 70]
IDENTITIES = [([54, 70], [0, 6, 55, 67]), ([54, 1], [0, 6, 53])]     # p. 1, (3)
SYM_ACTION = {0: 67, 67: 0, 1: 53, 53: 1, 6: 55, 55: 6, 54: 70, 70: 54, 2: 2, 5: 5}
SYM_DOC = ((1, 3, 0, 2), (1, 3, 0, 2), 1)     # p. 1: rows and columns permuted by (1,3,0,2), then transposition
ORDER = [0, 2, 5, 6, 53, 54, 55, 67]          # p. 1, (4): coefficients q, p, j, c, e, d, b, a
NAMES = ["q", "p", "j", "c", "e", "d", "b", "a", "z"]
SOC_REP, SOC_USE = (1, 9, 22), (2, 5, 54)             # (1): p f2 + j f5 + d f54, pj >= d^2
FOUR_REP, FOUR_USE = (0, 6, 10, 54), (67, 55, 53, 54)  # (2): a f67 + b f55 + e f53 + d f54


def normalise(lam):
    """The document's 'nonnegative source normalization' on a dict {j: Fraction >= 0} over the
    ten generators.  Returns (kind, coefficients, used_symmetry) with kind 'minor' or 'normal'."""
    lam = {j: Fraction(lam.get(j, 0)) for j in DOC_SOURCES}
    s = min(lam[54], lam[70])
    lam[54] -= s
    lam[70] -= s
    for j in (0, 6, 55, 67):
        lam[j] += s
    if lam[54] == 0 and lam[70] == 0:
        return "minor", lam, False
    sym = False
    if lam[70] > 0:
        lam = {SYM_ACTION[j]: v for j, v in lam.items()}
        sym = True
    s = min(lam[54], lam[1])
    lam[54] -= s
    lam[1] -= s
    for j in (0, 6, 53):
        lam[j] += s
    if lam[54] == 0:
        return "minor", lam, sym
    return "normal", lam, sym


def main():
    numeric = "--no-numeric" not in sys.argv
    ck = K.Checker("Replay of Case T_27 = %s  (document final_referee_T26_uniform_20261022, its 'T26')" % (list(T),))
    K.print_versions(ck, numeric)
    K.check_document_identity(ck, DOC)
    K.check_foundations(ck)

    # (1) contact set
    K.check_contact_set(ck, NU, T, "T26")
    ck.info("scope: the argument uses only x in X and c(x) in N_T (gaps >= 0, gaps = 0 on T); further contacts, "
            "repeated eigenvalues and a_i + b_j = 0 are not excluded anywhere (see T27.md)")

    # (2) generators
    gens = K.check_generators(ck, T, DOC_SOURCES, DOC_NONMINOR)

    # (3) the two identities, the symmetry and the normalisation
    ck.section("Source identities, symmetry and normal form (document p. 1)")
    K.check_linear_identities(ck, IDENTITIES, "identities (3)")
    K.check_symmetry(ck, T, gens, SYM_ACTION, SYM_DOC, seed=2701)
    ck.check("EXACT", "the eight generators of the normal form (4) are linearly independent, so the coefficients "
             "(q,p,j,c,e,d,b,a) of a vector in their span are unique", K.rank([K.F_VALUES[j] for j in ORDER]) == 8)
    rng = random.Random(2702)
    ok = True
    n_minor = n_normal = n_sym = 0
    for _ in range(4000):
        mode = rng.random()
        lam = {j: (Fraction(rng.randint(0, 40), rng.randint(1, 6)) if rng.random() < mode else Fraction(0))
               for j in DOC_SOURCES}
        before = K.combo_values(lam)
        kind, out, sym = normalise(lam)
        after = K.combo_values(out)
        if sym:   # the normal form represents g.c: undo the symmetry on the 24 values
            gm = K.g_index_map(SYM_DOC)
            after = list(K.g_values(gm, after))     # the element is an involution on these vectors
        ok &= after == before and all(v >= 0 for v in out.values())
        if kind == "minor":
            ok &= out[54] == 0 and out[70] == 0
            n_minor += 1
        else:
            ok &= out[54] > 0 and out[70] == 0 and out[1] == 0
            n_normal += 1
        n_sym += sym
    ck.check("EXACT", "normalisation run on 4000 random nonnegative rational coefficient vectors (seed 2702): the "
             "cost is unchanged (up to the symmetry when it is applied), coefficients stay >= 0, and the result "
             "is a minor combination (%d) or the normal form (4) with d > 0 and no f1, f70 (%d); symmetry used "
             "%d times" % (n_minor, n_normal, n_sym), ok and n_minor > 0 and n_normal > 0 and n_sym > 0)
    gmS = K.g_index_map(SYM_DOC)
    ck.check("EXACT", "the symmetry element acts as an involution on the 24 permutation indices",
             all(gmS[gmS[s]] == s for s in range(24)))

    # (4) radial necessity
    K.check_radial_necessity(ck)

    # (5) the cubic and the Hessian determinant on the normal form
    ck.section("The cubic and the Hessian determinant on y0 + z kappa (document p. 2)")
    R = K.Ring(NAMES)
    q, p, j, c, e, d, b, a, z = R.gens
    rd = K.radial_data(R, ORDER, log=None)
    Tt = (c + d) * (d + q)
    A0 = (c + d + q) * (a + d) * (b + d) + c * q * (a + b + d)
    Cdoc = -(A0 + e * Tt)
    Qdoc = -(d * (d + e)) - p * j
    ck.check("EXACT", "C(y0 + z kappa) = C + z Q with C = -(A0 + e T), Q = -d(d+e) - pj, "
             "T = (c+d)(d+q), A0 = (c+d+q)(a+d)(b+d) + cq(a+b+d)   (polynomial identity)",
             rd["C0"] == Cdoc and rd["Q"] == Qdoc)
    cC = R.to_dict(rd["C0"])
    d3 = tuple(3 if n == "d" else 0 for n in NAMES)
    d2 = tuple(2 if n == "d" else 0 for n in NAMES)
    ck.check("EXACT", "all %d coefficients of C are -1 and C contains -d^3; all 3 coefficients of Q are -1 and Q "
             "contains -d^2: so C <= -d^3 < 0 and Q <= -d^2 < 0 on {coefficients >= 0, d > 0}, faces included"
             % len(cC), len(cC) == 19 and all(v == -1 for v in cC.values()) and cC.get(d3) == -1
             and sorted(R.to_dict(rd["Q"]).values()) == [-1, -1, -1] and R.to_dict(rd["Q"]).get(d2) == -1)
    ck.check("EXACT", "det Hess C(y0 + z kappa): %d terms, degree <= 4 in z (here z-degrees %s), homogeneous of "
             "degree 14; equal to python-flint's integer determinant at 6 random integer points"
             % (len(rd["det"]), sorted(rd["det_by_z"])),
             R.total_degrees(rd["det"]) == {14} and K.det_spot_check(R, ORDER, rd["det"], random.Random(2703)))
    E = rd["E"]
    ck.info("E := sum_k det_k (-C)^k Q^(4-k) = Q^4 det Hess C(y0 - (C/Q) kappa): %d terms, degrees %s"
            % (len(E), sorted(R.total_degrees(E))))

    # the compact expression (6) of the document
    I = c + d + q
    Rr = p * j
    m = p + j
    U = (a + d) * (b + d)
    V = U + d * e
    Ss = a + b + d + e
    Wc = j * (a + d) + p * (b + d) + Rr * 2
    Wq = p * (a + d) + j * (b + d) + Rr * 2
    alpha = I * (U - Rr) + c * q * (a + b + d)
    Kcap = a * b * (e - d) - d * d * (a + b + e + d)
    beta = -Kcap + (m + c + q) * V + (Rr + c * q) * Ss
    Hrem = (V * (c * (c + d) * Wc + q * (q + d) * Wq)
            + c * q * (c * q * m * Ss + c * (Ss * Wc + m * V) + q * (Ss * Wq + m * V)
                       + m * (a + b + d * 2) * V + Rr * d * Ss * 4 + a * b * (d + e) * m))
    G6 = alpha * beta + c * q * (d * d - Rr) * m * Ss + Hrem
    ck.check("EXACT", "(5) with (6): Q^4 det Hess C(y0 - (C/Q) kappa) = -4 Q^2 C^2 G6^2 as polynomials, where "
             "G6 := alpha beta + cq(d^2 - R) m S + H is the document's compact expression",
             E == (rd["Q"] * rd["Q"]) * (rd["C0"] * rd["C0"]) * (G6 * G6) * (-4))
    ck.check("EXACT", "G6 has %d monomials, is homogeneous of degree 6 and does not involve z "
             "(document: degree six, 401 monomials)" % len(G6),
             len(G6) == 401 and R.total_degrees(G6) == {6} and R.degrees(G6)[8] == 0)
    F = E // (rd["Q"] * rd["Q"] * rd["C0"] * rd["C0"] * (-4))
    root = K.poly_sqrt(R, F)
    ck.check("EXACT", "independent of (6): E / (-4 Q^2 C^2) is a polynomial and my own exact square root of it "
             "equals +-G6", root is not None and (root == G6 or root == -G6))

    # (6) the sign argument
    ck.section("Sign of G6 where both coefficient inequalities fail (document p. 2)")
    hd = R.to_dict(Hrem)
    ck.check("EXACT", "H has %d monomials, all with positive integer coefficients (document: every factor in H "
             "is nonnegative)" % len(hd), len(hd) == 156 and all(v > 0 for v in hd.values()))
    ck.check("EXACT", "U - d^2 = ab + ad + bd, alpha - I(U - R) = cq(a+b+d) and beta + K = (m+c+q)V + (R+cq)S have "
             "only nonnegative coefficients; I - d = c + q; V - d^2 >= 0 coefficientwise",
             all(v > 0 for v in R.coeffs(U - d * d)) and all(v > 0 for v in R.coeffs(alpha - I * (U - Rr)))
             and all(v > 0 for v in R.coeffs(beta + Kcap)) and (I - d) == c + q
             and all(v > 0 for v in R.coeffs(V - d * d)))
    ck.info("PROOF-READ: on {q,p,j,c,e,b,a >= 0, d > 0} with R = pj < d^2 and K < 0: I >= d > 0 and "
            "U - R >= d^2 - R > 0, so alpha > 0; beta >= -K > 0; cq(d^2 - R)mS >= 0; H >= 0.  Hence G6 >= alpha beta > 0.")
    ck.info("            No division occurs, so p = 0, j = 0, a = 0, b = 0, e = 0, q = 0, c = 0 are included.")
    rng = random.Random(2704)
    tested = 0
    worst = None
    ok = True
    while tested < 20000:
        mode = rng.randint(0, 3)
        pt = []
        for n in NAMES[:8]:
            if n == "d":
                pt.append(Fraction(rng.randint(1, 10 ** 4), 10 ** rng.randint(0, 4)))
            elif rng.random() < (0.0, 0.2, 0.4, 0.7)[mode]:
                pt.append(Fraction(0))
            else:
                pt.append(Fraction(rng.randint(0, 10 ** 4), 10 ** rng.randint(0, 4)))
        qv, pv, jv, cv, ev, dv, bv, av = pt
        Kv = av * bv * (ev - dv) - dv * dv * (av + bv + ev + dv)
        if not (pv * jv < dv * dv and Kv < 0):
            continue
        tested += 1
        ipt = [int(x * 10 ** 4) for x in pt]            # G6 is homogeneous: same sign at the integer point
        val = R.evaluate_int(G6, ipt + [0])
        ok &= val > 0
        ratio = Fraction(val, ipt[5] ** 6)
        worst = ratio if worst is None or ratio < worst else worst
    ck.check("EXACT", "falsification attempt, exact arithmetic: G6 > 0 at 20000 random rational points of the "
             "region {pj < d^2, K < 0, d > 0, others >= 0} (seed 2704; up to 70 per cent of the coefficients set "
             "to 0, values spread over eight orders of magnitude)", ok, "smallest G6/d^6 = %.3e" % float(worst))
    neg = []
    rng = random.Random(2705)
    for _ in range(20000):
        pt = [Fraction(rng.randint(0, 30), rng.randint(1, 5)) for _ in range(8)]
        if pt[5] == 0:
            continue
        val = R.evaluate_int(G6, [int(x * 60) for x in pt] + [0])
        if val <= 0:
            neg.append(pt)
    okneg = all((pt[1] * pt[2] >= pt[5] ** 2) or
                (pt[7] * pt[6] * (pt[4] - pt[5]) - pt[5] ** 2 * (pt[7] + pt[6] + pt[4] + pt[5]) >= 0) for pt in neg)
    ck.check("EXACT", "the sign statement is not vacuous: among 20000 random rational points of the orthant "
             "(seed 2705) G6 <= 0 occurs %d times, every time with pj >= d^2 or K >= 0" % len(neg),
             okneg and len(neg) > 0)

    # consistency with the two seven-contact notes (not used by the document)
    ck.section("Side checks against the two seven-contact notes in the same folder (not used by the proof)")
    Gj0 = R.split(G6, R.idx("j")).get(0, R.zero)
    Xx = (a + c + d) * (b + d + q)
    ck.check("EXACT", "G6 at j = 0 equals A0 (A0 - e(ab - T)) + p (a+c+d)(b+d+q)(A0 + eT), %d monomials "
             "(boundary10 note: 206)" % len(Gj0),
             Gj0 == A0 * (A0 - e * (a * b - Tt)) + p * Xx * (A0 + e * Tt) and len(Gj0) == 206)
    Ge0 = R.split(G6, R.idx("e")).get(0, R.zero)
    base = R.from_dict({k: v for k, v in R.to_dict(Ge0).items()
                        if all(k[R.idx(n)] == 0 for n in ("q", "c", "b", "a"))})
    ck.check("EXACT", "G6 at e = 0 has %d monomials (boundary16 note: 277) and at e = q = c = b = a = 0 equals "
             "-d^2 (p+d)(j+d)(pj - d^2)" % len(Ge0),
             len(Ge0) == 277 and base == -(d * d * (p + d) * (j + d) * (p * j - d * d)))

    # (7) the two cited inequalities on SO(4)
    ck.section("Cited inequalities on SO(4): orbit membership with the roles used (document p. 1, (1) and (2))")
    rngc = random.Random(2706)
    hit = K.orbit_check(ck, "SOC", SOC_REP, SOC_USE, sym_positions=(0, 1), expected_orbit=144)
    if hit:
        K.covariance_check(ck, "SOC transfer", hit[0], list(zip(hit[1], SOC_USE)), rngc)
    hit = K.orbit_check(ck, "four-coefficient", FOUR_REP, FOUR_USE, sym_positions=(0, 1), expected_orbit=288)
    if hit:
        K.covariance_check(ck, "four-coefficient transfer", hit[0], list(zip(hit[1], FOUR_USE)), rngc)
    ck.info("coefficient forms: SOC  a f_i + b f_j + d f_k >= 0 if ab >= d^2  ->  p f2 + j f5 + d f54 if pj >= d^2;")
    ck.info("  four-coefficient  a f_a + b f_b + c f_c + d f_d >= 0 if ab(c-d) >= d^2(a+b+c+d)  ->  "
            "a f67 + b f55 + e f53 + d f54 if K >= 0.  Both statements are CITED (T27.md, citations).")
    K.replay_SOC_representative(ck)
    K.replay_four_coefficient_note(ck)
    K.check_minor_leibniz(ck, [jx for jx in DOC_SOURCES if K.F_IS_MINOR[jx]])

    # comparisons
    K.crosscheck_with_audit_library(ck, [27])
    K.crosscheck_with_supplement_data(ck)
    rj = K.load_radial_json(RADIAL_JSON)
    if rj is None:
        ck.info("MOC4_SUPP not set or radial JSON missing: stored radial factors not compared")
    else:
        K.check_document_identity(ck, RADIAL_JSON)
        six = [f for f in rj["factors"] if f[0] == 6]
        ck.check("DATA", "the supplement's stored radial record (JSON read as data) has the same variable order, "
                 "C, Q, constant -4 and factors of (degree, multiplicity) (2,2), (3,2), (6,2); its degree-6 "
                 "factor equals +-G6 rebuilt here",
                 rj["variables"] == NAMES and rj["scale"] == 1 and rj["constant"] == "-4"
                 and rj["C"] == R.to_dict(rd["C0"]) and rj["Q"] == R.to_dict(rd["Q"])
                 and sorted((f[0], f[1]) for f in rj["factors"]) == [(2, 2), (3, 2), (6, 2)]
                 and len(six) == 1 and (six[0][2] == R.to_dict(G6) or six[0][2] == R.to_dict(-G6)))

    if numeric:
        numeric_part(ck, R, G6)
    ck.finish()


def numeric_part(ck, R, G6):
    """Floating-point attempts to falsify intermediate claims.  Evidence only."""
    import numpy as np
    ck.section("NUMERIC falsification attempts (floating point; evidence, not proof)")
    # (a) the two cited inequalities on their boundaries, on the generators used here
    rng = np.random.default_rng(27001)
    worst = float("inf")
    for k in range(12):
        tau = float(np.exp(rng.uniform(-3, 3)))
        vals = [tau * x + y + z / tau for x, y, z in zip(K.F_VALUES[2], K.F_VALUES[54], K.F_VALUES[5])]
        worst = min(worst, K.numeric_min_phi(vals, 27100 + k, starts=16))
    ck.check("NUMERIC", "SOC on its boundary pj = d^2 for (f2, f5; f54): min over SO(4) of Phi[tau f2 + f54 + f5/tau] "
             ">= -1e-9 for 12 values of tau in [e^-3, e^3] (seed 27001; 16 BFGS starts each)", worst >= -1e-9,
             "smallest value found %.2e" % worst)
    vals = [0.5 * x + y + 0.5 * z for x, y, z in zip(K.F_VALUES[2], K.F_VALUES[54], K.F_VALUES[5])]
    mn = K.numeric_min_phi(vals, 27150, starts=16)
    ck.check("NUMERIC", "the test has power: with pj = d^2/4 the minimum is negative", mn < -1e-3, "%.4f" % mn)
    worst = float("inf")
    for k in range(12):
        av, bv = np.exp(rng.uniform(-1, 3, size=2))
        if av * bv <= 1.05:
            av = 1.1 / bv + av
        ev = (av * bv + av + bv + 1.0) / (av * bv - 1.0)        # K = 0 with d = 1
        vals = [av * x + bv * y + ev * z + w for x, y, z, w in
                zip(K.F_VALUES[67], K.F_VALUES[55], K.F_VALUES[53], K.F_VALUES[54])]
        worst = min(worst, K.numeric_min_phi(vals, 27200 + k, starts=16))
    ck.check("NUMERIC", "four-coefficient on its boundary K = 0 for (f67, f55; f53; f54): min over SO(4) of "
             "Phi[a f67 + b f55 + e f53 + f54] >= -1e-9 at 12 random (a, b) with ab > 1 (seed 27001)",
             worst >= -1e-9, "smallest value found %.2e" % worst)
    vals = [2.0 * x + 2.0 * y + 1.0 * z + w for x, y, z, w in
            zip(K.F_VALUES[67], K.F_VALUES[55], K.F_VALUES[53], K.F_VALUES[54])]
    mn = K.numeric_min_phi(vals, 27250, starts=16)
    ck.check("NUMERIC", "the test has power: at (a, b, e, d) = (2, 2, 1, 1), where K < 0, the minimum is negative",
             mn < -1e-3, "%.4f" % mn)

    # (b) points of the stratum: G6 = 0 and Lemma 1
    pts = []
    for k, eps in enumerate((0.02, 0.2, 1.0)):
        pts += K.numeric_stratum_search(T, 27300 + k, 500, eps=eps)
    gm = K.g_index_map(SYM_DOC)
    Gd = R.to_dict(G6)

    def g6(v):
        """(value of G6, sum of the absolute values of its terms) at the float point v."""
        tot = 0.0
        scale = 0.0
        for ex, co in Gd.items():
            t = float(co)
            for i in range(8):
                if ex[i]:
                    t *= v[i] ** ex[i]
            tot += t
            scale += abs(t)
        return tot, scale

    n_minor = n_normal = 0
    worst_g = 0.0
    worst_res = 0.0
    lemma_ok = True
    n_soc = n_four = 0
    margins = []
    for pnt in pts:
        g = np.maximum(pnt["gaps"], 0.0)
        lam, res = K.numeric_cone_coefficients(g, DOC_SOURCES)
        worst_res = max(worst_res, res)
        kind, out, _sym = normalise({jx: max(float(v), 0.0) for jx, v in lam.items()})
        coef = np.array([float(out[jx]) for jx in ORDER])
        if kind == "minor" or coef[5] < 1e-6 * float(np.sum(coef)):
            n_minor += 1                    # no non-minor generator survives (up to rounding)
            continue
        n_normal += 1
        v = coef / coef[5]
        qv, pv, jv, cv, ev, dv, bv, av = v
        val, scale = g6(v)
        worst_g = max(worst_g, abs(val) / scale)
        soc = pv * jv - 1.0
        four = av * bv * (ev - 1.0) - (av + bv + ev + 1.0)
        margins.append(max(soc, four))
        n_soc += soc >= -1e-7
        n_four += four >= -1e-7
        lemma_ok &= (soc >= -1e-6) or (four >= -1e-6)
    ck.info("stratum search: 1500 starts (seeds 27300-27302, repulsion 0.02, 0.2, 1.0), %d points with residual < 1e-10; largest distance "
            "of a gap vector from the cone of the ten generators %.1e; %d points are minor combinations or are "
            "absorbed by f1, %d reach the normal form (4) with d > 0"
            % (len(pts), worst_res, n_minor, n_normal))
    if n_normal:
        ck.check("NUMERIC", "radial equation at stratum points in normal form (d = 1): |G6| <= 1e-6 times the sum "
                 "of the absolute values of its terms, at all %d points" % n_normal, worst_g <= 1e-6,
                 "largest relative |G6| = %.2e" % worst_g)
        ck.check("NUMERIC", "Lemma 1 at these points: pj >= d^2 or K >= 0 (tolerance 1e-6); pj >= d^2 at %d, "
                 "K >= 0 at %d of %d points" % (n_soc, n_four, n_normal), lemma_ok,
                 "smallest max(pj - d^2, K) = %.3e" % min(margins))
    else:
        ck.info("no stratum point in normal form with d > 0 was found: the radial equation and Lemma 1 are NOT "
                "tested by this search")


if __name__ == "__main__":
    main()
