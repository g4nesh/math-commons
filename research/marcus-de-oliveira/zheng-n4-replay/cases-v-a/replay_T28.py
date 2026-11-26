"""Independent exact replay of the proof of Case T_28 = {2,3,5,6,9,12} (group V) of
D. Zheng, "The Marcus-de Oliveira conjecture for n = 4".

Document replayed (read as a PDF, nothing executed):
    MOC4_SUPP/proofs/manin_transfer_20261021/blind_moc4_T27_uniform_20261022.pdf   (5 pages)
The document calls the set "T27" (its numbering is one less than the paper's).

Run:   MOC4_SUPP=/path/to/supplement  python replay_T28.py  [--no-numeric]
MOC4_SUPP is used only (i) to compare sha256 hashes of the files that were read and (ii) to
compare JSON data (P, the cubic, the 82 facet vectors, the stored radial factor G9) with the
objects rebuilt here.  Nothing from the supplement is executed and no step of the replay uses
its data: the 2124-term polynomial G9 is recomputed here from the cubic (A.2), and all forty
coefficient tables (16 for the constant term, 24 Bernstein tables) and the 302 Sturm chains are
rebuilt from the formulas printed in the PDF.  Without MOC4_SUPP the comparisons are skipped.
Exit status 0 iff every check passes.  Run time about 4 minutes (1 minute with --no-numeric).
Packages: python-flint (exact polynomials, Sturm chains), sympy (one copied replay of a cited
inequality), numpy + scipy (NUMERIC part only).  Report: T28.md in this directory.
"""
import os
import random
import sys
import time
from fractions import Fraction

sys.dont_write_bytecode = True          # nothing is written to disk by this script
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import replay_T27_common as K   # noqa: E402

NU = 28
DOC = "proofs/manin_transfer_20261021/blind_moc4_T27_uniform_20261022.pdf"
RADIAL_JSON = "proofs/manin_transfer_20261021/blind_moc4_T27_full_radial_20261022.json"
T = (2, 3, 5, 6, 9, 12)                       # the document's "T27"

# ---- statements of the document, restated as data -------------------------------------------
DOC_SOURCES = [0, 1, 2, 5, 7, 9, 39, 40, 41, 67, 68]                 # p. 1
DOC_NONMINOR = [40, 68]
IDENTITIES = [([40, 68], [0, 7, 41, 67]), ([40, 1], [0, 7, 39]), ([40, 9], [2, 7, 41])]
SYM_DOC = ((2, 0, 1, 3), (2, 0, 1, 3), 1)     # p. 1: rows and columns permuted by (2,0,1,3), then transposition
# The document only says that the symmetry exchanges 40 and 68.  The full action below is mine
# (computed) and is verified by check_symmetry.
SYM_ACTION = {0: 41, 41: 0, 1: 39, 39: 1, 2: 9, 9: 2, 5: 5, 7: 67, 67: 7, 40: 68, 68: 40}
ORDER = [7, 39, 5, 41, 2, 40, 0, 67]          # p. 1, (1): coefficients q, p, j, c, e, d, b, a
NAMES = ["q", "p", "j", "c", "e", "d", "b", "a", "z"]
KNAMES = ["a", "b", "e", "p", "q", "c"]       # variables of k0, k1, k2 (d = 1)
CORE_REP, CORE_USE = (59, 53, 55, 67, 54), (5, 2, 0, 67, 40)         # p. 2: "T24 core"
CORE_DOC = ((2, 1, 3, 0), (2, 1, 0, 3), 0)    # row (2,1,0,3), column (3,1,0,2) in the document's convention
FOUR_REP = (0, 6, 10, 54)
FOUR_USE_1, FOUR_USE_2 = (0, 67, 2, 40), (67, 41, 39, 40)            # p. 2
OLDH_REP, OLDH_USE = (12, 2, 70), (2, 5, 40)                         # (3): e on f2 (aux), j on f5 (main)
# p. 3 and p. 4: numbers of terms of the cleared tables
STRICT_COUNTS = {("lt", "lt"): (135, 182, 182, 342), ("lt", "gt"): (135, 188, 164, 291),
                 ("gt", "lt"): (135, 164, 188, 291), ("gt", "gt"): (130, 180, 180, 268)}
BERN_COUNTS = {("both_le1", "lt"): (23160, 30838, 30838), ("both_le1", "gt"): (62476, 82332, 82332),
               ("a_ge1", "lt"): (8928, 18266, 25267), ("a_ge1", "gt"): (12719, 25703, 35113),
               ("b_ge1", "lt"): (9857, 20025, 27337), ("b_ge1", "gt"): (33660, 58371, 72870),
               ("gt", "lt"): (4548, 7772, 10983), ("gt", "alloc"): (30504, 63441, 100982)}
STURM_GROUPS = {"middle": 105, "endpoint": 197}                      # fourth row only


def normalise(lam):
    """The document's 'positive-cone reduction' on a dict {j: Fraction >= 0} over the eleven
    generators.  Returns (kind, coefficients, used_symmetry), kind 'minor' or 'normal'."""
    lam = {j: Fraction(lam.get(j, 0)) for j in DOC_SOURCES}
    s = min(lam[40], lam[68])
    lam[40] -= s
    lam[68] -= s
    for j in (0, 7, 41, 67):
        lam[j] += s
    if lam[40] == 0 and lam[68] == 0:
        return "minor", lam, False
    sym = False
    if lam[68] > 0:
        lam = {SYM_ACTION[j]: v for j, v in lam.items()}
        sym = True
    for res, right in ((1, (0, 7, 39)), (9, (2, 7, 41))):
        s = min(lam[40], lam[res])
        lam[40] -= s
        lam[res] -= s
        for j in right:
            lam[j] += s
    if lam[40] == 0:
        return "minor", lam, sym
    return "normal", lam, sym


def doc_endpoint(a, b, c, e, p):
    """The document's capacity endpoint J for a normal form with d = 1 (floats or Fractions are
    not mixed here: floats only).  Returns 0.0 when one of the two four-coefficient tests
    applies (or e >= E_r on the allocation chart) and inf when e = 0."""
    import math
    Le = max((a * b - 1.0) / ((a + 1.0) * (b + 1.0)), 0.0)
    Lp = max((a * c - 1.0) / ((a + 1.0) * (c + 1.0)), 0.0)
    if Le * e >= 1.0 or Lp * p >= 1.0:
        return 0.0
    if e <= 0.0:
        return math.inf
    if a * b < 1.0:
        z = math.sqrt(a * b)
        if a <= 1.0 and b <= 1.0:
            return 1.0 + 1.0 / e
        big = a if a >= 1.0 else b
        return 1.0 / e + big * (1.0 - z) ** 2 / (big - z) ** 2
    if a * c < 1.0:
        return 1.0 / e - Le
    # allocation: (a+1) delta^2 + ((a+1) c + p) delta - p a c = 0
    A2, B2, C2 = a + 1.0, (a + 1.0) * c + p, -p * a * c
    delta = (-B2 + math.sqrt(B2 * B2 - 4.0 * A2 * C2)) / (2.0 * A2)
    r = 1.0 - delta
    if a * b - r <= 0.0:
        return r * r / e
    Er = r * (a + 1.0) * (b + r) / (a * b - r)
    if e >= Er:
        return 0.0
    return r * r / e - r * r / Er


def main():
    numeric = "--no-numeric" not in sys.argv
    ck = K.Checker("Replay of Case T_28 = %s  (document blind_moc4_T27_uniform_20261022, its 'T27')" % (list(T),))
    K.print_versions(ck, numeric)
    K.check_document_identity(ck, DOC)
    K.check_foundations(ck)

    # (1) contact set
    K.check_contact_set(ck, NU, T, "T27")
    ck.info("scope: the argument uses only x in X and c(x) in N_T (gaps >= 0, gaps = 0 on T); further contacts, "
            "repeated eigenvalues and a_i + b_j = 0 are not excluded anywhere (see T28.md)")

    # (2) generators
    gens = K.check_generators(ck, T, DOC_SOURCES, DOC_NONMINOR)

    # (3) identities, symmetry, normal form
    ck.section("Source identities, symmetry and normal form (document p. 1)")
    K.check_linear_identities(ck, IDENTITIES, "the three identities")
    hits = K.check_symmetry(ck, T, gens, SYM_ACTION, SYM_DOC, seed=2801)
    stab = [g for g, gm in zip(K.G_ELEMENTS, K.G_MAPS) if K.g_set(gm, T) == frozenset(T)]
    ck.check("EXACT", "the stabiliser of T in G has order 2, so this is the only non-trivial symmetry of the cone",
             len(stab) == 2 and SYM_DOC in stab and hits == [SYM_DOC])
    ck.check("EXACT", "the eight generators of the normal form (1) are linearly independent",
             K.rank([K.F_VALUES[j] for j in ORDER]) == 8)
    rng = random.Random(2802)
    gmS = K.g_index_map(SYM_DOC)
    ok = all(gmS[gmS[s]] == s for s in range(24))
    n_minor = n_normal = n_sym = 0
    for _ in range(4000):
        mode = rng.random()
        lam = {j: (Fraction(rng.randint(0, 40), rng.randint(1, 6)) if rng.random() < mode else Fraction(0))
               for j in DOC_SOURCES}
        before = K.combo_values(lam)
        kind, out, sym = normalise(lam)
        after = K.combo_values(out)
        if sym:
            after = list(K.g_values(gmS, after))
        ok &= after == before and all(v >= 0 for v in out.values())
        if kind == "minor":
            ok &= out[40] == 0 and out[68] == 0
            n_minor += 1
        else:
            ok &= out[40] > 0 and out[68] == 0 and out[1] == 0 and out[9] == 0
            n_normal += 1
        n_sym += sym
    ck.check("EXACT", "reduction run on 4000 random nonnegative rational coefficient vectors (seed 2802): the cost "
             "is unchanged (up to the involutive symmetry when it is applied), coefficients stay >= 0, and the "
             "result is a minor combination (%d) or the normal form (1) with d > 0 and no f1, f9, f68 (%d); "
             "symmetry used %d times" % (n_minor, n_normal, n_sym), ok and n_minor > 0 and n_normal > 0 and n_sym > 0)

    # (4) radial necessity
    K.check_radial_necessity(ck)

    # (5) cubic, Hessian determinant, G9
    ck.section("The cubic and the Hessian determinant on y0 + z kappa (document p. 2)")
    R = K.Ring(NAMES)
    q, p, j, c, e, d, b, a, z = R.gens
    rd = K.radial_data(R, ORDER, log=None)
    Qdoc = d * d + d * e + d * p + e * p - e * j
    cC = R.to_dict(rd["C0"])
    d3 = tuple(3 if n == "d" else 0 for n in NAMES)
    ck.check("EXACT", "C(y0 + z kappa) = C + z Q with Q = d^2 + de + dp + ep - ej; C has %d monomials, all with "
             "coefficient +1, among them d^3: so C >= d^3 > 0 on {coefficients >= 0, d > 0}, hence Q != 0 and "
             "z = -C/Q at every point of X with this normal form" % len(cC),
             rd["Q"] == Qdoc and len(cC) == 25 and all(v == 1 for v in cC.values()) and cC.get(d3) == 1)
    ck.check("EXACT", "det Hess C(y0 + z kappa): %d terms, z-degrees %s, homogeneous of degree 14; equal to "
             "python-flint's integer determinant at 6 random integer points"
             % (len(rd["det"]), sorted(rd["det_by_z"])),
             R.total_degrees(rd["det"]) == {14} and max(rd["det_by_z"]) <= 4
             and K.det_spot_check(R, ORDER, rd["det"], random.Random(2803)))
    E = rd["E"]
    F, rem = divmod(E, rd["Q"] * rd["Q"] * (-4))
    G9 = K.poly_sqrt(R, F) if len(rem) == 0 else None
    okG = G9 is not None
    if okG:
        d9 = tuple(9 if n == "d" else 0 for n in NAMES)
        lead = R.to_dict(G9).get(d9, 0)
        if lead < 0:
            G9 = -G9
            lead = -lead
    ck.check("EXACT", "Q^4 det Hess C(y0 - (C/Q) kappa) = -4 Q^2 G9^2 as polynomials (%d terms), with G9 := my own "
             "exact square root of E / (-4 Q^2); G9 has %d monomials, is homogeneous of degree 9, does not "
             "involve z, has degree 2 in j, and the coefficient of d^9 is normalised to +1"
             % (len(E), len(G9) if okG else -1),
             okG and E == rd["Q"] * rd["Q"] * G9 * G9 * (-4) and len(G9) == 2124 and R.total_degrees(G9) == {9}
             and R.degrees(G9)[8] == 0 and R.degrees(G9)[2] == 2 and lead == 1)
    ck.info("so every point of X with normal form (1) has G9 = 0 (Q != 0).  Degrees of G9 in (q,p,j,c,e,d,b,a): %s"
            % (R.degrees(G9)[:8],))

    # d = 1: G9 = k0 + k1 j + k2 j^2 in the variables (a, b, e, p, q, c)
    RK = K.Ring(KNAMES)
    kmap = {NAMES.index(n): KNAMES.index(n) for n in KNAMES}
    kd = [{}, {}, {}]
    for ex, co in R.to_dict(G9).items():
        ne = [0] * 6
        for i, kx in kmap.items():
            ne[kx] = ex[i]
        ne = tuple(ne)
        kd[ex[2]][ne] = kd[ex[2]].get(ne, 0) + co
    kk = [RK.from_dict(x) for x in kd]
    Cn = K.convert(cC, RK, [kmap.get(i) for i in range(9)], fixed={5: 1})
    ka, kb, ke, kp, kq, kc = RK.gens
    ck.info("d = 1: k0, k1, k2 have %d, %d, %d monomials in (a,b,e,p,q,c)" % tuple(len(x) for x in kk))

    # (6) cited inequalities: orbits and coefficient forms
    ck.section("Cited inequalities on SO(4): orbit membership with the roles used (document p. 2)")
    rngc = random.Random(2804)
    for tgt in (CORE_USE, (5, 2, 67, 0, 40)):
        ks = K.find_group_elements(CORE_REP, tgt)
        ck.check("EXACT", "T24-face core: some group element maps (f59, f53, f55, f67, f54) to %s simultaneously"
                 % (tuple("f%d" % x for x in tgt),), len(ks) > 0,
                 "%d elements, e.g. %s" % (len(ks), K.G_ELEMENTS[ks[0]] if ks else None))
    ck.check("EXACT", "the element printed in the document (row (2,1,0,3), column (3,1,0,2); in my convention "
             "rho' = row, rho = column^{-1}) is one of those for the first orientation",
             K.G_ELEMENTS.index(CORE_DOC) in K.find_group_elements(CORE_REP, CORE_USE)
             and K.inverse(CORE_DOC[0]) == (3, 1, 0, 2))
    K.covariance_check(ck, "T24-face core transfer", CORE_DOC, list(zip(CORE_REP, CORE_USE)), rngc)
    ck.check("EXACT", "the T24 face is the cone of the paper's T_25 = {3,4,7,8,10,12}: its generators are "
             "f0, f1, f2, f6, f53, f54, f55, f59, f67, f70, which contain the five core vectors",
             K.generators(K.PAPER_T[25]) == [0, 1, 2, 6, 53, 54, 55, 59, 67, 70])
    hit = K.orbit_check(ck, "old-H (bound (3))", OLDH_REP, OLDH_USE, expected_orbit=288)
    if hit:
        K.covariance_check(ck, "old-H transfer", hit[0], list(zip(hit[1], OLDH_USE)), rngc)
    for nm, tgt in (("four-coefficient, first", FOUR_USE_1), ("four-coefficient, second", FOUR_USE_2)):
        hit = K.orbit_check(ck, nm, FOUR_REP, tgt, sym_positions=(0, 1),
                            expected_orbit=288 if tgt == FOUR_USE_1 else None)
        if hit:
            K.covariance_check(ck, nm + " transfer", hit[0], list(zip(hit[1], tgt)), rngc)
    # coefficient forms
    RC = K.Ring(["a", "b", "c", "e", "p", "dl", "r"])
    xa, xb, xc, xe, xp, xdl, xr = RC.gens
    one = RC.one
    ck.check("EXACT", "four-coefficient conditions with d = 1: ab(e-1) - (a+b+e+1) = e(ab-1) - (a+1)(b+1) and "
             "ac(p-1) - (a+c+p+1) = p(ac-1) - (a+1)(c+1): for ab > 1 (ac > 1) they read e >= E_ab (p >= P_ac)",
             xa * xb * (xe - one) - (xa + xb + xe + one) == xe * (xa * xb - one) - (xa + one) * (xb + one)
             and xa * xc * (xp - one) - (xa + xc + xp + one) == xp * (xa * xc - one) - (xa + one) * (xc + one))
    ck.check("EXACT", "old-H condition t(c - d) >= d^2 with (t, c, d) = (e, j, 1) is j >= 1 + 1/e, bound (3): "
             "e(j - 1) - 1 = e (j - (1 + 1/e))",
             K.RF(RC, xe * (xp - one) - one) == K.RF(RC, xe) * (K.RF(RC, xp) - (K.RF(RC, 1) + K.RF(RC, 1) / xe)))
    ck.check("EXACT", "allocation, second piece with coefficients (delta a, c, p, delta) on (f67, f41, f39, f40): "
             "delta a c (p - delta) - delta^2 (delta a + c + p + delta) = delta [p(ac - delta) - delta (a+1)(c+delta)], "
             "which is 0 on (10)",
             xdl * xa * xc * (xp - xdl) - xdl * xdl * (xdl * xa + xc + xp + xdl)
             == xdl * (xp * (xa * xc - xdl) - xdl * (xa + one) * (xc + xdl)))
    ck.check("EXACT", "allocation, first piece with coefficients (r a, b, e, r) on (f67, f0, f2, f40): the "
             "four-coefficient condition r a b (e - r) - r^2 (r a + b + e + r) = r [e (ab - r) - r (a+1)(b+r)], "
             "i.e. e >= E_r = r(a+1)(b+r)/(ab-r)",
             xr * xa * xb * (xe - xr) - xr * xr * (xr * xa + xb + xe + xr)
             == xr * (xe * (xa * xb - xr) - xr * (xa + one) * (xb + xr)))
    fa, fb, fe, fr = (K.RF(RC, v) for v in (xa, xb, xe, xr))
    f1 = K.RF(RC, 1)
    Er = fr * (fa + 1) * (fb + fr) / (fa * fb - fr)
    scaled = fr * (f1 / (fe / fr) - (fa * (fb / fr) - 1) / ((fa + 1) * (fb / fr + 1)))
    ck.check("EXACT", "allocation, first piece: bound (6) for the coefficients (r a, b, e, j, r), i.e. for "
             "(a, b/r, e/r, j/r, 1) after division by r, reads j >= r [ r/e - (ab/r - 1)/((a+1)(b/r + 1)) ] "
             "= r^2/e - r^2/E_r, which is (11)", scaled == fr * fr / fe - fr * fr / Er)
    K.replay_four_coefficient_note(ck)
    K.check_minor_leibniz(ck, [jx for jx in DOC_SOURCES if K.F_IS_MINOR[jx]])

    # (7) strict positivity of k0 on the closed hard chamber
    ck.section("k0 = C G6 and the strict constant term (document pp. 2-3)")
    G6, rem = divmod(kk[0], Cn)
    ck.check("EXACT", "k0 = C G6 exactly, with C the cubic coefficient at d = 1 (%d monomials, all +1, constant 1) "
             "and G6 a polynomial with %d monomials of degree <= 1 in e and in p" % (len(Cn), len(G6)),
             len(rem) == 0 and Cn * G6 == kk[0] and len(Cn) == 25 and len(G6) == 196
             and RK.degrees(G6)[2] == 1 and RK.degrees(G6)[3] == 1
             and all(v == 1 for v in RK.coeffs(Cn)) and RK.to_dict(Cn).get((0,) * 6) == 1)
    RG = K.Ring(["a", "b", "c", "q"])
    ga, gb, gc, gq = RG.gens
    gg = {}
    for ex, co in RK.to_dict(G6).items():
        key = (ex[2], ex[3])
        gg.setdefault(key, {})[(ex[0], ex[1], ex[5], ex[4])] = co
    gg = {key: RG.from_dict(v) for key, v in gg.items()}
    g00, gE, gP, gEP = gg[(0, 0)], gg[(1, 0)], gg[(0, 1)], gg[(1, 1)]
    ck.check("EXACT", "g00 has %d monomials, all coefficients positive, constant term 1: g00 >= 1 on the orthant"
             % len(g00), all(v > 0 for v in RG.coeffs(g00)) and RG.to_dict(g00).get((0, 0, 0, 0)) == 1)
    RS = K.Ring(["A", "T", "U", "Q"])
    sA, sT, sU, sQ = RS.gens
    gone = RG.one
    all_ok = True
    for ab in ("lt", "gt"):
        for ac in ("lt", "gt"):
            le = ga * gb - gone if ab == "gt" else RG.zero
            lp = ga * gc - gone if ac == "gt" else RG.zero
            de = (ga + gone) * (gb + gone)
            dp = (ga + gone) * (gc + gone)
            polys = [g00, gE * de + g00 * le, gP * dp + g00 * lp,
                     gEP * de * dp + gE * de * lp + gP * dp * le + g00 * le * lp]
            bn, bd = (sT, sA * (sT + 1)) if ab == "lt" else (sT + 1, sA)
            cn, cd = (sU, sA * (sU + 1)) if ac == "lt" else (sU + 1, sA)
            counts = []
            okc = True
            for poly in polys:
                out = K.subst_cleared(RS, RG.to_dict(poly), [sA, bn, cn, sQ], [RS.one, bd, cd, RS.one],
                                      RG.degrees(poly))
                counts.append(len(out))
                okc &= len(out) > 0 and all(v >= 0 for v in RS.coeffs(out))
            okc &= all(len(x) > 0 and all(v > 0 for v in RS.coeffs(x)) for x in (bd, cd))
            all_ok &= ck.check("EXACT", "chamber ab %s 1, ac %s 1: the four functions (8), multiplied by (a+1)(b+1), "
                               "(a+1)(c+1) as needed and cleared of the (positive) chart denominators "
                               "in the chart a = A, b = %s, c = %s, have only nonnegative integer coefficients; "
                               "numbers of terms %s (document: %s)"
                               % ("<" if ab == "lt" else ">", "<" if ac == "lt" else ">",
                                  "T/(A(1+T))" if ab == "lt" else "(1+T)/A",
                                  "U/(A(1+U))" if ac == "lt" else "(1+U)/A",
                                  tuple(counts), STRICT_COUNTS[(ab, ac)]),
                               okc and tuple(counts) == STRICT_COUNTS[(ab, ac)])
    rngs = random.Random(2806)
    oks = True
    nedge = 0
    for _ in range(3000):
        av = Fraction(rngs.randint(0, 60), 12) if rngs.random() > 0.1 else Fraction(0)
        bv = Fraction(rngs.randint(0, 60), 12) if rngs.random() > 0.1 else Fraction(0)
        cv = Fraction(rngs.randint(0, 60), 12) if rngs.random() > 0.1 else Fraction(0)
        if av > 0 and rngs.random() < 0.15:
            bv = 1 / av                                   # ab = 1
            nedge += 1
        if av > 0 and rngs.random() < 0.15:
            cv = 1 / av                                   # ac = 1
            nedge += 1
        qv = Fraction(rngs.randint(0, 60), 12) if rngs.random() > 0.15 else Fraction(0)
        pt = [av, bv, cv, qv]
        v00, vE, vP, vEP = (RG.evaluate(x, pt) for x in (g00, gE, gP, gEP))
        Le = max((av * bv - 1) / ((av + 1) * (bv + 1)), Fraction(0))
        Lp = max((av * cv - 1) / ((av + 1) * (cv + 1)), Fraction(0))
        oks = oks and v00 >= 1 and vE + Le * v00 >= 0 and vP + Lp * v00 >= 0 \
            and vEP + Lp * vE + Le * vP + Le * Lp * v00 >= 0
    ck.check("EXACT", "falsification attempt in the original variables: the four functions (8) are >= 0 (and "
             "g00 >= 1) at 3000 random rational points (a, b, c, q) of the closed orthant (seed 2806; zeros and "
             "%d points with ab = 1 or ac = 1 included)" % nedge, oks)
    RL = K.Ring(["g00", "gE", "gP", "gEP", "Le", "Lp", "e1", "p1"])
    l00, lE, lP, lEP, lLe, lLp, le1, lp1 = RL.gens
    # e = e1/(1 + Le e1), p = p1/(1 + Lp p1): G6 (1 + Le e1)(1 + Lp p1) = sum of the four functions times 1, e1, p1, e1 p1
    lhs = (l00 * (lLe * le1 + 1) * (lLp * lp1 + 1) + lE * le1 * (lLp * lp1 + 1) + lP * lp1 * (lLe * le1 + 1)
           + lEP * le1 * lp1)
    rhs = l00 + (lE + lLe * l00) * le1 + (lP + lLp * l00) * lp1 + (lEP + lLp * lE + lLe * lP + lLe * lLp * l00) * le1 * lp1
    ck.check("EXACT", "with e = e'/(1 + Le e'), p = p'/(1 + Lp p'): (1 + Le e')(1 + Lp p') G6 = g00 + (gE + Le g00) e' + "
             "(gP + Lp g00) p' + (gEP + Lp gE + Le gP + Le Lp g00) e'p'   (identity in eight indeterminates)",
             lhs == rhs)
    ck.info("PROOF-READ: on the hard chamber Le e < 1, Lp p < 1, so e' = e/(1 - Le e) >= 0, p' = p/(1 - Lp p) >= 0 and "
            "1/(1 + Le e') = 1 - Le e; the four functions are >= 0 on a dense set of (a,b,c,q) >= 0 and continuous, "
            "hence G6 >= g00 (1 - Le e)(1 - Lp p) > 0 and k0 = C G6 > 0 (C >= 1 at d = 1).")
    k1e0 = RK.split(kk[1], 2).get(0, RK.zero)
    k2e0 = RK.split(kk[2], 2).get(0, RK.zero)
    ck.check("EXACT", "e = 0: k2 vanishes identically and k1 has %d monomials, all with nonnegative coefficients "
             "(document: 375); so G9 = k0 + k1 j > 0 for e = 0, j >= 0 on the hard chamber" % len(k1e0),
             len(k2e0) == 0 and len(k1e0) == 375 and all(v > 0 for v in RK.coeffs(k1e0)))

    # (8) the 24 Bernstein tables
    ck.section("Bernstein data k0, 2k0 + J k1, k0 + J k1 + J^2 k2 on the eight charts (document pp. 3-4)")
    RB = K.Ring(["t", "J", "k0", "k1", "k2"])
    bt, bJ, b0, b1, b2 = RB.gens
    omt = RB.one - bt
    # with j = t J:
    ck.check("EXACT", "Bernstein identity: k0 + k1 (tJ) + k2 (tJ)^2 = (1-t)^2 k0 + t(1-t)(2k0 + J k1) + t^2 (k0 + J k1 + J^2 k2)",
             b0 + b1 * bt * bJ + b2 * bt * bt * bJ * bJ
             == omt * omt * b0 + bt * omt * (b0 * 2 + bJ * b1) + bt * bt * (b0 + bJ * b1 + bJ * bJ * b2))
    okS = (K.sturm_positive_roots([-6, 1, 1])[0] == 1 and K.sturm_positive_roots([6, -5, 1])[0] == 2
           and K.sturm_positive_roots([1, -2, 1])[0] == 1 and K.sturm_positive_roots([1, 0, 1])[0] == 0
           and K.sturm_positive_roots([-2, -2, -1, 1, 1])[0] == 1 and K.sturm_positive_roots([2, -3, 2])[0] == 0
           and K.sturm_positive_roots([-6, 11, -6, 1])[0] == 3 and K.sturm_positive_roots([4, -8, 5, -1])[0] == 2)
    ck.check("EXACT", "self-test of my Sturm routine on eight polynomials with known numbers of distinct positive "
             "roots (1, 2, 1 (double root), 0, 1, 0, 3, 2 (one double))", okS)
    RO = K.Ring(["A", "T", "W", "E", "P", "Q"])
    oA, oT, oW, oE, oP, oQ = RO.gens
    o1 = RO.one
    mds = [RK.degrees(x) for x in kk]
    md = tuple(max(m[i] for m in mds) for i in range(6))
    ck.info("degrees of (k0, k1, k2) in (a,b,e,p,q,c): %s; every chart clears each variable's denominator to the "
            "maximal degree %s" % (mds, md))
    kdict = [RK.to_dict(x) for x in kk]
    sturm_total = {"middle": 0, "endpoint": 0}
    sturm_polys = []

    def three_tables(key, ns, ds, num, den, desc):
        t0 = time.time()
        sk = [K.subst_cleared(RO, kdict[n], ns, ds, md) for n in range(3)]
        tables = [("base", sk[0]), ("middle", sk[0] * den * 2 + sk[1] * num),
                  ("endpoint", sk[0] * den * den + sk[1] * num * den + sk[2] * num * num)]
        counts = []
        notes = []
        ok = True
        for kind, poly in tables:
            counts.append(len(poly))
            cs = RO.coeffs(poly)
            nneg = sum(1 for v in cs if v < 0)
            if nneg == 0:
                ok &= len(poly) > 0
                continue
            # group by the five variables other than E and certify each group in E
            groups = {}
            for i in range(len(poly)):
                ex = tuple(int(x) for x in poly.exponent_tuple(i))
                groups.setdefault(ex[:3] + ex[4:], {})[ex[3]] = int(poly.coefficient(i))
            bad = 0
            for rest, co in groups.items():
                if min(co.values()) >= 0:
                    continue
                low = min(co)
                f = [co.get(i + low, 0) for i in range(max(co) - low + 1)]
                roots, _ = K.sturm_positive_roots(f) if f[0] > 0 else (None, 0)
                if f[0] > 0 and roots == 0:
                    bad += 1
                    sturm_polys.append(f)
                else:
                    ok = False
            sturm_total[kind] += bad
            notes.append("%s: %d negative coefficients in %d groups by (A,T,W,P,Q); each group is E^low times a "
                         "polynomial in E with positive constant term and no root in (0, infinity) (Sturm)"
                         % (kind, nneg, bad))
            ok &= key == ("a_ge1", "gt") and bad == STURM_GROUPS[kind]
        posden = all(len(x) > 0 and all(v > 0 for v in RO.coeffs(x)) for x in list(ds) + [den, num])
        ck.check("EXACT", "chart %s: all chart denominators, den and num have only positive coefficients (so the "
                 "clearing factor is > 0 for A, T, W, E, P, Q > 0), and the three cleared tables (base, middle, "
                 "endpoint) are nonnegative there; numbers of terms %s (document: %s)%s"
                 % (desc, tuple(counts), BERN_COUNTS[key], "; " + "; ".join(notes) if notes else
                    "; all integer coefficients nonnegative"),
                 ok and posden and tuple(counts) == BERN_COUNTS[key], "%.1f s" % (time.time() - t0))

    F = lambda num, den=None: K.RF(RO, num, den)
    fA, fT, fW, fE, fP, fQ = (F(v) for v in RO.gens)
    f1 = F(1)
    zz = fT / (fT + 1)                                   # z = T/(1+T) = sqrt(ab) on the charts with ab < 1

    def chart(abk, ack):
        """Numerators and denominators used for the substitution, and the document's formulas
        (p. 3) as rational functions for comparison.  Returns (ns, ds, num, den, description, ok)."""
        if abk == "both_le1":           # a = (z^2 + A)/(1 + A), b = z^2/a, J = 1 + 1/e
            an = oT * oT + oA * (oT + 1) ** 2
            ad = (oA + 1) * (oT + 1) ** 2
            bn, bd = oT * oT * (oA + 1), an
            a_doc = (zz * zz + fA) / (fA + 1)
            b_doc = zz * zz / a_doc
            dsc = "ab < 1, a, b <= 1"
        elif abk == "a_ge1":            # a = 1 + A, b = z^2/(1 + A), J = 1/e + a(1-z)^2/(a-z)^2
            an, ad = oA + 1, o1
            bn, bd = oT * oT, (oT + 1) ** 2 * (oA + 1)
            a_doc = fA + 1
            b_doc = zz * zz / (fA + 1)
            dsc = "ab < 1, a >= 1"
        elif abk == "b_ge1":            # a = z^2/(1 + A), b = 1 + A, J = 1/e + b(1-z)^2/(b-z)^2
            bn, bd = oA + 1, o1
            an, ad = oT * oT, (oT + 1) ** 2 * (oA + 1)
            b_doc = fA + 1
            a_doc = zz * zz / (fA + 1)
            dsc = "ab < 1, b >= 1"
        else:                           # a = A, b = (1 + T)/A
            an, ad = oA, o1
            bn, bd = oT + 1, oA
            a_doc = fA
            b_doc = (fT + 1) / fA
            dsc = "ab > 1"
        fa, fb = F(an, ad), F(bn, bd)
        ok = fa == a_doc and fb == b_doc
        if abk == "gt":                 # e = E_ab E/(1+E), J = 1/(E_ab E) = 1/e - 1/E_ab
            en = (oA + 1) * (oA + 1 + oT) * oE
            ed = oA * oT * (oE + 1)
            num, den = oA * oT, (oA + 1) * (oA + 1 + oT) * oE
            Eab = (fa + 1) * (fb + 1) / (fa * fb - 1)
            fe = F(en, ed)
            ok = ok and fa * fb == fT + 1 and fe == Eab * fE / (fE + 1) and F(num, den) == f1 / (Eab * fE) \
                and F(num, den) == f1 / fe - f1 / Eab
        else:
            en, ed = oE, o1
            fe = fE
            ok = ok and fa * fb == zz * zz
            if abk == "both_le1":
                num, den = oE + 1, oE
                ok = ok and F(num, den) == f1 + f1 / fe
            else:
                ff = (oA + 1) * (oT + 1) - oT
                num, den = ff * ff + oE * (oA + 1), oE * ff * ff
                big = fa if abk == "a_ge1" else fb
                ok = ok and F(num, den) == f1 / fe + big * (f1 - zz) * (f1 - zz) / ((big - zz) * (big - zz))
        if ack == "lt":                 # c = W/(a(1+W)), p = P
            cn, cd = oW * ad, an * (oW + 1)
            pn, pd = oP, o1
            ok = ok and F(cn, cd) == fW / (fa * (fW + 1)) and fa * F(cn, cd) == fW / (fW + 1)
            dsc += ", ac < 1"
        else:                           # c = (1+W)/a, p = P_ac P/(1+P), P_ac = (a+1)(c+1)/(ac-1)
            cn, cd = (oW + 1) * ad, an
            pn = (an + ad) * (an + (oW + 1) * ad) * oP
            pd = ad * an * oW * (oP + 1)
            fc = F(cn, cd)
            Pac = (fa + 1) * (fc + 1) / (fa * fc - 1)
            ok = ok and fc == (fW + 1) / fa and fa * fc == fW + 1 and F(pn, pd) == Pac * fP / (fP + 1)
            dsc += ", ac > 1"
        return [an, bn, en, pn, oQ, cn], [ad, bd, ed, pd, o1, cd], num, den, dsc, ok

    for abk in ("both_le1", "a_ge1", "b_ge1", "gt"):
        for ack in ("lt", "gt"):
            if abk == "gt" and ack == "gt":
                continue
            ns, ds, num, den, dsc, okf = chart(abk, ack)
            ck.check("EXACT", "chart %s: the substituted expressions for a, b, c, e, p and the endpoint J = num/den "
                     "are the document's formulas (identities of rational functions in A, T, W, E, P)" % dsc, okf)
            three_tables((abk, ack), ns, ds, num, den, dsc)
    # allocation chart: a = A, b = (1+T)/A, c = (1+W)/A, delta = U/(1+U) (U is the ring variable P here),
    # e = E_r E/(1+E), q = Q
    oU = oP
    oD = oU + 1
    gap = (oT + 1) * oD - 1                               # Delta
    hh = (oT + 1) * oD + oA                               # H
    pn = oU * (oA + 1) * ((oW + 1) * oD + oA * oU)
    pd = oA * oD * ((oW + 1) * oD - oU)
    en = (oA + 1) * hh * oE
    ed = oA * oD * gap * (oE + 1)
    num, den = oA * gap, oD * (oA + 1) * hh * oE
    fU = fP
    fa, fb, fc = fA, (fT + 1) / fA, (fW + 1) / fA
    fdl = fU / (fU + 1)
    frr = f1 - fdl
    p_doc = fdl * (fa + 1) * (fc + fdl) / (fa * fc - fdl)
    Er = frr * (fa + 1) * (fb + frr) / (fa * fb - frr)
    e_doc = Er * fE / (fE + 1)
    ck.check("EXACT", "allocation chart (ab >= 1, ac >= 1): with a = A, b = (1+T)/A, c = (1+W)/A, delta = U/(1+U), "
             "r = 1 - delta, the substituted p, e and J = num/den equal delta(a+1)(c+delta)/(ac-delta), "
             "E_r E/(1+E) and r^2/e - r^2/E_r; also (1+W)D - U = 1 + W + WU, so the denominator of p is positive",
             F(pn, pd) == p_doc and F(en, ed) == e_doc and F(num, den) == frr * frr / e_doc - frr * frr / Er
             and (oW + 1) * oD - oU == oW * oU + oW + 1)
    three_tables(("gt", "alloc"), [oA, oT + 1, en, pn, oQ, oW + 1], [o1, oA, ed, pd, o1, oA],
                 num, den, "ab >= 1, ac >= 1 (allocation)")
    ck.check("EXACT", "Sturm certificates were needed only on the chart ab < 1, a >= 1, ac > 1: %d groups in the "
             "middle table and %d in the endpoint table (document: 105 and 197, 302 in all)"
             % (sturm_total["middle"], sturm_total["endpoint"]), sturm_total == STURM_GROUPS)
    import sympy
    xs = sympy.Symbol("x")
    oks = True
    degs = set()
    for f in sturm_polys:
        pol = sympy.Poly(list(reversed(f)), xs)
        degs.add(pol.degree())
        oks = oks and pol.count_roots(0, sympy.oo) == 0 and f[0] > 0 and K.descartes_bound(f) % 2 == 0
    ck.check("EXACT", "second opinion on the %d univariate polynomials (degrees %s): sympy's exact real-root count "
             "on [0, infinity) is 0 for each, and each has an even number of coefficient sign changes"
             % (len(sturm_polys), sorted(degs)), oks and len(sturm_polys) == 302)
    ck.info("PROOF-READ: every table is (a Bernstein datum) x (product of the chart denominators to the degrees %s) "
            "x (den or den^2), and these factors are positive when A, T, W, E, P, Q > 0.  So the three data are "
            ">= 0 on the image of the open orthant, and on its closure in {e > 0} by continuity of J; see T28.md "
            "for the boundaries ab = 1, ac = 1, a = 0 and (ab, p) = (1, 0)." % (md,))

    # (9) exact sampling of the chambers in the original variables
    exact_sampling(ck, R, G9)

    # comparisons
    K.crosscheck_with_audit_library(ck, [25, 28])
    K.crosscheck_with_supplement_data(ck)
    rj = K.load_radial_json(RADIAL_JSON)
    if rj is None:
        ck.info("MOC4_SUPP not set or radial JSON missing: stored radial factors not compared")
    else:
        K.check_document_identity(ck, RADIAL_JSON)
        nine = [f for f in rj["factors"] if f[0] == 9]
        ck.check("DATA", "the supplement's stored radial record (JSON read as data) has the same variable order, "
                 "C, Q, constant -4 and factors of (degree, multiplicity) (2,2), (9,2); its degree-9 factor "
                 "(2124 terms) equals +-G9 computed here",
                 rj["variables"] == NAMES and rj["scale"] == 1 and rj["constant"] == "-4"
                 and rj["C"] == R.to_dict(rd["C0"]) and rj["Q"] in (R.to_dict(rd["Q"]), R.to_dict(-rd["Q"]))
                 and sorted((f[0], f[1]) for f in rj["factors"]) == [(2, 2), (9, 2)]
                 and len(nine) == 1 and (nine[0][2] == R.to_dict(G9) or nine[0][2] == R.to_dict(-G9)))

    if numeric:
        numeric_part(ck, R, G9)
    ck.finish()


def exact_sampling(ck, R, G9):
    """Falsification attempt in exact rational arithmetic, in the ORIGINAL variables (no chart):
    random points of every chamber, including the closed boundaries that the charts do not reach.
    At each point: k0 > 0, 2k0 + J k1 >= 0, k0 + J k1 + J^2 k2 >= 0, and G9 > 0 at five j in [0, J)."""
    ck.section("Exact rational sampling of the chambers in the original variables (falsification attempt)")
    parts = R.split(G9, 2)                               # homogeneous K_i of degree 9 - i in (q,p,c,e,d,b,a)
    Kh = [parts.get(i, R.zero) for i in range(3)]
    rng = random.Random(2805)

    def rnd(lo=0, hi=6, den=12, zero=0.15):
        if rng.random() < zero:
            return Fraction(0)
        return Fraction(rng.randint(lo * den, hi * den), den)

    def pos(hi=6, den=12):
        return Fraction(rng.randint(1, hi * den), den)

    def test(av, bv, cv, ev, pv, qv, J):
        """Bernstein data and G9 below the endpoint at the rational point (d = 1)."""
        pt = [qv, pv, Fraction(0), cv, ev, Fraction(1), bv, av]
        from math import gcd
        N = 1
        for x in pt + [J]:
            N = N * x.denominator // gcd(N, x.denominator)
        ipt = [int(x * N) for x in pt] + [0]
        k = [Fraction(R.evaluate_int(Kh[i], ipt), N ** (9 - i)) for i in range(3)]
        ok = k[0] > 0 and 2 * k[0] + J * k[1] >= 0 and k[0] + J * k[1] + J * J * k[2] >= 0
        for _ in range(5):
            jv = J * Fraction(rng.randint(0, 999), 1000)
            ok = ok and k[0] + k[1] * jv + k[2] * jv * jv > 0
        if k[0] + 4 * J * k[1] + 16 * J * J * k[2] < 0:
            power[0] += 1
        return ok

    tallies = {}
    power = [0]

    def run(name, n, gen):
        ok = True
        cnt = 0
        while cnt < n:
            out = gen()
            if out is None:
                continue
            cnt += 1
            ok = ok and test(*out)
        tallies[name] = cnt
        return ok

    def cp_for(av, edge=False):
        """(c, p) in the hard chamber for a given a > 0: ac < 1, ac = 1 or ac > 1 with p < P_ac."""
        mode = rng.randint(0, 5)
        if av == 0:
            return rnd(), rnd()
        if mode <= 1:
            cv = Fraction(rng.randint(0, 99), 100) / av
            return cv, rnd()
        if mode == 2:
            return 1 / av, rnd(hi=30)
        cv = (1 + pos()) / av
        Pac = (av + 1) * (cv + 1) / (av * cv - 1)
        return cv, Pac * Fraction(rng.randint(0, 999), 1000)

    def g_both_le1():
        av = Fraction(rng.randint(0, 12), 12)
        bv = Fraction(rng.randint(0, 12), 12)
        if av * bv >= 1:
            return None
        ev = pos()
        cv, pv = cp_for(av)
        return av, bv, cv, ev, pv, rnd(), 1 + 1 / ev

    def g_one_ge1(swap):
        def gen():
            zv = Fraction(rng.randint(0, 19), 20)
            big = 1 + rnd(zero=0.2)
            small = zv * zv / big
            av, bv = (small, big) if swap else (big, small)
            ev = pos()
            cv, pv = cp_for(av)
            return av, bv, cv, ev, pv, rnd(), 1 / ev + big * (1 - zv) ** 2 / (big - zv) ** 2
        return gen

    def g_gt_lt():
        av = pos()
        bv = (1 + rnd(zero=0.25)) / av                    # ab >= 1, with ab = 1 a quarter of the time
        cv = Fraction(rng.randint(0, 99), 100) / av       # ac < 1
        if av * bv == 1:
            ev = pos(hi=30)
            return av, bv, cv, ev, rnd(), rnd(), 1 / ev
        Eab = (av + 1) * (bv + 1) / (av * bv - 1)
        ev = Eab * Fraction(rng.randint(1, 999), 1000)
        return av, bv, cv, ev, rnd(), rnd(), 1 / ev - 1 / Eab

    def g_alloc():
        av = pos()
        bv = (1 + rnd(zero=0.25)) / av
        cv = (1 + rnd(zero=0.25)) / av
        dl = Fraction(rng.randint(0, 99), 100) if rng.random() > 0.2 else Fraction(0)
        pv = dl * (av + 1) * (cv + dl) / (av * cv - dl)
        r = 1 - dl
        if av * bv == r:                                   # ab = 1 and p = 0: endpoint 1/e
            ev = pos(hi=30)
            return av, bv, cv, ev, pv, rnd(), 1 / ev
        Er = r * (av + 1) * (bv + r) / (av * bv - r)
        ev = Er * Fraction(rng.randint(1, 999), 1000)
        return av, bv, cv, ev, pv, rnd(), r * r / ev - r * r / Er

    res = [run("ab<1, a,b<=1", 1500, g_both_le1), run("ab<1, a>=1", 1500, g_one_ge1(False)),
           run("ab<1, b>=1", 1500, g_one_ge1(True)), run("ab>=1, ac<1", 1500, g_gt_lt),
           run("ab>=1, ac>=1 (allocation)", 1500, g_alloc)]
    ck.check("EXACT", "at %d random rational points of the five chamber types (seed 2805; zeros, a = 1, b = 1, "
             "ab = 1, ac = 1, delta = 0 and the point ab = 1, p = 0 included): k0 > 0, 2k0 + J k1 >= 0, "
             "k0 + J k1 + J^2 k2 >= 0 with the document's endpoint J, and G9 > 0 at five random j in [0, J)"
             % sum(tallies.values()), all(res), str(tallies))
    ck.check("EXACT", "the test has power: with the endpoint 4J in place of J the third Bernstein datum is negative "
             "at %d of these points (the radial equation does have roots j > J)" % power[0], power[0] > 0)
    # e = 0 on the hard chamber
    ok = True
    for _ in range(1500):
        av, bv = rnd(), rnd()
        cv, pv = cp_for(av)
        pt = [rnd(), pv, rnd(hi=40), cv, Fraction(0), Fraction(1), bv, av]
        from math import gcd
        N = 1
        for x in pt:
            N = N * x.denominator // gcd(N, x.denominator)
        ok = ok and R.evaluate_int(G9, [int(x * N) for x in pt] + [0]) > 0
    ck.check("EXACT", "G9 > 0 at 1500 random rational points with e = 0, j >= 0 and (a, c, p) in the hard chamber", ok)


def numeric_part(ck, R, G9):
    """Floating-point attempts to falsify intermediate claims.  Evidence only."""
    import numpy as np
    ck.section("NUMERIC falsification attempts (floating point; evidence, not proof)")
    FV = K.F_VALUES

    def cost(av, bv, cv, ev, pv, jv, dv=1.0):
        # q f7 + p f39 + j f5 + c f41 + e f2 + d f40 + b f0 + a f67 with q = 0
        w = {39: pv, 5: jv, 41: cv, 2: ev, 40: dv, 0: bv, 67: av}
        return [sum(w[k] * FV[k][s] for k in w) for s in range(24)]

    rng = np.random.default_rng(28001)
    tests = []
    # (3) old-H: j = 1 + 1/e
    for _ in range(6):
        ev = float(np.exp(rng.uniform(-2, 2)))
        tests.append(("(3)", cost(0, 0, 0, ev, 0, 1 + 1 / ev)))
    # (4), (5) parallel
    for k in range(8):
        zv = float(rng.uniform(0.0, 0.95))
        big = float(1 + np.exp(rng.uniform(-3, 2)))
        small = zv * zv / big
        ev = float(np.exp(rng.uniform(-2, 2)))
        J = 1 / ev + big * (1 - zv) ** 2 / (big - zv) ** 2
        av, bv = (big, small) if k % 2 else (small, big)
        tests.append(("(4)/(5)", cost(av, bv, 0, ev, 0, J)))
    # (6) coherent and the double-SOC endpoint at ab = 1
    for k in range(8):
        av = float(np.exp(rng.uniform(-1.5, 1.5)))
        bv = (1.0 if k < 2 else 1 + float(np.exp(rng.uniform(-3, 2)))) / av
        if av * bv > 1:
            Eab = (av + 1) * (bv + 1) / (av * bv - 1)
            ev = Eab * float(rng.uniform(0.05, 0.98))
            J = 1 / ev - 1 / Eab
        else:
            ev = float(np.exp(rng.uniform(-2, 2)))
            J = 1 / ev
        tests.append(("(6)", cost(av, bv, 0, ev, 0, J)))
    # the two four-coefficient boundaries
    for k in range(6):
        av = float(np.exp(rng.uniform(-1, 1.5)))
        ov = (1 + float(np.exp(rng.uniform(-2, 2)))) / av
        bound = (av + 1) * (ov + 1) / (av * ov - 1)
        tests.append(("four", cost(av, ov, 0, bound, 0, 0) if k % 2 else cost(av, 0, ov, 0, bound, 0)))
    # (11) allocation
    for k in range(10):
        av = float(np.exp(rng.uniform(-1.5, 1.5)))
        bv = (1 + float(np.exp(rng.uniform(-3, 2)))) / av
        cv = (1 + float(np.exp(rng.uniform(-3, 2)))) / av
        dl = float(rng.uniform(0.02, 0.98))
        pv = dl * (av + 1) * (cv + dl) / (av * cv - dl)
        r = 1 - dl
        Er = r * (av + 1) * (bv + r) / (av * bv - r)
        ev = Er * float(rng.uniform(0.05, 0.98))
        tests.append(("(11)", cost(av, bv, cv, ev, pv, r * r / ev - r * r / Er)))
    worst = {}
    for i, (name, vals) in enumerate(tests):
        worst[name] = min(worst.get(name, float("inf")), K.numeric_min_phi(vals, 28100 + i, starts=16))
    ck.check("NUMERIC", "the cited capacities on their boundaries j = J, on the generators used here (q = 0, unused "
             "coefficients 0): min over SO(4) of the Leibniz sum >= -1e-9 at %d coefficient points (seed 28001, "
             "16 BFGS starts each)" % len(tests), min(worst.values()) >= -1e-9,
             "; ".join("%s: %.1e" % kv for kv in sorted(worst.items())))
    mn = K.numeric_min_phi(cost(2.0, 0.1, 0, 1.0, 0, 0.0), 28200, starts=16)
    ck.check("NUMERIC", "the test has power: with j = 0 at (a, b, e) = (2, 0.1, 1) the minimum is negative",
             mn < -1e-3, "%.4f" % mn)

    # roots of j -> G9 against the endpoints, on a float sample of the hard chambers
    import math
    parts = R.split(G9, 2)
    Kn = []
    for i in range(3):
        dct = R.to_dict(parts.get(i, R.zero))
        ex = np.array([[e_[7], e_[6], e_[3], e_[4], e_[1], e_[0]] for e_ in dct], float)   # a, b, c, e, p, q; d = 1
        Kn.append((ex, np.array(list(dct.values()), float)))

    def kvals(X):
        out = []
        for ex, co in Kn:
            vals = np.zeros(len(X))
            for s0 in range(0, len(X), 400):
                blk = X[s0:s0 + 400]
                vals[s0:s0 + 400] = np.prod(blk[:, None, :] ** ex[None, :, :], axis=2) @ co
            out.append(vals)
        return out

    rng = np.random.default_rng(28002)
    X = np.exp(rng.uniform(-3, 3, size=(40000, 6)))
    for col in (2, 4, 5):                                   # c, p, q vanish at about a tenth of the points
        X[rng.random(40000) < 0.1, col] = 0.0
    k0v, k1v, k2v = kvals(X)
    stats = {}
    for n in range(len(X)):
        av, bv, cv, ev, pv, qv = X[n]
        J = doc_endpoint(av, bv, cv, ev, pv)
        if J == 0.0:
            key = "not hard (a four-coefficient test applies, or e >= E_r)"
        elif av * bv < 1:
            key = "ab<1, a,b<=1" if (av <= 1 and bv <= 1) else "ab<1, max(a,b)>=1"
        else:
            key = "ab>=1, ac<1" if av * cv < 1 else "allocation"
        st = stats.setdefault(key, [0, 0, math.inf])
        st[0] += 1
        A2, B2, C2 = k2v[n], k1v[n], k0v[n]
        roots = []
        if abs(A2) > 1e-12 * (abs(B2) + abs(C2) + 1.0):
            disc = B2 * B2 - 4 * A2 * C2
            if disc >= 0:
                roots = [(-B2 - math.sqrt(disc)) / (2 * A2), (-B2 + math.sqrt(disc)) / (2 * A2)]
        elif B2 != 0:
            roots = [-C2 / B2]
        roots = [x for x in roots if x >= 0]
        if roots:
            st[1] += 1
            if J > 0:
                st[2] = min(st[2], min(roots) / J)
    hard = {k: v for k, v in stats.items() if not k.startswith("not hard")}
    ck.check("NUMERIC", "float sample of 40000 coefficient points, log-uniform in [e^-3, e^3]^6 with zeros (seed "
             "28002): wherever the quadratic j -> G9 has a root j >= 0 in the hard chamber, the smallest such root "
             "is >= J (tolerance 1e-9 J); per chamber (points, points with a root j >= 0, smallest root/J): %s"
             % "; ".join("%s: %d, %d, %s" % (k, v[0], v[1], "-" if v[2] == math.inf else "%.3f" % v[2])
                         for k, v in sorted(stats.items())),
             all(v[2] >= 1 - 1e-9 for v in hard.values()) and sum(v[1] for v in hard.values()) > 0)

    # stratum points: G9 = 0 and j >= J
    pts = []
    for k, eps in enumerate((0.02, 0.2, 1.0)):
        pts += K.numeric_stratum_search(T, 28300 + k, 800, eps=eps)
    Gd = [(tuple(ex[:8]), float(co)) for ex, co in R.to_dict(G9).items()]

    def g9(v):
        tot = 0.0
        for ex, co in Gd:
            t = co
            for i in range(8):
                if ex[i]:
                    t *= v[i] ** ex[i]
            tot += t
        return tot

    n_minor = n_normal = 0
    worst_res = worst_g = 0.0
    margin = []
    kinds = {}
    for pnt in pts:
        g = np.maximum(pnt["gaps"], 0.0)
        lam, res = K.numeric_cone_coefficients(g, DOC_SOURCES)
        worst_res = max(worst_res, res)
        kind, out, _sym = normalise({jx: max(float(v), 0.0) for jx, v in lam.items()})
        coef = np.array([float(out[jx]) for jx in ORDER])
        if kind == "minor" or coef[5] < 1e-6 * float(np.sum(coef)):
            n_minor += 1
            continue
        n_normal += 1
        v = coef / coef[5]
        qv, pv, jv, cv, ev, dv, bv, av = v
        scale = sum(abs(co) * np.prod([v[i] ** ex[i] for i in range(8) if ex[i]]) for ex, co in Gd)
        worst_g = max(worst_g, abs(g9(v)) / scale)
        J = doc_endpoint(av, bv, cv, ev, pv)
        margin.append(jv - J if np.isfinite(J) else -np.inf)
        if J == 0.0:
            if av * cv > 1 and pv * (av * cv - 1) >= (av + 1) * (cv + 1):
                key = "not hard: p >= P_ac"
            elif av * bv > 1 and ev * (av * bv - 1) >= (av + 1) * (bv + 1):
                key = "not hard: e >= E_ab"
            else:
                key = "allocation with e >= E_r"
        else:
            key = "hard, ab<1" if av * bv < 1 else ("hard, ab>=1, ac<1" if av * cv < 1 else "hard, allocation")
        kinds[key] = kinds.get(key, 0) + 1
    ck.info("stratum search: 2400 starts (seeds 28300-28302, repulsion 0.02, 0.2, 1.0), %d points with residual "
            "< 1e-10; largest distance of a gap vector from the cone of the eleven generators %.1e; %d points are "
            "minor combinations after the reduction, %d reach the normal form (1) with d > 0"
            % (len(pts), worst_res, n_minor, n_normal))
    if n_normal:
        ck.check("NUMERIC", "radial equation at stratum points in normal form (d = 1): |G9| <= 1e-6 times the sum "
                 "of the absolute values of its terms, at all %d points" % n_normal, worst_g <= 1e-6,
                 "largest relative |G9| = %.2e" % worst_g)
        ck.check("NUMERIC", "the radial root lies at or above the document's endpoint: j >= J - 1e-6 at all %d "
                 "points (J = 0 where a four-coefficient test applies); where the points lie: %s"
                 % (n_normal, kinds), min(margin) >= -1e-6, "smallest j - J = %.3e" % min(margin))
        if not any(k.startswith("hard") for k in kinds):
            ck.info("none of these points lies in the hard chamber, so on points of X the endpoints J > 0 were not "
                    "exercised by this search; they are tested on the whole coefficient region by the exact sampling "
                    "and by the float sample above")
    else:
        ck.info("no stratum point in normal form with d > 0 was found: G9 = 0 and j >= J are NOT tested by this search")


if __name__ == "__main__":
    main()
