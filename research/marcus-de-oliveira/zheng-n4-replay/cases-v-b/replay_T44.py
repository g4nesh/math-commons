"""Independent exact replay of Case T_44 = {1,4,5,6,7,8} (group V) of D. Zheng's proof of the
Marcus-de Oliveira conjecture for n = 4.

Proof document replayed:
    MOC4_SUPP/proofs/manin_transfer_20261021/independent_moc4_T43_uniform_20261022.pdf
(5 pages; the document calls the set "T43").  Every polynomial is rebuilt here from the
definitions (the matrix P, the cubic C, the vectors f_j of the paper); nothing from the
supplement is executed or imported.  If MOC4_SUPP is set, two JSON files are read as data and
COMPARED with the rebuilt objects (the stored radial polynomials, and the term counts and
sha256 digests of the 33 tables); no check depends on them.

Run:   MOC4_SUPP=/path/to/supplement  python replay_T44.py
Exit status 0 iff every check passes.  Checks tagged EXACT use integer/rational/symbolic
arithmetic only and prove what they state; checks tagged NUMERIC are falsification attempts
(floating point, or exact evaluation at random sample points) and are evidence only.
"""
import math
import os
import random
import sys
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import replay_T44_common as M                                   # noqa: E402
from replay_T44_common import Ring, PR, RF, F_VALUES             # noqa: E402

DIR = "proofs/manin_transfer_20261021/"
DOC = DIR + "independent_moc4_T43_uniform_20261022.pdf"
RADIAL_JSON = DIR + "independent_moc4_T43_full_radial_20261022.json"
TABLES_JSON = DIR + "independent_moc4_T43_uniform_exact_20261022.json"
DOC_SHA256 = {
    DOC: "6c64d2905de0720936dc4f9c9c3cb2a6e82df32e93cd51b5336580ea0e65302c",
    RADIAL_JSON: "19d676bcf16ded5704632d2f8e3bfca8dea9209bb99bd3746023a22d1177f13a",
    TABLES_JSON: "ddbb812bec5c8e2be507594e170b3e2e09d4b66d94ab693b5771aca3781cdc25",
}

T = (1, 4, 5, 6, 7, 8)
ORDER = [3, 33, 1, 59, 5, 34, 6, 0]     # g = q f3 + p f33 + j f1 + c f59 + e f5 + d f34 + b f6 + a f0
NAMES = "q,p,j,c,e,d,b,a,h"
# the document's coordinates: (q,p,j,c,e,d,b,a) = (g14, g13, g16, g0, g10, g3 - g13, g15, g19)
FORMS = [{14: 1}, {13: 1}, {16: 1}, {0: 1}, {10: 1}, {3: 1, 13: -1}, {15: 1}, {19: 1}]

# the document's term counts (constant, middle, endpoint), pp. 3-4
DOC_COUNTS = {
    ("I-", False): (1002, 1764, 2544), ("I+", False): (2982, 6381, 11472), ("II", False): (4632, 10089, 17442),
    ("IIIa-", False): (3777, 4533, 4940), ("IIIa+", False): (4284, 7308, 10564),
    ("IIIb+", False): (4572, 8784, 14012), ("IIIb-", False): (7098, 13743, 21347),
    ("IIIa-", True): (3777, 4533, 4940), ("IIIa+", True): (4284, 7308, 10564),
    ("IIIb+", True): (4428, 8448, 13472), ("IIIb-", True): (7020, 13398, 20834),
}


def coords_of(vals):
    """The eight coordinates (q,p,j,c,e,d,b,a) of a gap vector (24 values) by the document's forms."""
    return [sum(cf * vals[s] for s, cf in f.items()) for f in FORMS]


def build_charts(Rc):
    """The seven charts of Section 4 of the document as pairs (numerator, denominator) over
    Z[A,T,V,U,Q,P] for a, b, e, j and for the endpoint J = N/D0.  The representation of each
    fraction (which common factors are kept) is the one of the document's verifier, learnt by
    reading it as text, so that the term counts can be compared; any other positive clearing
    would prove the same inequalities."""
    A, T_, V, U = Rc["A"], Rc["T"], Rc["V"], Rc["U"]
    one = Rc.one
    X = {}
    X["I-"] = dict(a=(A, one), b=(T_, A * (T_ + 1)), e=(V, one), j=(V + U, one), J=(one, V + U))
    en = (A + 1) * (A + T_ + 1) * V
    ed = A * T_ * (V + 1)
    X["I+"] = dict(a=(A, one), b=(T_ + 1, A), e=(en, ed), j=(en + U * ed, ed), J=(ed, en + U * ed))
    Ln = (A + 1) * (A + T_ + 1)
    Ld = A * T_
    X["II"] = dict(a=(A, one), b=(T_ + 1, A), e=(en, ed), j=(en * U, ed * (U + 1)),
                   J=(ed * Ln - en * Ld, en * Ln))
    la = (T_ * V, (T_ + 1) * (V + 1))
    lb = (T_, T_ + 1)
    X["IIIa-"] = dict(a=la, b=lb, e=(A, one), j=(A * U, (A + 1) * (U + 1)), J=(A + 1, A))
    nj = A * (U + 1) + A * A * U
    dj = (A + 1) * (U + 1)
    X["IIIa+"] = dict(a=la, b=lb, e=(A, one), j=(nj, dj), J=(dj, nj))
    nj = A * (T_ + 1) * (U + 1) + A * A * U
    dj = (A + T_ + 1) * (U + 1)
    X["IIIb+"] = dict(a=(V, (T_ + 1) * (V + 1)), b=(T_ + 1, one), e=(A, one), j=(nj, dj), J=(dj, nj))
    ft = (T_ + 1) * (V + 1) - V
    X["IIIb-"] = dict(a=(V * V, (V + 1) ** 2 * (T_ + 1)), b=(T_ + 1, one), e=(A, one),
                      j=(A * (T_ + 1) * U, (A + T_ + 1) * (U + 1)), J=(ft * ft + A * (T_ + 1), A * ft * ft))
    return X


def main():
    import sympy as sp
    import numpy
    import scipy
    import flint
    ck = M.Checker("Case T_44 = {1,4,5,6,7,8} (document 'T43'): independent exact replay")
    ck.info("python %s, sympy %s, python-flint %s (integer polynomials of the chart tables; integer determinants "
            "at test points), numpy %s (int64 in det_linear_forms; floats in NUMERIC checks), scipy %s (NUMERIC only)"
            % (sys.version.split()[0], sp.__version__, flint.__version__, numpy.__version__, scipy.__version__))

    # ------------------------------------------------------------------------------------
    ck.section("A. Provenance and foundations")
    for rel in (DOC, RADIAL_JSON, TABLES_JSON):
        M.check_document(ck, rel, DOC_SHA256[rel])
    M.check_foundations(ck)
    M.check_radial_necessity(ck)
    M.crosscheck_supplement_tables(ck)
    M.crosscheck_audit_library(ck, [44, 46, 25])

    # ------------------------------------------------------------------------------------
    ck.section("B. Contact set (document p. 1)")
    ck.check("EXACT", "the document's set {1,4,5,6,7,8} is the paper's T_44 (Appendix E) as a set of "
             "lexicographic indices: group element = identity; the document's number is the paper's minus one",
             tuple(sorted(T)) == M.PAPER_T_V[44])
    rep = M.orbit_rep(T)
    stab = [t for t in M.G_MAPS if set(t[s] for s in T) == set(T)]
    ck.check("EXACT", "orbit representative of T_44 under G is (0,1,2,7,8,13); T_46 and T_25 lie in other orbits; "
             "the stabiliser of T_44 in G is trivial (the document: 'no nontrivial stabilizer')",
             rep == (0, 1, 2, 7, 8, 13) and M.orbit_rep(M.PAPER_T_V[46]) != rep
             and M.orbit_rep(M.PAPER_T_V[25]) != rep and len(stab) == 1,
             "T_44 -> %s, T_46 -> %s, T_25 -> %s, stabiliser order %d"
             % (rep, M.orbit_rep(M.PAPER_T_V[46]), M.orbit_rep(M.PAPER_T_V[25]), len(stab)))

    # ------------------------------------------------------------------------------------
    ck.section("C. The contact cone and the coordinates (1) (document Section 1)")
    gens = M.generators(T)
    ck.check("EXACT", "the f_j of Appendix D vanishing on T are exactly f0..f6, f33, f34, f59 (ten vectors)",
             gens == [0, 1, 2, 3, 4, 5, 6, 33, 34, 59], "computed %s" % gens)
    ck.check("EXACT", "f34 = h_{4321,+} is the only non-minor generator; the nine others are minor vectors",
             [j for j in gens if not M.F_IS_MINOR[j]] == [34],
             ", ".join("f%d = %s" % (j, M.f_name(j)) for j in gens))
    null = M.nullspace([M.P[s] for s in T], 14)
    ck.check("EXACT", "the six contact equations have rank 6 on V_R; their solution space L has dimension 8",
             M.rank([M.P[s] for s in T]) == 6 and len(null) == 8)
    ok = M.rank([M.F_COORD[j] for j in ORDER]) == 8
    for v in null:
        vals = M.values(v)
        ok &= M.combo_values(list(zip(ORDER, coords_of(vals)))) == vals
    ck.check("EXACT", "identity (1) on all of L: every c in L equals q f3 + p f33 + j f1 + c f59 + e f5 + d f34 + b f6 "
             "+ a f0 with (q,p,j,c,e,d,b,a) = (g14, g13, g16, g0, g10, g3 - g13, g15, g19), g = Pc (checked on a basis "
             "of L; the eight vectors have rank 8)", ok)
    rows = {s: [F_VALUES[j][s] for j in ORDER] for s in range(24)}
    ck.check("EXACT", "gap identities: g3 = p + d, g9 = b + d, g18 = q + d + a, g22 = j + d + a (coordinates in the "
             "order q,p,j,c,e,d,b,a)",
             rows[3] == [0, 1, 0, 0, 0, 1, 0, 0] and rows[9] == [0, 0, 0, 0, 0, 1, 1, 0]
             and rows[18] == [1, 0, 0, 0, 0, 1, 0, 1] and rows[22] == [0, 0, 1, 0, 0, 1, 0, 1])
    ck.check("EXACT", "all 24 gaps are combinations of (q,p,j,c,e,d,b,a) with coefficients in {0,1,2}: for d >= 0 the "
             "cone condition Pc >= 0 is exactly q,p,j,c,e,b,a >= 0 (seven of the coordinates are gaps themselves)",
             all(x in (0, 1, 2) for s in range(24) for x in rows[s]))
    d_, rays = M.extreme_rays_of_contact_cone(T)
    ck.check("EXACT", "direct enumeration of the extreme rays of N_T = {Pc >= 0, (Pc)_T = 0} (dimension %d; no use of "
             "the 82-facet list) gives exactly the ten vectors f0..f6, f33, f34, f59" % d_,
             rays == sorted(F_VALUES[j] for j in gens), "%d rays" % len(rays))
    add = lambda js: [sum(F_VALUES[j][s] for j in js) for s in range(24)]
    ck.check("EXACT", "f2 + f34 = f0 + f6 + f33 and f4 + f34 = f1 + f3 + f6 + f33 (all 24 values)",
             add([2, 34]) == add([0, 6, 33]) and add([4, 34]) == add([1, 3, 6, 33]))
    # the rewriting for d = -m <= 0, m = x + y: linear identity in (q,p,j,c,e,b,a,x,y)
    sq, sp_, sj, sc, se, sb, sa, sx, sy = sp.symbols("q p j c e b a x y")
    lhs = M.combo_values([(3, sq), (33, sp_), (1, sj), (59, sc), (5, se), (34, -(sx + sy)), (6, sb), (0, sa)])
    rhs = M.combo_values([(3, sq - sy), (33, sp_ - sx - sy), (1, sj - sy), (59, sc), (5, se), (6, sb - sx - sy),
                          (0, sa - sx), (2, sx), (4, sy)])
    ck.check("EXACT", "for d = -(x+y): q f3 + p f33 + j f1 + c f59 + e f5 + d f34 + b f6 + a f0 = (q-y) f3 + (p-m) f33 "
             "+ (j-y) f1 + c f59 + e f5 + (b-m) f6 + (a-x) f0 + x f2 + y f4, m = x + y (identity in nine symbols, all "
             "24 values)", all(sp.expand(l - r) == 0 for l, r in zip(lhs, rhs)))
    rng = random.Random(4401)
    n_neg = n_pos = 0
    okneg = True
    for _ in range(4000):
        lam = {j: (Fraction(rng.randint(0, 12), rng.randint(1, 4)) if rng.random() < 0.7 else Fraction(0)) for j in gens}
        vals = M.combo_values(lam)
        q, p, j, c, e, d, b, a = coords_of(vals)
        okneg &= min(q, p, j, c, e, b, a) >= 0
        if d > 0:
            n_pos += 1
            continue
        n_neg += 1
        m = -d
        x = min(a, m)
        y = m - x
        new = {3: q - y, 33: p - m, 1: j - y, 59: c, 5: e, 6: b - m, 0: a - x, 2: x, 4: y}
        okneg &= all(v >= 0 for v in new.values()) and M.combo_values(new) == vals
    ck.check("EXACT", "the document's split for d <= 0 (m = -d, x = min(a,m), y = m - x), run on 4000 random rational "
             "points of N_T (random nonnegative combinations of the ten generators, seed 4401): in all %d cases with "
             "d <= 0 the nine new coefficients are >= 0 and reproduce the vector (%d cases had d > 0)"
             % (n_neg, n_pos), okneg and n_neg > 500 and n_pos > 500)
    ck.info("PROOF-READ: g3, g9 >= 0 give p, b >= m; g18, g22 >= 0 give a + q, a + j >= m; with x = min(a, m), "
            "y = m - x: if a >= m then y = 0, else y = m - a <= q, j.  So every c in N_T with d <= 0 is a nonnegative "
            "combination of the nine minor generators.")
    M.check_minor_leibniz(ck, [j for j in gens if j != 34])

    # ------------------------------------------------------------------------------------
    ck.section("D. The sufficient coefficient conditions (document Section 2, eqs. (2), (3))")
    rng = random.Random(20261126)
    M.UNIVERSAL["Lemma 1 of the T_25 document (core f59; f53; f55, f67; f54)"] = (
        (59, 53, 55, 67, 54), (2, 3), "j f59 + e f53 + b f55 + a f67 + d f54: bounds (5), (6), double SOC")
    M.check_transport(ck, "SOC", (1, 59, 34), rng)
    M.check_transport(ck, "old-H", (5, 59, 34), rng)
    M.check_transport(ck, "four-coefficient", (0, 6, 5, 34), rng)
    M.check_transport(ck, "Lemma 1 of the T_25 document (core f59; f53; f55, f67; f54)", (59, 5, 6, 0, 34), rng)
    M.check_transport(ck, "chain", (59, 5, 6, 34), rng)
    M.check_transport(ck, "chain", (59, 5, 0, 34), rng)
    six = M.find_transport((59, 53, 55, 67, 54, 2), (59, 5, 6, 0, 34, 1))
    six2 = M.find_transport((59, 53, 55, 67, 54, 2), (59, 5, 0, 6, 34, 1))
    ck.check("EXACT", "the six-profile core of the T_25 document (f59, f53, f55, f67, f54 and the SOC endpoint f2) is "
             "mapped to (f59, f5, f6, f0, f34, f1) by a group element, and to (f59, f5, f0, f6, f34, f1) (leaves "
             "exchanged) by another one", bool(six) and bool(six2),
             "first: %s; second: %s" % (six[0][0] if six else None, six2[0][0] if six2 else None))
    # the document's explicit action row = (2,0,3,1), column = (2,3,0,1), no transpose
    row, col = (2, 0, 3, 1), (2, 3, 0, 1)
    act = [M.PIDX[tuple(col[s[row[i]]] for i in range(4))] for s in M.PERMS]      # sigma -> col o sigma o row
    old6, new6 = (59, 53, 55, 67, 54, 2), (59, 5, 6, 0, 34, 1)
    pull = all([F_VALUES[o][act[i]] for i in range(24)] == list(F_VALUES[n]) for o, n in zip(old6, new6))
    push = all(all(F_VALUES[n][act[i]] == F_VALUES[o][i] for i in range(24)) for o, n in zip(old6, new6))
    ck.check("EXACT", "the document's stated action (row (2,0,3,1), column (2,3,0,1), no transpose) does map the six "
             "profiles (59,53,55,67,54,2) to (59,5,6,0,34,1): f_new(sigma) = f_old(col o sigma o row) for all 24 sigma",
             pull, "pull-back reading: %s; push-forward reading: %s" % (pull, push))
    a, b, c, e, j, w = sp.symbols("a b c e j w", positive=True)
    four = lambda A_, B_, C_, D_: A_ * B_ * (C_ - D_) - D_ ** 2 * (A_ + B_ + C_ + D_)
    Kc = (a * b - 1) * e - (a + 1) * (b + 1)
    Emax = (a + 1) * (b + 1) / (a * b - 1)
    ck.check("EXACT", "(2): SOC ab >= d^2 for j f1 + c f59 + f34 is cj >= 1; old-H t(c-d) >= d^2 with t = e (f5), "
             "c = c (f59) is e(c-1) >= 1, i.e. c >= 1 + 1/e; the four-coefficient condition for a f0 + b f6 + e f5 + f34 "
             "is K = (ab-1)e - (a+1)(b+1) >= 0, i.e. ab > 1 and e >= Emax = (a+1)(b+1)/(ab-1)",
             sp.expand(four(a, b, e, 1) - Kc) == 0 and sp.cancel(Kc - (a * b - 1) * (e - Emax)) == 0
             and sp.cancel(e * (c - 1) - 1 - e * (c - (1 + 1 / e))) == 0)
    ck.check("EXACT", "(3), first line: 1/e - 1/Emax = -K/(e (a+1)(b+1)) (the bound (5) of the T_25 document with "
             "d = 1: d^2/e - d(ab - d^2)/((a+d)(b+d)))",
             sp.cancel(1 / e - 1 / Emax + Kc / (e * (a + 1) * (b + 1))) == 0
             and sp.cancel(1 / e - (a * b - 1) / ((a + 1) * (b + 1)) - (1 / e - 1 / Emax)) == 0)
    par_b = b * (1 - w) ** 2 / (b - w) ** 2
    ck.check("EXACT", "(3), second line: with a = w^2/b, b(1-w)^2/(b-w)^2 = (1-w)^2/(sqrt b - sqrt a)^2, symmetric in "
             "the two leaves: a(1-w)^2/(a-w)^2 has the same value; at w = 0 it is 1/b (chain endpoint 1/e + 1/b); at "
             "b = 1 it is 1 (old-H endpoint 1 + 1/e); and 1/b - b(1-w)^2/(b-w)^2 = w(b-1)(2b - w - bw)/(b(b-w)^2) >= 0 "
             "for b >= 1 > w",
             sp.cancel((w ** 2 / b) * (1 - w) ** 2 / (w ** 2 / b - w) ** 2 - par_b) == 0
             and par_b.subs(w, 0) == 1 / b and sp.cancel(par_b.subs(b, 1) - 1) == 0
             and sp.cancel(1 / b - par_b - w * (b - 1) * (2 * b - w - b * w) / (b * (b - w) ** 2)) == 0)
    ck.info("CITED statements (see the report): SOC, old-H, four-coefficient, chain as in PROOF.md 7.6; Lemma 1 of the "
            "T_25 document as verified in ../cases-v-a/T25.md, steps 14-16: for j f59 + e f53 + b f55 + a f67 + d f54 "
            "(d = 1): Phi >= 0 on SO(4) if (5) ab > 1, e > 0, j >= 1/e - (ab-1)/((a+1)(b+1)); or (6) a != b, e > 0, "
            "j >= 1/e + (1-w)^2/(sqrt b - sqrt a)^2, w = sqrt(ab) (a = 0 or b = 0 allowed); or ab >= 1 and ej >= 1")

    # ------------------------------------------------------------------------------------
    ck.section("E. Radial necessity in the coordinates (document Section 2, eqs. (4), (5))")
    rd = M.radial_data(ck, ORDER, NAMES)
    R9, C0, Q, G = rd["R"], rd["C0"], rd["Q"], rd["G"]
    ck.info("C0 = " + R9.fmt(C0, 30))
    ck.info("Q  = " + R9.fmt(Q, 30))
    v9 = {nm: R9.var(nm) for nm in R9.names}
    mul = lambda *xs: __import__("functools").reduce(M.p_mul, xs)
    brace = M.p_add(M.p_add(mul(v9["e"], v9["d"]), mul(v9["e"], v9["c"]), -1), mul(v9["d"], v9["d"]))
    Rrem = M.p_add(C0, mul(v9["p"], brace), -1)
    ck.check("EXACT", "(4): C0 = p (e(d - c) + d^2) + R, where R has %d terms, all with coefficient +1, among them d^3 "
             "(the document's C is the paper's cubic restricted to g, same sign)" % len(Rrem),
             all(vv == 1 for vv in Rrem.values()) and Rrem.get(R9.key((0, 0, 0, 0, 0, 3, 0, 0, 0))) == 1
             and len(C0) == 17)
    ck.info("PROOF-READ: if e = 0 the brace is d^2; if e > 0 and old-H fails (c < d + d^2/e) then e(d-c) + d^2 > 0.  So "
            "C0 >= d^3 > 0 at every point with d > 0 at which old-H fails, and C0 + hQ = 0 forces Q != 0, h = -C0/Q.")
    ck.check("EXACT", "G has 1232 terms (the document's count for G_9), is homogeneous of degree 9, and its d^9 "
             "coefficient is +1 (the document's sign normalisation)",
             len(G) == 1232 and all(sum(R9.exps(k)) == 9 for k in G)
             and G.get(R9.key((0, 0, 0, 0, 0, 9, 0, 0, 0))) == 1)
    M.compare_radial_json(ck, rd, RADIAL_JSON, alias={"z": "h"})
    G1 = R9.subs_one(G, "d")
    KN = ["a", "b", "e", "j", "q", "p"]
    Gt = M.packed_terms(R9, G1, KN + ["c"])
    ks = M.split_terms(Gt, 6)
    ck.check("EXACT", "at d = 1: G_9 = k0 + k1 c + k2 c^2 (degree 2 in c; no cancellation of terms at d = 1)",
             set(ks) == {0, 1, 2} and len(G1) == 1232,
             "terms of k0, k1, k2: %d, %d, %d" % (len(ks[0]), len(ks[1]), len(ks[2])))
    k = [ks[0], ks[1], ks[2]]
    swapped = {(ex[1], ex[0]) + ex[2:]: cf for ex, cf in Gt.items()}
    ck.check("EXACT", "G_9 is NOT invariant under the exchange of the two leaves a <-> b (so the document is right to "
             "certify both orientations of the ordered-leaf charts separately)", swapped != Gt,
             "%d of the 1232 monomials have a different coefficient after the exchange"
             % sum(1 for ex, cf in Gt.items() if swapped.get(ex) != cf))

    # ------------------------------------------------------------------------------------
    ck.section("F. The strict constant coefficient (document Section 3, eqs. (6), (7)) and the layer e = j = 0 (10)")
    Rk = PR(KN)
    ka, kb, ke, kj, kq, kp = Rk.gens
    kk = [Rk.from_terms(x) for x in k]
    P2 = (ka + 1) * (kb + 1) + kq * (kb + ke + 1) + ke
    pref = (kb + 1) * (kp + 1) ** 2 * P2
    H4 = kk[0] // pref
    ck.check("EXACT", "k0 = (b+1)(p+1)^2 P2 H4 with P2 = (a+1)(b+1) + q(b+e+1) + e (exact division, verified by "
             "multiplying out); H4 has 32 terms", H4 * pref == kk[0] and len(H4) == 32)
    Kk = (ka * kb - 1) * ke - (ka + 1) * (kb + 1)
    Dk = H4 + (ka + kq + 1) * Kk
    Ddoc = (ka ** 2 * kb * kj + ka ** 2 * ke * kq + ka ** 2 * kj + ka * kb * kj * kq + ka * kb * kj * 2
            + ka * ke * kj * kq + ka * ke * kj + ka * ke * kq + ka * kj * kq + ka * kj * 2 + kb * kj * kq + kb * kj
            + ke * kj * kq + ke * kj + kj * kq + kj)
    ck.check("EXACT", "(6): H4 = (a+q+1)(-K) + D, where D is the document's list of 16 monomials with positive "
             "coefficients", Dk == Ddoc and len(Dk) == 16 and M.fl_nonneg(Dk))
    ck.info("PROOF-READ: P2 >= 1 and D >= 0 on the closed orthant, so k0 >= (b+1)(p+1)^2 P2 (a+q+1)(-K) > 0 whenever "
            "K < 0, i.e. whenever the four-coefficient condition fails; this is (7), valid with any zero coefficients "
            "and at ab = 1.")
    outer = (kb + 1) ** 2 * (kp + 1) * (ka + kq + 1) ** 2
    want = [outer * (ka + 1) * (kb + 1) * (kp + 1), outer * (ka + kb + 1) * (kp + 1), outer * kp]
    lay = [Rk.from_terms({ex: cf for ex, cf in x.items() if ex[2] == 0 and ex[3] == 0}) for x in k]
    ck.check("EXACT", "(10): at e = j = 0, G_9 = (b+1)^2 (p+1) (a+q+1)^2 [(a+1)(b+1)(p+1) + (a+b+1)(p+1) c + p c^2], "
             "which is >= 1 on the closed orthant", all(x == y for x, y in zip(lay, want)))

    # ------------------------------------------------------------------------------------
    ck.section("G. The seven charts and the 33 tables (document Section 4)")
    Rc = PR("A,T,V,U,Q,P")
    X = build_charts(Rc)
    cA, cT, cV, cU = Rc["A"], Rc["T"], Rc["V"], Rc["U"]
    f = lambda pair: RF(Rc, pair[0], pair[1])
    one = RF(Rc, Rc.one)
    rA, rT, rV, rU = RF(Rc, cA), RF(Rc, cT), RF(Rc, cV), RF(Rc, cU)
    val = {nm: {x: f(ch[x]) for x in "abejJ"} for nm, ch in X.items()}
    emax = lambda v: (v["a"] + 1) * (v["b"] + 1) / (v["a"] * v["b"] - 1)
    Kof = lambda v: (v["a"] * v["b"] - 1) * v["e"] - (v["a"] + 1) * (v["b"] + 1)
    v = val["I-"]
    ck.check("EXACT", "chart I-: a = A, b = T/(A(1+T)), e = V, j = V + U: 1 - ab = 1/(1+T) (so ab < 1), j - e = U >= 0, "
             "J = 1/j", (1 - v["a"] * v["b"]) == 1 / (1 + rT) and (v["j"] - v["e"]) == rU and v["J"] == 1 / v["j"])
    v = val["I+"]
    ck.check("EXACT", "chart I+: a = A, b = (1+T)/A: ab - 1 = T; e = Emax V/(1+V), so Emax - e = Emax/(1+V) and "
             "-K = T Emax/(1+V) > 0; j - e = U; J = 1/j",
             (v["a"] * v["b"] - 1) == rT and v["e"] == emax(v) * rV / (1 + rV) and (emax(v) - v["e"]) == emax(v) / (1 + rV)
             and (0 - Kof(v)) == rT * emax(v) / (1 + rV) and (v["j"] - v["e"]) == rU and v["J"] == 1 / v["j"])
    v = val["II"]
    ck.check("EXACT", "chart II: a, b, e as in I+; j = e U/(1+U), so e - j = e/(1+U) > 0; "
             "J = 1/e - 1/Emax = -K/(e (a+1)(b+1))",
             v["a"] == val["I+"]["a"] and v["b"] == val["I+"]["b"] and v["e"] == val["I+"]["e"]
             and v["j"] == v["e"] * rU / (1 + rU) and (v["e"] - v["j"]) == v["e"] / (1 + rU)
             and v["J"] == 1 / v["e"] - 1 / emax(v) and v["J"] == (0 - Kof(v)) / (v["e"] * (v["a"] + 1) * (v["b"] + 1)))
    for nm in ("IIIa-", "IIIa+"):
        v = val[nm]
        okc = (1 - v["b"]) == 1 / (1 + rT) and (v["b"] - v["a"]) == v["b"] / (1 + rV) and v["e"] == rA
        h0 = v["e"] / (1 + v["e"])
        if nm == "IIIa-":
            okc &= v["j"] == h0 * rU / (1 + rU) and (h0 - v["j"]) == h0 / (1 + rU) and v["J"] == 1 + 1 / v["e"]
            txt = "j = (e/(1+e)) U/(1+U) < e/(1+e); J = 1 + 1/e"
        else:
            okc &= (v["j"] - h0) == v["e"] ** 2 / (1 + v["e"]) * rU / (1 + rU) \
                and (v["e"] - v["j"]) == v["e"] ** 2 / ((1 + v["e"]) * (1 + rU)) and v["J"] == 1 / v["j"]
            txt = "j - e/(1+e) = (e^2/(1+e)) U/(1+U) >= 0, e - j = e^2/((1+e)(1+U)) > 0; J = 1/j"
        ck.check("EXACT", "chart %s: b = T/(1+T) < 1, a = b V/(1+V) (b - a = b/(1+V) >= 0), e = A; %s" % (nm, txt), okc)
    for nm in ("IIIb+", "IIIb-"):
        v = val[nm]
        okc = (v["b"] - 1) == rT and v["e"] == rA
        h1 = v["e"] * v["b"] / (v["e"] + v["b"])
        if nm == "IIIb+":
            okc &= (1 - v["a"] * v["b"]) == 1 / (1 + rV) \
                and (v["j"] - h1) == v["e"] ** 2 / (v["e"] + v["b"]) * rU / (1 + rU) \
                and (v["e"] - v["j"]) == v["e"] ** 2 / ((v["e"] + v["b"]) * (1 + rU)) and v["J"] == 1 / v["j"]
            txt = ("a = V/(b(1+V)) (ab = V/(1+V) < 1); j - eb/(e+b) = (e^2/(e+b)) U/(1+U) >= 0, "
                   "e - j = e^2/((e+b)(1+U)) > 0; J = 1/j")
        else:
            wv = rV / (1 + rV)
            okc &= (v["a"] * v["b"]) == wv ** 2 and v["j"] == h1 * rU / (1 + rU) and (h1 - v["j"]) == h1 / (1 + rU) \
                and (v["b"] - wv) == (1 + rT + rT * rV) / (1 + rV) \
                and v["J"] == 1 / v["e"] + v["b"] * (1 - wv) ** 2 / (v["b"] - wv) ** 2
            txt = ("a = V^2/(b(1+V)^2), i.e. sqrt(ab) = w = V/(1+V) < 1 <= b, b - w = (1+T+TV)/(1+V) > 0; "
                   "j = (eb/(e+b)) U/(1+U) < eb/(e+b); J = 1/e + b(1-w)^2/(b-w)^2")
        ck.check("EXACT", "chart %s: b = 1 + T >= 1, e = A; %s" % (nm, txt), okc)
    # inverse maps: every point of the stated region is in the image of the chart orthant
    sA, sT, sV, sU, sQ, sP = sp.symbols("A T V U Q P")
    syms = (sA, sT, sV, sU, sQ, sP)

    def to_sym(pair):
        def poly(pp):
            return sp.Add(*[cf * sp.Mul(*[s ** ee for s, ee in zip(syms, ex) if ee]) for ex, cf in PR.terms(pp).items()])
        return poly(pair[0]) / poly(pair[1])
    a, b, e, j, w = sp.symbols("a b e j w", positive=True)
    Em = (a + 1) * (b + 1) / (a * b - 1)
    inv = {
        "I-": ({sA: a, sT: a * b / (1 - a * b), sV: e, sU: j - e}, None),
        "I+": ({sA: a, sT: a * b - 1, sV: e / (Em - e), sU: j - e}, None),
        "II": ({sA: a, sT: a * b - 1, sV: e / (Em - e), sU: j / (e - j)}, None),
        "IIIa-": ({sT: b / (1 - b), sV: a / (b - a), sA: e, sU: j * (1 + e) / (e - j * (1 + e))}, None),
        "IIIa+": ({sT: b / (1 - b), sV: a / (b - a), sA: e, sU: (j * (1 + e) - e) / ((1 + e) * (e - j))}, None),
        "IIIb+": ({sT: b - 1, sV: a * b / (1 - a * b), sA: e, sU: (j * (e + b) - e * b) / ((e + b) * (e - j))}, None),
        "IIIb-": ({sT: b - 1, sV: w / (1 - w), sA: e, sU: j * (e + b) / (e * b - j * (e + b))}, w ** 2 / b),
    }
    okinv = True
    for nm, (sub, aval) in inv.items():
        tgt = {"a": a if aval is None else aval, "b": b, "e": e, "j": j}
        for x in "abej":
            okinv &= sp.cancel(to_sym(X[nm][x]).subs(sub, simultaneous=True) - tgt[x]) == 0
    ck.check("EXACT", "inverse maps (rational identities): I-: A = a, T = ab/(1-ab), V = e, U = j - e; I+, II: A = a, "
             "T = ab - 1, V = e/(Emax - e), U = j - e resp. j/(e - j); IIIa: T = b/(1-b), V = a/(b-a), A = e, "
             "U = j(1+e)/(e - j(1+e)) resp. (j(1+e) - e)/((1+e)(e-j)); IIIb: T = b - 1, A = e, V = ab/(1-ab) and "
             "U = (j(e+b) - eb)/((e+b)(e-j)), resp. V = w/(1-w), a = w^2/b and U = j(e+b)/(eb - j(e+b)).  Each is "
             ">= 0 on the region of its chart, so the chart orthant maps ONTO that region", okinv)
    ck.info("regions (PROOF-READ from the two checks above; q = Q, p = P free):  I-: a > 0, ab < 1, j >= e >= 0;  "
            "I+: ab > 1, 0 <= e < Emax, j >= e;  II: ab > 1, 0 <= j < e < Emax (or e = j = 0);  IIIa-: 0 <= a < b < 1 "
            "(or a = b = 0), 0 <= j < e/(1+e);  IIIa+: same leaves, e/(1+e) <= j < e;  IIIb+: b >= 1, ab < 1, "
            "eb/(e+b) <= j < e;  IIIb-: b >= 1, ab < 1, 0 <= j < eb/(e+b).  'Opposite orientation' = the same with "
            "a and b exchanged.")
    common = [max(max((ex[i] for ex in x), default=0) for x in k) for i in range(6)]
    stored = None
    root = os.environ.get("MOC4_SUPP")
    if root and os.path.isfile(os.path.join(root, TABLES_JSON)):
        import json
        with open(os.path.join(root, TABLES_JSON)) as fh:
            stored = {(r["chart"], bool(r["leaves_swapped"]), r["coefficient"]): (r["terms"], r["sha256"])
                      for r in json.load(fh)["tables"]}
    ntab = 0
    nhash = 0
    for nm in ("I-", "I+", "II", "IIIa-", "IIIa+", "IIIb+", "IIIb-"):
        for swap in ((False,) if nm in ("I-", "I+", "II") else (False, True)):
            ch = X[nm]
            pa, pb = (ch["b"], ch["a"]) if swap else (ch["a"], ch["b"])
            nums = [pa[0], pb[0], ch["e"][0], ch["j"][0], Rc["Q"], Rc["P"]]
            dens = [pa[1], pb[1], ch["e"][1], ch["j"][1], Rc.one, Rc.one]
            K_ = [Rc.compose(x, nums, dens, common) for x in k]
            N_, D0 = ch["J"]
            tabs = [("constant", K_[0]), ("middle", K_[0] * D0 * 2 + K_[1] * N_),
                    ("endpoint", K_[0] * D0 * D0 + K_[1] * N_ * D0 + K_[2] * N_ * N_)]
            st = [M.fl_stats(pol) for _, pol in tabs]
            okh = True
            if stored is not None:
                for lab, pol in tabs:
                    rec = stored.get((nm, swap, lab))
                    okh &= rec is not None and rec[0] == len(pol) and rec[1] == M.fl_digest(pol)
                    nhash += 1
            ntab += 3
            ck.check("EXACT", "chart %s%s: the numerators of k0, of 2 D0 k0 + N k1 and of D0^2 k0 + N D0 k1 + N^2 k2 "
                     "(common clearing multiplier prod den_v^deg_v, degrees (3,3,2,2,2,2) in (a,b,e,j,q,p)) have only "
                     "nonnegative coefficients; term counts %s = the document's%s"
                     % (nm, " (opposite orientation: a and b exchanged)" if swap else "", DOC_COUNTS[(nm, swap)],
                        "; sha256 digests equal the supplement's records" if stored is not None else ""),
                     all(s[1] == 0 and s[0] > 0 for s in st) and tuple(s[0] for s in st) == DOC_COUNTS[(nm, swap)]
                     and okh and common == [3, 3, 2, 2, 2, 2],
                     "terms %s, smallest coefficients %s" % ([s[0] for s in st], [s[2] for s in st]))
    ck.check("EXACT", "33 tables in all (9 + 12 + 12)%s" % ("; %d digests compared with the supplement's JSON (data, "
             "comparison only)" % nhash if stored is not None else ""), ntab == 33 and nhash in (0, 33))
    if stored is None:
        ck.info("MOC4_SUPP not set: the digests of the 33 tables are not compared with the supplement's records")
    x0, x1, x2, s_, Jv = sp.symbols("x0 x1 x2 s J")
    quad = lambda cc: x0 + x1 * cc + x2 * cc ** 2
    ck.check("EXACT", "Bernstein identity: k0 + k1 c + k2 c^2 = k0 (1-s)^2 + (k0 + J k1/2) 2s(1-s) + (k0 + k1 J + k2 J^2) "
             "s^2 at c = sJ; with J = N/D0, D0 > 0: the middle and last coefficients are (2 D0 k0 + N k1)/(2 D0) and "
             "(D0^2 k0 + N D0 k1 + N^2 k2)/D0^2",
             sp.expand(quad(s_ * Jv) - (x0 * (1 - s_) ** 2 + (x0 + Jv * x1 / 2) * 2 * s_ * (1 - s_)
                                        + quad(Jv) * s_ ** 2)) == 0)
    ck.info("PROOF-READ: at a point of a chart image with D0 > 0 and N > 0 the clearing multiplier is > 0, so the two "
            "numerators are >= 0 there and G_9(c) >= k0 (1 - c/J)^2 for 0 <= c <= J: this is (9); with (7) it is > 0 "
            "for c < J when K < 0.")

    # ------------------------------------------------------------------------------------
    ck.section("H. Closed layers and completeness of the cover (document Section 5)")
    ck.info("PROOF-READ (details in the report).  Domain: d = 1, the other seven coefficients >= 0, K < 0.  "
            "(i) e = j = 0: (10).  (ii) j >= e, j > 0: J = 1/j (SOC); the two numerators 2 j k0 + k1 and "
            "j^2 k0 + j k1 + k2 are polynomials in the original coefficients, >= 0 on the images of I- and I+, hence "
            "on their closure {j >= e, K <= 0} (a = 0 and ab = 1 included).  (iii) 0 <= j < e, ab >= 1: "
            "J = -K/(e(a+1)(b+1)), again polynomial data; chart II gives ab > 1, its closure adds ab = 1, where "
            "J = 1/e and the capacity is the double-SOC bound.  (iv) 0 <= j < e, ab < 1, max(a,b) <= 1: J = 1 + 1/e "
            "or 1/j (polynomial data); the two orientations of IIIa give a < b < 1 and b < a < 1 and a = b = 0, their "
            "closure adds a = b and max(a,b) = 1.  (v) 0 <= j < e, ab < 1, max(a,b) >= 1: both orientations of IIIb "
            "are onto, zero smaller leaf included (V = 0), where the parallel endpoint is the chain endpoint.")
    Bmid_I = kk[0] * kj * 2 + kk[1]
    Bend_I = kk[0] * kj * kj + kk[1] * kj + kk[2]
    lay0 = lambda pol: Rk.from_terms({ex: cf for ex, cf in PR.terms(pol).items() if ex[0] == 0})
    ck.check("EXACT", "own direct certificate on the layer a = 0, j >= e (which the document reaches by a limit in I-): "
             "at a = 0, e = V, j = V + U the polynomials 2 j k0 + k1 and j^2 k0 + j k1 + k2 in (b,V,U,q,p) have only "
             "nonnegative coefficients",
             all(M.fl_nonneg(Rk.compose(PR.terms(lay0(pol)), [Rk.zero, kb, ke, ke + kj, kq, kp]))
                 for pol in (Bmid_I, Bend_I)))
    # ab = 1: a = 1/b; j >= e with J = 1/j and j < e with J = 1/e
    R1 = PR("B,V,U,Q,P")
    B_, V1, U1, Q1, P1 = R1.gens
    cm = [3, 3, 2, 2, 2, 2]
    o1 = R1.one
    Ks = [R1.compose(x, [o1, B_, V1, V1 + U1, Q1, P1], [B_, o1, o1, o1, o1, o1], cm) for x in k]
    okab = M.fl_nonneg(Ks[0] * (V1 + U1) * 2 + Ks[1]) and M.fl_nonneg(Ks[0] * (V1 + U1) ** 2 + Ks[1] * (V1 + U1) + Ks[2])
    Ks = [R1.compose(x, [o1, B_, V1, V1 * U1, Q1, P1], [B_, o1, o1, U1 + 1, o1, o1], cm) for x in k]
    okab &= M.fl_nonneg(Ks[0] * V1 * 2 + Ks[1]) and M.fl_nonneg(Ks[0] * V1 * V1 + Ks[1] * V1 + Ks[2])
    ck.check("EXACT", "own direct certificates on the layer ab = 1 (a = 1/b, b = B > 0; the document reaches it by "
             "limits): with e = V, j = V + U, J = 1/j, and with e = V, j = V U/(1+U), J = 1/e, the cleared middle and "
             "endpoint numerators have only nonnegative coefficients", okab)
    # exact random test of the whole cover in the original coefficients
    cover_test(ck, k)

    # ------------------------------------------------------------------------------------
    numeric_part(ck, rd, k)
    ck.finish()


def region_of(a, b, e, j, s, t):
    """The document's region and endpoint J = N/D0 at a point of the domain d = 1 (a = s^2, b = t^2 so
    that w = sqrt(ab) = st is rational).  Returns (name, N, D0)."""
    K = (a * b - 1) * e - (a + 1) * (b + 1)
    if K >= 0:
        return "four", None, None
    if e == 0 and j == 0:
        return "layer", None, None
    if j >= e:
        return "I", Fraction(1), j
    if a * b >= 1:
        return "II", -K, e * (a + 1) * (b + 1)
    m = max(a, b)
    if m <= 1:
        if j * (1 + e) <= e:
            return "IIIa-", e + 1, e
        return "IIIa+", Fraction(1), j
    if j * (e + m) >= e * m:
        return "IIIb+", Fraction(1), j
    w = s * t
    return "IIIb-", (m - w) ** 2 + e * m * (1 - w) ** 2, e * (m - w) ** 2


def cover_test(ck, k):
    rng = random.Random(4402)
    special = [Fraction(0), Fraction(0), Fraction(1), Fraction(1), Fraction(1, 2), Fraction(2), Fraction(1, 3),
               Fraction(3), Fraction(1, 10), Fraction(10)]

    def rnd():
        if rng.random() < 0.4:
            return rng.choice(special)
        return Fraction(rng.randint(0, 40), rng.randint(1, 12))

    cnt = {}
    bad = 0
    nG = 0
    minratio = None
    for it in range(6000):
        s, t = rnd(), rnd()
        r = rng.random()
        if r < 0.12 and s > 0:
            t = 1 / s                                  # ab = 1
        elif r < 0.2:
            t = s                                      # a = b
        a, b = s * s, t * t
        e, q, p = rnd(), rnd(), rnd()
        r = rng.random()
        m = max(a, b)
        if r < 0.15:
            j = e
        elif r < 0.3:
            j = e / (1 + e)
        elif r < 0.45 and e + m > 0:
            j = e * m / (e + m)
        elif r < 0.7:
            j = e * Fraction(rng.randint(0, 20), 20)
        else:
            j = rnd()
        name, N, D0 = region_of(a, b, e, j, s, t)
        cnt[name] = cnt.get(name, 0) + 1
        pt = [a, b, e, j, q, p]
        k0, k1, k2 = M.eval_terms_common(k, pt)
        if it < 25 and [k0, k1, k2] != [M.eval_terms(x, pt) for x in k]:      # the two evaluators agree
            bad += 1
        if name == "four":
            continue
        if k0 <= 0:
            bad += 1
        if name == "layer":
            for _ in range(3):
                cv = rnd()
                nG += 1
                if k0 + k1 * cv + k2 * cv * cv <= 0:
                    bad += 1
            continue
        if not (N > 0 and D0 > 0 and 2 * D0 * k0 + N * k1 >= 0 and D0 * D0 * k0 + N * D0 * k1 + N * N * k2 >= 0):
            bad += 1
        Jv = N / D0
        for _ in range(3):
            sv = Fraction(rng.randint(0, 999), 1000)
            cv = sv * Jv
            gv = k0 + k1 * cv + k2 * cv * cv
            nG += 1
            if gv <= 0 or gv < k0 * (1 - sv) ** 2:
                bad += 1
            rr = gv / (k0 * (1 - sv) ** 2)
            minratio = rr if minratio is None or rr < minratio else minratio
    ck.check("NUMERIC", "exact evaluation at 6000 random rational points of the closed domain in the ORIGINAL "
             "coefficients (seed 4402; a = s^2, b = t^2; zeros, ab = 1, a = b, max(a,b) = 1, j = e, j = e/(1+e), "
             "j = e max(a,b)/(e + max(a,b)) all occur), with the region and endpoint J = N/D0 of the report's table: "
             "N > 0, D0 > 0, k0 > 0, both Bernstein numerators >= 0, and G_9(c) >= k0 (1 - c/J)^2 > 0 at three random "
             "c in [0, J) (%d values of G_9; on the layer e = j = 0: G_9 > 0 at three random c >= 0): no failure"
             % nG, bad == 0 and all(cnt.get(nm, 0) > 20 for nm in ("I", "II", "IIIa-", "IIIa+", "IIIb+", "IIIb-", "layer")),
             "points per region %s; min of G_9/(k0 (1 - c/J)^2) = %.4f" % (dict(sorted(cnt.items())), float(minratio)))


def numeric_part(ck, rd, k):
    """Floating-point attempts to falsify the intermediate claims."""
    import numpy as np
    ck.section("I. NUMERIC falsification attempts (evidence only)")
    fv = lambda idx: np.array(F_VALUES[idx], float)
    rng = random.Random(4403)
    lu = lambda: math.exp(rng.gauss(0.0, 1.0))
    worst = {}

    def cost(q=0.0, p=0.0, j=0.0, c=0.0, e=0.0, b=0.0, a=0.0):
        return q * fv(3) + p * fv(33) + j * fv(1) + c * fv(59) + e * fv(5) + fv(34) + b * fv(6) + a * fv(0)

    def run(name, vals, seed):
        mval = M.numeric_min_phi(vals, seed, 4000, 6)
        worst[name] = min(worst.get(name, 1e9), mval)

    for trial in range(12):
        extra = dict(q=lu() * (rng.random() < 0.5), p=lu() * (rng.random() < 0.5))
        jv, ev, bv, av = lu(), lu(), lu(), lu()
        run("SOC cj = 1", cost(j=jv, c=1 / jv, e=ev * (rng.random() < 0.5), a=av * (rng.random() < 0.5), **extra),
            44000 + trial)
        run("old-H c = 1 + 1/e", cost(c=1 + 1 / ev, e=ev, j=jv * (rng.random() < 0.3), **extra), 44100 + trial)
        a2 = lu()
        b2 = (1 + lu()) / a2
        Em = (a2 + 1) * (b2 + 1) / (a2 * b2 - 1)
        run("four-coefficient K = 0", cost(a=a2, b=b2, e=Em, **extra), 44200 + trial)
        e2 = Em * rng.uniform(0.05, 1.0)
        run("coherent c = 1/e - 1/Emax (ab > 1)", cost(a=a2, b=b2, e=e2, c=1 / e2 - 1 / Em, **extra), 44300 + trial)
        wv = rng.uniform(0.0, 0.98)
        m = 1 + lu()
        sm = wv * wv / m
        par = 1 / ev + m * (1 - wv) ** 2 / (m - wv) ** 2
        if trial % 2:
            run("parallel c = 1/e + m(1-w)^2/(m-w)^2 (ab < 1)", cost(a=sm, b=m, e=ev, c=par, **extra), 44400 + trial)
        else:
            run("parallel c = 1/e + m(1-w)^2/(m-w)^2 (ab < 1)", cost(a=m, b=sm, e=ev, c=par, **extra), 44400 + trial)
        run("chain c = 1/e + 1/b (a = 0)", cost(b=bv, e=ev, c=1 / ev + 1 / bv, **extra), 44500 + trial)
        run("chain c = 1/e + 1/a (b = 0)", cost(a=av, e=ev, c=1 / ev + 1 / av, **extra), 44600 + trial)
        run("double SOC ab = 1, c = 1/e", cost(a=av, b=1 / av, e=ev, c=1 / ev, **extra), 44700 + trial)
    for nm, wv in worst.items():
        ck.check("NUMERIC", "capacity %s at its threshold on this face (12 draws, coefficients exp(N(0,1)), q and p "
                 "present at random; 4000 random rotations + 6 BFGS descents each): min Phi found >= -1e-9" % nm,
                 wv > -1e-9, "min = %.3e" % wv)
    ck.info("controls (the search does detect negative values): min Phi found for f34 alone %.4f; for "
            "f5 + 1.5 f59 + f34 (c = 1.5 < 1 + 1/e = 2) %.4f; for 2 f0 + 2 f6 + f5 + c f59 + f34 with c = 0: %.4f, with "
            "c = 0.3 < 1/e - 1/Emax = 2/3: %.1e (the coherent and parallel capacities are sufficient conditions; "
            "they are not sharp and the document does not claim that they are)"
            % (M.numeric_min_phi(fv(34), 98, 4000, 6), M.numeric_min_phi(cost(e=1.0, c=1.5), 99, 4000, 6),
               M.numeric_min_phi(cost(a=2.0, b=2.0, e=1.0, c=0.0), 97, 4000, 6),
               M.numeric_min_phi(cost(a=2.0, b=2.0, e=1.0, c=0.3), 97, 4000, 6)))

    # genuine points of X: the necessary equations
    R9, G = rd["R"], rd["G"]
    eq = M.numeric_equality_points(T, 44005, 80)
    fit, cub, hes, qmin, gmax = M.numeric_radial_test(rd, ORDER, eq)
    ck.check("NUMERIC", "end-to-end test of the necessary equations on %d generic points x of X with "
             "Re(F_sigma - F_pi*) = 0 on T (seed 44005; eight labels apart by > 0.05 in normalised distance; the "
             "other gaps have either sign, so the coordinates have mixed signs): relative |C0 + H Q| < 1e-9 and "
             "relative |G_9| < 1e-9" % len(eq),
             len(eq) >= 40 and fit < 1e-9 and cub < 1e-9 and gmax < 1e-9,
             "fit %.1e, cubic %.1e, |QG| %.1e, |G| %.1e, smallest relative |Q| %.1e" % (fit, cub, hes, gmax, qmin))
    # the stratum: gaps >= 0, zero on T
    absG = {kk_: abs(vv) for kk_, vv in G.items()}
    pts = []
    for i, eps in enumerate((0.02, 0.1, 0.3, 0.005)):
        pts += M.numeric_stratum_search(T, 44100 + i, 300, eps=eps)
    dpos = dneg = cover = negok = 0
    dmax = 0.0
    rel = [0.0]
    seps, ncont = [], []
    which = {"SOC": 0, "old-H": 0, "four-coefficient": 0, "coherent": 0, "parallel": 0}
    for pnt in pts:
        g = pnt["gaps"]
        q, p, j, c, e, d, b, a = coords_of([float(x) for x in g])
        if d <= 1e-4:
            dneg += 1
            m = max(-d, 0.0)
            x = min(a, m)
            y = m - x
            negok += min(q - y, p - m, j - y, c, e, b - m, a - x) >= -1e-7
            continue
        dpos += 1
        dmax = max(dmax, d)
        seps.append(pnt["sep"])
        ncont.append(sum(1 for x_ in g if abs(x_) < 1e-7))
        co = dict(zip(R9.names[:8], [max(x_, 0.0) for x_ in (q, p, j, c, e, d, b, a)]))
        co["h"] = 0.0
        rel.append(abs(R9.evaluate(G, co)) / R9.evaluate(absG, co))
        a, b, c, e, j = [max(x_, 0.0) / d for x_ in (a, b, c, e, j)]
        tol = 1e-7
        hit = []
        if c * j >= 1 - tol:
            hit.append("SOC")
        if e > 0 and e * (c - 1) >= 1 - tol:
            hit.append("old-H")
        if (a * b - 1) * e - (a + 1) * (b + 1) >= -tol:
            hit.append("four-coefficient")
        if a * b > 1 and e > 0 and c >= 1 / e - (a * b - 1) / ((a + 1) * (b + 1)) - tol:
            hit.append("coherent")
        if a * b < 1 and e > 0 and abs(math.sqrt(a) - math.sqrt(b)) > 1e-9 and \
                c >= 1 / e + (1 - math.sqrt(a * b)) ** 2 / (math.sqrt(b) - math.sqrt(a)) ** 2 - tol:
            hit.append("parallel")
        for h_ in hit:
            which[h_] += 1
        cover += bool(hit)
    ck.info("stratum (gaps >= 0, zero on T, sum 1): %d points from two-stage LM searches (seeds 44100..44103, 300 "
            "starts each, repulsion 0.02, 0.1, 0.3, 0.005; residual < 1e-10); %d have d > 1e-4 (largest d = %.2e), %d "
            "have d <= 1e-4; among the points with d > 1e-4 the smallest normalised distance of two labels is %.2e and "
            "the numbers of contacts (|gap| < 1e-7) are %s"
            % (len(pts), dpos, dmax, dneg, min(seps or [float("nan")]), sorted(set(ncont))))
    ck.check("NUMERIC", "at every stratum point found with d <= 1e-4 (%d points) the document's split gives "
             "nonnegative minor coefficients (tolerance 1e-7)" % dneg, negok == dneg and dneg > 0)
    ck.check("NUMERIC", "at every stratum point found with d > 1e-4 (%d points): relative |G_9| < 1e-7 and at least one "
             "capacity of Section 2 holds%s" % (dpos, "" if dpos else "  [VOID: no such point was found]"),
             cover == dpos and max(rel) < 1e-7,
             "max relative |G_9| %.1e; number of points at which each capacity holds: %s" % (max(rel), which))


if __name__ == "__main__":
    main()
