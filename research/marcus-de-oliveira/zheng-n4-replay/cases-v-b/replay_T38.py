"""Independent exact replay of Case T_38 = {1,2,3,8,11,12} (group V) of D. Zheng's proof of the
Marcus-de Oliveira conjecture for n = 4.

Proof document replayed: MOC4_SUPP/proofs/manin_transfer_20261021/root_moc4_T37_uniform_20261022.pdf
(4 pages; the document calls the set "T37").  Every polynomial is rebuilt here from the
definitions (the matrix P, the cubic C, the vectors f_j of the paper); nothing from the
supplement is executed or imported.  If MOC4_SUPP is set, the hash of the PDF is compared and a
few JSON files are read as data and COMPARED with the rebuilt objects; no step depends on them.

Run:   MOC4_SUPP=/path/to/supplement  python replay_T38.py  [--no-numeric]
Exit status 0 iff every check passes.  Checks tagged EXACT use integer/rational/symbolic
arithmetic only; checks tagged NUMERIC are falsification attempts on random samples (floating
point, or exact evaluation at random rational points where the text says so).
"""
import math
import os
import random
import sys
from fractions import Fraction

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import replay_T38_common as M                                   # noqa: E402
from replay_T38_common import (Ring, FL, RF, p_add, p_mul, p_scale, p_nonneg, p_compose,   # noqa: E402
                               p_compose_fl, p_rename, p_divexact, F_VALUES, SH, MASK)

DIR = "proofs/manin_transfer_20261021/"
DOC = DIR + "root_moc4_T37_uniform_20261022.pdf"
JSON_RADIAL = DIR + "root_moc4_T37_full_radial_20261022.json"
DOC_SHA256 = {
    DOC: "871f5ba783021ecd0d222845e1ead374df89b75aea7bc90131d5661933e6df1d",
}

T = (1, 2, 3, 8, 11, 12)
ORDER = [1, 2, 15, 9, 10, 56, 0, 53]          # (1): q f1 + p f2 + j f15 + c f9 + e f10 + d f56 + b f0 + a f53
NAMES = "q,p,j,c,e,d,b,a,h"
EXCL = (17, 15, 4, 10, 9, 6, 20, 14)          # permutations at which the eight generators are the unit vectors
# zero faces (p. 1): coefficient, extra contact, the paper's number of the cited case, the document's
# group element (rho, rho', transpose) and the image of the paper's contact set printed in the document
ZERO_FACES = [
    ("q", 17, 46, ((1, 3, 0, 2), (1, 2, 0, 3), 1), (1, 2, 8, 11, 12, 17), "T45"),
    ("j", 4, 4, ((0, 1, 2, 3), (0, 3, 1, 2), 0), (1, 2, 3, 4, 8, 11), "T4"),
    ("e", 9, 8, ((0, 1, 3, 2), (0, 3, 2, 1), 0), (1, 2, 3, 8, 9, 11), "T7"),
    ("a", 14, 11, ((3, 1, 2, 0), (3, 2, 1, 0), 1), (1, 2, 3, 8, 12, 14), "T10"),
]
ONE = {0: 1}


def mul(*xs):
    out = {0: 1}
    for x in xs:
        out = p_mul(out, x)
    return out


def P1(x):
    return p_add(x, ONE)


# --------------------------------------------------------------------------------------------
# the seventeen charts of pp. 2-4, as rational functions of (T, J, U, A, C, B, Q)
# --------------------------------------------------------------------------------------------
OUT_NAMES = "T,J,U,A,C,B,Q"
SRC_NAMES = "q,p,j,c,e,b,a"


def build_charts(E):
    """Returns the list of charts: dict(name, region, vals (RF for q,p,j,c,e,b,a), strict
    (names of the parameters that are > 0 on the domain), terms (the document's count), cert)."""
    Tp, Jp, Up, Ap, Cp, Bp, Qp = [RF(E, E.var(n)) for n in OUT_NAMES.split(",")]
    one = RF(E, 1)
    ch = []

    def add(name, region, vals, strict, terms, cert="nonneg", u0=None, bkind="free"):
        ch.append(dict(name=name, region=region, vals=vals, strict=strict, terms=terms, cert=cert, u0=u0,
                       bkind=bkind))
    # both centres below one: p = T/(1+T), j = J/(1+J)
    p_, j_ = Tp / (1 + Tp), Jp / (1 + Jp)
    add("A1", "p<1, j<1, u0<1", [Qp, p_, j_, Cp, (Up / (1 + Up)) / j_, Bp, Ap], "JUAQ", 3186, "square", "lt")
    add("A2", "p<1, j<1, u0=1", [Qp, p_, j_, Cp, one / j_, Bp, Ap], "JAQ", 1309, "square", "eq")
    e_ = (1 + Up) / j_
    add("A3", "p<1, j<1, u0>1, j<=b<e/U", [Qp, p_, j_, Cp, e_, j_ + (e_ / Up - j_) * Bp / (1 + Bp), Ap],
        "JUAQ", 12522, "nonneg", "gt", "upper")
    add("A4", "p<1, j<1, u0>1, b<j", [Qp, p_, j_, Cp, e_, j_ * Bp / (1 + Bp), Ap], "JUAQ", 7134, "factor", "gt",
        "lower")
    # j > 1: j = 1+J, p = T/((1+T) j)
    j_ = 1 + Jp
    p_ = Tp / ((1 + Tp) * j_)
    for nm, u0_, st, terms, u0 in (("B1", Up / (1 + Up), "JUAQ", 12416, "lt"), ("B2", one, "JAQ", 2564, "eq")):
        e_ = u0_ / j_
        a_ = (1 / (j_ - 1) - e_) * Ap / (1 + Ap)
        add(nm, "j>1, u0%s1" % ("<" if u0 == "lt" else "="), [Qp, p_, j_, Cp, e_, Bp, a_], st, terms, "nonneg", u0)
    D = 1 + Up + Ap
    e_ = (1 + Up / (Jp * D)) / j_
    a_ = Ap / (Jp * D * j_)
    b_ = (Jp * D + Up) * Bp / (j_ * Up * (1 + Bp))
    add("B3", "j>1, u0>1", [Qp, p_, j_, Cp, e_, b_, a_], "JUAQ", 31706, "TB-squares", "gt", "cap")
    # p > 1: p = 1+J, j = T/((1+T) p), c = C/(J(1+C))
    p_ = 1 + Jp
    j_ = Tp / ((1 + Tp) * p_)
    c_ = Cp / (Jp * (1 + Cp))
    for nm, u0_, st, terms, u0 in (("C1", Up / (1 + Up), "TJUAQ", 6354, "lt"), ("C2", one, "TJAQ", 2943, "eq"),
                                  ("C3", 1 + Up, "TJUAQ", 23198, "gt")):
        e_ = u0_ / j_
        b_ = e_ * Bp / ((u0_ - 1) * (1 + Bp)) if u0 == "gt" else Bp
        add(nm, "p>1, u0%s1" % {"lt": "<", "eq": "=", "gt": ">"}[u0], [Qp, p_, j_, c_, e_, b_, Ap], st, terms,
            "nonneg", u0, "cap" if u0 == "gt" else "free")
    # p = 1, j < 1
    p_, j_ = one, Jp / (1 + Jp)
    for nm, u0_, st, terms, u0 in (("D1", Up / (1 + Up), "JUAQ", 1077, "lt"), ("D2", one, "JAQ", 450, "eq"),
                                  ("D3", 1 + Up, "JUAQ", 2843, "gt")):
        e_ = u0_ / j_
        b_ = e_ * Bp / ((u0_ - 1) * (1 + Bp)) if u0 == "gt" else Bp
        add(nm, "p=1, j<1, u0%s1" % {"lt": "<", "eq": "=", "gt": ">"}[u0], [Qp, p_, j_, Cp, e_, b_, Ap], st, terms,
            "nonneg", u0, "cap" if u0 == "gt" else "free")
    # j = 1, p < 1
    p_, j_ = Tp / (1 + Tp), one
    add("E1", "j=1, p<1, u0<1", [Qp, p_, j_, Cp, Up / (1 + Up), Bp, Ap], "UAQ", 833, "nonneg", "lt")
    add("E2", "j=1, p<1, u0=1", [Qp, p_, j_, Cp, one, Bp, Ap], "AQ", 296, "nonneg", "eq")
    add("E3", "j=1, p<1, u0>1, 1<=b<e/U", [Qp, p_, j_, Cp, 1 + Up, 1 + Bp / (Up * (1 + Bp)), Ap], "UAQ", 1359,
        "nonneg", "gt", "upper")
    add("E4", "j=1, p<1, u0>1, b<1", [Qp, p_, j_, Cp, 1 + Up, Bp / (1 + Bp), Ap], "UAQ", None, "lead", "gt", "lower")
    return ch


def rf_sign(Ro, rf, strict):
    """'pos' if rf > 0 on the chart domain (parameters in `strict` > 0, the others >= 0) by the
    coefficient criterion: denominator with positive coefficients and a monomial in strict
    variables only, numerator with nonnegative coefficients and such a monomial; 'zero' if rf = 0;
    else 'unknown'."""
    E = rf.E
    if len(rf.n) == 0:
        return "zero"
    mask_names = [nm for nm in Ro.names if nm not in strict]

    def strict_mono(P_):
        for i in range(len(P_)):
            ex = P_.exponent_tuple(i)
            if all(int(ex[Ro.pos[nm]]) == 0 for nm in mask_names):
                return True
        return False
    okd = all(rf.d.coefficient(i) > 0 for i in range(len(rf.d))) and strict_mono(rf.d)
    okn = all(rf.n.coefficient(i) >= 0 for i in range(len(rf.n))) and strict_mono(rf.n)
    return "pos" if (okd and okn) else "unknown"


def main():
    numeric = "--no-numeric" not in sys.argv
    import sympy as sp
    import numpy
    import scipy
    import flint
    ck = M.Checker("Case T_38 = {1,2,3,8,11,12} (document 'T37'): independent exact replay")
    ck.info("python %s, sympy %s, python-flint %s (fmpz_mpoly arithmetic in the charts; integer determinants at "
            "test points), numpy %s (int64 in det_linear_forms; floats in NUMERIC checks), scipy %s (NUMERIC only)"
            % (sys.version.split()[0], sp.__version__, flint.__version__, numpy.__version__, scipy.__version__))

    # ------------------------------------------------------------------------------------
    ck.section("A. Provenance and foundations")
    M.check_document(ck, DOC, DOC_SHA256[DOC])
    M.check_foundations(ck)
    M.check_radial_necessity(ck)
    M.crosscheck_supplement_tables(ck)
    M.crosscheck_audit_library(ck, [38, 39, 4, 8, 11, 46])

    # ------------------------------------------------------------------------------------
    ck.section("B. Contact set (p. 1)")
    ck.check("EXACT", "the document's set {1,2,3,8,11,12} is the paper's T_38 (Appendix E) as a set of "
             "lexicographic indices: group element = identity; the document's number is the paper's minus one",
             tuple(sorted(T)) == M.PAPER_T_V[38])
    r38, r39 = M.orbit_rep(T), M.orbit_rep(M.PAPER_T_V[39])
    ck.check("EXACT", "orbit representative of T_38 under G is (0,1,2,6,11,22); T_39 lies in another orbit",
             r38 == (0, 1, 2, 6, 11, 22) and r39 != r38, "T_38 -> %s, T_39 -> %s" % (r38, r39))

    # ------------------------------------------------------------------------------------
    ck.section("C. The contact cone (p. 1)")
    gens = M.generators(T)
    ck.check("EXACT", "the f_j of Appendix D vanishing on T are exactly f0, f1, f2, f9, f10, f15, f53, f56 "
             "(the document's list); the only non-minor one is f56 = h_{4132,+}",
             gens == sorted(ORDER) and [j for j in gens if not M.F_IS_MINOR[j]] == [56]
             and M.f_name(56) == "h_{4132,+}", ", ".join("f%d = %s" % (j, M.f_name(j)) for j in gens))
    rk = M.rank([M.P[s] for s in T])
    ck.check("EXACT", "the six contact equations have rank 6 in V_R = Q^14, so their solution space has dimension 8; "
             "the eight generators are linearly independent, hence span it",
             rk == 6 and M.rank([M.F_COORD[j] for j in ORDER]) == 8)
    ck.check("EXACT", "the values of (f1, f2, f15, f9, f10, f56, f0, f53) at the permutations "
             "(17, 15, 4, 10, 9, 6, 20, 14) form the identity matrix: the coefficients (q,p,j,c,e,d,b,a) of a "
             "vector of the solution space are its gaps at these permutations, so N_T is the simplicial cone "
             "{q,...,a >= 0}",
             [[F_VALUES[j][s] for j in ORDER] for s in EXCL] == [[int(i == k) for k in range(8)] for i in range(8)])
    d, rays = M.extreme_rays_of_contact_cone(T)
    ck.check("EXACT", "independently: direct enumeration of the extreme rays of N_T = {Pc >= 0, (Pc)_T = 0} "
             "(dimension %d; no use of the 82-facet list) gives exactly these eight vectors" % d,
             d == 8 and rays == sorted(F_VALUES[j] for j in ORDER))
    M.check_minor_leibniz(ck, [j for j in ORDER if j != 56])

    # ------------------------------------------------------------------------------------
    ck.section("D. The four zero faces q = 0, j = 0, e = 0, a = 0 (p. 1): reduction to cited cases")
    for nm, extra, k, g, image, docname in ZERO_FACES:
        gm = M.g_index_map(g)
        img = tuple(sorted(gm[s] for s in M.PAPER_T_V[k]))
        S = set(T) | {extra}
        hits = [g2 for g2, gm2 in M.G_EXPLICIT if set(gm2[s] for s in M.PAPER_T_V[k]) <= S]
        subs = [sub for sub in __import__("itertools").combinations(sorted(S), 6)
                if M.orbit_rep(sub) == M.orbit_rep(M.PAPER_T_V[k])]
        ck.check("EXACT", "%s = 0 is the gap at permutation %d (generator f%d); the document's element "
                 "(row %s, col %s, tr %d) maps the paper's T_%d = %s (the document's '%s') onto %s, which is "
                 "contained in T + {%d}; the covariance Phi_{g.c}(g.M) = +-Phi_c(M) holds as a monomial identity"
                 % (nm, extra, ORDER[NAMES.split(",").index(nm)], g[0], g[1], g[2], k, list(M.PAPER_T_V[k]),
                    docname, list(img), extra),
                 EXCL[NAMES.split(",").index(nm)] == extra and img == image and set(img) <= S
                 and g in hits and M.generic_covariance(g),
                 "%d group element(s) do this; six-subsets of T + {%d} in the orbit of T_%d: %s"
                 % (len(hits), extra, k, subs))
    ck.info("PROOF-READ: if the coefficient vanishes, c(x) lies in N_{T + {extra}}, a subset of N_{g.T_k}; "
            "g^-1.x lies in X, its gap vector g^-1.c(x) lies in N_{T_k}, and Phi_{c(x)}(O) = Phi_{g^-1.c(x)}(O'') "
            "with O'' in SO(4); so Case T_k gives Phi >= 0.  Cited: T_46 (q), T_4 (j), T_8 (e), T_11 (a).")

    # ------------------------------------------------------------------------------------
    ck.section("E. The sufficient conditions (2)-(5) (pp. 1-2): orbit membership and algebra")
    rng = random.Random(38002)
    M.check_transport(ck, "SOC", (2, 15, 56), rng)
    M.check_transport(ck, "old-H", (9, 2, 56), rng)
    M.check_transport(ck, "old-H", (10, 15, 56), rng)
    M.check_transport(ck, "old-H", (53, 15, 56), rng)
    M.check_transport(ck, "chain", (15, 0, 10, 56), rng)
    a, b, c, e, j, p, dd = sp.symbols("a b c e j p d", positive=True)
    d1, d2 = dd * e / (a + e), dd * a / (a + e)
    ck.check("EXACT", "coefficient forms: SOC ab >= d^2 on (f2, f15; f56) is pj >= d^2, (2); old-H t(c-d) >= d^2 "
             "with t-role f9 (coefficient c), c-role f2 (coefficient p) is c(p-d) >= d^2, i.e. p >= 1 + 1/c at "
             "d = 1, (3); chain a >= d^2/p + d^2/b with centre f15, leaves f0, f10 is j >= 1/b + 1/e at d = 1, (5)",
             sp.cancel(c * (p - 1) - 1 - c * (p - (1 + 1 / c))) == 0)
    ck.check("EXACT", "(4): split d = d1 + d2 with d1 : d2 = e : a between the old-H members (f10, f15; f56) and "
             "(f53, f15; f56); they need the coefficients d1 + d1^2/e and d2 + d2^2/a of f15, whose sum is "
             "d + d^2/(a+e): so j >= d + d^2/(a+e) (j >= 1 + 1/(a+e) at d = 1) suffices, and its failure is "
             "d S + d^2 - j S > 0 with S = a + e",
             sp.cancel(d1 + d1 ** 2 / e + d2 + d2 ** 2 / a - (dd + dd ** 2 / (a + e))) == 0
             and sp.cancel(d1 + d2 - dd) == 0
             and sp.cancel((a + e) * (dd + dd ** 2 / (a + e) - j) - (dd * (a + e) + dd ** 2 - j * (a + e))) == 0)

    # ------------------------------------------------------------------------------------
    ck.section("F. Radial necessity in the generators (p. 2)")
    rd = M.radial_data(ck, ORDER, NAMES)
    R9, C0, Q0, G = rd["R"], rd["C0"], rd["Q"], rd["G"]
    va, ve, vd, vj = [R9.var(x) for x in "aedj"]
    S_ = p_add(va, ve)
    guard = p_add(p_mul(S_, C0), mul(va, ve, p_add(p_add(p_mul(vd, S_), mul(vd, vd)), p_mul(vj, S_), -1)), -1)
    ck.check("EXACT", "C0 has %d terms; with S = a + e, P := S C0 - a e (d S + d^2 - j S) has %d terms, all "
             "coefficients nonnegative, and contains a d^3 + e d^3 with coefficient 1.  So where (4) fails "
             "(d S + d^2 - j S > 0) and a, e, d > 0: S C0 >= (a+e) d^3 > 0, hence C0 > 0 and C0 + hQ = 0 forces "
             "Q != 0" % (len(C0), len(guard)),
             p_nonneg(guard) and guard.get(R9.key((0, 0, 0, 0, 0, 3, 0, 1, 0))) == 1
             and guard.get(R9.key((0, 0, 0, 0, 1, 3, 0, 0, 0))) == 1)
    d9 = R9.key((0, 0, 0, 0, 0, 9, 0, 0, 0))
    if G.get(d9, 0) < 0:
        G = p_scale(G, -1)
    ck.check("EXACT", "G_9 (sign fixed by coefficient +1 of d^9) has 1154 terms and is homogeneous of degree 9 "
             "(the document's count); every point of X in the hard region has G_9 = 0",
             len(G) == 1154 and G.get(d9) == 1 and all(sum(R9.exps(k)) == 9 for k in G))
    M.compare_radial_json(ck, rd, JSON_RADIAL)

    # ------------------------------------------------------------------------------------
    ck.section("G. The seventeen charts (pp. 2-4): definitions, regions, term counts, sign certificates")
    Rs = Ring(SRC_NAMES)
    G1 = p_rename(R9.subs_one(G, "d"), R9, Rs)
    ck.check("EXACT", "at d = 1, G_9 is a polynomial in (q,p,j,c,e,b,a) with 1154 terms and degrees %s"
             % ([Rs.deg(G1, n) for n in Rs.names],), len(G1) == 1154)
    Ro = Ring(OUT_NAMES)
    E = FL(Ro)
    charts = build_charts(E)
    one = RF(E, 1)
    polys = {}
    for chd in charts:
        nm, vals, strict = chd["name"], chd["vals"], set(chd["strict"])
        qv, pv, jv, cv, ev, bv, av = vals
        u0 = jv * ev
        # (i) the chart maps its parameter domain into the stated part of the hard region
        okreg = all(rf_sign(Ro, x, strict) == "pos" for x in (qv, jv, ev, av))          # q, j, e, a > 0
        okreg &= all(len(x.n) == 0 or x.num_nonneg() for x in (pv, cv, bv))             # p, c, b >= 0
        fails = [one - pv * jv,                      # (2) fails: pj < 1
                 one - cv * (pv - 1),                # (3) fails: c(p-1) < 1
                 one - (jv - 1) * (av + ev),         # (4) fails: (j-1)(a+e) < 1
                 ev - bv * (u0 - 1)]                 # (5) fails: b(je-1) < e
        okreg &= all(rf_sign(Ro, x, strict) == "pos" for x in fails)
        reg = chd["region"]
        U_, J_, T_ = RF(E, E.var("U")), RF(E, E.var("J")), RF(E, E.var("T"))
        okid = True
        if chd["u0"] == "lt":
            okid &= u0 == U_ / (1 + U_)
        elif chd["u0"] == "eq":
            okid &= u0 == one
        else:
            okid &= u0 == 1 + U_ or nm == "B3"
        if nm[0] == "A":
            okid &= pv == T_ / (1 + T_) and jv == J_ / (1 + J_)
        elif nm[0] == "B":
            okid &= jv == 1 + J_ and pv * jv == T_ / (1 + T_)
        elif nm[0] == "C":
            okid &= pv == 1 + J_ and pv * jv == T_ / (1 + T_) and one / J_ - cv == one / (J_ * (1 + RF(E, E.var("C"))))
        elif nm[0] == "D":
            okid &= pv == one and jv == J_ / (1 + J_)
        else:
            okid &= jv == one and pv == T_ / (1 + T_)
        Bv, Av = RF(E, E.var("B")), RF(E, E.var("A"))
        if chd["bkind"] == "cap":                    # b = cap B/(1+B), cap = e/(u0-1)
            cap = ev / (u0 - 1)
            okid &= bv == cap * Bv / (1 + Bv) and rf_sign(Ro, cap, strict) == "pos"
        elif chd["bkind"] == "upper":                # b = lo + (cap - lo) B/(1+B), lo = j
            cap = ev / (u0 - 1)
            okid &= bv - jv == (cap - jv) * Bv / (1 + Bv) and rf_sign(Ro, cap - jv, strict) == "pos"
        elif chd["bkind"] == "lower":                # b = j B/(1+B)
            okid &= bv == jv * Bv / (1 + Bv)
        else:
            okid &= bv == Bv
        if nm in ("B1", "B2"):                       # a = (1/(j-1) - e) A/(1+A)
            room = one / (jv - 1) - ev
            okid &= av == room * Av / (1 + Av) and rf_sign(Ro, room, strict) == "pos"
        elif nm == "B3":                             # je - 1 = U/(JD), ja = A/(JD), 1 - J(je - 1 + ja) = 1/D
            D_ = 1 + U_ + Av
            okid &= u0 - 1 == U_ / (J_ * D_) and jv * av == Av / (J_ * D_) \
                and one - J_ * (u0 - 1 + jv * av) == one / D_
        else:
            okid &= av == Av
        okden = all(rf_sign(Ro, RF(E, x.d), strict) == "pos" for x in vals)
        ck.check("EXACT", "chart %s (%s): the substituted q, p, j, c, e, b, a are the document's formulas; on the "
                 "domain {%s > 0, the other parameters >= 0} all denominators are > 0, q, j, e, a > 0, "
                 "p, c, b >= 0 and the four conditions (2)-(5) fail (each difference is a quotient of polynomials "
                 "with nonnegative coefficients having a monomial in the strictly positive parameters)"
                 % (nm, reg, ", ".join(sorted(strict))), okreg and okid and okden)
        # (ii) the cleared polynomial
        subs = {s_: v.pair() for s_, v in zip(Rs.names, vals)}
        out, degs = p_compose_fl(G1, Rs, E, subs)
        polys[nm] = out
        chd["degs"] = degs
    # cross-check of the flint composition against the pure Python composition
    okx = True
    for nm in ("A1", "D3", "E3"):
        chd = next(c_ for c_ in charts if c_["name"] == nm)
        subs = {s_: (E.from_fl(v.n), E.from_fl(v.d)) for s_, v in zip(Rs.names, chd["vals"])}
        ref, _ = p_compose(G1, Rs, Ro, subs)
        okx &= ref == polys[nm]
    ck.check("EXACT", "the python-flint composition agrees with the pure Python composition p_compose of PART 2 "
             "on the charts A1, D3, E3 (same numerators, same denominators, same clearing)", okx)
    ck.info("clearing convention (learned by reading the supplement's verifier as text, not executed): every "
            "variable x of G_9 is replaced by its reduced fraction n_x/d_x and the result is multiplied by "
            "prod_x d_x^(deg_x G_9); degrees of G_9 in (q,p,j,c,e,b,a): %s" % ([Rs.deg(G1, n) for n in Rs.names],))

    key = Ro.key
    Tv, Jv, Uv, Av_, Cv, Bv_, Qv = [Ro.var(n) for n in Ro.names]
    rems = {}                       # chart -> coefficientwise nonnegative remainder (a positive multiple of it)
    for chd in charts:
        nm, strict = chd["name"], chd["strict"]
        poly = polys[nm]
        nneg = sum(1 for v in poly.values() if v < 0)
        if chd["cert"] == "nonneg":
            rems[nm] = poly
            w = M.strict_witness(Ro, poly, strict)
            ck.check("EXACT", "chart %s: the cleared polynomial has %d terms (document: %d), all coefficients "
                     "nonnegative, and a monomial in the strictly positive parameters only: G_9 > 0 on the chart"
                     % (nm, len(poly), chd["terms"]),
                     len(poly) == chd["terms"] and nneg == 0 and w is not None,
                     "strict monomial %s * %s" % (poly.get(w), M.mono_str(Ro, w)) if w is not None else "")
        elif chd["cert"] == "square":
            Bm1 = p_add(Bv_, ONE, -1)
            sq = mul(Jv, Jv, P1(Jv), Cv, Cv, Qv, Qv, Bm1, Bm1)
            if nm == "A1":
                sq = mul(sq, Uv, Uv)
            rem2 = p_add(p_scale(poly, 2), sq, -1)          # 2 (poly - sq/2)
            rems[nm] = rem2
            wkey = key((0, 0, 2, 0, 0, 0, 0)) if nm == "A1" else 0
            ck.check("EXACT", "chart %s: the cleared polynomial has %d terms (document: %d) and exactly two "
                     "negative coefficients; poly - (1/2) J^2 (J+1) %sC^2 Q^2 (B-1)^2 has nonnegative coefficients "
                     "and contains %s with coefficient 1: G_9 > 0 on the chart"
                     % (nm, len(poly), chd["terms"], "U^2 " if nm == "A1" else "",
                        "U^2" if nm == "A1" else "the constant term"),
                     len(poly) == chd["terms"] and nneg == 2 and p_nonneg(rem2) and rem2.get(wkey) == 2,
                     "negative terms: %s" % sorted((M.mono_str(Ro, k), v) for k, v in poly.items() if v < 0))
        elif chd["cert"] == "TB-squares":
            direction = key((1, 0, 0, 0, 0, 1, 0))
            rem, sqs = M.absorb_binomial_squares(poly, direction)
            ok = rem is not None
            recon = {}
            if ok:
                recon = dict(rem)
                for k, t_, m_, r_ in sqs:
                    ok &= m_ * m_ == 4 * t_ * r_ and t_ > 0 and r_ > 0
                    for kk, vv in ((k - direction, t_), (k, -m_), (k + direction, r_)):
                        recon[kk] = recon.get(kk, 0) + vv
                recon = {k: v for k, v in recon.items() if v}
                ok &= recon == {k: Fraction(v) for k, v in poly.items()}
                ok &= all(v >= 0 for v in rem.values())
            w = M.strict_witness(Ro, rem, strict) if ok else None
            if ok:
                rems[nm] = rem
            ck.check("EXACT", "chart %s: the cleared polynomial has %d terms (document: %d) and %d negative "
                     "coefficients (document: twenty); each is absorbed by a quadratic "
                     "x/(TB) (t - m TB + r (TB)^2) with m^2 = 4 t r (so it is t x/(TB) (1 - sqrt(r/t) TB)^2 >= 0); "
                     "the polynomial equals the remainder plus the twenty quadratics, the remainder has "
                     "nonnegative coefficients and contains 3 J U^4 A^2 (the document's example)"
                     % (nm, len(poly), chd["terms"], nneg),
                     ok and len(poly) == chd["terms"] and nneg == 20 and len(sqs) == 20 and w is not None
                     and rem.get(key((0, 1, 4, 2, 0, 0, 0))) == 3)
            if ok:
                ck.info("the twenty allocations (own; the supplement's JSON with the author's allocations is not "
                        "public), as (monomial x, t, m, r): %s"
                        % "; ".join("(%s, %s, %s, %s)" % (M.mono_str(Ro, k), t_, m_, r_) for k, t_, m_, r_ in sqs))
    # the lower-b chart A4 and the layer E4
    Z = polys["A4"]
    sA, sJ = SH * Ro.pos["A"], SH * Ro.pos["J"]
    Z0 = {k: v for k, v in Z.items() if not (k >> sA) & MASK}
    ZA = {k: v for k, v in Z.items() if (k >> sA) & MASK}
    ck.check("EXACT", "chart A4: the cleared polynomial Z has %d terms (document: 7134); its part Z_A of positive "
             "degree in A has nonnegative coefficients (%d terms); its A^0 part Z_0 has %d terms"
             % (len(Z), len(ZA), len(Z0)), len(Z) == 7134 and p_nonneg(ZA) and len(Z0) == 2330)
    q1 = p_divexact(Z0, P1(Jv))
    okf = q1 is not None
    P5 = R11 = F11 = None
    sqfac = None
    if okf:
        g = None
        for _, part in Ro.split(q1, "T").items():
            pf = E.to_fl(part)
            g = pf if g is None else g.gcd(pf)
        P5 = E.from_fl(g)
        if max(P5.values()) < 0:
            P5 = p_scale(P5, -1)
        F11 = p_divexact(q1, P5)
        okf = F11 is not None
    if okf:
        lin = p_add(p_add(p_add(mul(Tv, Jv, Bv_), mul(Tv, Bv_)), Tv), Jv, -1)       # TJB + TB + T - J
        sqfac = mul(Uv, Cv, Qv, P1(Jv), P1(Jv), lin, lin)
        R11 = p_add(F11, sqfac, -1)
        recon = p_add(ZA, mul(P1(Jv), P5, p_add(R11, sqfac)))
        okf = recon == Z
    ck.check("EXACT", "chart A4: Z = Z_A + (J+1) P5 [R11 + U C Q (J+1)^2 (TJB + TB + T - J)^2] (verified by "
             "multiplying out), where P5 (found here as the gcd of the T-coefficients of Z_0/(J+1)) has 40 terms "
             "and R11 has 309 terms (the document's counts), all coefficients of P5 and R11 are positive, P5 "
             "contains 2 J^2 and R11 contains 5 J^4.  Hence Z >= 10 J^6 (J+1) > 0 for J > 0 and all "
             "T, U, A, C, B, Q >= 0",
             okf and len(P5) == 40 and len(R11) == 309 and all(v > 0 for v in P5.values())
             and all(v > 0 for v in R11.values()) and P5.get(key((0, 2, 0, 0, 0, 0, 0))) == 2
             and R11.get(key((0, 4, 0, 0, 0, 0, 0))) == 5 and Ro.tdeg(P5) == 5 and Ro.tdeg(F11) == 11,
             "P5 has total degree 5 and does not involve T, A; the cofactor F11 = R11 + square has total degree "
             "11 and %d terms, %d of them negative; R11 has total degree %d"
             % (len(F11 or {}), sum(1 for v in (F11 or {}).values() if v < 0), Ro.tdeg(R11 or {})))
    Z1 = polys["E4"]

    def lead(pp, deg):
        assert Ro.deg(pp, "J") <= deg
        return {k - (deg << sJ): v for k, v in pp.items() if ((k >> sJ) & MASK) == deg}
    okl = okf
    if okf:
        Pl, Rl, Al = lead(P5, 2), lead(R11, 4), lead(ZA, 7)
        tb1 = p_add(mul(Tv, Bv_), ONE, -1)
        okl = lead(Z, 7) == Z1 and Z1 == p_add(Al, mul(Pl, p_add(Rl, mul(Uv, Cv, Qv, tb1, tb1)))) \
            and p_nonneg(Al) and p_nonneg(Pl) and p_nonneg(Rl) and Pl.get(0) == 2 and Rl.get(0) == 5
    ck.check("EXACT", "chart E4 (j = 1, e = 1+U, b = B/(1+B)): the cleared polynomial Z1 (%d terms) equals the "
             "coefficient of J^7 in Z, and Z1 = A_lead + P_lead [R_lead + U C Q (TB - 1)^2] with A_lead, P_lead, "
             "R_lead the J-leading coefficients of Z_A, P5, R11; all three have nonnegative coefficients and "
             "P_lead, R_lead have constant terms 2, 5.  Hence Z1 >= 10 > 0 for all T, U, A, C, B, Q >= 0"
             % len(Z1), okl)
    ck.info("PROOF-READ (coverage): a hard point has q, j, e, a > 0, p, c, b >= 0 (d = 1) and (2)-(5) all fail.  "
            "(2) fails: pj < 1, so exactly one of: p<1 & j<1; j>1 (then p < 1/j); p>1 (then j < 1/p); p = 1 (j<1); "
            "j = 1 (p<1).  (3) can hold only for p > 1 and fails there iff c < 1/(p-1); (4) can hold only for j > 1 "
            "and fails there iff a + e < 1/(j-1); (5) can hold only for u0 = je > 1 and fails there iff "
            "b < e/(u0-1).  Each chart has an explicit inverse on its part of the region (tested on random points "
            "in section H), so the seventeen charts cover the hard region; in each, G_9 > 0.  With section F this "
            "contradicts G_9 = 0: no point of X lies in the hard region, so one of (2)-(5) holds, or a "
            "coefficient among q, j, e, a, d vanishes (sections C, D).")

    # ------------------------------------------------------------------------------------
    ck.section("I. Observation (not in the document): the four zero faces need no cited case")
    Cplus = p_add(C0, mul(R9.var("j"), ve, va))
    ck.check("EXACT", "C0 + j e a has only positive coefficients and contains d^3: the only negative term of C0 is "
             "-j e a.  So C0 >= d^3 > 0 whenever j e a = 0 and d > 0; with section F, C0 > 0 at every point with "
             "d > 0 at which (4) fails or one of j, e, a vanishes",
             all(v > 0 for v in Cplus.values()) and Cplus.get(R9.key((0, 0, 0, 0, 0, 3, 0, 0, 0))) == 1
             and len(Cplus) == len(C0) - 1)
    Gj0 = R9.subs_zero(G, ["j"])
    ck.check("EXACT", "face j = 0: G_9 restricted to j = 0 has %d terms, all coefficients nonnegative, and "
             "contains d^9 with coefficient 1.  So G_9 >= d^9 > 0 there, while C0 > 0 forces G_9 = 0: no point "
             "of X with c(x) in N_T has j = 0 < d.  The citation of T_4 is not needed" % len(Gj0),
             p_nonneg(Gj0) and Gj0.get(d9) == 1)
    okq = oke = True
    wq, we = [], []
    for chd in charts:
        nm, strict = chd["name"], chd["strict"]
        # the chart denominators stay positive at Q = 0, A = 0 (and at U = 0 in the charts with u0 < 1)
        st_qa = set(strict) - set("QA")
        okq &= all(rf_sign(Ro, RF(E, x.d), st_qa) == "pos" for x in chd["vals"])
        if chd["u0"] == "lt":
            oke &= all(rf_sign(Ro, RF(E, x.d), st_qa - set("U")) == "pos" for x in chd["vals"])
        if nm in ("A4", "E4"):
            continue                # Z >= 10 J^6 (J+1) and Z1 >= 10 hold for all T, U, A, C, B, Q >= 0
        w1 = M.strict_witness(Ro, rems[nm], strict.replace("Q", "").replace("A", ""))
        okq &= w1 is not None
        wq.append("%s: %s" % (nm, M.mono_str(Ro, w1) if w1 is not None else "-"))
        if chd["u0"] == "lt":
            w2 = M.strict_witness(Ro, rems[nm], strict.replace("Q", "").replace("A", "").replace("U", ""))
            oke &= w2 is not None
            we.append("%s: %s" % (nm, M.mono_str(Ro, w2) if w2 is not None else "-"))
    ck.check("EXACT", "faces q = 0 and a = 0 (separately or together): the charts are defined at Q = 0 and at "
             "A = 0 (all denominators stay positive) and still parametrize the hard points with q = 0 or a = 0 "
             "(same inverses; tested on random points in section H); in each of the fifteen direct "
             "charts the nonnegative remainder has a positive monomial in the strictly positive parameters "
             "other than Q and A (and the bounds of A4, E4 do not use Q or A).  So G_9 > 0 on the hard region "
             "also for q = 0 or a = 0: the citations of T_46 and T_11 are not needed", okq, "; ".join(wq))
    ck.check("EXACT", "face e = 0 (u0 = 0): the five charts with u0 < 1 are defined at U = 0 (all denominators "
             "stay positive), where e = 0, and "
             "parametrize the hard points with e = 0 (there (5) is unavailable and (4) reads j < 1 + 1/a); "
             "their remainders have a positive monomial free of U, Q and A.  So G_9 > 0 on the hard region also "
             "for e = 0 (alone or together with q = 0, a = 0): the citation of T_8 is not needed", oke,
             "; ".join(we))

    if numeric:
        numeric_part(ck, rd, G, G1, Rs, Ro, E, charts, polys, rems)
    ck.finish()


# ==========================================================================================
def classify(pt):
    """pt: dict of Fractions q,p,j,c,e,b,a (d = 1), a hard point.  Returns (chart name, parameters)
    by the explicit inverse of the chart, all in exact arithmetic."""
    q, p, j, c, e, b, a = [pt[n] for n in "qpjceba"]
    u0 = j * e
    par = dict(Q=q, A=a, C=c, B=b, U=Fraction(0), T=Fraction(0), J=Fraction(0))
    if u0 < 1:
        par["U"] = u0 / (1 - u0)
        us = "lt"
    elif u0 == 1:
        us = "eq"
    else:
        par["U"] = u0 - 1
        us = "gt"
    idx = {"lt": "1", "eq": "2", "gt": "3"}[us]
    if p < 1 and j < 1:
        par["T"], par["J"] = p / (1 - p), j / (1 - j)
        if us != "gt":
            return "A" + idx, par
        cap = e / (u0 - 1)
        if b >= j:
            par["B"] = (b - j) / (cap - b)
            return "A3", par
        par["B"] = b / (j - b)
        return "A4", par
    if j > 1:
        J = j - 1
        par["J"], par["T"] = J, p * j / (1 - p * j)
        if us != "gt":
            room = 1 / J - e
            par["A"] = a / (room - a)
            return "B" + idx, par
        x, y = u0 - 1, j * a
        Dv = 1 / (1 - J * (x + y))
        par["U"], par["A"] = x * J * Dv, y * J * Dv
        cap = e / (u0 - 1)
        par["B"] = b / (cap - b)
        return "B3", par
    if p > 1:
        J = p - 1
        par["J"], par["T"] = J, p * j / (1 - p * j)
        par["C"] = c * J / (1 - c * J)
        if us == "gt":
            cap = e / (u0 - 1)
            par["B"] = b / (cap - b)
        return "C" + idx, par
    if p == 1:
        par["J"] = j / (1 - j)
        if us == "gt":
            cap = e / (u0 - 1)
            par["B"] = b / (cap - b)
        return "D" + idx, par
    assert j == 1 and p < 1
    par["T"] = p / (1 - p)
    if us != "gt":
        return "E" + idx, par
    if b >= 1:
        par["B"] = (b - 1) / (1 / par["U"] - (b - 1))
        return "E3", par
    par["B"] = b / (1 - b)
    return "E4", par


def numeric_part(ck, rd, G, G1, Rs, Ro, E, charts, polys, rems):
    import numpy as np
    ck.section("H. NUMERIC falsification attempts (evidence only)")
    f = lambda idx: np.array(F_VALUES[idx], float)
    rng = random.Random(38003)
    ex = lambda: math.exp(rng.gauss(0.0, 1.0))
    worst = {}

    def upd(name, vec, seed):
        worst[name] = min(worst.get(name, 1e9), M.numeric_min_phi(vec, seed, 6000, 8))
    for trial in range(10):
        a, b, c, e, p = [ex() for _ in range(5)]
        upd("SOC (f2, f15; f56), pj = 1", p * f(2) + (1 / p) * f(15) + f(56), 100 + trial)
        upd("old-H (f9, f2; f56), p = 1 + 1/c", c * f(9) + (1 + 1 / c) * f(2) + f(56), 200 + trial)
        upd("old-H (f10, f15; f56), j = 1 + 1/e", e * f(10) + (1 + 1 / e) * f(15) + f(56), 300 + trial)
        upd("old-H (f53, f15; f56), j = 1 + 1/a", a * f(53) + (1 + 1 / a) * f(15) + f(56), 400 + trial)
        upd("(4): e f10 + a f53 + (1 + 1/(a+e)) f15 + f56", e * f(10) + a * f(53) + (1 + 1 / (a + e)) * f(15) + f(56),
            500 + trial)
        upd("chain (f15; f0, f10; f56), j = 1/e + 1/b", (1 / e + 1 / b) * f(15) + b * f(0) + e * f(10) + f(56),
            600 + trial)
    for nm, wv in worst.items():
        ck.check("NUMERIC", "%s (coefficients exp(N(0,1)), 10 draws, seed 38003; 6000 random rotations + 8 BFGS "
                 "descents each): min Phi found >= -1e-9" % nm, wv > -1e-9, "min = %.3e" % wv)
    ck.info("controls (the search does detect negative values): min Phi found for f56 alone %.4f; for "
            "f2 + 0.5 f15 + f56 (pj = 0.5 < 1) %.4f; for f9 + f2 + f56 (p = 1 < 1 + 1/c = 2) %.4f; for "
            "f0 + f10 + 0.2 f15 + f56 (j = 0.2 < 1/e + 1/b = 2) %.4f"
            % (M.numeric_min_phi(f(56), 96, 6000, 8), M.numeric_min_phi(f(2) + 0.5 * f(15) + f(56), 97, 6000, 8),
               M.numeric_min_phi(f(9) + f(2) + f(56), 98, 6000, 8),
               M.numeric_min_phi(f(0) + f(10) + 0.2 * f(15) + f(56), 99, 6000, 8)))

    # ---- exact rational sampling of the hard region in the ORIGINAL coefficients
    R9 = rd["R"]
    cd = {c_["name"]: {s_: (E.from_fl(v.n), E.from_fl(v.d)) for s_, v in zip(Rs.names, c_["vals"])} for c_ in charts}
    strict = {c_["name"]: set(c_["strict"]) for c_ in charts}
    rng = random.Random(38004)

    def pos(hi=3, den=20):
        return Fraction(rng.randint(1, hi * den), den)

    def frac(lo=0):
        return Fraction(rng.randint(lo, 99), 100)

    def gen():
        """A random rational point; the first draw only steers the distribution towards one of the
        five ranges of (p, j).  Whether the point is hard is decided afterwards by hard()."""
        kind = rng.choice("AABBCCDE")
        if kind == "A":
            p = Fraction(0) if rng.random() < 0.1 else frac()
            j = frac(1)
        elif kind == "B":
            j = 1 + pos(2)
            p = (Fraction(0) if rng.random() < 0.1 else frac()) / j
        elif kind == "C":
            p = 1 + pos(2)
            j = frac(1) / p
        elif kind == "D":
            p, j = Fraction(1), frac(1)
        else:
            p, j = (Fraction(0) if rng.random() < 0.1 else frac()), Fraction(1)
        r = rng.random()
        u0 = Fraction(1) if r < 0.15 else (frac(1) if r < 0.5 else 1 + pos(2))
        e = u0 / j
        a = pos(4)
        if j > 1 and rng.random() < 0.9:                       # pull a + e below 1/(j-1) when possible
            room = 1 / (j - 1) - e
            if room > 0:
                a = room * frac(1)
        r = rng.random()
        b = Fraction(0) if r < 0.1 else (j if r < 0.18 else (Fraction(1) if r < 0.24 else pos(4)))
        if u0 > 1 and rng.random() < 0.7:
            b = e / (u0 - 1) * frac()
        r = rng.random()
        c = Fraction(0) if r < 0.15 else pos(4)
        if p > 1 and rng.random() < 0.9:
            c = frac() / (p - 1)
        q = pos(4)
        return dict(q=q, p=p, j=j, c=c, e=e, b=b, a=a)

    def hard(pt):
        """All of (2)-(5) fail; a condition with a zero denominator is unavailable."""
        q, p, j, c, e, b, a = [pt[n] for n in "qpjceba"]
        return not (p * j >= 1 or (c > 0 and p >= 1 + 1 / c) or (a + e > 0 and j >= 1 + 1 / (a + e))
                    or (b > 0 and e > 0 and j >= 1 / e + 1 / b))

    def sample(n_target, faces):
        """Random hard points (faces: some of q, a, e are set to 0).  For each: the explicit inverse
        of the chart, the domain of the parameters, the chart image, and the sign of G_9."""
        counts = {c_["name"]: 0 for c_ in charts}
        bad_inv = bad_dom = bad_g = 0
        n_hard = tries = 0
        minval = None
        while n_hard < n_target and tries < 100 * n_target:
            tries += 1
            pt = gen()
            if faces:
                zero = [n_ for n_ in "qae" if rng.random() < 0.45] or [rng.choice("qae")]
                for n_ in zero:
                    pt[n_] = Fraction(0)
            if not hard(pt):
                continue
            n_hard += 1
            nm, par = classify(pt)
            counts[nm] += 1
            need = strict[nm] - set("QAU") if faces else strict[nm]
            if any(par[v] <= 0 for v in need) or any(par[v] < 0 for v in par):
                bad_dom += 1
                continue
            for s_ in Rs.names:
                num, den = cd[nm][s_]
                if Fraction(Ro.evaluate(num, par)) != pt[s_] * Fraction(Ro.evaluate(den, par)):
                    bad_inv += 1
                    break
            # G_9 at the point, through the homogeneous polynomial at an integer point
            N = 1
            for v in pt.values():
                N = N * v.denominator // math.gcd(N, v.denominator)
            ipt = {n_: int(pt[n_] * N) for n_ in "qpjceba"}
            ipt["d"] = N
            ipt["h"] = 0
            val = R9.evaluate(G, ipt)
            if val <= 0:
                bad_g += 1
            fv = float(Fraction(val, N ** 9))
            minval = fv if minval is None else min(minval, fv)
        return n_hard, counts, minval, bad_dom, bad_inv, bad_g
    n_hard, counts, minval, bad_dom, bad_inv, bad_g = sample(4000, False)
    ck.check("NUMERIC", "exact evaluation at %d random rational hard points in the original coefficients (seed 38004; "
             "q, j, e, a > 0, p, c, b >= 0 with p = 0, p = 1, j = 1, je = 1, b = 0, b = j, b = 1, c = 0 occurring; "
             "all of (2)-(5) fail): the explicit inverse of the chart selected by the case distinction gives "
             "parameters in the chart's domain, the chart maps them back to the point, and G_9 > 0"
             % n_hard, n_hard == 4000 and bad_inv == 0 and bad_dom == 0 and bad_g == 0
             and all(v > 0 for v in counts.values()),
             "points per chart: %s; smallest G_9 = %.3e; violations: domain %d, inverse %d, sign %d"
             % (counts, minval, bad_dom, bad_inv, bad_g))
    n_hard, counts, minval, bad_dom, bad_inv, bad_g = sample(1500, True)
    ck.check("NUMERIC", "the observation of section I, exact evaluation at %d random rational points with j > 0, at "
             "least one of q, a, e equal to 0, and all available conditions among (2)-(5) failing: the same "
             "inverses give parameters >= 0 (J > 0, and T > 0 in the charts C), the chart maps them back to the "
             "point, and G_9 > 0" % n_hard,
             n_hard == 1500 and bad_inv == 0 and bad_dom == 0 and bad_g == 0,
             "points per chart: %s; smallest G_9 = %.3e; violations: domain %d, inverse %d, sign %d"
             % ({k: v for k, v in counts.items() if v}, minval, bad_dom, bad_inv, bad_g))
    okj = True
    for _ in range(1000):
        pt = {n_: Fraction(rng.randint(0, 40), rng.randint(1, 8)) for n_ in "qpcedba"}
        pt["d"] += Fraction(1, 8)
        pt["j"] = pt["h"] = Fraction(0)
        okj &= R9.evaluate(G, pt) >= pt["d"] ** 9
    ck.check("NUMERIC", "the observation of section I, exact evaluation at 1000 random rational points with j = 0 < d: "
             "G_9 >= d^9", okj)

    # ---- the cleared polynomials at random points of the closed chart domains
    rng = random.Random(38006)
    okv = True
    npts = 0
    for chd in charts:
        nm, st = chd["name"], set(chd["strict"])
        poly = polys[nm]
        w = M.strict_witness(Ro, rems[nm], chd["strict"]) if nm in rems else None
        for _ in range(30):
            par = {v: (rng.randint(1, 4) if v in st else rng.choice((0, 0, 1, 2, 3))) for v in Ro.names}
            val = Ro.evaluate(poly, par)
            npts += 1
            if nm == "A4":
                okv &= val >= 10 * par["J"] ** 6 * (par["J"] + 1)
            elif nm == "E4":
                okv &= val >= 10
            else:
                scale = 2 if chd["cert"] == "square" else 1
                okv &= scale * val >= rems[nm][w] * Ro.evaluate({w: 1}, par) > 0
    ck.check("NUMERIC", "exact evaluation of the seventeen cleared polynomials at %d random integer points of "
             "their closed parameter domains (seed 38006; T, C, B = 0 occur): each is >= its strict monomial "
             "(> 0), Z >= 10 J^6 (J+1) on A4 and Z1 >= 10 on E4" % npts, okv)

    # ---- the guard of section F at random points
    rng = random.Random(38005)
    C0 = rd["C0"]
    okg = True
    for _ in range(2000):
        pt = {n_: Fraction(rng.randint(0, 40), rng.randint(1, 8)) for n_ in "qpjcedba"}
        pt["a"] += Fraction(1, 8)
        pt["e"] += Fraction(1, 8)
        pt["d"] += Fraction(1, 8)
        pt["h"] = Fraction(0)
        S = pt["a"] + pt["e"]
        if pt["d"] * S + pt["d"] ** 2 - pt["j"] * S > 0:
            okg &= R9.evaluate(C0, pt) > 0
    ck.check("NUMERIC", "exact evaluation at 2000 random rational points (seed 38005): C0 > 0 whenever a, e, d > 0 "
             "and (4) fails", okg)

    # ---- genuine points of X
    eq = M.numeric_equality_points(T, 38009, 80)
    fit, cub, hes, qmin, gmax = M.numeric_radial_test(rd, ORDER, eq)
    ck.check("NUMERIC", "end-to-end test of the necessary equations on %d generic points x of X with "
             "Re(F_sigma - F_pi*) = 0 on T (seed 38009; labels apart by > 0.05; the other gaps have either sign): "
             "relative |C0 + H Q| < 1e-9 and relative |G_9| < 1e-9" % len(eq),
             len(eq) >= 40 and fit < 1e-9 and cub < 1e-9 and gmax < 1e-9,
             "fit %.1e, cubic %.1e, |QG| %.1e, |G| %.1e, smallest relative |Q| %.1e" % (fit, cub, hes, gmax, qmin))
    pts = []
    for k, eps in enumerate((0.02, 0.1, 0.3)):
        pts += M.numeric_stratum_search(T, 38100 + k, 250, eps=eps)
    absG = {k: abs(v) for k, v in G.items()}
    sel = []
    for p_ in pts:
        co = np.array([p_["gaps"][s] for s in EXCL], float)
        if co[5] > 1e-4:
            sel.append((p_, co))
    relmax = 0.0
    concl = {"SOC (2)": 0, "old-H (3)": 0, "(4)": 0, "chain (5)": 0, "a zero face q, j, e or a (tolerance 1e-6)": 0}
    okc = True
    n_strict = 0
    faces = {"q": 0, "j": 0, "e": 0, "a": 0}
    for p_, co in sel:
        v = co / co[5]
        val = dict(zip(R9.names[:8], [float(x) for x in v]))
        val["h"] = 0.0
        relmax = max(relmax, abs(R9.evaluate(G, val)) / R9.evaluate(absG, {k: abs(x) for k, x in val.items()}))
        q, p, j, c, e, d_, b, a = v
        hit = False
        if p * j >= 1 - 1e-7:
            concl["SOC (2)"] += 1
            hit = True
        if c > 1e-9 and p >= 1 + 1 / c - 1e-7:
            concl["old-H (3)"] += 1
            hit = True
        if a + e > 1e-9 and j >= 1 + 1 / (a + e) - 1e-7:
            concl["(4)"] += 1
            hit = True
        if b > 1e-9 and e > 1e-9 and j >= 1 / e + 1 / b - 1e-7:
            concl["chain (5)"] += 1
            hit = True
        for nm_, x_ in (("q", q), ("j", j), ("e", e), ("a", a)):
            faces[nm_] += int(x_ < 1e-6)
        if min(q, j, e, a) < 1e-6:
            concl["a zero face q, j, e or a (tolerance 1e-6)"] += 1
        else:
            n_strict += 1
        okc &= hit or min(q, j, e, a) < 1e-6
    ck.info("stratum (gaps >= 0, zero on T, sum 1): %d points from two-stage LM searches (seeds 38100..38102, 250 "
            "starts each, repulsion 0.02, 0.1, 0.3; residual < 1e-10); %d have d > 1e-4 (largest d = %.2e); among "
            "these the smallest normalised distance of two labels is %.2e, %d have q, j, e, a > 1e-6 d (numbers of "
            "points with a vanishing coefficient: %s), and the "
            "numbers of contacts (|gap| < 1e-7) are %s"
            % (len(pts), len(sel), max([co[5] for _, co in sel] or [float("nan")]),
               min([p_["sep"] for p_, _ in sel] or [float("nan")]), n_strict, faces,
               sorted(set(sum(1 for x_ in p_["gaps"] if abs(x_) < 1e-7) for p_, _ in sel))))
    ck.check("NUMERIC", "at every stratum point found with d > 1e-4 (%d points): relative |G_9| < 1e-6, and one of "
             "(2)-(5) holds or the point lies on a zero face%s"
             % (len(sel), "" if sel else "  [VOID: no such point was found]"),
             okc and relmax < 1e-6, "max relative |G_9| %.1e; which conditions hold: %s" % (relmax, concl))


if __name__ == "__main__":
    main()
