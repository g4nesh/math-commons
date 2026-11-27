"""Independent exact replay of Case T_39 = {2,3,5,7,9,12} (group V) of D. Zheng's proof of the
Marcus-de Oliveira conjecture for n = 4.

Proof document replayed: MOC4_SUPP/proofs/manin_transfer_20261021/blind_moc4_T38_uniform_20261022.pdf
(6 pages; the document calls the set "T38").  The supplement contains no program and no data
file for this document, so every polynomial and every table is rebuilt here from the
definitions (the matrix P, the cubic C, the vectors f_j of the paper) and from the formulas
printed in the PDF.  Nothing from the supplement is executed or imported.  If MOC4_SUPP is
set, the hash of the PDF is compared and the supplement's stored P, cubic and facet table are
compared with the rebuilt objects; no step depends on them.

Run:   MOC4_SUPP=/path/to/supplement  python replay_T39.py  [--no-numeric]
Exit status 0 iff every check passes.  Checks tagged EXACT use integer/rational/symbolic
arithmetic only; checks tagged NUMERIC are falsification attempts on random samples (floating
point, or exact evaluation at random rational points where the text says so).
"""
import math
import os
import random
import sys
from fractions import Fraction
from math import comb

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import replay_T38_common as M                                   # noqa: E402
from replay_T38_common import (Ring, FL, RF, p_add, p_mul, p_scale, p_nonneg, p_compose,   # noqa: E402
                               p_compose_fl, p_rename, p_divexact, F_VALUES, SH, MASK)

DIR = "proofs/manin_transfer_20261021/"
DOC = DIR + "blind_moc4_T38_uniform_20261022.pdf"
DOC_SHA256 = {
    DOC: "11a6c8fda813c4036f8e28f2733dcac6ddc53d1f7b44f179802ee8c6e20a396d",
}

T = (2, 3, 5, 7, 9, 12)
GENS = [0, 1, 2, 5, 7, 12, 39, 40, 67, 70]
ORDER40 = [39, 12, 5, 7, 2, 40, 67, 0]       # (1): q f39 + p f12 + j f5 + c f7 + e f2 + d f40 + b f67 + a f0
NAMES40 = "q,p,j,c,e,d,b,a,z"
ORDER70 = [0, 1, 2, 12, 39, 67, 70, 5]       # (2): p f0 + r f1 + a f2 + b f12 + c f39 + q f67 + d f70 + s f5
NAMES70 = "p,r,a,b,c,q,d,s,z"
IDENTITIES = [([40, 70], [0, 5, 12, 39, 67]), ([1, 40], [0, 7, 39]), ([7, 70], [1, 5, 12, 67])]
CORE_REP = (59, 53, 55, 67, 54)              # the "T24 core": j f59 + e f53 + b f55 + a f67 + d f54
CORE_USE = (5, 2, 67, 0, 40)                 # its image printed on p. 2
CORE_SWAP = (5, 2, 0, 67, 40)                # the same statement with the roles of a and b exchanged
# numbers of terms printed in the document
STRICT_COUNTS = {("lt", "lt"): (87, 118, 108, 239), ("lt", "gt"): (87, 124, 72, 190),
                 ("gt", "lt"): (72, 97, 99, 195), ("gt", "gt"): (72, 114, 93, 172)}
BERN_COUNTS = {("both_le1", "lt"): (17336, 23174, 23174), ("both_le1", "gt"): (46488, 62640, 62640),
               ("a_ge1", "lt"): (6722, 13776, 18895), ("a_ge1", "gt"): (9177, 19517, 26634),
               ("b_ge1", "lt"): (7143, 14455, 19835), ("b_ge1", "gt"): (24225, 42780, 53640),
               ("gt", "lt"): (3022, 5126, 7341), ("gt", "alloc"): (19910, 41864, 67816)}
STURM_GROUPS = {"middle": 308, "endpoint": 491}
KKT_COUNTS = {"HH": ([6372, 7776, 7757], [9408, 11466, 13524, 13502], [0, 15624, 18921, 21801, 21863]),
              "HC": ([6084, 7380, 7379], [7728, 9366, 11004, 11003], [0, 11514, 13575, 15591, 15599]),
              "CH": ([5628, 6738, 6737], [6601, 8694, 10094, 10093], [0, 9810, 12631, 14351, 14359]),
              "CC": ([6559, 8176, 8071], [6831, 8760, 10557, 10557], [0, 8850, 10959, 12936, 12936]),
              "CCswap": ([6212, 7663, 7588], [6138, 8281, 9887, 9887], [0, 7929, 10363, 12124, 12124])}
ONE = {0: 1}
TWO = {0: 2}


def mul(*xs):
    out = {0: 1}
    for x in xs:
        out = p_mul(out, x)
    return out


def P1(x):
    return p_add(x, ONE)


def reduce_to_one_nonminor(lam):
    """The document's reduction (p. 1).  lam: dict generator -> nonnegative Fraction.  Returns
    (kind, coefficients) with kind in {'minor', '40', '70'}; the vector is unchanged."""
    lam = {j: Fraction(lam.get(j, 0)) for j in GENS}

    def use(lhs, rhs, amount):
        for j in lhs:
            lam[j] -= amount
        for j in rhs:
            lam[j] += amount
    use(*IDENTITIES[0], min(lam[40], lam[70]))
    if lam[40] > 0:
        use(*IDENTITIES[1], min(lam[1], lam[40]))
        return ("40" if lam[40] > 0 else "minor"), lam
    if lam[70] > 0:
        use(*IDENTITIES[2], min(lam[7], lam[70]))
        return ("70" if lam[70] > 0 else "minor"), lam
    return "minor", lam


def functional_on_cone(targets):
    """w in Q^24, zero on T, with w . f_j = targets.get(j, 0) for all ten generators."""
    others = [s for s in range(24) if s not in T]
    A = [[Fraction(F_VALUES[j][s]) for s in others] + [Fraction(targets.get(j, 0))] for j in GENS]
    Rm, piv = M.rref(A)
    assert len(others) not in piv, "inconsistent targets"
    w = [Fraction(0)] * 24
    for r, pc in enumerate(piv):
        w[others[pc]] = Rm[r][-1]
    assert all(sum(w[s] * F_VALUES[j][s] for s in range(24)) == targets.get(j, 0) for j in GENS)
    return w


def check_roles(ck, label, rep, target, rng, trials=4):
    """EXACT: some group element maps f_rep[m] to f_target[m] for all m simultaneously, with the
    covariance Phi_{g.f}(O'') = Phi_f(O) (monomial identity for a generic matrix, and exact
    values on random rational O in SO(4))."""
    hits = [g for g, gm in M.G_EXPLICIT if all(M.g_facet(gm, a) == b for a, b in zip(rep, target))]
    ok = bool(hits)
    detail = ""
    if ok:
        g = hits[0]
        ok = M.generic_covariance(g)
        for _ in range(trials):
            O = M.random_rational_SO4(rng)
            O2 = M.g_matrix(g, O)
            if M.det4(O2) == -1:
                O2 = [[-row[0]] + list(row[1:]) for row in O2]
            ok &= M.det4(O) == 1 and M.det4(O2) == 1
            for a_, b_ in zip(rep, target):
                ok &= M.leibniz(F_VALUES[a_], O) == M.leibniz(F_VALUES[b_], O2)
        detail = "%d group elements, e.g. (rho=%s, rho'=%s, transpose=%d)" % (len(hits), g[0], g[1], g[2])
    ck.check("EXACT", "%s: a group element maps (%s) to (%s) simultaneously, and Leibniz sums transfer exactly"
             % (label, ", ".join("f%d" % x for x in rep), ", ".join("f%d" % x for x in target)), ok, detail)
    return hits


def main():
    numeric = "--no-numeric" not in sys.argv
    import sympy as sp
    import numpy
    import scipy
    import flint
    ck = M.Checker("Case T_39 = {2,3,5,7,9,12} (document 'T38'): independent exact replay")
    ck.info("python %s, sympy %s, python-flint %s (fmpz_mpoly arithmetic in the charts; integer determinants at "
            "test points), numpy %s (int64 in det_linear_forms; floats in NUMERIC checks), scipy %s (NUMERIC only)"
            % (sys.version.split()[0], sp.__version__, flint.__version__, numpy.__version__, scipy.__version__))

    # ------------------------------------------------------------------------------------
    ck.section("A. Provenance and foundations")
    M.check_document(ck, DOC, DOC_SHA256[DOC])
    M.check_foundations(ck)
    M.check_radial_necessity(ck)
    M.crosscheck_supplement_tables(ck)
    M.crosscheck_audit_library(ck, [38, 39, 25, 28, 29])

    # ------------------------------------------------------------------------------------
    ck.section("B. Contact set (p. 1)")
    ck.check("EXACT", "the document's set {2,3,5,7,9,12} is the paper's T_39 (Appendix E) as a set of "
             "lexicographic indices: group element = identity; the document's number is the paper's minus one",
             tuple(sorted(T)) == M.PAPER_T_V[39])
    r39, r38 = M.orbit_rep(T), M.orbit_rep(M.PAPER_T_V[38])
    stab = [g for g, gm in M.G_EXPLICIT if set(gm[s] for s in T) == set(T)]
    ck.check("EXACT", "orbit representative of T_39 under G is (0,1,2,6,11,23); T_38 lies in another orbit; the "
             "stabiliser of T_39 in G is trivial (the two non-minor generators are not exchanged by a symmetry)",
             r39 == (0, 1, 2, 6, 11, 23) and r38 != r39 and len(stab) == 1,
             "T_39 -> %s, T_38 -> %s" % (r39, r38))

    # ------------------------------------------------------------------------------------
    ck.section("C. Generators and the reduction to two cones with one non-minor generator (p. 1)")
    gens = M.generators(T)
    ck.check("EXACT", "the f_j vanishing on T are exactly f0, f1, f2, f5, f7, f12, f39, f40, f67, f70",
             gens == GENS, "computed %s" % gens)
    nm = [j for j in gens if not M.F_IS_MINOR[j]]
    ck.check("EXACT", "the non-minor generators are f40 = h_{4213,+} and f70 = h_{4231,+}",
             nm == [40, 70] and M.f_name(40) == "h_{4213,+}" and M.f_name(70) == "h_{4231,+}",
             ", ".join("f%d = %s" % (j, M.f_name(j)) for j in GENS))
    d, rays = M.extreme_rays_of_contact_cone(T)
    ck.check("EXACT", "direct enumeration of the extreme rays of N_T = {Pc >= 0, (Pc)_T = 0} (dimension %d; no "
             "use of the 82-facet list) gives exactly these ten vectors, so every c in N_T is a nonnegative "
             "combination of them" % d, d == 8 and rays == sorted(F_VALUES[j] for j in GENS), "%d rays" % len(rays))
    ck.check("EXACT", "the three identities of p. 1 hold on all 24 values (P is injective, so they hold in V_R): "
             "f40+f70 = f0+f5+f12+f39+f67, f1+f40 = f0+f7+f39, f7+f70 = f1+f5+f12+f67",
             all(M.combo_values([(j, 1) for j in l]) == M.combo_values([(j, 1) for j in r]) for l, r in IDENTITIES))
    for nm_, order in (("(1)", ORDER40), ("(2)", ORDER70)):
        ck.check("EXACT", "the eight generators of %s, %s, are linearly independent and their common zero set "
                 "is exactly T" % (nm_, order),
                 M.rank([M.F_COORD[j] for j in order]) == 8
                 and tuple(s for s in range(24) if all(F_VALUES[j][s] == 0 for j in order)) == T)
    w40 = functional_on_cone({40: 1, 70: -1, 1: -1})
    w70 = functional_on_cone({70: 1, 40: -1, 7: -1})
    ck.check("EXACT", "delta40 := (coefficient of f40) - (coefficients of f70, f1) and delta70 := (coefficient of "
             "f70) - (coefficients of f40, f7) are well defined linear functionals on N_T (independent of the "
             "representation): the vectors w with w . f_j = the stated values on all ten generators exist",
             all(sum(w40[s] * F_VALUES[j][s] for s in range(24)) == {40: 1, 70: -1, 1: -1}.get(j, 0) for j in GENS)
             and all(sum(w70[s] * F_VALUES[j][s] for s in range(24)) == {70: 1, 40: -1, 7: -1}.get(j, 0) for j in GENS)
             and all(w40[s] == 0 == w70[s] for s in T),
             "delta40 = %s, delta70 = %s (as combinations of gaps)"
             % ({s: str(v) for s, v in enumerate(w40) if v}, {s: str(v) for s, v in enumerate(w70) if v}))
    rng = random.Random(39001)
    ok = True
    kinds = {"minor": 0, "40": 0, "70": 0}
    for _ in range(4000):
        lam = {j: (Fraction(0) if rng.random() < 0.3 else Fraction(rng.randint(0, 12), rng.randint(1, 6)))
               for j in GENS}
        vec = M.combo_values(list(lam.items()))
        kind, out = reduce_to_one_nonminor(lam)
        kinds[kind] += 1
        ok &= M.combo_values(list(out.items())) == vec and all(v >= 0 for v in out.values())
        d40 = sum(w40[s] * vec[s] for s in range(24))
        d70 = sum(w70[s] * vec[s] for s in range(24))
        supp = {j for j, v in out.items() if v}
        if kind == "40":
            ok &= supp <= set(ORDER40) and out[40] == d40 > 0
        elif kind == "70":
            ok &= supp <= set(ORDER70) and out[70] == d70 > 0
        else:
            ok &= out[40] == 0 == out[70] and d40 <= 0 and d70 <= 0
    ck.check("EXACT", "the reduction of p. 1 (cancel the common part of f40, f70; then absorb f40 against f1, or "
             "f70 against f7) on 4000 random rational coefficient vectors with zeros (seed 39001): the vector is "
             "unchanged, all coefficients stay >= 0, and the result is a sum of minors (iff delta40, delta70 <= 0), "
             "or of the form (1) with d = delta40 > 0, or of the form (2) with d = delta70 > 0",
             ok and min(kinds.values()) > 0, "outcomes %s" % kinds)
    M.check_minor_leibniz(ck, [j for j in GENS if j not in (40, 70)])

    # ------------------------------------------------------------------------------------
    ck.section("D. The capacities used (pp. 2-4): orbit membership with roles, and algebra")
    rng = random.Random(39002)
    check_roles(ck, "T24-face core (Lemma 1 of the document for the paper's T_25)", CORE_REP, CORE_USE, rng)
    check_roles(ck, "the same core with the roles of f55 and f67 exchanged", CORE_REP, CORE_SWAP, rng)
    ck.check("EXACT", "the five core vectors f59, f53, f55, f67, f54 are generators of the cone of the paper's "
             "T_25 = {3,4,7,8,10,12} (the document's 'T24')",
             set(CORE_REP) <= set(M.generators(M.PAPER_T_V[25])))
    M.check_transport(ck, "old-H", (2, 5, 40), rng)
    M.check_transport(ck, "four-coefficient", (0, 67, 2, 40), rng)
    M.check_transport(ck, "four-coefficient", (0, 7, 12, 40), rng)
    M.check_transport(ck, "old-H", (12, 2, 70), rng)
    M.check_transport(ck, "old-H", (39, 2, 70), rng)
    M.check_transport(ck, "chain", (2, 0, 12, 70), rng)
    M.check_transport(ck, "chain", (2, 39, 67, 70), rng)
    M.check_transport(ck, "SOC", (2, 5, 70), rng)
    a, b, c, e, j, p, q = sp.symbols("a b c e j p q", positive=True)
    dl, r0, al, be = sp.symbols("delta r0 alpha beta", positive=True)
    four = lambda A_, B_, C_, D_: A_ * B_ * (C_ - D_) - D_ ** 2 * (A_ + B_ + C_ + D_)
    ck.check("EXACT", "(1), d = 1: four-coefficient with roles (f0, f67; f2; f40) = (a, b; e; 1) is "
             "e(ab-1) >= (a+1)(b+1), i.e. e >= E_ab for ab > 1; with roles (f0, f7; f12; f40) = (a, c; p; 1) it is "
             "p(ac-1) >= (a+1)(c+1), i.e. p >= P_ac for ac > 1; old-H t(c-d) >= d^2 with t-role f2 (coefficient e), "
             "c-role f5 (coefficient j) is j >= 1 + 1/e",
             sp.expand(four(a, b, e, 1) - (e * (a * b - 1) - (a + 1) * (b + 1))) == 0
             and sp.expand(four(a, c, p, 1) - (p * (a * c - 1) - (a + 1) * (c + 1))) == 0
             and sp.cancel(e * (j - 1) - 1 - e * (j - (1 + 1 / e))) == 0)
    par = lambda m_, w_: m_ * (1 - w_) ** 2 / (m_ - w_) ** 2
    ck.check("EXACT", "the parallel endpoint is symmetric: with a = alpha^2, b = beta^2, w = alpha beta, "
             "a(1-w)^2/(a-w)^2 = b(1-w)^2/(b-w)^2 = (1-w)^2/(alpha-beta)^2",
             sp.cancel(par(al ** 2, al * be) - (1 - al * be) ** 2 / (al - be) ** 2) == 0
             and sp.cancel(par(be ** 2, al * be) - (1 - al * be) ** 2 / (al - be) ** 2) == 0)
    Er = r0 * (a + 1) * (b + r0) / (a * b - r0)
    coh = lambda j_, e_, b_, a_, d_: d_ ** 2 / e_ - d_ * (a_ * b_ - d_ ** 2) / ((a_ + d_) * (b_ + d_))
    ck.check("EXACT", "allocation (eighth domain): the piece (delta a, c, p, delta) on (f0, f7; f12; f40) has "
             "four(delta a, c, p, delta) = delta [p(ac - delta) - delta (a+1)(c+delta)], zero on the document's "
             "p(delta); the piece (r0 a, b, e, r0) on (f0, f67; f2; f40) has four = r0 [e(ab - r0) - r0 (a+1)(b+r0)], "
             "i.e. e >= E_r0; the coherent bound d^2/e - d(ab-d^2)/((a+d)(b+d)) at (a, b, d) -> (r0 a, b, r0) equals "
             "r0^2/e - r0^2/E_r0; and at d = 1 it equals 1/e - 1/E_ab",
             sp.expand(four(dl * a, c, p, dl) - dl * (p * (a * c - dl) - dl * (a + 1) * (c + dl))) == 0
             and sp.expand(four(r0 * a, b, e, r0) - r0 * (e * (a * b - r0) - r0 * (a + 1) * (b + r0))) == 0
             and sp.cancel(coh(j, e, b, r0 * a, r0) - (r0 ** 2 / e - r0 ** 2 / Er)) == 0
             and sp.cancel(coh(j, e, b, a, 1) - (1 / e - (a * b - 1) / ((a + 1) * (b + 1)))) == 0)
    # H/chain resource lemma for (2)
    b1, b2, d1, d2, x, y, t1, t2, tt = sp.symbols("beta1 beta2 delta1 delta2 x y t1 t2 t", positive=True)
    e1, e2 = b1 - d1, b2 - d2
    bH, bC, cH, cC = b * e1 / b1, b * d1 / b1, c * e2 / b2, c * d2 / b2
    ck.check("EXACT", "resource lemma for (2), two old-H inequalities (f12, f2; f70), (f39, f2; f70) with "
             "t-coefficients t1, t2 and non-minor shares t t1/(t1+t2), t t2/(t1+t2): the required coefficient of "
             "f2 is t + t^2/(t1+t2)",
             sp.cancel(sum(dd_ + dd_ ** 2 / ti for dd_, ti in ((tt * t1 / (t1 + t2), t1), (tt * t2 / (t1 + t2), t2)))
                       - (tt + tt ** 2 / (t1 + t2))) == 0)
    ck.check("EXACT", "resource lemma, Cauchy step: eta1^2/x + eta2^2/y - (eta1+eta2)^2/(x+y) = "
             "(eta1 y - eta2 x)^2/(xy(x+y)) >= 0",
             sp.cancel(d1 ** 2 / x + d2 ** 2 / y - (d1 + d2) ** 2 / (x + y)
                       - (d1 * y - d2 * x) ** 2 / (x * y * (x + y))) == 0)
    total = (e1 + e2) + e1 ** 2 / bH + e2 ** 2 / cH + d1 ** 2 / p + d1 ** 2 / bC + d2 ** 2 / cC + d2 ** 2 / q
    rhsJ = 1 - d1 - d2 + b1 ** 2 / b + b2 ** 2 / c + d1 ** 2 / p + d2 ** 2 / q
    ck.check("EXACT", "resource lemma, bookkeeping: with b split as b eta1/beta1 (old-H) + b delta1/beta1 (chain "
             "with leaf p) and c split likewise (chain with leaf q), eta_i = beta_i - delta_i, the sum of the "
             "three requirements (eta1+eta2) + eta1^2/b_H + eta2^2/c_H + delta1^2/p + delta1^2/b_ch + "
             "delta2^2/c_ch + delta2^2/q equals 1 - delta1 - delta2 + beta1^2/b + beta2^2/c + delta1^2/p + "
             "delta2^2/q when beta1 + beta2 = 1; the shares of b, c and of the non-minor coefficient add up",
             sp.cancel((total - rhsJ).subs(b2, 1 - b1)) == 0 and sp.cancel(bH + bC - b) == 0
             and sp.cancel(cH + cC - c) == 0 and sp.expand((e1 + e2 + d1 + d2 - 1).subs(b2, 1 - b1)) == 0)

    # ------------------------------------------------------------------------------------
    ck.section("E. Radial necessity in the generators of (1) and (2) (pp. 1-2)")
    rd40 = M.radial_data(ck, ORDER40, NAMES40)
    R40, C40, Q40, G40 = rd40["R"], rd40["C0"], rd40["Q"], rd40["G"]
    ck.check("EXACT", "(1): C0 has 21 terms, every coefficient is positive, and the coefficient of d^3 is 1, so "
             "C0 >= d^3 > 0 and C0 + zQ = 0 forces Q != 0",
             len(C40) == 21 and all(v > 0 for v in C40.values())
             and C40.get(R40.key((0, 0, 0, 0, 0, 3, 0, 0, 0))) == 1)
    ck.check("EXACT", "(1): Q = d^2 + de + (d+e)(p+q) - ej (the document's formula)",
             Q40 == R40.from_sympy(sp.sympify("d**2 + d*e + (d+e)*(p+q) - e*j")))
    d9 = R40.key((0, 0, 0, 0, 0, 9, 0, 0, 0))
    if G40.get(d9, 0) < 0:
        G40 = p_scale(G40, -1)
    ck.check("EXACT", "(1): G_40 (sign fixed by coefficient +1 of d^9) has 1455 terms and degree 9 (the "
             "document's count)", len(G40) == 1455 and G40.get(d9) == 1 and all(sum(R40.exps(k)) == 9 for k in G40))
    rd70 = M.radial_data(ck, ORDER70, NAMES70, seed=20261126)
    R70, C70, Q70, G70 = rd70["R"], rd70["C0"], rd70["Q"], rd70["G"]
    T0s = "b*(d+r)*(d+q) + c*d*(d+p+r) + d*(d+p+r)*(d+q)"
    ck.check("EXACT", "(2): C0 = bc(a - U) - T0 and Q = (a - U) M - dU - as with M = b + c, U = d + r, "
             "T0 = bU(d+q) + cd(d+p+r) + d(d+p+r)(d+q) (the document's formulas)",
             C70 == R70.from_sympy(sp.sympify("b*c*(a - d - r) - (%s)" % T0s))
             and Q70 == R70.from_sympy(sp.sympify("(a - d - r)*(b + c) - d*(d + r) - a*s")))
    vb, vc, vd, va, vr = [R70.var(x) for x in "bcdar"]
    M_ = p_add(vb, vc)
    T0 = R70.from_sympy(sp.sympify(T0s))
    rest = p_add(p_mul(M_, T0), mul(vb, vc, vd, vd), -1)
    lhs = p_scale(p_mul(M_, C70), -1)
    rhs = p_add(p_add(mul(vb, vc, p_add(p_add(p_mul(M_, vd), mul(vd, vd)), p_mul(va, M_), -1)), mul(vr, M_, vb, vc)),
                rest)
    ck.check("EXACT", "(2): -M C0 = bc(Md + d^2 - aM) + rMbc + (M T0 - bc d^2), and M T0 - bc d^2 has %d terms, all "
             "coefficients nonnegative, among them b d^3 + c d^3 (at d = 1 this is the document's identity).  So "
             "if M > 0 and a < d + d^2/M then C0 < 0; and for b = c = 0, C0 = -d(d+p+r)(d+q) < 0.  In both cases "
             "C0 + zQ = 0 forces Q != 0" % len(rest),
             lhs == rhs and p_nonneg(rest) and rest.get(R70.key((0, 0, 0, 1, 0, 0, 3, 0, 0))) == 1
             and rest.get(R70.key((0, 0, 0, 0, 1, 0, 3, 0, 0))) == 1
             and R70.subs_zero(C70, ["b", "c"]) == R70.from_sympy(sp.sympify("-d*(d+p+r)*(d+q)")))
    d9 = R70.key((0, 0, 0, 0, 0, 0, 9, 0, 0))
    if G70.get(d9, 0) < 0:
        G70 = p_scale(G70, -1)
    ck.check("EXACT", "(2): G_70 (sign fixed by coefficient +1 of d^9) has 1113 terms and degree 9 (the "
             "document's count)", len(G70) == 1113 and G70.get(d9) == 1 and all(sum(R70.exps(k)) == 9 for k in G70))
    ck.info("the supplement contains no stored radial polynomials for this document; nothing to compare")

    b40 = branch40(ck, rd40, G40)
    b70 = branch70(ck, rd70, G70)
    if numeric:
        numeric_part(ck, rd40, G40, rd70, G70, w40, w70, b40, b70)
    ck.finish()


# ==========================================================================================
def branch40(ck, rd, G40):
    import sympy as sp
    ck.section("F. Branch (1), non-minor f40: the strict constant coefficient (p. 2)")
    R9 = rd["R"]
    Rs = Ring("q,p,j,c,e,b,a")
    G1 = p_rename(R9.subs_one(G40, "d"), R9, Rs)
    ks = Rs.split(G1, "j")
    ck.check("EXACT", "at d = 1: G_40 = k0 + k1 j + k2 j^2", set(ks) == {0, 1, 2} and len(G1) == 1455,
             "terms of k0, k1, k2: %d, %d, %d" % (len(ks[0]), len(ks[1]), len(ks[2])))
    k0, k1, k2 = ks[0], ks[1], ks[2]
    C1 = p_rename(R9.subs_one(rd["C0"], "d"), R9, Rs)
    G6 = p_divexact(k0, C1)
    ck.check("EXACT", "k0 = C G6 with C the cubic coefficient at d = 1 (exact division, re-multiplied); G6 has 128 "
             "terms and degree <= 1 in e and in p",
             G6 is not None and p_mul(G6, C1) == k0 and len(G6) == 128 and Rs.deg(G6, "e") == 1
             and Rs.deg(G6, "p") == 1 and C1.get(0) == 1)
    Rg = Ring("a,b,c,q")
    se, sp_ = SH * Rs.pos["e"], SH * Rs.pos["p"]
    gg = {}
    for k, cf in G6.items():
        ee, pp = (k >> se) & MASK, (k >> sp_) & MASK
        gg.setdefault((ee, pp), {})[k - (ee << se) - (pp << sp_)] = cf
    gg = {key: p_rename(v, Rs, Rg) for key, v in gg.items()}
    g00, gE, gP, gEP = gg[(0, 0)], gg[(1, 0)], gg[(0, 1)], gg[(1, 1)]
    ck.check("EXACT", "G6 = g00 + gE e + gP p + gEP ep; g00 has %d terms, all coefficients positive, constant term 1"
             % len(g00), all(v > 0 for v in g00.values()) and g00.get(0) == 1)
    ga, gb, gc = [Rg.var(n) for n in "abc"]
    RS = Ring("A,T,W,Q")
    ES = FL(RS)
    sA, sT, sW = [RS.var(n) for n in "ATW"]
    for ab in ("lt", "gt"):
        for ac in ("lt", "gt"):
            le = p_add(mul(ga, gb), ONE, -1) if ab == "gt" else {}
            lp = p_add(mul(ga, gc), ONE, -1) if ac == "gt" else {}
            de, dp = mul(P1(ga), P1(gb)), mul(P1(ga), P1(gc))
            polys = [g00, p_add(mul(gE, de), mul(g00, le)), p_add(mul(gP, dp), mul(g00, lp)),
                     p_add(p_add(p_add(mul(gEP, de, dp), mul(gE, de, lp)), mul(gP, dp, le)), mul(g00, le, lp))]
            bn, bd = (sT, mul(sA, P1(sT))) if ab == "lt" else (P1(sT), sA)
            cn, cd = (sW, mul(sA, P1(sW))) if ac == "lt" else (P1(sW), sA)
            counts, okc = [], True
            for poly in polys:
                out, _ = p_compose_fl(poly, Rg, ES, {"a": (sA, ONE), "b": (bn, bd), "c": (cn, cd)}, rename={"q": "Q"})
                counts.append(len(out))
                okc &= len(out) > 0 and p_nonneg(out)
            ck.check("EXACT", "chamber ab %s 1, ac %s 1: the four functions g00, gE + Le g00, gP + Lp g00, "
                     "gEP + Lp gE + Le gP + Le Lp g00, multiplied by (a+1)(b+1), (a+1)(c+1) as needed and cleared "
                     "of the chart denominators in the chart a = A, b = %s, c = %s, have only nonnegative "
                     "coefficients; numbers of terms %s (document: %s)"
                     % ("<" if ab == "lt" else ">", "<" if ac == "lt" else ">",
                        "T/(A(1+T))" if ab == "lt" else "(1+T)/A", "W/(A(1+W))" if ac == "lt" else "(1+W)/A",
                        tuple(counts), STRICT_COUNTS[(ab, ac)]), okc and tuple(counts) == STRICT_COUNTS[(ab, ac)])
    l00, lE, lP, lEP, lLe, lLp, le1, lp1 = sp.symbols("g00 gE gP gEP Le Lp e1 p1")
    eo, po = le1 / (1 + lLe * le1), lp1 / (1 + lLp * lp1)
    ck.check("EXACT", "with e = e'/(1 + Le e'), p = p'/(1 + Lp p'): (1 + Le e')(1 + Lp p') G6 = g00 + (gE + Le g00) e' "
             "+ (gP + Lp g00) p' + (gEP + Lp gE + Le gP + Le Lp g00) e'p'; and 1/(1 + Le e') = 1 - Le e",
             sp.cancel((1 + lLe * le1) * (1 + lLp * lp1) * (l00 + lE * eo + lP * po + lEP * eo * po)
                       - (l00 + (lE + lLe * l00) * le1 + (lP + lLp * l00) * lp1
                          + (lEP + lLp * lE + lLe * lP + lLe * lLp * l00) * le1 * lp1)) == 0
             and sp.cancel(1 / (1 + lLe * le1) - (1 - lLe * eo)) == 0)
    ck.info("PROOF-READ: on the hard chamber Le e < 1, Lp p < 1 put e' = e/(1 - Le e) >= 0, p' = p/(1 - Lp p) >= 0.  "
            "The four functions are >= 0 on the dense set {a > 0, ab != 1, ac != 1} of the closed orthant and "
            "continuous (Le, Lp are continuous), hence >= 0 everywhere; so G6 >= g00 (1 - Le e)(1 - Lp p) > 0 and "
            "k0 = C G6 >= (1 - Le e)(1 - Lp p) > 0 (C >= 1, g00 >= 1 at d = 1).")
    k1e0, k2e0 = Rs.subs_zero(k1, ["e"]), Rs.subs_zero(k2, ["e"])
    ck.check("EXACT", "e = 0: k2 vanishes identically and k1 has %d terms, all coefficients nonnegative (document: "
             "234); so G_40 = k0 + k1 j > 0 for e = 0, j >= 0 on the hard chamber" % len(k1e0),
             len(k2e0) == 0 and len(k1e0) == 234 and p_nonneg(k1e0))

    ck.section("G. Branch (1): the Bernstein data on the eight domains (pp. 3-4)")
    x0, x1, x2, tt, JJ = sp.symbols("k0 k1 k2 t J")
    ck.check("EXACT", "Bernstein identity: k0 + k1 (tJ) + k2 (tJ)^2 = (1-t)^2 k0 + t(1-t)(2 k0 + J k1) + "
             "t^2 (k0 + J k1 + J^2 k2)",
             sp.expand(x0 + x1 * tt * JJ + x2 * (tt * JJ) ** 2 - ((1 - tt) ** 2 * x0 + tt * (1 - tt) * (2 * x0 + JJ * x1)
                                                                + tt ** 2 * (x0 + JJ * x1 + JJ ** 2 * x2))) == 0)
    RO = Ring("A,T,W,E,P,Q")
    EO = FL(RO)
    oA, oT, oW, oE, oP, oQ = EO.gens
    o1 = EO.one
    md = {n_: max(Rs.deg(kk, n_) for kk in (k0, k1, k2)) for n_ in "abepqc"}
    ck.info("degrees of (k0, k1, k2): every chart clears the denominator of each variable to the maximal degree "
            "%s (the convention that reproduced the counts of the document for the paper's T_28; the present "
            "document has no public verifier)" % (md,))
    F_ = lambda n_, d_=None: RF(EO, n_, d_)
    fA, fT, fW, fE, fP = (F_(v) for v in (oA, oT, oW, oE, oP))
    zz = fT / (fT + 1)
    sturm_total = {"middle": 0, "endpoint": 0}
    sturm_polys = []

    def chart(abk, ack):
        if abk == "both_le1":           # a = (z^2 + A)/(1 + A), b = z^2/a, J = 1 + 1/e
            an, ad = oT * oT + oA * (oT + o1) ** 2, (oA + o1) * (oT + o1) ** 2
            bn, bd = oT * oT * (oA + o1), an
            a_doc = (zz * zz + fA) / (fA + 1)
            b_doc = zz * zz / a_doc
            dsc = "ab < 1, a, b <= 1"
        elif abk == "a_ge1":            # a = 1 + A, b = z^2/(1 + A), J = 1/e + a(1-z)^2/(a-z)^2
            an, ad = oA + o1, o1
            bn, bd = oT * oT, (oT + o1) ** 2 * (oA + o1)
            a_doc, b_doc = fA + 1, zz * zz / (fA + 1)
            dsc = "ab < 1, a >= 1"
        elif abk == "b_ge1":            # a = z^2/(1 + A), b = 1 + A, J = 1/e + b(1-z)^2/(b-z)^2
            bn, bd = oA + o1, o1
            an, ad = oT * oT, (oT + o1) ** 2 * (oA + o1)
            b_doc, a_doc = fA + 1, zz * zz / (fA + 1)
            dsc = "ab < 1, b >= 1"
        else:                           # a = A, b = (1 + T)/A
            an, ad = oA, o1
            bn, bd = oT + o1, oA
            a_doc, b_doc = fA, (fT + 1) / fA
            dsc = "ab > 1"
        fa, fb = F_(an, ad), F_(bn, bd)
        ok = fa == a_doc and fb == b_doc
        if abk == "gt":                 # e = E_ab E/(1+E), J = 1/(E_ab E) = 1/e - 1/E_ab
            en, ed = (oA + o1) * (oA + o1 + oT) * oE, oA * oT * (oE + o1)
            num, den = oA * oT, (oA + o1) * (oA + o1 + oT) * oE
            Eab = (fa + 1) * (fb + 1) / (fa * fb - 1)
            fe = F_(en, ed)
            ok = ok and fa * fb == fT + 1 and fe == Eab * fE / (fE + 1) and F_(num, den) == 1 / (Eab * fE) \
                and F_(num, den) == 1 / fe - 1 / Eab
        else:
            en, ed = oE, o1
            fe = fE
            ok = ok and fa * fb == zz * zz
            if abk == "both_le1":
                num, den = oE + o1, oE
                ok = ok and F_(num, den) == 1 + 1 / fe and (1 - fa) == (1 - zz * zz) / (1 + fA) \
                    and (1 - fb) == fA * (1 - zz * zz) / (zz * zz + fA)
            else:
                ff = (oA + o1) * (oT + o1) - oT
                num, den = ff * ff + oE * (oA + o1), oE * ff * ff
                big = fa if abk == "a_ge1" else fb
                ok = ok and F_(num, den) == 1 / fe + big * (1 - zz) * (1 - zz) / ((big - zz) * (big - zz))
        if ack == "lt":                 # c = W/(a(1+W)), p = P
            cn, cd = oW * ad, an * (oW + o1)
            pn, pd = oP, o1
            ok = ok and F_(cn, cd) == fW / (fa * (fW + 1))
            dsc += ", ac < 1"
        else:                           # c = (1+W)/a, p = P_ac P/(1+P), P_ac = (a+1)(c+1)/(ac-1)
            cn, cd = (oW + o1) * ad, an
            pn = (an + ad) * (an + (oW + o1) * ad) * oP
            pd = ad * an * oW * (oP + o1)
            fc = F_(cn, cd)
            Pac = (fa + 1) * (fc + 1) / (fa * fc - 1)
            ok = ok and fc == (fW + 1) / fa and F_(pn, pd) == Pac * fP / (fP + 1)
            dsc += ", ac > 1"
        return dict(a=(an, ad), b=(bn, bd), e=(en, ed), p=(pn, pd), c=(cn, cd)), num, den, dsc, ok

    def pos_coeffs(P_):
        return len(P_) > 0 and all(P_.coefficient(i) > 0 for i in range(len(P_)))

    def three_tables(key, subs, num, den, dsc):
        sk = [p_compose_fl(kk, Rs, EO, subs, degs=dict(md), rename={"q": "Q"}, as_fl=True)[0] for kk in (k0, k1, k2)]
        tabs = [("base", sk[0]), ("middle", sk[0] * den * EO.fmpz(2) + sk[1] * num),
                ("endpoint", sk[0] * den * den + sk[1] * num * den + sk[2] * num * num)]
        counts, notes, ok = [], [], True
        for kind, tab in tabs:
            n_, neg = M.fl_len_neg(tab)
            counts.append(n_)
            if neg == 0:
                ok &= n_ > 0
                continue
            poly = EO.from_fl(tab)
            ngr, bad = M.grouped_univariate(RO, poly, "E")
            good = 0
            for _, cf in bad:
                if M.uni_positive(cf):
                    good += 1
                    low = min(cf)
                    sturm_polys.append([cf.get(i + low, 0) for i in range(max(cf) - low + 1)])
            ok &= good == len(bad) and key == ("a_ge1", "gt") and len(bad) == STURM_GROUPS[kind]
            sturm_total[kind] += len(bad)
            notes.append("%s: %d negative coefficients in %d of %d groups by (A,T,W,P,Q); each group is E^low times "
                         "a polynomial in E with positive constant term and no root in (0, infinity) (own Sturm "
                         "sequences)" % (kind, neg, len(bad), ngr))
        posden = all(pos_coeffs(x_) for pr in subs.values() for x_ in pr) and pos_coeffs(num) and pos_coeffs(den)
        ck.check("EXACT", "chart %s: all chart numerators and denominators, num and den have only positive "
                 "coefficients (so the clearing factor is > 0 for A, T, W, E, P, Q > 0), and the three cleared "
                 "tables (base, middle, endpoint) are nonnegative there; numbers of terms %s (document: %s)%s"
                 % (dsc, tuple(counts), BERN_COUNTS[key],
                    "; " + "; ".join(notes) if notes else "; all coefficients nonnegative"),
                 ok and posden and tuple(counts) == BERN_COUNTS[key])
        return sk

    keep = {}
    for abk in ("both_le1", "a_ge1", "b_ge1", "gt"):
        for ack in ("lt", "gt"):
            if abk == "gt" and ack == "gt":
                continue
            subs, num, den, dsc, okf = chart(abk, ack)
            ck.check("EXACT", "chart %s: the substituted expressions for a, b, c, e, p and the endpoint J = num/den "
                     "are the document's formulas (identities of rational functions in A, T, W, E, P)" % dsc, okf)
            sk = three_tables((abk, ack), subs, num, den, dsc)
            if (abk, ack) == ("gt", "lt"):
                keep = dict(subs=subs, sk0=sk[0])
    # allocation chart: a = A, b = (1+T)/A, c = (1+W)/A, delta = U/(1+U) (U is the ring variable P here)
    oU = oP
    oD = oU + o1
    gap = (oT + o1) * oD - o1                               # Delta
    hh = (oT + o1) * oD + oA                                # H
    pn = oU * (oA + o1) * ((oW + o1) * oD + oA * oU)
    pd = oA * oD * ((oW + o1) * oD - oU)
    en, ed = (oA + o1) * hh * oE, oA * oD * gap * (oE + o1)
    num, den = oA * gap, oD * (oA + o1) * hh * oE
    fU = fP
    fa, fb, fc = fA, (fT + 1) / fA, (fW + 1) / fA
    fdl = fU / (fU + 1)
    frr = 1 - fdl
    p_doc = fdl * (fa + 1) * (fc + fdl) / (fa * fc - fdl)
    Er = frr * (fa + 1) * (fb + frr) / (fa * fb - frr)
    e_doc = Er * fE / (fE + 1)
    ck.check("EXACT", "allocation chart (ab >= 1, ac >= 1): with a = A, b = (1+T)/A, c = (1+W)/A, delta = U/(1+U), "
             "r0 = 1 - delta, the substituted p, e and J = num/den equal delta(a+1)(c+delta)/(ac-delta), "
             "E_r0 E/(1+E) and r0^2/e - r0^2/E_r0 (the document's p, e, J in terms of D, Delta, H); also "
             "(1+W)D - U = 1 + W + WU, so the denominator of p has positive coefficients",
             F_(pn, pd) == p_doc and F_(en, ed) == e_doc and F_(num, den) == frr * frr / e_doc - frr * frr / Er
             and (oW + o1) * oD - oU == oW * oU + oW + o1)
    three_tables(("gt", "alloc"), dict(a=(oA, o1), b=(oT + o1, oA), e=(en, ed), p=(pn, pd), c=(oW + o1, oA)),
                 num, den, "ab >= 1, ac >= 1 (allocation)")
    ck.check("EXACT", "sign certificates in E were needed only on the chart ab < 1, a >= 1, ac > 1: %d groups in "
             "the middle table and %d in the endpoint table (document: 308 and 491, 799 in all)"
             % (sturm_total["middle"], sturm_total["endpoint"]), sturm_total == STURM_GROUPS)
    xs = sp.Symbol("x")
    oks, degs = True, set()
    for f in sturm_polys:
        pol = sp.Poly(list(reversed(f)), xs)
        degs.add(pol.degree())
        oks = oks and f[0] > 0 and pol.count_roots(0, sp.oo) == 0
    ck.check("EXACT", "second opinion on the %d univariate polynomials (degrees %s): positive constant term, and "
             "sympy's exact real-root count on [0, infinity) is 0 for each" % (len(sturm_polys), sorted(degs)),
             oks and len(sturm_polys) == 799)
    # cross-check of the flint composition against the pure Python composition (seventh chart, k0)
    subs = {n_: (EO.from_fl(x_), EO.from_fl(y_)) for n_, (x_, y_) in keep["subs"].items()}
    ref, _ = p_compose(k0, Rs, RO, subs, degs=dict(md), rename={"q": "Q"})
    ck.check("EXACT", "the python-flint composition agrees with the pure Python composition p_compose of PART 2 "
             "on the base table of the seventh chart (%d terms)" % len(ref), ref == EO.from_fl(keep["sk0"]))
    ck.info("PROOF-READ: every table is a Bernstein datum times positive factors on the open chart, so the three "
            "data k0, 2k0 + J k1, k0 + J k1 + J^2 k2 are >= 0 on the image of the open chart and, by continuity "
            "of J on {e > 0}, on its closure; with k0 > 0 (section F) the Bernstein identity gives G_40 > 0 for "
            "0 <= j < J.  See T39.md for the boundaries.")
    return dict(Rs=Rs, G1=G1, ks=(k0, k1, k2), Rg=Rg, gs=(g00, gE, gP, gEP), C1=C1)


# ==========================================================================================
def kkt_charts(Rc):
    """The five rational charts of p. 5 for the positive quadruples (p, b, c, q) at the minimiser
    of the resource energy: substitutions and J = nn/dd."""
    L_, T_, U_, W_ = [Rc.var(x) for x in "LTUW"]
    out = {}
    for case in ("HH", "HC", "CH", "CC", "CCswap"):
        if case == "HH":
            pn, pd = mul(TWO, T_, U_), mul(P1(T_), P1(U_))
            bn, bd = T_, mul(P1(T_), L_)
            cn, cd = ONE, mul(P1(T_), L_)
            qn, qd = mul(TWO, W_), mul(P1(T_), P1(W_))
            nn = p_add(p_add(mul(TWO, P1(L_), P1(T_), P1(U_), P1(W_)), mul(T_, U_, P1(W_)), -1), mul(W_, P1(U_)), -1)
            dd = mul(TWO, P1(T_), P1(U_), P1(W_))
            states = "HH"
        elif case in ("HC", "CH"):
            pn, pd = mul(TWO, T_, U_), mul(P1(T_), P1(U_))
            bn, bd = T_, mul(P1(T_), L_)
            cn, cd = mul(TWO, P1(W_)), mul(P1(T_), p_add(mul(TWO, L_, P1(W_)), W_))
            qn, qd = mul(TWO, P1(W_)), P1(T_)
            nn = p_add(mul(p_add(mul(TWO, L_), ONE), P1(T_), P1(U_)), T_)
            dd = mul(TWO, P1(T_), P1(U_))
            states = "HC"
            if case == "CH":
                pn, pd, qn, qd, bn, bd, cn, cd = qn, qd, pn, pd, cn, cd, bn, bd
                states = "CH"
        else:
            V_ = p_add(U_, W_)
            pn, pd = mul(TWO, T_, P1(U_)), P1(T_)
            bn, bd = T_, mul(P1(T_), L_)
            cn, cd = mul(TWO, P1(U_), P1(V_)), mul(P1(T_), p_add(mul(TWO, L_, P1(U_), P1(V_)), W_))
            qn, qd = mul(TWO, P1(V_)), P1(T_)
            nn = p_add(mul(TWO, L_, P1(U_)), ONE)
            dd = mul(TWO, P1(U_))
            states = "CC"
            if case == "CCswap":
                pn, pd, qn, qd, bn, bd, cn, cd = qn, qd, pn, pd, cn, cd, bn, bd
        out[case] = (dict(p=(pn, pd), b=(bn, bd), c=(cn, cd), q=(qn, qd)), nn, dd, states,
                     case in ("CH", "CCswap"))
    return out


def branch70(ck, rd, G70):
    import sympy as sp
    ck.section("H. Branch (2), non-minor f70: the resource charts and the SOC peeling (pp. 4-5)")
    R9 = rd["R"]
    Rs = Ring("p,r,a,b,c,q,s")
    G1 = p_rename(R9.subs_one(G70, "d"), R9, Rs)
    Ls = Rs.split(G1, "s")
    ck.check("EXACT", "at d = 1: G_70 = L0 + s L1 + s^2 L2 with each L_i of degree 2 in a",
             set(Ls) == {0, 1, 2} and len(G1) == 1113 and all(Rs.deg(Ls[i], "a") == 2 for i in range(3)),
             "terms of L0, L1, L2: %d, %d, %d" % (len(Ls[0]), len(Ls[1]), len(Ls[2])))
    L0, L1, L2 = Ls[0], Ls[1], Ls[2]
    va = Rs.var("a")
    B1 = p_add(L1, p_scale(p_mul(va, L0), 2))
    B2 = p_add(p_add(L2, p_mul(va, L1)), mul(va, va, L0))
    L00 = Rs.split(L0, "a")[0]
    ck.check("EXACT", "the a^0 coefficient of L0 has 231 terms, all coefficients positive, constant term 1; "
             "B1 := L1 + 2a L0 has degree 3 in a, B2 := L2 + a L1 + a^2 L0 has degree 4 in a and vanishes at a = 0",
             len(L00) == 231 and all(v > 0 for v in L00.values()) and L00.get(0) == 1
             and Rs.deg(B1, "a") == 3 and Rs.deg(B2, "a") == 4 and not Rs.split(B2, "a").get(0))
    Rt = Ring("p,r,a,b,c,q,S")
    vS, vat = Rt.var("S"), Rt.var("a")
    lhs, _ = p_compose(G1, Rs, Rt, {"s": (vS, p_add(ONE, p_mul(vat, vS)))})
    ren = lambda X: p_rename(X, Rs, Rt)
    rhs = p_add(p_add(ren(L0), p_mul(vS, ren(B1))), mul(vS, vS, ren(B2)))
    ck.check("EXACT", "(1 + aS)^2 G_70(a, s = S/(1 + aS)) = L0 + S B1 + S^2 B2 (polynomial identity); for as < 1 "
             "and S = s/(1 - as) >= 0 one has 1 + aS = 1/(1 - as), so G_70 = (1 - as)^2 (L0 + S B1 + S^2 B2)",
             lhs == rhs)
    Rc = Ring("L,T,U,W,r")
    E = FL(Rc)
    sL, sT, sU, sW = sp.symbols("L T U W")
    theta = sT / (1 + sT)
    charts = kkt_charts(Rc)

    def nonneg_rat(expr):
        n_, d_ = sp.fraction(sp.cancel(expr))
        cn, cd = sp.Poly(n_, sL, sT, sU, sW).coeffs(), sp.Poly(d_, sL, sT, sU, sW).coeffs()
        return (all(x >= 0 for x in cn) and all(x > 0 for x in cd)) or \
            (all(x <= 0 for x in cn) and all(x < 0 for x in cd))

    doc_forms = {
        "HH": dict(b=theta / sL, c=(1 - theta) / sL, p=2 * theta * sU / (1 + sU), q=2 * (1 - theta) * sW / (1 + sW)),
        "HC": dict(b=theta / sL, p=2 * theta * sU / (1 + sU), q=2 * (1 - theta) * (1 + sW),
                   c=2 * (1 - theta) * (1 + sW) / (2 * sL * (1 + sW) + sW)),
        "CC": dict(p=2 * theta * (1 + sU), q=2 * (1 - theta) * (1 + sU + sW), b=theta / sL,
                   c=2 * (1 - theta) * (1 + sU) * (1 + sU + sW) / (2 * sL * (1 + sU) * (1 + sU + sW) + sW)),
    }
    doc_J = {"HH": lambda f: 1 + sL - (f["p"] + f["q"]) / 4, "HC": lambda f: sL + sp.Rational(1, 2) + theta / (2 * (1 + sU)),
             "CC": lambda f: sL + 1 / (2 * (1 + sU))}
    swap = {"p": "q", "q": "p", "b": "c", "c": "b"}
    for case, (subs, nn, dd, states, swapped) in charts.items():
        val = {x: Rc.to_sympy(subs[x][0]) / Rc.to_sympy(subs[x][1]) for x in "pbcq"}
        base = {"CH": "HC", "CCswap": "CC"}.get(case, case)
        form = doc_forms[base]
        okform = all(sp.cancel(val[swap[x] if swapped else x] - form[x]) == 0 for x in "pbcq")
        Jc = Rc.to_sympy(nn) / Rc.to_sympy(dd)
        okform &= sp.cancel(Jc - doc_J[base](form)) == 0
        pp, bb, cc, qq = val["p"], val["b"], val["c"], val["q"]
        beta = [theta, 1 - theta]
        if swapped:
            beta.reverse()
        delta = [pp / 2 if states[0] == "H" else beta[0], qq / 2 if states[1] == "H" else beta[1]]
        energy = 1 - delta[0] - delta[1] + beta[0] ** 2 / bb + beta[1] ** 2 / cc + delta[0] ** 2 / pp + delta[1] ** 2 / qq
        feas = all(nonneg_rat(x) for x in (delta[0], delta[1], beta[0] - delta[0], beta[1] - delta[1], Jc))
        der = [2 * beta[0] / bb, 2 * beta[1] / cc]
        for i, res in enumerate((pp, qq)):
            if states[i] == "C":
                der[i] += 2 * beta[i] / res - 1
        leaf = True
        for i, res in enumerate((pp, qq)):
            leaf &= nonneg_rat(beta[i] - res / 2) if states[i] == "H" else nonneg_rat(res / 2 - beta[i])
        ck.check("EXACT", "chart %s: the substitutions and J are the document's formulas%s; with theta = T/(1+T), "
                 "(beta1, beta2) = %s and delta_i = %s, the point (beta, delta) is feasible (0 <= delta_i <= beta_i "
                 "as rational functions with nonnegative coefficients), the resource energy equals the chart's "
                 "J = n/d0 > 0, the leaf conditions hold, and d/dbeta1 of the reduced function vanishes"
                 % (case, " with (p, b) and (q, c) exchanged" if swapped else "",
                    "(1-theta, theta)" if swapped else "(theta, 1-theta)",
                    "(" + ", ".join("resource_i/2" if st == "H" else "beta_i" for st in states) + ")"),
                 okform and feas and leaf and sp.cancel(energy - Jc) == 0 and sp.cancel(der[0] - der[1]) == 0)
    Rk = Ring("p,r,b,c,q")
    keep = None
    for case, (subs, nn, dd, states, swapped) in charts.items():
        nnf, ddf = E.to_fl(nn), E.to_fl(dd)
        rows, ok = [], True
        for X, m in ((L0, 2), (B1, 3), (B2, 4)):
            Xk = {k: p_rename(v, Rs, Rk) for k, v in Rs.split(X, "a").items()}
            degs = {}
            for v in Xk.values():
                for nm in subs:
                    degs[nm] = max(degs.get(nm, 0), Rk.deg(v, nm))
            comp = {k: p_compose_fl(v, Rk, E, subs, degs=degs, as_fl=True)[0] for k, v in Xk.items()}
            if case == "CH" and m == 2:
                keep = (Xk[0], dict(degs), comp[0], subs)
            tabs = []
            for i in range(m + 1):
                tot = E.zero
                for k, Xt in comp.items():
                    if k <= i:
                        tot = tot + Xt * nnf ** k * ddf ** (m - k) * E.fmpz(comb(m - k, i - k))
                tabs.append(tot)
            st = [M.fl_len_neg(t_) for t_ in tabs]
            rows.append([s_[0] for s_ in st])
            ok &= all(s_[1] == 0 for s_ in st)
        ck.check("EXACT", "chart %s: with a = J Y/(1+Y), J = n/d0, and one common denominator for every power of a "
                 "within each polynomial, the coefficients of Y^0..Y^m in (d0 (1+Y))^m X, X = L0 (m = 2), B1 (m = 3), "
                 "B2 (m = 4), are polynomials in (L, T, U, W, r) with nonnegative coefficients; numbers of terms "
                 "%s (document: %s)" % (case, rows, list(KKT_COUNTS[case])),
                 ok and rows == list(KKT_COUNTS[case]))
    Xk0, degs, compfl, subs = keep
    ref, _ = p_compose(Xk0, Rk, Rc, subs, degs=degs)
    ck.check("EXACT", "the python-flint composition agrees with the pure Python composition p_compose of PART 2 "
             "on the a^0 coefficient of L0 in chart CH (%d terms)" % len(ref), ref == E.from_fl(compfl))
    ck.info("PROOF-READ: for positive (p, b, c, q) in a chart and 0 <= a < J (Y = a/(J - a) >= 0, 1/(1+Y) = 1 - a/J): "
            "(1+Y)^2 d0^2 Delta L0 >= its Y^0 coefficient = d0^2 Delta L0(a = 0), so L0 >= L0(a = 0)(1 - a/J)^2 >= "
            "(1 - a/J)^2; B1, B2 >= 0; hence for as < 1: G_70 = (1 - as)^2 (L0 + S B1 + S^2 B2) >= "
            "(1 - as)^2 (1 - a/J)^2 > 0.  r >= 0 is a chart variable.  So a point of X of the form (2) has as >= 1 "
            "(SOC) or a >= J (resource lemma).")
    return dict(Rs=Rs, G1=G1, L00=L00, LB=(L0, B1, B2))


# ==========================================================================================
def energy_min(p, b, c, q, n=4000):
    """NUMERIC: min over beta1 in [0,1] of 1 - delta1 - delta2 + beta1^2/b + beta2^2/c + delta1^2/p +
    delta2^2/q with the optimal delta_i = min(resource/2, beta_i); zero coefficients follow the
    document's conventions (0/0 = 0, x^2/0 = +infinity).  Grid + golden-section refinement of a
    convex function.  Returns (J, beta1)."""
    def h(res, beta):
        if res <= 0:
            return 0.0
        dlt = min(res / 2, beta)
        return dlt * dlt / res - dlt

    def sq(beta, res):
        if beta == 0:
            return 0.0
        return float("inf") if res <= 0 else beta * beta / res

    def phi(x):
        return 1 + sq(x, b) + sq(1 - x, c) + h(p, x) + h(q, 1 - x)
    xs = [i / n for i in range(n + 1)]
    vals = [phi(x) for x in xs]
    i = min(range(n + 1), key=lambda k: vals[k])
    lo, hi = xs[max(i - 1, 0)], xs[min(i + 1, n)]
    for _ in range(80):
        m1, m2 = lo + (hi - lo) * 0.381966, lo + (hi - lo) * 0.618034
        if phi(m1) <= phi(m2):
            hi = m2
        else:
            lo = m1
    x = (lo + hi) / 2
    return min([(phi(x), x), (vals[i], xs[i]), (phi(0.0), 0.0), (phi(1.0), 1.0)])


def doc_endpoint(a, b, c, e, p):
    """NUMERIC (floats): the document's endpoint J for the form (1) with d = 1.  Returns 0.0 when
    one of the two four-coefficient tests applies (or e >= E_r0 on the allocation domain) and inf
    when e = 0."""
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
    A2, B2, C2 = a + 1.0, (a + 1.0) * c + p, -p * a * c
    delta = (-B2 + math.sqrt(B2 * B2 - 4.0 * A2 * C2)) / (2.0 * A2)
    r = 1.0 - delta
    if a * b - r <= 0.0:
        return r * r / e
    Er = r * (a + 1.0) * (b + r) / (a * b - r)
    if e >= Er:
        return 0.0
    return r * r / e - r * r / Er


def exact_sampling40(ck, R9, G40):
    """Exact rational sampling of the chambers of branch (1) in the ORIGINAL coefficients (no
    chart), including the closed boundaries that the open charts do not reach."""
    parts = R9.split(G40, "j")                             # K_i homogeneous of degree 9 - i
    Kh = [parts.get(i, {}) for i in range(3)]
    rng = random.Random(39005)
    names = R9.names

    def rnd(lo=0, hi=6, den=12, zero=0.15):
        if rng.random() < zero:
            return Fraction(0)
        return Fraction(rng.randint(lo * den, hi * den), den)

    def pos(hi=6, den=12):
        return Fraction(rng.randint(1, hi * den), den)
    power = [0]

    def test(av, bv, cv, ev, pv, qv, J):
        pt = dict(q=qv, p=pv, j=Fraction(0), c=cv, e=ev, d=Fraction(1), b=bv, a=av, z=Fraction(0))
        N = 1
        for x_ in list(pt.values()) + [J]:
            N = N * x_.denominator // math.gcd(N, x_.denominator)
        ipt = {n_: int(pt[n_] * N) for n_ in names}
        k = [Fraction(R9.evaluate(Kh[i], ipt), N ** (9 - i)) for i in range(3)]
        ok = k[0] > 0 and 2 * k[0] + J * k[1] >= 0 and k[0] + J * k[1] + J * J * k[2] >= 0
        for _ in range(5):
            jv = J * Fraction(rng.randint(0, 999), 1000)
            ok = ok and k[0] + k[1] * jv + k[2] * jv * jv > 0
        if k[0] + 4 * J * k[1] + 16 * J * J * k[2] < 0:
            power[0] += 1
        return ok
    tallies = {}

    def run(name, n, gen):
        ok, cnt = True, 0
        while cnt < n:
            out = gen()
            if out is None:
                continue
            cnt += 1
            ok = ok and test(*out)
        tallies[name] = cnt
        return ok

    def cp_for(av):
        """(c, p) in the hard chamber for a given a: ac < 1, ac = 1, or ac > 1 with p < P_ac."""
        mode = rng.randint(0, 5)
        if av == 0:
            return rnd(), rnd()
        if mode <= 1:
            return Fraction(rng.randint(0, 99), 100) / av, rnd()
        if mode == 2:
            return 1 / av, rnd(hi=30)
        cv = (1 + pos()) / av
        Pac = (av + 1) * (cv + 1) / (av * cv - 1)
        return cv, Pac * Fraction(rng.randint(0, 999), 1000)

    def g_both_le1():
        av, bv = Fraction(rng.randint(0, 12), 12), Fraction(rng.randint(0, 12), 12)
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
    res = [run("ab<1, a,b<=1", 700, g_both_le1), run("ab<1, a>=1", 700, g_one_ge1(False)),
           run("ab<1, b>=1", 700, g_one_ge1(True)), run("ab>=1, ac<1", 700, g_gt_lt),
           run("ab>=1, ac>=1 (allocation)", 700, g_alloc)]
    ck.check("NUMERIC", "branch (1), exact evaluation at %d random rational points of the five chamber types in "
             "the original coefficients (seed 39005; zeros, a = 1, b = 1, ab = 1, ac = 1, delta = 0 and the point "
             "ab = 1, p = 0 included): k0 > 0, 2k0 + J k1 >= 0, k0 + J k1 + J^2 k2 >= 0 with the document's endpoint "
             "J, and G_40 > 0 at five random j in [0, J)" % sum(tallies.values()), all(res), str(tallies))
    ck.check("NUMERIC", "the test has power: with the endpoint 4J in place of J the third Bernstein datum is "
             "negative at %d of these points (the radial equation does have roots j > J)" % power[0], power[0] > 0)
    ok = True
    for _ in range(800):
        av, bv = rnd(), rnd()
        cv, pv = cp_for(av)
        pt = dict(q=rnd(), p=pv, j=rnd(hi=40), c=cv, e=Fraction(0), d=Fraction(1), b=bv, a=av, z=Fraction(0))
        N = 1
        for x_ in pt.values():
            N = N * x_.denominator // math.gcd(N, x_.denominator)
        ok = ok and R9.evaluate(G40, {n_: int(pt[n_] * N) for n_ in names}) > 0
    ck.check("NUMERIC", "branch (1), exact evaluation: G_40 > 0 at 800 random rational points with e = 0, j >= 0 "
             "and (a, c, p) in the hard chamber", ok)


def numeric_part(ck, rd40, G40, rd70, G70, w40, w70, b40, b70):
    import numpy as np
    ck.section("I. NUMERIC falsification attempts (evidence only)")
    f = lambda idx: np.array(F_VALUES[idx], float)
    rng = random.Random(39003)
    ex = lambda: math.exp(rng.gauss(0.0, 1.0))
    worst = {}

    def upd(name, vec, seed):
        worst[name] = min(worst.get(name, 1e9), M.numeric_min_phi(vec, seed, 6000, 8))
    for trial in range(10):
        a, b, c, e, p, q = [ex() for _ in range(6)]
        # branch (1): j f5, e f2, b f67, a f0, p f12, c f7, d f40
        upd("old-H (f2, f5; f40), j = 1 + 1/e", e * f(2) + (1 + 1 / e) * f(5) + f(40), 100 + trial)
        z = rng.uniform(0.05, 0.95)
        big = 1 + ex()
        small = z * z / big
        Jp = 1 / e + big * (1 - z) ** 2 / (big - z) ** 2
        upd("parallel, a >= 1 > ab", big * f(0) + small * f(67) + e * f(2) + Jp * f(5) + f(40), 200 + trial)
        upd("parallel, b >= 1 > ab", small * f(0) + big * f(67) + e * f(2) + Jp * f(5) + f(40), 300 + trial)
        bb_ = (1 + ex()) / a
        Eab = (a + 1) * (bb_ + 1) / (a * bb_ - 1)
        e2 = Eab * rng.uniform(0.05, 0.95)
        upd("coherent, ab > 1, j = 1/e - 1/E_ab", a * f(0) + bb_ * f(67) + e2 * f(2) + (1 / e2 - 1 / Eab) * f(5) + f(40),
            400 + trial)
        upd("double SOC, ab = 1, j = 1/e", a * f(0) + (1 / a) * f(67) + e * f(2) + (1 / e) * f(5) + f(40), 500 + trial)
        upd("four-coefficient (f0, f67; f2; f40), e = E_ab", a * f(0) + bb_ * f(67) + Eab * f(2) + f(40), 600 + trial)
        cc_ = (1 + ex()) / a
        Pac = (a + 1) * (cc_ + 1) / (a * cc_ - 1)
        upd("four-coefficient (f0, f7; f12; f40), p = P_ac", a * f(0) + cc_ * f(7) + Pac * f(12) + f(40), 700 + trial)
        dl = rng.uniform(0.05, 0.95)
        pv = dl * (a + 1) * (cc_ + dl) / (a * cc_ - dl)
        r0 = 1 - dl
        Er = r0 * (a + 1) * (bb_ + r0) / (a * bb_ - r0)
        e3 = Er * rng.uniform(0.05, 0.95)
        upd("allocation (eighth domain), j = r0^2/e - r0^2/E_r0",
            a * f(0) + bb_ * f(67) + cc_ * f(7) + pv * f(12) + e3 * f(2) + (r0 * r0 / e3 - r0 * r0 / Er) * f(5) + f(40),
            800 + trial)
        # branch (2): p f0, a f2, b f12, c f39, q f67, s f5, d f70
        upd("old-H (f12, f2; f70), a = 1 + 1/b", b * f(12) + (1 + 1 / b) * f(2) + f(70), 900 + trial)
        upd("old-H (f39, f2; f70), a = 1 + 1/c", c * f(39) + (1 + 1 / c) * f(2) + f(70), 1000 + trial)
        upd("chain (f2; f0, f12; f70), a = 1/p + 1/b", (1 / p + 1 / b) * f(2) + p * f(0) + b * f(12) + f(70), 1100 + trial)
        upd("chain (f2; f39, f67; f70), a = 1/c + 1/q", (1 / c + 1 / q) * f(2) + c * f(39) + q * f(67) + f(70),
            1200 + trial)
        upd("SOC (f2, f5; f70), as = 1", a * f(2) + (1 / a) * f(5) + f(70), 1300 + trial)
        beta1 = rng.random()
        dl1, dl2 = rng.random() * beta1, rng.random() * (1 - beta1)
        aa = 1 - dl1 - dl2 + beta1 ** 2 / b + (1 - beta1) ** 2 / c + dl1 ** 2 / p + dl2 ** 2 / q
        upd("resource lemma at a random feasible (beta, delta)",
            p * f(0) + b * f(12) + c * f(39) + q * f(67) + aa * f(2) + f(70), 1400 + trial)
        Jm, _ = energy_min(p, b, c, q)
        upd("resource lemma at the minimum J(p,b,c,q)",
            p * f(0) + b * f(12) + c * f(39) + q * f(67) + Jm * f(2) + f(70), 1500 + trial)
    for nm, wv in worst.items():
        ck.check("NUMERIC", "%s (10 random coefficient draws, seed 39003; 6000 random rotations + 8 BFGS descents "
                 "each): min Phi found >= -1e-9" % nm, wv > -1e-9, "min = %.3e" % wv)
    ck.info("controls (the search does detect negative values): min Phi found for f40 alone %.4f, for f70 alone "
            "%.4f, for f2 + f5 + f40 (j = 1 < 1 + 1/e = 2) %.4f, for 2 f0 + 0.1 f67 + f2 + f40 (j = 0) %.4f, for "
            "f12 + 1.2 f2 + f70 (a = 1.2 < 1 + 1/b = 2) %.4f, for f0 + f12 + 0.2 f2 + f70 (a = 0.2 < 1/p + 1/b = 2) %.4f"
            % (M.numeric_min_phi(f(40), 96, 6000, 8), M.numeric_min_phi(f(70), 97, 6000, 8),
               M.numeric_min_phi(f(2) + f(5) + f(40), 98, 6000, 8),
               M.numeric_min_phi(2 * f(0) + 0.1 * f(67) + f(2) + f(40), 99, 6000, 8),
               M.numeric_min_phi(f(12) + 1.2 * f(2) + f(70), 99, 6000, 8),
               M.numeric_min_phi(f(0) + f(12) + 0.2 * f(2) + f(70), 99, 6000, 8)))

    # ---- branch (1): the four coefficient functions and the strict bound for k0, original coefficients
    Rg, (g00, gE, gP, gEP), Rs1, C1, k0 = b40["Rg"], b40["gs"], b40["Rs"], b40["C1"], b40["ks"][0]
    rngs = random.Random(39004)
    oks = okk = True
    nedge = 0
    for _ in range(2500):
        rv = lambda: Fraction(rngs.randint(0, 60), 12) if rngs.random() > 0.12 else Fraction(0)
        av, bv, cv, qv = rv(), rv(), rv(), rv()
        if av > 0 and rngs.random() < 0.15:
            bv = 1 / av                                   # ab = 1
            nedge += 1
        if av > 0 and rngs.random() < 0.15:
            cv = 1 / av                                   # ac = 1
            nedge += 1
        pt = dict(a=av, b=bv, c=cv, q=qv)
        v00, vE, vP, vEP = (Rg.evaluate(x_, pt) for x_ in (g00, gE, gP, gEP))
        Le = max((av * bv - 1) / ((av + 1) * (bv + 1)), Fraction(0))
        Lp = max((av * cv - 1) / ((av + 1) * (cv + 1)), Fraction(0))
        oks = oks and v00 >= 1 and vE + Le * v00 >= 0 and vP + Lp * v00 >= 0 \
            and vEP + Lp * vE + Le * vP + Le * Lp * v00 >= 0
        ev = Fraction(rngs.randint(0, 999), 1000) / Le if Le > 0 else rv()
        pv = Fraction(rngs.randint(0, 999), 1000) / Lp if Lp > 0 else rv()
        full = dict(q=qv, p=pv, j=Fraction(0), c=cv, e=ev, b=bv, a=av)
        okk = okk and Rs1.evaluate(k0, full) >= Rs1.evaluate(C1, full) * v00 * (1 - Le * ev) * (1 - Lp * pv) > 0
    ck.check("NUMERIC", "branch (1), exact evaluation at 2500 random rational points of the closed orthant (seed "
             "39004; zeros, and %d points with ab = 1 or ac = 1): the four coefficient functions are >= 0 and "
             "g00 >= 1; and with (e, p) random in the hard chamber, k0 >= C g00 (1 - Le e)(1 - Lp p) > 0" % nedge,
             oks and okk)

    # ---- branch (1): exact sampling in the original coefficients
    exact_sampling40(ck, rd40["R"], G40)

    # ---- branch (2): floats, J from a direct one-dimensional minimisation (no charts)
    Rs2, G2, L00 = b70["Rs"], b70["G1"], b70["L00"]
    L0_, B1_, B2_ = b70["LB"]
    absB1 = {k: abs(v) for k, v in B1_.items()}
    absB2 = {k: abs(v) for k, v in B2_.items()}
    rng = random.Random(39006)
    n8 = bad8 = ninf = badL = 0
    minratio = minL = float("inf")
    minB = 0.0
    for _ in range(3000):
        def rv():
            return 0.0 if rng.random() < 0.12 else math.exp(rng.gauss(0.0, 1.5))
        p, b, c, q, r, s = rv(), rv(), rv(), rv(), rv(), rv()
        Jm, _ = energy_min(p, b, c, q)
        if Jm == float("inf"):
            ninf += 1
            av = math.exp(rng.gauss(0.0, 2.0))
            fa_ = 0.0
        else:
            fa_ = rng.random() * 0.98
            av = fa_ * Jm
        fs_ = rng.random() * 0.98
        s = fs_ / av if av > 0 else s                       # as < 1
        ptf = dict(p=p, r=r, a=av, b=b, c=c, q=q, s=s)
        val = Rs2.evaluate(G2, ptf)
        low = (1 - fs_ if av > 0 else 1.0) ** 2 * (1 - fa_) ** 2
        n8 += 1
        if not val >= low * (1 - 1e-9):
            bad8 += 1
        minratio = min(minratio, val / low)
        l0 = Rs2.evaluate(L0_, ptf)
        rb1 = Rs2.evaluate(B1_, ptf) / max(Rs2.evaluate(absB1, ptf), 1e-300)
        rb2 = Rs2.evaluate(B2_, ptf) / max(Rs2.evaluate(absB2, ptf), 1e-300)
        minL = min(minL, l0 / (1 - fa_) ** 2)
        minB = min(minB, rb1, rb2)
        if not (l0 >= (1 - fa_) ** 2 * (1 - 1e-9) and rb1 >= -1e-12 and rb2 >= -1e-12):
            badL += 1
    ck.check("NUMERIC", "branch (2) in floating point at %d random points (seed 39006; p, b, c, q, r = 0 with "
             "probability 0.12 each, else exp(N(0,1.5)); J computed by direct minimisation of the resource energy, "
             "not through the charts; %d points with J = infinity, a arbitrary; a < J and as < 1): "
             "G_70 >= (1 - as)^2 (1 - a/J)^2 (relative tolerance 1e-9)" % (n8, ninf), bad8 == 0,
             "min ratio %.6f" % minratio)
    ck.check("NUMERIC", "branch (2), the three intermediate inequalities at the same %d points: L0 >= (1 - a/J)^2 "
             "(relative tolerance 1e-9), B1 >= 0 and B2 >= 0 (tolerance 1e-12 relative to the sum of the absolute "
             "values of the terms)" % n8, badL == 0,
             "smallest L0/(1 - a/J)^2 = %.6f; smallest relative value of B1, B2 = %.1e" % (minL, minB))
    # chart coverage
    rng = random.Random(39007)
    Rc = Ring("L,T,U,W,r")
    charts = kkt_charts(Rc)
    cover = {}
    worst_err, okcov = 0.0, True
    for _ in range(400):
        p, b, c, q = [math.exp(rng.gauss(0.0, 1.3)) for _ in range(4)]
        Jm, th = energy_min(p, b, c, q, n=20000)
        st = ("H" if p / 2 < th else "C") + ("H" if q / 2 < 1 - th else "C")
        if st == "HH":
            case, L, T_, U, W = "HH", th / b, th / (1 - th), p / (2 * th - p), q / (2 * (1 - th) - q)
        elif st == "HC":
            case, L, T_, U, W = "HC", th / b, th / (1 - th), p / (2 * th - p), q / (2 * (1 - th)) - 1
        elif st == "CH":
            case, L, T_, U, W = "CH", (1 - th) / c, (1 - th) / th, q / (2 * (1 - th) - q), p / (2 * th) - 1
        else:
            up, uq = p / (2 * th) - 1, q / (2 * (1 - th)) - 1
            if uq >= up:
                case, L, T_, U, W = "CC", th / b, th / (1 - th), up, uq - up
            else:
                case, L, T_, U, W = "CCswap", (1 - th) / c, (1 - th) / th, uq, up - uq
        cover[case] = cover.get(case, 0) + 1
        subs, nn, dd, states, swapped = charts[case]
        pt = dict(L=L, T=T_, U=U, W=W, r=0.0)
        img = [Rc.evaluate(subs[x][0], pt) / Rc.evaluate(subs[x][1], pt) for x in "pbcq"]
        Jc = Rc.evaluate(nn, pt) / Rc.evaluate(dd, pt)
        err = max(abs(u_ - v_) / v_ for u_, v_ in zip(img + [Jc], [p, b, c, q, Jm]))
        worst_err = max(worst_err, err)
        okcov &= min(L, T_, U, W) >= -1e-6 and err < 1e-4
    ck.check("NUMERIC", "coverage of the five charts: for 400 random positive (p,b,c,q) (seed 39007) the minimiser "
             "theta of the resource energy (computed directly) determines a chart and parameters L, T > 0, "
             "U, W >= 0 whose image is (p,b,c,q) and whose J equals the minimum (relative error < 1e-4)",
             okcov and len(cover) == 5, "charts hit: %s; worst relative error %.1e" % (cover, worst_err))

    # ---- genuine points of X
    eq = M.numeric_equality_points(T, 39009, 80)
    rd40n, rd70n = dict(rd40, G=G40), dict(rd70, G=G70)
    for nm, rd, order in (("(1)", rd40n, ORDER40), ("(2)", rd70n, ORDER70)):
        fit, cub, hes, qmin, gmax = M.numeric_radial_test(rd, order, eq)
        ck.check("NUMERIC", "end-to-end test of the necessary equations in the basis of %s on %d generic points x of "
                 "X with Re(F_sigma - F_pi*) = 0 on T (seed 39009; labels apart by > 0.05; the other gaps have "
                 "either sign, so the eight coordinates have mixed signs): relative |C0 + z Q| < 1e-9 and "
                 "relative |G| < 1e-9" % (nm, len(eq)),
                 len(eq) >= 40 and fit < 1e-9 and cub < 1e-9 and gmax < 1e-9,
                 "fit %.1e, cubic %.1e, |QG| %.1e, |G| %.1e, smallest relative |Q| %.1e" % (fit, cub, hes, gmax, qmin))
    pts = []
    for k, eps in enumerate((0.02, 0.1, 0.3)):
        pts += M.numeric_stratum_search(T, 39100 + k, 250, eps=eps)
    w40f, w70f = np.array([float(x) for x in w40]), np.array([float(x) for x in w70])
    for nm, wv, rd, order in (("delta40", w40f, rd40n, ORDER40), ("delta70", w70f, rd70n, ORDER70)):
        sel = [p_ for p_ in pts if float(wv @ p_["gaps"]) > 1e-4]
        R = rd["R"]
        A = np.array([F_VALUES[j] for j in order], float).T
        absG = {k: abs(v) for k, v in rd["G"].items()}
        relmax, okc, details, fitmax = 0.0, True, {}, 0.0
        for p_ in sel:
            co, *_ = np.linalg.lstsq(A, p_["gaps"], rcond=None)
            fitmax = max(fitmax, float(np.max(np.abs(A @ co - p_["gaps"]))))
            okc &= float(np.min(co)) > -1e-8
            dpos = order.index(40 if nm == "delta40" else 70)
            v = np.maximum(co, 0.0) / co[dpos]
            val = dict(zip(R.names[:8], [float(x) for x in v]))
            val["z"] = 0.0
            relmax = max(relmax, abs(R.evaluate(rd["G"], val)) / R.evaluate(absG, val))
            if nm == "delta40":
                J = doc_endpoint(val["a"], val["b"], val["c"], val["e"], val["p"])
                if J == 0.0:
                    av_, bv_, cv_, ev_, pv_ = (val[n_] for n_ in "abcep")
                    if av_ * cv_ > 1 and pv_ * (av_ * cv_ - 1) >= (av_ + 1) * (cv_ + 1):
                        key = "four-coefficient test p >= P_ac"
                    elif av_ * bv_ > 1 and ev_ * (av_ * bv_ - 1) >= (av_ + 1) * (bv_ + 1):
                        key = "four-coefficient test e >= E_ab"
                    else:
                        key = "allocation domain with e >= E_r0"
                    hit = True
                else:
                    hit = val["j"] >= J - 1e-6
                    key = "j >= J" if hit else "NONE"
            else:
                Jm, _ = energy_min(val["p"], val["b"], val["c"], val["q"])
                hs, ha = val["a"] * val["s"] >= 1 - 1e-7, val["a"] >= Jm - 1e-7
                hit = hs or ha
                key = ("as >= 1" if hs else "") + (" and " if hs and ha else "") + ("a >= J" if ha else "") or "NONE"
            details[key] = details.get(key, 0) + 1
            okc &= hit
        ck.info("stratum (gaps >= 0, zero on T, sum 1): %d points from two-stage LM searches (seeds 39100..39102, "
                "250 starts each, repulsion 0.02, 0.1, 0.3; residual < 1e-10); largest %s = %.2e; %d points with "
                "%s > 1e-4; among these the smallest normalised distance of two labels is %.2e and the numbers of "
                "contacts (|gap| < 1e-7) are %s"
                % (len(pts), nm, max([float(wv @ p_["gaps"]) for p_ in pts] or [float("nan")]), len(sel), nm,
                   min([p_["sep"] for p_ in sel] or [float("nan")]),
                   sorted(set(sum(1 for x_ in p_["gaps"] if abs(x_) < 1e-7) for p_ in sel))))
        ck.check("NUMERIC", "at every stratum point found with %s > 1e-4 (%d points; branch %s): the gap vector is a "
                 "nonnegative combination of the eight generators of the branch, relative |G| < 1e-6, and the "
                 "branch's conclusion holds%s"
                 % (nm, len(sel), "(1)" if nm == "delta40" else "(2)", "" if sel else "  [VOID: no such point was found]"),
                 okc and relmax < 1e-6 and fitmax < 1e-8,
                 "max relative |G| %.1e; fit %.1e; conclusions: %s" % (relmax, fitmax, details))


if __name__ == "__main__":
    main()
