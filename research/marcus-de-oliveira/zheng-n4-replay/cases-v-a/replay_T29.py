"""Independent exact replay of Case T_29 = {0,4,5,6,7,8} (group V) of D. Zheng's proof of the
Marcus-de Oliveira conjecture for n = 4.

Proof document replayed: MOC4_SUPP/proofs/manin_transfer_20261021/blind_moc4_T28_uniform_20261022.pdf
(5 pages; the document calls the set "T28").  Every polynomial is rebuilt here from the
definitions (the matrix P, the cubic C, the vectors f_j of the paper); nothing from the
supplement is executed or imported.  If MOC4_SUPP is set, a few JSON files are read as data and
COMPARED with the rebuilt objects.

Run:   MOC4_SUPP=/path/to/supplement  python replay_T29.py
Exit status 0 iff every check passes.  Checks tagged EXACT use integer/rational/symbolic
arithmetic only; checks tagged NUMERIC are floating-point falsification attempts.
"""
import math
import os
import random
import sys
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import replay_T29_common as M                                   # noqa: E402
from replay_T29_common import (Ring, p_add, p_mul, p_scale, p_nonneg, p_min, p_compose,   # noqa: E402
                               p_rename, p_divexact, F_VALUES, SH, MASK)

DIR = "proofs/manin_transfer_20261021/"
DOC = DIR + "blind_moc4_T28_uniform_20261022.pdf"
JSON34 = DIR + "blind_moc4_T28_bad34_full_radial_20261022.json"
JSON46 = DIR + "blind_moc4_T28_bad46_full_radial_20261022.json"
DOC_SHA256 = {
    DOC: "973f17a6585606dc3dd0c69e18cdd8bd2be4970ae3bb0540a5c9536327fa89bf",
    JSON34: "bb84ddd5a2cf9f31c9750cf2c7d10d31e4709790cce42c22c40395c7caee00c2",
    JSON46: "bd8b3f1997257ee87b436acee3e70e9f371a4387c776e6e91274728ad272b3b4",
}

T = (0, 4, 5, 6, 7, 8)
GENS = [0, 1, 2, 3, 4, 5, 6, 33, 34, 45, 46]
ORDER34 = [3, 6, 1, 33, 45, 34, 0, 5]        # (1): q f3 + p f6 + j f1 + c f33 + e f45 + d f34 + b f0 + a f5
NAMES34 = "q,p,j,c,e,d,b,a,h"
ORDER46 = [4, 3, 2, 33, 5, 6, 46, 45]        # (2): p f4 + r f3 + a f2 + b f33 + c f5 + q f6 + d f46 + s f45
NAMES46 = "p,r,a,b,c,q,d,s,h"
IDENTITIES = [([34, 46], [3, 6, 6, 33, 45]), ([2, 34], [0, 6, 33]), ([4, 34], [1, 3, 6, 33]),
              ([0, 46], [2, 3, 6, 45]), ([1, 46], [4, 6, 45])]
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
    """The document's reduction.  lam: dict generator -> nonnegative Fraction.  Returns
    (kind, coefficients) with kind in {'minor', '34', '46'}; the vector is unchanged."""
    lam = {j: Fraction(lam.get(j, 0)) for j in GENS}

    def use(lhs, rhs, amount):
        for j in lhs:
            lam[j] -= amount
        for j in rhs:
            lam[j] += amount
    use(*IDENTITIES[0], min(lam[34], lam[46]))
    if lam[34] > 0:
        use(*IDENTITIES[1], min(lam[2], lam[34]))
        use(*IDENTITIES[2], min(lam[4], lam[34]))
        return ("34" if lam[34] > 0 else "minor"), lam
    if lam[46] > 0:
        use(*IDENTITIES[3], min(lam[0], lam[46]))
        use(*IDENTITIES[4], min(lam[1], lam[46]))
        return ("46" if lam[46] > 0 else "minor"), lam
    return "minor", lam


def functional_on_cone(targets):
    """w in Q^24, zero on T, with w . f_j = targets.get(j, 0) for all eleven generators."""
    others = [s for s in range(24) if s not in T]
    A = [[Fraction(F_VALUES[j][s]) for s in others] + [Fraction(targets.get(j, 0))] for j in GENS]
    Rm, piv = M.rref(A)
    assert len(others) not in piv, "inconsistent targets"
    w = [Fraction(0)] * 24
    for r, pc in enumerate(piv):
        w[others[pc]] = Rm[r][-1]
    assert all(sum(w[s] * F_VALUES[j][s] for s in range(24)) == targets.get(j, 0) for j in GENS)
    return w


def main():
    import sympy as sp
    import numpy
    import scipy
    ck = M.Checker("Case T_29 = {0,4,5,6,7,8} (document 'T28'): independent exact replay")
    try:
        import flint
        fl = flint.__version__
    except ImportError:
        fl = "absent"
    ck.info("python %s, sympy %s, python-flint %s (integer determinants at test points only), numpy %s "
            "(int64 in det_linear_forms; floats in NUMERIC checks), scipy %s (NUMERIC only)"
            % (sys.version.split()[0], sp.__version__, fl, numpy.__version__, scipy.__version__))

    # ------------------------------------------------------------------------------------
    ck.section("A. Provenance and foundations")
    for rel in (DOC, JSON34, JSON46):
        M.check_document(ck, rel, DOC_SHA256[rel])
    M.check_foundations(ck)
    M.check_radial_necessity(ck)
    M.crosscheck_supplement_tables(ck)
    M.crosscheck_audit_library(ck, [29, 31])

    # ------------------------------------------------------------------------------------
    ck.section("B. Contact set")
    ck.check("EXACT", "the document's set {0,4,5,6,7,8} is the paper's T_29 (Appendix E) as a set of "
             "lexicographic indices: group element = identity", tuple(sorted(T)) == M.PAPER_T_V[29])
    r29, r31 = M.orbit_rep(T), M.orbit_rep(M.PAPER_T_V[31])
    ck.check("EXACT", "orbit representative of T_29 under G is (0,1,2,6,10,11); T_31 lies in another orbit",
             r29 == (0, 1, 2, 6, 10, 11) and r31 != r29, "T_29 -> %s, T_31 -> %s" % (r29, r31))

    # ------------------------------------------------------------------------------------
    ck.section("C. Generators and the reduction to two cones with one non-minor generator (p. 1)")
    gens = M.generators(T)
    ck.check("EXACT", "the f_j vanishing on T are exactly f0..f6, f33, f34, f45, f46", gens == GENS, "computed %s" % gens)
    nm = [j for j in gens if not M.F_IS_MINOR[j]]
    ck.check("EXACT", "the non-minor generators are f34 = h_{4321,+} and f46 = h_{3241,-}", nm == [34, 46],
             ", ".join("f%d = %s" % (j, M.f_name(j)) for j in GENS))
    d, rays = M.extreme_rays_of_contact_cone(T)
    ck.check("EXACT", "direct enumeration of the extreme rays of N_T = {Pc >= 0, (Pc)_T = 0} (dimension %d; no "
             "use of the 82-facet list) gives exactly these eleven vectors, so every c in N_T is a nonnegative "
             "combination of them" % d, rays == sorted(F_VALUES[j] for j in GENS), "%d rays" % len(rays))
    ck.check("EXACT", "the five identities of p. 1 hold on all 24 values (P is injective, so they hold in V_R): "
             "f34+f46 = f3+2f6+f33+f45, f2+f34 = f0+f6+f33, f4+f34 = f1+f3+f6+f33, f0+f46 = f2+f3+f6+f45, "
             "f1+f46 = f4+f6+f45",
             all(M.combo_values([(j, 1) for j in l]) == M.combo_values([(j, 1) for j in r]) for l, r in IDENTITIES))
    for nm_, order in (("(1)", ORDER34), ("(2)", ORDER46)):
        ck.check("EXACT", "the eight generators of %s, %s, are linearly independent and their common zero set "
                 "is exactly T" % (nm_, order),
                 M.rank([M.F_COORD[j] for j in order]) == 8
                 and tuple(s for s in range(24) if all(F_VALUES[j][s] == 0 for j in order)) == T)
    w34 = functional_on_cone({34: 1, 46: -1, 2: -1, 4: -1})
    w46 = functional_on_cone({46: 1, 34: -1, 0: -1, 1: -1})
    ck.check("EXACT", "delta34 := (coefficient of f34) - (coefficients of f46, f2, f4) and delta46 := (coefficient "
             "of f46) - (coefficients of f34, f0, f1) are well defined linear functionals on N_T (independent of "
             "the representation): delta34 = gap_2 - gap_13, delta46 = gap_12 - gap_2 - gap_16 - gap_18",
             {s: v for s, v in enumerate(w34) if v} == {2: 1, 13: -1}
             and {s: v for s, v in enumerate(w46) if v} == {12: 1, 2: -1, 16: -1, 18: -1})
    rng = random.Random(29001)
    ok = True
    kinds = {"minor": 0, "34": 0, "46": 0}
    for _ in range(4000):
        lam = {j: (Fraction(0) if rng.random() < 0.3 else Fraction(rng.randint(0, 12), rng.randint(1, 6)))
               for j in GENS}
        vec = M.combo_values(list(lam.items()))
        kind, out = reduce_to_one_nonminor(lam)
        kinds[kind] += 1
        ok &= M.combo_values(list(out.items())) == vec and all(v >= 0 for v in out.values())
        d34 = sum(w34[s] * vec[s] for s in range(24))
        d46 = sum(w46[s] * vec[s] for s in range(24))
        supp = {j for j, v in out.items() if v}
        if kind == "34":
            ok &= supp <= set(ORDER34) and out[34] == d34 > 0
        elif kind == "46":
            ok &= supp <= set(ORDER46) and out[46] == d46 > 0
        else:
            ok &= out[34] == 0 == out[46] and d34 <= 0 and d46 <= 0
    ck.check("EXACT", "the reduction of p. 1 (cancel the common part of f34, f46; then absorb f34 against f2, f4, "
             "or f46 against f0, f1) on 4000 random rational coefficient vectors with zeros (seed 29001): the "
             "vector is unchanged, all coefficients stay >= 0, and the result is a sum of minors (iff "
             "delta34, delta46 <= 0), or of the form (1) with d = delta34 > 0, or of the form (2) with "
             "d = delta46 > 0", ok, "outcomes %s" % kinds)
    M.check_minor_leibniz(ck, [j for j in GENS if j not in (34, 46)])

    # ------------------------------------------------------------------------------------
    ck.section("D. The coefficient conditions (3) and Lemma 2: orbit membership and algebra")
    rng = random.Random(29002)
    M.check_transport(ck, "old-H", (45, 1, 34), rng)
    M.check_transport(ck, "chain", (1, 0, 45, 34), rng)
    M.check_transport(ck, "four-coefficient", (6, 0, 5, 34), rng)
    M.check_transport(ck, "old-H", (33, 2, 46), rng)
    M.check_transport(ck, "old-H", (5, 2, 46), rng)
    M.check_transport(ck, "chain", (2, 4, 33, 46), rng)
    M.check_transport(ck, "chain", (2, 5, 6, 46), rng)
    a, b, c, e, j, p, q = sp.symbols("a b c e j p q", positive=True)
    four = lambda A_, B_, C_, D_: A_ * B_ * (C_ - D_) - D_ ** 2 * (A_ + B_ + C_ + D_)
    ck.check("EXACT", "(3) for (1) with d = 1: old-H t(c-d) >= d^2 with t-role f45 (coefficient e), c-role f1 "
             "(coefficient j) is e(j-1) >= 1, i.e. j >= 1 + 1/e; chain with centre f1, leaves f0, f45 is "
             "j >= 1/b + 1/e; four-coefficient with roles (f6, f0, f5; f34) = (p, b, a; 1) is "
             "bp(a-1) >= b + p + a + 1",
             sp.expand(four(p, b, a, 1) - (b * p * (a - 1) - (b + p + a + 1))) == 0
             and sp.simplify(e * (j - 1) - 1 - e * (j - 1 - 1 / e)) == 0)
    ck.check("EXACT", "bp(a-1) - (b+p+a+1) = a(bp-1) - (b+1)(p+1): the four-coefficient test fails iff L a < 1, "
             "L = max((bp-1)/((b+1)(p+1)), 0)",
             sp.expand(b * p * (a - 1) - (b + p + a + 1) - (a * (b * p - 1) - (b + 1) * (p + 1))) == 0)
    # Lemma 2
    b1, b2, d1, d2, x, y, t1, t2, tt = sp.symbols("beta1 beta2 delta1 delta2 x y t1 t2 t", positive=True)
    e1, e2 = b1 - d1, b2 - d2
    bH, bC, cH, cC = b * e1 / b1, b * d1 / b1, c * e2 / b2, c * d2 / b2
    ck.check("EXACT", "Lemma 2, two old-H inequalities (f33, f2; f46), (f5, f2; f46) with t-coefficients t1, t2 and "
             "non-minor shares t t1/(t1+t2), t t2/(t1+t2): the required coefficient of f2 is t + t^2/(t1+t2)",
             sp.cancel(sum(dd + dd ** 2 / ti for dd, ti in ((tt * t1 / (t1 + t2), t1), (tt * t2 / (t1 + t2), t2)))
                       - (tt + tt ** 2 / (t1 + t2))) == 0)
    ck.check("EXACT", "Lemma 2, Cauchy step: eta1^2/x + eta2^2/y - (eta1+eta2)^2/(x+y) = (eta1 y - eta2 x)^2/(xy(x+y)) >= 0",
             sp.cancel(d1 ** 2 / x + d2 ** 2 / y - (d1 + d2) ** 2 / (x + y)
                       - (d1 * y - d2 * x) ** 2 / (x * y * (x + y))) == 0)
    total = (e1 + e2) + e1 ** 2 / bH + e2 ** 2 / cH + d1 ** 2 / p + d1 ** 2 / bC + d2 ** 2 / cC + d2 ** 2 / q
    rhs6 = 1 - d1 - d2 + b1 ** 2 / b + b2 ** 2 / c + d1 ** 2 / p + d2 ** 2 / q
    ck.check("EXACT", "Lemma 2, bookkeeping: with b split as b eta1/beta1 (old-H) + b delta1/beta1 (chain 1) and c "
             "split likewise, (eta1+eta2) + eta1^2/b_H + eta2^2/c_H + delta1^2/p + delta1^2/b_ch + delta2^2/c_ch "
             "+ delta2^2/q equals the right side of (6) when beta1 + beta2 = 1; the shares of b, c and of the "
             "non-minor coefficient add up",
             sp.cancel((total - rhs6).subs(b2, 1 - b1)) == 0 and sp.cancel(bH + bC - b) == 0
             and sp.cancel(cH + cC - c) == 0 and sp.expand((e1 + e2 + d1 + d2 - 1).subs(b2, 1 - b1)) == 0)

    # ------------------------------------------------------------------------------------
    ck.section("E. Radial necessity in the generators of (1) and (2) (pp. 1-2)")
    rd34 = M.radial_data(ck, ORDER34, NAMES34)
    R34, C34, Q34, G34 = rd34["R"], rd34["C0"], rd34["Q"], rd34["G"]
    ck.check("EXACT", "(1): every coefficient of C0 is positive (%d terms) and the coefficient of d^3 is 1, so "
             "C0 >= d^3 > 0 and C0 + hQ = 0 forces Q != 0" % len(C34),
             all(v > 0 for v in C34.values()) and C34.get(R34.key((0, 0, 0, 0, 0, 3, 0, 0, 0))) == 1)
    Qdoc = R34.from_sympy(sp.sympify("a*d + a*e + a*q + c*d + c*e + c*q + d**2 + d*e + d*q - e*j + e*q"))
    ck.check("EXACT", "(1): Q = ad + ae + aq + cd + ce + cq + d^2 + de + dq - ej + eq (the document's formula)",
             Q34 == Qdoc)
    ck.check("EXACT", "(1): G has 1376 terms and degree 9 (the document's count for G_34)",
             len(G34) == 1376 and all(sum(R34.exps(k)) == 9 for k in G34))
    M.compare_radial_json(ck, rd34, JSON34)
    rd46 = M.radial_data(ck, ORDER46, NAMES46, seed=20261126)
    R46, C46, Q46, G46 = rd46["R"], rd46["C0"], rd46["Q"], rd46["G"]
    Qdoc = R46.from_sympy(sp.sympify("-(b+c)*a + (d+r+s)*(b+c) + d*(d+r+s) + r*s"))
    ck.check("EXACT", "(2): Q = -M a + U M + d U + r s with M = b + c, U = d + r + s (the document's formula)",
             Q46 == Qdoc)
    vb, vc = R46.var("b"), R46.var("c")
    guard = p_add(p_mul(p_add(vb, vc), C46), p_mul(p_mul(vb, vc), Q46), -1)
    ck.check("EXACT", "(2): C0 has one negative coefficient (-abc); M C0 - bc Q has 44 terms, all coefficients "
             "nonnegative, and contains b d^3 + c d^3.  So if M > 0 and Q = 0 then C0 > 0, and if M = 0 then "
             "Q = dU + rs > 0: in both cases C0 + hQ = 0 forces Q != 0",
             sorted(v for v in C46.values() if v < 0) == [-1]
             and C46.get(R46.key((0, 0, 1, 1, 1, 0, 0, 0, 0))) == -1
             and len(guard) == 44 and p_nonneg(guard)
             and guard.get(R46.key((0, 0, 0, 1, 0, 0, 3, 0, 0))) == 1
             and guard.get(R46.key((0, 0, 0, 0, 1, 0, 3, 0, 0))) == 1
             and R46.subs_zero(Q46, ["b", "c"]) == R46.from_sympy(sp.sympify("d*(d+r+s) + r*s")))
    ck.check("EXACT", "(2): G has 1549 terms and degree 9 (the document's count for G_46)",
             len(G46) == 1549 and all(sum(R46.exps(k)) == 9 for k in G46))
    M.compare_radial_json(ck, rd46, JSON46, alias={"sa": "s"})

    branch34(ck, rd34)
    branch46(ck, rd46)
    numeric_part(ck, rd34, rd46, w34, w46)
    ck.finish()


# ==========================================================================================
def branch34(ck, rd):
    import sympy as sp
    ck.section("F. Branch (1), non-minor f34: strict guards (p. 2)")
    R9 = rd["R"]
    Rs = Ring("q,p,j,c,e,b,a")
    G1 = p_rename(R9.subs_one(rd["G"], "d"), R9, Rs)
    ks = Rs.split(G1, "j")
    ck.check("EXACT", "at d = 1: G_34 = k0 + k1 j + k2 j^2", set(ks) == {0, 1, 2} and len(G1) == 1376,
             "terms of k0, k1, k2: %d, %d, %d" % (len(ks[0]), len(ks[1]), len(ks[2])))
    k0, k1, k2 = ks[0], ks[1], ks[2]
    va, vb, ve, vp = [Rs.var(x) for x in "abep"]
    F8 = p_divexact(k0, P1(vp))
    ck.check("EXACT", "k0 = (p+1) F8 with F8 of total degree 8 and degree 2 in a (exact division, re-multiplied)",
             F8 is not None and p_mul(F8, P1(vp)) == k0 and Rs.tdeg(F8) == 8 and Rs.deg(F8, "a") == 2,
             "F8 has %d terms" % len(F8 or {}))
    f8 = Rs.split(F8, "a")
    ck.check("EXACT", "the a^0 coefficient of F8 has 204 terms, all coefficients positive, constant term 1",
             len(f8[0]) == 204 and all(v > 0 for v in f8[0].values()) and f8[0].get(0) == 1)
    F0, F1, F2, L, a_, ap = sp.symbols("F0 F1 F2 L a aprime")
    Fq = lambda z: F0 + F1 * z + F2 * z ** 2
    ck.check("EXACT", "identity (4): (1 + L a')^2 F(a'/(1 + L a')) = F0 + (F1 + 2 L F0) a' + (F2 + L F1 + L^2 F0) a'^2; "
             "and a = a'/(1 + L a') is equivalent to a' = a/(1 - L a), 1 + L a' = 1/(1 - L a)",
             sp.cancel((1 + L * ap) ** 2 * Fq(ap / (1 + L * ap))
                       - (F0 + (F1 + 2 * L * F0) * ap + (F2 + L * F1 + L ** 2 * F0) * ap ** 2)) == 0
             and sp.cancel((a_ / (1 - L * a_)) / (1 + L * a_ / (1 - L * a_)) - a_) == 0
             and sp.cancel(1 + L * a_ / (1 - L * a_) - 1 / (1 - L * a_)) == 0)
    cut = p_add(p_add(k0, p_mul(vb, k1)), p_mul(mul(vb, vb), k2))
    ck.check("EXACT", "the cut value G_34(j = b) = k0 + b k1 + b^2 k2 has degree 2 in a", Rs.deg(cut, "a") == 2,
             "%d terms" % len(cut))
    Rc5 = Ring("B,T,E,Q,C")
    B5, T5 = Rc5.var("B"), Rc5.var("T")
    ren = {"e": "E", "q": "Q", "c": "C"}
    doc_sizes = {("base", "lt"): (291, 542, 556), ("base", "gt"): (276, 467, 426),
                 ("cut", "lt"): (598, 869, 846), ("cut", "gt"): (736, 1131, 1377)}
    sB, sT = sp.symbols("B T")
    for label, poly in (("base", F8), ("cut", cut)):
        f = Rs.split(poly, "a")
        f0, f1, f2 = f.get(0, {}), f.get(1, {}), f.get(2, {})
        if label == "cut":
            hh, _ = p_compose(f0, Rs, Rc5, {"b": (B5, P1(B5))}, rename=dict(ren, p="T"))
            b1 = Rs.subs_one(f0, "b")
            ck.check("EXACT", "the a^0 coefficient of G_34(b): with b = B/(1+B) its numerator has 460 terms, "
                     "nonnegative coefficients and constant term 1 (so it is > 0 for 0 <= b < 1); at b = 1 it "
                     "has 92 terms, nonnegative coefficients and constant term 16 (> 0 at b = 1)",
                     len(hh) == 460 and p_nonneg(hh) and hh.get(0) == 1
                     and len(b1) == 92 and p_nonneg(b1) and b1.get(0) == 16)
        for bp in ("lt", "gt"):
            ln = p_add(p_mul(vb, vp), ONE, -1) if bp == "gt" else {}
            ld = mul(P1(vb), P1(vp))
            coeff = [f0, p_add(p_mul(f1, ld), p_scale(p_mul(ln, f0), 2)),
                     p_add(p_add(p_mul(f2, mul(ld, ld)), p_mul(mul(ln, ld), f1)), p_mul(mul(ln, ln), f0))]
            bn, bd = (B5, ONE) if label == "base" else (B5, P1(B5))
            pn, pd = (p_mul(T5, bd), p_mul(bn, P1(T5))) if bp == "lt" else (p_mul(P1(T5), bd), bn)
            sb_, sp_ = Rc5.to_sympy(bn) / Rc5.to_sympy(bd), Rc5.to_sympy(pn) / Rc5.to_sympy(pd)
            okc = sp.cancel(sb_ * sp_ - (sT / (1 + sT) if bp == "lt" else 1 + sT)) == 0
            sizes, okn = [], True
            for h_ in coeff:
                out, _ = p_compose(h_, Rs, Rc5, {"b": (bn, bd), "p": (pn, pd)}, rename=ren)
                sizes.append(len(out))
                okn &= p_nonneg(out)
            ck.check("EXACT", "%s, chart b = %s, bp = %s: the three transformed coefficient functions of (4), "
                     "F0, (b+1)(p+1) F1 + 2 l F0, ((b+1)(p+1))^2 F2 + l (b+1)(p+1) F1 + l^2 F0 with l = %s, "
                     "have nonnegative coefficients (document: %s terms)"
                     % ("F8" if label == "base" else "G_34(b)", "B" if label == "base" else "B/(1+B)",
                        "T/(1+T) < 1" if bp == "lt" else "1+T > 1", "0" if bp == "lt" else "bp - 1",
                        doc_sizes[(label, bp)]),
                     okc and okn and tuple(sizes) == doc_sizes[(label, bp)], "terms %s" % sizes)
    ck.info("PROOF-READ: for L a < 1, F(a) = [F0 + ..a' + ..a'^2] (1 - L a)^2 >= F0 (1 - L a)^2 with a' = a/(1-La) >= 0. "
            "So k0 = (p+1) F8 >= (p+1)(1 - L a)^2 > 0, and G_34(b) >= (its a^0 coefficient)(1 - L a)^2 > 0 for "
            "0 <= b <= 1: the guard (5).  The charts need b > 0 and bp != 1; the transformed coefficient "
            "functions are polynomials in (b, p, e, q, c) (after multiplying by (b+1)(p+1) > 0), so their "
            "signs extend to b = 0 and bp = 1 by continuity, and the two a^0 coefficients are checked to be "
            "strictly positive on the closed sets directly.")
    k2e0 = Rs.subs_zero(k2, ["e"])
    k1e0 = Rs.subs_zero(k1, ["e"])
    ck.check("EXACT", "e = 0: k2 vanishes identically and k1 has nonnegative coefficients, so G_34 >= k0 > 0 at "
             "e = 0 (given L a < 1)", not k2e0 and p_nonneg(k1e0), "k1 at e = 0: %d terms" % len(k1e0))

    ck.section("G. Branch (1): interval certificates (pp. 2-3)")
    Rc = Ring("B,T,A,E,Q,C")
    B_, T_, A_ = Rc.var("B"), Rc.var("T"), Rc.var("A")
    sB, sT, sA = sp.symbols("B T A")
    doc = {("lt", "lt"): (1799, 1799, 2488, 2128), ("lt", "gt"): (3555, 3555, 5043, 4590),
           ("gt", "lt"): (2837, 3201), ("gt", "gt"): (4167, 4663)}
    doc_groups = {("gt", "lt", "endpoint"): 19, ("gt", "gt", "middle"): 24, ("gt", "gt", "endpoint"): 42}
    total_groups = 0
    for bcase in ("lt", "gt"):
        for bp in ("lt", "gt"):
            bn, bd = (B_, P1(B_)) if bcase == "lt" else (P1(B_), ONE)
            pn, pd = (p_mul(T_, bd), p_mul(bn, P1(T_))) if bp == "lt" else (p_mul(P1(T_), bd), bn)
            if bp == "lt":
                an, ad = A_, ONE
            else:
                an, ad = mul(p_add(bn, bd), p_add(pn, pd), A_), mul(bd, pd, T_, P1(A_))
            sb_, sp_, sa_ = [Rc.to_sympy(n_) / Rc.to_sympy(d_) for n_, d_ in ((bn, bd), (pn, pd), (an, ad))]
            okc = sp.cancel(sb_ - (sB / (1 + sB) if bcase == "lt" else 1 + sB)) == 0
            okc &= sp.cancel(sb_ * sp_ - (sT / (1 + sT) if bp == "lt" else 1 + sT)) == 0
            if bp == "gt":
                okc &= sp.cancel(sa_ * (sb_ * sp_ - 1) / ((sb_ + 1) * (sp_ + 1)) - sA / (1 + sA)) == 0
            subs = {"a": (an, ad), "b": (bn, bd), "p": (pn, pd)}
            if bcase == "lt":
                ps = [("left middle 2 k0 + b k1", p_add(p_scale(k0, 2), p_mul(vb, k1))),
                      ("cut G_34(b)", cut),
                      ("right middle 2e k0 + (be+e+1) k1 + 2b(e+1) k2",
                       p_add(p_add(p_scale(p_mul(ve, k0), 2), p_mul(p_add(p_add(mul(vb, ve), ve), ONE), k1)),
                             p_scale(p_mul(mul(vb, P1(ve)), k2), 2))),
                      ("endpoint e^2 k0 + (e+1)e k1 + (e+1)^2 k2",
                       p_add(p_add(p_mul(mul(ve, ve), k0), p_mul(mul(P1(ve), ve), k1)),
                             p_mul(mul(P1(ve), P1(ve)), k2)))]
            else:
                nn, dd = p_add(vb, ve), mul(vb, ve)
                ps = [("middle 2be k0 + (b+e) k1", p_add(p_scale(p_mul(dd, k0), 2), p_mul(nn, k1))),
                      ("endpoint (be)^2 k0 + (b+e)be k1 + (b+e)^2 k2",
                       p_add(p_add(p_mul(mul(dd, dd), k0), p_mul(mul(nn, dd), k1)), p_mul(mul(nn, nn), k2)))]
            for idx, (label, poly) in enumerate(ps):
                out, _ = p_compose(poly, Rs, Rc, subs, rename={"e": "E", "q": "Q", "c": "C"})
                want = doc[(bcase, bp)][idx]
                chart = "b = %s, bp = %s%s" % ("B/(1+B) <= 1" if bcase == "lt" else "1+B >= 1",
                                               "T/(1+T), a = A" if bp == "lt" else "1+T, L a = A/(1+A)", "")
                if p_nonneg(out):
                    ck.check("EXACT", "chart %s: numerator of the %s has nonnegative coefficients (document: %d terms)"
                             % (chart, label, want), okc and len(out) == want, "%d terms" % len(out))
                else:
                    ngr, bad = M.grouped_univariate(Rc, out, "E")
                    okg = all(M.uni_nonneg(cf) and M.uni_sympy_nonneg(cf) for _, cf in bad)
                    total_groups += len(bad)
                    key = (bcase, bp, label.split()[0])
                    ck.check("EXACT", "chart %s: numerator of the %s (%d terms; document: %d) is not coefficientwise "
                             "nonnegative; grouped by the monomial in (B,T,A,Q,C), %d of its %d coefficient "
                             "polynomials in E have a negative coefficient (document: %d) and each of them is "
                             ">= 0 on [0, inf) (own squarefree decomposition + Sturm; cross-checked with sympy)"
                             % (chart, label, len(out), want, len(bad), ngr, doc_groups.get(key, -1)),
                             okc and len(out) == want and okg and len(bad) == doc_groups.get(key, -1),
                             "max degree in E: %d" % max(max(cf) for _, cf in bad))
    ck.check("EXACT", "the univariate groups number 19 + 24 + 42 = 85 (the document's total)", total_groups == 85)
    x0, x1, x2, s_, Jv_, l_ = sp.symbols("x0 x1 x2 s J l")
    quad = lambda jj: x0 + x1 * jj + x2 * jj ** 2
    id2 = sp.expand(quad(l_ + s_ * (Jv_ - l_)) - (quad(l_) * (1 - s_) ** 2
                                                  + (x0 + x1 * (l_ + Jv_) / 2 + x2 * l_ * Jv_) * 2 * s_ * (1 - s_)
                                                  + quad(Jv_) * s_ ** 2))
    eb, ee = sp.symbols("b e", positive=True)
    ck.check("EXACT", "Bernstein form on [l, J]: coefficients G(l), k0 + k1 (l+J)/2 + k2 l J, G(J); with l = 0, "
             "J = b: middle (2 k0 + b k1)/2; with l = b, J = 1 + 1/e: middle [2e k0 + (be+e+1) k1 + 2b(e+1) k2]/(2e), "
             "end [e^2 k0 + e(e+1) k1 + (e+1)^2 k2]/e^2; with l = 0, J = 1/e + 1/b: middle [2be k0 + (b+e) k1]/(2be)",
             id2 == 0
             and sp.cancel((x0 + x1 * (eb + 1 + 1 / ee) / 2 + x2 * eb * (1 + 1 / ee)) * 2 * ee
                           - (2 * ee * x0 + (eb * ee + ee + 1) * x1 + 2 * eb * (ee + 1) * x2)) == 0
             and sp.cancel(quad(1 + 1 / ee) * ee ** 2 - (ee ** 2 * x0 + ee * (ee + 1) * x1 + (ee + 1) ** 2 * x2)) == 0
             and sp.cancel((x0 + x1 * (1 / ee + 1 / eb) / 2) * 2 * eb * ee - (2 * eb * ee * x0 + (eb + ee) * x1)) == 0)
    ck.info("PROOF-READ: b <= 1, e > 0, La < 1: on [0,b], G_34 >= k0 (1-s)^2 + G_34(b) s^2 > 0; on [b, 1+1/e), "
            "G_34 >= G_34(b)(1-s)^2 > 0; so G_34 = 0 forces j >= 1 + 1/e (old-H).  b >= 1, e > 0, La < 1: on "
            "[0, 1/e+1/b), G_34 >= k0 (1-s)^2 > 0; so G_34 = 0 forces j >= 1/e + 1/b (chain).  The charts are "
            "open; the numerators are polynomials in the coefficients times positive denominators, so the "
            "nonnegativity extends to b = 0, b = 1, bp = 1 by continuity, and the strict guards hold on the "
            "closed sets.")
    return Rs, G1, k0, k1, k2, cut


# ==========================================================================================
def charts46(Rc):
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


def branch46(ck, rd):
    import sympy as sp
    ck.section("H. Branch (2), non-minor f46: the allocation charts (pp. 3-4)")
    R9 = rd["R"]
    Rs = Ring("p,r,a,b,c,q,s")
    G1 = p_rename(R9.subs_one(rd["G"], "d"), R9, Rs)
    ks = Rs.split(G1, "a")
    ck.check("EXACT", "at d = 1: G_46 = k0 + k1 a + k2 a^2", set(ks) == {0, 1, 2} and len(G1) == 1549,
             "terms of k0, k1, k2: %d, %d, %d" % (len(ks[0]), len(ks[1]), len(ks[2])))
    k0, k1, k2 = ks[0], ks[1], ks[2]
    ck.check("EXACT", "k0 has 612 terms, all coefficients positive, constant term 1 (so k0 >= 1)",
             len(k0) == 612 and all(v > 0 for v in k0.values()) and k0.get(0) == 1)
    Rc = Ring("L,T,U,W,s,r")
    sL, sT, sU, sW = sp.symbols("L T U W")
    theta = sT / (1 + sT)
    charts = charts46(Rc)

    def nonneg_rat(expr):
        n_, d_ = sp.fraction(sp.cancel(expr))
        cn, cd = sp.Poly(n_, sL, sT, sU, sW).coeffs(), sp.Poly(d_, sL, sT, sU, sW).coeffs()
        return (all(x >= 0 for x in cn) and all(x > 0 for x in cd)) or \
            (all(x <= 0 for x in cn) and all(x < 0 for x in cd))

    for case, (subs, nn, dd, states, swapped) in charts.items():
        pp, bb, cc, qq = [Rc.to_sympy(subs[x][0]) / Rc.to_sympy(subs[x][1]) for x in "pbcq"]
        beta = [theta, 1 - theta]
        if swapped:
            beta.reverse()
        delta = [pp / 2 if states[0] == "H" else beta[0], qq / 2 if states[1] == "H" else beta[1]]
        Jc = Rc.to_sympy(nn) / Rc.to_sympy(dd)
        energy = 1 - delta[0] - delta[1] + beta[0] ** 2 / bb + beta[1] ** 2 / cc + delta[0] ** 2 / pp + delta[1] ** 2 / qq
        feas = all(nonneg_rat(x) for x in (delta[0], delta[1], beta[0] - delta[0], beta[1] - delta[1], Jc))
        der = [2 * beta[0] / bb, 2 * beta[1] / cc]
        for i, res in enumerate((pp, qq)):
            if states[i] == "C":
                der[i] += 2 * beta[i] / res - 1
        # leaf conditions: H: resource/2 < beta ; C: beta <= resource/2
        leaf = True
        for i, res in enumerate((pp, qq)):
            leaf &= nonneg_rat(beta[i] - res / 2) if states[i] == "H" else nonneg_rat(res / 2 - beta[i])
        ck.check("EXACT", "chart %s: with theta = T/(1+T), (beta1, beta2) = %s and delta_i = %s, the point "
                 "(beta, delta) is feasible (0 <= delta_i <= beta_i as rational functions with nonnegative "
                 "coefficients), the right side of (6) equals the chart's J = n/d0 > 0 (energy identity), the "
                 "leaf conditions hold, and d/dbeta1 of the reduced function vanishes (stationarity)"
                 % (case, "(1-theta, theta)" if swapped else "(theta, 1-theta)",
                    "(" + ", ".join("resource_i/2" if st == "H" else "beta_i" for st in states) + ")"),
                 feas and leaf and sp.cancel(energy - Jc) == 0 and sp.cancel(der[0] - der[1]) == 0)
    ck.info("PROOF-READ (coverage): for positive (p,b,c,q) the function of beta1 obtained after delta_i = "
            "min(resource_i/2, beta_i) is C^1 and strictly convex on [0,1] with negative derivative at 0 and "
            "positive derivative at 1, so it has a stationary point theta in (0,1).  With T = theta/(1-theta), "
            "L = theta/b (or (1-theta)/c in the exchanged charts) and U, W read off from p and q, the "
            "stationarity equation is exactly the chart's formula for the remaining coefficient; the leaf "
            "types (H: resource < 2 beta_i; C: resource >= 2 beta_i) and, in CC, the order of p/(2 beta1) and "
            "q/(2 beta2) select one of the five charts with U, W >= 0.  What the proof uses is only: J_chart is "
            "the value of the right side of (6) at a feasible point (so a >= J_chart gives the support "
            "inequality by Lemma 2), and every positive quadruple lies in the image of a chart.")

    ck.section("I. Branch (2): sign certificates in the five charts and the bound (8) (p. 4)")
    sS, sR, sUu = SH * Rc.pos["s"], SH * Rc.pos["r"], SH * Rc.pos["U"]

    def has_aux(key):
        return bool(((key >> sS) & MASK) or ((key >> sR) & MASK))

    doc_exc = {"HH": 18, "HC": 21, "CC": 26, "CCswap": 38}
    for case, (subs, nn, dd, states, swapped) in charts.items():
        degs = {}
        for kk in (k0, k1, k2):
            for nm in subs:
                degs[nm] = max(degs.get(nm, 0), Rs.deg(kk, nm))
        sk = [p_compose(kk, Rs, Rc, subs, degs=degs)[0] for kk in (k0, k1, k2)]
        P0 = p_mul(mul(dd, dd), sk[0])
        P1_ = p_mul(dd, p_add(p_scale(p_mul(dd, sk[0]), 2), p_mul(nn, sk[1])))
        P2 = p_add(p_add(p_mul(mul(dd, dd), sk[0]), p_mul(mul(nn, dd), sk[1])), p_mul(mul(nn, nn), sk[2]))
        polys = [P0, P1_, P2]
        aux0 = all(c >= 0 for P_ in polys for k, c in P_.items() if not has_aux(k))
        ck.check("EXACT", "chart %s: with a = J Y/(1+Y) and one common denominator, the numerator of G_46 is "
                 "P0 + P1 Y + P2 Y^2 as in (7); the part of P0, P1, P2 free of r and s has nonnegative "
                 "coefficients" % case, aux0, "terms of P0, P1, P2: %s" % [len(P_) for P_ in polys])
        if case in ("HH", "HC"):
            um = MASK << sUu
            groups = {}
            for i, P_ in enumerate(polys):
                for k, c in P_.items():
                    groups.setdefault(k & ~um, [{}, {}, {}])[i][(k >> sUu) & MASK] = c
            exc = 0
            ok = True
            kinds = {"B >= 0": 0, "4AC - B^2 >= 0": 0}
            for key, (cc_, bb_, aa_) in groups.items():
                if all(v >= 0 for col in (cc_, bb_, aa_) for v in col.values()):
                    continue
                exc += 1
                ok &= has_aux(key) and M.uni_nonneg(aa_) and M.uni_nonneg(cc_)
                ok &= M.uni_sympy_nonneg(aa_) and M.uni_sympy_nonneg(cc_)
                if M.uni_nonneg(bb_):
                    kinds["B >= 0"] += 1
                    ok &= M.uni_sympy_nonneg(bb_)
                else:
                    disc = M.uni_sub(M.uni_mul({k: 4 * v for k, v in aa_.items()}, cc_), M.uni_mul(bb_, bb_))
                    ok &= M.uni_nonneg(disc) and M.uni_sympy_nonneg(disc)
                    kinds["4AC - B^2 >= 0"] += 1
            ck.check("EXACT", "chart %s: grouped by the monomial in (L,T,W,s,r), the numerator is a sum of "
                     "monomial * (A(U) Y^2 + B(U) Y + C(U)); %d of the %d groups have a negative coefficient "
                     "(document: %d); in each of them A >= 0, C >= 0 and either B >= 0 or 4AC - B^2 >= 0 on "
                     "U >= 0 (own Sturm certificates, cross-checked with sympy), and each carries a positive "
                     "power of r or s" % (case, exc, len(groups), doc_exc[case]),
                     ok and exc == doc_exc[case], str(kinds))
        elif case == "CH":
            ck.check("EXACT", "chart CH: P0, P1, P2 have nonnegative coefficients", all(p_nonneg(P_) for P_ in polys))
        else:
            ngr, bad = M.grouped_univariate(Rc, P2, "L")
            okg = all(M.uni_nonneg(cf) and M.uni_sympy_nonneg(cf) and has_aux(key) for key, cf in bad)
            ck.check("EXACT", "chart %s: P0 and P1 have nonnegative coefficients; P2, grouped by the monomial in "
                     "(T,U,W,s,r), has %d coefficient polynomials in L with a negative coefficient out of %d "
                     "(document: %d), each >= 0 on [0, inf) and each carrying a positive power of r or s"
                     % (case, len(bad), ngr, doc_exc[case]),
                     p_nonneg(P0) and p_nonneg(P1_) and okg and len(bad) == doc_exc[case])
    ck.info("PROOF-READ: (1+Y)^2 d0^2 Delta G_46 = P0 + P1 Y + P2 Y^2 >= (part free of r, s) >= P0 at r = s = 0 "
            "= d0^2 Delta k0(p,b,c,q,0,0); with 1/(1+Y) = 1 - a/J this is (8): G_46(a,r,s) >= k0(p,b,c,q,0,0) "
            "(1 - a/J)^2 >= (1 - a/J)^2 > 0 for 0 <= a < J, all r, s >= 0 and positive p, b, c, q.")
    return Rs, G1, k0


# ==========================================================================================
def rhs6_min(p, b, c, q, n=4000):
    """NUMERIC: min over beta1 in [0,1] of the right side of (6) with the optimal delta_i =
    min(resource/2, beta_i); zero coefficients follow the document's conventions.  Grid +
    golden-section refinement of a convex function.  Returns (J, beta1)."""
    def h(res, beta):
        if res <= 0:
            return 0.0                              # delta must be 0
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
    cand = min([(phi(x), x), (vals[i], xs[i]), (phi(0.0), 0.0), (phi(1.0), 1.0)])
    return cand


def numeric_part(ck, rd34, rd46, w34, w46):
    import numpy as np
    ck.section("J. NUMERIC falsification attempts (evidence only)")
    f = lambda idx: np.array(F_VALUES[idx], float)
    rng = random.Random(29003)
    ex = lambda: math.exp(rng.gauss(0.0, 1.0))
    worst = {}

    def upd(name, vec, seed):
        worst[name] = min(worst.get(name, 1e9), M.numeric_min_phi(vec, seed, 3000, 5))
    for trial in range(10):
        a, b, c, e, p, q = [ex() for _ in range(6)]
        upd("old-H (f45, f1; f34)", e * f(45) + (1 + 1 / e) * f(1) + f(34), 100 + trial)
        upd("chain (f1; f0, f45; f34)", (1 / e + 1 / b) * f(1) + b * f(0) + e * f(45) + f(34), 200 + trial)
        if b * p > 1:
            upd("four-coefficient (f6, f0, f5; f34)",
                p * f(6) + b * f(0) + (b + 1) * (p + 1) / (b * p - 1) * f(5) + f(34), 300 + trial)
        upd("old-H (f33, f2; f46)", b * f(33) + (1 + 1 / b) * f(2) + f(46), 400 + trial)
        upd("old-H (f5, f2; f46)", c * f(5) + (1 + 1 / c) * f(2) + f(46), 500 + trial)
        upd("chain (f2; f4, f33; f46)", (1 / p + 1 / b) * f(2) + p * f(4) + b * f(33) + f(46), 600 + trial)
        upd("chain (f2; f5, f6; f46)", (1 / c + 1 / q) * f(2) + c * f(5) + q * f(6) + f(46), 700 + trial)
        beta1 = rng.random()
        dl1, dl2 = rng.random() * beta1, rng.random() * (1 - beta1)
        aa = 1 - dl1 - dl2 + beta1 ** 2 / b + (1 - beta1) ** 2 / c + dl1 ** 2 / p + dl2 ** 2 / q
        upd("Lemma 2 at a random feasible (beta, delta)",
            p * f(4) + b * f(33) + c * f(5) + q * f(6) + aa * f(2) + f(46), 800 + trial)
        Jm, _ = rhs6_min(p, b, c, q)
        upd("Lemma 2 at the minimum J(p,b,c,q)",
            p * f(4) + b * f(33) + c * f(5) + q * f(6) + Jm * f(2) + f(46), 900 + trial)
    for nm, wv in worst.items():
        ck.check("NUMERIC", "%s at its threshold (coefficients exp(N(0,1)), 10 draws, seed 29003; 3000 random "
                 "rotations + 5 BFGS descents each): min Phi found >= -1e-9" % nm, wv > -1e-9, "min = %.3e" % wv)
    ck.info("controls (the search does detect negative values): min Phi found for f34 alone %.4f, for f46 alone "
            "%.4f, for f45 + f1 + f34 (j = 1 < 1 + 1/e = 2) %.4f, for f33 + 1.2 f2 + f46 (a = 1.2 < 1 + 1/b = 2) %.4f, "
            "for f4 + f33 + 0.2 f2 + f46 (a = 0.2 < 1/p + 1/b = 2) %.4f"
            % (M.numeric_min_phi(f(34), 96, 3000, 5), M.numeric_min_phi(f(46), 97, 3000, 5),
               M.numeric_min_phi(f(45) + f(1) + f(34), 98, 3000, 5),
               M.numeric_min_phi(f(33) + 1.2 * f(2) + f(46), 99, 6000, 8),
               M.numeric_min_phi(f(4) + f(33) + 0.2 * f(2) + f(46), 99, 6000, 8)))

    # ---- branch (1): exact rational sampling
    R9 = rd34["R"]
    Rs = Ring("q,p,j,c,e,b,a")
    G1 = p_rename(R9.subs_one(rd34["G"], "d"), R9, Rs)
    k0 = Rs.split(G1, "j")[0]
    rng = random.Random(29004)
    vals = [Fraction(0), Fraction(0), Fraction(1), Fraction(1, 2), Fraction(2), Fraction(1, 9), Fraction(7),
            Fraction(1, 40), Fraction(30)]

    def rnd():
        return rng.choice(vals) if rng.random() < 0.45 else Fraction(rng.randint(0, 60), rng.randint(1, 12))
    n_g = bad_g = n_k = bad_k = n_cut = bad_cut = 0
    for _ in range(5000):
        b, p, e, q, c = rnd(), rnd(), rnd(), rnd(), rnd()
        Lv = max((b * p - 1) / ((b + 1) * (p + 1)), Fraction(0))
        if Lv > 0:
            a = Fraction(rng.randint(0, 999), 1000) / Lv
        else:
            a = rnd()
        pt = dict(q=q, p=p, j=Fraction(0), c=c, e=e, b=b, a=a)
        k0v = Rs.evaluate(k0, pt)
        n_k += 1
        bad_k += k0v < (p + 1) * (1 - Lv * a) ** 2
        if b <= 1:
            pt["j"] = b
            n_cut += 1
            bad_cut += Rs.evaluate(G1, pt) <= 0
        s = Fraction(rng.randint(0, 999), 1000)
        if e == 0:
            jv = rnd()
        elif b <= 1:
            jv = s * (1 + 1 / e)
        else:
            jv = s * (1 / e + 1 / b)
        pt["j"] = jv
        n_g += 1
        bad_g += Rs.evaluate(G1, pt) <= 0
    ck.check("NUMERIC", "branch (1), exact evaluation at %d random rational points with L a < 1 (seed 29004, zeros "
             "and boundary values included): k0 >= (p+1)(1 - L a)^2 always; G_34(j = b) > 0 at the %d points "
             "with b <= 1; G_34 > 0 at j = s J, s in [0,1), J = 1 + 1/e (b <= 1) or 1/e + 1/b (b >= 1), and for "
             "arbitrary j >= 0 when e = 0" % (n_k, n_cut), bad_k == 0 and bad_cut == 0 and bad_g == 0,
             "violations: %d, %d, %d" % (bad_k, bad_cut, bad_g))

    # ---- branch (2): floats, J from a direct one-dimensional minimisation (no charts)
    R9b = rd46["R"]
    Rs2 = Ring("p,r,a,b,c,q,s")
    G2 = p_rename(R9b.subs_one(rd46["G"], "d"), R9b, Rs2)
    k0b = Rs2.split(G2, "a")[0]
    rng = random.Random(29005)
    n8 = bad8 = ninf = 0
    minratio = float("inf")
    for _ in range(2500):
        def rv():
            u = rng.random()
            if u < 0.12:
                return 0.0
            return math.exp(rng.gauss(0.0, 1.5))
        p, b, c, q, r, s = rv(), rv(), rv(), rv(), rv(), rv()
        Jm, _ = rhs6_min(p, b, c, q)
        if Jm == float("inf"):
            ninf += 1
            av = math.exp(rng.gauss(0.0, 2.0))
            frac = 0.0
        else:
            frac = rng.random() * 0.98
            av = frac * Jm
        val = Rs2.evaluate(G2, dict(p=p, r=r, a=av, b=b, c=c, q=q, s=s))
        low = Rs2.evaluate(k0b, dict(p=p, r=0.0, a=0.0, b=b, c=c, q=q, s=0.0)) * (1 - frac) ** 2
        n8 += 1
        if not val >= low * (1 - 1e-9):
            bad8 += 1
        minratio = min(minratio, val / low)
    ck.check("NUMERIC", "branch (2), (8) in floating point at %d random points (seed 29005; p, b, c, q, r, s = 0 "
             "with probability 0.12 each, else exp(N(0,1.5)); J computed by direct minimisation of the right "
             "side of (6), not through the charts; %d points with J = infinity, a arbitrary): "
             "G_46(a,r,s) >= k0(p,b,c,q,0,0)(1 - a/J)^2 (relative tolerance 1e-9)" % (n8, ninf), bad8 == 0,
             "min ratio %.6f" % minratio)
    # chart coverage
    rng = random.Random(29006)
    Rc = Ring("L,T,U,W,s,r")
    charts = charts46(Rc)
    cover = Counter_()
    worst_err = 0.0
    okcov = True
    for _ in range(400):
        p, b, c, q = [math.exp(rng.gauss(0.0, 1.3)) for _ in range(4)]
        Jm, th = rhs6_min(p, b, c, q, n=20000)
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
        cover[case] += 1
        subs, nn, dd, states, swapped = charts[case]
        pt = dict(L=L, T=T_, U=U, W=W, s=0.0, r=0.0)
        img = [Rc.evaluate(subs[x][0], pt) / Rc.evaluate(subs[x][1], pt) for x in "pbcq"]
        Jc = Rc.evaluate(nn, pt) / Rc.evaluate(dd, pt)
        err = max(abs(u_ - v_) / v_ for u_, v_ in zip(img + [Jc], [p, b, c, q, Jm]))
        worst_err = max(worst_err, err)
        okcov &= min(L, T_, U, W) >= -1e-6 and err < 1e-4
    ck.check("NUMERIC", "coverage of the five charts: for 400 random positive (p,b,c,q) (seed 29006) the "
             "minimiser theta of the right side of (6) (computed directly) determines a chart and parameters "
             "L, T > 0, U, W >= 0 whose image is (p,b,c,q) and whose J equals the minimum (relative error "
             "< 1e-4)", okcov, "charts hit: %s; worst relative error %.1e" % (dict(cover), worst_err))

    # ---- genuine points of X
    eq = M.numeric_equality_points(T, 29009, 80)
    for nm, rd, order in (("(1)", rd34, ORDER34), ("(2)", rd46, ORDER46)):
        fit, cub, hes, qmin, gmax = M.numeric_radial_test(rd, order, eq)
        ck.check("NUMERIC", "end-to-end test of the necessary equations in the basis of %s on %d generic points x of "
                 "X with Re(F_sigma - F_pi*) = 0 on T (seed 29009; labels apart by > 0.05; the other gaps have "
                 "either sign, so the eight coordinates have mixed signs): relative |C0 + H Q| < 1e-9 and "
                 "relative |G| < 1e-9" % (nm, len(eq)),
                 len(eq) >= 40 and fit < 1e-9 and cub < 1e-9 and gmax < 1e-9,
                 "fit %.1e, cubic %.1e, |QG| %.1e, |G| %.1e, smallest relative |Q| %.1e" % (fit, cub, hes, gmax, qmin))
    pts = []
    for k, eps in enumerate((0.02, 0.1, 0.3, 0.005)):
        pts += M.numeric_stratum_search(T, 29100 + k, 300, eps=eps)
    w34f, w46f = np.array([float(x) for x in w34]), np.array([float(x) for x in w46])
    for nm, wv in (("delta34", w34f), ("delta46", w46f)):
        sel = [p_ for p_ in pts if float(wv @ p_["gaps"]) > 1e-4]
        pos = [p_["gaps"] for p_ in sel]
        res = [stratum_ok(g, rd34, rd46) for g in pos]
        ck.info("stratum (gaps >= 0, zero on T, sum 1): %d points from two-stage LM searches (seeds 29100..29103, "
                "300 starts each, repulsion 0.02, 0.1, 0.3, 0.005; residual < 1e-10); largest %s = %.2e; %d points "
                "with %s > 1e-4; among these the smallest normalised distance of two labels is %.2e and the numbers "
                "of contacts (|gap| < 1e-7) are %s"
                % (len(pts), nm, max([float(wv @ p_["gaps"]) for p_ in pts] or [float("nan")]), len(pos), nm,
                   min([p_["sep"] for p_ in sel] or [float("nan")]),
                   sorted(sum(1 for x_ in p_["gaps"] if abs(x_) < 1e-7) for p_ in sel)))
        ck.check("NUMERIC", "at every stratum point found with %s > 1e-4 (%d points; branch %s): relative |G| < 1e-7 "
                 "in the basis of the branch and the branch's conclusion holds (%s)%s"
                 % (nm, len(pos), "(1)" if nm == "delta34" else "(2)",
                    "one of (3)" if nm == "delta34" else "a >= J(p,b,c,q)",
                    "" if pos else "  [VOID: no such point was found]"),
                 all(r_[0] for r_ in res) and max([r_[1] for r_ in res] or [0.0]) < 1e-7,
                 "max relative |G| %.1e; details %s" % (max([r_[1] for r_ in res] or [0.0]), [r_[2] for r_ in res]))


def stratum_ok(g, rd34, rd46):
    """At a stratum point with a positive reduced non-minor coefficient: (conclusion of the branch
    holds, relative |G| in the basis of the branch, short description)."""
    import numpy as np
    from scipy.optimize import nnls
    A = np.array([F_VALUES[j] for j in GENS], float).T
    x, res = nnls(A, np.maximum(g, 0.0))
    if res > 1e-7:
        return False, 1.0, "not in the cone"
    kind, out = reduce_to_one_nonminor({j: Fraction(float(v)).limit_denominator(10 ** 13) for j, v in zip(GENS, x)})
    if kind == "minor":
        return False, 1.0, "minor"
    rd, order, names = (rd34, ORDER34, "qpjcedba") if kind == "34" else (rd46, ORDER46, "prabcqds")
    R = rd["R"]
    val = {n_: float(out[j_]) for n_, j_ in zip(R.names[:8], order)}
    val[R.names[8]] = 0.0
    relg = abs(R.evaluate(rd["G"], val)) / R.evaluate({k: abs(v) for k, v in rd["G"].items()}, val)
    co = {n_: float(out[j_]) / float(out[int(kind)]) for n_, j_ in zip(names, order)}
    if kind == "34":
        a, b, e, j, p = co["a"], co["b"], co["e"], co["j"], co["p"]
        if b * p * (a - 1) >= b + p + a + 1 - 1e-7:
            return True, relg, "four-coefficient"
        if e > 0 and j >= 1 + 1 / e - 1e-7:
            return True, relg, "old-H"
        if e > 0 and b > 0 and j >= 1 / e + 1 / b - 1e-7:
            return True, relg, "chain"
        return False, relg, "none of (3)"
    Jm, _ = rhs6_min(co["p"], co["b"], co["c"], co["q"])
    return co["a"] >= Jm - 1e-7, relg, "a/J = %.3f" % (co["a"] / Jm)


class Counter_(dict):
    def __missing__(self, key):
        return 0


if __name__ == "__main__":
    main()
