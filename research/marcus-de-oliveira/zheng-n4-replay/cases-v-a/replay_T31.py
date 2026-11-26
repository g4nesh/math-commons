"""Independent exact replay of Case T_31 = {0,1,5,6,9,12} (group V) of D. Zheng's proof of the
Marcus-de Oliveira conjecture for n = 4.

Proof document replayed: MOC4_SUPP/proofs/manin_transfer_20261021/root_moc4_T30_uniform_20261022.pdf
(4 pages; the document calls the set "T30").  Every polynomial is rebuilt here from the
definitions (the matrix P, the cubic C, the vectors f_j of the paper); nothing from the
supplement is executed or imported.  If MOC4_SUPP is set, three JSON files are read as data and
COMPARED with the rebuilt objects.

Run:   MOC4_SUPP=/path/to/supplement  python replay_T31.py
Exit status 0 iff every check passes.  Checks tagged EXACT use integer/rational/symbolic
arithmetic only; checks tagged NUMERIC are floating-point falsification attempts.
"""
import os
import random
import sys
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import replay_T29_common as M                                   # noqa: E402
from replay_T29_common import (Ring, p_add, p_mul, p_scale, p_nonneg, p_min, p_compose,   # noqa: E402
                               p_rename, p_divexact, F_VALUES)

DOC = "proofs/manin_transfer_20261021/root_moc4_T30_uniform_20261022.pdf"
RADIAL_JSON = "proofs/manin_transfer_20261021/root_moc4_T30_full_radial_20261022.json"

T = (0, 1, 5, 6, 9, 12)
ORDER = [2, 29, 1, 9, 5, 35, 0, 7]                  # g = q f2 + p f29 + j f1 + c f9 + e f5 + d f35 + b f0 + a f7
NAMES = "q,p,j,c,e,d,b,a,h"
EXCL = (15, 4, 17, 7, 11, 3, 18, 14)                # permutation indices dual to the eight generators


def mul(*xs):
    out = {0: 1}
    for x in xs:
        out = p_mul(out, x)
    return out


def main():
    import sympy as sp
    ck = M.Checker("Case T_31 = {0,1,5,6,9,12} (document 'T30'): independent exact replay")
    import numpy
    import scipy
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
    M.check_document(ck, DOC, DOC_SHA256[DOC])
    M.check_document(ck, RADIAL_JSON, DOC_SHA256[RADIAL_JSON])
    M.check_foundations(ck)
    M.check_radial_necessity(ck)
    M.crosscheck_supplement_tables(ck)
    M.crosscheck_audit_library(ck, [29, 31])

    # ------------------------------------------------------------------------------------
    ck.section("B. Contact set")
    ck.check("EXACT", "the document's set {0,1,5,6,9,12} is the paper's T_31 (Appendix E) as a set of "
             "lexicographic indices: group element = identity", tuple(sorted(T)) == M.PAPER_T_V[31])
    r31, r29 = M.orbit_rep(T), M.orbit_rep(M.PAPER_T_V[29])
    ck.check("EXACT", "orbit representative of T_31 under G is (0,1,2,6,10,15); T_29 lies in another orbit",
             r31 == (0, 1, 2, 6, 10, 15) and r29 != r31, "T_31 -> %s, T_29 -> %s" % (r31, r29))

    # ------------------------------------------------------------------------------------
    ck.section("C. Generators of the contact cone (document, Section 1, eq. (1))")
    gens = M.generators(T)
    ck.check("EXACT", "the f_j vanishing on T are exactly f2, f29, f1, f9, f5, f35, f0, f7", gens == sorted(ORDER),
             "computed %s" % gens)
    nm = [j for j in gens if not M.F_IS_MINOR[j]]
    ck.check("EXACT", "f35 = h_{4312,-} is the only non-minor generator", nm == [35],
             ", ".join("f%d = %s" % (j, M.f_name(j)) for j in ORDER))
    rkT = M.rank([M.P[s] for s in T])
    null = M.nullspace([M.P[s] for s in T], 14)
    ck.check("EXACT", "the six contact equations have rank 6 on V_R; their solution space has dimension 8",
             rkT == 6 and len(null) == 8)
    ident = [[F_VALUES[j][s] for j in ORDER] for s in EXCL]
    ck.check("EXACT", "the values of the eight generators at the permutation indices (15,4,17,7,11,3,18,14) "
             "form the identity matrix, so the eight generators are a basis of the solution space and the "
             "coefficient of each generator in c is the gap of c at the corresponding index",
             ident == [[int(i == k) for k in range(8)] for i in range(8)]
             and M.rank([M.F_COORD[j] for j in ORDER]) == 8)
    d, rays = M.extreme_rays_of_contact_cone(T)
    ck.check("EXACT", "direct enumeration of the extreme rays of N_T = {Pc >= 0, (Pc)_T = 0} (dimension %d; "
             "no use of the 82-facet list) gives exactly these eight vectors" % d,
             rays == sorted(F_VALUES[j] for j in ORDER), "%d rays" % len(rays))
    M.check_minor_leibniz(ck, [j for j in ORDER if j != 35])

    # ------------------------------------------------------------------------------------
    ck.section("D. The three sufficient coefficient conditions (2), (3), (4) and Lemma 1")
    rng = random.Random(20261126)
    M.check_transport(ck, "four-coefficient", (7, 0, 5, 35), rng)
    M.check_transport(ck, "old-H", (9, 1, 35), rng)
    M.check_transport(ck, "chain", (1, 7, 9, 35), rng)
    ck.check("EXACT", "the document's base triple (f3, f5, f54) for old-H lies in the orbit of the audit's "
             "representative (f12, f2, f70) with the same roles", bool(M.find_transport((12, 2, 70), (3, 5, 54))))
    a, b, c, e, j, t, d_ = sp.symbols("a b c e j t d", positive=True)
    # four-coefficient: coefficients (a_role, b_role, c_role, d_role) = (a, b, e, 1) on (f7, f0, f5, f35)
    four = lambda A_, B_, C_, D_: A_ * B_ * (C_ - D_) - D_ ** 2 * (A_ + B_ + C_ + D_)
    ck.check("EXACT", "(2): the four-coefficient condition for a f7 + b f0 + e f5 + 1 f35 is "
             "ab(e-1) - (a+b+e+1) >= 0", sp.expand(four(a, b, e, 1) - (a * b * (e - 1) - (a + b + e + 1))) == 0)
    ck.check("EXACT", "(3): the old-H condition t(c-d) >= d^2 for c f9 + j f1 + 1 f35 (t-role f9, c-role f1) is "
             "c(j-1) >= 1, i.e. j >= 1 + 1/c for c > 0",
             sp.simplify((c * (j - 1) - 1) - c * (j - (1 + 1 / c))) == 0)
    chain = lambda centre, leaf1, leaf2, D_: centre - D_ ** 2 / leaf1 - D_ ** 2 / leaf2
    ck.check("EXACT", "(4): the chain condition a_centre >= d^2/p + d^2/b (centre f1 with coefficient j, leaves "
             "f7, f9 with coefficients a, c, non-minor coefficient 1) is j >= 1/a + 1/c",
             sp.cancel(chain(j, a, c, 1) - (j - 1 / a - 1 / c)) == 0)
    E = t * (a + 1) * (b + t) / (a * b - t)
    J = (1 - t) / a + (1 - t) ** 2 / c
    ck.check("EXACT", "Lemma 1, four-coefficient share (t a) f7 + b f0 + e f5 + t f35: "
             "four(ta, b, e, t) = t (ab - t) (e - E(t)) with E(t) = t(a+1)(b+t)/(ab-t)",
             sp.cancel(four(t * a, b, e, t) - t * (a * b - t) * (e - E)) == 0)
    ck.check("EXACT", "Lemma 1, chain share j f1 + (1-t)a f7 + c f9 + (1-t) f35: "
             "(1-t)^2/((1-t)a) + (1-t)^2/c = J(t) = (1-t)/a + (1-t)^2/c",
             sp.cancel((1 - t) ** 2 / ((1 - t) * a) + (1 - t) ** 2 / c - J) == 0)
    ck.check("EXACT", "Lemma 1: the two shares add up to the coefficients of (1): t a + (1-t) a = a on f7 and "
             "t + (1-t) = 1 on f35; b, e are used only by the first share and c, j only by the second; "
             "q f2 and p f29 are unused minor terms", sp.expand(t * a + (1 - t) * a - a) == 0)
    ck.check("EXACT", "E'(t) = (a+1)(a b^2 + 2abt - t^2)/(ab-t)^2",
             sp.cancel(sp.diff(E, t) - (a + 1) * (a * b ** 2 + 2 * a * b * t - t ** 2) / (a * b - t) ** 2) == 0)
    ck.check("EXACT", "a b^2 + 2abt - t^2 = a b^2 + t(ab - t) + abt (so E' > 0 for 0 <= t < ab); E(0) = 0; "
             "E(1) = (a+1)(b+1)/(ab-1) = Emax",
             sp.expand(a * b ** 2 + 2 * a * b * t - t ** 2 - (a * b ** 2 + t * (a * b - t) + a * b * t)) == 0
             and E.subs(t, 0) == 0
             and sp.cancel(E.subs(t, 1) - (a + 1) * (b + 1) / (a * b - 1)) == 0)
    ck.check("EXACT", "ab(e-1) - (a+b+e+1) = e(ab-1) - (a+1)(b+1): so (2) never holds for ab <= 1, and for "
             "ab > 1 it holds iff e >= Emax, i.e. iff L_e e >= 1 with L_e = (ab-1)/((a+1)(b+1))",
             sp.expand(a * b * (e - 1) - (a + b + e + 1) - (e * (a * b - 1) - (a + 1) * (b + 1))) == 0)

    # ------------------------------------------------------------------------------------
    ck.section("E. Radial necessity in the generators (document, Section 2)")
    rd = M.radial_data(ck, ORDER, NAMES)
    R9, C0, Q, G = rd["R"], rd["C0"], rd["Q"], rd["G"]
    ck.check("EXACT", "every coefficient of C0 is negative (%d terms) and the coefficient of d^3 is -1: "
             "C0 <= -d^3 < 0 for d > 0 on the closed orthant, hence Q != 0 at a point with C0 + hQ = 0"
             % len(C0), all(v < 0 for v in C0.values()) and C0.get(R9.key((0, 0, 0, 0, 0, 3, 0, 0, 0))) == -1)
    ck.info("C0 = " + R9.fmt(C0, 30))
    ck.info("Q  = " + R9.fmt(Q, 30))
    ck.check("EXACT", "G has 1183 terms and degree 9 (the document's count for G_9)",
             len(G) == 1183 and R9.tdeg(G) == 9 and all(sum(R9.exps(k)) == 9 for k in G))
    M.compare_radial_json(ck, rd, RADIAL_JSON)
    Rs = Ring("q,p,j,c,e,b,a")
    G1 = p_rename(R9.subs_one(G, "d"), R9, Rs)
    ks = Rs.split(G1, "j")
    ck.check("EXACT", "at d = 1: G_9 = k0 + j k1 + j^2 k2 (degree 2 in j; no cancellation of terms at d = 1)",
             set(ks) == {0, 1, 2} and len(G1) == 1183,
             "terms of k0, k1, k2: %d, %d, %d" % (len(ks[0]), len(ks[1]), len(ks[2])))
    k0, k1, k2 = ks[0], ks[1], ks[2]
    ck.check("EXACT", "k2 has only positive coefficients", all(v > 0 for v in k2.values()))

    # ------------------------------------------------------------------------------------
    ck.section("F. The strict constant coefficient: k0 >= B (1 - L_e e)^2, eq. (5)")
    va, vb = Rs.var("a"), Rs.var("b")
    one = {0: 1}
    B = p_mul(Rs.pow(p_add(va, one), 3), Rs.pow(p_add(vb, one), 3))
    ck.check("EXACT", "k0 at c = e = q = p = 0 equals B = (a+1)^3 (b+1)^3",
             Rs.subs_zero(k0, ["c", "e", "q", "p"]) == B)
    ck.check("EXACT", "k0 at e = 0, minus B, has nonnegative coefficients",
             p_nonneg(p_add(Rs.subs_zero(k0, ["e"]), B, -1)))
    Rc = Ring("A,T,V,C,Q,P")
    A_, T_, V_, C_ = Rc.var("A"), Rc.var("T"), Rc.var("V"), Rc.var("C")
    ren = {"c": "C", "q": "Q", "p": "P"}
    sA, sT, sV, sC = sp.symbols("A T V C")          # the same (assumption-free) symbols as Rc.to_sympy
    one_c = {0: 1}

    def P1(x):
        return p_add(x, one_c)

    def rat(pair):
        return Rc.to_sympy(pair[0]) / Rc.to_sympy(pair[1])

    def emax_e(an, ad, bn, bd):
        """e = Emax V/(1+V) with Emax = (a+1)(b+1)/(ab-1): numerator and denominator."""
        return mul(p_add(an, ad), p_add(bn, bd), V_), p_mul(p_add(mul(an, bn), mul(ad, bd), -1), P1(V_))

    doc_counts5 = {(False, False): 2854, (True, False): 2875, (False, True): 4982, (True, True): 5047}
    for swap in (False, True):
        for high in (False, True):
            an = mul(P1(T_), P1(A_)) if high else mul(T_, P1(A_))
            ad = one_c if high else P1(T_)
            bn = P1(T_) if high else T_
            bd = P1(A_) if high else mul(P1(T_), P1(A_))
            if swap:
                an, bn, ad, bd = bn, an, bd, ad
            sa_, sb_ = rat((an, ad)), rat((bn, bd))
            z = (1 + sT) if high else sT / (1 + sT)
            okchart = sp.cancel(sa_ * sb_ - z ** 2) == 0 and \
                sp.cancel((sb_ / sa_ if swap else sa_ / sb_) - (1 + sA) ** 2) == 0
            if high:
                en, ed = emax_e(an, ad, bn, bd)
                se = rat((en, ed))
                emax = (sa_ + 1) * (sb_ + 1) / (sa_ * sb_ - 1)
                Le = (sa_ * sb_ - 1) / ((sa_ + 1) * (sb_ + 1))
                okchart &= sp.cancel(se - emax * sV / (1 + sV)) == 0 and sp.cancel(1 - Le * se - 1 / (1 + sV)) == 0
            else:
                en, ed = V_, one_c
            subs = {"a": (an, ad), "b": (bn, bd), "e": (en, ed)}
            numk, degs = p_compose(k0, Rs, Rc, subs, rename=ren)
            numB, _ = p_compose(B, Rs, Rc, subs, degs=degs, rename=ren)
            if high:
                # k0 - B/(1+V)^2: the numerator of B is divisible by (1+V)^2 (e-degree 2 of k0)
                nb = p_divexact(p_divexact(numB, P1(V_)), P1(V_))
                diff = p_add(numk, nb, -1)
            else:
                diff = p_add(numk, numB, -1)
            name = "ab %s 1, %s" % (">" if high else "<", "a <= b" if swap else "a >= b")
            ck.check("EXACT", "(5) chart %s: sqrt(ab) = %s, ratio (1+A)^2%s; the numerator of "
                     "k0 - B(1 - L_e e)^2 has nonnegative coefficients (document: %d terms)"
                     % (name, "1+T" if high else "T/(1+T)", ", e = Emax V/(1+V) so 1 - L_e e = 1/(1+V)" if high
                        else ", e = V, L_e = 0", doc_counts5[(swap, high)]),
                     okchart and p_nonneg(diff) and len(diff) == doc_counts5[(swap, high)],
                     "%d terms, min coefficient %d" % (len(diff), p_min(diff)))
    ck.info("PROOF-READ: the four charts cover {a, b > 0, ab != 1, e >= 0, L_e e < 1}; both sides of (5) are "
            "continuous in (a,b,c,e,q,p) on the closed orthant (L_e = max(.,0) is continuous), so (5) holds "
            "on the closure, i.e. for all nonnegative coefficients with L_e e <= 1; it is > 0 iff L_e e < 1, "
            "i.e. iff the four-coefficient condition (2) fails.")

    # ------------------------------------------------------------------------------------
    ck.section("G. Seven domains: Bernstein coefficients of G_9 on [0, J) (document, Section 3)")

    def common_k(subs):
        degs = {}
        for kk in (k0, k1, k2):
            for nm_ in subs:
                degs[nm_] = max(degs.get(nm_, 0), Rs.deg(kk, nm_))
        return [p_compose(kk, Rs, Rc, subs, degs=degs, rename=ren)[0] for kk in (k0, k1, k2)], degs

    def bern(K, N, D):
        mid = p_add(p_scale(p_mul(D, K[0]), 2), p_mul(N, K[1]))
        end = p_add(p_add(p_mul(mul(D, D), K[0]), p_mul(mul(N, D), K[1])), p_mul(mul(N, N), K[2]))
        return mid, end

    def nonneg_rat(expr):
        """numerator and denominator of a rational function of A, T, V, C have coefficients >= 0."""
        n_, d_ = sp.fraction(sp.cancel(expr))
        syms = (sA, sT, sV, sC)
        cn, cd = sp.Poly(n_, *syms).coeffs(), sp.Poly(d_, *syms).coeffs()
        return (all(co >= 0 for co in cn) and all(co > 0 for co in cd)) or \
            (all(co <= 0 for co in cn) and all(co < 0 for co in cd))

    def lemma1_e(an, ad, bn, bd, tn, td):
        """E(t) = t(a+1)(b+t)/(ab-t) for a = an/ad, b = bn/bd, t = tn/td, as (numerator, denominator);
        the denominator contains the factor (an bn td - tn ad bd) = numerator of ab - t."""
        return (mul(tn, p_add(an, ad), p_add(mul(bn, td), mul(tn, bd))),
                p_mul(td, p_add(mul(an, bn, td), mul(tn, ad, bd), -1)))

    sa_s, sb_s, sc_s, st_s = sp.symbols("a b c t", positive=True)
    E_s = st_s * (sa_s + 1) * (sb_s + st_s) / (sa_s * sb_s - st_s)
    J_s = (1 - st_s) / sa_s + (1 - st_s) ** 2 / sc_s
    doms = {}
    # 1: 0 <= b <= a <= 1
    doms[1] = dict(a=(A_, P1(A_)), b=(mul(A_, T_), mul(P1(A_), P1(T_))), e=(V_, one_c), kind="oldH",
                   region=lambda a_, b_: [a_ - b_, 1 - a_])
    # 2: 0 <= a <= b <= 1
    doms[2] = dict(a=(mul(A_, T_), mul(P1(A_), P1(T_))), b=(A_, P1(A_)), e=(V_, one_c), kind="oldH",
                   region=lambda a_, b_: [b_ - a_, 1 - b_])
    # 3: ab < 1, a >= 1; t = ab V/(1+V), e = E(t), J = J(t)
    an3, ad3, bn3, bd3 = P1(A_), one_c, T_, mul(P1(T_), P1(A_))
    tn3, td3 = mul(an3, bn3, V_), mul(ad3, bd3, P1(V_))
    e3 = (mul(T_, V_, P1(an3), p_add(one_c, p_mul(P1(an3), V_))), mul(P1(T_), an3, P1(V_)))
    doms[3] = dict(a=(an3, ad3), b=(bn3, bd3), t=(tn3, td3), kind="alloc", e_doc=e3,
                   region=lambda a_, b_: [a_ - 1, 1 - a_ * b_])
    # 4: ab < 1, b >= 1
    doms[4] = dict(a=(T_, mul(P1(T_), P1(A_))), b=(P1(A_), one_c), e=(V_, one_c), kind="oldH",
                   region=lambda a_, b_: [b_ - 1, 1 - a_ * b_])
    # 5: ab > 1, a <= 1
    a5, b5 = (A_, P1(A_)), (mul(P1(T_), P1(A_)), A_)
    doms[5] = dict(a=a5, b=b5, e=emax_e(a5[0], a5[1], b5[0], b5[1]), kind="oldH", emax=True,
                   region=lambda a_, b_: [1 - a_, a_ * b_ - 1])
    # 6: ab > 1, 1 <= a <= b
    a6, b6 = (p_add(P1(A_), T_), P1(A_)), (mul(P1(T_), P1(T_), P1(A_)), p_add(P1(A_), T_))
    doms[6] = dict(a=a6, b=b6, e=emax_e(a6[0], a6[1], b6[0], b6[1]), kind="chain", emax=True,
                   region=lambda a_, b_: [a_ - 1, b_ - a_, a_ * b_ - 1])
    # 7: ab > 1, a >= b; t = V/(1+V)
    an7, ad7, bn7, bd7 = mul(P1(T_), P1(A_)), one_c, P1(T_), P1(A_)
    e7 = (mul(V_, p_add(an7, ad7), p_add(mul(bn7, P1(V_)), mul(bd7, V_))),
          p_mul(P1(V_), p_add(mul(an7, bn7, P1(V_)), mul(ad7, bd7, V_), -1)))
    doms[7] = dict(a=(an7, ad7), b=(bn7, bd7), t=(V_, P1(V_)), kind="alloc", e_doc=e7,
                   region=lambda a_, b_: [a_ - b_, a_ * b_ - 1])
    doc_counts = {2: (2886, 2940), 3: (18189, 33322), 4: (2469, 2472), 5: (4217, 4221),
                  6: (15906, 18183), 7: (21350, 34392)}
    Ksave = {}
    for dn in range(1, 8):
        dm = doms[dn]
        an, ad = dm["a"]
        bn, bd = dm["b"]
        sa_, sb_ = rat(dm["a"]), rat(dm["b"])
        okreg = all(nonneg_rat(x) for x in dm["region"](sa_, sb_))
        if dm["kind"] == "alloc":
            tn, td = dm["t"]
            en, ed = dm["e_doc"]                      # the document's reduced form of E(t)
            st_ = rat((tn, td))
            # 0 <= t, t < 1 and t < ab
            okreg &= nonneg_rat(1 - st_) and nonneg_rat(sa_ * sb_ - st_) and nonneg_rat(st_)
            # the unreduced closed form of E(t) agrees with it
            en0, ed0 = lemma1_e(an, ad, bn, bd, tn, td)
            okreg &= sp.cancel(rat((en0, ed0)) - rat((en, ed))) == 0
            se = rat((en, ed))
            okreg &= sp.cancel(se - E_s.subs({sa_s: sa_, sb_s: sb_, st_s: st_})) == 0
            # J(t) = (1-t)/a + (1-t)^2/c = S (S an + C td ad ... ) in the document's scaling:
            S_ = p_add(td, tn, -1)                                 # numerator of 1 - t over td
            Nn = p_mul(S_, p_add(mul(C_, td, ad), mul(an, S_)))
            Dd = mul(an, C_, td, td)
            okreg &= sp.cancel(rat((Nn, Dd)) - J_s.subs({sa_s: sa_, sc_s: sC, st_s: st_})) == 0
            if dn == 7:
                # document scaling: N = C(1+V) + an, D = an C (1+V)^2  (here S_ = 1, td = 1+V)
                Nn, Dd = p_add(mul(C_, P1(V_)), an), mul(an, C_, P1(V_), P1(V_))
                okreg &= sp.cancel(rat((Nn, Dd)) - J_s.subs({sa_s: sa_, sc_s: sC, st_s: st_})) == 0
            if dn == 3:
                # document scaling: rD = (1+T)(1+V), rN = rD - TV; N = rN (C rD + an rN), D = an C rD^2
                rD = mul(P1(T_), P1(V_))
                rN = p_add(rD, mul(T_, V_), -1)
                Nn, Dd = p_mul(rN, p_add(mul(C_, rD), mul(an, rN))), mul(an, C_, rD, rD)
                okreg &= sp.cancel(rat((Nn, Dd)) - J_s.subs({sa_s: sa_, sc_s: sC, st_s: st_})) == 0
            dm["e"] = (en, ed)
            jtxt = "t = %s, e = E(t), J = J(t)" % ("ab V/(1+V)" if dn == 3 else "V/(1+V)")
        else:
            en, ed = dm["e"]
            se = rat((en, ed))
            if dm.get("emax"):
                okreg &= sp.cancel(se - (sa_ + 1) * (sb_ + 1) / (sa_ * sb_ - 1) * sV / (1 + sV)) == 0
            if dm["kind"] == "oldH":
                Nn, Dd = P1(C_), C_
                okreg &= sp.cancel(rat((Nn, Dd)) - (1 + 1 / sC)) == 0
                jtxt = "e = %s, J = 1 + 1/c" % ("Emax V/(1+V)" if dm.get("emax") else "V")
            else:
                Nn, Dd = p_add(mul(C_, ad), an), mul(an, C_)
                okreg &= sp.cancel(rat((Nn, Dd)) - (1 / sa_ + 1 / sC)) == 0
                jtxt = "e = Emax V/(1+V), J = 1/a + 1/c"
        K, degs = common_k({"a": dm["a"], "b": dm["b"], "e": dm["e"]})
        Ksave[dn] = (K, Nn, Dd, degs)
        ck.check("EXACT", "domain %d: chart a = %s, b = %s maps the orthant into the stated region "
                 "(differences have nonnegative numerators and positive denominators); %s as rational "
                 "functions of the chart variables" % (dn, sp.factor(sa_), sp.factor(sb_), jtxt), okreg)
        if dn == 1:
            continue
        mid, end = bern(K, Nn, Dd)
        ck.check("EXACT", "domain %d: Bernstein middle numerator 2 D k0 + N k1 has nonnegative coefficients "
                 "(document: %d terms)" % (dn, doc_counts[dn][0]),
                 p_nonneg(mid) and len(mid) == doc_counts[dn][0],
                 "%d terms, min coefficient %d" % (len(mid), p_min(mid)))
        if dn != 3:
            ck.check("EXACT", "domain %d: Bernstein endpoint numerator D^2 k0 + N D k1 + N^2 k2 has nonnegative "
                     "coefficients (document: %d terms)" % (dn, doc_counts[dn][1]),
                     p_nonneg(end) and len(end) == doc_counts[dn][1],
                     "%d terms, min coefficient %d" % (len(end), p_min(end)))
        else:
            ngroups, bad = M.grouped_univariate(Rc, end, "C")
            okg = all(M.uni_positive(cf) for _, cf in bad)
            oks = all(M.uni_sympy_nonneg(cf) for _, cf in bad)
            ck.check("EXACT", "domain 3: the endpoint numerator (%d terms; document: 33322) is NOT coefficientwise "
                     "nonnegative; grouped by the monomial in (A,T,V,Q,P) it has %d coefficient polynomials in C, "
                     "of which %d have a negative coefficient (document: 21); each of these, divided by its "
                     "lowest power of C, has a positive constant term and no root in (0, inf) (own Sturm "
                     "sequences over Q; cross-checked with sympy factor_list/count_roots), so the numerator "
                     "is >= 0 for C >= 0" % (len(end), ngroups, len(bad)),
                     len(end) == doc_counts[3][1] and len(bad) == 21 and okg and oks,
                     "degrees in C of the 21 polynomials: %s" % sorted(max(cf) - min(cf) for _, cf in bad))
    ck.info("PROOF-READ: with J = N/D, D > 0 and s = j/J, D^2 Delta G_9(j) = D^2 k0~ (1-s)^2 + D (2 D k0~ + N k1~)"
            " s(1-s) + (D^2 k0~ + N D k1~ + N^2 k2~) s^2 (k~ = numerators over the common denominator Delta > 0);"
            " so nonnegative middle and endpoint numerators give G_9(j) >= k0 (1 - j/J)^2 > 0 on [0, J), eq. (6).")
    x0, x1, x2, s_, Jv_, l_ = sp.symbols("x0 x1 x2 s J l")
    quad = lambda jj: x0 + x1 * jj + x2 * jj ** 2
    id1 = sp.expand(quad(s_ * Jv_) - (x0 * (1 - s_) ** 2 + (x0 + Jv_ * x1 / 2) * 2 * s_ * (1 - s_)
                                      + quad(Jv_) * s_ ** 2))
    id2 = sp.expand(quad(l_ + s_ * (Jv_ - l_)) - (quad(l_) * (1 - s_) ** 2
                                                  + (x0 + x1 * (l_ + Jv_) / 2 + x2 * l_ * Jv_) * 2 * s_ * (1 - s_)
                                                  + quad(Jv_) * s_ ** 2))
    ck.check("EXACT", "Bernstein identities used: on [0,J], k0 + k1 j + k2 j^2 = k0 (1-s)^2 + (k0 + J k1/2) 2s(1-s)"
             " + G(J) s^2 at j = sJ; on [l,J] the middle coefficient is k0 + k1 (l+J)/2 + k2 l J",
             id1 == 0 and id2 == 0)

    # ---- domain 1: the cut l = a - b
    K, Nn, Dd, degs1 = Ksave[1]
    an, ad = doms[1]["a"]
    bn, bd = doms[1]["b"]
    ln_, ld_ = A_, mul(P1(A_), P1(T_))                           # l = a - b = A/((1+A)(1+T))
    ck.check("EXACT", "domain 1: l = a - b = A/((1+A)(1+T)), 0 <= l < 1 < J = 1 + 1/c",
             sp.cancel(rat((ln_, ld_)) - (rat(doms[1]["a"]) - rat(doms[1]["b"]))) == 0
             and nonneg_rat(1 - rat((ln_, ld_))))
    sqr = Rs.pow(p_add(p_add(va, vb, -1), Rs.var("j"), -1), 2)
    grp = Rs.group(G1, ["a", "b", "j"])
    ck.check("EXACT", "the document's remark behind the cut: the coefficients of q c^2 e^2, q^2 c e^2 and q^2 c^2 e in "
             "G_9 (as polynomials in a, b, j) are exactly (a - b - j)^2",
             all(grp.get(Rs.key(ex)) == sqr for ex in ((1, 0, 0, 2, 2, 0, 0), (2, 0, 0, 1, 2, 0, 0), (2, 0, 0, 2, 1, 0, 0))))
    whole_mid, _ = bern(K, Nn, Dd)
    ck.info("domain 1: on the whole interval [0, J] the middle numerator has a negative coefficient "
            "(min %d), so a cut is needed, as the document says" % p_min(whole_mid))
    leftmid = p_add(p_scale(p_mul(ld_, K[0]), 2), p_mul(ln_, K[1]))
    cut = p_add(p_add(p_mul(mul(ld_, ld_), K[0]), p_mul(mul(ln_, ld_), K[1])), p_mul(mul(ln_, ln_), K[2]))
    rightmid = p_add(p_add(p_scale(p_mul(mul(ld_, Dd), K[0]), 2),
                           p_mul(p_add(mul(ln_, Dd), mul(Nn, ld_)), K[1])), p_scale(p_mul(mul(ln_, Nn), K[2]), 2))
    endp = p_add(p_add(p_mul(mul(Dd, Dd), K[0]), p_mul(mul(Nn, Dd), K[1])), p_mul(mul(Nn, Nn), K[2]))
    NB, _ = p_compose(B, Rs, Rc, {"a": doms[1]["a"], "b": doms[1]["b"], "e": doms[1]["e"]}, degs=degs1, rename=ren)
    cutB = p_add(cut, p_mul(mul(ld_, ld_), NB), -1)
    ck.check("EXACT", "domain 1, own (reduced) scaling: the numerators of k0, of the middle coefficient on [0,l], "
             "of G_9(l), of the middle coefficient on [l,J], of G_9(J), and of G_9(l) - B all have "
             "nonnegative coefficients",
             all(p_nonneg(x) for x in (K[0], leftmid, cut, rightmid, endp, cutB)),
             "terms: %s" % [len(x) for x in (K[0], leftmid, cut, rightmid, endp, cutB)])
    # the document's scaling (read from its verifier's text): D = C ad bd, L = C (an bd - bn ad), H = (C+1) ad bd
    Dv = mul(C_, ad, bd)
    Lv = p_mul(C_, p_add(mul(an, bd), mul(bn, ad), -1))
    Hv = mul(P1(C_), ad, bd)
    G0 = p_mul(mul(Dv, Dv), K[0])
    GL = p_add(p_add(G0, p_mul(mul(Lv, Dv), K[1])), p_mul(mul(Lv, Lv), K[2]))
    d_left = p_add(p_scale(G0, 2), p_mul(mul(Lv, Dv), K[1]))
    d_right = p_add(p_scale(GL, 2), p_mul(p_add(Hv, Lv, -1), p_add(p_mul(Dv, K[1]), p_scale(p_mul(Lv, K[2]), 2))))
    d_end = p_add(p_add(G0, p_mul(mul(Hv, Dv), K[1])), p_mul(mul(Hv, Hv), K[2]))
    cut_src = p_add(p_add(k0, p_mul(p_add(va, vb, -1), k1)), p_mul(Rs.pow(p_add(va, vb, -1), 2), k2))
    d_cutB, _ = p_compose(p_add(cut_src, B, -1), Rs, Rc,
                          {"a": doms[1]["a"], "b": doms[1]["b"], "e": doms[1]["e"]}, rename=ren)
    cnt = [len(x) for x in (d_left, GL, d_right, d_end, d_cutB)]
    ck.check("EXACT", "domain 1, the document's scaling: the four numerators and G_9(a-b) - B have nonnegative "
             "coefficients with the document's term counts 4913, 4913, 6696, 6921 and 3331",
             all(p_nonneg(x) for x in (d_left, GL, d_right, d_end, d_cutB)) and cnt == [4913, 4913, 6696, 6921, 3331],
             "terms: %s" % cnt)
    ck.info("PROOF-READ: on [0,l]: G_9 >= k0 (1-s)^2 > 0 for s < 1 and G_9(l) >= B > 0 (eq. (7)); on [l,J): "
            "G_9 >= G_9(l)(1-s)^2 > 0.  If l = 0 only the second interval is used.")

    # ------------------------------------------------------------------------------------
    ck.section("H. Closed layers (document, Section 4)")
    Ga0 = Rs.subs_zero(G1, ["a"])
    ck.check("EXACT", "a = 0: every coefficient of G_9 (as a polynomial in b,c,e,q,p,j) is nonnegative and the "
             "constant term is 1, so G_9 >= 1 > 0 at a = 0: no spectrally realizable point has a = 0 < d",
             p_nonneg(Ga0) and Ga0.get(0) == 1, "%d terms" % len(Ga0))
    ck.check("EXACT", "c = 0: k1 and k2 have nonnegative coefficients at c = 0, so G_9 >= k0 > 0 there by (5) "
             "unless the four-coefficient condition holds",
             p_nonneg(Rs.subs_zero(k1, ["c"])) and p_nonneg(Rs.subs_zero(k2, ["c"])))
    # b = 0, a >= 1: the document passes to the limit b -> 0 in domain 3.  Direct check of the limit statement:
    Kb, _ = common_k({"a": (P1(A_), one_c), "b": ({}, one_c), "e": (V_, one_c)})
    Nb_, Db_ = p_add(C_, P1(A_)), mul(P1(A_), C_)                 # J = 1/a + 1/c at a = 1 + A
    midb, endb = bern(Kb, Nb_, Db_)
    ngb, badb = M.grouped_univariate(Rc, endb, "C")
    ck.check("EXACT", "b = 0, a = 1 + A >= 1, e = V >= 0 arbitrary, J = 1/a + 1/c (own direct certificate for the "
             "layer that the document reaches by the limit b -> 0 in domain 3): the middle numerator has "
             "nonnegative coefficients; the endpoint numerator is nonnegative for C >= 0 after grouping by the "
             "monomial in (A,V,Q,P): %d of its %d coefficient polynomials in C have a negative coefficient and "
             "each is >= 0 on [0, inf) (Sturm)" % (len(badb), ngb),
             p_nonneg(midb) and all(M.uni_nonneg(cf) and M.uni_sympy_nonneg(cf) for _, cf in badb),
             "terms %d, %d; the exceptional polynomials: %s" % (len(midb), len(endb), [cf for _, cf in badb]))
    ck.info("PROOF-READ: boundaries between domains (ab = 1, a = b, a = 1, b = 1) are in the closures of the "
            "charts; the Bernstein numerators with J = 1 + 1/c or 1/a + 1/c are polynomials in the original "
            "coefficients, so they stay >= 0 on closures, and strictness comes from (5), (7).  ab = 1 with "
            "a >= b is the face T = 0 of the chart of domain 7.")

    # ------------------------------------------------------------------------------------
    numeric_part(ck, Rs, G1, k0, k1, k2, B, rd)
    ck.finish()


def numeric_part(ck, Rs, G1, k0, k1, k2, B, rd):
    """Floating-point / random-sample attempts to falsify the intermediate claims."""
    import numpy as np
    ck.section("I. NUMERIC falsification attempts (evidence only)")
    # (i) the transported universal inequalities at their thresholds
    rng = random.Random(31001)
    worst = {}
    f = lambda idx: np.array(F_VALUES[idx], float)
    for trial in range(12):
        a, b, c = [math_exp(rng) for _ in range(3)]
        e_thr = None
        if a * b > 1:
            e_thr = (a + 1) * (b + 1) / (a * b - 1)
            v = a * f(7) + b * f(0) + e_thr * f(5) + f(35)
            worst["four-coefficient"] = min(worst.get("four-coefficient", 1e9),
                                            M.numeric_min_phi(v, 5000 + trial, 4000, 6))
        v = c * f(9) + (1 + 1 / c) * f(1) + f(35)
        worst["old-H"] = min(worst.get("old-H", 1e9), M.numeric_min_phi(v, 6000 + trial, 4000, 6))
        v = (1 / a + 1 / c) * f(1) + a * f(7) + c * f(9) + f(35)
        worst["chain"] = min(worst.get("chain", 1e9), M.numeric_min_phi(v, 7000 + trial, 4000, 6))
        tmax = min(a * b, 1.0)
        t = rng.random() * tmax * 0.98
        Et = t * (a + 1) * (b + t) / (a * b - t)
        Jt = (1 - t) / a + (1 - t) ** 2 / c
        v = a * f(7) + b * f(0) + Et * f(5) + Jt * f(1) + c * f(9) + f(35)
        worst["Lemma 1"] = min(worst.get("Lemma 1", 1e9), M.numeric_min_phi(v, 8000 + trial, 4000, 6))
    for nm, w in worst.items():
        ck.check("NUMERIC", "%s at its threshold (coefficients exp(N(0,1)), 12 draws, seeds 5000.., 4000 random "
                 "rotations + 6 BFGS descents each): min Phi found >= -1e-9" % nm, w > -1e-9, "min = %.3e" % w)
    ck.info("controls (the search does detect negative values): min Phi found for f35 alone %.4f, for "
            "f9 + f1 + f35 (j = 1 < 1 + 1/c = 2) %.4f"
            % (M.numeric_min_phi(f(35), 98, 4000, 6), M.numeric_min_phi(f(9) + f(1) + f(35), 99, 4000, 6)))

    # (ii) exact rational sampling of the polynomial claims, including boundary layers
    rng = random.Random(31002)
    vals = [Fraction(0), Fraction(0), Fraction(1), Fraction(1, 2), Fraction(2), Fraction(1, 7), Fraction(9),
            Fraction(1, 30), Fraction(40)]

    def rnd():
        r = rng.random()
        if r < 0.45:
            return rng.choice(vals)
        return Fraction(rng.randint(0, 60), rng.randint(1, 12))

    def ev(poly, pt):
        return Rs.evaluate(poly, pt)

    n5 = bad5 = nG = badG = nlay = badlay = 0
    per_dom = {i: 0 for i in range(1, 8)}
    minratio = None
    for _ in range(6000):
        a, b, c, q, p = rnd(), rnd(), rnd(), rnd(), rnd()
        ab = a * b
        # e: either free or through the allocation parameter
        if ab > 1 and a >= b and rng.random() < 0.8:
            t = Fraction(rng.randint(0, 99), 100)
            e = t * (a + 1) * (b + t) / (ab - t)
        elif ab < 1 and a >= 1 and rng.random() < 0.8 and ab > 0:
            t = ab * Fraction(rng.randint(0, 99), 100)
            e = t * (a + 1) * (b + t) / (ab - t)
        else:
            t = None
            e = rnd()
        pt = dict(a=a, b=b, c=c, e=e, q=q, p=p, j=Fraction(0))
        Bv = (a + 1) ** 3 * (b + 1) ** 3
        Le = max((ab - 1) / ((a + 1) * (b + 1)), Fraction(0))
        four_holds = ab > 1 and e * (ab - 1) >= (a + 1) * (b + 1)
        k0v = ev(k0, pt)
        if Le * e <= 1:
            n5 += 1
            if k0v < Bv * (1 - Le * e) ** 2:
                bad5 += 1
        if a == 0:
            pt["j"] = rnd()
            nlay += 1
            badlay += ev(G1, pt) < 1
            continue
        if four_holds:
            continue
        if c == 0:
            pt["j"] = rnd()
            nlay += 1
            badlay += not (ev(G1, pt) >= k0v > 0)
            continue
        # threshold of the document's domain
        if ab > 1:
            if a <= 1:
                dom, Jv = 5, 1 + 1 / c
            elif a <= b:
                dom, Jv = 6, 1 / a + 1 / c
            else:
                dom = 7
                if t is None:
                    continue
                Jv = (1 - t) / a + (1 - t) ** 2 / c
        elif ab < 1:
            if a >= 1:
                dom = 3
                if t is None:
                    if b == 0:
                        Jv = 1 / a + 1 / c
                    else:
                        continue
                else:
                    Jv = (1 - t) / a + (1 - t) ** 2 / c
            elif b >= 1:
                dom, Jv = 4, 1 + 1 / c
            else:
                dom, Jv = (1 if b <= a else 2), 1 + 1 / c
        else:
            dom, Jv = 0, 1 + 1 / c if a <= 1 else None
            if Jv is None:
                continue
        s = Fraction(rng.randint(0, 999), 1000)
        jv = s * Jv
        pt["j"] = jv
        gv = ev(G1, pt)
        nG += 1
        if dom in per_dom:
            per_dom[dom] += 1
        bound = k0v * (1 - s) ** 2
        if gv <= 0 or (dom != 1 and gv < bound):
            badG += 1
        if bound > 0:
            r = gv / bound
            minratio = r if minratio is None or r < minratio else minratio
    ck.check("NUMERIC", "(5) on %d random rational points with L_e e <= 1 (seed 31002; zeros and boundary values "
             "included; exact evaluation): no violation of k0 >= B(1 - L_e e)^2" % n5, bad5 == 0)
    ck.check("NUMERIC", "layers a = 0 (G_9 >= 1) and c = 0 with (2) failing (G_9 >= k0 > 0), j >= 0 arbitrary: "
             "%d random rational points, exact evaluation, no violation" % nlay, badlay == 0 and nlay > 500)
    ck.check("NUMERIC", "G_9 > 0 (and (6) G_9 >= k0 (1 - j/J)^2 outside domain 1) at %d random rational points "
             "with (2) failing, a, c > 0 and j = sJ, s in [0,1), J the document's threshold of the domain "
             "(points per domain %s; exact evaluation)" % (nG, per_dom), badG == 0 and nG > 1500,
             "min of G_9 / (k0 (1-s)^2) = %.4f" % float(minratio))

    # (iii) genuine points of X: the necessary equations, and the conclusion on the stratum
    R9, G = rd["R"], rd["G"]
    eq = M.numeric_equality_points(T, 31005, 80)
    fit, cub, hes, qmin, gmax = M.numeric_radial_test(rd, ORDER, eq)
    ck.check("NUMERIC", "end-to-end test of the necessary equations on %d generic points x of X with "
             "Re(F_sigma - F_pi*) = 0 on T (seed 31005; eight labels apart by > 0.05 in normalised distance; the "
             "other gaps have either sign, so the coordinates q,...,a have mixed signs): relative "
             "|C0 + H Q| < 1e-9 and relative |G| < 1e-9" % len(eq),
             len(eq) >= 40 and fit < 1e-9 and cub < 1e-9 and gmax < 1e-9,
             "fit %.1e, cubic %.1e, |QG| %.1e, |G| %.1e, smallest relative |Q| %.1e" % (fit, cub, hes, gmax, qmin))
    absG = {k: abs(v) for k, v in G.items()}
    pts = []
    for k, eps in enumerate((0.02, 0.1, 0.3, 0.005)):
        pts += M.numeric_stratum_search(T, 31100 + k, 300, eps=eps)
    dmax = max([float(p_["gaps"][EXCL[5]]) for p_ in pts] or [float("nan")])
    cover = dpos = 0
    seps, ncont = [], []
    rel = [0.0]
    which = {"(2) four-coefficient": 0, "(3) old-H": 0, "(4) chain": 0, "Lemma 1": 0}
    for p_ in pts:
        g = p_["gaps"]
        co = dict(zip(R9.names[:8], [float(g[s]) for s in EXCL]))
        if co["d"] <= 1e-4:
            continue
        dpos += 1
        seps.append(p_["sep"])
        ncont.append(sum(1 for x_ in g if abs(x_) < 1e-7))
        co["h"] = 0.0
        rel.append(abs(R9.evaluate(G, co)) / R9.evaluate(absG, {k_: abs(v) for k_, v in co.items()}))
        a, b, c, e, j = [co[x] / co["d"] for x in "abcej"]
        c2 = a * b > 1 and e * (a * b - 1) >= (a + 1) * (b + 1) - 1e-7
        c3 = c > 0 and j >= 1 + 1 / c - 1e-7
        c4 = a > 0 and c > 0 and j >= 1 / a + 1 / c - 1e-7
        which["(2) four-coefficient"] += c2
        which["(3) old-H"] += c3
        which["(4) chain"] += c4
        ok = c2 or c3 or c4
        if not ok and a > 0 and b > 0 and c > 0:
            for tt in np.linspace(0, min(a * b, 1.0), 2001)[:-1]:
                if e >= tt * (a + 1) * (b + tt) / (a * b - tt) - 1e-7 and j >= (1 - tt) / a + (1 - tt) ** 2 / c - 1e-7:
                    ok = True
                    which["Lemma 1"] += 1
                    break
        cover += ok
    ck.info("stratum (gaps >= 0, zero on T, sum 1): %d points from two-stage LM searches (seeds 31100..31103, 300 "
            "starts each, repulsion 0.02, 0.1, 0.3, 0.005; residual < 1e-10); %d of them have d > 1e-4; "
            "largest d = %.2e; among these %d points the smallest normalised distance of two labels is %.2e and "
            "the numbers of contacts (|gap| < 1e-7) are %s"
            % (len(pts), dpos, dmax, dpos, min(seps or [float("nan")]), sorted(ncont)))
    ck.check("NUMERIC", "at every stratum point found with d > 1e-4 (%d points): relative |G_9| < 1e-7 and one of "
             "(2), (3), (4), Lemma 1 holds%s" % (dpos, "" if dpos else "  [VOID: no such point was found]"),
             cover == dpos and max(rel) < 1e-7, "max relative |G| %.1e; number of points at which each condition holds (Lemma 1 is tried only "
             "if the other three fail): %s" % (max(rel), which))


def math_exp(rng):
    import math
    return math.exp(rng.gauss(0.0, 1.0))


DOC_SHA256 = {
    DOC: "6b69d7659f15cbe53ecc73b215f0122da6cb51d9e288197090e2cc8534fd0043",
    RADIAL_JSON: "b95de3e71b5927ee21bebfa1a5d4826151abaec7509a2becc175ef7d877f88ad",
}

if __name__ == "__main__":
    main()
