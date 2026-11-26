#!/usr/bin/env python3
"""Independent exact replay of the proof of contact case T_1 = {0,1,2,3,6,7}.

Proof replayed: supplement, proofs/manin_transfer_20261021/independent_moc4_six_contact_reduction_20261022.pdf,
Section 5.1 (pp. 6-7).  Every check is exact (integers, rationals, symbolic); nothing from the supplement is
executed or imported.  The supplement is read only in the last section, as data, to compare lists, and only
when the environment variable MOC4_SUPP points to it.

Run:  MOC4_SUPP=/path/to/supplement python replay_T1.py        exit code 0 iff every check passes.
Packages: sympy (plus the standard library).  Shared objects: t12_common.py in this directory.

Notation of the document: A_i = a_i, C_j = -b_j, F_sigma = lam * prod_i (A_i - C_sigma(i)), D_rs = F_r - F_s,
permutations indexed 0..23 in lexicographic order.
"""
import os
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sympy as sp

import t12_common as c

T = frozenset({0, 1, 2, 3, 6, 7})        # the document's contact set, p. 6


def main():
    ck = c.Checker("Exact replay, contact case T_1 = {0,1,2,3,6,7} "
                   "(independent_moc4_six_contact_reduction_20261022.pdf, Section 5.1, pp. 6-7)")
    c.check_basic_objects(ck)

    # ------------------------------------------------------------------------------------------------
    ck.section("1. Identification of the contact set (task item 1)")
    ck.check("document's set equals the paper's T_1 (Appendix E) as a set of lexicographic indices",
             T == c.PAPER_T[1])
    ck.info("permutations of T_1 (0-based one-line): %s" % ", ".join("%d:%s" % (k, c.oneline(c.PERMS[k]))
                                                                    for k in sorted(T)))
    ident = c.G_ELEMENTS.index((c.IDENT, c.IDENT, 0))
    ck.check("the group element relating them is the identity (rho = rho' = id, no transposition)",
             c.g_image(ident, T) == c.PAPER_T[1])
    orb = c.orbit_of(T)
    ck.check("orbit of T_1 under the group of order 1152 has 576 elements; T_1 is its lexicographic minimum",
             len(orb) == 576 and c.orbit_min(T) == tuple(sorted(T)), "orbit size %d" % len(orb))
    stab = [k for k in range(1152) if c.g_image(k, T) == T]
    ck.info("stabiliser of T_1: %d elements: %s" % (len(stab), ", ".join(
        "(rho=%s, rho'=%s, eps=%d)" % (c.oneline(c.G_ELEMENTS[k][0]), c.oneline(c.G_ELEMENTS[k][1]),
                                       c.G_ELEMENTS[k][2]) for k in stab)))

    # ------------------------------------------------------------------------------------------------
    ck.section("2. Generators of the contact cone N_{T_1} (task item 2)")
    gens = c.generators(T)
    nonminor = [j for j in gens if not c.F_IS_MINOR[j]]
    ck.check("f_j vanishing on T_1: j in {0,1,2,3,4,5,6,7,8,15,16,17,18}",
             gens == [0, 1, 2, 3, 4, 5, 6, 7, 8, 15, 16, 17, 18], str(gens))
    ck.check("exactly two of them are not minor vectors: f_17 = h_{3421,+} and f_18 = h_{3412,-}",
             nonminor == [17, 18] and c.F_NAME[17] == "h3421+" and c.F_NAME[18] == "h3412-")
    ck.info("generators: %s" % ", ".join("f%d=%s" % (j, c.F_NAME[j]) for j in gens))
    d, rays = c.contact_cone_extreme_rays(T)
    ck.check("direct extreme-ray enumeration of N_{T_1} (dimension %d) returns exactly these 13 profiles "
             "(no use of the 82-facet theorem)" % d, set(rays) == set(c.F_PROF[j] for j in gens))

    c.check_entry_faces(ck)

    # ------------------------------------------------------------------------------------------------
    ck.section("3. Spectral identities of the proof (task item 3)")
    A = sp.symbols("A0:4")
    C = sp.symbols("C0:4")
    lam = sp.Symbol("lam")
    F = c.F_products(A, C, lam)

    def D(r, s):
        return F[r] - F[s]

    doc = {
        (0, 1): lam * (A[0] - C[0]) * (A[1] - C[1]) * (A[2] - A[3]) * (C[2] - C[3]),
        (2, 3): lam * (A[0] - C[0]) * (A[1] - C[2]) * (A[2] - A[3]) * (C[1] - C[3]),
        (0, 2): lam * (A[0] - C[0]) * (A[3] - C[3]) * (A[1] - A[2]) * (C[1] - C[2]),
        (1, 3): lam * (A[0] - C[0]) * (A[2] - C[3]) * (A[1] - A[3]) * (C[1] - C[2]),
        (0, 6): lam * (A[2] - C[2]) * (A[3] - C[3]) * (A[0] - A[1]) * (C[0] - C[1]),
        (1, 7): lam * (A[2] - C[3]) * (A[3] - C[2]) * (A[0] - A[1]) * (C[0] - C[1]),
    }
    ck.check("the six differences D_01, D_23, D_02, D_13, D_06, D_17 factor as products of lam, two factors of "
             "contact products and two within-spectrum differences (expanded identities)",
             all(sp.expand(D(r, s) - v) == 0 for (r, s), v in doc.items()))
    ck.check("all six pairs have both endpoints in T_1, so each difference has real part 0 on the contact stratum",
             all(r in T and s in T for (r, s) in doc))
    # which hypotheses make each difference nonzero
    cross_in_T = set((i, c.PERMS[k][i]) for k in T for i in range(4))
    ok, usedA, usedC = True, set(), set()
    for (r, s) in doc:
        cross, wA, wC, has_lam, left, coeff = c.classify_factors(D(r, s), A, C, lam)
        ok = ok and has_lam and not left and abs(coeff) == 1 and cross <= cross_in_T \
            and len(cross) == 2 and len(wA) == 1 and len(wC) == 1
        usedA |= wA
        usedC |= wC
    ck.check("own factorisation: every cross factor A_i - C_j of the six differences divides a product F_sigma with "
             "sigma in T_1; hence all six are nonzero when lam != 0, the six contact products are nonzero and the "
             "points within each spectrum are distinct", ok)
    ck.info("within-spectrum differences actually used: A-pairs %s, C-pairs %s" % (sorted(usedA), sorted(usedC)))
    ck.info("cross pairs (i,j) with A_i - C_j dividing a contact product: %s" % sorted(cross_in_T))
    ck.info("cross pairs occurring in no contact product (collision allowed in the main branch): %s"
            % sorted(set((i, j) for i in range(4) for j in range(4)) - cross_in_T))

    r1n, r1d = (A[2] - C[2]) * (A[3] - C[3]), (A[2] - C[3]) * (A[3] - C[2])
    r2n, r2d = (A[3] - C[3]) * (A[1] - A[2]), (A[2] - C[3]) * (A[1] - A[3])
    r3n, r3d = (A[1] - C[1]) * (C[2] - C[3]), (A[1] - C[2]) * (C[1] - C[3])
    ck.check("ratio 1: D_06 / D_17 = (A2-C2)(A3-C3) / ((A2-C3)(A3-C2))  (cross-multiplied identity)",
             sp.expand(D(0, 6) * r1d - D(1, 7) * r1n) == 0)
    ck.check("ratio 2: D_02 / D_13 = (A3-C3)(A1-A2) / ((A2-C3)(A1-A3))",
             sp.expand(D(0, 2) * r2d - D(1, 3) * r2n) == 0)
    ck.check("ratio 3: D_01 / D_23 = (A1-C1)(C2-C3) / ((A1-C2)(C1-C3))",
             sp.expand(D(0, 1) * r3d - D(2, 3) * r3n) == 0)
    ck.check("the four points A2, A3, C2, C3 are pairwise distinct under the hypotheses: the four cross pairs "
             "(2,2),(3,3),(2,3),(3,2) divide contact products",
             {(2, 2), (3, 3), (2, 3), (3, 2)} <= cross_in_T)

    # real cross-ratio <=> concyclic (determinant criterion), as an exact identity
    xs = sp.symbols("x1:5", real=True)
    ys = sp.symbols("y1:5", real=True)
    z = [xs[k] + sp.I * ys[k] for k in range(4)]
    num = (z[0] - z[2]) * (z[1] - z[3])
    den = (z[0] - z[3]) * (z[1] - z[2])
    imcr = sp.im(sp.expand(num * sp.conjugate(den)))
    circ = sp.Matrix([[xs[k] ** 2 + ys[k] ** 2, xs[k], ys[k], 1] for k in range(4)]).det()
    ck.check("Im[(z1-z3)(z2-z4) conj((z1-z4)(z2-z3))] = det[|z_k|^2, Re z_k, Im z_k, 1] (rows k = 1..4): a "
             "cross-ratio with nonzero denominator is real iff the four points lie on a circle or a line",
             sp.expand(imcr - circ) == 0)

    # ------------------------------------------------------------------------------------------------
    ck.section("4. The projective change of coordinate (task item 4: hypotheses of the normalisation)")
    al, be, ga, de, a_, c_ = sp.symbols("alpha beta gamma delta a c")
    ck.check("two-point identity: (alpha a + beta)(gamma c + delta) - (alpha c + beta)(gamma a + delta) "
             "= (alpha delta - beta gamma)(a - c)",
             sp.expand((al * a_ + be) * (ga * c_ + de) - (al * c_ + be) * (ga * a_ + de)
                       - (al * de - be * ga) * (a_ - c_)) == 0)
    # exact instance: all 24 products scale by one common factor
    Aq = [sp.Rational(3, 2) + sp.I, sp.Rational(-1, 3) + 2 * sp.I, 2 - sp.I / 2, sp.Rational(5, 7) + 3 * sp.I]
    Cq = [-1 + sp.I / 3, sp.Rational(1, 4) - sp.I, sp.Rational(-7, 5) - 2 * sp.I, 3 + sp.I / 5]
    mob = (sp.Rational(2) + sp.I, sp.Rational(1, 3) - sp.I, sp.Rational(1, 2) + sp.I / 4, 1 - 2 * sp.I)

    def M(w):
        return (mob[0] * w + mob[1]) / (mob[2] * w + mob[3])

    F0 = [sp.nsimplify(sp.expand(v)) for v in c.F_products(Aq, Cq, 1)]
    F1 = c.F_products([M(w) for w in Aq], [M(w) for w in Cq], 1)
    common = (mob[0] * mob[3] - mob[1] * mob[2]) ** 4
    for w in Aq + Cq:
        common = common / (mob[2] * w + mob[3])
    ck.check("exact Gaussian-rational instance: after a Moebius map all 24 products are multiplied by the same "
             "number (alpha delta - beta gamma)^4 / prod(gamma w + delta), w over the eight points",
             all(sp.simplify(F1[k] - common * F0[k]) == 0 for k in range(24)) and sp.simplify(common) != 0)
    ck.info("consequence (read, not computed): replacing lam by lam/common leaves all 24 numbers F_sigma unchanged, "
            "so the gap vector, the contact set and the Leibniz sum are unchanged; a circle through four of the "
            "points has infinitely many points, so the point sent to infinity can avoid all eight labelled points.")

    # ------------------------------------------------------------------------------------------------
    ck.section("5. The normalised coordinates A2 = x, A3 = y, C2 = t, C3 = 0 (x, y, t real)")
    x, y, t, u, v, p, q, l, m = sp.symbols("x y t u v p q l m", real=True)
    alpha, beta, LL = u + sp.I * v, p + sp.I * q, l + sp.I * m
    sub2 = {A[1]: alpha, A[2]: x, A[3]: y, C[3]: 0}
    n1, d1 = r2n.subs(sub2), r2d.subs(sub2)
    ck.check("ratio 2 becomes y(A1 - x) / (x(A1 - y))", sp.expand(n1 - y * (alpha - x)) == 0
             and sp.expand(d1 - x * (alpha - y)) == 0)
    ck.check("its imaginary cross-multiplied numerator is x y (x - y) Im(A1)",
             sp.expand(sp.im(sp.expand(n1 * sp.conjugate(d1))) - x * y * (x - y) * v) == 0)
    ck.check("x, y, t != 0, x != y and the denominator of ratio 2 are nonzero under the hypotheses "
             "(pairs (2,3), (3,3) divide contact products; the rest are within-spectrum differences)",
             {(2, 3), (3, 3)} <= cross_in_T)
    sub3 = {A[1]: u, C[1]: beta, C[2]: t, C[3]: 0}
    n2, d2 = r3n.subs(sub3), r3d.subs(sub3)
    ck.check("ratio 3 becomes t(u - C1) / ((u - t) C1) once A1 = u is real",
             sp.expand(n2 - t * (u - beta)) == 0 and sp.expand(d2 - (u - t) * beta) == 0)
    ck.check("its imaginary cross-multiplied numerator is -t u (u - t) Im(C1)",
             sp.expand(sp.im(sp.expand(n2 * sp.conjugate(d2))) + t * u * (u - t) * q) == 0)
    ck.check("u - t = A1 - C2 is nonzero under the hypotheses (pair (1,2) divides a contact product), "
             "so the numerator vanishes iff u Im(C1) = 0", (1, 2) in cross_in_T)
    ck.check("at u = 0 ratio 3 is identically 1, so it carries no information: the collision A1 = C3 needs "
             "the separate argument", sp.expand((n2 - d2).subs(u, 0)) == 0)
    ck.check("the collision A1 = C3 is compatible with the hypotheses: the pair (1,3) divides no contact product",
             (1, 3) not in cross_in_T)

    ck.section("6. The exceptional collision A1 = C3 (= 0 in the normalised coordinate)")
    special = {A[1]: 0, A[2]: x, A[3]: y, C[1]: beta, C[2]: t, C[3]: 0}
    Lsym = lam * (A[0] - C[0])
    ck.check("D_01 = -L C1 t (x - y) with L = lam (A0 - C0)",
             sp.expand(D(0, 1).subs(special) + Lsym * beta * t * (x - y)) == 0)
    ck.check("D_02 = -L x y (C1 - t)", sp.expand(D(0, 2).subs(special) + Lsym * x * y * (beta - t)) == 0)
    e1 = sp.re(sp.expand(LL * beta))
    e2 = sp.re(sp.expand(LL * (beta - t)))
    ck.check("Re(L C1) - Re(L (C1 - t)) = t Re(L)", sp.expand(e1 - e2 - t * l) == 0)
    ck.check("if Re L = 0 then Re(L C1) = -Im(L) Im(C1)", sp.expand(e1.subs(l, 0) + m * q) == 0)
    ck.info("read: Re D_01 = Re D_02 = 0 with t (x - y) != 0 and x y != 0 give Re(L C1) = Re(L (C1 - t)) = 0, "
            "hence Re L = 0 (t != 0), hence Im C1 = 0 (L != 0).")

    ck.section("7. Conclusion of the main branch: six real points force the entry face sigma(0) = 0")
    face = sorted(c.entry_face(0, 0))
    ck.check("the permutations with sigma(0) = 0 are the indices 0..5", face == [0, 1, 2, 3, 4, 5])
    six = {A[1], A[2], A[3], C[1], C[2], C[3]}
    okp = True
    for k in face:
        quo = sp.cancel(F[k] / (lam * (A[0] - C[0])))
        okp = okp and quo.free_symbols <= six and sp.Poly(quo, *sorted(six, key=str)).coeffs() != [] \
            and all(co.is_integer for co in sp.Poly(quo, *sorted(six, key=str)).coeffs())
    ck.check("for these six sigma, F_sigma / (lam (A0 - C0)) is a polynomial with integer coefficients in "
             "A1, A2, A3, C1, C2, C3: real when the six points are real", okp)
    ck.check("D_01 = L (A1 - C1)(A2 - A3)(C2 - C3): L times a real number, nonzero because D_01 != 0",
             sp.expand(D(0, 1) - Lsym * (A[1] - C[1]) * (A[2] - A[3]) * (C[2] - C[3])) == 0)
    ck.check("T_1 together with the face sigma(0)=0 is {0,..,7}; the f_j vanishing on it are minor vectors only",
             sorted(T | set(face)) == list(range(8)) and all(c.F_IS_MINOR[j] for j in c.generators(T | set(face))),
             "generators of {0..7}: %s" % c.generators(T | set(face)))
    ck.info("read: Re D_01 = 0 gives Re L = 0, so Re F_sigma = 0 for sigma(0) = 0; index 0 is a contact, so the "
            "minimum is 0 and indices 0..5 are contacts; the gap vector lies in N_{sigma(0)=0}, generated by 15 "
            "minor vectors (section 'Entry faces'), and Lemma 6.2 gives Phi >= 0 on SO(4).")

    # exact example: the stratum is not empty
    Ae = [sp.Rational(158, 85) - sp.Rational(14, 85) * sp.I, sp.Integer(-3), sp.Integer(-6), sp.Integer(6)]
    Ce = [2 + 4 * sp.I, sp.Integer(2), sp.Integer(3), sp.Integer(4)]
    le = sp.Rational(-59, 246) - sp.Rational(1, 123) * sp.I
    ge = [sp.re(sp.expand(val)) for val in c.F_products(Ae, Ce, le)]
    ck.check("exact example (A1..A3, C1..C3 real, distinct, no zero product): all 24 real parts are >= 0 and "
             "vanish exactly on {0,..,7}",
             min(ge) == 0 and [k for k in range(24) if ge[k] == 0] == list(range(8)),
             "gaps at 8..23: %s" % ", ".join(str(val) for val in ge[8:]))
    from fractions import Fraction
    gq = [Fraction(int(sp.numer(val)), int(sp.denom(val))) for val in ge]
    g8 = c.generators(range(8))
    xdec = c.exact_cone_decomposition(gq, [c.F_PROF[j] for j in g8])
    ck.check("exact minor decomposition of this gap vector with nonnegative rational coefficients on the minor "
             "vectors f_0..f_8 that vanish on {0,..,7}", xdec is not None,
             "" if xdec is None else ", ".join("%s f%d" % (xdec[n], j) for n, j in enumerate(g8) if xdec[n] != 0))

    # ------------------------------------------------------------------------------------------------
    ck.section("7b. Remark (not in the document): the other four cross collisions allowed by the hypotheses")
    ck.info("In the main branch the six points A1..A3, C1..C3 are real, Re(lam (A0 - C0)) = 0, and Re F_6 = 0 "
            "(index 6 is a contact and the common real part is 0).")
    a0, cr, ci, ar, ai, c0, pr = sp.symbols("a0 cr ci ar ai c0 p", real=True)
    # F_6 = lam (A0 - C1)(A1 - C0)(A2 - C2)(A3 - C3); (A2 - C2)(A3 - C3) is real and nonzero
    ck.check("F_6 = lam (A0 - C1)(A1 - C0)(A2 - C2)(A3 - C3)",
             sp.expand(F[6] - lam * (A[0] - C[1]) * (A[1] - C[0]) * (A[2] - C[2]) * (A[3] - C[3])) == 0)
    # (i) A0 real (this is the case when A0 = C2 or A0 = C3): A0 - C1 is real and nonzero, so Re(lam (A1 - C0)) = 0
    C0c = cr + sp.I * ci
    eA = sp.re(sp.expand(LL * (a0 - C0c)))
    eB = sp.re(sp.expand(LL * (u - C0c)))
    ck.check("A0 real: Re(lam (A0 - C0)) - Re(lam (A1 - C0)) = Re(lam) (A0 - A1), and for Re(lam) = 0 the first "
             "is Im(lam) Im(C0)", sp.expand(eA - eB - l * (a0 - u)) == 0 and sp.expand(eA.subs(l, 0) - m * ci) == 0)
    # (ii) C0 real (this is the case when A2 = C0 or A3 = C0): A1 - C0 is real and nonzero, so Re(lam (A0 - C1)) = 0
    A0c = ar + sp.I * ai
    eA2 = sp.re(sp.expand(LL * (A0c - c0)))
    eB2 = sp.re(sp.expand(LL * (A0c - pr)))
    ck.check("C0 real: Re(lam (A0 - C0)) - Re(lam (A0 - C1)) = Re(lam) (C1 - C0), and for Re(lam) = 0 the first "
             "is -Im(lam) Im(A0)", sp.expand(eA2 - eB2 - l * (pr - c0)) == 0
             and sp.expand(eA2.subs(l, 0) + m * ai) == 0)
    ck.info("read: A0 != A1 and C0 != C1 in the main branch, so in both cases lam is purely imaginary and all eight "
            "points are real; then every F_sigma is purely imaginary and the gap vector is 0.  So among the five "
            "cross collisions that the hypotheses allow, only A1 = C3 occurs with a nonzero gap vector.")

    # ------------------------------------------------------------------------------------------------
    ck.section("8. Zero contact products (second paragraph of the proof)")
    okz = True
    for i in range(4):
        for j in range(4):
            E = c.entry_face(i, j)
            okz = okz and all(sp.expand(F[k].subs(A[i], C[j])) == 0 for k in E) \
                and all(sp.expand(F[k].subs(A[i], C[j])) != 0 for k in range(24) if k not in E)
    ck.check("for every (i,j): A_i = C_j makes exactly the six products with sigma(i) = j vanish identically", okz)
    ck.info("read: if a contact product is 0 then the minimum of Re F is 0 and those six permutations are contacts; "
            "the gap vector lies in the cone of that entry face, generated by minor vectors.")

    c.check_lemma_6_2(ck)

    # ------------------------------------------------------------------------------------------------
    ck.section("9. Repeated eigenvalues inside one spectrum: the sub-case the document delegates to the "
               "one-circle theorem.  Own elementary argument (not in the document)")
    ck.info("If A_i = A_j then F_sigma = F_{sigma o (i j)} for all sigma; if C_i = C_j then F_sigma = F_{(i j) o sigma}.")
    ck.info("So the gap vector g and the contact set are invariant under that involution s, and the contact set "
            "contains T' = T_1 u s(T_1).")
    closed_by_enlargement, need_sym = [], []
    okinv = True
    for kind in ("A", "C"):
        for i in range(4):
            for j in range(i + 1, 4):
                if kind == "A":
                    smap = [c.PIDX[tuple(pp[{i: j, j: i}.get(k, k)] for k in range(4))] for pp in c.PERMS]
                    subs = {A[j]: A[i]}
                else:
                    smap = [c.PIDX[tuple({i: j, j: i}.get(val, val) for val in pp)] for pp in c.PERMS]
                    subs = {C[j]: C[i]}
                okinv = okinv and all(sp.expand((F[k] - F[smap[k]]).subs(subs)) == 0 for k in range(24))
                Tp = set(T) | set(smap[k] for k in T)
                if all(c.F_IS_MINOR[n] for n in c.generators(Tp)):
                    closed_by_enlargement.append("%s%d=%s%d" % (kind, i, kind, j))
                else:
                    need_sym.append((kind, i, j, smap, Tp))
    ck.check("invariance F_sigma = F_{s(sigma)} under each of the 12 coincidences A_i = A_j, C_i = C_j", okinv)
    ck.check("11 of the 12 coincidences enlarge the contact set to a set T' on which only minor vectors vanish",
             len(closed_by_enlargement) == 11, ", ".join(closed_by_enlargement))
    ck.check("the remaining coincidence is A2 = A3, whose involution maps T_1 to itself",
             len(need_sym) == 1 and need_sym[0][:3] == ("A", 2, 3) and need_sym[0][4] == set(T))
    kind, i, j, smap, Tp = need_sym[0]

    def act(f):
        out = [0] * 24
        for k in range(24):
            out[smap[k]] = f[k]
        return tuple(out)

    img = {n: c.PROF_INDEX.get(act(c.F_PROF[n])) for n in gens}
    ck.check("that involution permutes the 13 generators of N_{T_1}: minor vectors go to minor vectors and "
             "f_17 <-> f_18", all(img[n] in gens for n in gens) and img[17] == 18 and img[18] == 17
             and all(c.F_IS_MINOR[img[n]] for n in gens if c.F_IS_MINOR[n]),
             ", ".join("%d->%d" % (n, img[n]) for n in gens))
    lhs = c.vadd(c.F_PROF[17], c.F_PROF[18])
    rhs = (0,) * 24
    for n in (1, 4, 8, 15, 16):
        rhs = c.vadd(rhs, c.F_PROF[n])
    ck.check("profile identity f_17 + f_18 = f_1 + f_4 + f_8 + f_15 + f_16 (all 24 entries)", lhs == rhs)
    ck.info("read: for A2 = A3, g = s.g = (g + s.g)/2 = sum_j a_j (f_j + s.f_j)/2 with a_j >= 0; every summand is a "
            "nonnegative combination of minor vectors by the two checks above.  Hence every repeated-eigenvalue "
            "sub-case of T_1 has a minor decomposition, without the one-circle theorem.")

    c.crosscheck_with_audit_lib(ck)
    c.compare_with_supplement_data(ck, T)
    return ck.finish()


if __name__ == "__main__":
    sys.exit(main())
