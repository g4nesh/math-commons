#!/usr/bin/env python3
"""Independent exact replay of the proof of contact case T_2 = {0,1,2,3,6,8}.

Proof replayed: supplement, proofs/manin_transfer_20261021/independent_moc4_six_contact_reduction_20261022.pdf,
Section 5.2 (pp. 7-8; the document calls this set "T_1").  Every check is exact (integers, rationals,
symbolic); nothing from the supplement is executed or imported.  The supplement is read only in the last
section, as data, to compare lists, and only when the environment variable MOC4_SUPP points to it.

Run:  MOC4_SUPP=/path/to/supplement python replay_T2.py        exit code 0 iff every check passes.
Packages: sympy (plus the standard library).  Shared objects: t12_common.py in this directory.

Notation of the document: A_i = a_i, C_j = -b_j, F_sigma = lam * prod_i (A_i - C_sigma(i)), D_rs = F_r - F_s,
g_sigma = Re(F_sigma) - min, permutations indexed 0..23 in lexicographic order.
"""
import os
import sys
from fractions import Fraction

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sympy as sp

import t12_common as c

T = frozenset({0, 1, 2, 3, 6, 8})                       # the document's contact set, p. 7
J_DOC = [0, 1, 2, 3, 4, 5, 6, 9, 10, 15, 16, 17]        # the document's list of generators, p. 8


def main():
    ck = c.Checker("Exact replay, contact case T_2 = {0,1,2,3,6,8} "
                   "(independent_moc4_six_contact_reduction_20261022.pdf, Section 5.2, pp. 7-8)")
    c.check_basic_objects(ck)

    # ------------------------------------------------------------------------------------------------
    ck.section("1. Identification of the contact set (task item 1)")
    ck.check("document's set equals the paper's T_2 (Appendix E) as a set of lexicographic indices",
             T == c.PAPER_T[2])
    ck.info("permutations of T_2 (0-based one-line): %s" % ", ".join("%d:%s" % (k, c.oneline(c.PERMS[k]))
                                                                    for k in sorted(T)))
    ident = c.G_ELEMENTS.index((c.IDENT, c.IDENT, 0))
    ck.check("the group element relating them is the identity (rho = rho' = id, no transposition)",
             c.g_image(ident, T) == c.PAPER_T[2])
    orb = c.orbit_of(T)
    ck.check("orbit of T_2 under the group of order 1152 has 576 elements; T_2 is its lexicographic minimum",
             len(orb) == 576 and c.orbit_min(T) == tuple(sorted(T)), "orbit size %d" % len(orb))
    ck.check("T_2 is not in the orbit of T_1", frozenset(c.PAPER_T[1]) not in orb)

    # ------------------------------------------------------------------------------------------------
    ck.section("2. Generators of the contact cone N_{T_2} (task item 2)")
    gens = c.generators(T)
    ck.check("f_j vanishing on T_2 = the document's list J = {0,1,2,3,4,5,6,9,10,15,16,17}", gens == J_DOC, str(gens))
    ck.check("all are minor vectors except f_17 = h_{3421,+}",
             [j for j in gens if not c.F_IS_MINOR[j]] == [17] and c.F_NAME[17] == "h3421+")
    ck.info("generators: %s" % ", ".join("f%d=%s" % (j, c.F_NAME[j]) for j in gens))
    d, rays = c.contact_cone_extreme_rays(T)
    ck.check("direct extreme-ray enumeration of N_{T_2} (dimension %d) returns exactly these 12 profiles "
             "(no use of the 82-facet theorem)" % d, set(rays) == set(c.F_PROF[j] for j in gens))

    # ------------------------------------------------------------------------------------------------
    ck.section("3. The spectral identity (6) and the equality g_7 = g_9 (task item 3)")
    A = sp.symbols("A0:4")
    C = sp.symbols("C0:4")
    lam = sp.Symbol("lam")
    F = c.F_products(A, C, lam)

    def D(r, s):
        return F[r] - F[s]

    d02 = lam * (A[0] - C[0]) * (A[3] - C[3]) * (A[1] - A[2]) * (C[1] - C[2])
    ck.check("D_02 = lam (A0-C0)(A3-C3)(A1-A2)(C1-C2)", sp.expand(D(0, 2) - d02) == 0)
    ck.check("identity (6): D_79 D_02 = D_68 D_13 as a polynomial identity in A, C, lam",
             sp.expand(D(7, 9) * D(0, 2) - D(6, 8) * D(1, 3)) == 0)
    own = {
        (6, 8): lam * (A[0] - C[1]) * (A[3] - C[3]) * (A[1] - A[2]) * (C[0] - C[2]),
        (1, 3): lam * (A[0] - C[0]) * (A[2] - C[3]) * (A[1] - A[3]) * (C[1] - C[2]),
        (7, 9): lam * (A[0] - C[1]) * (A[2] - C[3]) * (A[1] - A[3]) * (C[0] - C[2]),
    }
    ck.check("own factorisations of D_68, D_13, D_79 (not displayed in the document)",
             all(sp.expand(D(r, s) - val) == 0 for (r, s), val in own.items()))
    ck.check("the endpoints of D_02, D_68, D_13 lie in T_2 (so these three have real part 0 on the stratum); "
             "7 and 9 do not", {0, 2, 6, 8, 1, 3} <= T and 7 not in T and 9 not in T)
    cross_in_T = set((i, c.PERMS[k][i]) for k in T for i in range(4))
    cross, wA, wC, has_lam, left, coeff = c.classify_factors(D(0, 2), A, C, lam)
    ck.check("D_02 != 0 needs exactly: lam != 0, A0 != C0, A3 != C3 (both factors of the contact product F_0), "
             "A1 != A2, C1 != C2",
             has_lam and not left and cross == {(0, 0), (3, 3)} and wA == {(1, 2)} and wC == {(1, 2)}
             and cross <= set((i, c.PERMS[0][i]) for i in range(4)))
    ok68 = all(c.classify_factors(D(r, s), A, C, lam)[0] <= cross_in_T for (r, s) in ((6, 8), (1, 3)))
    ck.check("D_68 and D_13 are also nonzero under the document's hypotheses (their cross factors divide contact "
             "products); this is stated in the document but not needed", ok68)
    d1, d2, d3 = sp.symbols("d1 d2 d3", real=True)
    ck.check("if D_02 = i d1 (d1 != 0), D_68 = i d2, D_13 = i d3 with d1, d2, d3 real, then D_79 = D_68 D_13 / D_02 "
             "= i d2 d3 / d1 has real part 0",
             sp.re(sp.simplify((sp.I * d2) * (sp.I * d3) / (sp.I * d1))) == 0)
    ck.info("read: hence Re F_7 = Re F_9, i.e. g_7 = g_9, whenever D_02 != 0.  Nothing is claimed about 7 or 9 "
            "being contacts.")

    # ------------------------------------------------------------------------------------------------
    ck.section("4. The two integer identities on p. 8")
    lhs = c.vadd(c.F_PROF[17], c.F_PROF[9])
    rhs = c.vadd(c.vadd(c.F_PROF[4], c.F_PROF[10]), c.F_PROF[16])
    ck.check("f_17 + f_9 = f_4 + f_10 + f_16 (all 24 entries)", lhs == rhs)
    ell = {j: c.F_PROF[j][7] - c.F_PROF[j][9] for j in J_DOC}
    ck.check("l(g) = g_7 - g_9 takes the values l(f_6) = -1, l(f_9) = 1, l(f_17) = -1 and 0 on the other f_j, j in J",
             {j: val for j, val in ell.items() if val} == {6: -1, 9: 1, 17: -1}, str(ell))

    ck.section("5. The coefficient transfer (formal coefficients)")
    # coefficient vectors over the basis a_j, j in J_DOC
    idx = {j: n for n, j in enumerate(J_DOC)}

    def e(j):
        v = [0] * len(J_DOC)
        v[idx[j]] = 1
        return v

    beta = {j: e(j) for j in J_DOC if j != 17}
    beta[9] = e(6)
    for j in (4, 10, 16):
        beta[j] = [a + b for a, b in zip(beta[j], e(17))]
    ck.check("each new coefficient beta_j is a sum of old coefficients a_k with coefficients in {0, 1}: "
             "beta >= 0 whenever a >= 0", all(set(v) <= {0, 1} for v in beta.values()))
    ok = True
    for s in range(24):
        orig = [c.F_PROF[j][s] for j in J_DOC]                                    # coefficient of a_j in g_s
        conv = [sum(beta[j][n] * c.F_PROF[j][s] for j in beta) for n in range(len(J_DOC))]
        defect = [(e(9)[n] - e(6)[n] - e(17)[n]) * c.F_PROF[9][s] for n in range(len(J_DOC))]
        ok = ok and [o - cv for o, cv in zip(orig, conv)] == defect
    ck.check("sum_J a_j f_j - sum_{J minus 17} beta_j f_j = (a_9 - a_6 - a_17) f_9 identically in the a_j "
             "(all 24 entries)", ok)
    ck.check("l(sum_J a_j f_j) = a_9 - a_6 - a_17, so g_7 = g_9 is the equation a_9 = a_6 + a_17",
             [ell[j] for j in J_DOC] == [e(9)[n] - e(6)[n] - e(17)[n] for n in range(len(J_DOC))])

    ck.section("6. Independent route to the same conclusion: the slice N_{T_2} n {g_7 = g_9}")
    lvec = [0] * 24
    lvec[7], lvec[9] = 1, -1
    d_s, rays_s = c.contact_cone_extreme_rays(T, extra_eqs=[lvec])
    expected = set(c.F_PROF[j] for j in J_DOC if j not in (6, 9, 17)) | {c.vadd(c.F_PROF[6], c.F_PROF[9])}
    ck.check("direct extreme-ray enumeration of the slice (dimension %d): its %d extreme rays are the nine minor "
             "generators with l = 0 and f_6 + f_9; every one is a sum of minor vectors"
             % (d_s, len(rays_s)), set(rays_s) == expected and len(rays_s) == 10)

    c.check_entry_faces(ck)

    # ------------------------------------------------------------------------------------------------
    ck.section("7. Zero contact products (first sentence of the proof)")
    okz = True
    for i in range(4):
        for j in range(4):
            E = c.entry_face(i, j)
            okz = okz and all(sp.expand(F[k].subs(A[i], C[j])) == 0 for k in E) \
                and all(sp.expand(F[k].subs(A[i], C[j])) != 0 for k in range(24) if k not in E)
    ck.check("for every (i,j): A_i = C_j makes exactly the six products with sigma(i) = j vanish identically", okz)
    ck.info("read: if a contact product is 0 the minimum of Re F is 0 and those six permutations are contacts; "
            "the gap vector lies in the cone of that entry face, generated by minor vectors.")
    ck.info("cross pairs (i,j) with A_i - C_j dividing a contact product: %s" % sorted(cross_in_T))

    c.check_lemma_6_2(ck)

    # ------------------------------------------------------------------------------------------------
    ck.section("8. Repeated eigenvalues inside one spectrum: the sub-case the document delegates to the "
               "one-circle theorem.  Own elementary argument (not in the document)")
    ck.info("If A_i = A_j then F_sigma = F_{sigma o (i j)}; if C_i = C_j then F_sigma = F_{(i j) o sigma}.  The contact "
            "set is invariant under that involution s and contains T' = T_2 u s(T_2).")
    okinv, closed, rows = True, [], []
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
                    closed.append("%s%d=%s%d" % (kind, i, kind, j))
                rows.append("%s%d=%s%d: T' = %s" % (kind, i, kind, j, sorted(Tp)))
    ck.check("invariance F_sigma = F_{s(sigma)} under each of the 12 coincidences A_i = A_j, C_i = C_j", okinv)
    ck.check("all 12 coincidences enlarge the contact set to a set T' on which only minor vectors vanish "
             "(f_17 does not vanish on T')", len(closed) == 12, ", ".join(closed))
    for r in rows:
        ck.info(r)
    ck.info("read: so every repeated-eigenvalue sub-case of T_2 has a minor decomposition, without the one-circle "
            "theorem.  In particular A1 = A2 gives T' containing the entry face {0,..,5} and C1 = C2 gives "
            "T' = {0,1,2,3,6,8,12,14}; these are the only two coincidences in which D_02 = 0.")

    # ------------------------------------------------------------------------------------------------
    ck.section("9. Exact example: the stratum with contact set exactly T_2 is not empty")
    R, I = sp.Rational, sp.I
    Ae = [R(54791, 58393) + R(260077, 58393) * I, R(-357, 53) + R(57, 53) * I, -4 + 4 * I, 1 - I]
    Ce = [-2 + I, R(-391, 183) + R(132, 61) * I, 2 + 4 * I, 3 * I]
    le = R(-295600505, 25595661) - R(100299065, 8531887) * I
    Fe = [sp.expand(val) for val in c.F_products(Ae, Ce, le)]
    re = [sp.re(val) for val in Fe]
    mn = min(re)
    ge = [val - mn for val in re]
    ck.check("all 24 gaps are >= 0 and vanish exactly on T_2 = {0,1,2,3,6,8}",
             min(ge) == 0 and [k for k in range(24) if ge[k] == 0] == sorted(T))
    ck.check("the eight points are pairwise distinct inside each spectrum and D_02 != 0 at this point",
             len(set(Ae)) == 4 and len(set(Ce)) == 4 and sp.expand(Fe[0] - Fe[2]) != 0)
    ck.check("g_7 = g_9 > 0 at this point", ge[7] == ge[9] and ge[7] > 0, "g_7 = g_9 = %s" % ge[7])
    gq = [Fraction(int(sp.numer(val)), int(sp.denom(val))) for val in ge]
    minor_gens = [j for j in J_DOC if j != 17]
    x = c.exact_cone_decomposition(gq, [c.F_PROF[j] for j in minor_gens])
    ck.check("exact minor decomposition of this gap vector with nonnegative rational coefficients on the 11 minor "
             "generators of N_{T_2}", x is not None,
             "" if x is None else ", ".join("%s f%d" % (x[n], j) for n, j in enumerate(minor_gens) if x[n] != 0))

    c.crosscheck_with_audit_lib(ck)
    c.compare_with_supplement_data(ck, T)
    return ck.finish()


if __name__ == "__main__":
    sys.exit(main())
