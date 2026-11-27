"""Independent exact replay of Case T_88 = {1,5,6,10,13,14} (group V) of D. Zheng's proof of the
Marcus-de Oliveira conjecture for n = 4.

Proof document replayed: MOC4_SUPP/proofs/manin_transfer_20261021/final_referee_inventory87_uniform_20261023.pdf
(4 pages; the document calls the set "inventory 87").  Every polynomial is rebuilt here from the
definitions (the matrix P, the cubic C, the vectors f_j of the paper); nothing from the supplement is
executed or imported.

The script also re-establishes, with own code, every finite assertion made by the supplement's
"portable check" MOC4_SUPP/checks/fresh_T88_audit.py (read as text only).  Checks that answer one of
its assertions carry a tag [PC n]; the numbering is that of the report T88.md.

If MOC4_SUPP is set, the document's certificate JSON is read as DATA and COMPARED with the rebuilt
objects (C, Q, G, the lift, the group elements, the hashes of the ten tables).  No step depends on it.

Run:   MOC4_SUPP=/path/to/supplement  python replay_T88.py
Exit status 0 iff every check passes.  Checks tagged EXACT use integer/rational/symbolic
arithmetic only; checks tagged NUMERIC are falsification attempts (floating point, or exact
evaluation at sampled points).
"""
import json
import os
import random
import sys
import time
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import replay_T72_common as M                                   # noqa: E402
from replay_T72_common import F_VALUES, Ring                    # noqa: E402

DIR = "proofs/manin_transfer_20261021/"
DOC = DIR + "final_referee_inventory87_uniform_20261023.pdf"
CERT_JSON = DIR + "final_referee_inventory87_uniform_exact_20261023.json"
PORTABLE = "checks/fresh_T88_audit.py"
SHA = {
    DOC: "63542699a4d7e955768d597d9ab04d7bedaa273463e81e2e5bd104f9c1cbd0b8",
    CERT_JSON: "1602e6b3745900572cfced9f7c6bcc292bfdf4e8a9844bd56ddb3234d6a2a3a6",
    PORTABLE: "5aa37a0064868df5990d4638147673d549d7291ebcfd1339cec27271a63d606f",
}

T = (1, 5, 6, 10, 13, 14)                     # the paper's T_88 = the document's set
GT = ((0, 1, 2, 3), (0, 1, 2, 3), 1)          # sigma -> sigma^-1 (the document's "spectral transpose")
ORDER = [0, 1, 6, 10, 27, 36, 59, 30]         # g = q f0 + a f1 + v f6 + b f10 + w f27 + t f36 + c f59 + d f30
NAMES = "q,a,v,b,w,t,c,d,h"
NATIVE = ("q", "a", "v", "b", "w", "t", "c")
EXCL = (20, 16, 11, 7, 12, 3, 0)              # permutation indices dual to the first seven generators
DEGS = (3, 2, 3, 2, 3, 3, 2)                  # the document's clearing degrees (11), d = 1
PIVOT_ROWS = (0, 1, 2, 3, 4, 6, 7, 8, 9, 10, 12, 13, 16, 18)
DOC_COUNTS = {"balanced": 20241, "dominant": 27372, "solo": 450}
DOC_TOTAL = 185823


def load_json(relpath):
    root = os.environ.get("MOC4_SUPP")
    if not root:
        return None
    path = os.path.join(root, relpath)
    if not os.path.isfile(path):
        return None
    with open(path) as fh:
        return json.load(fh)


def main():
    import sympy as sp
    import numpy
    import scipy
    import flint
    ck = M.Checker("Case T_88 = {1,5,6,10,13,14} (document 'inventory 87'): independent exact replay")
    ck.info("python %s, sympy %s, python-flint %s (exact integer polynomials for the chart tables and the "
            "fraction-free elimination; integer determinants at test points), numpy %s (int64 in det_linear_forms; "
            "floats in NUMERIC checks), scipy %s (NUMERIC only)"
            % (sys.version.split()[0], sp.__version__, flint.__version__, numpy.__version__, scipy.__version__))

    # ------------------------------------------------------------------------------------
    ck.section("A. Provenance and foundations")
    ck.info("elapsed %.0f s" % (time.time() - ck.t0))
    for rel in (DOC, CERT_JSON, PORTABLE):
        M.check_document(ck, rel, SHA[rel])
    M.check_foundations(ck)
    ck.info("[PC 1] is the first foundations check (every row of P is a polynomial identity; rank 14); "
            "[PC 11] is the check 'P kappa = (1,...,1)'; the vanishing of the fourteen partial derivatives of C on "
            "products (document p. 3) is the fifth foundations check")
    Msub = sp.Matrix([M.P[r] for r in PIVOT_ROWS])
    ck.check("EXACT", "[PC 2] the rows 0,1,2,3,4,6,7,8,9,10,12,13,16,18 of P form a matrix of determinant 1 "
             "(so its inverse is an integer matrix of determinant 1)", Msub.det() == 1, "det = %s" % Msub.det())
    M.check_radial_necessity(ck)
    M.crosscheck_supplement_tables(ck)
    M.crosscheck_audit_library(ck, [72, 88])

    # ------------------------------------------------------------------------------------
    ck.section("B. Contact set and the transpose (document pp. 1-2)")
    ck.info("elapsed %.0f s" % (time.time() - ck.t0))
    ck.check("EXACT", "the document's set {1,5,6,10,13,14} is the paper's T_88 (Appendix E) as a set of "
             "lexicographic indices: group element = identity", tuple(sorted(T)) == M.PAPER_T[88])
    rT = M.orbit_rep(T)
    ck.check("EXACT", "orbit representative of T_88 under G; T_72 lies in another orbit",
             rT == (0, 3, 7, 8, 19, 20) and M.orbit_rep(M.PAPER_T[72]) != rT, "T_88 -> %s" % (rT,))
    gmt = M.g_index_map(GT)
    ck.check("EXACT", "[PC 10a] sigma -> sigma^-1 maps T onto itself",
             tuple(sorted(gmt[j] for j in T)) == T, "index map %s" % (list(gmt),))
    M.check_group_on_X(ck, [("the transpose", GT)])

    # ------------------------------------------------------------------------------------
    ck.section("C. The contact kernel and its signed generator representation (document pp. 1-2, eq. (1)-(3))")
    ck.info("elapsed %.0f s" % (time.time() - ck.t0))
    gens = M.generators(T)
    ck.check("EXACT", "the f_j of Appendix D vanishing on T are exactly f0, f1, f6, f10, f27, f30, f34, f36, f59 "
             "(the document's nine profiles)", gens == [0, 1, 6, 10, 27, 30, 34, 36, 59],
             ", ".join("f%d = %s" % (j, M.f_name(j)) for j in gens))
    ck.check("EXACT", "[PC 3] these nine value vectors lie in the image of P with integer coordinates; "
             "f30 = h_{4321,-} and f34 = h_{4321,+} are the two that are not minor vectors",
             [j for j in gens if not M.F_IS_MINOR[j]] == [30, 34]
             and all(x.denominator == 1 for j in gens for x in M.F_COORD[j])
             and all(M.values(M.F_COORD[j]) == list(F_VALUES[j]) for j in gens))
    rkT = M.rank([M.P[s] for s in T])
    null = M.nullspace([M.P[s] for s in T], 14)
    ck.check("EXACT", "[PC 4] the six rows of P at T have rank 6 (solution space of dimension 8) and the eight vectors "
             "f0, f1, f6, f10, f27, f36, f59, f30 are linearly independent",
             rkT == 6 and len(null) == 8 and M.rank([M.F_COORD[j] for j in ORDER]) == 8)
    VL = [[F_VALUES[j][s] for j in ORDER] for s in range(24)]          # 24 x 8 value matrix
    ck.check("EXACT", "[PC 5] the eight vectors vanish at the six indices of T",
             all(VL[s][k] == 0 for s in T for k in range(8)))
    unit = lambda *ks: [int(k in ks) for k in range(8)]
    ck.check("EXACT", "[PC 6] rows 20, 16, 11, 7, 12, 3, 0 of the value matrix are the unit vectors e_0..e_6: for "
             "g = q f0 + a f1 + v f6 + b f10 + w f27 + t f36 + c f59 + d f30 the coefficients q, a, v, b, w, t, c "
             "are the gaps g20, g16, g11, g7, g12, g3, g0",
             [VL[s] for s in EXCL] == [unit(k) for k in range(7)])
    ck.check("EXACT", "[PC 7] row 19 minus row 20 is e_7: d = g19 - g20 (so the form (1) is unique and d may have "
             "either sign)", [x - y for x, y in zip(VL[19], VL[20])] == unit(7))
    ck.check("EXACT", "[PC 8] (3): f30 + f34 = f0 + f6 + f27 + f36 on all 24 permutations; the nine vectors have "
             "rank 8, so this is their only linear relation up to a factor",
             [x + y for x, y in zip(F_VALUES[30], F_VALUES[34])]
             == [sum(z) for z in zip(F_VALUES[0], F_VALUES[6], F_VALUES[27], F_VALUES[36])]
             and M.rank([M.F_COORD[j] for j in gens]) == 8)
    ck.check("EXACT", "[PC 9] rows 19, 15, 4, 8 are e_0 + e_7, e_2 + e_7, e_4 + e_7, e_5 + e_7: g19 = q + d, "
             "g15 = v + d, g4 = w + d, g8 = t + d", [VL[19], VL[15], VL[4], VL[8]]
             == [unit(0, 7), unit(2, 7), unit(4, 7), unit(5, 7)])
    sq, sa, sv, sb, sw, st, sc, sd = sp.symbols("q a v b w t c d")
    lhs = M.combo_values(list(zip(ORDER, (sq, sa, sv, sb, sw, st, sc, sd))))
    rhs = M.combo_values(list(zip([0, 1, 6, 10, 27, 36, 59, 34], (sq + sd, sa, sv + sd, sb, sw + sd, st + sd, sc, -sd))))
    ck.check("EXACT", "signed side: q f0 + a f1 + v f6 + b f10 + w f27 + t f36 + c f59 + d f30 = (q+d) f0 + a f1 + "
             "(v+d) f6 + b f10 + (w+d) f27 + (t+d) f36 + c f59 + (-d) f34 on all 24 permutations; for d < 0 the "
             "coefficients q+d, v+d, w+d, t+d are the gaps g19, g15, g4, g8, so all eight coefficients on the "
             "right are >= 0", all(sp.expand(x - y) == 0 for x, y in zip(lhs, rhs)))
    img = {j: M.g_facet(gmt, j) for j in gens}
    ck.check("EXACT", "[PC 10b] the transpose exchanges f30 and f34 and permutes the seven minor vectors: "
             "f0 <-> f6, f27 <-> f36, and f1, f10, f59 are fixed", img == {0: 6, 6: 0, 27: 36, 36: 27, 1: 1, 10: 10,
                                                                          59: 59, 30: 34, 34: 30}, str(img))
    # the transposed vector of the rewritten source is of the form (1) with d' = -d > 0
    tr = M.combo_values([(img[j], cf) for j, cf in
                         zip([0, 1, 6, 10, 27, 36, 59, 34], (sq + sd, sa, sv + sd, sb, sw + sd, st + sd, sc, -sd))])
    form1 = M.combo_values(list(zip(ORDER, (sv + sd, sa, sq + sd, sb, st + sd, sw + sd, sc, -sd))))
    perm_ok = True
    for s in range(24):                        # (g.c) has the value of c at sigma on the permutation g.sigma
        perm_ok &= sp.expand(tr[gmt[s]] - lhs[s]) == 0
    ck.check("EXACT", "for d < 0 the transposed gap vector is of the form (1) with coefficients "
             "(q', a', v', b', w', t', c', d') = (v+d, a, q+d, b, t+d, w+d, c, -d), all >= 0 and d' > 0: so the case "
             "d < 0 reduces to d > 0 (T is preserved, and by section B the transpose preserves X and the set of "
             "values of the Leibniz sum)", perm_ok and all(sp.expand(x - y) == 0 for x, y in zip(tr, form1)))
    dim, rays = M.extreme_rays_of_contact_cone(T)
    ck.check("EXACT", "direct enumeration of the extreme rays of N_T = {Pc >= 0, (Pc)_T = 0} (dimension %d; no use "
             "of the 82-facet list) gives exactly the nine vectors above; the cone is not simplicial" % dim,
             rays == sorted(F_VALUES[j] for j in gens), "%d rays" % len(rays))
    M.check_minor_leibniz(ck, [0, 1, 6, 10, 27, 36, 59])

    # ------------------------------------------------------------------------------------
    ck.section("D. The three SOC atoms and the capacity H3 (document p. 2, eq. (4)-(6))")
    ck.info("elapsed %.0f s" % (time.time() - ck.t0))
    rng = random.Random(20261188)
    for tgt in ((1, 10, 30), (1, 59, 30), (10, 59, 30)):
        M.check_transport(ck, "SOC", tgt, rng)
    ck.check("EXACT", "the document's reference tuple (f1, f9; f22) is the audit's SOC representative",
             M.UNIVERSAL["SOC"][0] == (1, 9, 22))
    M.check_H3_allocation(ck, 1, 10, 59, 30)
    M.replay_SOC_representative(ck)

    # ------------------------------------------------------------------------------------
    ck.section("E. Radial necessity in the generators (document pp. 2-3, eq. (7)-(9))")
    ck.info("elapsed %.0f s" % (time.time() - ck.t0))
    rd = M.radial_data(ck, ORDER, NAMES)
    R9, C0, Q, G = rd["R"], rd["C0"], rd["Q"], rd["G"]
    syms9 = sp.symbols(NAMES)
    q9, a9, v9, b9, w9, t9, c9, d9, h9 = syms9
    aux = (q9, v9, w9, t9)
    import itertools
    e1 = sum(aux)
    e2 = sum(x * y for x, y in itertools.combinations(aux, 2))
    e3 = sum(x * y * z for x, y, z in itertools.combinations(aux, 3))
    Cdoc = R9.from_sympy(a9 * b9 * c9 + e3 + d9 * e2 + d9 ** 2 * e1 + d9 ** 3)
    Qdoc = R9.from_sympy(a9 * b9 + a9 * c9 + b9 * c9 + d9 ** 2)
    ck.check("EXACT", "[PC 12] (7): C(g + h kappa) = C0 + h Q with C0 = abc + e3(q,v,w,t) + d e2(q,v,w,t) + "
             "d^2 e1(q,v,w,t) + d^3 and Q = ab + ac + bc + d^2 (same sign as the paper's cubic)",
             C0 == Cdoc and Q == Qdoc)
    kd3 = R9.key((0, 0, 0, 0, 0, 0, 0, 3, 0))
    kd2 = R9.key((0, 0, 0, 0, 0, 0, 0, 2, 0))
    ck.check("EXACT", "all %d coefficients of C0 and all %d coefficients of Q are +1; C0 contains d^3 and Q contains "
             "d^2: C0 >= d^3 > 0 and Q >= d^2 > 0 on the closed orthant with d > 0, so h = -C0/Q"
             % (len(C0), len(Q)), all(x == 1 for x in C0.values()) and all(x == 1 for x in Q.values())
             and C0.get(kd3) == 1 and Q.get(kd2) == 1)
    kd9 = R9.key((0, 0, 0, 0, 0, 0, 0, 9, 0))
    if G.get(kd9, 0) < 0:
        G = M.p_scale(G, -1)
    degs8 = [R9.deg(G, nm) for nm in R9.names[:8]]
    ck.check("EXACT", "(8), (11): G has 2285 terms, is homogeneous of degree 9, does not involve h, its d^9 coefficient "
             "is +1, and its degrees in (q, a, v, b, w, t, c, d) are (3, 2, 3, 2, 3, 3, 2, 9)",
             len(G) == 2285 and all(sum(R9.exps(k)) == 9 and R9.exps(k)[8] == 0 for k in G)
             and G.get(kd9) == 1 and degs8 == [3, 2, 3, 2, 3, 3, 2, 9], "degrees %s" % degs8)
    ck.info("G has %d negative coefficients in the native variables, so its sign is not visible without the charts"
            % sum(1 for cf in G.values() if cf < 0))

    def permuted(p, swaps):
        idx = list(range(9))
        for x, y in swaps:
            i, j = R9.pos[x], R9.pos[y]
            idx[i], idx[j] = idx[j], idx[i]
        out = {}
        for k, cf in p.items():
            e = R9.exps(k)
            out[R9.key(tuple(e[idx[i]] for i in range(9)))] = cf
        return out
    sym_t = permuted(G, [("q", "v"), ("w", "t")]) == G
    sym_abc = [permuted(G, [pr]) == G for pr in (("a", "b"), ("a", "c"), ("b", "c"))]
    ck.check("EXACT", "G is invariant under q <-> v, w <-> t (the action of the transpose on the coefficients) and "
             "is NOT invariant under any transposition of a, b, c alone; this is why the document computes the "
             "dominant and solo charts for each of a, b, c separately",
             sym_t and not any(sym_abc), "a<->b %s, a<->c %s, b<->c %s" % tuple(sym_abc))
    M.check_bareiss(ck, rd, ORDER)
    ck.info("[PC 14] is answered by the two determinant computations above and the identity N = -4 Q^2 G^2, here as "
            "an identity of homogeneous polynomials in eight variables (the portable check works at d = 1)")
    cert = load_json(CERT_JSON)
    if cert is None:
        ck.info("MOC4_SUPP not set (or file absent): the certificate JSON is not compared")
    else:
        def load(rows):
            return {R9.key(tuple(int(x) for x in e)): int(cf) for e, cf in rows}
        grows = cert["G2285"]
        fmt_ok = (len(grows) == 2285 and all(len(e) == 9 and e[-1] == 0 and sum(e) == 9 for e, _ in grows)
                  and len({tuple(e) for e, _ in grows}) == 2285
                  and [max(e[i] for e, _ in grows) for i in range(8)] == [3, 2, 3, 2, 3, 3, 2, 9])
        ck.check("EXACT", "[PC 13] supplement data (JSON, comparison only): its axes are q,a,v,b,w,t,c,d,height; its G "
                 "list has 2285 distinct exponent vectors of length 9 with last entry 0 and sum 9, maximal exponents "
                 "(3,2,3,2,3,3,2,9), d^9 coefficient 1; it EQUALS the G computed here; its C and Q equal C0 and Q",
                 fmt_ok and cert["native_axes"] == list(R9.names[:8]) + ["height"] and load(grows) == G
                 and load(cert["C"]) == C0 and load(cert["Q"]) == Q)
        lift_ok = [[int(x) for x in row] for row in cert["source_lift14"]] == \
            [[int(M.F_COORD[j][i]) for j in ORDER] for i in range(14)]
        okc = lift_ok and cert["source_rows24"] == VL and cert["contacts"] == list(T) and cert["source8"] == ORDER \
            and cert["source9_crosscheck"] == gens and cert["actual_bad_swap"]["action24"] == list(gmt)
        want = {"ab": (1, 10, 30), "ac": (1, 59, 30), "bc": (10, 59, 30)}
        for nm, rec in cert["SOC_atoms"].items():
            ge = (tuple(rec["row"]), tuple(rec["column"]), int(bool(rec["transpose"])))
            ac = rec["action24"]
            okc = okc and list(M.g_index_map(ge)) == ac and tuple(rec["target"]) == want.get(nm) \
                and tuple(rec["mother"]) == (1, 9, 22)
            # the certificate's convention: the vector with value f_i(ac[j]) at j is the target vector
            okc = okc and all(tuple(F_VALUES[i][ac[j]] for j in range(24)) == F_VALUES[tg]
                              for i, tg in zip(rec["mother"], rec["target"]))
        ck.check("EXACT", "supplement data (JSON, comparison only): its 14 x 8 lift and its 24 x 8 value matrix equal "
                 "the coordinates and values of the eight generators computed here; its transpose is the index map "
                 "of sigma -> sigma^-1; its three recorded index maps are elements of G that carry (f1, f9; f22) to "
                 "(f1, f10; f30), (f1, f59; f30), (f10, f59; f30) (pull-back convention of the certificate)",
                 okc and set(cert["SOC_atoms"]) == set(want))

    # ------------------------------------------------------------------------------------
    ck.section("F. Charts of the region H3 < 1 (document pp. 3-4, eq. (10), (12), (13))")
    ck.info("elapsed %.0f s" % (time.time() - ck.t0))
    CH = M.ChartRing("Q0,V,W,T,X,Y,Z")
    Q0, V, W, Tt, X, Y, Z = CH.syms
    Dn, L = 1 + X + Y + Z, 1 + X
    zero = sp.Integer(0)
    Mexp = (1 + Y * L ** 2) / (L * (1 + Y))
    leaves = [Z / (L * (1 + Z)), 1 / (L * (1 + Z))]
    tables = {"balanced": dict(abc=[(Y + Z) / Dn, (X + Z) / Dn, (X + Y) / Dn], opts=[set()], mark=({}, 1),
                               expected=DOC_COUNTS["balanced"], kind="balanced")}
    for i, lab in enumerate("abc"):
        others = [j for j in range(3) if j != i]
        abc = [zero] * 3
        abc[i] = X
        tables["solo_" + lab] = dict(abc=abc, opts=[set()], mark=({}, 1), expected=DOC_COUNTS["solo"], kind="solo",
                                     main=i)
        for rev in (False, True):
            abc = [None] * 3
            abc[i] = Mexp
            ls = leaves[::-1] if rev else leaves
            abc[others[0]], abc[others[1]] = ls
            tables["dominant_" + lab + ("_reverse" if rev else "")] = dict(
                abc=abc, opts=[{"X"}], mark=({"X": 1}, 8), expected=DOC_COUNTS["dominant"], kind="dominant", main=i,
                rev=rev)
    ck.check("EXACT", "ten tables are set up: the balanced chart, six dominant charts (each of a, b, c dominant, "
             "both orders of the two leaves) and three solo charts", len(tables) == 10)

    def sign_ok(expr, opts, strict):
        n, d = CH.reduced(expr)
        tn, td = CH.terms(CH.poly(n)), CH.terms(CH.poly(d))
        ok = all(cf >= 0 for cf in tn.values()) and all(cf > 0 for cf in td.values()) and bool(td)
        ok = ok and all(M.support_monomial(td, CH.names, o) for o in opts)
        if strict:
            ok = ok and all(M.support_monomial(tn, CH.names, o) for o in opts)
        return ok

    bad = []
    nmargins = 0
    for label, tb in tables.items():
        a_, b_, c_ = tb["abc"]
        opts = tb["opts"]
        req = [(x, False) for x in (a_, b_, c_)]
        if tb["kind"] == "balanced":
            req += [(b_ + c_ - a_, False), (a_ + c_ - b_, False), (a_ + b_ - c_, False), (2 - a_ - b_ - c_, True)]
        elif tb["kind"] == "dominant":
            m_ = tb["abc"][tb["main"]]
            s_ = a_ + b_ + c_ - m_
            req += [(m_ - s_, False), (1 - m_ * s_, True), (s_, True)]
        else:
            m_ = tb["abc"][tb["main"]]
            okz = all(x == 0 for k, x in enumerate(tb["abc"]) if k != tb["main"])
            if not okz:
                bad.append(label)
        nmargins += len(req)
        if not all(sign_ok(e, opts, s) for e, s in req):
            bad.append(label)
    ck.check("EXACT", "each of the 10 charts maps its parameter domain into the region it is meant to cover: a, b, c "
             ">= 0 and the defining margins (balanced: the three triangle inequalities and 2 - (a+b+c) > 0; dominant: "
             "M - S >= 0, 1 - M S > 0 and S > 0 for S = sum of the two leaves; solo: both leaves 0) are reduced "
             "fractions whose numerators have coefficients >= 0, with a surviving monomial where strictness is "
             "claimed, and whose denominators have positive coefficients and a surviving monomial (%d margins); "
             "the dominant charts need X > 0, the others no strict hypothesis" % nmargins, not bad,
             "failing: %s" % bad if bad else "")
    na, nb, nc, nm, ns, nl1, nl2 = sp.symbols("a b c m s l1 l2")
    xx, yy, zz = (nb + nc - na) / 2, (na + nc - nb) / 2, (na + nb - nc) / 2
    den = 1 - xx - yy - zz
    inv = {X: xx / den, Y: yy / den, Z: zz / den}
    ok_bal = all(sp.cancel(e.subs(inv, simultaneous=True) - tg) == 0
                 for e, tg in zip(tables["balanced"]["abc"], (na, nb, nc)))
    ck.check("EXACT", "[PC 16a] inverse of the balanced chart (10): X, Y, Z = (x, y, z)/(1 - x - y - z) with x, y, z "
             "the triangle half-differences returns (a, b, c); 1 - x - y - z = 1 - (a+b+c)/2", ok_bal
             and sp.expand(den - (1 - (na + nb + nc) / 2)) == 0)
    S_ = nl1 + nl2
    invd = {X: 1 / S_ - 1, Y: (nm * S_ - S_ ** 2) / (1 - nm * S_)}
    ok_dom = sp.cancel(Mexp.subs(invd, simultaneous=True) - nm) == 0
    ok_dom &= sp.cancel(leaves[0].subs({**invd, Z: nl1 / nl2}, simultaneous=True) - nl1) == 0
    ok_dom &= sp.cancel(leaves[1].subs({**invd, Z: nl1 / nl2}, simultaneous=True) - nl2) == 0
    ok_dom &= sp.cancel(sum(leaves) - 1 / L) == 0
    ok_dom &= sp.cancel(1 - Mexp / L - (1 - 1 / L ** 2) / (1 + Y)) == 0
    ck.check("EXACT", "[PC 16b] inverse of the dominant charts (12), (13): X = 1/S - 1, Y = (M S - S^2)/(1 - M S), "
             "Z = L1/L2 returns (M, L1, L2); the two leaves add up to 1/l with l = 1 + X; "
             "1 - M/l = (1 - 1/l^2)/(1 + Y), which is 1 - M S", ok_dom)

    # ------------------------------------------------------------------------------------
    ck.section("G. The ten positive tables (document pp. 3-4)")
    ck.info("elapsed %.0f s" % (time.time() - ck.t0))
    R7 = Ring(",".join(NATIVE))
    G1p = M.p_rename(R9.subs_one(G, "d"), R9, R7)               # G at d = 1, packed, native variables
    G1 = M.packed_to_terms(R7, G1p)
    ck.check("EXACT", "G at d = 1 still has 2285 terms and degrees (3,2,3,2,3,3,2) = the clearing degrees (11)",
             len(G1) == 2285 and tuple(max(e[i] for e in G1) for i in range(7)) == DEGS)
    cert_tables = cert["tables"] if cert is not None else None
    results = {}
    subs_of = {}
    rngt = random.Random(20261189)
    t0 = time.time()
    total = 0
    for label, tb in tables.items():
        a_, b_, c_ = tb["abc"]
        subs = [Q0, a_, V, b_, W, Tt, c_]
        subs_of[label] = subs
        red = [CH.reduced(e) for e in subs]
        pairs = [(CH.poly(n), CH.poly(d)) for n, d in red]
        pterms = [(CH.terms(n), CH.terms(d)) for n, d in pairs]
        tab = CH.terms(CH.compose_cleared(G1, pairs, DEGS, order=[1, 3, 6, 0, 2, 4, 5]))
        results[label] = (tab, pterms)
        total += len(tab)
        zero7 = (0,) * 7
        den_ok = all(all(cf > 0 for cf in td.values()) and td.get(zero7, 0) > 0 for _, td in pterms)
        mk, mcf = tb["mark"]
        mex = tuple(mk.get(nm_, 0) for nm_ in CH.names)
        mark_ok = tab.get(mex) == mcf and all(all(mex[i] == 0 or CH.names[i] in o for i in range(7))
                                              for o in tb["opts"])
        ev_ok = True
        for trial in range(3):
            pt = [rngt.randint(1, 4) for _ in range(7)]
            if trial == 2:                                   # a boundary point: zeros where the chart allows them
                o = sorted(tb["opts"][0])
                pt = [(rngt.randint(1, 3) if nm_ in o else rngt.choice([0, 0, 1, 2])) for nm_ in CH.names]
            nat, scale = [], 1
            for (tn, td), m in zip(pterms, DEGS):
                dv = M.eval_terms(td, pt)
                ev_ok &= dv > 0
                nat.append(Fraction(M.eval_terms(tn, pt), dv))
                scale *= dv ** m
            gv = M.eval_terms(G1, nat)
            tv = M.eval_terms(tab, pt)
            ev_ok &= gv * scale == tv and tv > 0
        stored = ""
        dig_ok = True
        if cert_tables is not None:
            ct_ = cert_tables.get(label)
            dig_ok = ct_ is not None and ct_["terms"] == len(tab) and ct_["sha256"] == M.table_digest(tab)
            stored = "; equals the supplement's stored table (sha256 of the canonical serialisation, comparison only)"
        ck.check("EXACT", "[PC 15] table %-18s: %6d terms (document: %d), every coefficient > 0 (min %d); marker %s "
                 "present with coefficient %d and involving only variables that are > 0 on the chart %s; all "
                 "denominators have positive coefficients and a positive constant term; table = G * prod D_i^m_i "
                 "checked by exact evaluation at 3 integer parameter points%s"
                 % (label, len(tab), tb["expected"], min(tab.values()),
                    "*".join("%s^%d" % kv for kv in sorted(mk.items())) or "1", mcf,
                    [sorted(o) for o in tb["opts"]], stored),
                 len(tab) == tb["expected"] and min(tab.values()) > 0 and den_ok and mark_ok and ev_ok and dig_ok)
    ck.check("EXACT", "[PC 17] the ten tables have %d coefficients in total (document: %d), all positive"
             % (total, DOC_TOTAL), total == DOC_TOTAL, "%.1f s for all tables" % (time.time() - t0))
    Rc = Ring(",".join(CH.names))
    t0 = time.time()
    ok2 = True
    done = []
    for label in ("solo_a", "solo_b", "solo_c", "balanced", "dominant_b_reverse"):
        subs = {}
        for nm_, e in zip(NATIVE, subs_of[label]):
            n, d = CH.reduced(e)
            subs[nm_] = (Rc.from_sympy(n), Rc.from_sympy(d))
        num, _ = M.p_compose(G1p, R7, Rc, subs, degs=dict(zip(NATIVE, DEGS)))
        ok2 &= {Rc.exps(k): cf for k, cf in num.items()} == results[label][0]
        done.append(label)
    ck.check("EXACT", "second implementation (pure Python integer arithmetic, own p_compose) reproduces the tables %s "
             "term by term" % done, ok2, "%.1f s" % (time.time() - t0))

    # ------------------------------------------------------------------------------------
    ck.section("H. NUMERIC falsification attempts (evidence only)")
    ck.info("elapsed %.0f s" % (time.time() - ck.t0))
    import numpy as np
    rngn = np.random.default_rng(8801)
    worst = 0.0
    seed = 88000
    for i in range(10):
        a, b, c = [float(x) for x in np.exp(rngn.normal(size=3))]
        if i % 5 == 1:
            a = b + c                     # the balanced/dominant interface
        if i % 5 == 2:
            b = 0.0                       # a zero leaf
        if i % 5 == 3:
            c = a + b
        if i % 5 == 4:
            a = 0.0
        vals = M.combo_values({1: a, 10: b, 59: c, 30: M.numeric_H3(a, b, c)})
        worst = min(worst, M.numeric_min_phi(vals, seed + i, 20000, 24))
    ck.check("NUMERIC", "SOC sum at d = H3: smallest value of Phi found over SO(4) for a f1 + b f10 + c f59 + H3 f30 at "
             "10 coefficient draws (20000 random rotations and 24 BFGS descents each, seeds 88000-88009) is >= -1e-9",
             worst >= -1e-9, "min %.3e" % worst)
    c1 = M.numeric_min_phi(F_VALUES[30], 88900, 20000, 24)
    c2 = M.numeric_min_phi(M.combo_values({1: 1.0, 10: 1.0, 30: 2.0}), 88901, 20000, 24)
    c3 = M.numeric_min_phi(M.combo_values({10: 1.0, 59: 1.0, 30: 2.0}), 88902, 20000, 24)
    ck.check("NUMERIC", "controls with the same search budget (it does find negative values where the SOC condition "
             "fails for a single pair): Phi_{f30} alone, f1 + f10 + 2 f30, f10 + f59 + 2 f30",
             c1 < -1e-3 and c2 < -1e-3 and c3 < -1e-3, "min %.4f, %.4f, %.4f" % (c1, c2, c3))

    # G > 0 on the failure region, sampled in the NATIVE coefficients, and the chart cover
    rq = random.Random(88100)

    def rfrac(zero_p=0.25, hi=6, den=12):
        return Fraction(0) if rq.random() < zero_p else Fraction(rq.randint(1, hi), rq.randint(1, den))

    redmaps = {label: results[label][1] for label in tables}
    labs = "abc"

    def classify(a, b, c):
        """(a, b, c) with H3 < 1 -> (table label, X, Y, Z) by the document's inverse formulas."""
        m3 = (a, b, c)
        s, m = a + b + c, max(m3)
        if 2 * m <= s:
            dd = 1 - s / 2
            return "balanced", (b + c - a) / 2 / dd, (a + c - b) / 2 / dd, (a + b - c) / 2 / dd
        i = m3.index(m)
        l1, l2 = [m3[j] for j in range(3) if j != i]
        S = l1 + l2
        if S == 0:
            return "solo_" + labs[i], m, Fraction(0), Fraction(0)
        Xp, Yp = 1 / S - 1, (m * S - S * S) / (1 - m * S)
        if l2 > 0:
            return "dominant_" + labs[i], Xp, Yp, l1 / l2
        return "dominant_" + labs[i] + "_reverse", Xp, Yp, l2 / l1

    counts = {label: 0 for label in tables}
    npts = 0
    gmin = None
    ok_pos = ok_cover = True
    while npts < 4000:
        mode = rq.randrange(8)
        m3 = [rfrac(), rfrac(), rfrac()]
        if mode < 3:
            i = rq.randrange(3)
            m3[i] = sum(m3) - m3[i]                    # an interface: one coefficient equals the sum of the others
        elif mode == 3:
            i = rq.randrange(3)
            m3 = [Fraction(0)] * 3
            m3[i] = rfrac(0.2, 9, 3)                   # two mains zero (solo), any size
        elif mode == 4:
            i = rq.randrange(3)
            m3[i] = m3[i] * 6                          # a large dominant coefficient
            for j in range(3):
                if j != i:
                    m3[j] = m3[j] / 12
        a, b, c = m3
        if not M.H3_below_one(a, b, c):
            continue
        q, v, w, t = rfrac(0.3, 12, 6), rfrac(0.3, 12, 6), rfrac(0.3, 12, 6), rfrac(0.3, 12, 6)
        npts += 1
        gv = M.eval_terms(G1, (q, a, v, b, w, t, c))
        ok_pos &= gv > 0
        gmin = gv if gmin is None else min(gmin, gv)
        label, Xp, Yp, Zp = classify(a, b, c)
        counts[label] += 1
        par = dict(Q0=q, V=v, W=w, T=t, X=Xp, Y=Yp, Z=Zp)
        pv = tuple(par[nm_] for nm_ in CH.names)
        okp = all(x >= 0 for x in pv) and all(par[nm_] > 0 for nm_ in tables[label]["opts"][0])
        for (tn, td), target in zip(redmaps[label], (q, a, v, b, w, t, c)):
            dv = M.eval_terms(td, pv)
            okp = okp and dv > 0 and Fraction(M.eval_terms(tn, pv), dv) == target
        ok_cover &= okp
    ck.check("NUMERIC", "G > 0 at %d random rational points of the region {H3(a,b,c) < 1} with d = 1 and q, v, w, t "
             ">= 0, sampled in the native coefficients with zeros, the balanced/dominant interfaces, two vanishing "
             "mains and large dominant coefficients over-represented (exact evaluation, seed 88100)" % npts,
             ok_pos, "smallest value %s = %.4g" % (gmin, float(gmin)))
    ck.check("NUMERIC", "cover: each of these points is the image of an admissible parameter point of one of the 10 "
             "charts under the document's inverse formulas: parameters >= 0, X > 0 in the dominant charts, the "
             "chart returns the point exactly", ok_cover and all(n > 0 for n in counts.values()),
             "points per table: %s" % counts)

    # the necessary equations on real points of X (contact equations only)
    pts = M.numeric_equality_points(T, 88200, 60)
    fit, cub, hes, qmin, gmax = M.numeric_radial_test(rd, ORDER, pts)
    ck.check("NUMERIC", "%d generic points x = alpha Psi(a,b) of X whose real parts satisfy the five contact equations "
             "of T (seed 88200; all eight labels apart): gap vector in the span of the eight generators, "
             "C0 + hQ = 0 and G = 0 up to rounding" % len(pts),
             len(pts) >= 40 and fit < 1e-8 and cub < 1e-8 and gmax < 1e-8,
             "fit %.1e, relative |C0 + hQ| <= %.1e, relative |G| <= %.1e, min relative |Q| %.1e" % (fit, cub, gmax, qmin))

    # the conclusion on the stratum itself
    found = []
    for sd in (88300, 88301, 88302, 88303):
        found += M.numeric_stratum_search(T, sd, 100)
    absG = {k: abs(cf) for k, cf in G.items()}
    nzero = 0
    info = []
    ok_h4 = True
    for pnt in found:
        g = pnt["gaps"]
        co = dict(q=g[20], a=g[16], v=g[11], b=g[7], w=g[12], t=g[3], c=g[0], d=g[19] - g[20], h=0.0)
        fitv = float(np.max(np.abs(np.array(M.combo_values(list(zip(ORDER, [co[n] for n in R9.names[:8]]))),
                                            float) - g)))
        ok_h4 &= fitv < 1e-7
        d = co["d"]
        if abs(d) <= 1e-6:
            nzero += 1
            continue
        if d < 0:                                  # the transposed source of section C
            co = dict(q=co["v"] + d, a=co["a"], v=co["q"] + d, b=co["b"], w=co["t"] + d, t=co["w"] + d,
                      c=co["c"], d=-d, h=0.0)
            ok_h4 &= min(co["q"], co["v"], co["w"], co["t"]) > -1e-7
        relG = abs(R9.evaluate(G, co)) / R9.evaluate(absG, {k: abs(x) for k, x in co.items()})
        a, b, c = [max(co[n], 0.0) / co["d"] for n in ("a", "b", "c")]
        h3 = M.numeric_H3(a, b, c)
        ok_h4 &= relG < 1e-6 and h3 >= 1 - 1e-9
        info.append((d, relG, h3, pnt["sep"]))
    npos = sum(1 for x in info if x[0] > 0)
    ck.check("NUMERIC", "stratum points (all gaps >= 0, zero on T, sum of gaps 1) from 400 Levenberg-Marquardt starts "
             "(seeds 88300-88303): %d found, %d with |d| <= 1e-6, %d with d > 1e-6, %d with d < -1e-6; at the points "
             "with d != 0 (for d < 0 after the transpose of section C) G = 0 up to rounding (relative size "
             "< 1e-6) and the SOC condition |d| <= H3 holds" % (len(found), nzero, npos, len(info) - npos), ok_h4 and len(found) > 0,
             "; ".join("d=%.3g relG=%.1e H3/|d|=%.3g sep=%.2g" % x for x in sorted(info, reverse=True)))
    ck.finish()


if __name__ == "__main__":
    main()
