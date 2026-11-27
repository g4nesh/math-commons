"""Independent exact replay of Case T_72 = {1,3,4,7,8,12} (group V) of D. Zheng's proof of the
Marcus-de Oliveira conjecture for n = 4.

Proof document replayed: MOC4_SUPP/proofs/manin_transfer_20261021/final_referee_inventory71_uniform_20261023.pdf
(5 pages; the document calls the set "inventory 71" and works on the equivalent set
T' = {0,1,2,10,14,17}).  Every polynomial is rebuilt here from the definitions (the matrix P, the
cubic C, the vectors f_j of the paper); nothing from the supplement is executed or imported.

The script also re-establishes, with own code, every finite assertion made by the supplement's
"portable check" MOC4_SUPP/checks/fresh_T72_audit.py (read as text only).  Checks that answer one of
its assertions carry a tag [PC n]; the numbering is that of the report T72.md.

If MOC4_SUPP is set, the document's certificate JSON is read as DATA and COMPARED with the rebuilt
objects (C, Q, G9, the determinant, the hashes of the seventeen tables).  No step depends on it.

Run:   MOC4_SUPP=/path/to/supplement  python replay_T72.py
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
DOC = DIR + "final_referee_inventory71_uniform_20261023.pdf"
CERT_JSON = DIR + "final_referee_inventory71_uniform_exact_20261023.json"
PORTABLE = "checks/fresh_T72_audit.py"
SHA = {
    DOC: "5d82cd3d2c8e1fe44919a5e4c48f3177b3a7be0abbc9a8553b8c175a6e291411",
    CERT_JSON: "edb8fc9755abf4d95eb213655bec8b7350564bfd931a3bd9ad584a266a737510",
    PORTABLE: "057ea2c726d9ce2bd5133bd8566c398f92310e2c4841a79145f5729e78b98176",
}

T = (1, 3, 4, 7, 8, 12)                       # the paper's T_72
TW = (0, 1, 2, 10, 14, 17)                    # the document's working set T'
G0 = ((0, 2, 1, 3), (0, 2, 3, 1), 0)          # sigma -> rho o sigma o chi, the document's transport T' -> T
ORDER = [0, 2, 12, 13, 15, 21, 23, 26]        # g = q f0 + a f2 + b f12 + v f13 + c f15 + w f21 + t f23 + d f26
NAMES = "q,a,b,v,c,w,t,d,h"
NATIVE = ("q", "a", "b", "v", "c", "w", "t")
EXCL = (23, 15, 8, 7, 4, 16, 9)               # permutation indices dual to the first seven generators
DEGS = (3, 2, 2, 3, 2, 3, 2)                  # the document's clearing degrees (17), d = 1
PIVOT_ROWS = (0, 1, 2, 3, 4, 6, 7, 8, 9, 10, 12, 13, 16, 18)

# document, table (18) and the two counts after (19)
DOC_COUNTS = {
    "balanced": (28693, 10012, 93570), "dominant_a": (27186, 10665, 70122),
    "dominant_c_small": (36773, 13175, 70650), "dominant_c_large": (38081, 12700, 79550),
    "solo_c_small": (703, 271, 1674), "c0_cube": 1293, "c0_a_big": 1680,
}
DOC_TOTAL = 496798


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
    ck = M.Checker("Case T_72 = {1,3,4,7,8,12} (document 'inventory 71'): independent exact replay")
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
            "[PC 11] is the check 'P kappa = (1,...,1)'")
    Msub = sp.Matrix([M.P[r] for r in PIVOT_ROWS])
    ck.check("EXACT", "[PC 2] the rows 0,1,2,3,4,6,7,8,9,10,12,13,16,18 of P form a matrix of determinant 1 "
             "(so its inverse is an integer matrix of determinant 1)", Msub.det() == 1, "det = %s" % Msub.det())
    M.check_radial_necessity(ck)
    M.crosscheck_supplement_tables(ck)
    M.crosscheck_audit_library(ck, [72, 88])

    # ------------------------------------------------------------------------------------
    ck.section("B. Contact set and the transport T' -> T (document p. 1)")
    ck.info("elapsed %.0f s" % (time.time() - ck.t0))
    ck.check("EXACT", "the document's set {1,3,4,7,8,12} is the paper's T_72 (Appendix E) as a set of "
             "lexicographic indices", tuple(sorted(T)) == M.PAPER_T[72])
    gm0 = M.g_index_map(G0)
    ck.check("EXACT", "[PC 10] sigma -> rho o sigma o chi with rho = (0,2,1,3), chi = (0,2,3,1) (no transpose) maps "
             "T' = {0,1,2,10,14,17} onto T", tuple(sorted(gm0[j] for j in TW)) == T,
             "index map %s" % (list(gm0),))
    rT, rW = M.orbit_rep(T), M.orbit_rep(TW)
    ck.check("EXACT", "T and T' have the same orbit representative under G; T_88 lies in another orbit",
             rT == rW and M.orbit_rep(M.PAPER_T[88]) != rT, "T_72 -> %s, T_88 -> %s" % (rT, M.orbit_rep(M.PAPER_T[88])))
    M.check_group_on_X(ck, [("the transport g0", G0), ("its inverse", M.g_inverse(G0))])
    ck.check("EXACT", "g0 maps the f_j vanishing on T' onto the f_j vanishing on T",
             sorted(M.g_facet(gm0, j) for j in M.generators(TW)) == M.generators(T),
             "generators of T: %s" % M.generators(T))

    # ------------------------------------------------------------------------------------
    ck.section("C. The contact kernel of T' and its signed generator representation (document p. 1, eq. (1), (2))")
    ck.info("elapsed %.0f s" % (time.time() - ck.t0))
    gens = M.generators(TW)
    ck.check("EXACT", "the f_j of Appendix D vanishing on T' are exactly f0, f2, f10, f12, f13, f15, f21, f23, f26 "
             "(the document's 'nine-source list')", gens == [0, 2, 10, 12, 13, 15, 21, 23, 26],
             ", ".join("f%d = %s" % (j, M.f_name(j)) for j in gens))
    ck.check("EXACT", "[PC 3] these nine value vectors lie in the image of P with integer coordinates; f26 = h_{4132,-} "
             "is the only one that is not a minor vector",
             [j for j in gens if not M.F_IS_MINOR[j]] == [26]
             and all(x.denominator == 1 for j in gens for x in M.F_COORD[j])
             and all(M.values(M.F_COORD[j]) == list(F_VALUES[j]) for j in gens))
    rkT = M.rank([M.P[s] for s in TW])
    null = M.nullspace([M.P[s] for s in TW], 14)
    ck.check("EXACT", "[PC 4] the six rows of P at T' have rank 6 (solution space of dimension 8) and the eight vectors "
             "f0, f2, f12, f13, f15, f21, f23, f26 are linearly independent",
             rkT == 6 and len(null) == 8 and M.rank([M.F_COORD[j] for j in ORDER]) == 8)
    VL = [[F_VALUES[j][s] for j in ORDER] for s in range(24)]          # 24 x 8 value matrix
    ck.check("EXACT", "[PC 5] the eight vectors vanish at the six indices of T'",
             all(VL[s][k] == 0 for s in TW for k in range(8)))
    unit = lambda *ks: [int(k in ks) for k in range(8)]
    ck.check("EXACT", "[PC 6] rows 23, 15, 8, 7, 4, 16, 9 of the value matrix are the unit vectors e_0..e_6: for "
             "g = q f0 + a f2 + b f12 + v f13 + c f15 + w f21 + t f23 + d f26 the coefficients q, a, b, v, c, w, t "
             "are the gaps g23, g15, g8, g7, g4, g16, g9",
             [VL[s] for s in EXCL] == [unit(k) for k in range(7)])
    ck.check("EXACT", "[PC 7] row 11 minus row 9 is e_7: d = g11 - g9 (so the form (1) is unique and d may have "
             "either sign)", [x - y for x, y in zip(VL[11], VL[9])] == unit(7))
    ck.check("EXACT", "[PC 8] f10 + f26 = f0 + f13 + f23 on all 24 permutations; the nine vectors have rank 8, so "
             "this is their only linear relation up to a factor",
             [x + y for x, y in zip(F_VALUES[10], F_VALUES[26])]
             == [x + y + z for x, y, z in zip(F_VALUES[0], F_VALUES[13], F_VALUES[23])]
             and M.rank([M.F_COORD[j] for j in gens]) == 8)
    ck.check("EXACT", "[PC 9] rows 20, 12, 11 are e_0 + e_7, e_3 + e_7, e_6 + e_7: g20 = q + d, g12 = v + d, "
             "g11 = t + d", [VL[20], VL[12], VL[11]] == [unit(0, 7), unit(3, 7), unit(6, 7)])
    sq, sa, sb, sv, sc, sw, st, sd = sp.symbols("q a b v c w t d")
    lhs = M.combo_values(list(zip(ORDER, (sq, sa, sb, sv, sc, sw, st, sd))))
    rhs = M.combo_values(list(zip([0, 2, 12, 13, 15, 21, 23, 10], (sq + sd, sa, sb, sv + sd, sc, sw, st + sd, -sd))))
    ck.check("EXACT", "signed side: q f0 + a f2 + b f12 + v f13 + c f15 + w f21 + t f23 + d f26 = (q+d) f0 + a f2 + "
             "b f12 + (v+d) f13 + c f15 + w f21 + (t+d) f23 + (-d) f10 on all 24 permutations; for d <= 0 the eight "
             "coefficients on the right are the gaps g20, g15, g8, g12, g4, g16, g11 and -d, all >= 0, on minor "
             "vectors only", all(sp.expand(x - y) == 0 for x, y in zip(lhs, rhs))
             and all(M.F_IS_MINOR[j] for j in (0, 2, 12, 13, 15, 21, 23, 10)))
    dim, rays = M.extreme_rays_of_contact_cone(TW)
    ck.check("EXACT", "direct enumeration of the extreme rays of N_T' = {Pc >= 0, (Pc)_T' = 0} (dimension %d; no use "
             "of the 82-facet list) gives exactly the nine vectors above; so the cone is not simplicial, and the "
             "document's split by the sign of d covers it" % dim,
             rays == sorted(F_VALUES[j] for j in gens), "%d rays" % len(rays))
    M.check_minor_leibniz(ck, [0, 2, 10, 12, 13, 15, 21, 23])

    # ------------------------------------------------------------------------------------
    ck.section("D. The sufficient coefficient conditions (document pp. 1-2, eq. (3)-(7))")
    ck.info("elapsed %.0f s" % (time.time() - ck.t0))
    rng = random.Random(20261172)
    for tgt in ((2, 12, 26), (2, 15, 26), (12, 15, 26)):
        M.check_transport(ck, "SOC", tgt, rng)
    M.check_transport(ck, "old-H", (23, 15, 26), rng)
    M.check_transport(ck, "chain", (15, 21, 23, 26), rng)
    ck.check("EXACT", "the document's own reference tuples lie in the orbits of the audit's representatives with the "
             "same roles: old-H (f3, f15; f26) ~ (f12, f2; f70); SOC (f1, f9; f22) and chain (f2; f0, f12; f70) are "
             "the representatives themselves",
             bool(M.find_transport((12, 2, 70), (3, 15, 26))) and M.UNIVERSAL["SOC"][0] == (1, 9, 22)
             and M.UNIVERSAL["chain"][0] == (2, 0, 12, 70))
    ck.check("EXACT", "the roles of old-H matter and are the document's: (f15, f23; f26), with f15 in the t-role and "
             "f23 in the c-role, is NOT in the orbit of (f12, f2; f70); nor is a chain tuple with centre f21 or f23",
             not M.find_transport((12, 2, 70), (15, 23, 26))
             and not M.find_transport((2, 0, 12, 70), (21, 15, 23, 26), (1, 2))
             and not M.find_transport((2, 0, 12, 70), (23, 15, 21, 26), (1, 2)))
    c_, t_, w_, d_ = sp.symbols("c t w d", positive=True)
    ck.check("EXACT", "(4): the old-H condition t(c - d) >= d^2 for t f23 + c f15 + d f26 (t-role f23, c-role f15) "
             "is t(c - 1) >= 1 at d = 1; (5): the chain condition (centre f15 with coefficient c, leaves f21, f23 "
             "with coefficients w, t) is c >= d^2/w + d^2/t, and c - 1/w - 1/t = (w(ct - 1) - t)/(wt): for ct > 1 it "
             "holds iff w >= t/(ct - 1), and for ct <= 1 it never holds; this is (7)",
             sp.cancel((c_ - 1 / w_ - 1 / t_) - (w_ * (c_ * t_ - 1) - t_) / (w_ * t_)) == 0
             and sp.expand((t_ * (c_ - d_) - d_ ** 2).subs(d_, 1) - (t_ * (c_ - 1) - 1)) == 0)
    M.check_H3_allocation(ck, 2, 12, 15, 26)
    M.replay_SOC_representative(ck)

    # ------------------------------------------------------------------------------------
    ck.section("E. Radial necessity in the generators (document p. 2, eq. (8)-(10))")
    ck.info("elapsed %.0f s" % (time.time() - ck.t0))
    rd = M.radial_data(ck, ORDER, NAMES)
    R9, C0, Q, G, Ddet = rd["R"], rd["C0"], rd["Q"], rd["G"], rd["D"]
    syms9 = sp.symbols(NAMES)
    q9, a9, b9, v9, c9, w9, t9, d9, h9 = syms9
    Cdoc = R9.from_sympy((q9 + d9) * (v9 + d9) * (w9 + t9 + d9) + a9 * b9 * c9)
    Qdoc = R9.from_sympy(a9 * b9 + a9 * c9 + b9 * c9 + d9 ** 2 - t9 * (c9 - d9))
    ck.check("EXACT", "[PC 12] (8): C(g + h kappa) = C0 + h Q with C0 = (q+d)(v+d)(w+t+d) + abc and "
             "Q = ab + ac + bc + d^2 - t(c - d) (same sign as the paper's cubic)", C0 == Cdoc and Q == Qdoc)
    kd3 = R9.key((0, 0, 0, 0, 0, 0, 0, 3, 0))
    ck.check("EXACT", "all %d coefficients of C0 are +1 and C0 contains d^3: C0 >= d^3 > 0 on the closed orthant with "
             "d > 0; hence C0 + hQ = 0 forces Q != 0 and h = -C0/Q" % len(C0),
             all(v == 1 for v in C0.values()) and C0.get(kd3) == 1)
    ck.check("EXACT", "Q - (d^2 - t(c - d)) = ab + ac + bc has nonnegative coefficients: Q > 0 whenever "
             "t(c - d) < d^2 (the document's 'no Q = 0 hole')",
             M.p_nonneg(M.p_add(Q, R9.from_sympy(d9 ** 2 - t9 * (c9 - d9)), -1)))
    kd9 = R9.key((0, 0, 0, 0, 0, 0, 0, 9, 0))
    if G.get(kd9, 0) < 0:
        G = M.p_scale(G, -1)
    degs8 = [R9.deg(G, nm) for nm in R9.names[:8]]
    ck.check("EXACT", "(9): G9 has 1615 terms, is homogeneous of degree 9, does not involve h, its d^9 coefficient is "
             "+1, and its degrees in (q, a, b, v, c, w, t, d) are (3, 2, 2, 3, 2, 3, 2, 9)",
             len(G) == 1615 and all(sum(R9.exps(k)) == 9 and R9.exps(k)[8] == 0 for k in G)
             and G.get(kd9) == 1 and degs8 == [3, 2, 2, 3, 2, 3, 2, 9], "degrees %s" % degs8)
    ck.info("G9 has %d negative coefficients in the native variables, so its sign is not visible without the charts"
            % sum(1 for cf in G.values() if cf < 0))

    def swapped(p):
        out = {}
        for k, cf in p.items():
            e = R9.exps(k)
            out[R9.key((e[3], e[2], e[1], e[0]) + e[4:])] = cf
        return out
    ck.check("EXACT", "[PC 13] [PC 15] (10): C0, Q and G9 are invariant under q <-> v, a <-> b",
             swapped(C0) == C0 and swapped(Q) == Q and swapped(G) == G)
    M.check_bareiss(ck, rd, ORDER)
    ck.info("[PC 16] is answered by the two determinant computations above and the identity N = -4 Q^2 G^2, here as "
            "an identity of homogeneous polynomials in eight variables (the portable check works at d = 1)")
    cert = load_json(CERT_JSON)
    if cert is None:
        ck.info("MOC4_SUPP not set (or file absent): the certificate JSON is not compared")
    else:
        def load(rows):
            return {R9.key(tuple(int(x) for x in e)): int(cf) for e, cf in rows}
        g9rows = cert["G9"]
        fmt_ok = (len(g9rows) == 1615 and all(len(e) == 9 and e[-1] == 0 and sum(e) == 9 for e, _ in g9rows)
                  and len({tuple(e) for e, _ in g9rows}) == 1615
                  and [max(e[i] for e, _ in g9rows) for i in range(8)] == [3, 2, 2, 3, 2, 3, 2, 9])
        ck.check("EXACT", "[PC 14] supplement data (JSON, comparison only): its G9 list has 1615 distinct exponent "
                 "vectors of length 9 with last entry 0 and sum 9, maximal exponents (3,2,2,3,2,3,2,9), d^9 "
                 "coefficient 1; it EQUALS the G9 computed here; its C, Q and its stored determinant equal C0, Q "
                 "and det Hess C(g + h kappa) computed here (variables q,a,b,v,c,w,t,d,z)",
                 fmt_ok and cert["variables"] == list(R9.names[:8]) + ["z"] and load(g9rows) == G
                 and load(cert["C"]) == C0 and load(cert["Q"]) == Q and load(cert["original14_det"]) == Ddet)
        okc = cert["transport"]["action24"] == list(gm0) and cert["working_contacts"] == list(TW) \
            and cert["contacts"] == list(T) and cert["canonical_source"] == ORDER
        want = {"SOCab": (2, 12, 26), "SOCac": (2, 15, 26), "SOCbc": (12, 15, 26), "H": (23, 15, 26),
                "chain": (15, 21, 23, 26)}
        for nm, rec in cert["capacities"].items():
            ge = (tuple(rec["row"]), tuple(rec["column"]), int(bool(rec["transpose"])))
            ac = rec["action24"]
            okc = okc and list(M.g_index_map(ge)) == ac and tuple(rec["target"]) == want.get(nm)
            # the certificate's convention: the vector with value f_i(ac[j]) at j is the target vector
            okc = okc and all(tuple(F_VALUES[i][ac[j]] for j in range(24)) == F_VALUES[tg]
                              for i, tg in zip(rec["mother"], rec["target"]))
        ck.check("EXACT", "supplement data (JSON, comparison only): its transport is the index map of g0, and its "
                 "five recorded group elements (three SOC pairs, old-H, chain) are elements of G that carry the "
                 "reference tuples (f1,f9;f22), (f3,f15;f26), (f2;f0,f12;f70) to the tuples used in section D "
                 "(pull-back convention of the certificate)", okc and set(cert["capacities"]) == set(want))

    # ------------------------------------------------------------------------------------
    ck.section("F. Charts of the region where every sufficient condition fails (document pp. 3-4, eq. (11)-(16), (19))")
    ck.info("elapsed %.0f s" % (time.time() - ck.t0))
    CH = M.ChartRing("Q0,V,W,T,X,Y,Z")
    Q0, V, W, Tt, X, Y, Z = CH.syms
    Dn, L = 1 + X + Y + Z, 1 + X
    zero, one = sp.Integer(0), sp.Integer(1)
    MAIN = {      # name: ((a, b, c), alternatives for 'these variables are > 0', kind)
        "balanced": (((Y + Z) / Dn, (X + Z) / Dn, (X + Y) / Dn), [{"X"}, {"Y"}], "balanced"),
        "dominant_a": (((1 + Y * L ** 2) / (L * (1 + Y)), Z / (L * (1 + Z)), 1 / (L * (1 + Z))), [{"X"}], "dom_a"),
        "dominant_c_small": ((Z / (L * (1 + Y) * (1 + Z)), 1 / (L * (1 + Y) * (1 + Z)), 1 / L), [{"X"}, {"Y"}], "dom_c"),
        "dominant_c_large": ((Y * Z / (L * (1 + Y) * (1 + Z)), Y / (L * (1 + Y) * (1 + Z)), L), [{"X"}], "dom_c"),
        "solo_c_small": ((zero, zero, 1 / L), [set()], "dom_c"),
    }
    # document's strict markers: {table: [({variable: exponent}, coefficient), ...]}
    MARK = {
        "balanced": ([({"X": 2}, 1), ({"Y": 2}, 1)], [({}, 1)], [({"T": 3, "X": 3}, 1), ({"T": 3, "Y": 3}, 1)]),
        "dominant_a": ([({"X": 1}, 8)], [({}, 8)], [({"T": 3}, 8)]),
        "dominant_c_small": ([({"X": 1}, 8), ({"Y": 1}, 4)], [({}, 8)], [({"T": 3}, 8)]),
        "dominant_c_large": ([({}, 2)], [({}, 10)], [({"T": 5}, 3)]),
        "solo_c_small": ([({}, 2)], [({}, 10)], [({"T": 3}, 10)]),
    }
    tables = {}          # label -> dict(subs=[7 sympy expressions], opts, marks, expected, kind, branch)
    for name, (abc, opts, kind) in MAIN.items():
        a_, b_, c_e = abc
        for bi, branch in enumerate(("low", "equal", "high")):
            o2 = [set(o) for o in opts]
            if branch == "low":
                te, we = Tt / (c_e * (1 + Tt)), W
            elif branch == "equal":
                te, we = 1 / c_e, W
            else:
                o2 = [o | {"T"} for o in o2]
                tau = Tt / (X * (1 + Tt)) if name == "dominant_c_large" else Tt
                te = (1 + tau) / c_e
                we = te * W / (tau * (1 + W))
            tables[name + "_" + branch] = dict(subs=[Q0, a_, b_, V, c_e, we, te], opts=o2, marks=MARK[name][bi],
                                               expected=DOC_COUNTS[name][bi], kind=kind, branch=branch, main=name)
    tables["c0_cube"] = dict(subs=[Q0, X / (1 + X), Y / (1 + Y), V, zero, W, Tt], opts=[set()], marks=[({}, 1)],
                             expected=DOC_COUNTS["c0_cube"], kind="c0_cube", branch="c0", main="c0_cube")
    tables["c0_a_big"] = dict(subs=[Q0, 1 + X, Y / ((1 + X) * (1 + Y)), V, zero, W, Tt], opts=[set()],
                              marks=[({}, 2)], expected=DOC_COUNTS["c0_a_big"], kind="c0_a_big", branch="c0",
                              main="c0_a_big")
    ck.check("EXACT", "seventeen tables are set up: five main charts times three chain branches, and two charts of "
             "the face c = 0", len(tables) == 17)

    def sign_ok(expr, opts, strict):
        """expr = N/D reduced; N has coefficients >= 0, D has coefficients > 0 and is > 0 on the chart domain;
        if strict, N is > 0 on the chart domain (a monomial survives for each alternative)."""
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
        q_, a_, b_, v_, c_e, we, te = tb["subs"]
        opts, kind, branch = tb["opts"], tb["kind"], tb["branch"]
        req = [(x, False) for x in (a_, b_, c_e, we, te)]                 # the coefficients are >= 0
        if kind == "balanced":
            req += [(b_ + c_e - a_, False), (a_ + c_e - b_, False), (a_ + b_ - c_e, False),
                    (2 - a_ - b_ - c_e, True), (c_e, True)]
        elif kind == "dom_a":
            req += [(a_ - b_ - c_e, False), (1 - a_ * (b_ + c_e), True), (c_e, True)]
        elif kind == "dom_c":
            req += [(c_e - a_ - b_, False), (1 - c_e * (a_ + b_), True), (c_e, True)]
            req += [(c_e - 1, True)] if tb["main"] == "dominant_c_large" else [(1 - c_e, False)]
        elif kind == "c0_cube":
            req += [(1 - a_, True), (1 - b_, True), (1 - a_ * b_, True)]
        elif kind == "c0_a_big":
            req += [(a_ - 1, False), (1 - a_ * b_, True)]
        req += [(1 - te * (c_e - 1), True)]                              # old-H fails
        if branch == "low":
            req += [(1 - c_e * te, True)]
        elif branch == "high":
            req += [(c_e * te - 1, True), (te - we * (c_e * te - 1), True)]   # chain applicable and failing
        okl = all(sign_ok(e, opts, s) for e, s in req)
        if branch == "equal":
            okl = okl and sp.cancel(c_e * te - 1) == 0
        if branch == "c0":
            okl = okl and c_e == 0
        nmargins += len(req)
        if not okl:
            bad.append(label)
    ck.check("EXACT", "each of the 17 charts maps its parameter domain into the region it is meant to cover: the "
             "coefficients a, b, c, w, t are >= 0, and the defining margins (triangle inequalities and 2 - (a+b+c) > 0, "
             "or dominance and 1 - m(s-m) > 0; c > 0 resp. c = 0; 1 - t(c-1) > 0; 1 - ct > 0, ct = 1, or ct - 1 > 0 "
             "with t - w(ct-1) > 0) are reduced fractions whose numerators have coefficients >= 0 with a surviving "
             "monomial where strictness is claimed and whose denominators have positive coefficients and a "
             "surviving monomial (%d margins)" % nmargins, not bad, "failing: %s" % bad if bad else "")

    # inverses: every point of each region is the image of a point of the chart domain
    na, nb, nc, nt, nw = sp.symbols("a b c t w")

    def back(exprs, inv, targets):
        return all(sp.cancel(e.subs(inv, simultaneous=True) - tg) == 0 for e, tg in zip(exprs, targets))
    xx, yy, zz = (nb + nc - na) / 2, (na + nc - nb) / 2, (na + nb - nc) / 2
    den = 1 - xx - yy - zz
    ok_inv = back(MAIN["balanced"][0], {X: xx / den, Y: yy / den, Z: zz / den}, (na, nb, nc))
    Sa = nb + nc
    ok_inv &= back(MAIN["dominant_a"][0], {X: 1 / Sa - 1, Y: (na - Sa) / (1 / Sa - na), Z: nb / nc}, (na, nb, nc))
    Sc = na + nb
    ok_inv &= back(MAIN["dominant_c_small"][0], {X: 1 / nc - 1, Y: nc / Sc - 1, Z: na / nb}, (na, nb, nc))
    ok_inv &= back(MAIN["dominant_c_large"][0], {X: nc - 1, Y: nc * Sc / (1 - nc * Sc), Z: na / nb}, (na, nb, nc))
    ok_inv &= back(MAIN["solo_c_small"][0], {X: 1 / nc - 1}, (0, 0, nc))
    ck.check("EXACT", "[PC 18a] inverses of the main charts (11)-(14) and of the solo chart: substituting "
             "X = x/(1-x-y-z) etc. (x, y, z the triangle half-differences); X = 1/S - 1, Y = (a-S)/(1/S - a), "
             "Z = b/c (S = b+c); X = 1/c - 1, Y = c/S - 1, Z = a/b (S = a+b); X = c - 1, Y = cS/(1 - cS), Z = a/b "
             "returns (a, b, c)", ok_inv)
    ct = nc * nt
    ok_ch = sp.cancel((Tt / (nc * (1 + Tt))).subs(Tt, ct / (1 - ct)) - nt) == 0
    ok_ch &= sp.cancel(((1 + Tt) / nc).subs(Tt, ct - 1) - nt) == 0
    ok_ch &= sp.cancel((nt * W / ((ct - 1) * (1 + W))).subs(W, nw / (nt / (ct - 1) - nw)) - nw) == 0
    Xc = nc - 1
    tau_e = Tt / (Xc * (1 + Tt))
    ok_ch &= sp.cancel(tau_e.subs(Tt, Xc * (ct - 1) / (1 - Xc * (ct - 1))) - (ct - 1)) == 0
    ok_ch &= sp.cancel((1 / Xc - 1 / nc) - 1 / (Xc * nc)) == 0 and sp.cancel((nc / Xc - 1) - 1 / Xc) == 0
    ck.check("EXACT", "[PC 18b] inverses of the chain interface (15), (16): T = ct/(1-ct) gives t for ct < 1; "
             "T = ct - 1 gives t for ct > 1; W = w/(t/(ct-1) - w) gives w; for c > 1 (X = c - 1), "
             "T = X(ct-1)/(1 - X(ct-1)) gives tau = ct - 1, and tau ranges over (0, 1/X) = (0, c/X - 1), i.e. "
             "1/c < t < 1/X", ok_ch)
    ok_c0 = back([X / (1 + X), Y / (1 + Y)], {X: na / (1 - na), Y: nb / (1 - nb)}, (na, nb))
    ok_c0 &= back([1 + X, Y / ((1 + X) * (1 + Y))], {X: na - 1, Y: na * nb / (1 - na * nb)}, (na, nb))
    ck.check("EXACT", "inverses of the two charts (19) of the face c = 0: X = a/(1-a), Y = b/(1-b) for a, b < 1; "
             "X = a - 1, Y = ab/(1-ab) for a >= 1, ab < 1", ok_c0)

    # ------------------------------------------------------------------------------------
    ck.section("G. The seventeen positive tables (document p. 4, eq. (17)-(19), and p. 5)")
    ck.info("elapsed %.0f s" % (time.time() - ck.t0))
    R7 = Ring(",".join(NATIVE))
    G1p = M.p_rename(R9.subs_one(G, "d"), R9, R7)               # G9 at d = 1, packed, native variables
    G1 = M.packed_to_terms(R7, G1p)                             # {exponent tuple: coeff}
    ck.check("EXACT", "G9 at d = 1 still has 1615 terms and degrees (3,2,2,3,2,3,2) = the clearing degrees (17)",
             len(G1) == 1615 and tuple(max(e[i] for e in G1) for i in range(7)) == DEGS)
    cert_tables = cert["tables"] if cert is not None else None
    results = {}
    rngt = random.Random(20261173)
    t0 = time.time()
    total = 0
    for label, tb in tables.items():
        red = [CH.reduced(e) for e in tb["subs"]]
        pairs = [(CH.poly(n), CH.poly(d)) for n, d in red]
        pterms = [(CH.terms(n), CH.terms(d)) for n, d in pairs]
        tab = CH.terms(CH.compose_cleared(G1, pairs, DEGS, order=[5, 6, 1, 2, 4, 0, 3]))
        results[label] = (tab, pterms)
        total += len(tab)
        den_ok = all(all(cf > 0 for cf in td.values()) and all(M.support_monomial(td, CH.names, o) for o in tb["opts"])
                     for _, td in pterms)
        marks = [(tuple(mk.get(nm, 0) for nm in CH.names), cf) for mk, cf in tb["marks"]]
        mark_ok = all(tab.get(ex) == cf for ex, cf in marks) and \
            all(any(all(ex[i] == 0 or CH.names[i] in o for i in range(7)) for ex, _ in marks) for o in tb["opts"])
        # exact evaluation at random integer parameter points of the chart domain: table = G9 * prod D^m
        ev_ok = True
        for trial in range(3):
            lo = 1 if trial < 2 else 0
            pt = [rngt.randint(lo, 4) for _ in range(7)]
            if trial == 2:                                   # a boundary point: zeros where the chart allows them
                o = sorted(tb["opts"][0])
                pt = [(rngt.randint(1, 3) if nm in o else rngt.choice([0, 0, 1, 2])) for nm in CH.names]
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
        ck.check("EXACT", "[PC 17] table %-22s: %6d terms (document: %d), every coefficient > 0 (min %d); markers %s "
                 "present with these coefficients, one of them surviving in each case of the strict hypothesis %s; "
                 "all seven denominators have positive coefficients and are > 0 on the chart domain; "
                 "table = G9 * prod D_i^m_i checked by exact evaluation at 3 integer parameter points%s"
                 % (label, len(tab), tb["expected"], min(tab.values()),
                    [("*".join("%s^%d" % kv for kv in sorted(mk.items())) or "1", cf) for mk, cf in tb["marks"]],
                    [sorted(o) for o in tb["opts"]], stored),
                 len(tab) == tb["expected"] and min(tab.values()) > 0 and den_ok and mark_ok and ev_ok and dig_ok)
    ck.check("EXACT", "[PC 19] the seventeen tables have %d coefficients in total (document: %d), all positive"
             % (total, DOC_TOTAL), total == DOC_TOTAL, "%.1f s for all tables" % (time.time() - t0))
    # second implementation (pure Python integers, p_compose) for six of the tables
    Rc = Ring(",".join(CH.names))
    t0 = time.time()
    ok2 = True
    done = []
    for label in ("solo_c_small_low", "solo_c_small_equal", "solo_c_small_high", "c0_cube", "c0_a_big",
                  "balanced_equal"):
        tb = tables[label]
        subs = {}
        for nm, e in zip(NATIVE, tb["subs"]):
            n, d = CH.reduced(e)
            subs[nm] = (Rc.from_sympy(n), Rc.from_sympy(d))
        num, _ = M.p_compose(G1p, R7, Rc, subs, degs=dict(zip(NATIVE, DEGS)))
        ok2 &= {Rc.exps(k): cf for k, cf in num.items()} == results[label][0]
        done.append(label)
    ck.check("EXACT", "second implementation (pure Python integer arithmetic, own p_compose) reproduces the tables %s "
             "term by term" % done, ok2, "%.1f s" % (time.time() - t0))

    # ------------------------------------------------------------------------------------
    ck.section("H. NUMERIC falsification attempts (evidence only)")
    ck.info("elapsed %.0f s" % (time.time() - ck.t0))
    # H1. the cited inequalities at their thresholds
    import numpy as np
    rngn = np.random.default_rng(7201)
    worst = {"SOC sum at d = H3": 0.0, "old-H at t(c-1) = 1": 0.0, "chain at c = 1/w + 1/t": 0.0}
    seed = 72000
    for i in range(8):
        a, b, c = [float(x) for x in np.exp(rngn.normal(size=3))]
        if i % 4 == 1:
            a = b + c                     # the balanced/dominant interface
        if i % 4 == 2:
            b = 0.0                       # a zero leaf
        if i % 4 == 3:
            a = b = 0.0                   # H3 = 0
        h3 = M.numeric_H3(a, b, c)
        vals = M.combo_values({2: a, 12: b, 15: c, 26: h3})
        worst["SOC sum at d = H3"] = min(worst["SOC sum at d = H3"], M.numeric_min_phi(vals, seed, 20000, 24))
        t, w = [float(x) for x in np.exp(rngn.normal(size=2))]
        vals = M.combo_values({23: t, 15: 1 + 1 / t, 26: 1.0})
        worst["old-H at t(c-1) = 1"] = min(worst["old-H at t(c-1) = 1"], M.numeric_min_phi(vals, seed + 1, 20000, 24))
        vals = M.combo_values({15: 1 / w + 1 / t, 21: w, 23: t, 26: 1.0})
        worst["chain at c = 1/w + 1/t"] = min(worst["chain at c = 1/w + 1/t"],
                                              M.numeric_min_phi(vals, seed + 2, 20000, 24))
        seed += 10
    for k, v in worst.items():
        ck.check("NUMERIC", "%s: smallest value of Phi found over SO(4) for 8 coefficient draws (20000 random "
                 "rotations and 24 BFGS descents each, seeds 72000..) is >= -1e-9" % k, v >= -1e-9, "min %.3e" % v)
    c1 = M.numeric_min_phi(F_VALUES[26], 72900, 20000, 24)
    c2 = M.numeric_min_phi(M.combo_values({2: 1.0, 12: 1.0, 26: 2.0}), 72901, 20000, 24)
    c3 = M.numeric_min_phi(M.combo_values({23: 1.0, 15: 1.0, 26: 1.0}), 72902, 20000, 24)
    ck.check("NUMERIC", "controls with the same search budget (it does find negative values in simple cases where no "
             "condition holds): Phi_{f26} alone, f2 + f12 + 2 f26 (SOC pair, ab < d^2), f23 + f15 + f26 (old-H pair, "
             "t(c-1) < 1)",
             c1 < -1e-3 and c2 < -1e-3 and c3 < -1e-3, "min %.4f, %.4f, %.4f" % (c1, c2, c3))

    # H2. G9 > 0 on the failure region, sampled in the NATIVE coefficients, and the chart cover
    rq = random.Random(72100)

    def rfrac(zero_p=0.25, hi=12, den=6):
        return Fraction(0) if rq.random() < zero_p else Fraction(rq.randint(1, hi), rq.randint(1, den))

    def failure(a, b, c, w, t):
        if not M.H3_below_one(a, b, c):
            return False
        if not t * (c - 1) < 1:
            return False
        return c * t <= 1 or w * (c * t - 1) < t

    redmaps = {label: [(CH.terms(CH.poly(n)), CH.terms(CH.poly(d))) for n, d in
                       (CH.reduced(e) for e in tb["subs"])] for label, tb in tables.items()}

    def classify(q, a, b, v, c, w, t):
        """native point of the failure region -> (table label, chart parameters, the point after the swap (10))."""
        def swp():
            nonlocal q, a, b, v
            q, v, a, b = v, q, b, a
        Xp = Yp = Zp = Tp = Fraction(0)
        if c == 0:
            if a < 1 and b < 1:
                label, Xp, Yp = "c0_cube", a / (1 - a), b / (1 - b)
            else:
                if a < 1:
                    swp()
                label, Xp, Yp = "c0_a_big", a - 1, a * b / (1 - a * b)
            return label, dict(Q0=q, V=v, W=w, T=t, X=Xp, Y=Yp, Z=Zp), (q, a, b, v, c, w, t)
        s, m = a + b + c, max(a, b, c)
        if 2 * m <= s:
            main_, dd = "balanced", 1 - s / 2
            Xp, Yp, Zp = (b + c - a) / 2 / dd, (a + c - b) / 2 / dd, (a + b - c) / 2 / dd
        elif c == m:
            S = a + b
            if c <= 1:
                if S > 0:
                    if b == 0:
                        swp()
                    main_, Xp, Yp, Zp = "dominant_c_small", 1 / c - 1, c / S - 1, a / b
                else:
                    main_, Xp = "solo_c_small", 1 / c - 1
            else:
                main_, Xp = "dominant_c_large", c - 1
                if S > 0:
                    if b == 0:
                        swp()
                    Yp, Zp = c * S / (1 - c * S), a / b
        else:
            if b == m:
                swp()
            S = b + c
            main_, Xp, Yp, Zp = "dominant_a", 1 / S - 1, (a - S) / (1 / S - a), b / c
        ctv = c * t
        if ctv < 1:
            br, Tp, Wp = "low", ctv / (1 - ctv), w
        elif ctv == 1:
            br, Wp = "equal", w
        else:
            br, Wp = "high", w / (t / (ctv - 1) - w)
            Tp = ctv - 1 if main_ != "dominant_c_large" else Xp * (ctv - 1) / (1 - Xp * (ctv - 1))
        return main_ + "_" + br, dict(Q0=q, V=v, W=Wp, T=Tp, X=Xp, Y=Yp, Z=Zp), (q, a, b, v, c, w, t)

    counts = {label: 0 for label in tables}
    npts = 0
    gmin = None
    ok_pos = ok_cover = True
    while npts < 4000:
        mode = rq.randrange(8)
        a, b, c = rfrac(0.25, 6, 12), rfrac(0.25, 6, 12), rfrac(0.25, 6, 12)
        if mode == 0:
            c = Fraction(0)
        elif mode == 1:
            c = Fraction(1)
        elif mode == 2:
            a = b + c                                  # balanced / dominant interface
        elif mode == 3:
            c = a + b                                  # the other interface
        elif mode == 4:
            c = 1 + rfrac(0.0)                         # c > 1
            a, b = a / 20, b / 20
        sc_ = Fraction(rq.randint(1, 8), 8)
        if mode not in (1, 4):
            a, b, c = a * sc_, b * sc_, c * sc_
        t = rfrac()
        if c > 0 and rq.random() < 0.25:
            t = 1 / c                                  # the interface ct = 1
        w = rfrac()
        if c * t > 1 and rq.random() < 0.5:
            w = t / (c * t - 1) * Fraction(rq.randint(0, 15), 16)
        q, v = rfrac(0.3), rfrac(0.3)
        if not failure(a, b, c, w, t):
            continue
        npts += 1
        gv = M.eval_terms(G1, (q, a, b, v, c, w, t))
        ok_pos &= gv > 0
        gmin = gv if gmin is None else min(gmin, gv)
        label, par, pt = classify(q, a, b, v, c, w, t)
        counts[label] += 1
        pv = tuple(par[nm] for nm in CH.names)
        okp = all(x >= 0 for x in pv) and any(all(par[nm] > 0 for nm in o) for o in tables[label]["opts"])
        for (tn, td), target in zip(redmaps[label], pt):
            dv = M.eval_terms(td, pv)
            okp = okp and dv > 0 and Fraction(M.eval_terms(tn, pv), dv) == target
        okp = okp and M.eval_terms(G1, pt) == gv        # the swap (10) does not change G9
        ok_cover &= okp
    ck.check("NUMERIC", "G9 > 0 at %d random rational points of the failure region {H3 < 1, t(c-1) < 1, ct <= 1 or "
             "w < t/(ct-1)} with d = 1, sampled in the native coefficients with zeros and the interfaces c = 0, "
             "c = 1, c > 1, a = b + c, c = a + b, ct = 1 over-represented (exact evaluation, seed 72100)" % npts,
             ok_pos, "smallest value %s = %.4g" % (gmin, float(gmin)))
    ck.check("NUMERIC", "cover: each of these points is the image of an admissible parameter point of one of the 17 "
             "charts under the document's inverse formulas (after the swap (10) where the document uses it): "
             "parameters >= 0, strict hypotheses met, the chart returns the point exactly",
             ok_cover and all(n > 0 for n in counts.values()), "points per table: %s" % counts)

    # H3. the necessary equations on real points of X (contact equations only)
    pts = M.numeric_equality_points(TW, 72200, 60)
    fit, cub, hes, qmin, gmax = M.numeric_radial_test(rd, ORDER, pts)
    ck.check("NUMERIC", "%d generic points x = alpha Psi(a,b) of X whose real parts satisfy the five contact equations "
             "of T' (seed 72200; all eight labels apart): gap vector in the span of the eight generators, "
             "C0 + hQ = 0 and G9 = 0 up to rounding" % len(pts),
             len(pts) >= 40 and fit < 1e-8 and cub < 1e-8 and gmax < 1e-8,
             "fit %.1e, relative |C0 + hQ| <= %.1e, relative |G9| <= %.1e, min relative |Q| %.1e" % (fit, cub, gmax, qmin))

    # H4. the conclusion on the stratum itself
    found = []
    for sd in (72300, 72301, 72302, 72303):
        found += M.numeric_stratum_search(TW, sd, 100)
    absG = {k: abs(cf) for k, cf in G.items()}
    npos = nneg = 0
    info = []
    ok_h4 = True
    for pnt in found:
        g = pnt["gaps"]
        co = dict(q=g[23], a=g[15], b=g[8], v=g[7], c=g[4], w=g[16], t=g[9], d=g[11] - g[9], h=0.0)
        fitv = float(np.max(np.abs(np.array(M.combo_values(list(zip(ORDER, [co[n] for n in R9.names[:8]]))),
                                            float) - g)))
        ok_h4 &= fitv < 1e-7
        if co["d"] <= 1e-6:
            nneg += 1
            continue
        npos += 1
        d = co["d"]
        relG = abs(R9.evaluate(G, co)) / R9.evaluate(absG, {k: abs(x) for k, x in co.items()})
        a, b, c, w, t = [max(co[n], 0.0) / d for n in ("a", "b", "c", "w", "t")]
        h3 = M.numeric_H3(a, b, c)
        soc = h3 >= 1 - 1e-9
        oldh = t * (c - 1) >= 1 - 1e-9
        chain = w > 0 and t > 0 and c >= 1 / w + 1 / t - 1e-9
        ok_h4 &= relG < 1e-6 and (soc or oldh or chain)
        info.append((d, relG, h3, soc, oldh, chain, pnt["sep"]))
    ck.check("NUMERIC", "stratum points (all gaps >= 0, zero on T', sum of gaps 1) from 400 Levenberg-Marquardt "
             "starts (seeds 72300-72303): %d found, %d with d <= 1e-6 (minor vectors only), %d with d > 1e-6; at "
             "the latter G9 = 0 up to rounding (relative size < 1e-6) and at least one of the sufficient conditions "
             "(3), (4), (5) holds"
             % (len(found), nneg, npos), ok_h4 and len(found) > 0,
             "; ".join("d=%.3g relG=%.1e H3/d=%.3g SOC=%s oldH=%s chain=%s sep=%.2g" % x for x in sorted(info, reverse=True)))
    ck.finish()


if __name__ == "__main__":
    main()
