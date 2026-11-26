"""Certificate replay (auditor cert-T39) for independent_moc4_T39_high_j_uniform_20261022.

Imports the shared CORE (own code) from the high_p replay in this directory.  Then replays
the eight high-j certificates: four under the chart (5) for j > 1 and four on the wall j = 1,
including the explicit base part.  A final section checks the T_40 cone data used by the
assembly document (generators, the three rewriting identities, the transpose symmetry) and
the four-region cover.  Supplement JSON is read as DATA only.

usage: MOC4_SUPP=<supplement root> python -I -B this_file.py
"""
import importlib.util
import os
import sys
import time

import flint

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location(
    "hp", os.path.join(HERE, "independent_moc4_T39_high_p_uniform_20261022__certificates.py"))
hp = importlib.util.module_from_spec(spec)
spec.loader.exec_module(hp)
check, all_negative, all_nonneg = hp.check, hp.all_negative, hp.all_nonneg
L = hp.L
q, p, j, c, e, b, a = hp.q, hp.p, hp.j, hp.c, hp.e, hp.b, hp.a
IQ, IP, IJ, IC, IE, IB, IA, IH = range(8)

HCTX = flint.fmpz_mpoly_ctx.get(["J", "P", "V", "L", "M", "q", "b"])
hJ, hP, hV, hL, hM, hq, hb = HCTX.gens()
WCTX = flint.fmpz_mpoly_ctx.get(["q", "P", "e", "b", "a"])
wq, wP, we, wb, wa = WCTX.gens()


def highj_clear(pol):
    """Chart (5): j = 1+J, p = P/(j(1+P)), S = V/(J(1+V)), a = S L/(L+M), e = S M/(L+M).
    Multiplier [j(1+P)]^dp [J(1+V)]^D (L+M)^D with D = max total degree in (a, e)."""
    d = pol.to_dict()
    assert all(ex[IC] == 0 and ex[IH] == 0 for ex in d)
    D = max(int(ex[IE] + ex[IA]) for ex in d)
    dp = hp.deg(pol, IP)
    one = HCTX.constant(1)
    jj = one + hJ
    out = HCTX.constant(0)
    for ex, co in d.items():
        s = int(ex[IE] + ex[IA])
        t = (HCTX.constant(int(co)) * hq ** int(ex[IQ]) * hb ** int(ex[IB]) * jj ** int(ex[IJ])
             * hP ** int(ex[IP]) * (jj * (one + hP)) ** (dp - int(ex[IP]))
             * hV ** s * (hJ * (one + hV)) ** (D - s)
             * hL ** int(ex[IA]) * hM ** int(ex[IE]) * (hL + hM) ** (D - s))
        out += t
    return out, D, dp


def wall_clear(pol):
    """j = 1, p = P/(1+P); multiplier (1+P)^dp."""
    d = pol.to_dict()
    dp = hp.deg(pol, IP)
    one = WCTX.constant(1)
    out = WCTX.constant(0)
    for ex, co in d.items():
        assert ex[IC] == 0 and ex[IH] == 0
        out += (WCTX.constant(int(co)) * wq ** int(ex[IQ]) * wP ** int(ex[IP]) * (one + wP) ** (dp - int(ex[IP]))
                * we ** int(ex[IE]) * wb ** int(ex[IB]) * wa ** int(ex[IA]))   # j -> 1
    return out


def cone_and_cover():
    print("---- T_40 cone data used by the assembly document, and the four-region cover ----")
    T40 = [0, 1, 5, 9, 14, 18]
    gens = L.generators(T40)
    check("generators of N_{T_40} (Lemma 4.2 via library) = {1,2,5,9,11,12,28,29,33,35}",
          gens == [1, 2, 5, 9, 11, 12, 28, 29, 33, 35], str(gens))
    check("non-minor generators are exactly f28, f35",
          [g for g in gens if not L.F_IS_MINOR[g]] == [28, 35])
    V = L.F_VALUES
    add = lambda *js: [sum(V[k][s] for k in js) for s in range(24)]
    check("identity f28 + f12 = f1 + f2 + f11 + f29 (all 24 values)", add(28, 12) == add(1, 2, 11, 29))
    check("identity f35 + f11 = f1 + f5 + f12 + f33", add(35, 11) == add(1, 5, 12, 33))
    check("identity f28 + f35 = 2 f1 + f2 + f5 + f29 + f33", add(28, 35) == add(1, 1, 2, 5, 29, 33))
    g = (L.PERMS[0], L.PERMS[0], 1)      # sigma -> sigma^{-1}
    check("inversion sigma -> sigma^{-1} maps T_40 onto itself",
          L.g_act_set(g, T40) == frozenset(T40))
    pairs = {2: 5, 1: 1, 11: 12, 9: 9, 5: 2, 28: 35, 29: 33, 33: 29}
    ok = all(L.g_act_values(g, list(V[k])) == list(V[pairs[k]]) for k in pairs)
    check("inversion maps g28-source (f2,f1,f11,f9,f5,f28,f29,f33) to g35-source (f5,f1,f12,f9,f2,f35,f33,f29) termwise", ok)
    # exhaustive cover of the closed nonnegative cone (d = 1) by the four regions, checked
    # on a fine rational grid including all walls, plus the logical proof in the report.
    from fractions import Fraction as Fr
    grid = [Fr(k, 4) for k in range(0, 13)] + [Fr(1, 3), Fr(3, 2) + Fr(1, 100), Fr(99, 100), Fr(101, 100)]
    miss = 0
    for pv in grid:
        for jv in grid:
            for S in grid:          # S = a + e
                soc = pv * jv >= 1
                oldH = S > 0 and jv >= 1 + 1 / S
                hj = jv >= 1 and pv * jv < 1 and (jv - 1) * S < 1
                lo = pv < 1 and jv < 1
                hpch = pv >= 1 and pv * jv < 1
                if not (soc or oldH or hj or lo or hpch):
                    miss += 1
    check("four-region cover: every grid point (p, j, a+e) lies in SOC, old-H, high-j, low-p,j or high-p", miss == 0)


def main():
    R = hp.core(full=True)
    A_, B_, N_, Nj, Q, k0, k1, k2, Fa, Fb = (R[k] for k in
                                            ["A", "B", "N", "Nj", "Q", "gamma", "beta", "alpha", "Fa", "Fb"])
    S = a + e
    J = j - 1
    check("N - (1 + S + ea) has nonnegative coefficients (S = a + e)", all_nonneg(N_ - (1 + S + e * a)))
    check("-Q = (q+b+1+S-p)c + 1 - JS + (q+b)(1+S) + pj (J = j - 1)",
          -Q == (q + b + 1 + S - p) * c + 1 - J * S + (q + b) * (1 + S) + p * j)
    data = hp.load_supp_json("independent_moc4_T39_high_j_uniform_exact_20261022.json")
    stored = {cc["name"]: cc["coefficients"] for cc in data["strict_certificates"]} if data else {}
    polys = [("gamma", k0, 3648, 228),
             ("AB_beta_minus_Nj_alpha", A_ * B_ * k1 - Nj * k2, 14147, 491),
             ("AB_Fa_plus_N_alpha", A_ * B_ * Fa + N_ * k2, 14123, 479),
             ("AB_Fb_plus_N_alpha", A_ * B_ * Fb + N_ * k2, 14123, 479)]
    outs = {}
    for name, pol, want, wantw in polys:
        out, D, dp = highj_clear(pol)
        base = HCTX.from_dict({ex: co for ex, co in out.to_dict().items() if ex[2] == 0 and ex[5] == 0 and ex[6] == 0})
        jj = 1 + hJ
        expected = -jj * (jj + 1) * (jj * (1 + hP) + hP) * hJ ** D * (hL + hM) ** D
        check("high-j certificate %s: D=%d, deg_p=%d; %d terms (claimed %d), all coefficients < 0; "
              "q=b=V=0 part = -j(j+1)[j(1+P)+P] J^D (L+M)^D" % (name, D, dp, len(out), want),
              len(out) == want and all_negative(out) and base == expected and dp == 2)
        w = wall_clear(pol)
        cst = int(w.to_dict().get((0,) * 5, 0))
        check("j=1 wall certificate %s: %d terms (claimed %d), all coefficients < 0, constant %d (claimed -2)"
              % (name, len(w), wantw, cst),
              len(w) == wantw and all_negative(w) and cst == -2)
        outs["high_j_" + name] = out
        outs["j_equals_one_" + name] = w
    if data:
        print("INFO stored certificate names: " + ", ".join(stored))
        for name in ["gamma", "AB_beta_minus_Nj_alpha"]:
            for pre in ["high_j_", "j_equals_one_"]:
                hp.compare_with_data(pre + name, outs[pre + name], None, stored.get(pre + name))
        for pre in ["high_j_", "j_equals_one_"]:
            hp.compare_factor_pair(pre + "F+- certificates",
                                   [outs[pre + "AB_Fa_plus_N_alpha"], outs[pre + "AB_Fb_plus_N_alpha"]],
                                   [stored.get(pre + "AB_Fa_plus_N_alpha"), stored.get(pre + "AB_Fb_plus_N_alpha")])
        hp.compare_factor_pair("octics F+, F-", [Fa, Fb], [data.get("Fa_terms"), data.get("Fb_terms")])
        hp.compare_factor_pair("factors L, M of r1", [R["L"], R["M"]],
                               [data.get("U16_terms"), data.get("V16_terms")])
    # chart (5) inverse, checked symbolically on the rational function level
    import sympy as sp
    jv, pv, av, ev = sp.symbols("j p a e", positive=True)
    Jv = jv - 1
    Sv = av + ev
    Pv = pv * jv / (1 - pv * jv)
    Vv = Jv * Sv / (1 - Jv * Sv)
    p_back = Pv / (jv * (1 + Pv))
    S_back = Vv / (Jv * (1 + Vv))
    a_back = S_back * av / (av + ev)
    e_back = S_back * ev / (av + ev)
    check("chart (5) inverse: P = pj/(1-pj), V = JS/(1-JS), L:M = a:e reproduce (p, S, a, e)",
          all(sp.simplify(x) == 0 for x in [p_back - pv, S_back - Sv, a_back - av, e_back - ev]))
    cone_and_cover()
    print("SUMMARY %d failures; elapsed %.1f s" % (len(hp.FAILS), time.monotonic() - hp.T0))
    return 1 if hp.FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
