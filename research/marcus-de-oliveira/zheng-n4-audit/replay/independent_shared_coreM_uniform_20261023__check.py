"""Independent replay for independent_shared_coreM_uniform_20261023 (seven-contact set M).

Own code only. Nothing from the supplement is executed, imported or read.
Builds everything from the shared audit library (../lib/moc4lib.py: the paper's P, f_j,
cubic C of (A.2), group G) and python-flint.

Usage:  python -I -B independent_shared_coreM_uniform_20261023__check.py [--numeric]
Prints one line per check; exits 1 if an EXACT check fails.
"""
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import moc4lib as L  # noqa: E402
import flint  # noqa: E402

FAIL = []


def check(name, ok, info=""):
    print(f"[{'PASS' if ok else 'FAIL'}] EXACT {name} {info}".rstrip(), flush=True)
    if not ok:
        FAIL.append(name)


T_M = [0, 1, 2, 7, 10, 16, 17]
NAMES = ("q", "p", "t", "b", "c", "v")       # then d
FIDX = [0, 2, 3, 12, 15, 23, 26]             # coefficients q,p,t,b,c,v,d


def group_image_maps():
    idx = {vv: j for j, vv in enumerate(L.F_VALUES)}
    maps = []
    for mp in L._G_MAPS:
        img = []
        for j in range(82):
            new = [None] * 24
            for s in range(24):
                new[mp[s]] = L.F_VALUES[j][s]
            img.append(idx[tuple(new)])
        maps.append(img)
    return maps


def structural():
    check("M is Table 5 set #13", L.AUX7[13] == frozenset(T_M))
    gens = L.generators(T_M)
    check("generators(M) = [0,2,3,8,12,15,23,26] (the document's zero-active list), only f26 non-minor",
          gens == [0, 2, 3, 8, 12, 15, 23, 26] and [j for j in gens if not L.F_IS_MINOR[j]] == [26])
    check("contact rows of P on M have rank 7", L.rank([L.P[s] for s in T_M]) == 7)
    check("basis (1) independent", L.rank([L.F_VECS[j] for j in FIDX]) == 7)
    V = L.F_VALUES
    check("f8 + f26 = f2 + f3 + f12 + f23", all(V[8][s] + V[26][s] == V[2][s] + V[3][s] + V[12][s] + V[23][s]
                                                for s in range(24)))
    ctx = flint.fmpz_mpoly_ctx.get(NAMES + ("d",), "deglex")
    q, p, t, b, c, v, d = ctx.gens()
    co = [q, p, t, b, c, v, d]
    g = [sum((co[k] * V[j][s] for k, j in enumerate(FIDX)), 0 * q) for s in range(24)]
    check("cost vanishes on M", all(g[s] == 0 for s in T_M))
    check("(q,p,t,b,c,v) = (g23,g15,g14,g8,g4,g9) and d = g6 - g8",
          [g[23], g[15], g[14], g[8], g[4], g[9], g[6] - g[8]] == [q, p, t, b, c, v, d])
    check("(2): g13=p+d, g12=t+d, g6=b+d, g3=v+d", [g[13], g[12], g[6], g[3]] == [p + d, t + d, b + d, v + d])
    alt = [sum((q * V[0][s], (p + d) * V[2][s], (t + d) * V[3][s], (b + d) * V[12][s], c * V[15][s],
                (v + d) * V[23][s], -d * V[8][s]), 0 * q) for s in range(24)]
    check("d<=0 rewrite as a minor combination (f8 minor)", alt == g and L.F_IS_MINOR[8])
    maps = group_image_maps()

    def exists(ref, tgt, swap=None):
        for img in maps:
            r = [img[j] for j in ref]
            if r == list(tgt):
                return True
            if swap:
                i, k = swap
                r[i], r[k] = r[k], r[i]
                if r == list(tgt):
                    return True
        return False
    for tgt in [(2, 12, 26), (2, 15, 26), (12, 15, 26)]:
        check(f"some g in G maps SOC reference (f1,f9;f22) to {tgt}", exists((1, 9, 22), tgt, (0, 1)))
    check("old-H (f3,f15;f26) maps to (f23,f15;f26)", exists((3, 15, 26), (23, 15, 26)))
    check("chain (center f2; leaves f0,f12; f70) maps to (center f15; leaves f0,f3; f26)",
          exists((2, 0, 12, 70), (15, 0, 3, 26), (1, 2)))


def hessian_forms():
    import sympy
    ms, Cs = L.cubic_sympy()
    H = [[None] * 14 for _ in range(14)]
    for i in range(14):
        for j in range(14):
            ex = sympy.expand(sympy.diff(Cs, ms[i], ms[j]))
            dct = {}
            if ex != 0:
                for mon, co in sympy.Poly(ex, *ms).terms():
                    dct[mon.index(1)] = int(co)
            H[i][j] = dct
    return H


def bareiss(Mx, one):
    A = [r[:] for r in Mx]
    n = len(A)
    sign, prev, nd = 1, one, 0
    for k in range(n - 1):
        if A[k][k].is_zero():
            sw = next((r for r in range(k + 1, n) if not A[r][k].is_zero()), None)
            if sw is None:
                return 0 * one, nd
            A[k], A[sw] = A[sw], A[k]
            sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                A[i][j] = (A[i][j] * A[k][k] - A[i][k] * A[k][j]) / prev
                nd += 1
        prev = A[k][k]
    return sign * A[n - 1][n - 1], nd


def cubic_and_G():
    ctx = flint.fmpq_mpoly_ctx.get(NAMES + ("d", "z"), "deglex")
    q, p, t, b, c, v, d, z = ctx.gens()
    co = [q, p, t, b, c, v, d]
    w = [sum((co[k] * flint.fmpq(L.F_VECS[j][i].numerator, L.F_VECS[j][i].denominator)
              for k, j in enumerate(FIDX)), 0 * q) + (z if i == 13 else 0) for i in range(14)]
    F = L.cubic_C(w)
    C = d * (d + v) * (d + q + t) + c * (p * b - t * v)
    Q = p * b + c * (p + b) + d * d - (t + v) * (c - d)
    S = t + v
    check("(5): paper cubic C(w(z)) = C + Q z", F == C + Q * z)
    check("Q = pb + c(p+b) + [d^2 - S(c-d)]", Q == p * b + c * (p + b) + (d * d - S * (c - d)))
    check("(7, first): S C = S d(d+v)(d+q) + d^2 t^2 + S c p b + t v [d^2 - S(c-d)]",
          S * C == S * d * (d + v) * (d + q) + d * d * t * t + S * c * p * b + t * v * (d * d - S * (c - d)))
    check("C at S=0 equals d^2(d+q) + cpb",
          C.subs({"t": 0, "v": 0}) == (d * d * (d + q) + c * p * b).subs({"t": 0, "v": 0}))
    t0 = time.time()
    cz = flint.fmpz_mpoly_ctx.get(NAMES + ("z",), "deglex")
    G = cz.gens()
    z = G[-1]
    one = cz.from_dict({(0,) * 7: 1})
    coef = list(G[:6]) + [one]
    w2 = []
    for i in range(14):
        s = 0 * z
        for k, j in enumerate(FIDX):
            fr = L.F_VECS[j][i] * 2
            assert fr.denominator == 1
            s += int(fr) * coef[k]
        if i == 13:
            s += 2 * z
        w2.append(s)
    H = hessian_forms()
    Mx = [[sum((cc * w2[k] for k, cc in H[i][j].items()), 0 * z) for j in range(14)] for i in range(14)]
    D, nd = bareiss(Mx, one)
    check("Bareiss elimination uses 819 exact divisions", nd == 819, f"({time.time() - t0:.1f}s)")
    Dd = D.to_dict()
    zdeg = max(m[6] for m in Dd)
    q, p, t, b, c, v = G[:6]
    d = one
    C = d * (d + v) * (d + q + t) + c * (p * b - t * v)
    Q = p * b + c * (p + b) + d * d - (t + v) * (c - d)
    R = 0 * one
    for k in range(zdeg + 1):
        Dk = cz.from_dict({m[:6] + (0,): cc for m, cc in Dd.items() if m[6] == k})
        R += Dk * (-C) ** k * Q ** (4 - k)
    den = -4 * (d + v) ** 2 * Q ** 2 * 2 ** 14
    qq, rr = divmod(R, den)
    check("Q^4 det H(w(-C/Q)) divisible by -4(d+v)^2 Q^2", rr.is_zero())
    Gp = qq.sqrt()
    check("(6): quotient is a perfect square G^2", Gp is not None and Gp * Gp == qq)
    if Gp.to_dict().get((0,) * 7, 0) < 0:
        Gp = -Gp
    Gd = {tuple(int(x) for x in m[:6]): int(cc) for m, cc in Gp.to_dict().items()}
    Bd = {m[:5]: cc for m, cc in Gd.items() if m[5] == 0}
    Ad = {m[:5]: cc for m, cc in Gd.items() if m[5] == 1}
    check("G: 273 terms, homogeneous degree 8, G(0,..,0,1) = 1, linear in v; B 137 terms, A 136 terms",
          len(Gd) == 273 and Gp.total_degree() == 8 and Gd.get((0,) * 6) == 1
          and max(m[5] for m in Gd) == 1 and len(Bd) == 137 and len(Ad) == 136)
    sw = {(m[0], m[3], m[2], m[1], m[4], m[5]): cc for m, cc in Gd.items()}
    check("(10): G invariant under p <-> b", sw == Gd)
    # (11) homogeneous c = 0 face
    hc = flint.fmpz_mpoly_ctx.get(NAMES + ("d",), "deglex")
    q, p, t, b, c, v, d = hc.gens()
    Gh = hc.from_dict({m + (8 - sum(m),): cc for m, cc in Gd.items() if m[4] == 0})
    check("(11): G|_{c=0} = d(d+v)(d+b)(d+p)[d(q+d)-pb](q+t+d)^2",
          Gh == d * (d + v) * (d + b) * (d + p) * (d * (q + d) - p * b) * (q + t + d) ** 2)
    return Gd, Bd, Ad


def tables(Gd, Bd, Ad):
    def cleared(Rd, nv, ctx, subs):
        md = [max(m[i] for m in Rd) for i in range(nv)]
        one = ctx.from_dict({(0,) * ctx.nvars(): 1})
        pw = []
        for i, (N, D) in enumerate(subs):
            Np, Dp = [one], [one]
            for _ in range(md[i]):
                Np.append(Np[-1] * N)
                Dp.append(Dp[-1] * D)
            pw.append([Np[k] * Dp[md[i] - k] for k in range(md[i] + 1)])
        tot = 0 * one
        for m, cc in Rd.items():
            tt = one * cc
            for i in range(nv):
                tt = tt * pw[i][m[i]]
            tot += tt
        return tot
    total = 0
    cx = flint.fmpz_mpoly_ctx.get(("P", "B", "C", "T", "Q"), "deglex")

    def rep(label, Pl, marker, nterms, mcoef):
        nonlocal total
        cs = Pl.coeffs()
        e = tuple(marker.get(n, 0) for n in cx.names())
        got = int(Pl.to_dict().get(e, 0))
        total += len(cs)
        check(f"table {label}: {len(cs)} terms, all coefficients > 0, marker {marker} coefficient {mcoef}",
              all(x > 0 for x in cs) and len(cs) == nterms and got == mcoef, f"(got {got})")

    P_, B_, C_, T_, Q_ = cx.gens()
    o = cx.from_dict({(0,) * 5: 1})
    mains = {"(15)": [(P_, o + P_), (B_, o + B_), (C_, o + C_)],
             "(16)": [(o + P_, o), (B_, (o + P_) * (o + B_ + C_)), (C_, (o + P_) * (o + B_ + C_))]}
    exp = {("(15)", "ct<1"): ((378, {"C": 2}, 1), (383, {"C": 2}, 1)),
           ("(15)", "ct>1"): ((924, {"C": 3, "T": 3}, 1), (889, {"C": 3, "T": 3}, 1)),
           ("(15)", "ct=1"): ((147, {}, 1), (147, {}, 1)),
           ("(16)", "ct<1"): ((1535, {"C": 2}, 2), (1540, {"C": 2}, 2)),
           ("(16)", "ct>1"): ((5490, {"C": 3, "T": 3}, 2), (5486, {"C": 3, "T": 3}, 2)),
           ("(16)", "ct=1"): ((671, {}, 2), (671, {}, 2))}
    for mn, (ps, bs, cs) in mains.items():
        cN, cD = cs
        chains = {"ct<1": ((Q_, o), (T_ * cD, cN * (o + T_))),
                  "ct>1": (((o + T_) * cD * Q_, cN * T_ * (o + Q_)), ((o + T_) * cD, cN)),
                  "ct=1": ((Q_, o), (cD, cN))}
        for key, (qs, ts) in chains.items():
            subs = [qs, ps, ts, bs, cs]
            (nb, mb, cb), (na, ma, ca) = exp[(mn, key)]
            rep(f"{mn} {key} B", cleared(Bd, 5, cx, subs), mb, nb, cb)
            rep(f"{mn} {key} A", cleared(Ad, 5, cx, subs), ma, na, ca)
    # sector c > 1
    den = (o + C_) * (o + P_ + B_)
    ps, bs, cs = (P_, den), (B_, den), (o + C_, o)
    rows = {"ct<1": ((Q_, o), (T_, (o + C_) * (o + T_))),
            "ct>1": (((C_ * (o + T_) + T_) * Q_, (o + C_) * T_ * (o + Q_)),
                     (C_ * (o + T_) + T_, C_ * (o + T_) * (o + C_))),
            "ct=1": ((Q_, o), (o, o + C_))}
    exp2 = {"ct<1": ((1421, {}, 2), (2400, {}, 3)), "ct>1": ((3000, {"T": 5}, 3), (4320, {"C": 2, "T": 4}, 3)),
            "ct=1": ((480, {}, 10), (600, {}, 6))}
    for key, (qs, ts) in rows.items():
        tN, tD = ts
        VN, VD = tD - C_ * tN, C_ * tD            # V = 1/C - t
        check(f"c>1 {key}: V numerator and denominator have positive coefficients",
              all(x > 0 for x in VN.coeffs()) and all(x > 0 for x in VD.coeffs()))
        subs = [qs, ps, ts, bs, cs]
        (nb, mb, cb), (ng, mg, cg) = exp2[key]
        rep(f"c>1 {key} B", cleared(Bd, 5, cx, subs), mb, nb, cb)
        rep(f"c>1 {key} G(V)", cleared(Gd, 6, cx, subs + [(VN, VD)]), mg, ng, cg)
    # c = 1 wall
    ps, bs, cs = (P_, o + P_), (B_, o + B_), (o, o)
    rows = {"t<1": ((Q_, o), (T_, o + T_)), "t>1": (((o + T_) * Q_, T_ * (o + Q_)), (o + T_, o)),
            "t=1": ((Q_, o), (o, o))}
    exp3 = {"t<1": ((100, {}, 2), (105, {}, 3)), "t>1": ((159, {"T": 3}, 10), (124, {"T": 3}, 6)),
            "t=1": ((35, {}, 10), (35, {}, 6))}
    for key, (qs, ts) in rows.items():
        subs = [qs, ps, ts, bs, cs]
        (nb, mb, cb), (na, ma, ca) = exp3[key]
        rep(f"c=1 {key} B", cleared(Bd, 5, cx, subs), mb, nb, cb)
        rep(f"c=1 {key} A", cleared(Ad, 5, cx, subs), ma, na, ca)
    check("24 tables with 31040 coefficients in total", total == 31040, f"(got {total})")


def numeric():
    import numpy as np
    from scipy.optimize import minimize
    from scipy.linalg import expm
    P = np.array(L.P, float)
    SIG = np.array(L.PERMS)
    SG = np.array(L.SGN, float)
    FV = [np.array([float(x) for x in f]) for f in L.F_VECS]

    def W(O):
        return SG * O[0, SIG[:, 0]] * O[1, SIG[:, 1]] * O[2, SIG[:, 2]] * O[3, SIG[:, 3]]

    def sk(t):
        a, b, c, d, e, f = t
        return np.array([[0, a, b, c], [-a, 0, d, e], [-b, -d, 0, f], [-c, -e, -f, 0]])
    rng = np.random.default_rng(2)

    def minphi(coefs, restarts=12):
        vals = P @ sum(k * FV[j] for j, k in coefs.items())
        best = np.inf
        for _ in range(restarts):
            r = minimize(lambda t: float(vals @ W(expm(sk(t)))), rng.uniform(-np.pi, np.pi, 6), method="BFGS")
            best = min(best, r.fun)
        return best
    w = min(minphi({0: qq, 3: tt, 15: 1 / qq + 1 / tt, 26: 1.0}) for qq in [0.2, 1, 5] for tt in [0.2, 1, 5])
    print(f"[INFO] NUMERIC chain (center f15; leaves f0,f3; f26) at c = 1/q + 1/t, 9 grid points, 12 restarts, "
          f"seed 2: min Phi = {w:.3e}")


if __name__ == "__main__":
    structural()
    Gd, Bd, Ad = cubic_and_G()
    tables(Gd, Bd, Ad)
    if "--numeric" in sys.argv:
        numeric()
    print("SUMMARY:", "all EXACT checks passed" if not FAIL else f"FAILED: {FAIL}")
    sys.exit(1 if FAIL else 0)
