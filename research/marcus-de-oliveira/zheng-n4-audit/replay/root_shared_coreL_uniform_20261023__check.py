"""Independent replay for root_shared_coreL_uniform_20261023 (seven-contact set L).

Own code only. Nothing from the supplement is executed, imported or read.
Builds everything from the shared audit library (../lib/moc4lib.py: the paper's P, f_j,
cubic C of (A.2), group G) and python-flint.

Usage:  python -I -B root_shared_coreL_uniform_20261023__check.py [--numeric]
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


T_L = [0, 1, 2, 7, 10, 16, 23]
NAMES = ("a", "b", "q", "c", "e", "v")       # d = 1 after scaling; d kept where needed
FIDX = [2, 3, 11, 12, 15, 23, 26]            # coefficients a,b,q,c,e,v,d


# ---------------------------------------------------------------- A/B: contact set, basis, gaps
def structural():
    check("L is Table 5 set #14", L.AUX7[14] == frozenset(T_L))
    gens = L.generators(T_L)
    check("generators(L) = {2,3,8,11,12,15,23,26}, only f26 non-minor",
          gens == [2, 3, 8, 11, 12, 15, 23, 26] and [j for j in gens if not L.F_IS_MINOR[j]] == [26])
    check("contact rows of P on L have rank 7", L.rank([L.P[s] for s in T_L]) == 7)
    check("the seven listed profiles are independent", L.rank([L.F_VECS[j] for j in FIDX]) == 7)
    V = L.F_VALUES
    check("f8 + f26 = f2 + f3 + f12 + f23", all(V[8][s] + V[26][s] == V[2][s] + V[3][s] + V[12][s] + V[23][s]
                                                for s in range(24)))
    # gaps of the cost a f2 + b f3 + q f11 + c f12 + e f15 + v f23 + d f26
    ctx = flint.fmpz_mpoly_ctx.get(NAMES + ("d",), "deglex")
    a, b, q, c, e, v, d = ctx.gens()
    co = [a, b, q, c, e, v, d]
    g = [sum((co[k] * V[j][s] for k, j in enumerate(FIDX)), 0 * a) for s in range(24)]
    check("cost vanishes on L", all(g[s] == 0 for s in T_L))
    check("coefficients = gaps (g15,g14,g17,g8,g4,g9) = (a,b,q,c,e,v)",
          [g[15], g[14], g[17], g[8], g[4], g[9]] == [a, b, q, c, e, v])
    check("g13=a+d, g20=b+d, g22=c+d, g3=v+d", [g[13], g[20], g[22], g[3]] == [a + d, b + d, c + d, v + d])
    # d <= 0 decomposition: cost = (a+d)f2+(b+d)f3+q f11+(c+d)f12+e f15+(v+d)f23+(-d)f8
    alt = [sum(((a + d) * V[2][s], (b + d) * V[3][s], q * V[11][s], (c + d) * V[12][s], e * V[15][s],
                (v + d) * V[23][s], -d * V[8][s]), 0 * a) for s in range(24)]
    check("d<=0 rewrite as minor combination (f8 minor)", alt == g and L.F_IS_MINOR[8])
    # group images of the cited atoms
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
    check("some g in G maps old-H (f3,f15;f26) to (f23,f15;f26)", exists((3, 15, 26), (23, 15, 26)))


# ---------------------------------------------------------------- cubic, C, Q, identity (9)
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


def cubic_and_G8():
    # homogeneous check of F(w0 + z kappa) = C + z Q
    ctx = flint.fmpq_mpoly_ctx.get(NAMES + ("d", "z"), "deglex")
    a, b, q, c, e, v, d, z = ctx.gens()
    co = [a, b, q, c, e, v, d]
    w = [sum((co[k] * flint.fmpq(L.F_VECS[j][i].numerator, L.F_VECS[j][i].denominator)
              for k, j in enumerate(FIDX)), 0 * a) + (z if i == 13 else 0) for i in range(14)]
    F = L.cubic_C(w)
    C = (b + d) * (q + d) * (v + d) + a * c * e - b * e * v
    Q = a * c + a * e + b * q - b * e + b * d - q * c + q * v + q * d + c * e - e * v + v * d + d * d
    check("(6)-(7): paper cubic C(w0 + z kappa) = C + z Q", F == C + z * Q)
    # d = 1 identities
    one = ctx.from_dict({(0,) * 8: 1})
    C1 = (b + 1) * (q + 1) * (v + 1) + a * c * e - b * e * v
    xi, Pp = b + v, b * v
    check("C = ace + (q+1)(xi+1) + P(q+1-e) at d=1", C1 == a * c * e + (q + 1) * (xi + 1) + Pp * (q + 1 - e))
    lhs = 4 * xi * (C1 - 1 - flint.fmpq(3, 4) * xi - q * (xi + 1 + Pp) - a * c * e)
    check("(9): 4 xi [C - 1 - 3xi/4 - q(xi+1+P) - ace] = 4P(1 - xi(e-1)) + (b-v)^2",
          lhs == 4 * Pp * (1 - xi * (e - 1)) + (b - v) ** 2)
    del one
    # Hessian determinant along the kappa line, d = 1, integer scaling w2 = 2w
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
    a, b, q, c, e, v = G[:6]
    d = one
    C = (b + d) * (q + d) * (v + d) + a * c * e - b * e * v
    Q = a * c + a * e + b * q - b * e + b * d - q * c + q * v + q * d + c * e - e * v + v * d + d * d
    R = 0 * one
    for k in range(zdeg + 1):
        Dk = cz.from_dict({m[:6] + (0,): cc for m, cc in Dd.items() if m[6] == k})
        R += Dk * (-C) ** k * Q ** (4 - k)
    den = -4 * (a + 1) ** 2 * Q ** 2 * 2 ** 14     # det H(2w) = 2^14 det H(w)
    qq, rr = divmod(R, den)
    check("Q^4 det H(w0 - C kappa/Q) divisible by -4(a+d)^2 Q^2", rr.is_zero())
    G8 = qq.sqrt()
    okc = G8 is not None and G8 * G8 == qq
    check("(10): quotient is a perfect square G8^2", okc)
    if G8.to_dict().get((0,) * 7, 0) < 0:
        G8 = -G8
    Gd = {tuple(int(x) for x in m[:6]): int(cc) for m, cc in G8.to_dict().items()}
    md = [max(m[i] for m in Gd) for i in range(6)]
    check("G8: 327 terms, [d^8]G8 = 1, degrees (a,b,q,c,e,v) = (1,2,2,2,2,2), total 8",
          len(Gd) == 327 and Gd.get((0,) * 6) == 1 and md == [1, 2, 2, 2, 2, 2] and G8.total_degree() == 8)
    sw = {(m[0], m[5], m[2], m[3], m[4], m[1]): cc for m, cc in Gd.items()}
    check("G8 invariant under b <-> v", sw == Gd)
    return Gd


# ---------------------------------------------------------------- the 12 coefficient tables
def tables(Gd):
    md = [max(m[i] for m in Gd) for i in range(6)]

    def cleared(ctx, subs):
        one = ctx.from_dict({(0,) * ctx.nvars(): 1})
        pw = []
        for i, (N, D) in enumerate(subs):
            Np, Dp = [one], [one]
            for _ in range(md[i]):
                Np.append(Np[-1] * N)
                Dp.append(Dp[-1] * D)
            pw.append([Np[k] * Dp[md[i] - k] for k in range(md[i] + 1)])
        tot = 0 * one
        for m, cc in Gd.items():
            t = one * cc
            for i in range(6):
                t = t * pw[i][m[i]]
            tot += t
        return tot

    total = 0

    def rep(label, P, ctx, markers, nterms):
        nonlocal total
        cs = P.coeffs()
        dct = P.to_dict()
        mk = [int(dct.get(tuple(mm.get(n, 0) for n in ctx.names()), 0)) for mm, _ in markers]
        ok = all(x >= 0 for x in cs) and len(cs) == nterms and mk == [w for _, w in markers]
        total += len(cs)
        check(f"table {label}: {len(cs)} terms, all coefficients >= 0, markers {markers}", ok,
              f"(got markers {mk})")

    cx = flint.fmpz_mpoly_ctx.get(("X", "Y", "Z", "b", "q", "v"), "deglex")
    X, Y, Z, b, q, v = cx.gens()
    o = cx.from_dict({(0,) * 6: 1})
    D = o + X + Y + Z
    rep("balanced (12)", cleared(cx, [(Y + Z, D), (b, o), (q, o), (X + Z, D), (X + Y, D), (v, o)]), cx,
        [({}, 1)], 1493)
    mn, mdn = o + Y * (o + X) ** 2, (o + X) * (o + Y)
    s1, s2 = (Z, (o + X) * (o + Z)), (o, (o + X) * (o + Z))
    for main, nt, mks in (("a", 1411, (4, 4)), ("c", 2113, (8, 4))):
        for k, (u, w) in enumerate([(s1, s2), (s2, s1)]):
            subs = ([(mn, mdn), (b, o), (q, o), u, w, (v, o)] if main == "a"
                    else [u, (b, o), (q, o), (mn, mdn), w, (v, o)])
            rep(f"dominant {main} (13) order {k}", cleared(cx, subs), cx, [({"X": 1}, mks[k])], nt)
    c4 = flint.fmpz_mpoly_ctx.get(("m", "b", "q", "v"), "deglex")
    m_, b4, q4, v4 = c4.gens()
    o4 = c4.from_dict({(0,) * 4: 1})
    zr = 0 * o4
    rep("solo a (c=e=0)", cleared(c4, [(m_, o4), (b4, o4), (q4, o4), (zr, o4), (zr, o4), (v4, o4)]), c4,
        [({}, 1)], 27)
    rep("solo c (a=e=0)", cleared(c4, [(zr, o4), (b4, o4), (q4, o4), (m_, o4), (zr, o4), (v4, o4)]), c4,
        [({}, 1)], 72)
    ce = flint.fmpz_mpoly_ctx.get(("E", "T", "Z", "b", "q", "v"), "deglex")
    E, T, Z, b, q, v = ce.gens()
    o = ce.from_dict({(0,) * 6: 1})
    den = (o + E) * (o + T) * (o + Z)
    s1, s2 = (Z, den), (o, den)
    for k, (u, w) in enumerate([(s1, s2), (s2, s1)]):
        rep(f"dominant e<=1 (14) order {k}", cleared(ce, [u, (b, o), (q, o), w, (o, o + E), (v, o)]), ce,
            [({"E": 1}, (8, 4)[k]), ({"T": 1}, (4, 2)[k])], 2572)
    rep("e<=1 with a=c=0", cleared(ce, [(0 * o, o), (b, o), (q, o), (0 * o, o), (o, o + E), (v, o)]), ce,
        [({}, 2)], 80)
    cf = flint.fmpz_mpoly_ctx.get(("E", "T", "Z", "B", "V", "q"), "deglex")
    E, T, Z, B, V, q = cf.gens()
    o = cf.from_dict({(0,) * 6: 1})
    Sd = (o + E) * (o + T) * (o + Z)
    u1, u2 = (T * Z, Sd), (T, Sd)
    bb, vv = (B * V, E * (o + B) * (o + V)), (B, E * (o + B) * (o + V))
    for k, (u, w) in enumerate([(u1, u2), (u2, u1)]):
        rep(f"dominant e>1 (15)-(16) order {k}", cleared(cf, [u, bb, (q, o), w, (o + E, o), vv]), cf,
            [({"E": 4}, 2)], 7860)
    check("total number of coefficients = 29584", total == 29584, f"(got {total})")


# ---------------------------------------------------------------- NUMERIC (optional)
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

    rng = np.random.default_rng(1)

    def minphi(coefs, restarts=12):
        vals = P @ sum(k * FV[j] for j, k in coefs.items())
        best = np.inf
        for _ in range(restarts):
            r = minimize(lambda t: float(vals @ W(expm(sk(t)))), rng.uniform(-np.pi, np.pi, 6), method="BFGS")
            best = min(best, r.fun)
        return best
    grid = [0.05, 0.2, 0.5, 1, 2, 5, 20]
    for A, B in [(2, 12), (2, 15), (12, 15)]:
        w = min(minphi({A: s, B: 1 / s, 26: 1.0}) for s in grid)
        print(f"[INFO] NUMERIC SOC atom (f{A},f{B};f26) at xy = delta^2, s-grid {grid}, 12 restarts, seed 1: "
              f"min Phi = {w:.3e}")
    for aux in (3, 23):
        w = min(minphi({aux: s, 15: 1 + 1 / s, 26: 1.0}) for s in grid)
        print(f"[INFO] NUMERIC old-H atom (f{aux},f15;f26) at b(e-delta) = delta^2: min Phi = {w:.3e}")


if __name__ == "__main__":
    structural()
    Gd = cubic_and_G8()
    tables(Gd)
    if "--numeric" in sys.argv:
        numeric()
    print("SUMMARY:", "all EXACT checks passed" if not FAIL else f"FAILED: {FAIL}")
    sys.exit(1 if FAIL else 0)
