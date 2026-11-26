"""Second-pass independent replay for final_referee_inventory71_uniform_20261023 (paper case T_72).

Own code only.  Nothing from the supplement is executed, imported or read.  Everything is
rebuilt from the paper's definitions through the audit library lib/moc4lib.py.

Usage:  python -I -B <this file> exact   [CACHE_DIR]     (exact checks, ~7-8 minutes)
        python -I -B <this file> numeric [CACHE_DIR]     (floating-point tests; needs the
                                                         G9 cache written by 'exact')
CACHE_DIR defaults to the system temporary directory.  Prints one line per check and exits
nonzero if an EXACT check fails.
"""
import sys, os, time, pickle, tempfile, importlib.util, math

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("moc4lib", os.path.join(HERE, "..", "lib", "moc4lib.py"))
L = importlib.util.module_from_spec(spec); spec.loader.exec_module(L)

MODE = sys.argv[1] if len(sys.argv) > 1 else "exact"
CACHE = sys.argv[2] if len(sys.argv) > 2 else tempfile.gettempdir()
G9_CACHE = os.path.join(CACHE, "inv71_second_G9.pkl")
FAIL = []
T0 = time.time()


def check(name, ok, info=""):
    print(("PASS " if ok else "FAIL ") + name + (" | " + info if info else ""), flush=True)
    if not ok:
        FAIL.append(name)


T = frozenset({1, 3, 4, 7, 8, 12})
TP = frozenset({0, 1, 2, 10, 14, 17})
BASIS = [0, 2, 12, 13, 15, 21, 23, 26]          # q a b v c w t d
READ_AT = [23, 15, 8, 7, 4, 16, 9]
VAL = L.F_VALUES


def exact_part():
    import flint, sympy
    # ---------------- A. contact set and transport ----------------
    check("T equals paper T_72 (Appendix E)", L.CONTACT_SETS[72] == T)
    g = ((0, 2, 1, 3), (0, 2, 3, 1), 0)          # sigma -> rho o sigma o chi
    check("literal action rho=0213, chi=0231 maps T' onto T", L.g_act_set(g, TP) == T)
    idx = {VAL[j]: j for j in range(82)}
    gf = lambda gg, j: idx.get(tuple(L.g_act_values(gg, list(VAL[j]))))
    check("that action permutes the 82 generators", sorted(gf(g, j) for j in range(82)) == list(range(82)))
    # ---------------- generators, relation, representation ----------------
    gens = L.generators(TP)
    check("generators(T') = [0,2,10,12,13,15,21,23,26]", gens == [0, 2, 10, 12, 13, 15, 21, 23, 26], str(gens))
    check("all but f26 are minor vectors; f26 = h4132,- is not",
          [L.F_IS_MINOR[j] for j in gens] == [True] * 8 + [False] and L.f_name(26) == "h4132,-")
    check("relation f10 + f26 = f0 + f13 + f23",
          all(VAL[10][s] + VAL[26][s] == VAL[0][s] + VAL[13][s] + VAL[23][s] for s in range(24)))
    check("contact rows of P have rank 6", L.rank([L.P[s] for s in TP]) == 6)
    check("the 8 basis profiles are independent", L.rank([list(L.F_VECS[j]) for j in BASIS]) == 8)
    ok = all(VAL[k][p] == (1 if k == j else 0) for j, p in zip(BASIS[:7], READ_AT) for k in BASIS)
    check("coefficients q,a,b,v,c,w,t are the gaps at 23,15,8,7,4,16,9", ok)
    lin = lambda p: {k: VAL[k][p] for k in BASIS if VAL[k][p]}
    check("(2): g20 = q+d, g12 = v+d, g11 = t+d, g9 = t",
          lin(20) == {0: 1, 26: 1} and lin(12) == {13: 1, 26: 1} and lin(11) == {23: 1, 26: 1} and lin(9) == {23: 1})
    # ---------------- reference-atom images ----------------
    def has_image(src, tgt):
        return any(tuple(gf(gg, j) for j in src) == tgt for gg in L.G_ELEMENTS)
    check("SOC (f1,f9;f22) has images (f2,f12;f26), (f2,f15;f26), (f12,f15;f26)",
          all(has_image((1, 9, 22), t) for t in [(2, 12, 26), (2, 15, 26), (12, 15, 26)]))
    check("old-H (f3,f15;f26) -> (f23,f15;f26) with f3->f23, f15->f15", has_image((3, 15, 26), (23, 15, 26)))
    check("chain (f2; f0,f12; f70) -> (f15; f21,f23; f26) role by role", has_image((2, 0, 12, 70), (15, 21, 23, 26)))
    # ---------------- (8): cubic on the lift ----------------
    names = ('q', 'a', 'b', 'v', 'c', 'w', 't', 'd', 'z')
    ctx = flint.fmpz_mpoly_ctx.get(names, 'degrevlex')
    q, a, b, v, c, w, t, d, z = ctx.gens()
    zero = ctx.from_dict({})
    coef = [q, a, b, v, c, w, t, d]
    m = []
    for i in range(14):
        s = zero
        for k, j in enumerate(BASIS):
            f = L.F_VECS[j][i]
            assert f.denominator == 1
            s = s + coef[k] * int(f.numerator)
        m.append(s)
    mz = m[:13] + [m[13] + z]
    F = L.cubic_C(mz)
    Cc = (q + d) * (v + d) * (w + t + d) + a * b * c
    Qc = a * b + a * c + b * c + d * d - t * (c - d)
    check("(8): C(g + z kappa) = C + z Q with the stated C, Q (paper cubic (A.2))", F == Cc + z * Qc)
    # ---------------- (9): Hessian determinant ----------------
    ms, Cs = L.cubic_sympy()
    def tolin(e):
        s = zero
        if e == 0:
            return s
        for mon, cf in sympy.Poly(e, *ms).terms():
            k = [i for i, e2 in enumerate(mon) if e2]
            assert len(k) == 1 and mon[k[0]] == 1
            s = s + mz[k[0]] * int(cf)
        return s
    A = [[tolin(sympy.diff(Cs, ms[i], ms[j])) for j in range(14)] for i in range(14)]
    zpos = sorted((i, j) for i in range(14) for j in range(14) if A[i][j].degrees()[8] > 0)
    check("z enters the Hessian only at (0,6),(1,5),(5,1),(6,0) (H_kappa of (3.10))",
          zpos == [(0, 6), (1, 5), (5, 1), (6, 0)])
    one = ctx.from_dict({(0,) * 9: 1})
    M = [row[:] for row in A]; n = 14; prev = one; sign = 1
    for k in range(n - 1):
        if M[k][k] == 0:
            r = next(r for r in range(k + 1, n) if M[r][k] != 0)
            M[k], M[r] = M[r], M[k]; sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                num = M[i][j] * M[k][k] - M[i][k] * M[k][j]
                M[i][j] = num / prev if num != 0 else zero     # exact division (raises if inexact)
        prev = M[k][k]
    Dz = sign * M[n - 1][n - 1]
    check("det Hess C(g + z kappa) computed by fraction-free elimination", True,
          f"{len(Dz)} terms, degree {Dz.total_degree()}, z-degree {Dz.degrees()[8]}, {time.time()-T0:.0f}s")
    parts = {}
    for mon, cf in Dz.to_dict().items():
        parts.setdefault(mon[8], {})[tuple(mon[:8]) + (0,)] = cf
    R = zero
    for k, p in parts.items():
        R = R + ctx.from_dict(p) * (-Cc) ** k * Qc ** (4 - k)
    S = R / (-4 * Qc * Qc)
    G9 = S.sqrt()
    if G9.coefficient(G9.monoms().index((0,) * 7 + (9, 0))) < 0:
        G9 = -G9
    check("(9): Q^4 det H(g - C kappa/Q) = -4 Q^2 G9^2 for a polynomial G9", R == -4 * Qc * Qc * G9 * G9)
    check("G9: 1615 terms, homogeneous of degree 9, [d^9] = 1, degree profile (3,2,2,3,2,3,2,9)",
          len(G9) == 1615 and all(sum(mn) == 9 for mn in G9.monoms())
          and G9.coefficient(G9.monoms().index((0,) * 7 + (9, 0))) == 1
          and tuple(G9.degrees()[:8]) == (3, 2, 2, 3, 2, 3, 2, 9))
    check("(10): G9(q,a,b,v,..) = G9(v,b,a,q,..); C and Q also symmetric",
          G9.compose(v, b, a, q, c, w, t, d, z) == G9 and Cc.compose(v, b, a, q, c, w, t, d, z) == Cc
          and Qc.compose(v, b, a, q, c, w, t, d, z) == Qc)
    Gd = {tuple(int(e) for e in mn[:8]): int(cf) for mn, cf in zip(G9.monoms(), G9.coeffs())}
    pickle.dump(Gd, open(G9_CACHE, "wb"))
    # ---------------- (17)-(19): positivity tables ----------------
    MUL = (3, 2, 2, 3, 2, 3, 2)
    hn = []
    for nm in ('q', 'a', 'b', 'v', 'c', 'w', 't'):
        hn += [nm + 'N', nm + 'D']
    hctx = flint.fmpz_mpoly_ctx.get(tuple(hn), 'degrevlex')
    hd = {}
    for mn, cf in Gd.items():
        e = []
        for i in range(7):
            e += [mn[i], MUL[i] - mn[i]]
        assert min(e) >= 0
        hd[tuple(e)] = cf
    GH = hctx.from_dict(hd)
    cn = ('X', 'Y', 'Z', 'T', 'W', 'Q0', 'V')
    cctx = flint.fmpz_mpoly_ctx.get(cn, 'degrevlex')
    X, Y, Z, TT, W, Q0, V = cctx.gens()
    o = cctx.from_dict({(0,) * 7: 1}); zc = cctx.from_dict({}); Lx = 1 + X
    mains = {
        '(11)': dict(a=(Y + Z, 1 + X + Y + Z), b=(X + Z, 1 + X + Y + Z), c=(X + Y, 1 + X + Y + Z)),
        '(12)': dict(a=(1 + Y * Lx * Lx, Lx * (1 + Y)), b=(Z, Lx * (1 + Z)), c=(o, Lx * (1 + Z))),
        '(13)': dict(a=(Z, Lx * (1 + Y) * (1 + Z)), b=(o, Lx * (1 + Y) * (1 + Z)), c=(o, Lx)),
        '(14)': dict(a=(Y * Z, Lx * (1 + Y) * (1 + Z)), b=(Y, Lx * (1 + Y) * (1 + Z)), c=(Lx, o)),
        'solo': dict(a=(zc, o), b=(zc, o), c=(o, 1 + X)),
    }
    expected = {('(11)', 'low'): 28693, ('(11)', 'eq'): 10012, ('(11)', 'high'): 93570,
                ('(12)', 'low'): 27186, ('(12)', 'eq'): 10665, ('(12)', 'high'): 70122,
                ('(13)', 'low'): 36773, ('(13)', 'eq'): 13175, ('(13)', 'high'): 70650,
                ('(14)', 'low'): 38081, ('(14)', 'eq'): 12700, ('(14)', 'high'): 79550,
                ('solo', 'low'): 703, ('solo', 'eq'): 271, ('solo', 'high'): 1674}
    # strict markers: list of alternative monomials, at least one of which is positive on the domain
    e7 = lambda **kw: tuple(kw.get(nm, 0) for nm in cn)
    markers = {('(11)', 'low'): [e7(X=2), e7(Y=2)], ('(11)', 'high'): [e7(T=3, X=3), e7(T=3, Y=3)],
               ('(12)', 'low'): [e7(X=1)], ('(12)', 'high'): [e7(T=3)],
               ('(13)', 'low'): [e7(X=1), e7(Y=1)], ('(13)', 'high'): [e7(T=3)],
               ('(14)', 'high'): [e7(T=5)], ('solo', 'high'): [e7(T=3)]}
    total = 0
    for key, mm in mains.items():
        cN, cD = mm['c']
        for row in ('low', 'eq', 'high'):
            sub = dict(q=(Q0, o), v=(V, o), **mm)
            if row == 'low':
                sub.update(t=(TT * cD, cN * (1 + TT)), w=(W, o))
            elif row == 'eq':
                sub.update(t=(cD, cN), w=(W, o))
            elif key == '(14)':                       # (16), X = c - 1 > 0
                num = X * (1 + TT) + TT
                sub.update(t=(num, X * Lx * (1 + TT)), w=(num * W, Lx * TT * (1 + W)))
            else:
                sub.update(t=((1 + TT) * cD, cN), w=((1 + TT) * cD * W, cN * TT * (1 + W)))
            args = []
            for nm in ('q', 'a', 'b', 'v', 'c', 'w', 't'):
                args += list(sub[nm])
            Pp = GH.compose(*args, ctx=cctx)
            total += len(Pp)
            mons = Pp.monoms(); cfs = [int(x) for x in Pp.coeffs()]
            allpos = all(x > 0 for x in cfs)
            mk = markers.get((key, row), [(0,) * 7])
            mset = set(mons)
            check(f"table {key} {row}: {len(Pp)} terms (paper {expected[(key, row)]}), all coefficients > 0, strict marker present",
                  len(Pp) == expected[(key, row)] and allpos and all(x in mset for x in mk))
    for name, ab, exp in (('(19a)', dict(a=(X, 1 + X), b=(Y, 1 + Y)), 1293),
                          ('(19b)', dict(a=(1 + X, o), b=(Y, (1 + X) * (1 + Y))), 1680)):
        sub = dict(q=(Q0, o), v=(V, o), c=(zc, o), w=(W, o), t=(TT, o), **ab)
        args = []
        for nm in ('q', 'a', 'b', 'v', 'c', 'w', 't'):
            args += list(sub[nm])
        Pp = GH.compose(*args, ctx=cctx)
        total += len(Pp)
        cfs = [int(x) for x in Pp.coeffs()]
        check(f"table {name}: {len(Pp)} terms (paper {exp}), all coefficients > 0, constant term > 0",
              len(Pp) == exp and all(x > 0 for x in cfs) and (0,) * 7 in set(Pp.monoms()))
    check("total coefficients = 496798", total == 496798, str(total))
    # ---------------- chart inverses (coverage) ----------------
    Xs, Ys, Zs, Ts, Ws = sympy.symbols('X Y Z T W', positive=True)
    aS, bS, cS, tS, wS = sympy.symbols('a b c t w', positive=True)
    simp = lambda e: sympy.simplify(sympy.together(e)) == 0
    # (11): forward then inverse
    D = 1 + Xs + Ys + Zs
    fa, fb, fc = (Ys + Zs) / D, (Xs + Zs) / D, (Xs + Ys) / D
    s = fa + fb + fc; den = 1 - s / 2
    ok11 = all(simp(e) for e in [(fb + fc - fa) / 2 / den - Xs, (fa + fc - fb) / 2 / den - Ys, (fa + fb - fc) / 2 / den - Zs])
    # (12)
    Lc = 1 + Xs
    fa, fb, fc = (1 + Ys * Lc ** 2) / (Lc * (1 + Ys)), Zs / (Lc * (1 + Zs)), 1 / (Lc * (1 + Zs))
    S = fb + fc
    ok12 = all(simp(e) for e in [1 / S - 1 - Xs, (fa - S) / (1 / S - fa) - Ys, fb / fc - Zs])
    # (13)
    fc, fa, fb = 1 / Lc, Zs / (Lc * (1 + Ys) * (1 + Zs)), 1 / (Lc * (1 + Ys) * (1 + Zs))
    S = fa + fb
    ok13 = all(simp(e) for e in [1 / fc - 1 - Xs, fc / S - 1 - Ys, fa / fb - Zs])
    # (14)
    fc, fa, fb = Lc, Ys * Zs / (Lc * (1 + Ys) * (1 + Zs)), Ys / (Lc * (1 + Ys) * (1 + Zs))
    S = fa + fb
    ok14 = all(simp(e) for e in [fc * S / (1 - fc * S) - Ys, fa / fb - Zs, fc * S - Ys / (1 + Ys)])
    # (15) high and (16)
    ft = (1 + Ts) / cS; fw = ft * Ws / (Ts * (1 + Ws))
    ok15 = simp(cS * ft - 1 - Ts) and simp(fw / (ft / (cS * ft - 1) - fw) - Ws)
    Xc = cS - 1
    tau = Ts / (Xc * (1 + Ts)); ft = (1 + tau) / cS; fw = ft * Ws / (tau * (1 + Ws))
    ok16 = (simp(Xc * (cS * ft - 1) / (1 - Xc * (cS * ft - 1)) - Ts) and simp(fw / (ft / (cS * ft - 1) - fw) - Ws)
            and simp(ft - (Xc * (1 + Ts) + Ts) / (Xc * cS * (1 + Ts))) and simp(1 / Xc - ft - (1 - Xc * tau) / (Xc * cS)))
    # (19b)
    fa, fb = 1 + Xs, Ys / ((1 + Xs) * (1 + Ys))
    ok19 = simp(fa * fb / (1 - fa * fb) - Ys)
    check("chart inverses (11),(12),(13),(14),(15),(16),(19) recover the chart parameters",
          ok11 and ok12 and ok13 and ok14 and ok15 and ok16 and ok19)


def numeric_part():
    import numpy as np
    from scipy.optimize import minimize
    Gd = pickle.load(open(G9_CACHE, "rb"))
    Pn = L.P_numpy()
    val = {j: Pn @ np.array([float(x) for x in L.F_VECS[j]]) for j in range(82)}
    phi = lambda j, O: float(val[j] @ L.leibniz_weights_np(O))
    def minor(j, O):
        k, I, J = L.F_TABLE[j]
        return float(np.linalg.det(O[np.ix_([int(x) - 1 for x in I], [int(x) - 1 for x in J])]))
    forms = {}
    for (i, j) in ((2, 12), (2, 15), (12, 15)):
        forms[f"SOC image: Phi_f26 + 2|m{i} m{j}|"] = (lambda O, i=i, j=j: phi(26, O) + 2 * abs(minor(i, O) * minor(j, O)))
    forms["old-H image: Phi_f26 + m15^2 + 2|m23 m15|"] = lambda O: phi(26, O) + minor(15, O) ** 2 + 2 * abs(minor(23, O) * minor(15, O))
    forms["chain image: Phi_f26 + 2|m15|(|m21|+|m23|)"] = lambda O: phi(26, O) + 2 * abs(minor(15, O)) * (abs(minor(21, O)) + abs(minor(23, O)))
    for name, f in forms.items():
        best = np.inf; rng = np.random.default_rng(1)
        for _ in range(100):
            fun = lambda th: f(L.so4_from_params(th))
            r = minimize(fun, rng.uniform(-np.pi, np.pi, 6), method="BFGS")
            r2 = minimize(fun, r.x, method="Nelder-Mead", options=dict(xatol=1e-12, fatol=1e-14, maxiter=4000))
            best = min(best, r.fun, r2.fun)
        print(f"NUMERIC {name} >= 0 on SO(4): min ~ {best:.2e} (100 BFGS+NM restarts, seed 1)", flush=True)
    best = np.inf; rng = np.random.default_rng(2)
    for _ in range(30):
        r = minimize(lambda th: phi(26, L.so4_from_params(th)), rng.uniform(-np.pi, np.pi, 6), method="BFGS"); best = min(best, r.fun)
    print(f"NUMERIC Phi_f26 alone takes negative values: min ~ {best:.6f}", flush=True)
    def G9eval(x):
        s = sa = 0.0
        for mn, cf in Gd.items():
            tt = cf * np.prod([x[i] ** int(mn[i]) for i in range(8)]); s += tt; sa += abs(tt)
        return s, sa
    worst = np.inf; nd = 0; nsafe = 0; worstG = 0.0; nf = 0
    for seed in range(30):
        for margin in (0.0, 1e-2):
            r = L.random_point_on_contact_stratum(TP, seed=seed, tries=20, margin=margin)
            if r is None:
                continue
            nf += 1
            gv = r['gaps']
            x = np.array([gv[k] for k in READ_AT] + [gv[11] - gv[9]])
            qq, aa, bb, vv, cc, ww, tt, dd = x
            Q = aa * bb + aa * cc + bb * cc + dd * dd - tt * (cc - dd)
            if margin > 0 and abs(Q) > 1e-6:
                Gv, Ga = G9eval(x); worstG = max(worstG, abs(Gv) / Ga)
            if dd > 1e-6:
                nd += 1
                A_, B_, C_ = [max(u, 0) / dd for u in (aa, bb, cc)]
                s = A_ + B_ + C_; mx = max(A_, B_, C_)
                H3 = s / 2 if 2 * mx <= s else math.sqrt(mx * (s - mx))
                t_, w_ = max(tt, 0) / dd, max(ww, 0) / dd
                if H3 >= 1 - 1e-9 or t_ * (C_ - 1) >= 1 - 1e-9 or (w_ > 0 and t_ > 0 and C_ >= 1 / w_ + 1 / t_ - 1e-9):
                    nsafe += 1
            cvec = np.linalg.lstsq(Pn, gv, rcond=None)[0]
            val_, _ = L.min_phi_over_SO4(cvec, restarts=10, seed=seed)
            worst = min(worst, val_)
    print(f"NUMERIC stratum points with contacts containing T' (seeds 0-29, margins 0 and 1e-2): {nf} found; "
          f"worst min Phi over SO(4) = {worst:.2e}; points with d > 0: {nd}, of which satisfying (3), (4) or (5): {nsafe}; "
          f"max |G9|/sum|terms| at exact-contact points = {worstG:.1e}", flush=True)


if MODE == "exact":
    exact_part()
else:
    numeric_part()
print(f"done ({MODE}) in {time.time()-T0:.0f}s; EXACT failures: {len(FAIL)}")
sys.exit(1 if FAIL else 0)
