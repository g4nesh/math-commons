"""Own-code replay for the Zheng n=4 audit, document blind_inventory32_native_whole_20261023.
Exact checks run by default; floating-point falsification runs with --numeric.
Supplement data (JSON only, loaded as data) are located through $MOC4_SUPP or --supp PATH.
Nothing from the supplement is imported or executed.  Exit code 1 if an EXACT check fails."""
import sys, os, json, itertools, random, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import moc4lib as L
import sympy as sp
args = sys.argv[1:]
SUPP = os.environ.get("MOC4_SUPP")
if "--supp" in args: SUPP = args[args.index("--supp") + 1]
NUMERIC = "--numeric" in args
PARTS = []
KJ = None

def part_contacts():

    STAGE6 = {1: [1,2,3,4,5,6,7,15,18,19,20,21,22,23,25,27,28,29,31,34,35,36,37,39,40,42,43,44,72,88],
              2: [10,26,46,90], 3: [9,11,12,47], 4: [8,13,17,49], 5: [14,16,24,33,38,70],
              6: [52,57,58,65,68,71,76,78], 7: [45,51,66,69,75,81], 8: [48,50,53,59,62,67,73,79,86],
              9: [30,41,54,55,60,61,63,74,77,80,83,84], 10: [56,64,82,85,87], 11: [32], 12: [89]}
    STAGE7 = {1: [1,2,5,7,8,9,10,11,12,13,14,15,16,17,21,22,24], 2: [3,4,23,26], 4: [6,18], 5: [20,27], 7: [19,25]}
    st6 = {k: s for s, ks in STAGE6.items() for k in ks}
    st7 = {k: s for s, ks in STAGE7.items() for k in ks}
    assert len(st6) == 90
    rep6 = {L.orbit_rep(T): k for k, T in L.CONTACT_SETS.items()}
    rep7 = {L.orbit_rep(T): k for k, T in L.AUX7.items()}

    def trivial(S):
        return all(L.F_IS_MINOR[j] for j in L.generators(S))

    def classify7(T, j):
        """all ways the 7-set T+{j} is covered"""
        U = set(T) | {j}
        out = []
        if all(L.F_IS_MINOR[k] for k in L.generators(U)):
            out.append(("all-generators-minor-7", None, None))
        r = L.orbit_rep(U)
        if r in rep7:
            out.append(("AUX7", rep7[r], st7[rep7[r]]))
        for S in itertools.combinations(sorted(U), 6):
            if set(S) == set(T):
                continue
            if trivial(S):
                out.append(("trivial6", S, None))
            else:
                k = rep6.get(L.orbit_rep(S))
                out.append(("T", k, st6.get(k), S))
        return out

    def report(name, T, k):
        T = frozenset(T)
        print("==", name, sorted(T))
        print("paper T_%d =" % k, sorted(L.CONTACT_SETS[k]), "find_g:", L.find_g(L.CONTACT_SETS[k], T))
        gens = L.generators(T)
        print("generators:", gens, "nonminor:", [j for j in gens if not L.F_IS_MINOR[j]],
              [L.f_name(j) for j in gens])
        for j in gens:
            print("  f%d %s values" % (j, L.f_name(j)), L.F_VALUES[j])
        for j in range(24):
            if j in T: continue
            nm = [L.F_VALUES[g][j] for g in gens if not L.F_IS_MINOR[g]]
            print(" new index", j, "nonminor values", nm, "->", classify7(T, j) if not all(v > 0 for v in nm) else "pure minor")


    res = [True]
    report("doc32", {0,3,6,7,11,12}, 33)
    return True

print('==== part contacts', flush=True)
PARTS.append(('contacts', part_contacts()))

def part_cone():
    import sympy as sp

    def check(name, T, basis, names, gapclaims, idents):
        T = sorted(T)
        rows = [L.P[t] for t in T]
        print(name, "rank of contact rows", L.rank(rows))
        syms = sp.symbols(names)
        vec = [sum(s * L.F_VALUES[j][k] for s, j in zip(syms, basis)) for k in range(24)]
        # basis vanishes on T and independent
        print(" basis vanish on T:", all(L.F_VALUES[j][t] == 0 for j in basis for t in T),
              " rank of basis profiles:", L.rank([list(L.F_VALUES[j]) for j in basis]))
        ok = True
        for k, expr in gapclaims.items():
            e = sp.sympify(expr, locals=dict(zip(names, syms)))
            good = sp.expand(vec[k] - e) == 0
            ok &= good
            print("  gap", k, "claimed", expr, "computed", vec[k], "OK" if good else "FAIL")
        for k in range(24):
            if k not in gapclaims and k not in T:
                print("  gap", k, "(not stated)", vec[k])
        for lhs, rhs in idents:
            a = [sum(L.F_VALUES[j][k] for j in lhs) for k in range(24)]
            b = [sum(L.F_VALUES[j][k] for j in rhs) for k in range(24)]
            print("  identity", lhs, "=", rhs, "OK" if a == b else "FAIL")
            ok &= (a == b)
        # nonminor values at stated indices
        return ok

    ok2 = check("doc32", [0,3,6,7,11,12], [0,1,2,7,8,15,29,39], "x0 x1 a b c f l g".split(),
                {13:"a",9:"c",2:"l",5:"f",1:"g",14:"b+c",10:"b+l",20:"x0+b",16:"x1+b",17:"x1+g",18:"x0+f",21:"x0+a",23:"x0+x1"},
                [((18,0),(1,2,7,15)), ((40,1),(0,7,39)), ((30,7),(0,1,8,29))])
    vv = [(L.F_VALUES[j][8], L.F_VALUES[j][19]) for j in (18, 30, 40)]
    print(" f18,f30,f40 at 8 and 19:", vv)
    res = [ok2, vv == [(1, 1)] * 3]
    return all(res)

print('==== part cone', flush=True)
PARTS.append(('cone', part_cone()))

def part_chart():
    import sympy as sp
    from sympy.polys.fields import field

    res = []
    def kchk(name, ok):
        res.append(bool(ok)); print(("OK  " if ok else "FAIL"), name)

    # ---- forward direction in original spectral variables
    K1, Zo, bo, co, do, eo = field("Z,b,c,d,e", sp.QQ)
    def Fgen(k, A, B, K):
        p = K(1)
        for i, j in enumerate(L.PERMS[k]):
            if A[i] is not None:
                p = p * (A[i] + B[j])
        return p
    A1 = [None, K1(0), K1(1), Zo]; B1 = [bo, co, do, eo]
    F = lambda k: Fgen(k, A1, B1, K1)
    E1 = F(0) - F(6)
    kchk("E = (c-b)(1+d)(Z+e)", E1 == (co - bo) * (1 + do) * (Zo + eo))
    kchk("rE", F(7) - F(6) == bo * (eo - do) * (Zo - 1))
    kchk("sE", F(11) - F(6) == (1 + do) * Zo * (eo - bo))
    kchk("uE", F(12) - F(6) == bo * (co - do) * (Zo + eo))
    print("perm one-line 0,3,6,7,11,12:", [L.perm_to_oneline(L.PERMS[k]) for k in (0, 3, 6, 7, 11, 12)])
    rf = (F(7) - F(6)) / E1; sf = (F(11) - F(6)) / E1; uf = (F(12) - F(6)) / E1
    Vf = (co - bo) / (bo - do)
    T0f = (rf - sf) * Zo + sf
    N0f = uf * Zo * (Zo - 1) - rf * Zo**2
    N1f = uf * sf * (Zo - 1)**2 - rf * Zo**2
    Yf = (N0f + Vf * N1f) / T0f
    kchk("b(V+1)T0 = N0 + V N1 (original variables)", bo * (Vf + 1) * T0f == N0f + Vf * N1f)
    kchk("inverse: b = Y/(V+1)", Yf / (Vf + 1) == bo)
    kchk("inverse: d = Y/(uV) - 1", Yf / (uf * Vf) - 1 == do)
    kchk("inverse: c = V + (u-1)Y/u", Vf + (uf - 1) * Yf / uf == co)
    kchk("inverse: e = Z[b+s(c-b)]/(Z-s(c-b))", Zo * (bo + sf * (co - bo)) / (Zo - sf * (co - bo)) == eo)

    # ---- chart variables
    K, Z, V, r, s, u, h = field("Z,V,r,s,u,h", sp.QQ)
    T0 = (r - s) * Z + s; N0 = u * Z * (Z - 1) - r * Z**2; N1 = u * s * (Z - 1)**2 - r * Z**2
    Y = (N0 + V * N1) / T0
    b = Y / (V + 1); d = Y / (u * V) - 1; c = V + (u - 1) * Y / u
    e = Z * (b + s * (c - b)) / (Z - s * (c - b))
    A = [None, K(0), K(1), Z]; B = [b, c, d, e]
    Fc = {k: Fgen(k, A, B, K) for k in range(24)}
    Ec = Fc[0] - Fc[6]
    H = {k: (Fc[k] - Fc[6]) / Ec for k in range(24)}
    for k, v in [(6, 0), (0, 1), (7, r), (11, s), (12, u)]:
        kchk("H%d = %s" % (k, v), H[k] - v == 0)
    kchk("chart -> V recovered", (c - b) / (b - d) == V)
    G3 = H[3] - h
    C41 = G3.numer; D3 = G3.denom
    Cs = sp.Poly(C41.as_expr(), *sp.symbols("Z V r s u h"))
    print("H3-h numerator terms:", len(Cs.terms()), " bidegree (Z,V):", sp.Poly(C41.as_expr(), sp.Symbol("Z"), sp.Symbol("V")).degree(sp.Symbol("Z")), sp.Poly(C41.as_expr(), sp.Symbol("Z"), sp.Symbol("V")).degree(sp.Symbol("V")))
    print("H3-h denominator:", sp.factor(D3.as_expr()))
    kchk("H8 = -u/Y", H[8] == -u / Y)
    kchk("H19 = Y[(r-su)Z+su]/(uZ(Z-1))", H[19] == Y * ((r - s * u) * Z + s * u) / (u * Z * (Z - 1)))
    # image curve: substitute V = (Y T0 - N0)/N1
    K2, Z2, Y2, r2, s2, u2, h2 = field("Z,Y,r,s,u,h", sp.QQ)
    T02 = (r2 - s2) * Z2 + s2; N02 = u2 * Z2 * (Z2 - 1) - r2 * Z2**2; N12 = u2 * s2 * (Z2 - 1)**2 - r2 * Z2**2
    Vsub = (Y2 * T02 - N02) / N12
    G3img = K2(sp.sympify(G3.as_expr()).subs({sp.Symbol("V"): sp.Symbol("Vtmp")}).subs(sp.Symbol("Vtmp"), Vsub.as_expr()))
    P75 = G3img.numer
    Pe = sp.expand(P75.as_expr())
    Zs, Ys = sp.symbols("Z Y")
    PP = sp.Poly(Pe, Zs, Ys)
    print("P75 terms:", len(sp.Poly(Pe, *sp.symbols("Z Y r s u h")).terms()), " bidegree (Z,Y):", PP.degree(Zs), PP.degree(Ys))
    rr, ss, uu, hh = sp.symbols("r s u h")
    P0 = sp.factor(Pe.subs(Ys, 0))
    target = Zs * rr * uu**2 * (Zs - 1) * ((ss - 1) * Zs - ss) * ((rr - ss * uu) * Zs + uu * (ss + 1))
    ratio = sp.cancel(P0 / target)
    print("P75(Z,0) =", P0, " ratio to doc's (10):", ratio)
    const = ratio  # normalisation constant
    kchk("P75(Z,0) proportional to doc (10) by a constant", ratio.free_symbols == set())
    print("[Z^4] P75(Z,0) / const:", sp.factor(sp.Poly(Pe.subs(Ys, 0), Zs).coeff_monomial(Zs**4) / const))
    for pt, nm in [((0, 0), "(11)"), ((1, 0), "(12)")]:
        PZ = sp.factor(sp.diff(Pe, Zs).subs({Zs: pt[0], Ys: pt[1]}) / const)
        PY = sp.factor(sp.diff(Pe, Ys).subs({Zs: pt[0], Ys: pt[1]}) / const)
        print("  derivatives at", pt, nm, ":", PZ, ",", PY)
    # divisibility: substitute Y = (N0 + V N1)/T0 into P75 -> numerator divisible by C41
    P75K = K(sp.sympify(Pe).subs(Ys, sp.Symbol("Ytmp")).subs(sp.Symbol("Ytmp"), Y.as_expr()))
    q = P75K.numer
    Cnum = C41
    Rq = sp.rem(sp.Poly(q.as_expr(), *sp.symbols("V Z r s u h")), sp.Poly(Cnum.as_expr(), *sp.symbols("V Z r s u h")))
    print("remainder of (P75 o Y) numerator by C41 (as polys in V first):", Rq.is_zero)
    quo = sp.div(sp.Poly(q.as_expr(), *sp.symbols("Z V r s u h")), sp.Poly(Cnum.as_expr(), *sp.symbols("Z V r s u h")))
    kchk("C41 divides numerator of P75(Z, Y(Z,V))", quo[1].is_zero)
    # discriminant witness
    Vs = sp.Symbol("V")
    C0 = sp.expand(C41.as_expr()).subs({rr: 2, ss: 3, uu: 5, hh: 7})
    disc = sp.factor(sp.discriminant(sp.Poly(C0, Vs)))
    print("Disc_V C41 at (2,3,5,7):", disc)
    dp = sp.Poly(sp.discriminant(sp.Poly(C0, Vs)), Zs)
    print("  degree", dp.degree(), " gcd with derivative:", sp.gcd(dp, dp.diff(Zs)))
    cond9 = rr*ss*uu*(ss-1)*(ss+1)*(rr-ss*uu)*(rr+uu)*(rr*ss-uu)*(hh-rr-ss-1)*(hh-rr-uu-1)
    print("  (9) at (2,3,5,7):", cond9.subs({rr: 2, ss: 3, uu: 5, hh: 7}))
    # universal constant stratum
    coeffs = sp.Poly(sp.expand(C41.as_expr()), Zs, Vs).coeffs()
    lam = sp.Symbol("lam")
    G = sp.groebner(list(coeffs) + [lam * rr * ss * uu - 1], lam, rr, ss, uu, hh, order="lex")
    print("coefficient ideal saturated: ", list(G))
    H22 = sp.cancel(H[22].as_expr().subs({rr: 1, ss: -1, uu: -1, hh: 1}))
    print("H22 on (13):", H22)
    print("PART chart ALLOK", all(res))
    return all(res)

print('==== part chart', flush=True)
PARTS.append(('chart', part_chart()))

def numeric_stratum(T_arg, s0, s1, margin_arg):
    import numpy as np
    T = [int(v) for v in T_arg.split(",")]
    seeds = range(int(s0), int(s1))
    margin = float(margin_arg)
    Pn = L.P_numpy()
    worst = None; found = 0
    t0 = time.time()
    for sd in seeds:
        pt = L.random_point_on_contact_stratum(T, seed=sd, tries=30, margin=margin)
        if pt is None:
            print("seed", sd, "no point"); continue
        found += 1
        m = L.matching_vector_np(pt["a"], pt["b"], pt["alpha"])
        Fps = (Pn @ m)[pt["pi_star"]]
        c = m.real.copy(); c[13] -= Fps.real
        g = Pn @ c
        val, O = L.min_phi_over_SO4(c, restarts=25, seed=sd)
        # normalise by sum of gaps
        sc = val / g.sum()
        print("seed", sd, "contact", pt["contact"], "min gap off T %.3g" % min(g[j] for j in range(24) if j not in T),
              "min Phi/sum(gaps) = %.3e" % sc, flush=True)
        if worst is None or sc < worst[0]:
            worst = (sc, sd)
    print("T", T, "points", found, "worst normalised min Phi", worst, "time %.0fs" % (time.time() - t0))

def numeric_chart(seed_arg, n_arg):
    """NUMERIC falsification for T33 = {0,3,6,7,11,12} using the document's own chart."""
    import numpy as np
    from scipy.optimize import minimize, brentq
    T = [0, 3, 6, 7, 11, 12]; OFF = [j for j in range(24) if j not in T]
    PERMS = np.array(L.PERMS)
    Pn = L.P_numpy(); Pp = np.linalg.pinv(Pn)

    def sigma_pts(Z, V, r, s, u):
        T0 = (r - s) * Z + s; N0 = u * Z * (Z - 1) - r * Z**2; N1 = u * s * (Z - 1)**2 - r * Z**2
        Y = (N0 + V * N1) / T0
        b = Y / (V + 1); d = Y / (u * V) - 1; c = V + (u - 1) * Y / u
        e = Z * (b + s * (c - b)) / (Z - s * (c - b))
        A = [None, 0.0, 1.0, Z]; B = [b, c, d, e]
        F = np.ones(24, complex)
        for k, p in enumerate(L.PERMS):
            for i in range(4):
                if A[i] is not None:
                    F[k] *= A[i] + B[p[i]]
        return F

    def H(x):
        r, s, u, Zr, Zi, Vr, Vi = x
        F = sigma_pts(Zr + 1j * Zi, Vr + 1j * Vi, r, s, u)
        return (F - F[6]) / (F[0] - F[6])

    def obj(x, eta):
        with np.errstate(all="ignore"):
            Hh = H(x)
        if not np.all(np.isfinite(Hh)): return 1e6
        g = -eta * Hh.imag
        scale = np.mean(np.abs(g[OFF])) + 1e-12
        return 1e4 * (Hh[3].imag / (1 + abs(Hh[3])))**2 - min(g[j] for j in OFF) / scale

    rng = np.random.default_rng(int(seed_arg))
    nwant = int(n_arg); got = 0; worst = None; tries = 0
    while got < nwant and tries < 400:
        tries += 1
        eta = rng.choice([-1.0, 1.0])
        x0 = rng.normal(size=7) * 1.5
        rres = minimize(obj, x0, args=(eta,), method="Nelder-Mead", options=dict(maxiter=4000, xatol=1e-12, fatol=1e-14))
        x = rres.x
        with np.errstate(all="ignore"):
            Hh = H(x)
        if not np.all(np.isfinite(Hh)): continue
        g = -eta * Hh.imag
        if abs(Hh[3].imag) > 1e-8 * (1 + abs(Hh[3])) or min(g[j] for j in OFF) <= 1e-6 * np.mean(np.abs(g[OFF])):
            continue
        if abs(x[3 + 1]) < 1e-6:   # Z real: one-circle case, skip
            continue
        F = sigma_pts(x[3] + 1j * x[4], x[5] + 1j * x[6], x[0], x[1], x[2])
        alpha = 1j * eta / (F[0] - F[6])
        m = Pp @ (alpha * F)
        c = m.real.copy(); c[13] -= (alpha * F[6]).real
        gg = Pn @ c
        val, O = L.min_phi_over_SO4(c, restarts=20, seed=tries)
        sc = val / gg.sum()
        got += 1
        print("pt", got, "r,s,u=%.3f,%.3f,%.3f Z=%.3f%+.3fi" % (x[0], x[1], x[2], x[3], x[4]),
              "contacts(|g|<1e-9)", [j for j in range(24) if abs(gg[j]) < 1e-9 * gg.sum()],
              "min off-T gap/sum %.2e" % (min(gg[j] for j in OFF) / gg.sum()), "min Phi/sum %.3e" % sc, flush=True)
        if worst is None or sc < worst: worst = sc
    print("points", got, "tries", tries, "worst normalised min Phi", worst)

if NUMERIC:
    numeric_stratum('0,3,6,7,11,12', '200', '206', '0.01')
    numeric_chart('1', '3')

print("==== summary")
for nm, ok in PARTS: print("part", nm, "OK" if ok else "FAIL")
ok = all(ok for _, ok in PARTS)
print("EXACT CHECKS", "ALL PASS" if ok else "SOME FAILED")
sys.exit(0 if ok else 1)
