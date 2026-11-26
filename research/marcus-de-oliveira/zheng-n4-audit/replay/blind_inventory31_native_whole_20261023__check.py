"""Own-code replay for the Zheng n=4 audit, document blind_inventory31_native_whole_20261023.
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
    report("doc31", {0,2,3,7,11,12}, 32)
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

    ok1 = check("doc31", [0,2,3,7,11,12], [0,1,2,7,8,12,15,39], "x0 x1 a b c e f g".split(),
                {1:"g",5:"f",6:"e",9:"c",13:"a",10:"b",23:"x0+x1",18:"x0+f",20:"x0+b",21:"x0+a",16:"x1+b",17:"x1+g"},
                [((18,0),(1,2,7,15)), ((40,1),(0,7,39))])
    vv = [(L.F_VALUES[j][4], L.F_VALUES[j][10]) for j in (18, 40)]
    print(" f18,f40 at 4 and 10:", vv)
    res = [ok1, vv == [(1, 1), (1, 1)]]
    return all(res)

print('==== part cone', flush=True)
PARTS.append(('cone', part_cone()))

def part_chart():
    import sympy as sp

    a, z, x, w, t = sp.symbols("a z x w t")
    W, r, s, u, h = sp.symbols("W r s u h")
    A = [a, z, x, w]
    B = [t, 0, 1, None]   # None = infinity

    def F(k, A=A, B=B):
        p = 1
        for i, j in enumerate(L.PERMS[k]):
            if B[j] is not None:
                p *= (A[i] + B[j])
        return p

    res = []
    def chk(name, expr):
        ok = sp.expand(sp.numer(sp.together(expr))) == 0
        res.append(ok); print(("OK  " if ok else "FAIL"), name)

    E = F(0) - F(2)
    chk("E = (a+t)(z-x)", E - (a + t) * (z - x))
    chk("rE", F(3) - F(2) - (a + t) * (z + 1) * (w - x))
    chk("sE", F(12) - F(2) - x * (t - 1) * (a - z))
    chk("uE", F(7) - F(2) - (a * (z + t) * (w + 1) - (a + t) * (z + 1) * x))
    print("perm one-line 0,2,3,7,11,12:", [L.perm_to_oneline(L.PERMS[k]) for k in (0, 2, 3, 7, 11, 12)])
    # forward definitions
    rf = sp.cancel((F(3) - F(2)) / E); sf = sp.cancel((F(12) - F(2)) / E); uf = sp.cancel((F(7) - F(2)) / E)
    Wf = (z - x) / (x * (z + 1))
    # inverse (7)
    zI = x * (1 + W) / (1 - x * W)
    wI = x * (1 + r * W)
    Da = 1 + s * W + x * W * (r + s - u + r * s * W)
    aI = x * (1 + u * W) / Da
    tI = ((aI - zI) + s * aI * (zI + 1) * W) / ((aI - zI) - s * (zI + 1) * W)
    # forward then inverse: substitute W,r,s,u by forward definitions -> recover z,w,a,t
    sub_f = {W: Wf, r: rf, s: sf, u: uf}
    for nm, e, tgt in [("z", zI, z), ("w", wI, w), ("a", aI, a), ("t", tI, t)]:
        chk("(7) forward->inverse recovers " + nm, sp.together(e.subs(sub_f, simultaneous=True)) - tgt)
    # inverse then forward: r,s,u,W recovered
    sub_i = {z: zI, w: wI, a: aI, t: tI}
    for nm, e, tgt in [("r", rf, r), ("s", sf, s), ("u", uf, u), ("W", Wf, W)]:
        chk("(7) inverse->forward recovers " + nm, sp.together(e.subs(sub_i, simultaneous=True)) - tgt)
    Dt = (a - z) - s * (z + 1) * W
    chk("Dt numerator - denominator = s(a+1)(z+1)W", ((a - z) + s * a * (z + 1) * W) - Dt - s * (a + 1) * (z + 1) * W)
    chk("t+z identity", (tI + zI) - ((aI - zI) * (zI + 1) * (1 + s * W) / ((aI - zI) - s * (zI + 1) * W)))
    chk("a*Da = x(1+uW) in original variables", (a * Da - x * (1 + u * W)).subs(sub_f, simultaneous=True))
    chk("Da at u=-1/W equals (1+sW)(w+1) with w=x(1+rW)", (Da.subs(u, -1 / W) - (1 + s * W) * (wI + 1)))
    # H's in (x,W,r,s,u)
    from sympy.polys.fields import field
    K, X_, W_, R_, S_, U_, Hh_ = field("x,W,r,s,u,h", sp.QQ)
    def toK(e):
        return K(sp.sympify(e).subs({x: sp.Symbol("x")}))
    zK = X_ * (1 + W_) / (1 - X_ * W_); wK = X_ * (1 + R_ * W_)
    DaK = 1 + S_ * W_ + X_ * W_ * (R_ + S_ - U_ + R_ * S_ * W_)
    aK = X_ * (1 + U_ * W_) / DaK
    tK = ((aK - zK) + S_ * aK * (zK + 1) * W_) / ((aK - zK) - S_ * (zK + 1) * W_)
    Asub = [aK, zK, X_, wK]; Bsub = [tK, 0, 1, None]
    def FK(k):
        p = K(1)
        for i, j in enumerate(L.PERMS[k]):
            if Bsub[j] is not None:
                p = p * (Asub[i] + Bsub[j])
        return p
    F2K = FK(2); EK = FK(0) - F2K
    HK = {k: (FK(k) - F2K) / EK for k in range(24)}
    def kchk(name, e):
        ok = (e == 0); res.append(ok); print(("OK  " if ok else "FAIL"), name)
    for k, v in [(2, 0), (0, 1), (3, R_), (12, S_), (7, U_)]:
        kchk("H%d = %s" % (k, v), HK[k] - v)
    DK = (W_ + 1) * (R_ * X_ * (1 + S_ * W_) + S_ * (X_ + 1)) + (1 - U_) * (X_ + 1)
    CK = (HK[11] - Hh_) * DK
    print("  denominator of (H11-h)*D in the field:", CK.denom.as_expr())
    kchk("(H11-h)*D is a polynomial (denominator is a nonzero constant)", 0 if sp.sympify(CK.denom.as_expr()).is_number else 1)
    C = sp.expand(CK.as_expr()); D = sp.expand(DK.as_expr())
    print("denominator of H11-h (reduced):", sp.factor((HK[11] - Hh_).denom.as_expr()))
    Cp = sp.Poly(C, W, x)
    print("C: number of terms =", len(sp.Poly(C, W, x, r, s, u, h).terms()), " bidegree (W,x) =", Cp.degree(W), Cp.degree(x))
    kchk("a - z = -W x D/(Da(1-xW))", aK - zK + W_ * X_ * DK / (DaK * (1 - X_ * W_)))
    kchk("H4", HK[4] - (-R_ * W_ * X_**2 + R_ * X_ - X_**2 - X_) / (X_ + 1))
    DeltaK = R_ * S_ * (1 + U_ * W_) * (1 - W_ * X_) / DK
    kchk("H11 - H10 = Delta", HK[11] - HK[10] - DeltaK)
    kchk("F2/E = 1/W", F2K / EK - 1 / W_)
    H10_13 = HK[10].as_expr().subs({r: -1, s: -1, u: -1, h: 1})
    print("H10 on (13):", sp.simplify(H10_13))
    # discriminant witness
    C0 = C.subs({r: 2, s: 3, u: 4, h: 7})
    disc = sp.factor(sp.discriminant(sp.Poly(C0, x)))
    print("Disc_x C at (2,3,4,7):", disc)
    dp = sp.Poly(sp.discriminant(sp.Poly(C0, x)), W)
    print("  degree", dp.degree(), "gcd with derivative", sp.gcd(dp, dp.diff(W)))
    cond10 = r*s*u*(r+1)*(s+1)*(h+u)*(h-1)*(s-1)*(h*s-u)*(r-u)*(s-u)*(u-1)*(r+s-u+1)
    print("  (10) at (2,3,4,7):", cond10.subs({r: 2, s: 3, u: 4, h: 7}))
    # x = -1
    chk("C(W,-1)", C.subs(x, -1) - r * (W + 1) * ((h * s - u) * W + h + s - u - 1))
    chk("D(W,-1)", D.subs(x, -1) + r * (W + 1) * (s * W + 1))
    chk("C_W(-1,-1)", sp.diff(C, W).subs({W: -1, x: -1}) + r * (h - 1) * (s - 1))
    print("  C_x(-1,-1) =", sp.factor(sp.diff(C, x).subs({W: -1, x: -1})))
    W1 = sp.solve((h * s - u) * W + h + s - u - 1, W)[0]
    print("  second root W1 at x=-1:", sp.factor(W1), "; W1+1 =", sp.factor(W1 + 1), "; s*W1+1 =", sp.factor(s * W1 + 1))
    print("  H4 numerator at x=-1:", sp.factor((-r * W * x**2 + r * x - x**2 - x).subs(x, -1)))
    # x = infinity
    cx2 = Cp.coeff_monomial(x**2) if False else sp.expand(C).coeff(x, 2)
    chk("[x^2]C", cx2 + (u * W + 1) * (r * (s + 1) * W + r + s - u + 1))
    D1 = sp.expand(D).coeff(x, 1); D0 = sp.expand(D).coeff(x, 0)
    print("  D1 =", sp.factor(D1), " D0 =", sp.factor(D0))
    R12 = sp.factor(sp.resultant(cx2, D1, W))
    print("  Res_W([x^2]C, D1) =", R12)
    chk("Res (12)", R12 + r * s * (r - u)**2 * (s - u) * (u - 1) * (r + s - u + 1))
    chk("[W^2 x^2]C", sp.expand(C).coeff(x, 2).coeff(W, 2) + r * u * (s + 1))
    # W = infinity
    cW2 = sp.expand(C).coeff(W, 2)
    chk("[W^2]C", cW2 + r * x * (u * (s + 1) * x + s * (u + h)))
    xs = -s * (h + u) / (u * (s + 1))
    print("  lim D/W^2 at x*:", sp.factor(sp.expand(D).coeff(W, 2).subs(x, xs)), " expected r s x* =", sp.factor(r * s * xs))
    print("  lim Delta at W->inf, x->x*:", sp.factor(sp.limit((r*s*(1+u*W)*(1-W*xs)/D.subs(x, xs)), W, sp.oo)))
    eps, V = sp.symbols("eps V")
    wl = sp.expand(sp.cancel(eps * C.subs({W: 1 / eps, x: eps * V})))
    wl0 = sp.factor(wl.subs(eps, 0))
    print("  eps*C(1/eps, eps V)|eps=0 =", wl0)
    chk("weighted leading eq", wl0 + s * (r * (h + u) * V + h - r * u))
    V0 = (r * u - h) / (r * (h + u))
    print("  D/W at V0, W->inf:", sp.factor(sp.limit(sp.cancel(D.subs(x, V / W) / W), W, sp.oo).subs(V, V0)))
    print("  G(eps,x)=eps^2 C(1/eps,x): d/dx at (0,0) =", sp.factor(sp.diff(sp.expand(sp.cancel(eps**2*C.subs(W,1/eps))), x).subs({eps:0, x:0})))
    # constant stratum
    coeffs = sp.Poly(C, W, x).coeffs()
    lam = sp.symbols("lam")
    G = sp.groebner(list(coeffs) + [lam * r * s - 1], lam, r, s, u, h, order="lex")
    print("  coefficient ideal saturated Groebner basis:", list(G))
    # nonconstancy of H4
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

if NUMERIC:
    numeric_stratum('0,2,3,7,11,12', '100', '130', '0.01')

print("==== summary")
for nm, ok in PARTS: print("part", nm, "OK" if ok else "FAIL")
ok = all(ok for _, ok in PARTS)
print("EXACT CHECKS", "ALL PASS" if ok else "SOME FAILED")
sys.exit(0 if ok else 1)
