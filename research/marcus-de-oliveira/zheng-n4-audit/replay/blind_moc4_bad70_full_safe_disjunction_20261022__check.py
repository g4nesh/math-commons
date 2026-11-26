"""Own-code replay for the Zheng n=4 audit, document blind_moc4_bad70_full_safe_disjunction_20261022.
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
    report("doc70", {2,3,5,9,10,12}, 34)
    return True

print('==== part contacts', flush=True)
PARTS.append(('contacts', part_contacts()))

def part_radial():
    """Exact replay for blind_moc4_bad70_full_safe_disjunction_20261022 (paper case T_34).
    Supplement JSON files are loaded as DATA only (json.load); nothing of theirs is executed."""
    import sympy as sp
    import flint

    res = []
    def chk(name, ok):
        res.append(bool(ok)); print(("OK  " if ok else "FAIL"), name, flush=True)

    T = [2, 3, 5, 9, 10, 12]
    GENS = [0, 1, 2, 12, 39, 67, 70, 72]
    chk("generators(T) == doc list (f0,f1,f2,f12,f39,f67,f70,f72)", L.generators(T) == GENS)
    chk("only f70 is non-minor", [j for j in GENS if not L.F_IS_MINOR[j]] == [70])
    print("names:", [L.f_name(j) for j in GENS])

    # ---------- restricted cubic (2), (3)
    p, r, a, b, c, q, d, s, H = sp.symbols("p r a b c q d s H")
    coef = [p, r, a, b, c, q, d, s]
    v = [sum(cf * L.F_VECS[j][k] for cf, j in zip(coef, GENS)) for k in range(14)]
    vals = [sp.expand(sum(L.P[t][k] * v[k] for k in range(14))) for t in range(24)]
    chk("P v = sum coef * profile and vanishes on T", all(vals[t] == 0 for t in T))
    vH = list(v); vH[13] = vH[13] + H
    CH = sp.expand(L.cubic_C(vH))
    Cv = sp.expand(CH.subs(H, 0)); Qv = sp.expand(sp.diff(CH, H))
    chk("C(v + H kappa) is affine in H", sp.diff(CH, H, 2) == 0)
    U = d + r + s; M = b + c
    Qdoc = (a - U) * M - d * U
    Cdoc = b * c * (a - U) - b * (d + r) * (d + q + s) - c * (d + s) * (d + p + r) - d * (d + p + r) * (d + q + s)
    chk("(2) Q = (a-U)M - dU", sp.expand(Qv - Qdoc) == 0)
    chk("(3) C formula", sp.expand(Cv - Cdoc) == 0)
    chk("M C - b c Q = bcdU - M[...] (Q=0 exclusion identity)",
        sp.expand(M * Cdoc - b * c * Qdoc - (b * c * d * U - M * (b * (d + r) * (d + q + s) + c * (d + s) * (d + p + r) + d * (d + p + r) * (d + q + s)))) == 0)

    # ---------- R9 from the supplement JSON (data) and the Hessian identity at random rational points
    R9j = json.load(open(os.path.join(SUPP, "data", "blind_moc4_next_singlebad_R9_restriction_20261022.json")))
    names = R9j["variables"]
    symmap = dict(p=p, r=r, a=a, b=b, c=c, q=q, d=d, s=s, h=sp.Symbol("h"))
    assert len(R9j["factors"]) == 1 and R9j["factors"][0]["power"] == 1 and str(R9j["unit"]) == "1"
    R9 = sp.Integer(0)
    for exps, cf in R9j["factors"][0]["terms"]:
        t = sp.Rational(cf)
        for nm, e in zip(names, exps):
            t *= symmap[nm] ** e
        R9 += t
    R9P = sp.Poly(R9, *coef)
    degs = {sum(m) for m in R9P.monoms()}
    print("R9: terms", len(R9P.terms()), "total degrees", degs, "uses h:", sp.Symbol("h") in R9.free_symbols)
    ms = sp.symbols("m0:14")
    Csym = sp.expand(L.cubic_C(ms))
    Hess = [[sp.diff(Csym, ms[i], ms[j]) for j in range(14)] for i in range(14)]
    Hf = [[sp.lambdify(ms, Hess[i][j]) for j in range(14)] for i in range(14)]
    vf = sp.lambdify(coef, v); Cf = sp.lambdify(coef, Cdoc); Qf = sp.lambdify(coef, Qdoc); Rf = sp.lambdify(coef, R9)
    rng = random.Random(2026)
    okall = True; npts = 12
    for it in range(npts):
        pt = [Fr(rng.randint(-40, 40), rng.randint(1, 9)) for _ in range(8)]
        Qx = Qf(*pt)
        if Qx == 0: continue
        m = [Fr(x) for x in vf(*pt)]
        m[13] -= Fr(Cf(*pt)) / Qx
        Hm = [[Fr(Hf[i][j](*m)) for j in range(14)] for i in range(14)]
        lhs = Qx**3 * L.det(Hm)
        rhs = -4 * Qx * Fr(Rf(*pt))**2
        okall &= (lhs == rhs)
    chk("Q^3 det Hess C(v - (C/Q) kappa) = -4 Q R9^2 at %d random rational points (exact arithmetic)" % npts, okall)
    # paper p.27 argument: Re x satisfies C = 0 and det Hess = 0 -- check at a Gaussian-rational point of X
    mx = L.matching_vector([L.QI(1, 2), L.QI(-3, 1), L.QI(2, -1), L.QI(0, 5)], [L.QI(4, 1), L.QI(-1, -1), L.QI(3, 3), L.QI(-2, 0)], L.QI(1, 1))
    u_re = [L.re_part(z) for z in mx]
    Hm = [[Fr(sp.Rational(Hess[i][j].subs(dict(zip(ms, u_re))))) for j in range(14)] for i in range(14)]
    chk("C(Re x) = 0 and det Hess C(Re x) = 0 at a Gaussian-rational matching point", L.cubic_C(u_re) == 0 and L.det(Hm) == 0)

    # ---------- Lemma 2 structural claims on k0, k1, k2 (d = 1)
    R1 = sp.Poly(R9.subs(d, 1), a)
    k = {i: sp.Poly(R1.coeff_monomial(a**i), p, r, b, c, q, s) for i in range(3)}
    print("deg_a R9 =", R1.degree())
    chk("k0 has only positive coefficients and constant term 1",
        all(cf > 0 for cf in k[0].coeffs()) and k[0].coeff_monomial(1) == 1)
    k2_0 = k[2].as_expr().subs({b: 0, c: 0}); k1_0 = sp.Poly(k[1].as_expr().subs({b: 0, c: 0}), p, r, q, s)
    chk("b=c=0: k2 = 0 and k1 has positive coefficients", sp.expand(k2_0) == 0 and all(cf > 0 for cf in k1_0.coeffs()))
    global KJ
    KJ = {i: str(k[i].as_expr()) for i in range(3)}
    print("PART radial ALLOK", all(res))
    return all(res)

print('==== part radial', flush=True)
PARTS.append(('radial', part_radial()))

def part_gram():
    """Exact replay of the (PB) and (CQ) Gram identities of the T_34 document.
    The cubic feature polynomials, the four signed variable permutations and the 8x8 matrix T are
    read from supplement JSON as DATA; the identities are checked against Leibniz sums built here
    from the paper's generator profiles (moc4lib), modulo the unit-sphere relations."""
    import flint
    import sympy as sp

    res = []
    def chk(name, ok):
        res.append(bool(ok)); print(("OK  " if ok else "FAIL"), name, flush=True)

    names = ("p0", "p1", "p2", "p3", "q0", "q1", "q2", "q3")
    ctx = flint.fmpq_mpoly_ctx.get(names, "lex")
    X = ctx.gens()
    w, a1, b1, c1 = X[:4]; v0, d1, e1, f1 = X[4:]
    LL = [[w, -a1, -b1, -c1], [a1, w, -c1, b1], [b1, c1, w, -a1], [c1, -b1, a1, w]]
    RR = [[v0, -d1, -e1, -f1], [d1, v0, f1, -e1], [e1, -f1, v0, d1], [f1, e1, -d1, v0]]
    A = [[sum((LL[i][k] * RR[k][j] for k in range(4)), ctx.from_dict({})) for j in range(4)] for i in range(4)]
    npp = sum(x * x for x in X[:4]); nqq = sum(x * x for x in X[4:])
    AtA = [[sum((A[k][i] * A[k][j] for k in range(4)), ctx.from_dict({})) for j in range(4)] for i in range(4)]
    chk("A^T A = |p|^2 |q|^2 I", all(AtA[i][j] == (npp * nqq if i == j else 0) for i in range(4) for j in range(4)))
    detA = L.det_leibniz4(A)
    chk("det A = (|p|^2 |q|^2)^2 (so A in SO(4) on the unit spheres)", detA == (npp * nqq) ** 2)

    Sp = 1 - X[1]**2 - X[2]**2 - X[3]**2
    Sq = 1 - X[5]**2 - X[6]**2 - X[7]**2
    def NF(poly):
        """normal form modulo (p0^2 - Sp, q0^2 - Sq) (a Groebner basis: coprime leading terms)"""
        out = ctx.from_dict({})
        cache = {}
        for exps, cf in poly.to_dict().items():
            e0, f0 = exps[0], exps[4]
            rest = list(exps); rest[0] = e0 % 2; rest[4] = f0 % 2
            key = (e0 // 2, f0 // 2)
            if key not in cache:
                cache[key] = Sp ** key[0] * Sq ** key[1]
            out += ctx.from_dict({tuple(rest): cf}) * cache[key]
        return out

    W = [[A[i][s[i]] for i in range(4)] for s in L.PERMS]
    mono = []
    for k, s in enumerate(L.PERMS):
        m = L.SGN[k] * A[0][s[0]] * A[1][s[1]] * A[2][s[2]] * A[3][s[3]]
        mono.append(m)
    def Phi(j):
        return NF(sum((L.F_VALUES[j][k] * mono[k] for k in range(24) if L.F_VALUES[j][k]), ctx.from_dict({})))
    GENS = [0, 1, 2, 12, 39, 67, 70, 72]
    PH = {j: Phi(j) for j in [0, 2, 12, 39, 67, 70]}
    print("Leibniz sums built", flush=True)

    data = json.load(open(os.path.join(SUPP, "data", "blind_moc4_next_singlebad_cubic_reduction_20261022.json")))
    chg = json.load(open(os.path.join(SUPP, "data", "blind_moc4_next_singlebad_budget_PSD8_20261022.json")))
    feat = data["features"]; left = data["left_four"]
    if isinstance(feat, str): feat = json.loads(feat.replace("'", '"'))
    if isinstance(left, str): left = json.loads(left)
    Tm = chg["T"]
    if isinstance(Tm, str): Tm = json.loads(Tm.replace("'", '"'))
    Tm = sp.Matrix([[sp.Rational(x) for x in row] for row in Tm])
    Ti = Tm.inv()
    loc = {nm: g for nm, g in zip(names, X)}
    def parse(txt):
        e = sp.sympify(txt)
        pe = sp.Poly(e, *sp.symbols(names))
        return ctx.from_dict({tuple(m): flint.fmpq(int(sp.Rational(cf).p), int(sp.Rational(cf).q)) for m, cf in pe.terms()})
    seeds = [[parse(t) for t in block] for block in feat]
    def acted(poly, g):
        out = {}
        for exps, cf in poly.to_dict().items():
            e = [0] * 8; val = cf
            for power, (dest, sg) in zip(exps, g):
                e[dest] += power
                if sg == -1 and power % 2: val = -val
            e = tuple(e); out[e] = out.get(e, 0) + val
        return ctx.from_dict(out)
    blocks = []
    for block in seeds:
        for g in left:
            vv = [acted(f, g) for f in block]
            blocks.append([sum((flint.fmpq(int(Ti[i, j].p), int(Ti[i, j].q)) * vv[j] for j in range(8) if Ti[i, j] != 0), ctx.from_dict({})) for i in range(8)])
    print("feature blocks:", len(blocks), flush=True)
    def G(i, j):
        return NF(sum((bl[i] * bl[j] for bl in blocks), ctx.from_dict({})))
    G00, G22, G33, G44, G55, G66 = G(0, 0), G(2, 2), G(3, 3), G(4, 4), G(5, 5), G(6, 6)
    G02, G23, G24, G25 = G(0, 2), G(2, 3), G(2, 4), G(2, 5)
    lam = 2
    chk("diag G00 = 2 Phi(f0)", G00 == lam * PH[0])
    chk("diag G22 = 2 Phi(f2)", G22 == lam * PH[2])
    chk("diag G33 = 2 Phi(f12)", G33 == lam * PH[12])
    chk("diag G44 = 2 Phi(f39)", G44 == lam * PH[39])
    chk("diag G55 = 2 Phi(f67)", G55 == lam * PH[67])
    chk("(PB) G66 - 2 G02 - 2 G23 = 2 Phi(f70)", G66 - 2 * G02 - 2 * G23 == lam * PH[70])
    chk("(CQ) G66 + 2 G24 - 2 G25 = 2 Phi(f70)", G66 + 2 * G24 - 2 * G25 == lam * PH[70])
    pp, aa, bb, cc, qq, dd = sp.symbols("p a b c q d", positive=True)
    B3 = sp.Matrix([[pp, -dd, 0], [-dd, aa, -dd], [0, -dd, bb]])
    C3 = sp.Matrix([[aa, dd, -dd], [dd, cc, 0], [-dd, 0, qq]])
    chk("Schur: det[[p,-d,0],[-d,a,-d],[0,-d,b]] = p b (a - d^2/p - d^2/b)", sp.expand(B3.det() - pp * bb * (aa - dd**2 / pp - dd**2 / bb)) == 0)
    chk("Schur: det[[a,d,-d],[d,c,0],[-d,0,q]] = c q (a - d^2/c - d^2/q)", sp.expand(C3.det() - cc * qq * (aa - dd**2 / cc - dd**2 / qq)) == 0)
    print("PART gram ALLOK", all(res))
    return all(res)

print('==== part gram', flush=True)
PARTS.append(('gram', part_gram()))

def part_bernstein():
    """Exact replay of the Bernstein coefficient-sign certificate of Lemma 2 (T_34 document).
    Input: k0, k1, k2 of R9 (d = 1), produced by d70.py from the R9 polynomial."""
    import flint, sympy as sp

    kj = KJ
    src = sp.symbols("p r b c q s")
    kp = {int(i): sp.Poly(sp.sympify(e), *src) for i, e in kj.items()}
    deg = [max(kp[i].degree(x) for i in range(3)) for x in src]
    ctx = flint.fmpz_mpoly_ctx.get(("M", "T", "u", "W", "r", "s"), "lex")
    M, T, u, W, r, s = ctx.gens()
    one = ctx.from_dict({(0,) * 6: 1})

    def frac(n, d=None): return (n, one if d is None else d)

    def region_subs(region, order):
        # returns dict var -> (N, D), and (nN, nD), (eN, eD) for A = n/e
        if region in ("I", "II"):
            z = frac(M * (1 + T + M * T), (1 + M) * (1 + T))
            oz = frac(M * M, (1 + M) * (1 + T))                  # M - z
            P0 = frac(1 + T + M * T, (1 + M) * T)
            if region == "I":
                P = frac(1 + T + M * T + u * (1 + M) * T, (1 + M) * T)
            else:
                P = frac((1 + T + M * T) * u, (1 + M) * T * (1 + u))
        else:
            z = frac(M * T, (1 + M) * (1 + T))
            oz = frac(M * (1 + M + M * T), (1 + M) * (1 + T))
            P = frac(u)
        other = (P[0] * W, P[1] * (1 + W))
        if order == "p>=q":
            sub = {"p": P, "q": other, "b": z, "c": oz}
        else:
            sub = {"q": P, "p": other, "c": z, "b": oz}
        sub["r"] = frac(r); sub["s"] = frac(s)
        if region == "I":
            n = (P[0] * z[1] + z[0] * P[1], P[1] * z[1]); e = (P[0] * z[0], P[1] * z[1])
            # A = n/e = (P+z)/(Pz); write n/e = (P0' z1 + z0 P1) / (P0' z0)  (common factor P1 z1 cancels)
            n = (P[0] * z[1] + z[0] * P[1], one); e = (P[0] * z[0], one)
        else:
            n = (M + 1, one); e = (M, one)
        return sub, n, e

    def substitute(poly, sub):
        names = ["p", "r", "b", "c", "q", "s"]
        pw = {}
        out = ctx.from_dict({})
        for mon, cf in poly.terms():
            t = ctx.from_dict({(0,) * 6: int(cf)})
            for x, ex, dx in zip(names, mon, deg):
                N, D = sub[x]
                key = (x, ex)
                if key not in pw:
                    pw[key] = N ** ex * D ** (dx - ex)
                t = t * pw[key]
            out += t
        return out

    allok = True
    for region in ("I", "II", "III"):
        for order in ("p>=q", "q>=p"):
            sub, n, e = region_subs(region, order)
            K = [substitute(kp[i], sub) for i in range(3)]
            nN, eN = n[0], e[0]
            mid = 2 * K[0] * eN + K[1] * nN
            end = K[0] * eN * eN + K[1] * nN * eN + K[2] * nN * nN
            for nm, pol in (("middle", mid), ("endpoint", end)):
                cfs = [int(c) for c in pol.coeffs()]
                neg = sum(1 for c in cfs if c < 0)
                ok = neg == 0 and len(cfs) > 0
                allok &= ok
                print(region, order, nm, "terms", len(cfs), "negative coefficients", neg, "OK" if ok else "FAIL", flush=True)
    # denominators used: products of (1+M),(1+T),(1+W),(1+u),T and powers -- all positive on the regions
    print("PART bernstein ALLOK", allok)
    return allok

print('==== part bernstein', flush=True)
PARTS.append(('bernstein', part_bernstein()))

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

def numeric_H(seed_arg, n_arg):
    """NUMERIC test of the (H) sufficiency claim of the T_34 document: with p=q=r=s=0, d=1,
    a = 1 + 1/(b+c), is min over SO(4) of Phi(a f2 + b f12 + c f39 + f70) >= 0 ?"""
    import numpy as np
    rng = np.random.default_rng(int(seed_arg)); n = int(n_arg)
    FV = {j: np.array(L.F_VALUES[j], float) for j in (2, 12, 39, 70)}
    worst = None
    for it in range(n):
        b, c = rng.exponential(1.0, 2) * rng.choice([1e-2, 1, 1e2])
        if rng.random() < 0.2: c = 0.0
        a = 1 + 1 / (b + c)
        vals = a * FV[2] + b * FV[12] + c * FV[39] + FV[70]
        cvec = L.solve_P([float(x) for x in vals]) if False else np.linalg.lstsq(L.P_numpy(), vals, rcond=None)[0]
        val, O = L.min_phi_over_SO4(cvec, restarts=12, seed=it)
        print("b=%.4g c=%.4g a=%.4g  min Phi = %.3e" % (b, c, a, val), flush=True)
        if worst is None or val < worst[0]: worst = (val, b, c)
    print("worst", worst)

if NUMERIC:
    numeric_stratum('2,3,5,9,10,12', '300', '325', '0.01')
    numeric_H('7', '25')

print("==== summary")
for nm, ok in PARTS: print("part", nm, "OK" if ok else "FAIL")
ok = all(ok for _, ok in PARTS)
print("EXACT CHECKS", "ALL PASS" if ok else "SOME FAILED")
sys.exit(0 if ok else 1)
