# Auditor's own code (C24-T47). Nothing from the supplement is executed or imported.
# Usage: python -I -B <this file> <path to the audit lib directory (contains moc4lib.py)> [options]
# All objects are rebuilt from the paper's definitions through the shared library moc4lib.
import sys, os, time, random, itertools, pickle
from fractions import Fraction as Fr
LIBDIR = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'lib')
sys.path.insert(0, LIBDIR)
import moc4lib as L
import flint, sympy
engine = gpoint = sys.modules[__name__]
FAILS = []
def report(name, ok, info=""):
    print(("PASS " if ok else "FAIL ") + name + ((" : " + info) if info else ""), flush=True)
    if not ok: FAILS.append(name)
def build_G(L, rays, names, degbound, dpos=5):
    # variables other than d, with Kronecker radix
    vars_ = [i for i in range(8) if i != dpos]
    bases = [degbound[names[i]] + 1 for i in vars_]
    expo = {}; e = 1
    for i, bs in zip(vars_, bases):
        expo[i] = e; e *= bs
    t = flint.fmpz_poly([0, 1])
    def mono(i): return flint.fmpz_poly([0]*expo[i] + [1])
    X = [None]*8
    for i in vars_: X[i] = mono(i)
    X[dpos] = flint.fmpz_poly([1])
    def w_of(h):
        w = []
        for r in range(14):
            s = flint.fmpz_poly([0])
            for k, ray in enumerate(rays):
                v = L.F_VECS[ray][r]; assert v.denominator == 1
                if v: s += int(v) * X[k]
            if r == 13: s += h
            w.append(s)
        return w
    C0 = L.cubic_C(w_of(flint.fmpz_poly([0])))
    Q = L.cubic_C(w_of(flint.fmpz_poly([1]))) - C0
    ms, Cs = L.cubic_sympy()
    Hs = [[sympy.Poly(sympy.diff(Cs, ms[i], ms[k]), *ms) for k in range(14)] for i in range(14)]
    def detH(h0):
        w = w_of(flint.fmpz_poly([h0]))
        M = []
        for i in range(14):
            row = []
            for k in range(14):
                s = flint.fmpz_poly([0])
                for mon, co in Hs[i][k].terms():
                    idx = [z for z, tt in enumerate(mon) if tt]
                    s += int(co) * (w[idx[0]] if idx else 1)
                row.append(s)
            M.append(row)
        n = 14; prev = flint.fmpz_poly([1]); sign = 1
        for kk in range(n-1):
            if M[kk][kk] == 0:
                sw = next((r for r in range(kk+1, n) if M[r][kk] != 0), None)
                if sw is None: return flint.fmpz_poly([0])
                M[kk], M[sw] = M[sw], M[kk]; sign = -sign
            for i in range(kk+1, n):
                for j in range(kk+1, n):
                    num = M[i][j]*M[kk][kk] - M[i][kk]*M[kk][j]
                    qq, rr = divmod(num, prev); assert rr == 0
                    M[i][j] = qq
            prev = M[kk][kk]
        return sign * M[n-1][n-1]
    # det H(g + h kappa) has degree <= 4 in h (rank Hess(kappa) = 4): interpolate from h = 0..5
    vals = [detH(h0) for h0 in range(6)]
    # Lagrange over Q[t] -> coefficients D_k (integer polys after clearing)
    hs = sympy.Symbol('h')
    # D(h) = sum D_k h^k; solve Vandermonde with Fractions
    import itertools
    V = [[Fr(h0)**k for k in range(6)] for h0 in range(6)]
    Vi = L.inv(V)
    Dk = []
    for k in range(6):
        acc = flint.fmpq_poly([0])
        for h0 in range(6):
            if Vi[k][h0]: acc += flint.fmpq_poly(vals[h0]) * flint.fmpq(Vi[k][h0].numerator, Vi[k][h0].denominator)
        Dk.append(acc)
    assert Dk[5] == 0, "degree in h exceeds 4"
    num = flint.fmpq_poly([0])
    C0q, Qq = flint.fmpq_poly(C0), flint.fmpq_poly(Q)
    for k in range(5):
        num += Dk[k] * (-C0q)**k * Qq**(4-k)
    G2Q2 = -num / 4                       # = G^2 Q^2
    G2, rem = divmod(G2Q2, Qq*Qq); assert rem == 0
    den = G2.denom(); assert den == 1, den
    G2z = flint.fmpz_poly([int(c) for c in G2.numer().coeffs()])
    Gt = G2z.sqrt()
    assert Gt is not None and Gt*Gt == G2z, "G^2 not a perfect square"
    # decode
    ctx = flint.fmpz_mpoly_ctx.get(tuple(names[i] for i in vars_), 'deglex')
    terms = {}
    for E, co in enumerate(Gt.coeffs()):
        if co == 0: continue
        ex = []; r = E
        for bs in bases: ex.append(r % bs); r //= bs
        terms[tuple(ex)] = int(co)
    G = ctx.from_dict(terms)
    return G, ctx
def check_G_points(L, gpoint, rays, G, ctx, names, rng, n=5, dpos=5):
    """compare G^2 with the independent per-point exact computation (Fraction det)."""
    Hlin, vecs = gpoint.setup(L, rays)
    ok = True
    for _ in range(n):
        cf = [Fr(rng.randint(1, 40), rng.randint(1, 40)) for _ in range(8)]; cf[dpos] = Fr(1)
        g2, _, _ = gpoint.G2(L, Hlin, vecs, cf, Fr(1))
        vals = [c for i, c in enumerate(cf) if i != dpos]
        gv = evalF(G, vals)
        ok &= (gv**2 == g2)
    return ok

def evalF(G, vals):
    tot = Fr(0)
    for ex, co in G.to_dict().items():
        m = Fr(int(co))
        for v, k in zip(vals, ex):
            if k: m *= v**int(k)
        tot += m
    return tot
def coeffs_in(F, ctx, var):
    """F as dict k -> coefficient poly (in same ctx) of var^k"""
    i = ctx.names().index(var); out = {}
    for ex, co in F.to_dict().items():
        ex = list(ex); k = ex[i]; ex[i] = 0
        out.setdefault(k, {})[tuple(ex)] = co
    return {k: ctx.from_dict(d) for k, d in out.items()}
def substitute(F, ctx, subs, octx):
    """subs: var -> (num, den) polys in octx. Clears each variable with its own degree in F:
    returns sum_alpha F_alpha prod n_i^a_i d_i^(D_i - a_i) as poly in octx (positive multiple of F
    whenever all d_i > 0)."""
    names = ctx.names()
    degs = F.degrees()
    pw = {}
    for v, (n, d) in subs.items():
        D = degs[names.index(v)]
        pw[v] = ([n**k for k in range(D+1)], [d**k for k in range(D+1)], D)
    res = octx.from_dict({})
    for ex, co in F.to_dict().items():
        t = octx.constant(int(co))
        for v, k in zip(names, ex):
            if v in pw:
                n, d, D = pw[v]; t = t * n[k] * d[D - k]
            elif k:
                raise ValueError("unsubstituted variable " + v)
        res += t
    return res
def signs(P):
    cs = [int(c) for c in P.coeffs()]
    return len(cs), sum(1 for c in cs if c > 0), sum(1 for c in cs if c < 0)
def strict_ok(P, octx, always_pos, optional, constraints):
    """P has nonnegative coefficients (checked separately). For every assignment of the optional
    variables to {0, >0} allowed by constraints (list of sets: at least one member > 0), check that
    a monomial of P is supported on positive variables."""
    names = octx.names()
    mons = [tuple(names[i] for i, k in enumerate(ex) if k) for ex in P.to_dict()]
    for pat in itertools.product([0, 1], repeat=len(optional)):
        pos = set(always_pos) | {v for v, s in zip(optional, pat) if s}
        if any(not (set(c) & pos) for c in constraints): continue
        if not any(set(m) <= pos for m in mons): return False, sorted(pos)
    return True, None
def frac(expr, syms, octx):
    """sympy rational expression -> (num, den) flint polys in octx, den with positive content"""
    import sympy
    n, d = sympy.fraction(sympy.cancel(sympy.together(expr)))
    def conv(x):
        P = sympy.Poly(sympy.expand(x), *syms)
        return octx.from_dict({m: int(c) for m, c in P.terms()})
    return conv(n), conv(d)
def substitute_deg(F, ctx, subs, octx, degs):
    """like substitute, but clearing variable v with the given degree degs[v] >= deg_v F"""
    names = ctx.names()
    pw = {v: ([n**k for k in range(degs[v]+1)], [d**k for k in range(degs[v]+1)], degs[v]) for v, (n, d) in subs.items()}
    res = octx.from_dict({})
    for ex, co in F.to_dict().items():
        t = octx.constant(int(co))
        for v, k in zip(names, ex):
            if v in pw:
                n, d, D = pw[v]; t = t * n[k] * d[D - k]
            elif k: raise ValueError(v)
        res += t
    return res
def setup(L, rays):
    ms, Cs = L.cubic_sympy()
    Hs = [[sympy.Poly(sympy.diff(Cs, ms[i], ms[k]), *ms) for k in range(14)] for i in range(14)]
    # linear form coefficients
    Hlin = [[{(mon.index(1) if 1 in mon else -1): int(co) for mon, co in Hs[i][k].terms()} for k in range(14)] for i in range(14)]
    vecs = [[Fr(x) for x in L.F_VECS[r]] for r in rays]
    return Hlin, vecs
def w_of(vecs, coef, h):
    w = [sum(c*v[i] for c, v in zip(coef, vecs)) for i in range(14)]
    w[13] += h
    return w
def G2(L, Hlin, vecs, coef, dval):
    """returns (G^2, C0, Q) with G^2 = -Q^2 det H(w(g,-C0/Q)) / (4 d^2)"""
    C0 = L.cubic_C(w_of(vecs, coef, Fr(0)))
    Q = L.cubic_C(w_of(vecs, coef, Fr(1))) - C0
    h = -C0/Q
    w = w_of(vecs, coef, h)
    H = [[sum(co*(w[i] if i >= 0 else 1) for i, co in Hlin[r][k].items()) for k in range(14)] for r in range(14)]
    D = L.det(H)
    return -Q*Q*D/(4*dval*dval), C0, Q
def basics(ctx, G):
    q,p,j,c,e,b,a = ctx.gens()
    cc = coeffs_in(G, ctx, 'c')
    assert set(cc) <= {0,1,2}
    al, be, ga = cc.get(2, ctx.from_dict({})), cc.get(1, ctx.from_dict({})), cc.get(0, ctx.from_dict({}))
    H0 = p+j+1+a; K0 = p*j + b*(p+1+a) - p*a
    H = H0**2; J = -2*H0*K0 - 4*a*(p-b)*(p+1+a); K = K0**2
    Phi = ((p-c)*j + (b-c)*(p+1+a) - p*a)**2 - 4*a*c*(p-b)*(p+1+a)
    u0 = H*ga - al*K; u1 = H*be - al*J; v0 = K*be - ga*J
    report("Phi = H c^2 + J c + K with H=(p+j+1+a)^2, K=(pj+b(p+1+a)-pa)^2", Phi == H*c**2 + J*c + K)
    report("H G = alpha Phi + u0 + c u1 (polynomial identity)", H*G == al*Phi + u0 + c*u1)
    report("K G = gamma Phi + c (v0 - c u0) (polynomial identity)", K*G == ga*Phi + c*(v0 - c*u0))
    return dict(al=al, be=be, ga=ga, H=H, J=J, K=K, u0=u0, u1=u1, v0=v0, Phi=Phi)
def normalise(ctx, G):
    al = coeffs_in(G, ctx, 'c')[2]
    return G if al(*([1]*7)) > 0 else -G
NAMES = ['q','p','j','c','e','d','b','a']
def setup_source(which):
    """build the source, check contact set / generators / paper case, C0, Q, and G"""
    if which == 'T47':
        T = {0,4,5,6,9,12}; rays = [2,9,1,45,5,35,0,7]; db = dict(q=2,p=2,j=2,c=2,e=2,b=3,a=3); nterms, tdeg = 965, 8
    else:
        T = {1,3,5,6,7,8}; rays = [27,1,59,5,2,28,6,4]; db = dict(q=3,p=2,j=2,c=2,e=2,b=3,a=3); nterms, tdeg = 1556, 9
    gen = L.generators(T)
    if which == 'T47':
        report("supplement T47 contact set {0,4,5,6,9,12} equals paper T_65 (identity element)", L.CONTACT_SETS[65] == frozenset(T))
        report("generators(T_65) are exactly the 8 rays f2,f9,f1,f45,f5,f35,f0,f7", sorted(rays) == gen, str(gen))
    else:
        report("supplement T46 contact set {1,3,5,6,7,8} equals paper T_47 (identity element)", L.CONTACT_SETS[47] == frozenset(T))
        report("generators(T_47) = the 8 source rays plus f0, f3 (absorbed in the T46 assembly via f28+f0=f1+f2+f5+f27, f28+f3=f4+f5+f27)",
               gen == sorted(rays + [0, 3]) and
               [sum(L.F_VALUES[k][s] for k in (28, 0)) for s in range(24)] == [sum(L.F_VALUES[k][s] for k in (1, 2, 5, 27)) for s in range(24)] and
               [sum(L.F_VALUES[k][s] for k in (28, 3)) for s in range(24)] == [sum(L.F_VALUES[k][s] for k in (4, 5, 27)) for s in range(24)], str(gen))
    cx = flint.fmpz_mpoly_ctx.get(tuple(NAMES) + ('h',), 'deglex'); V = cx.gens(); h = V[-1]
    w = []
    for i in range(14):
        s = cx.from_dict({})
        for v, r in zip(V[:8], rays): s += int(L.F_VECS[r][i]) * v
        if i == 13: s += h
        w.append(s)
    Cw = L.cubic_C(w); C0 = Cw.subs({'h': 0}); Q = (Cw - C0) / h
    q,p,j,c,e,d,b,a = V[:8]
    report("C(w(g,h)) is affine in h", Cw.degrees()[-1] == 1)
    if which == 'T47':
        report("restored height constant C0 = -d{ab+ad+aq+bd+bp+d^2+de+dp+dq+ep+eq+pq} (as in the T47 negative-Phi note)",
               C0 == -d*(a*b+a*d+a*q+b*d+b*p+d*d+d*e+d*p+d*q+e*p+e*q+p*q))
        report("height slope Q = ce-cj-d^2-de-dp-dq-ep-eq+jp-pq", Q == c*e-c*j-d*d-d*e-d*p-d*q-e*p-e*q+j*p-p*q)
    else:
        M = (c+d)*(a+d+e); Lq = M + b*(a+d)
        report("restored height constant C0 = -(q+1)L - jM (d=1 homogenised: -(q+d)L - jM d)", C0 == -(q+d)*Lq - j*M)
        report("height slope Q = jp - j(c+e+1) - (c+1)(e+1) (homogenised)", Q == j*p - j*(c+e+d) - (c+d)*(e+d))
    t0 = time.time()
    G, ctx = build_G(L, rays, NAMES, db)
    report(f"radial polynomial G: {len(G)} terms, total degree {G.total_degree()}", len(G) == nterms and G.total_degree() == tdeg,
           f"{round(time.time()-t0,1)} s; per-variable degrees {G.degrees()}")
    report("G^2 agrees with an independent exact Fraction 14x14 determinant at 5 random rational points",
           check_G_points(L, sys.modules[__name__], rays, G, ctx, NAMES, random.Random(5)))
    return ctx, G
ctx, G = setup_source('T47')
# PCBE failure algebra
s_, b_, e_, j_ = sympy.symbols('s b e j', positive=True)
Kc = 1/j_ + 1/(b_+1); Bc = 1/e_ + 1/(b_+1)
lhs = (1 + Bc*s_**2)/(1-s_)**2 - 1/(b_+1)        # 1/j must exceed this at failure
rhs = (b_ + 2*s_ + (b_+1)*s_**2/e_)/((b_+1)*(1-s_)**2)
report("PCBE failure: K(1-s)^2 - 1 > B s^2 (s = sqrt(c/B) < 1) is equivalent to j < (b+1)(1-s)^2/(b+2s+(b+1)s^2/e)",
       sympy.simplify(lhs - rhs) == 0, "monotonicity part PROOF-READ: K(1-x)^2-1 decreasing, Bx^2 increasing on [0,1]")

print("SUMMARY:", "all exact checks passed" if not FAILS else "FAILED: " + "; ".join(FAILS))
sys.exit(1 if FAILS else 0)
