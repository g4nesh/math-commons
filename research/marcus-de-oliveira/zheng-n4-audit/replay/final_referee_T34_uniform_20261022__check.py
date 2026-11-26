"""Independent exact replay for supplement document final_referee_T34_uniform_20261022.
Own code only: rebuilds generators, profile identities, the radial cubic restriction C(c + h kappa)
and the full 14x14 Hessian of the paper's cubic (A.2) from moc4lib, and replays the positivity tables.
The degree-nine polynomial G is loaded as DATA from the supplement's JSON (never executed) and is
validated against the Hessian computed here by an exact randomized identity test.
Usage: MOC4_SUPP=<supplement root> python -I -B final_referee_T34_uniform_20261022__check.py   (or pass the root as argv[1])
"""
# helper block (to be embedded in replay scripts)
import os, sys, json, random, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
LIBDIR = os.environ.get("MOC4_LIB") or os.path.join(HERE, "..", "lib")
sys.path.insert(0, LIBDIR)
import moc4lib as L
import flint
import sympy

def cost_vec(gens, coeffs):
    """c = sum coeffs[k] * f_{gens[k]} in V_R (14 coords)."""
    c = [0]*14
    for j, x in zip(gens, coeffs):
        for r in range(14):
            c[r] = c[r] + x*L.F_VECS[j][r]
    return c

MS, CSYM = L.cubic_sympy()
T3 = [[[int(sympy.diff(CSYM, MS[i], MS[k], MS[r])) for r in range(14)] for k in range(14)] for i in range(14)]
def hess_at(m):
    return [[sum((T3[i][k][r]*m[r] for r in range(14) if T3[i][k][r]), Fr(0)) for k in range(14)] for i in range(14)]

def radial_CQ(gens, names):
    """Symbolic C(c + h kappa) = Cst + h Q, as flint polys in names (+h check)."""
    syms = sympy.symbols(names)
    h = sympy.Symbol('h')
    c = cost_vec(gens, syms)
    m = [c[r] + (h if r == 13 else 0) for r in range(14)]
    e = sympy.Poly(sympy.expand(L.cubic_C(m)), h, *syms)
    assert e.degree(h) <= 1, "C(c+h kappa) not linear in h"
    ctx = flint.fmpz_mpoly_ctx.get(tuple(names), 'lex')
    Cst, Q = {}, {}
    for mon, cf in e.terms():
        (Cst if mon[0] == 0 else Q)[tuple(mon[1:])] = int(cf)
    return ctx, ctx.from_dict(Cst), ctx.from_dict(Q)

def load_json_poly(ctx, terms, nv):
    return ctx.from_dict({tuple(t[0][:nv]): int(Fr(t[1])) for t in terms})

def evalp(p, vals):
    tot = Fr(0)
    for mon, cf in p.terms():
        t = Fr(int(cf))
        for v, e in zip(vals, mon):
            if e: t *= Fr(v)**int(e)
        tot += t
    return tot

def hessian_identity_test(gens, ctx, Cst, Q, G, npts=12, seed=1, R=10**6):
    """Randomized exact test of Q^2 det Hess C(c - (C*/Q) kappa) == -4 G^2 (Schwartz-Zippel)."""
    rng = random.Random(seed)
    ok = 0
    for _ in range(npts):
        while True:
            vals = [rng.randint(-R, R) for _ in range(ctx.nvars())]
            qv = evalp(Q, vals)
            if qv != 0: break
        cv = evalp(Cst, vals)
        h = -cv/qv
        c = cost_vec(gens, [Fr(v) for v in vals])
        m = [c[r] + (h if r == 13 else 0) for r in range(14)]
        D = L.det(hess_at(m))
        gv = evalp(G, vals)
        if qv*qv*D == -4*gv*gv: ok += 1
    return ok

def sym_to_flint(expr, tctx, tsyms):
    pe = sympy.Poly(sympy.expand(expr), *tsyms)
    return tctx.from_dict({tuple(m): int(c) for m, c in pe.terms()})

def clear_subst(P, pnames, subs, tnames, extra_den=None, degs=None):
    """P: fmpz_mpoly over pnames. subs: dict name -> sympy rational expr in tnames
    (names not in subs are mapped to the same-named target variable).
    Returns the cleared polynomial sum coef prod N_v^e D_v^(m_v - e) over target ctx."""
    tsyms = sympy.symbols(tnames)
    loc = dict(zip(tnames, tsyms))
    tctx = flint.fmpz_mpoly_ctx.get(tuple(tnames), 'lex')
    degs = P.degrees() if degs is None else degs
    hnames = tuple(pnames) + tuple('w_' + n for n in pnames)
    hctx = flint.fmpz_mpoly_ctx.get(hnames, 'lex')
    nv = len(pnames)
    H = hctx.from_dict({tuple(m) + tuple(int(degs[i]) - int(m[i]) for i in range(nv)): int(c) for m, c in P.terms()})
    Ns, Ds = [], []
    for n in pnames:
        ex = subs.get(n, None)
        if ex is None:
            ex = loc[n] if n in loc else sympy.Integer(0)
        ex = sympy.sympify(ex, locals=loc)
        num, den = sympy.fraction(sympy.together(ex))
        Ns.append(sym_to_flint(num, tctx, tsyms)); Ds.append(sym_to_flint(den, tctx, tsyms))
    R = H.compose(*(Ns + Ds), ctx=tctx)
    if extra_den is not None:
        R = R * sym_to_flint(sympy.sympify(extra_den, locals=loc), tctx, tsyms)
    return R, tctx

def negcount(R):
    cs = R.coeffs()
    return sum(1 for c in cs if c < 0), len(cs)

SUPP = os.environ.get("MOC4_SUPP") or (sys.argv[1] if len(sys.argv) > 1 else None)
if not SUPP:
    sys.exit("set MOC4_SUPP to the supplement root")
FAIL = []
def check(name, cond, info=""):
    print(("PASS " if cond else "FAIL ") + "EXACT " + name + ((" | " + str(info)) if info != "" else ""))
    if not cond: FAIL.append(name)
def numeric(name, info):
    print("NUMERIC " + name + " | " + str(info))
def finish():
    print("SUMMARY exact failures:", len(FAIL), FAIL)
    sys.exit(1 if FAIL else 0)
import numpy as np
DOC = 'final_referee_T34_uniform_20261022'
T = L.CONTACT_SETS[35]
check("doc contact set {2,3,4,6,9,12} equals paper T_35 literally", sorted(T) == [2,3,4,6,9,12])
gens = L.generators(T)
check("generators(T_35) == source 0,1,2,5,7,9,53,67,68", gens == [0,1,2,5,7,9,53,67,68], gens)
check("only nonminor generator is 68", [j for j in gens if not L.F_IS_MINOR[j]] == [68])
check("rank of source = 8", L.rank([L.F_VECS[j] for j in gens]) == 8)
V = L.F_VALUES
check("profile identity f68 + f2 = f0 + f9 + f67", all(V[68][s]+V[2][s] == V[0][s]+V[9][s]+V[67][s] for s in range(24)))
G1 = [67,5,53,9,1,68,7,0]; names = ['q','p','j','c','e','d','b','a']
check("representation (1) basis {67,5,53,9,1,68,7,0} independent (coefficients unique)", L.rank([L.F_VECS[j] for j in G1]) == 8)
ctx, Cst, Q = radial_CQ(G1, names)
dd = json.load(open(os.path.join(SUPP, 'proofs/manin_transfer_20261021/final_referee_T34_full_radial_20261022.json')))
check("C(c + h kappa) is linear in h; own C*, Q equal the JSON C_terms, Q_terms", load_json_poly(ctx, dd['C_terms'], 8) == Cst and load_json_poly(ctx, dd['Q_terms'], 8) == Q)
q,p,j,c,e,d,b,a = ctx.gens()
check("Q|d=1 == c + ce + cj - cp + 1 + e + j + ej", Q.subs({'d':1}) == c + c*e + c*j - c*p + 1 + e + j + e*j)
check("C* has only positive coefficients and contains d^3", all(x > 0 for x in Cst.coeffs()) and any(m == (0,0,0,0,0,3,0,0) for m, x in Cst.terms()))
G = load_json_poly(ctx, dd['radial_factors'][1]['terms'], 8)
check("JSON radial factor 0 is +-Q; G has 1475 terms, degree 9, d^9 coefficient +1", load_json_poly(ctx, dd['radial_factors'][0]['terms'], 8) in (Q, -Q) and len(G) == 1475 and G.total_degree() == 9 and [int(x) for m, x in G.terms() if m == (0,0,0,0,0,9,0,0)] == [1])
k_ok = hessian_identity_test(G1, ctx, Cst, Q, G, npts=12, seed=1)
check("(2) Q^2 det Hess C(c - (C*/Q) kappa) == -4 G^2 at 12 random integer points in [-1e6,1e6]^8 (exact arithmetic; Schwartz-Zippel)", k_ok == 12, k_ok)
Gd = G.subs({'d':1})
check("(7) G00 = G|p=e=j=0 == (1+b)(c+1+b)^2(a+1)^3(q+1)^3", Gd.subs({'p':0,'e':0,'j':0}) == (1+b)*(c+1+b)**2*(a+1)**3*(q+1)**3)
check("k0|b=0 == (a+1)(q+1)(e+1+a)^2(c+1)^2(q+j+1)^2", Gd.subs({'p':0,'b':0}) == (a+1)*(q+1)*(e+1+a)**2*(c+1)**2*(q+j+1)**2)
check("(9) baseline G|c=p=e=j=0 == (1+b)^3(a+1)^3(q+1)^3", Gd.subs({'c':0,'p':0,'e':0,'j':0}) == (1+b)**3*(a+1)**3*(q+1)**3)
check("G|c=b=0 == (qp+q+qa+p+pa+1+a)(e+1+a)^2(q+j+1)^2", Gd.subs({'c':0,'b':0}) == (q*p+q+q*a+p+p*a+1+a)*(e+1+a)**2*(q+j+1)**2)
pn = ['q','p','j','c','e','b','a']; c7 = flint.fmpz_mpoly_ctx.get(tuple(pn), 'lex')
G7 = c7.from_dict({tuple(m[:5]) + tuple(m[6:]): int(cf) for m, cf in Gd.terms()})
check("G has degree 2 in p", G7.degrees()[1] == 2)
k = [c7.from_dict({tuple(m[:1])+(0,)+tuple(m[2:]): int(cf) for m, cf in G7.terms() if m[1] == i}) for i in range(3)]
def charts(bexpr):
    lo_a = {'a': f'X/(({bexpr})*(1+X))', 'e': 'E'}
    hi_a = {'a': f'(1+X)/({bexpr})', 'e': f'((1+X)/({bexpr})+1)*(({bexpr})+1)/X*E/(1+E)'}
    lo_q = {'q': f'Y/(({bexpr})*(1+Y))', 'j': 'J'}
    hi_q = {'q': f'(1+Y)/({bexpr})', 'j': f'((1+Y)/({bexpr})+1)*(({bexpr})+1)/Y*J/(1+J)'}
    out = []
    for na, A in (('ab<1', lo_a), ('ab>1', hi_a)):
        for nq, Qd in (('bq<1', lo_q), ('bq>1', hi_q)):
            s = {}; s.update(A); s.update(Qd); s['b'] = bexpr; s['c'] = 'C'
            out.append((na + ',' + nq, s))
    return out
exp_k0 = [2227, 3268, 3268, 4945]; exp_c0 = [1419, 2110, 2110, 3211]
for (nm, s), ex in zip(charts('B'), exp_k0):
    R, _ = clear_subst(k[0], pn, s, ['X','Y','E','J','B','C'])
    nc = negcount(R); check(f"(8) k0 chart {nm}: cleared polynomial has no negative coefficient (doc count {ex})", nc[0] == 0, nc)
for (nm, s), ex in zip(charts('B'), exp_c0):
    s = dict(s); s['c'] = '0'; s['p'] = 'P'
    R, _ = clear_subst(G7, pn, s, ['X','Y','E','J','B','P'])
    nc = negcount(R); check(f"(9) G|c=0 chart {nm}: no negative coefficient (doc count {ex})", nc[0] == 0, nc)
cC = c7.gens()[3]
Mpol = 2*cC*k[0] + (cC+1)*k[1]; Epol = cC**2*k[0] + (cC+1)*cC*k[1] + (cC+1)**2*k[2]
exp_b = [3612, 7230, 7230, 12258]
for (nm, s), ex in zip(charts('B/(1+B)'), exp_b):
    for lab, Pp in (('middle', Mpol), ('endpoint', Epol)):
        R, _ = clear_subst(Pp, pn, s, ['X','Y','E','J','B','C'])
        nc = negcount(R); check(f"b<=1 Bernstein {lab} numerator, P*=1+1/c, chart {nm}: no negative coefficient (doc count {ex})", nc[0] == 0, nc)
# allocation identities
dl, aS, bS, ES = sympy.symbols('delta a b E', positive=True)
Ea = dl*(bS+1)*(aS+dl)/(aS*bS-dl)
check("(13) a(delta b)(E_a(delta)-delta) - delta^2(a + delta b + E_a(delta) + delta) == 0", sympy.simplify(aS*dl*bS*(Ea-dl) - dl**2*(aS+dl*bS+Ea+dl)) == 0)
check("(11) derivative of E_a == (b+1)(a^2 b + 2ab delta - delta^2)/(ab-delta)^2", sympy.simplify(sympy.diff(Ea, dl) - (bS+1)*(aS**2*bS+2*aS*bS*dl-dl**2)/(aS*bS-dl)**2) == 0)
check("E_a(1) == (a+1)(b+1)/(ab-1) (the four-coefficient threshold)", sympy.simplify(Ea.subs(dl, 1) - (aS+1)*(bS+1)/(aS*bS-1)) == 0)
B_, U_, V_, X_, Y_, C_ = sympy.symbols('B U V X Y C')
Dd = 1+U_+V_; bb = 1+B_; d1 = U_/Dd; d2 = V_/Dd; rr = 1/Dd
aa = (U_+X_*Dd)/(bb*Dd); qq = (V_+Y_*Dd)/(bb*Dd)
check("(15) chart: X = ab - delta1, Y = bq - delta2", sympy.simplify(aa*bb - d1 - X_) == 0 and sympy.simplify(qq*bb - d2 - Y_) == 0)
e15 = U_*(bb+1)*(X_*Dd+U_*(1+bb))/(bb*Dd**2*X_)
check("(15) e == E_a(delta1) in the chart", sympy.simplify(Ea.subs({dl: d1, aS: aa, bS: bb}) - e15) == 0)
check("(16) r^2/c + r/b == (b + C D)/(C b D^2)", sympy.simplify(rr**2/C_ + rr/bb - (bb + C_*Dd)/(C_*bb*Dd**2)) == 0)
s = {'b': '1+B', 'a': '(U+X*(1+U+V))/((1+B)*(1+U+V))', 'q': '(V+Y*(1+U+V))/((1+B)*(1+U+V))', 'c': 'C',
     'e': 'U*(2+B)*(X*(1+U+V)+U*(2+B))/((1+B)*(1+U+V)**2*X)', 'j': 'V*(2+B)*(Y*(1+U+V)+V*(2+B))/((1+B)*(1+U+V)**2*Y)'}
tn = ['B','U','V','X','Y','C']; tsy = sympy.symbols(tn); B_,U_,V_,X_,Y_,C_ = tsy
Dd = 1+U_+V_
tctx = flint.fmpz_mpoly_ctx.get(tuple(tn), 'lex')
Nf = sym_to_flint((1+B_) + C_*Dd, tctx, tsy); Df = sym_to_flint(sympy.expand(C_*(1+B_)*Dd**2), tctx, tsy)
deg = G7.degrees(); pad = c7.from_dict({tuple(int(deg[t]) if t != 1 else 0 for t in range(7)): 1})
Rpad, _ = clear_subst(pad, pn, s, tn)
ks = []
for i in range(3):
    R, _ = clear_subst(k[i] + pad, pn, s, tn); ks.append(R - Rpad)
for lab, R, ex in (('base', ks[0], 117110), ('middle', 2*Df*ks[0] + Nf*ks[1], 209489), ('endpoint', Df**2*ks[0] + Nf*Df*ks[1] + Nf**2*ks[2], 289848)):
    nc = negcount(R); check(f"b>=1 allocation chart (15)-(16) Bernstein {lab}: no negative coefficient (doc count {ex})", nc[0] == 0, nc)

if os.environ.get("MOC4_NUMERIC", "1") == "1":
    import numpy as np
    nrng = np.random.default_rng(2026)
    def piece(dct):
        cc = [0.0]*14
        for jx, x in dct.items():
            for rr_ in range(14): cc[rr_] += float(x)*float(L.F_VECS[jx][rr_])
        return cc
    def worst(cases, restarts=6):
        wv = 1e9
        for dct in cases:
            v, _ = L.min_phi_over_SO4(piece(dct), restarts=restarts, seed=int(nrng.integers(1e6)))
            wv = min(wv, v/sum(abs(x) for x in dct.values()))
        return wv
    RR = lambda lo=0.05, hi=5: float(np.exp(nrng.uniform(np.log(lo), np.log(hi))))
    def loadGnum(fn, nms):
        cx = flint.fmpz_mpoly_ctx.get(tuple(nms), 'lex')
        dj = json.load(open(os.path.join(SUPP, 'proofs/manin_transfer_20261021/' + fn)))
        Gx = load_json_poly(cx, dj['radial_factors'][1]['terms'], 8)
        E_ = np.array([np.array(m, float) for m, cf in Gx.terms()]); C_ = np.array([float(cf) for m, cf in Gx.terms()])
        sg = 1.0 if [cf for m, cf in Gx.terms() if m[nms.index('d')] == 9][0] > 0 else -1.0
        return lambda x: sg*float(C_ @ np.prod(np.asarray(x, float)[None, :]**E_, axis=1))
    def main_stmt(nu):
        wv = 1e9; found = 0
        for sd in range(3):
            pt = L.random_point_on_contact_stratum(sorted(L.CONTACT_SETS[nu]), seed=sd, tries=20)
            if pt is None: continue
            found += 1
            m = L.matching_vector_np(pt['a'], pt['b'], pt['alpha']); Fv = L.P_numpy() @ m
            cc = np.real(m).copy(); cc[13] -= np.real(Fv[pt['pi_star']])
            v, _ = L.min_phi_over_SO4(cc, restarts=10, seed=sd)
            wv = min(wv, v/np.abs(L.P_numpy() @ cc).sum())
        numeric("main statement T_%d: random_point_on_contact_stratum seeds 0-2, min_phi_over_SO4 restarts 10" % nu, "points %d/3, worst normalized min Phi %.2e" % (found, wv))
    n = 6
    numeric("(3) old-H at equality c(p-1)=1, piece f68 + c f9 + p f5, 6 samples x 6 restarts", "worst min Phi/l1 %.2e" % worst([{68:1, 9:cc, 5:1+1/cc} for cc in [RR() for _ in range(n)]]))
    numeric("(4) chain at equality p=1/b+1/c, piece f68 + p f5 + b f7 + c f9", "worst %.2e" % worst([{68:1, 5:1/bb+1/cc, 7:bb, 9:cc} for bb, cc in [(RR(), RR()) for _ in range(n)]]))
    cs = []
    for _ in range(n):
        aa = RR(); bb = (1+RR())/aa; cs.append({68:1, 7:bb, 0:aa, 1:(aa+1)*(bb+1)/(aa*bb-1)})
    numeric("(5) four-coefficient at equality, piece f68 + b f7 + a f0 + e f1", "worst %.2e" % worst(cs))
    cs = []
    for _ in range(n):
        qq = RR(); bb = (1+RR())/qq; cs.append({68:1, 7:bb, 67:qq, 53:(qq+1)*(bb+1)/(qq*bb-1)})
    numeric("(6) four-coefficient at equality, piece f68 + b f7 + q f67 + j f53", "worst %.2e" % worst(cs))
    G34 = loadGnum('final_referee_T34_full_radial_20261022.json', ['q','p','j','c','e','d','b','a'])
    from scipy.optimize import brentq
    def Ea_(x, a_, b_): return x*(b_+1)*(a_+x)/(a_*b_-x)
    wv = 1e9; cnt = 0; N = 20000
    while cnt < N:
        q_, j_, c_, e_, b_, a_ = [RR(0.02, 20) for _ in range(6)]
        if a_*b_*(e_-1) >= a_+b_+e_+1 or b_*q_*(j_-1) >= b_+q_+j_+1: continue
        if b_ <= 1: Ps = 1 + 1/c_
        else:
            d1 = brentq(lambda x: Ea_(x, a_, b_) - e_, 0, min(1, a_*b_)*(1-1e-15))
            d2 = brentq(lambda x: Ea_(x, q_, b_) - j_, 0, min(1, q_*b_)*(1-1e-15))
            r_ = 1 - d1 - d2
            if r_ <= 0: continue
            Ps = r_*r_/c_ + r_/b_
        cnt += 1
        wv = min(wv, G34([q_, Ps*nrng.uniform(0, 1), j_, c_, e_, 1, b_, a_]))
    numeric("G > 0 below the selected threshold P* (old-H for b<=1, allocation (14) for b>=1) on %d random hard points (log-uniform in [0.02,20], seed 2026)" % N, "min G = %.3e" % wv)
    main_stmt(35)

finish()
