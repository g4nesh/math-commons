"""Independent exact replay for supplement document independent_moc4_T36_uniform_20261022.
Own code only: rebuilds generators, profile identities, the radial cubic restriction C(c + h kappa)
and the full 14x14 Hessian of the paper's cubic (A.2) from moc4lib, and replays the positivity tables.
The degree-nine polynomial G is loaded as DATA from the supplement's JSON (never executed) and is
validated against the Hessian computed here by an exact randomized identity test.
Usage: MOC4_SUPP=<supplement root> python -I -B independent_moc4_T36_uniform_20261022__check.py   (or pass the root as argv[1])
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
T = L.CONTACT_SETS[37]
check("doc contact set {0,5,6,7,9,12} equals paper T_37 literally", sorted(T) == [0,5,6,7,9,12])
gens = L.generators(T)
check("generators(T_37) == {0,1,2,5,7,29,35,39,40,45}, nonminor 35, 40, rank 8", gens == [0,1,2,5,7,29,35,39,40,45] and [j for j in gens if not L.F_IS_MINOR[j]] == [35,40] and L.rank([L.F_VECS[j] for j in gens]) == 8, gens)
NS = L.nullspace([L.P[s] for s in sorted(T)])
check("contact nullspace has dimension 8", len(NS) == 8)
ok = {}
for c in NS:
    g = L.values(c); d35 = g[3]-g[1]; d40 = g[14]-g[16]
    ok.setdefault('delta35 + delta40 == -g17', []).append(d35 + d40 == -g[17])
    reps = {'(1)': {0: g[18], 1: -d40, 2: g[13], 5: g[11], 7: g[14], 29: g[2], 39: -d35, 45: g[3]},
            '(2)': {2: g[15], 29: g[4], 1: g[17], 45: g[1], 5: g[11], 35: d35, 0: g[18], 7: g[14]},
            '(3)': {0: g[23], 39: g[17], 5: g[11], 2: g[13], 29: g[2], 7: g[16], 40: d40, 45: g[3]}}
    for nm, rep in reps.items():
        ok.setdefault('identity ' + nm + ' g = sum coef f_j as linear operators on the contact nullspace', []).append([sum(cf*L.F_VALUES[j][s] for j, cf in rep.items()) for s in range(24)] == list(g))
for kx, v in ok.items(): check(kx, all(v))
check("(1) uses only minor generators 0,1,2,5,7,29,39,45", all(L.F_IS_MINOR[j] for j in (0,1,2,5,7,29,39,45)))
V = L.F_VALUES
def ident(lhs, rhs): return all(sum(V[j][s] for j in lhs) == sum(V[j][s] for j in rhs) for s in range(24))
check("f35 + f40 = f0 + f7 + f2 + f29 + f45", ident([35,40],[0,7,2,29,45]))
check("f35 + f39 = f1 + f2 + f29 + f45", ident([35,39],[1,2,29,45]))
check("f40 + f1 = f0 + f7 + f39", ident([40,1],[0,7,39]))
GG = {}
for tag, B, names, fn in (('35', [2,29,1,45,5,35,0,7], ['q','p','j','c','e','d','b','a'], 'independent_moc4_T36_bad35_full_radial_20261022.json'),
                          ('40', [0,39,5,2,29,7,40,45], ['p','r','a','b','c','q','d','s'], 'independent_moc4_T36_bad40_full_radial_20261022.json')):
    ctx, Cst, Q = radial_CQ(B, names)
    dd = json.load(open(os.path.join(SUPP, 'proofs/manin_transfer_20261021/' + fn)))
    check(f"({tag}) own C, Q equal JSON C_terms, Q_terms", load_json_poly(ctx, dd['C_terms'], 8) == Cst and load_json_poly(ctx, dd['Q_terms'], 8) == Q)
    G = -load_json_poly(ctx, dd['radial_factors'][1]['terms'], 8)
    di = names.index('d')
    check(f"({tag}) radial factor0 is +-Q; minus stored G has d^9 coefficient +1", load_json_poly(ctx, dd['radial_factors'][0]['terms'], 8) in (Q, -Q) and [int(x) for m, x in G.terms() if m[di] == 9] == [1], len(G))
    kk = hessian_identity_test(B, ctx, Cst, Q, G, npts=10, seed=11)
    check(f"({tag}) (4) Q^2 det Hess C(c - (C/Q)kappa) == -4 G^2 at 10 random integer points (exact; Schwartz-Zippel)", kk == 10, kk)
    GG[tag] = (ctx, G, Cst)
# branch 35
ctx, G, Cst = GG['35']
q,p,j,c,e,d,b,a = ctx.gens()
check("(35) all coefficients of C negative, including -d^3", all(x < 0 for x in Cst.coeffs()) and any(m[5] == 3 for m, x in Cst.terms()))
check("(35) exchange symmetry G(a,b,q,p) == G(b,a,p,q)", G.compose(p, q, j, c, e, d, a, b) == G)
pn = ['q','p','j','c','e','b','a']; c7 = flint.fmpz_mpoly_ctx.get(tuple(pn), 'lex')
Gd = G.subs({'d':1}); G7 = c7.from_dict({tuple(m[:5]) + tuple(m[6:]): int(cf) for m, cf in Gd.terms()})
check("(35) G degree 2 in c", G7.degrees()[3] == 2)
q,p,j,c,e,b,a = c7.gens()
k = [c7.from_dict({tuple(m[:3]) + (0,) + tuple(m[4:]): int(cf) for m, cf in G7.terms() if m[3] == i}) for i in range(3)]
Cb = (a+1)*(b+1) + a*q + b*p + (1+e)*(p*q+p+q) + e
R90, rem = divmod(k[0], Cb)
check("(35) k0 == Cbar * R90 exactly", rem == 0, len(R90))
check("(35) G|e=j=0 == [(1+a+b)c + (a+1)(b+1)](p+a+1)^2(q+b+1)^2", G7.subs({'e':0,'j':0}) == ((1+a+b)*c + (a+1)*(b+1))*(p+a+1)**2*(q+b+1)**2)
K = (a*b-1)*e - (a+1)*(b+1)
Diff = R90 - (a+1)*(b+1)*(-K)
Emax = '((A0+1)*(B0+1)/(A0*B0-1))'
for nm, (aa, bb) in (('a=zR,b=z/R', ('z*R','z/R')), ('exchanged', ('z/R','z*R'))):
    for zc, zz, ee, ex in (('ab<1', 'T/(1+T)', 'V', 280), ('ab>1', '1+T', None, 404)):
        s = {'a': aa.replace('z', f'({zz})').replace('R','(1+A)'), 'b': bb.replace('z', f'({zz})').replace('R','(1+A)')}
        if ee is None: ee = Emax.replace('A0', f"({s['a']})").replace('B0', f"({s['b']})") + '*V/(1+V)'
        s['e'] = ee
        s = {kx: str(sympy.cancel(sympy.sympify(v))) for kx, v in s.items()}
        R, _ = clear_subst(Diff, pn, s, ['A','T','V','q','p','j'])
        nc = negcount(R); check(f"(8) R90 - (a+1)(b+1)(-K), chart {nm}, {zc}: no negative coefficient (doc count {ex})", nc[0] == 0, nc)
A,T,V,U = sympy.symbols('A T V U')
b0 = T/(1+T); b1 = 1+T
def EM(a_, b_): return (a_+1)*(b_+1)/(a_*b_-1)
aI = A; bIm = T/(A*(1+T)); bIp = (1+T)/A; eIp = EM(aI, bIp)*V/(1+V); aIIIa = b0*V/(1+V); w = V/(1+V)
rows = [('I-', aI, bIm, V, V+U, None, (969,1704,2364)),
        ('I+', aI, bIp, eIp, eIp+U, None, (2643,5745,10406)),
        ('II', aI, bIp, eIp, eIp*U/(1+U), 1/eIp - 1/EM(aI, bIp), (4414,9616,16507)),
        ('IIIa-', aIIIa, b0, A, A*U/((1+A)*(1+U)), 1 + 1/A, (3779,4535,4535)),
        ('IIIa+', aIIIa, b0, A, A/(1+A) + A**2*U/((1+A)*(1+U)), None, (4284,7308,9756)),
        ('IIIb+', V/(b1*(1+V)), b1, A, A*b1/(A+b1) + A**2*U/((A+b1)*(1+U)), None, (4428,8448,12396)),
        ('IIIb-', V**2/(b1*(1+V)**2), b1, A, A*b1*U/((A+b1)*(1+U)), 1/A + b1*(1-w)**2/(b1-w)**2, (7028,13410,18874))]
tn = ['A','T','V','U','q','p']; tsy = sympy.symbols(tn); tctx = flint.fmpz_mpoly_ctx.get(tuple(tn), 'lex')
deg = G7.degrees(); pad = c7.from_dict({tuple(int(deg[t]) if t != 3 else 0 for t in range(7)): 1})
for nm, aa, bb, ee, jj, Jt, ex in rows:
    if Jt is None: Jt = 1/jj
    Jt = sympy.cancel(sympy.together(Jt)); N, Dn = sympy.fraction(Jt)
    s = {kx: str(sympy.cancel(sympy.together(v))) for kx, v in dict(a=aa, b=bb, e=ee, j=jj).items()}
    Rpad, _ = clear_subst(pad, pn, s, tn); ks = []
    for i in range(3):
        R, _ = clear_subst(k[i] + pad, pn, s, tn); ks.append(R - Rpad)
    Nf = sym_to_flint(sympy.expand(N), tctx, tsy); Df = sym_to_flint(sympy.expand(Dn), tctx, tsy)
    posND = all(x > 0 for x in Nf.coeffs()) and all(x > 0 for x in Df.coeffs())
    res = [negcount(ks[0]), negcount(2*Df*ks[0] + Nf*ks[1]), negcount(Df**2*ks[0] + Nf*Df*ks[1] + Nf**2*ks[2])]
    check(f"seven intervals row {nm}: J = N/D with positive N, D; base/middle/endpoint Bernstein numerators have no negative coefficient (doc counts {ex})", posND and all(r[0] == 0 for r in res), res)
# branch 40
ctx, G, Cst = GG['40']
p,r,a,b,c,q,d,s = ctx.gens()
M = b + c; U_ = d + r; T0 = b*(d+q) + c*(d+p) + (d+p)*(d+q)
check("(5) C == bc(U - a) + U T0", Cst == b*c*(U_-a) + U_*T0)
last = M*U_*T0 - d**2*b*c
check("(5) M C == bc(dM + d^2 - aM) + rMbc + (MUT0 - d^2 bc), last parenthesis nonnegative coefficients with b d^3, c d^3", M*Cst == b*c*(d*M + d**2 - a*M) + r*M*b*c + last and all(x >= 0 for x in last.coeffs()))
check("(5) C|b=c=0 == (d+r)(d+p)(d+q)", Cst.subs({'b':0,'c':0}) == (d+r)*(d+p)*(d+q))
pn6 = ['p','r','a','b','c','q']; c6 = flint.fmpz_mpoly_ctx.get(tuple(pn6), 'lex')
Gd = G.subs({'d':1})
check("(40) G degree 2 in s", Gd.degrees()[7] == 2)
Ls = [c6.from_dict({tuple(m[:6]): int(cf) for m, cf in Gd.terms() if m[7] == i}) for i in range(3)]
p,r,a,b,c,q = c6.gens()
L0, L1, L2 = Ls
B1 = L1 + 2*a*L0; B2 = L2 + a*L1 + a**2*L0
L00 = L0.subs({'a':0})
check("(12) L0|a=0 has only positive coefficients and constant term 1", all(x > 0 for x in L00.coeffs()) and L00.subs({'p':0,'r':0,'b':0,'c':0,'q':0}) == 1)
Ssym, A6, s6 = sympy.symbols('S a s')
l0, l1, l2 = sympy.symbols('l0 l1 l2')
check("(13) (1+aS)^2 G(S/(1+aS)) == L0 + B1 S + B2 S^2 (symbolic in L0, L1, L2)", sympy.expand(sympy.cancel((1+A6*Ssym)**2*(l0 + l1*Ssym/(1+A6*Ssym) + l2*(Ssym/(1+A6*Ssym))**2)) - (l0 + (l1 + 2*A6*l0)*Ssym + (l2 + A6*l1 + A6**2*l0)*Ssym**2)) == 0)
Tt, Lv, Uu, Ww, Yy = sympy.symbols('T L U W Y')
th = Tt/(1+Tt)
def energy(bb, cc, pp, qq, d1, d2, be1):
    be2 = 1 - be1
    return 1 - d1 - d2 + be1**2/bb + be2**2/cc + d1**2/pp + d2**2/qq
charts = {}
pHH = 2*th*Uu/(1+Uu); qHH = 2*(1-th)*Ww/(1+Ww)
charts['HH'] = (dict(b=th/Lv, c=(1-th)/Lv, p=pHH, q=qHH), 1 + Lv - (pHH+qHH)/4, (pHH/2, qHH/2))
charts['HC'] = (dict(b=th/Lv, p=pHH, q=2*(1-th)*(1+Ww), c=2*(1-th)*(1+Ww)/(2*Lv*(1+Ww)+Ww)), Lv + sympy.Rational(1,2) + th/(2*(1+Uu)), (pHH/2, 1-th))
qCH = 2*(1-th)*Uu/(1+Uu)
charts['CH'] = (dict(c=(1-th)/Lv, q=qCH, p=2*th*(1+Ww), b=2*th*(1+Ww)/(2*Lv*(1+Ww)+Ww)), Lv + sympy.Rational(1,2) + (1-th)/(2*(1+Uu)), (th, qCH/2))
Vv = Uu + Ww
charts['CC'] = (dict(p=2*th*(1+Uu), q=2*(1-th)*(1+Vv), b=th/Lv, c=2*(1-th)*(1+Uu)*(1+Vv)/(2*Lv*(1+Uu)*(1+Vv)+Ww)), Lv + 1/(2*(1+Uu)), (th, 1-th))
charts["CC'"] = (dict(q=2*(1-th)*(1+Uu), p=2*th*(1+Vv), c=(1-th)/Lv, b=2*th*(1+Uu)*(1+Vv)/(2*Lv*(1+Uu)*(1+Vv)+Ww)), Lv + 1/(2*(1+Uu)), (th, 1-th))
expc = {'HH': ((6372,7776,7757),(9408,11466,13524,13502),(8832,15960,18979,21801,21863)),
        'HC': ((6084,7380,7379),(7728,9366,11004,11003),(5456,11550,13581,15591,15599)),
        'CH': ((6084,7380,7379),(7728,9366,11004,11003),(5456,11550,13581,15591,15599)),
        'CC': ((6645,8394,8259),(7029,8958,10887,10887),(4252,9048,11157,13266,13266)),
        "CC'": ((6645,8394,8259),(7029,8958,10887,10887),(4252,9048,11157,13266,13266))}
tn = ['T','L','U','W','Y','r']
deg = Gd.degrees()
for nm, (sb, J, (d1, d2)) in charts.items():
    en = sympy.simplify(energy(sb['b'], sb['c'], sb['p'], sb['q'], d1, d2, th) - J)
    feas = [sympy.simplify(th - d1), sympy.simplify(1 - th - d2)]
    check(f"chart {nm}: J equals the energy (10) of the feasible allocation beta=(theta,1-theta), delta as in the chart, 0<=delta_i<=beta_i", en == 0 and all(sympy.Poly(sympy.numer(sympy.together(f)), Tt, Uu, Ww).coeffs() == [] or all(x >= 0 for x in sympy.Poly(sympy.numer(sympy.together(f)), Tt, Uu, Ww).coeffs()) for f in feas), feas)
    Jc = sympy.cancel(sympy.together(J)); nN, dD = sympy.fraction(Jc)
    sub = {kx: str(sympy.cancel(sympy.together(v))) for kx, v in sb.items()}
    sub['a'] = str(sympy.cancel(Jc*Yy/(1+Yy)))
    for (lab, Pp, m), ex in zip((('L0', L0, 2), ('B1', B1, 3), ('B2', B2, 4)), expc[nm]):
        dg = (int(deg[0]), int(deg[1]), m, int(deg[3]), int(deg[4]), int(deg[5]))
        R, _ = clear_subst(Pp, pn6, sub, tn, degs=dg)
        cnt = []
        for l in range(m+1):
            cl = [cf for mon, cf in R.terms() if mon[4] == l]
            cnt.append((sum(1 for x in cl if x < 0), len(cl)))
        check(f"chart {nm} {lab}: all Y^l coefficient polynomials (l=0..{m}) nonnegative; doc counts {ex}", all(cn[0] == 0 for cn in cnt) and tuple(cn[1] for cn in cnt) == ex, cnt)

VAL = {tuple(L.F_VALUES[kk]): kk for kk in range(82)}
def img(g, jx): return VAL.get(tuple(int(x) for x in L.g_act_values(g, list(L.F_VALUES[jx]))))
for nm, src, tgt in (("core action (67,55,53,54,59,2) -> (7,0,5,35,45,1)", (67,55,53,54,59,2), (7,0,5,35,45,1)),
                     ("transposed core action (4,2,33,5,6,46,45) -> (0,5,2,29,7,40,39)", (4,2,33,5,6,46,45), (0,5,2,29,7,40,39)),
                     ("SOC action (1,9,22) -> (5,45,40)", (1,9,22), (5,45,40)),
                     ("stabilizer (0,7,2,29) -> (7,0,29,2) fixing 1,5,35,45", (0,7,2,29,1,5,35,45), (7,0,29,2,1,5,35,45))):
    sols = [g for g in L.G_ELEMENTS if all(img(g, s_) == t_ for s_, t_ in zip(src, tgt))]
    check("some g in G realizes the " + nm, len(sols) > 0, "%d elements, e.g. (rho, rho2, eps) = %s" % (len(sols), sols[0] if sols else None))

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
    numeric("(6) SOC cj=1, piece f35 + j f1 + c f45", "worst %.2e" % worst([{35:1, 1:jj, 45:1/jj} for jj in [RR() for _ in range(n)]]))
    numeric("(6) old-H c=1+1/e, piece f35 + c f45 + e f5", "worst %.2e" % worst([{35:1, 45:1+1/ee, 5:ee} for ee in [RR() for _ in range(n)]]))
    cs = []; cs2 = []; cs3 = []
    for _ in range(n):
        aa = RR(); bb = (1+RR())/aa; cs.append({35:1, 0:bb, 7:aa, 5:(aa+1)*(bb+1)/(aa*bb-1)})
        ee = RR(0.05, (aa+1)*(bb+1)/(aa*bb-1)); cs2.append({35:1, 0:bb, 7:aa, 5:ee, 45:1/ee - (aa*bb-1)/((aa+1)*(bb+1))})
        bb = RR(); aa = min(RR(0.01, 1/bb) if bb >= 1 else RR(0.01, bb), bb, 0.999/bb); w_ = (aa*bb)**0.5; ee = RR()
        cs3.append({35:1, 0:bb, 7:aa, 5:ee, 45:1/ee + bb*(1-w_)**2/(bb-w_)**2})
    numeric("(6) four-coefficient K=0, piece f35 + e f5 + b f0 + a f7", "worst %.2e" % worst(cs))
    numeric("(7) coherent bound ab>1 at equality, piece f35 + c f45 + e f5 + b f0 + a f7", "worst %.2e" % worst(cs2))
    numeric("(7) coherent bound ab<1, b>=a at equality", "worst %.2e" % worst(cs3))
    numeric("SOC as=1, piece f40 + a f5 + s f45", "worst %.2e" % worst([{40:1, 5:aa, 45:1/aa} for aa in [RR() for _ in range(n)]]))
    def Jtrue(pp, bb, cc, qq):
        be = np.linspace(1e-7, 1-1e-7, 40001)
        d1 = np.minimum(pp/2, be); d2 = np.minimum(qq/2, 1-be)
        return float((1 - d1 - d2 + be**2/bb + (1-be)**2/cc + d1**2/pp + d2**2/qq).min())
    numeric("energy threshold a=J (numerically minimized (10)), piece f40 + a f5 + p f0 + b f2 + c f29 + q f7", "worst %.2e" % worst([{40:1, 5:Jtrue(pp,bb,cc,qq), 0:pp, 2:bb, 29:cc, 7:qq} for pp, bb, cc, qq in [(RR(), RR(), RR(), RR()) for _ in range(n)]]))
    G40n = loadGnum('independent_moc4_T36_bad40_full_radial_20261022.json', ['p','r','a','b','c','q','d','s'])
    N = 20000; wv = 1e9
    for _ in range(N):
        pp, bb, cc, qq, rr_ = RR(0.02, 20), RR(0.02, 20), RR(0.02, 20), RR(0.02, 20), RR(0.001, 50)
        J = Jtrue(pp, bb, cc, qq); a_ = J*nrng.uniform(0, 1); s_ = min(nrng.uniform(0, 1)/max(a_, 1e-9), 1e6)
        wv = min(wv, G40n([pp, rr_, a_, bb, cc, qq, 1, s_]) - (1-a_*s_)**2*(1-a_/J)**2)
    numeric("(14) G40 - (1-as)^2(1-a/J)^2 with J the numerically minimized energy, %d random points (tests tables and KKT chart coverage together)" % N, "min %.3e (must be >= 0)" % wv)
    main_stmt(37)

finish()
