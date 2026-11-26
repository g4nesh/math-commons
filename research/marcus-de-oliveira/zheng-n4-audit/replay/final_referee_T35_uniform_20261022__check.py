"""Independent exact replay for supplement document final_referee_T35_uniform_20261022.
Own code only: rebuilds generators, profile identities, the radial cubic restriction C(c + h kappa)
and the full 14x14 Hessian of the paper's cubic (A.2) from moc4lib, and replays the positivity tables.
The degree-nine polynomial G is loaded as DATA from the supplement's JSON (never executed) and is
validated against the Hessian computed here by an exact randomized identity test.
Usage: MOC4_SUPP=<supplement root> python -I -B final_referee_T35_uniform_20261022__check.py   (or pass the root as argv[1])
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
T = L.CONTACT_SETS[36]
check("doc contact set {2,3,4,6,11,12} equals paper T_36 literally", sorted(T) == [2,3,4,6,11,12])
gens = L.generators(T)
check("generators(T_36) == source 0,1,2,7,8,9,10,53,67,68,69", gens == [0,1,2,7,8,9,10,53,67,68,69], gens)
check("nonminor generators are exactly 68, 69", [j for j in gens if not L.F_IS_MINOR[j]] == [68, 69])
check("rank of source = 8", L.rank([L.F_VECS[j] for j in gens]) == 8)
V = L.F_VALUES
def ident(lhs, rhs): return all(sum(V[j][s] for j in lhs) == sum(V[j][s] for j in rhs) for s in range(24))
check("f68 + f69 = f1 + f9 + f10 + 2 f67", ident([68,69],[1,9,10,67,67]))
check("f68 + f2 = f0 + f9 + f67", ident([68,2],[0,9,67]))
check("f68 + f8 = f7 + f10 + f67", ident([68,8],[7,10,67]))
check("f69 + f0 = f1 + f2 + f10 + f67", ident([69,0],[1,2,10,67]))
check("f69 + f7 = f1 + f8 + f9 + f67", ident([69,7],[1,8,9,67]))
names = ['q','p','j','c','e','d','b','a']
pn = ['q','p','j','c','e','b','a']; c7 = flint.fmpz_mpoly_ctx.get(tuple(pn), 'lex')
GG = {}
for tag, B, fn in (('68', [0,1,7,9,53,68,67,10], 'final_referee_T35_bad68_full_radial_20261022.json'),
                   ('69', [1,2,8,9,53,69,67,10], 'final_referee_T35_bad69_full_radial_20261022.json')):
    check(f"representation ({tag}) basis independent", L.rank([L.F_VECS[j] for j in B]) == 8)
    ctx, Cst, Q = radial_CQ(B, names)
    dd = json.load(open(os.path.join(SUPP, 'proofs/manin_transfer_20261021/' + fn)))
    check(f"({tag}) own C*, Q equal JSON C_terms, Q_terms", load_json_poly(ctx, dd['C_terms'], 8) == Cst and load_json_poly(ctx, dd['Q_terms'], 8) == Q)
    G = -load_json_poly(ctx, dd['radial_factors'][1]['terms'], 8)
    check(f"({tag}) radial factor0 is +-Q; minus stored G has d^9 coefficient +1; term count", load_json_poly(ctx, dd['radial_factors'][0]['terms'], 8) in (Q, -Q) and [int(x) for m, x in G.terms() if m[5] == 9] == [1], len(G))
    kk = hessian_identity_test(B, ctx, Cst, Q, G, npts=10, seed=7)
    check(f"({tag}) (2) Q^2 det Hess C(c - (C*/Q)kappa) == -4 G^2 at 10 random integer points (exact; Schwartz-Zippel)", kk == 10, kk)
    q,p,j,c,e,d,b,a = ctx.gens()
    Gd = G.subs({'d':1}); Cd = Cst.subs({'d':1})
    GG[tag] = c7.from_dict({tuple(m[:5]) + tuple(m[6:]): int(cf) for m, cf in Gd.terms()})
    if tag == '68':
        check("(68) C* and Q have only positive coefficients, containing d^3 and d^2", all(x > 0 for x in Cst.coeffs()) and all(x > 0 for x in Q.coeffs()) and any(m[5] == 3 for m, x in Cst.terms()) and any(m[5] == 2 for m, x in Q.terms()))
        check("(5) baseline G68|p=e=a=c=0 == (j+1)^3(b+1)^3(q+1)^3", Gd.subs({'p':0,'e':0,'a':0,'c':0}) == (j+1)**3*(b+1)**3*(q+1)**3)
        check("G68 degrees in (q,p,j,c,e,b,a) == (3,2,3,2,2,3,2)", tuple(GG[tag].degrees()) == (3,2,3,2,2,3,2), GG[tag].degrees())
    else:
        negs = [(m, int(x)) for m, x in Cst.terms() if x < 0]
        check("(69) only negative monomial of C* is -a e j", negs == [((0,0,1,0,1,0,0,1), -1)], negs)
        check("(1) C* - (1 + a + e + ae(1-j)) has nonnegative coefficients (d=1)", all(x >= 0 for x in (Cd - (1 + a + e + a*e - a*e*j)).coeffs()))
        B0 = (e+1)*(p+1)*((j+1)*(1-p*j) + e*(j**2+j+1-j*p))
        check("(9)/(10) baseline G69|a=c=b=q=0 == B0", Gd.subs({'a':0,'c':0,'b':0,'q':0}) == B0)
        check("G69 degrees in (q,p,j,c,e,b,a) <= (2,2,2,2,2,3,2)", all(int(x) <= y for x, y in zip(GG[tag].degrees(), (2,2,2,2,2,3,2))), GG[tag].degrees())
# branch 68 tables
G68 = GG['68']
tn = ['J','X','Y','P','E','A','C']; tsy = sympy.symbols(tn); J,X,Y,P,E,A,C = tsy
tctx = flint.fmpz_mpoly_ctx.get(tuple(tn), 'lex')
F = lambda ex: sym_to_flint(sympy.expand(ex), tctx, tsy)
qlo = {'q': 'X/(J*(1+X))', 'p': 'P'}; qhi = {'q': '(1+X)/J', 'p': '(1+X+J)*(J+1)*P/(J*X*(1+P))'}
blo = {'b': 'Y/(J*(1+Y))', 'e': 'E'}; bhi = {'b': '(1+Y)/J', 'e': '(1+Y+J)*(J+1)*E/(J*Y*(1+E))'}
Hq = J**2*(X+1) - X
expZ = {'00': (6190, 5893), '01': (8949, 8802), '10': (8705, 8439), '11': (12375, 12197)}
for ch, s1, s2 in (('00', qlo, blo), ('01', qlo, bhi), ('10', qhi, blo), ('11', qhi, bhi)):
    s = {}; s.update(s1); s.update(s2); s.update({'j': 'J', 'a': 'A', 'c': 'C'})
    Z, _ = clear_subst(G68, pn, s, tn, degs=(3,2,3,2,2,3,2))
    if ch == '00':
        H = (J*Y+J+Y)*Hq - E*J*X*(Y+1); M = A*C*P*J**2*(X+1)*(Y+1)*(A*C+A*P+C*P)
    elif ch == '01':
        H = (J+Y+1)*Hq; M = A*C*P*J**4*Y**2*(X+1)*(A*C+A*P+C*P)
    else:
        H = J**2 - X - 1
        if ch == '10': F0 = J**2*(J+1)*(J+X+1)*(J*Y+J+Y)**2; Y0 = Y+1
        else: F0 = J**4*Y**2*(J+1)*(J+X+1)*(J+Y+1)**2; Y0 = 1
        M = F0*Y0*P*A*C*(P*A*(C*J*X + (J+1)**2*(X+1)) + P*C*(J+X+1)**2 + A*C*J*X + P*(J+X+1)*J*(X+1))
    check(f"(68) chamber {ch}: Z term count == {expZ[ch][0]}", len(Z) == expZ[ch][0], len(Z))
    check(f"(68) chamber {ch}: multiplier M has nonnegative coefficients", all(x >= 0 for x in F(M).coeffs()))
    R = Z - F(M*H**2)
    if ch == '00':
        Nset = [(5,3,0,1,0,1,2), (7,3,0,1,0,2,1), (7,3,1,1,0,2,1)]
        neg = sorted((tuple(m), int(cf)) for m, cf in R.terms() if cf < 0)
        check("(68) chamber 00: Z - M H^2 has exactly the three negative monomials of (7), each -1", neg == sorted((nu, -1) for nu in Nset), neg)
        corr = 0
        for nu in Nset:
            mon = list(nu); mon[0] -= 1
            corr = corr + F(sympy.Mul(*[v**kx for v, kx in zip(tsy, mon)])*(J-1)**2)
        R = 2*R - corr
    nc = negcount(R)
    check(f"(68) chamber {ch}: remainder has no negative coefficient (doc count {expZ[ch][1]})", nc[0] == 0, nc)
# branch 69 tables
G69 = GG['69']
tn = ['A','C','B','Q','t','x','u']; tsy = sympy.symbols(tn); A,C,B,Q,t,x,u = tsy
tctx = flint.fmpz_mpoly_ctx.get(tuple(tn), 'lex')
F = lambda ex: sym_to_flint(sympy.expand(ex), tctx, tsy)
h = t/(1+t)
rows = [('both low v<1', t/(1+t), x/(1+x), u/(1+u), A, C, (3226, 3199)),
        ('both low v>1', None, (h+x)/(1+x), 1+u, A, C, (41018, 41018)),
        ('j>1 v<1', None, 1+x, u/(1+u), 'H', C, (14406, 14406)),
        ('j>1 v>1', None, 1+x, 1+u/(x*(1+u)), 'H', C, (26220, 26220)),
        ('p>1 v<1', 1+x, h/(1+x), u/(1+u), A, C/(x*(1+C)), (6760, 6760)),
        ('p>1 v>1', 1+x, h/(1+x), 1+u, A, C/(x*(1+C)), (22059, 22059))]
pS, jS, eS, aS, cS, bS, vS = sympy.symbols('p j e a c b v', positive=True)
for nm, pp, jj, vv, aa, cc, ex in rows:
    if pp is None: pp = h/jj
    ee = vv/jj
    if aa == 'H': aa = (1/x - ee)*A/(1+A)
    bb = B if 'v<1' in nm else ee*B/((vv-1)*(1+B))
    s = {kx: sympy.cancel(sympy.together(v)) for kx, v in dict(q=Q, p=pp, j=jj, c=cc, e=ee, b=bb, a=aa).items()}
    Z, _ = clear_subst(G69, pn, s, tn, degs=(2,2,2,2,2,3,2))
    R = Z
    if nm == 'both low v<1':
        R = Z - F(x**2*(t-x)**2*(B+1)**2*A*C*Q*(A*C+A*Q+C*Q))
    nc = negcount(R)
    check(f"(69) row '{nm}': cleared polynomial (minus square (12) in row 1) has no negative coefficient; doc counts {ex}", nc[0] == 0 and len(Z) == ex[0], (len(Z), nc))
    # 1 - L a etc. in the chart
    La = sympy.Piecewise((0, True)) if True else None
# chart factors: 1 - La a = 1/(1+A) in the j>1 rows; 1 - Lc c = 1/(1+C) in the p>1 rows; 1 - Lb b = 1/(1+B) in v>1 rows
x_, A_, C_, B_, u_, e_ = sympy.symbols('x A C B u e', positive=True)
jj = 1+x_; ee = e_; aa = (1/x_ - ee)*A_/(1+A_); La = (jj-1)/(1-ee*(jj-1))
check("chart (11): 1 - La a == 1/(1+A) for j = 1+x", sympy.simplify(1 - La*aa - 1/(1+A_)) == 0)
pp = 1+x_; cc = C_/(x_*(1+C_))
check("chart (11): 1 - Lc c == 1/(1+C) for p = 1+x", sympy.simplify(1 - (pp-1)*cc - 1/(1+C_)) == 0)
vv = sympy.Symbol('v', positive=True); jv = sympy.Symbol('jv', positive=True); ev = vv/jv; bb = ev*B_/((vv-1)*(1+B_))
check("chart (11): 1 - Lb b == 1/(1+B) for v = je > 1", sympy.simplify(1 - (jv*ev-1)/ev*bb - 1/(1+B_)) == 0)

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
    cs = []
    for _ in range(n):
        qq = RR(); jj = (1+RR())/qq; cs.append({68:1, 0:qq, 7:jj, 1:(qq+1)*(jj+1)/(qq*jj-1)})
    numeric("(3a) four-coefficient at equality, piece f68 + q f0 + j f7 + p f1", "worst %.2e" % worst(cs))
    cs = []
    for _ in range(n):
        jj = RR(); bb = (1+RR())/jj; cs.append({68:1, 7:jj, 67:bb, 53:(jj+1)*(bb+1)/(jj*bb-1)})
    numeric("(3b) four-coefficient at equality, piece f68 + j f7 + b f67 + e f53", "worst %.2e" % worst(cs))
    numeric("(4) SOC pj=1, piece f69 + p f2 + j f8", "worst %.2e" % worst([{69:1, 2:pp, 8:1/pp} for pp in [RR() for _ in range(n)]]))
    numeric("(4) old-H c(p-1)=1, piece f69 + c f9 + p f2", "worst %.2e" % worst([{69:1, 9:cc, 2:1+1/cc} for cc in [RR() for _ in range(n)]]))
    numeric("(4) H-sum (a+e)(j-1)=1, piece f69 + a f10 + e f53 + j f8", "worst %.2e" % worst([{69:1, 10:aa, 53:ee, 8:1+1/(aa+ee)} for aa, ee in [(RR(), RR()) for _ in range(n)]]))
    numeric("(4) chain j=1/e+1/b, piece f69 + j f8 + e f53 + b f67", "worst %.2e" % worst([{69:1, 8:1/ee+1/bb, 53:ee, 67:bb} for ee, bb in [(RR(), RR()) for _ in range(n)]]))
    N = 20000
    G68n = loadGnum('final_referee_T35_bad68_full_radial_20261022.json', ['q','p','j','c','e','d','b','a'])
    wv = 1e9
    for _ in range(N):
        q_, j_, b_, a_, c_ = [RR(0.02, 20) for _ in range(5)]
        Lp = max((q_*j_-1)/((q_+1)*(j_+1)), 0); Le = max((j_*b_-1)/((b_+1)*(j_+1)), 0)
        p_ = nrng.uniform(0, 1)/Lp if Lp > 0 else RR(0.02, 20); e_ = nrng.uniform(0, 1)/Le if Le > 0 else RR(0.02, 20)
        wv = min(wv, G68n([q_, p_, j_, c_, e_, 1, b_, a_])/((j_+1)**3*(b_+1)**3*(q_+1)**3*(1-Lp*p_)**2*(1-Le*e_)**2) - 1)
    numeric("(5) G68 / bound - 1 on %d random hard points" % N, "min %.3e (must be >= 0)" % wv)
    G69n = loadGnum('final_referee_T35_bad69_full_radial_20261022.json', ['q','p','j','c','e','d','b','a'])
    wv = 1e9; cnt = 0
    while cnt < N:
        q_, p_, j_, c_, e_, b_, a_ = [RR(0.02, 20) for _ in range(7)]
        if p_*j_ >= 1 or c_*(p_-1) >= 1 or (a_+e_)*(j_-1) >= 1 or j_ >= 1/e_ + 1/b_: continue
        cnt += 1
        La = (j_-1)/(1-e_*(j_-1)) if j_ > 1 else 0; Lc = max(p_-1, 0); Lb = max(j_*e_-1, 0)/e_
        B0 = (e_+1)*(p_+1)*((j_+1)*(1-p_*j_) + e_*(j_*j_+j_+1-j_*p_))
        wv = min(wv, G69n([q_, p_, j_, c_, e_, 1, b_, a_])/(B0*(1-La*a_)**2*(1-Lc*c_)**2*(1-Lb*b_)**3) - 1)
    numeric("(10) G69 / bound - 1 on %d random hard points" % N, "min %.3e (must be >= 0)" % wv)
    main_stmt(36)

finish()
