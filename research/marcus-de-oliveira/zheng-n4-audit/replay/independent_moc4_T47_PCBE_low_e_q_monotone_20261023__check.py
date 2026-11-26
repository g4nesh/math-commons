# Replay for document independent_moc4_T47_PCBE_low_e_q_monotone_20261023
# Auditor's own exact replay (Zheng n=4 audit). Nothing from the supplement is executed or read.
# The T47 (= paper T_65 contact set {0,4,5,6,9,12}) radial polynomial G is rebuilt from the paper's
# cubic C (A.2) and the Appendix D generators via the shared library moc4lib (../lib).
import sys, os, time, collections
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.environ.get('MOC4LIB', os.path.join(HERE, '..', 'lib')))
import moc4lib as L
import flint
FAILS = []
def check(name, cond, info=''):
    print(('PASS ' if cond else 'FAIL ') + name + ((' :: ' + str(info)) if info != '' else ''), flush=True)
    if not cond: FAILS.append(name)
def info(name, val):
    print('INFO ' + name + ' :: ' + str(val), flush=True)
def finish():
    print('SUMMARY: %d EXACT check(s) failed' % len(FAILS), FAILS, flush=True)
    sys.exit(1 if FAILS else 0)
NAMES = ['q','p','j','c','e','d','b','a','h']
FACETS = [2,9,1,45,5,35,0,7]          # q p j c e d b a (paper Appendix D numbering)
def build(L):
    import sympy as sp
    ctx = flint.fmpz_mpoly_ctx.get(NAMES, 'lex')
    g = ctx.gens(); h = g[8]
    m = []
    for k in range(14):
        s = ctx.from_dict({})
        for f, cf in zip(FACETS, g[:8]):
            x = L.F_VECS[f][k]; assert x.denominator == 1; s += int(x) * cf
        m.append(s)
    Cc = L.cubic_C(m)
    Qh = m[0]*m[6] - m[1]*m[5]
    mh = m[:13] + [m[13] + h]
    assert L.cubic_C(mh) - Cc - h*Qh == 0
    ms, Cs = L.cubic_sympy()
    Hm = []
    for r in range(14):
        row = []
        for s_ in range(14):
            ex = sp.expand(sp.diff(Cs, ms[r], ms[s_]))
            out = ctx.from_dict({})
            for k in range(14):
                cf = ex.coeff(ms[k])
                if cf != 0: out += int(cf) * mh[k]
            row.append(out)
        Hm.append(row)
    cur = {0: ctx.from_dict({(0,)*9: 1})}
    for r in range(14):
        nxt = {}
        for S, val in cur.items():
            for col in range(14):
                if S >> col & 1: continue
                ent = Hm[r][col]
                if ent == 0: continue
                term = val*ent
                if bin(S >> (col+1)).count('1') % 2: term = -term
                T = S | (1 << col)
                nxt[T] = nxt.get(T, 0) + term
        cur = nxt
    D = cur[(1 << 14) - 1]
    Dk = {}
    for exps, cf in D.to_dict().items():
        ex = list(exps); k = ex[8]; ex[8] = 0
        Dk.setdefault(k, {})[tuple(ex)] = cf
    Dk = {k: ctx.from_dict(v) for k, v in Dk.items()}
    N = sum((Dk[k] * (-Cc)**k * Qh**(4-k) for k in Dk), ctx.from_dict({}))
    quo, rem = divmod(N, Qh**2)
    assert rem == 0
    dvar = g[5]
    G2, rem2 = divmod(quo, -4 * dvar**2)
    assert rem2 == 0
    G = G2.sqrt()          # raises if not a perfect square
    return ctx, g, m, Cc, Qh, sorted(Dk), G

def evalq(P, vals):
    from fractions import Fraction as Fr
    tot = Fr(0)
    for ex, cf in P.to_dict().items():
        r = Fr(int(cf))
        for v, k in zip(vals, ex):
            k = int(k)
            if k: r *= Fr(v)**k
        tot += r
    return tot

def clear(P, chart, cctx):
    """chart: list of 9 (num, den) pairs of fmpz_mpoly in cctx for (q,p,j,c,e,d,b,a,h).
    returns prod den_i^{deg_i P} * P(num/den)."""
    degs = [int(x) for x in P.degrees()]
    pw = {}
    def pwr(i, which, k):
        key = (i, which, k)
        if key not in pw:
            base = chart[i][which]
            pw[key] = base**k
        return pw[key]
    tot = cctx.from_dict({})
    for ex, cf in P.to_dict().items():
        t = cctx.constant(int(cf))
        for i, k in enumerate(ex):
            k = int(k)
            if degs[i] == 0: continue
            if k: t *= pwr(i, 0, k)
            if degs[i] - k: t *= pwr(i, 1, degs[i] - k)
        tot += t
    return tot

class RF:
    """rational function num/den over fmpz_mpoly, kept reduced."""
    def __init__(self, n, d=None):
        if d is None: d = n.context().constant(1)
        g = n.gcd(d)
        if g != 1:
            n = n // g; d = d // g
        # normalize sign of den: leading coefficient positive
        lc = list(d.to_dict().values())[0] if len(d) else 1
        if lc < 0: n, d = -n, -d
        self.n, self.d = n, d
    def __add__(s, o):
        o = o if isinstance(o, RF) else RF(s.n.context().constant(int(o)) if isinstance(o,int) else o)
        return RF(s.n*o.d + o.n*s.d, s.d*o.d)
    __radd__ = __add__
    def __neg__(s): return RF(-s.n, s.d)
    def __sub__(s, o):
        o = o if isinstance(o, RF) else RF(s.n.context().constant(int(o)) if isinstance(o,int) else o)
        return s + (-o)
    def __rsub__(s, o): return (-s) + o
    def __mul__(s, o):
        o = o if isinstance(o, RF) else RF(s.n.context().constant(int(o)) if isinstance(o,int) else o)
        return RF(s.n*o.n, s.d*o.d)
    __rmul__ = __mul__
    def __truediv__(s, o):
        o = o if isinstance(o, RF) else RF(s.n.context().constant(int(o)) if isinstance(o,int) else o)
        return RF(s.n*o.d, s.d*o.n)
    def __rtruediv__(s, o):
        o = o if isinstance(o, RF) else RF(s.n.context().constant(int(o)) if isinstance(o,int) else o)
        return o / s
    def __pow__(s, k): return RF(s.n**k, s.d**k)

def nonneg(P):
    """(all coefficients >= 0, number of terms, min coefficient)"""
    cs = [int(v) for v in P.to_dict().values()]
    return (len(cs) > 0 and min(cs) >= 0), len(cs), (min(cs) if cs else None)

def has_mono_free_of(P, idxs):
    """a positive monomial whose exponents at positions idxs are all zero"""
    for ex, cf in P.to_dict().items():
        if cf > 0 and all(int(ex[i]) == 0 for i in idxs): return True
    return False

def subs_rf(P, vals, cctx):
    """clear P at RF values: returns (cleared numerator, product of dens^deg); P(vals) = cleared/denprod"""
    chart = [(v.n, v.d) for v in vals]
    num = clear(P, chart, cctx)
    degs = [int(x) for x in P.degrees()]
    den = cctx.constant(1)
    for (n_, d_), k in zip(chart, degs):
        if k: den *= d_**k
    return num, den
def conic(ctx, g, G):
    q,p,j,c,e,d,b,a,h = g
    blocks = collections.defaultdict(dict)
    for ex, cf in G.to_dict().items():
        ex = list(ex); ke, kq = int(ex[4]), int(ex[0]); ex[4]=0; ex[0]=0
        blocks[(ke,kq)][tuple(ex)] = cf
    Bk = {k: ctx.from_dict(v) for k, v in blocks.items()}
    z = ctx.constant(0)
    e0 = sum((Bk.get((0,k), z) * q**k for k in range(3)), z)
    e1 = sum((Bk.get((1,k), z) * q**k for k in range(3)), z)
    K, r0 = divmod(e0, (q+b+d)**2); assert r0 == 0 and int(K.degrees()[0]) == 0
    Bq, r1 = divmod(e1, (q+b+d)); assert r1 == 0 and int(Bq.degrees()[0]) <= 1
    bd = collections.defaultdict(dict)
    for ex, cf in Bq.to_dict().items():
        ex = list(ex); kq = int(ex[0]); ex[0] = 0; bd[kq][tuple(ex)] = cf
    B0 = ctx.from_dict(bd.get(0, {})); B1 = ctx.from_dict(bd.get(1, {}))
    A0, A1, A2 = Bk.get((2,0), z), Bk.get((2,1), z), Bk.get((2,2), z)
    # reassemble
    assert G == e**2*(A0 + A1*q + A2*q**2) + e*(q+b+d)*(B0 + B1*q) + (q+b+d)**2*K
    return dict(A0=A0, A1=A1, A2=A2, B0=B0, B1=B1, K=K)

def setup():
    t0 = time.time()
    ctx, g, m, Cc, Qh, hdeg, G = build(L)
    q,p,j,c,e,d,b,a,h = g
    check('T65 generators are the eight facets (2,9,1,45,5,35,0,7) of the document',
          sorted(L.generators(L.CONTACT_SETS[65])) == sorted(FACETS))
    check('the eight generators are linearly independent (coefficients of a source are unique)',
          L.rank([L.F_VECS[f] for f in FACETS]) == 8)
    check('only f35 is non-minor', [L.F_IS_MINOR[f] for f in FACETS] == [True]*5 + [False] + [True]*2)
    target = -d*(q*p + q*e + q*d + q*a + p*e + p*d + p*b + e*d + d**2 + d*b + d*a + b*a)
    check('C(c) = -d(qp+qe+qd+qa+pe+pd+pb+ed+d^2+db+da+ba) (height-equation constant; < 0 when d>0, others >= 0)', Cc == target)
    check('C(c + h kappa) = C(c) + h (m0 m6 - m1 m5) (linear in h)', True)
    check('det Hess C(c + h kappa) has h-degree <= 3', max(hdeg) <= 3, hdeg)
    check('at h = -C(c)/Q: det Hess C = -4 d^2 G^2 / Q^2 with G a polynomial (965 terms, multidegree (2,2,2,2,2,8,3,3))',
          len(G) == 965 and tuple(int(x) for x in G.degrees()[:8]) == (2,2,2,2,2,8,3,3), len(G))
    info('G built in %.1fs' % (time.time()-t0), 'sign fixed so that G > 0 at (q,p,j,c,e,d,b,a)=(0,2,1/2,1/1000,2/3,1,2,1/10): G = %s' % float(evalq(G, [0,2,0.5,0.001,2/3.,1,2,0.1,0])))
    return ctx, g, G

import sympy as sp
ctx, g, G = setup()
# side fact: failure of (2) at L = b, low e  <=>  s = sqrt(c/B) in (0,1) and j < J*
bs, es, ss, js = sp.symbols('b e s j', positive=True)
Bq = 1/es + 1/(bs+1); Kq = 1/js + 1/(bs+1); cq = ss**2*Bq
Js = (bs+1)*(1-ss)**2/(bs + 2*ss + (bs+1)*ss**2/es)
check('chart origin: (c+1)/(1-s)^2 - 1/(b+1) == 1/J*  (so (c+1)/K < (1-s)^2 <=> j < J*)', sp.simplify((cq+1)/(1-ss)**2 - 1/(bs+1) - 1/Js) == 0)
check('chart origin: low e (e < b(b+1)) gives 1/e + 1/(b+1) > 1/b, so B = 1/e + 1/(b+1)', sp.simplify((1/es + 1/(bs+1) - 1/bs) - (bs*(bs+1) - es)/(es*bs*(bs+1))) == 0)
cc = flint.fmpz_mpoly_ctx.get(['A','T','U','V','S','P','Q'], 'lex')
A,T,U,V,S,P,Q = [RF(x) for x in cc.gens()]
one = RF(cc.constant(1)); zero = RF(cc.constant(0))
b = 1 + A; p = b + P; r = V/(1+V); a = r**2/b; s = S/(1+S)
e = b*(b+1)*U/(1+U); c = s**2*(1/e + 1/(b+1))
Jst = (b+1)*(1-s)**2/(b + 2*s + (b+1)*s**2/e); j = Jst*T/(1+T)
vals = [Q, p, j, c, e, one, b, a, zero]
num, den = subs_rf(G, vals, cc)
check('canonical clearing denominator has nonnegative coefficients and is free of Q', nonneg(den)[0] and int(den.degrees()[6]) == 0 and has_mono_free_of(den, [5]))
bl = collections.defaultdict(dict)
for ex, cf in num.to_dict().items():
    ex = list(ex); k = int(ex[6]); ex[6] = 0; bl[k][tuple(ex)] = cf
g0, g1, g2 = [cc.from_dict(bl.get(k, {})) for k in (0, 1, 2)]
check('cleared G is quadratic in Q with (40647, 41139, 38301) terms', max(bl) == 2 and (len(g0), len(g1), len(g2)) == (40647, 41139, 38301), (len(g0), len(g1), len(g2)))
o1 = nonneg(g1)
check('g1: every coefficient a positive integer; a positive monomial free of P', o1[0] and o1[2] > 0 and has_mono_free_of(g1, [5]), o1[1:])
Phi = ((p - c)*j + (b - c)*(p + 1 + a) - p*a)**2 - 4*a*c*(p - b)*(p + 1 + a)
check('reduced denominator D_Phi has nonnegative coefficients (positive on the chart)', nonneg(Phi.d)[0] and has_mono_free_of(Phi.d, [5]))
Ag, Ug, Vg = cc.gens()[0], cc.gens()[2], cc.gens()[3]
R = g2 - 2*Ag*Ug**2*Vg**2*Phi.n
oR = nonneg(R)
check('g2 = R + 2 A U^2 V^2 N_Phi with R of 38304 strictly positive coefficients, a positive monomial free of P', oR[0] and oR[2] > 0 and oR[1] == 38304 and has_mono_free_of(R, [5]), oR[1:])
finish()
