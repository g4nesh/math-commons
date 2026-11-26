#!/usr/bin/env python3
# Audit replay for supplement document independent_moc4_T47_IIIb_plus_unswapped_uniform_20261022 (auditor C22-T47).
# Own code only: nothing from the supplement is executed, imported or read.
# Rebuilds the source vector c(g) from the paper's f_j (lib moc4lib), the cubic C of (A.2),
# the height equation, the 14x14 Hessian determinant and the radial factor G, then replays
# the document's exact claims. Prints PASS/FAIL per check; exits 1 if an EXACT check fails.
# Run: python -I -B independent_moc4_T47_IIIb_plus_unswapped_uniform_20261022__check.py   (MOC4_LIB may point to the lib directory)
# ---- core (own code): rebuild c(g), the cubic C, the height equation and the radial octic G ----
import sys, os, time, json
AUD_LIB = os.environ.get('MOC4_LIB') or os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'lib')
sys.path.insert(0, AUD_LIB)
import moc4lib as L
import sympy as sp
import flint
FAIL = []
def check(name, ok, info=''):
    print(('PASS ' if ok else 'FAIL ') + name + ((' : ' + str(info)) if info != '' else ''), flush=True)
    if not ok: FAIL.append(name)
NAMES = ['q','p','j','c','e','d','b','a']
def build_frame(gens):
    """Return ctx8, vars, C0 (constant of C(c+h kappa)), Q (h-coefficient), G (radial octic), prefactor, Hmat-builder."""
    names9 = NAMES + ['h']
    ctx9 = flint.fmpz_mpoly_ctx.get(tuple(names9), 'deglex')
    V9 = dict(zip(names9, ctx9.gens()))
    cvec = []
    for k in range(14):
        s = ctx9.from_dict({})
        for nm, j in zip(NAMES, gens):
            f = L.F_VECS[j][k]; assert f.denominator == 1
            s += int(f.numerator) * V9[nm]
        cvec.append(s)
    w = [cvec[k] + (V9['h'] if k == 13 else 0) for k in range(14)]
    ms, Cexpr = L.cubic_sympy()
    def lin(e):
        e = sp.Poly(e, *ms); out = ctx9.from_dict({})
        for mon, co in e.terms():
            t = ctx9.from_dict({tuple([0]*9): int(co)})
            for i, ex in enumerate(mon):
                if ex: t *= w[i]**ex
            out += t
        return out
    H = [[lin(sp.diff(Cexpr, ms[i], ms[k])) for k in range(14)] for i in range(14)]
    Cw = L.cubic_C(w); Cd = Cw.to_dict()
    ctx8 = flint.fmpz_mpoly_ctx.get(tuple(NAMES), 'deglex')
    to8 = lambda dd: ctx8.from_dict({m[:8]: v for m, v in dd.items()})
    hdeg = max(m[8] for m in Cd)
    C0 = to8({m: v for m, v in Cd.items() if m[8] == 0}); Q = to8({m: v for m, v in Cd.items() if m[8] == 1})
    def bareiss(M, ctx, nv):
        M = [row[:] for row in M]; n = len(M); sign = 1; prev = ctx.from_dict({tuple([0]*nv): 1})
        for k in range(n-1):
            cand = [r for r in range(k, n) if not M[r][k].is_zero()]
            if not cand: return ctx.from_dict({})
            r = min(cand, key=lambda r: len(M[r][k]))
            if r != k: M[k], M[r] = M[r], M[k]; sign = -sign
            piv = M[k][k]
            for i in range(k+1, n):
                for jj in range(k+1, n):
                    M[i][jj] = (M[i][jj]*piv - M[i][k]*M[k][jj]) / prev
                M[i][k] = ctx.from_dict({})
            prev = piv
        return sign * M[n-1][n-1]
    D = bareiss(H, ctx9, 9); Dd = D.to_dict(); dh = max(m[8] for m in Dd)
    coef = [to8({m: v for m, v in Dd.items() if m[8] == k}) for k in range(dh+1)]
    tot = ctx8.from_dict({})
    for k in range(dh+1):
        tot += coef[k] * (-C0)**k * Q**(dh-k)
    num = tot / Q**(dh-2); assert num * Q**(dh-2) == tot      # = Q^2 det H(w(g,-C0/Q))
    fs = num.factor_squarefree()
    pf = ctx8.from_dict({tuple([0]*8): 1})
    for f, e in fs[1]:
        if e == 1: pf *= f
    G2 = num / (-4*pf); assert G2 * (-4*pf) == num
    G = G2.sqrt(); assert G*G == G2
    mc = G.term_content(); G = G / mc; pf = pf * mc * mc        # strip monomial factors into the prefactor
    assert -4 * pf * G * G == num
    return dict(ctx=ctx8, V=dict(zip(NAMES, ctx8.gens())), C0=C0, Q=Q, hdeg=hdeg, G=G, pref=pf, H=H, ctx9=ctx9, bareiss=bareiss, to8=to8)
def coeffs_in(F, var, ctx):
    i = NAMES.index(var); dd = F.to_dict(); n = max(m[i] for m in dd)
    return [ctx.from_dict({tuple(x if t != i else 0 for t, x in enumerate(m)): v for m, v in dd.items() if m[i] == k}) for k in range(n+1)]
def chartsub(F, subs, cctx):
    """F(x_i = num_i/den_i) * prod den_i^(deg_i F); subs: name -> (num, den) in cctx."""
    dd = F.to_dict(); degs = [max(m[i] for m in dd) for i in range(8)]
    one = cctx.from_dict({tuple([0]*cctx.nvars()): 1})
    pw = {}
    def P_(i, kind, e):
        key = (i, kind, e)
        if key not in pw:
            base = subs[NAMES[i]][0 if kind == 'n' else 1]
            pw[key] = base**e if e else one
        return pw[key]
    out = cctx.from_dict({})
    for m, v in dd.items():
        t = one * int(v)
        for i in range(8):
            t = t * P_(i, 'n', m[i]) * P_(i, 'd', degs[i]-m[i])
        out += t
    return out
def coeff_stats(F):
    vals = [int(v) for v in F.coeffs()]
    return len(vals), sum(1 for v in vals if v < 0), (min(vals) if vals else None)
def has_marker(F, zero_vars, cctx):
    """True if F has a positive monomial free of every variable in zero_vars; returns count."""
    idx = [cctx.names().index(z) for z in zero_vars]
    n = sum(1 for m, v in F.to_dict().items() if all(m[i] == 0 for i in idx) and int(v) > 0)
    return n

def get_frame(gens):
    t0 = time.time(); FR = build_frame(gens)
    print('INFO radial octic rebuilt from the cubic in %.1f s' % (time.time()-t0), flush=True)
    check('height equation C(c(g)+h kappa) is affine in h', FR['hdeg'] == 1)
    return FR

# ---- document-specific part: independent_moc4_T47_IIIb_plus_unswapped_uniform_20261022 ----
FR = get_frame([2,9,1,45,5,35,0,7])
ctx, Vv, G = FR['ctx'], FR['V'], FR['G']
check('G has 965 terms and total degree 8', len(G) == 965 and G.total_degree() == 8)
cctx = flint.fmpz_mpoly_ctx.get(('A','P','T','U','V','S','Q'), 'deglex')
A_,P_,T_,U_,V_,S_,Q_ = cctx.gens(); one = cctx.from_dict({(0,)*7: 1})
j_ = 1 + T_; b_ = j_ + A_
# a = theta(V)/b, p = theta(P)/T, q = Q ; A0 = a + b = (V + b^2(1+V))/(b(1+V))
A0n = V_ + b_*b_*(1+V_); A0d = b_*(1+V_)
# stationary range: e = A0^2/j + U, c = theta(S)/(j + A0)
e_n = A0n*A0n + U_*j_*A0d*A0d; e_d = j_*A0d*A0d
c_n = S_*A0d; c_d = (1+S_)*(j_*A0d + A0n)
subs = {'j': (j_, one), 'b': (b_, one), 'a': (V_, b_*(1+V_)), 'p': (P_, (1+P_)*T_), 'q': (Q_, one),
        'e': (e_n, e_d), 'c': (c_n, c_d), 'd': (one, one)}
t0 = time.time(); Gs = chartsub(G, subs, cctx); n, neg, mn = coeff_stats(Gs)
check('stationary chart: cleared G has 134717 terms, 24 negative (document)', (n, neg) == (134717, 24), 'terms=%d negative=%d min=%s (%.0f s)' % (n, neg, mn, time.time()-t0))
del Gs
# stationary-range two-guard certificate: the 11 positive rational weights are DATA transcribed from the
# document's companion verifier (read, not executed); the identity and the signs are recomputed here.
CERT = [('64/3', 0, (0,0,1,2,5,1,2)), ('12', 0, (0,0,1,2,6,1,2)), ('96/7', 0, (0,1,0,1,6,1,2)), ('96/7', 0, (0,1,0,2,6,1,1)),
        ('96/7', 0, (0,1,1,2,6,1,2)), ('13/2', 2, (0,0,0,1,3,0,2)), ('6/7', 2, (0,0,0,1,4,0,2)), ('13/2', 2, (0,0,0,2,3,0,1)),
        ('6/7', 2, (0,0,0,2,4,0,1)), ('13/2', 2, (0,0,1,2,3,0,2)), ('6/7', 2, (0,0,1,2,4,0,2))]
from fractions import Fraction as Fr
class RFn:
    def __init__(s, n, d=None):
        s.n = n; s.d = d if d is not None else one
        g = s.n.gcd(s.d)
        if not g.is_one(): s.n = s.n / g; s.d = s.d / g
        if int(s.d.leading_coefficient()) < 0: s.n = -s.n; s.d = -s.d
    def __add__(s, o): o = o if isinstance(o, RFn) else RFn(o*one); return RFn(s.n*o.d + o.n*s.d, s.d*o.d)
    __radd__ = __add__
    def __sub__(s, o): o = o if isinstance(o, RFn) else RFn(o*one); return RFn(s.n*o.d - o.n*s.d, s.d*o.d)
    def __mul__(s, o): o = o if isinstance(o, RFn) else RFn(o*one); return RFn(s.n*o.n, s.d*o.d)
    __rmul__ = __mul__
    def __truediv__(s, o): o = o if isinstance(o, RFn) else RFn(o*one); return RFn(s.n*o.d, s.d*o.n)
Ar,Pr,Tr,Ur,Vr,Sr,Qr = [RFn(x) for x in cctx.gens()]
jr = Tr + 1; br = jr + Ar; ar = Vr/((Vr+1)*br); Fr_ = ar + br; pr = Pr/(Tr*(Pr+1)); er = Fr_*Fr_/jr + Ur; cr = (Sr/(Sr+1))/(Fr_+jr)
fpc = ar*ar - pr*cr
rr = (pr - cr)*jr + (br - cr)*(pr + 1 + ar) - pr*ar
fphi = rr*rr - 4*ar*cr*(pr - br)*(pr + 1 + ar)
check('guard denominators have only positive coefficients', all(int(x) > 0 for x in fpc.d.coeffs()) and all(int(x) > 0 for x in fphi.d.coeffs()))
Gst = chartsub(G, subs, cctx)
Mpc = cctx.from_dict({}); Mph = cctx.from_dict({})
for w, i, ex in CERT:
    ww = Fr(w); assert ww > 0 and (42*ww).denominator == 1
    t = cctx.from_dict({ex: int(42*ww)})
    if i == 0: Mpc += t
    else: Mph += t
R0 = 42*Gst - Mpc*fpc.n - Mph*fphi.n
nR, negR, _ = coeff_stats(R0)
check('stationary certificate: R0 = 42 Gtilde - Mpc F_pc - MPhi F_Phi has only positive coefficients', negR == 0, 'R0 terms=%d (document 134706), negative=%d' % (nR, negR))
check('R0 has positive monomials with U = Q = 0', has_marker(R0, ['U','Q'], cctx) > 0, has_marker(R0, ['U','Q'], cctx))
check('the two multipliers have only positive coefficients (5 and 6 monomials)', (len(Mpc), len(Mph)) == (5, 6))

del Gst, R0
# cross range: R = 1 + (A0/j - 1) theta(U), e = j R^2, c = C*(A0,j,e) theta(S), C* = (A0+e)/(A0 j (1+R)^2)
# R = (j(1+U) + (A0 - j) U)/(j(1+U)) = Rn/Rd with A0 = A0n/A0d
Rn = j_*A0d*(1+U_) + (A0n - j_*A0d)*U_; Rd = j_*A0d*(1+U_)
e2n = j_*Rn*Rn; e2d = Rd*Rd
# C* = (A0 + e)/(A0 j (1+R)^2) = (A0n e2d + e2n A0d) Rd^2 / (A0d e2d * A0n j (Rd+Rn)^2) ; simplify: e2d = Rd^2
c2n = (A0n*e2d + e2n*A0d) * S_
c2d = A0d * A0n * j_ * (Rd + Rn)**2 * (1+S_)
subs2 = dict(subs); subs2['e'] = (e2n, e2d); subs2['c'] = (c2n, c2d)
sp_ok = True
# sanity (sympy): C* formula equals (A0+e)/(A0 (sqrt j + sqrt e)^2) on e = j R^2
As,Vs,Ts,Us = sp.symbols('A V T U', positive=True)
js = 1+Ts; bs = js+As; a_s = Vs/(1+Vs)/bs; A0s = a_s + bs
Rs = 1 + (A0s/js - 1)*Us/(1+Us); es = js*Rs**2
check('cross chart: e=jR^2 gives sqrt(je) = jR and the C* branch formula', sp.simplify((A0s+es)/(A0s*(sp.sqrt(js)+sp.sqrt(es))**2) - (A0s+es)/(A0s*js*(1+Rs)**2)) == 0)
check('cross chart: R -> A0/j as U -> oo (je -> A0^2)', sp.simplify(sp.limit(Rs, Us, sp.oo) - A0s/js) == 0)
if os.environ.get('MOC4_SKIP_CROSS') != '1':
    t0 = time.time(); Gc = chartsub(G, subs2, cctx); n2, neg2, mn2 = coeff_stats(Gc)
    mk = has_marker(Gc, ['U','Q'], cctx)
    check('cross chart: cleared G has only positive coefficients', neg2 == 0, 'terms=%d negative=%d markers(U=Q=0)=%d (%.0f s); document: 1232529 terms' % (n2, neg2, mk, time.time()-t0))

# ---- NUMERIC test of the stationary-range claim (pc <= a^2, Phi >= 0 => G > 0) ----
# ---- NUMERIC sampler of spectrally realizable points on a contact stratum (own code, lib helpers) ----
import numpy as np
def evalf(F, pt):
    """evaluate an 8-variable fmpz_mpoly at a float point dict"""
    x = [pt[n] for n in NAMES]
    s = 0.0
    for m, v in F.to_dict().items():
        t = float(int(v))
        for i, ex in enumerate(m):
            if ex: t *= x[i]**int(ex)
        s += t
    return s
def sample_realizable(gens, nwant, seed0, maxseeds=None, margin=0.0):
    T = frozenset(range(24))
    for jg in gens: T &= L.F_ZERO[jg]
    M = np.array([[float(L.F_VECS[jg][k]) for jg in gens] for k in range(14)])
    out = []; s = seed0; tried = 0
    while len(out) < nwant and (maxseeds is None or tried < maxseeds):
        tried += 1
        r = L.random_point_on_contact_stratum(T, seed=s, tries=6, tol=1e-10, margin=margin); s += 1
        if r is None: continue
        m = L.matching_vector_np(r['a'], r['b'], r['alpha'])
        F = L.F_values_np(r['a'], r['b'], r['alpha'])
        cvec = np.real(m).copy(); cvec[13] -= np.real(F[r['pi_star']])
        lam, *_ = np.linalg.lstsq(M, cvec, rcond=None)
        resid = np.max(np.abs(M @ lam - cvec))
        if resid > 1e-7 or lam.min() < -1e-8: continue
        dsc = lam[NAMES.index('d')]
        if dsc <= 1e-9: continue
        pt = {n: float(v)/dsc for n, v in zip(NAMES, lam)}
        out.append(pt)
    return out, tried
def rel(F, pt):
    """|F(pt)| divided by the sum of |terms| (a scale-free size of the value)"""
    x = [pt[n] for n in NAMES]; s = 0.0; sa = 0.0
    for m, v in F.to_dict().items():
        t = float(int(v))
        for i, ex in enumerate(m):
            if ex: t *= x[i]**int(ex)
        s += t; sa += abs(t)
    return s / sa if sa else 0.0

rng = np.random.default_rng(11); tested = 0; mn = np.inf
th = lambda X: X/(1+X)
while tested < 3000:
    A, P, T, V, S, U, Qv = np.exp(rng.uniform(-4, 4, 7))
    j_ = 1+T; b_ = j_+A; a_ = th(V)/b_; p_ = th(P)/T; A0 = a_+b_
    e_ = A0*A0/j_ + U; c_ = th(S)/(j_+A0)
    pt = dict(q=Qv, p=p_, j=j_, c=c_, e=e_, d=1.0, b=b_, a=a_)
    Ph = ((p_-c_)*j_ + (b_-c_)*(p_+1+a_) - p_*a_)**2 - 4*a_*c_*(p_-b_)*(p_+1+a_)
    if p_*c_ > a_*a_ or Ph < 0: continue
    tested += 1; mn = min(mn, rel(G, pt))
print('NUMERIC stationary range, pc<=a^2, Phi>=0: %d log-uniform chart samples (seed 11): min G/scale = %.3e (positive means G>0)' % (tested, mn), flush=True)

print('SUMMARY: %d EXACT check(s) failed' % len(FAIL) if FAIL else 'SUMMARY: all EXACT checks passed')
sys.exit(1 if FAIL else 0)
