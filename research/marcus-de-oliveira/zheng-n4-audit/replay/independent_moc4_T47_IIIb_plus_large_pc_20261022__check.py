#!/usr/bin/env python3
# Audit replay for supplement document independent_moc4_T47_IIIb_plus_large_pc_20261022 (auditor C22-T47).
# Own code only: nothing from the supplement is executed, imported or read.
# Rebuilds the source vector c(g) from the paper's f_j (lib moc4lib), the cubic C of (A.2),
# the height equation, the 14x14 Hessian determinant and the radial factor G, then replays
# the document's exact claims. Prints PASS/FAIL per check; exits 1 if an EXACT check fails.
# Run: python -I -B independent_moc4_T47_IIIb_plus_large_pc_20261022__check.py   (MOC4_LIB may point to the lib directory)
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

# ---- document-specific part: independent_moc4_T47_IIIb_plus_large_pc_20261022 ----
FR = get_frame([2,9,1,45,5,35,0,7])
ctx, Vv, G = FR['ctx'], FR['V'], FR['G']
check('G has 965 terms and total degree 8', len(G) == 965 and G.total_degree() == 8, (len(G), G.total_degree()))
cctx = flint.fmpz_mpoly_ctx.get(('A','P','T','U','V','S','Q'), 'deglex')
A_,P_,T_,U_,V_,S_,Q_ = cctx.gens(); one = cctx.from_dict({(0,)*7: 1})
iV = cctx.names().index('V')
def certify(label, subs, zero_vars):
    Gt = chartsub(G, subs, cctx)
    n, neg, mn = coeff_stats(Gt)
    dd = Gt.to_dict()
    negs = [m for m, v in dd.items() if int(v) < 0]
    check('%s: negative coefficients are all -2 with V-exponent >= 2' % label,
          all(int(dd[m]) == -2 and m[iV] >= 2 for m in negs), 'terms=%d negatives=%d' % (n, len(negs)))
    W = cctx.from_dict({tuple(x - (2 if t == iV else 0) for t, x in enumerate(m)): 1 for m in negs})
    R = Gt - (V_*V_ - 1)**2 * W
    nR, negR, _ = coeff_stats(R)
    check('%s: R = Gtilde - (V^2-1)^2 W has only nonnegative coefficients' % label, negR == 0, 'R terms=%d, W terms=%d' % (nR, len(W)))
    mk = has_marker(R, zero_vars, cctx)
    check('%s: R has a positive monomial free of %s' % (label, ','.join(zero_vars)), mk > 0, mk)
    return n, len(negs), nR
# chart 1: 1 < j < b
j1 = 1 + T_; b1 = j1 + A_
den_a = b1*(1+V_)
subs1 = {'j': (j1, one), 'b': (b1, one), 'a': (V_, den_a), 'p': (P_, (1+P_)*T_), 'q': (Q_, one),
         'e': (j1*(1+U_)*A_ + j1*j1*U_, (1+U_)*A_), 'c': ((b1*(1+V_) - j1*V_)*S_, j1*b1*(1+V_)*(1+S_)), 'd': (one, one)}
certify('chart 1 (1<j<b)', subs1, ['U','Q'])
# chart 2: j >= b
b2 = 1 + A_; j2 = b2 + T_
subs2 = {'b': (b2, one), 'j': (j2, one), 'a': (V_, (1+V_)*j2), 'p': (P_, (1+P_)*(A_+T_)), 'q': (Q_, one),
         'e': (j2 + U_, one), 'c': (S_, j2*(1+V_)*(1+S_)), 'd': (one, one)}
certify('chart 2 (j>=b)', subs2, ['T','U','Q'])
# coverage of the charts (sympy, exact)
At,Pt,Tt,Ut,Vt,St = sp.symbols('A P T U V S', positive=True)
th = lambda X: X/(1+X)
j = 1+Tt; b = j+At; a = Vt/(b*(1+Vt)); e = j + j**2*th(Ut)/At; Lc = 1/j - a
check('chart 1: e upper limit = bj/(b-j)', sp.simplify(sp.limit(e, Ut, sp.oo) - b*j/(b-j)) == 0)
check('chart 1: p upper limit = 1/(j-1)', sp.simplify(sp.limit(th(Pt)/Tt, Pt, sp.oo) - 1/(j-1)) == 0)
check('chart 1: c numerator/denominator consistent with c = (1/j - a) S/(1+S)', sp.simplify(Lc*th(St) - (b*(1+Vt)-j*Vt)*St/(j*b*(1+Vt)*(1+St))) == 0)
check('chart 1: 1/j - a > 0 identically (b(1+V) - jV = b + AV... )', sp.expand((b*(1+Vt) - j*Vt) - (b + At*Vt)) == 0)
j = 1+At+Tt
check('chart 2: c = (1/j - a) S/(1+S) with a = V/((1+V) j)', sp.simplify((1/j - th(Vt)/j)*th(St) - St/(j*(1+Vt)*(1+St))) == 0)
# numeric probe of Lemma 1 (Gram criterion) on its boundary pc = a^2, j = d^2/(a+c): min Phi over SO(4)
import numpy as np
rng = np.random.default_rng(20261125)
worst = 1e9
for trial in range(12):
    a_, p_, d_ = rng.uniform(0.2, 2.0, 3); c_ = a_*a_/p_; j_ = d_*d_/(a_+c_)
    coef = dict(q=0.0, p=p_, j=j_, c=c_, e=0.0, d=d_, b=0.0, a=a_)
    cv = [sum(coef[nm]*float(L.F_VECS[f][k]) for nm, f in zip(NAMES, [2,9,1,45,5,35,0,7])) for k in range(14)]
    val = L.min_phi_over_SO4(cv, restarts=8, seed=trial)
    val = val[0] if isinstance(val, tuple) else val
    worst = min(worst, float(val))
print('NUMERIC Lemma 1 boundary pc=a^2, j=d^2/(a+c): 12 samples (seed 20261125), min over SO(4) of Phi >= %.3e' % worst)

print('SUMMARY: %d EXACT check(s) failed' % len(FAIL) if FAIL else 'SUMMARY: all EXACT checks passed')
sys.exit(1 if FAIL else 0)
