#!/usr/bin/env python3
# Audit replay for supplement document independent_moc4_T47_IIIb_plus_unswapped_low_j_20261022 (auditor C22-T47).
# Own code only: nothing from the supplement is executed, imported or read.
# Rebuilds the source vector c(g) from the paper's f_j (lib moc4lib), the cubic C of (A.2),
# the height equation, the 14x14 Hessian determinant and the radial factor G, then replays
# the document's exact claims. Prints PASS/FAIL per check; exits 1 if an EXACT check fails.
# Run: python -I -B independent_moc4_T47_IIIb_plus_unswapped_low_j_20261022__check.py   (MOC4_LIB may point to the lib directory)
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

# ---- document-specific part: independent_moc4_T47_IIIb_plus_unswapped_low_j_20261022 ----
FR = get_frame([2,9,1,45,5,35,0,7])
ctx, Vv, G = FR['ctx'], FR['V'], FR['G']
q,p,j,c,e,d,b,a = [Vv[n] for n in NAMES]
check('G has 965 terms and total degree 8', len(G) == 965 and G.total_degree() == 8, (len(G), G.total_degree()))
cc = coeffs_in(G, 'c', ctx)
check('G is quadratic in c', len(cc) == 3)
gam, bet, alp = cc
Phi = ((p-c)*j + (b-c)*(p+d+a) - p*a)**2 - 4*a*c*(p-b)*(p+d+a)
pc_ = coeffs_in(Phi, 'c', ctx); K_, J_, H_ = pc_
check('H = (p+j+d+a)^2', H_ == (p+j+d+a)**2)
check('K = [pj + b(p+d+a) - pa]^2', K_ == (p*j + b*(p+d+a) - p*a)**2)
U0 = H_*gam - alp*K_; U1 = H_*bet - alp*J_; V0 = K_*bet - gam*J_
check('identity H G = alpha Phi + U0 + c U1', H_*G == alp*Phi + U0 + c*U1)
check('identity K G = gamma Phi + c (V0 - c U0)', K_*G == gam*Phi + c*(V0 - c*U0))
# chart (d = 1)
cctx = flint.fmpz_mpoly_ctx.get(('B','V','P','T','W','U','Q'), 'deglex')
B_,V_,P_,T_,W_,U_,Q_ = cctx.gens(); one = cctx.from_dict({(0,)*7: 1})
bb = 1 + B_; Z = B_*T_ + bb*W_
subs = {'b': (bb, one), 'a': (V_, bb*(1+V_)), 'p': (P_, one), 'q': (Q_, one), 'j': (T_, T_+W_),
        'e': (T_*(Z*(1+U_) + T_*U_), (T_+W_)*Z*(1+U_)), 'd': (one, one), 'c': (one, one)}
# EJ = j U0 + U1 : treat as polynomial in (q,p,j,e,d,b,a); c absent from all four
EJ = j*U0 + U1
res = {}
for nm, F in [('alpha', alp), ('gamma', gam), ('V0', V0), ('EJ', EJ)]:
    Ft = chartsub(F, subs, cctx)
    n, neg, mn = coeff_stats(Ft)
    mk = has_marker(Ft, ['W','U','Q'], cctx)
    res[nm] = (n, neg, mk)
    check('low-j chart: cleared %s has only positive coefficients' % nm, neg == 0 and n > 0, 'terms=%d negative=%d markers(W=U=Q=0)=%d' % (n, neg, mk))
check('gamma and EJ have a positive monomial free of W,U,Q', res['gamma'][2] > 0 and res['EJ'][2] > 0)
# direct symbolic check of the e-map with sympy
Bs,Ts,Ws,Us = sp.symbols('B T W U', positive=True)
bs = 1+Bs; js = Ts/(Ts+Ws); Zs = Bs*Ts + bs*Ws
es = Ts*(Zs*(1+Us) + Ts*Us)/((Ts+Ws)*Zs*(1+Us))
check('e-map equals j + j^2 U/((b-j)(1+U)) (sympy)', sp.simplify(es - (js + js**2*Us/((bs-js)*(1+Us)))) == 0)
check('e-map upper limit U->oo equals bj/(b-j)', sp.simplify(sp.limit(es, Us, sp.oo) - bs*js/(bs-js)) == 0)

print('SUMMARY: %d EXACT check(s) failed' % len(FAIL) if FAIL else 'SUMMARY: all EXACT checks passed')
sys.exit(1 if FAIL else 0)
