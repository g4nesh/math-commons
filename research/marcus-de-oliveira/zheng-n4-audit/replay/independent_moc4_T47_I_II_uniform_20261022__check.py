#!/usr/bin/env python3
# Audit replay for supplement document independent_moc4_T47_I_II_uniform_20261022 (auditor C22-T47).
# Own code only: nothing from the supplement is executed, imported or read.
# Rebuilds the source vector c(g) from the paper's f_j (lib moc4lib), the cubic C of (A.2),
# the height equation, the 14x14 Hessian determinant and the radial factor G, then replays
# the document's exact claims. Prints PASS/FAIL per check; exits 1 if an EXACT check fails.
# Run: python -I -B independent_moc4_T47_I_II_uniform_20261022__check.py   (MOC4_LIB may point to the lib directory)
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

# ---- document-specific part: independent_moc4_T47_I_II_uniform_20261022 ----
# (a) atom tuples in the proved orbits (so4-atoms.md representatives), roles as printed
vidx = {tuple(v): k for k, v in enumerate(L.F_VALUES)}
fmaps = []
for g in L.G_ELEMENTS:
    fmaps.append([vidx[tuple(int(x) for x in L.g_act_values(g, list(L.F_VALUES[k])))] for k in range(82)])
def in_orbit(rep, tgt, unordered_groups):
    for fm in fmaps:
        img = [[fm[x] for x in grp] for grp in rep]
        if all((sorted(i) == sorted(t)) if u else (i == t) for i, t, u in zip(img, tgt, unordered_groups)):
            return True
    return False
check('SOC (45,1;35) in the orbit of (1,9;22), pair unordered', in_orbit([[1,9],[22]], [[45,1],[35]], [True, False]))
check('old-H (aux 9, main 1; 35) in the orbit of (3,15;26), ordered', in_orbit([[3],[15],[26]], [[9],[1],[35]], [False, False, False]))
check('old-H reversed roles (aux 1, main 9; 35) NOT in that orbit', not in_orbit([[3],[15],[26]], [[1],[9],[35]], [False, False, False]))
check('chain (centre 1; leaves 7,9; 35) in the orbit of (2;0,12;70)', in_orbit([[2],[0,12],[70]], [[1],[7,9],[35]], [False, True, False]))
check('four (0,7;5;35) in the orbit of (0,6;10;54), a,b unordered', in_orbit([[0,6],[10],[54]], [[0,7],[5],[35]], [True, False, False]))
check('T24 core profile map (59,53,55,67,54,2) -> (45,5,7,0,35,1) is realised by one group element',
      any([fm[x] for x in (59,53,55,67,54,2)] == [45,5,7,0,35,1] for fm in fmaps))
# capacity algebra: four-coefficient with (a<-b-role f0 coefficient b, f7 coefficient a, c<-e, d=1)
as_, bs, es = sp.symbols('a b e', positive=True)
check('four atom ab(c-d) >= d^2(a+b+c+d) with (a,b,c,d)=(b,a,e,1) is e >= (a+1)(b+1)/(ab-1) for ab>1',
      sp.simplify(sp.solve(sp.Eq(bs*as_*(es-1), (as_+bs+es+1)), es)[0] - (as_+1)*(bs+1)/(as_*bs-1)) == 0)
FR = get_frame([2,9,1,45,5,35,0,7])
ctx, Vv, G, C0, Qh = FR['ctx'], FR['V'], FR['G'], FR['C0'], FR['Q']
q,p,j,c,e,d,b,a = [Vv[n] for n in NAMES]
check('Q (d=1) = ce - cj - 1 - e - p - q - ep - eq + jp - pq', Qh == c*e - c*j - d*d - d*e - d*p - d*q - e*p - e*q + j*p - p*q)
# (b) sector I, ab <= 1 charts
cctx = flint.fmpz_mpoly_ctx.get(('A','P','T','U','V','S','Q'), 'deglex')
A_,P_,T_,U_,V_,S_,Q_ = cctx.gens(); one = cctx.from_dict({(0,)*7: 1})
for lab, (an, ad), (rn, rd) in [('a = A/(1+A)', (A_, 1+A_), (1+P_, P_)), ('a = 1+A', (1+A_, one), (1+A_+P_, P_*(1+A_)))]:
    # j = r theta(T), e = j theta(U), b = theta(V)/a, c = theta(S)/j
    jn, jd = rn*T_, rd*(1+T_)
    subs = {'a': (an, ad), 'p': (P_, one), 'j': (jn, jd), 'e': (jn*U_, jd*(1+U_)), 'b': (V_*ad, (1+V_)*an),
            'c': (S_*jd, (1+S_)*jn), 'q': (Q_, one), 'd': (one, one)}
    Gt = chartsub(G, subs, cctx); n, neg, mn = coeff_stats(Gt)
    check('sector I, ab<=1, %s: cleared G has only positive coefficients' % lab, neg == 0, 'terms=%d (document: 19536 / 43392)' % n)
# (c) ab = 1 wall, a = A/(1+A): b = 1/a, p = P, j = r theta(T), e = j + U, c = theta(S)/e
an, ad = A_, 1+A_; rn, rd = 1+P_, P_
jn, jd = rn*T_, rd*(1+T_)
en, ed = jn + U_*jd, jd
subs = {'a': (an, ad), 'b': (ad, an), 'p': (P_, one), 'j': (jn, jd), 'e': (en, ed), 'c': (S_*ed, (1+S_)*en), 'q': (Q_, one), 'd': (one, one)}
Gw = chartsub(G, subs, cctx); n, neg, mn = coeff_stats(Gw)
print('INFO ab=1 wall chart: terms=%d negative=%d (document: 9233 terms, 4 negative)' % (n, neg), flush=True)
iU = cctx.names().index('U')
fib = {}
for m, v in Gw.to_dict().items():
    key = tuple(x if t != iU else 0 for t, x in enumerate(m)); fib.setdefault(key, {})[m[iU]] = int(v)
Us = sp.Symbol('U')
bad = [k for k, f in fib.items() if min(f.values()) < 0]
facs = []
for k in bad:
    f = fib[k]; lo = min(f); poly = sp.Poly(sum(v*Us**(e_-lo) for e_, v in f.items()), Us)
    facs.append(poly.as_expr())
print('INFO fibers with a negative coefficient (after removing the power of U):', facs, flush=True)
target = [4*Us**2-Us+13, 8*Us**2-2*Us+6, 5*Us**2-Us+1, Us**3+7*Us**2-Us+2]
check('the negative fibers are, up to a positive integer factor, the four listed univariate polynomials',
      len(facs) == 4 and all(any(sp.simplify(fa/t).is_number and sp.simplify(fa/t) > 0 for t in target) for fa in facs))
check('each listed polynomial is positive on U >= 0 (positive at 0, every real root negative)', all(sp.Poly(t, Us).eval(0) > 0 and all(r < 0 for r in sp.real_roots(sp.Poly(t, Us))) for t in target))

print('SUMMARY: %d EXACT check(s) failed' % len(FAIL) if FAIL else 'SUMMARY: all EXACT checks passed')
sys.exit(1 if FAIL else 0)
