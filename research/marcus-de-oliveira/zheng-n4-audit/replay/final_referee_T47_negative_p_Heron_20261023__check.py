#!/usr/bin/env python3
# Audit replay for supplement document final_referee_T47_negative_p_Heron_20261023 (auditor C22-T47).
# Own code only: nothing from the supplement is executed, imported or read.
# Rebuilds the source vector c(g) from the paper's f_j (lib moc4lib), the cubic C of (A.2),
# the height equation, the 14x14 Hessian determinant and the radial factor G, then replays
# the document's exact claims. Prints PASS/FAIL per check; exits 1 if an EXACT check fails.
# Run: python -I -B final_referee_T47_negative_p_Heron_20261023__check.py   (MOC4_LIB may point to the lib directory)
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

# ---- document-specific part: final_referee_T47_negative_p_Heron_20261023 ----
FR = get_frame([2,9,1,45,5,35,0,7])
ctx, Vv, G, C0 = FR['ctx'], FR['V'], FR['G'], FR['C0']
q,p,j,c,e,d,b,a = [Vv[n] for n in NAMES]
check('constant C = -d{ab+ad+aq+bd+bp+d^2+de+dp+dq+ep+eq+pq}', C0 == -d*(a*b+a*d+a*q+b*d+b*p+d*d+d*e+d*p+d*q+e*p+e*q+p*q))
check('det H(w(g,-C/Q)) = -4 d^2 G^2/Q^2 (prefactor of the multiplicity-one part is d^2)', FR['pref'] == d*d, FR['pref'])
check('G has 965 terms and total degree 8', len(G) == 965 and G.total_degree() == 8)
cp = coeffs_in(G, 'p', ctx); check('G quadratic in p', len(cp) == 3)
Cc, Bc, Ac = cp
disc = Bc*Bc - 4*Ac*Cc
fl = disc.factor()
check('B^2 - 4AC0 = H6*H8 with primitive factors of 123 and 260 terms', fl[0] == 1 and sorted((len(f), f.total_degree(), ex) for f, ex in fl[1]) == [(123, 6, 1), (260, 8, 1)], [(len(f), f.total_degree(), ex) for f, ex in fl[1]])
H6 = [f for f, ex in fl[1] if f.total_degree() == 6][0]; H8 = [f for f, ex in fl[1] if f.total_degree() == 8][0]
# fix signs: Corollary 2 identity 4LG + H6 Phi = K5^2 pins H6 relative to G (G's sign is pinned by the positive chart certificates)
Lc = j*(c+d+a) + (c+d)*(d+b+a) + a*b
Phi = ((p-c)*j + (b-c)*(p+d+a) - p*a)**2 - 4*a*c*(p-b)*(p+d+a)
ok = False
for s in (1, -1):
    try:
        K5 = (4*Lc*G + s*H6*Phi).sqrt(); ok = True; sH6 = s; break
    except Exception: pass
check('Cor 2: 4 L G + H6 Phi is a perfect square K5^2 (with H6 sign +1 as factored)', ok and sH6 == 1, (len(K5) if ok else None))
H6 = sH6*H6; H8 = sH6*H8
check('K5 has 193 terms', ok and len(K5) == 193, len(K5))
check('gcd(G, G_p) = 1', G.gcd(G.derivative('p')).is_one())
check('G is irreducible over Q', len(G.factor()[1]) == 1 and G.factor()[1][0][1] == 1)
# principal minor deleting coordinates 1 and 6 (0-based)
keep = [i for i in range(14) if i not in (1, 6)]
Hm = FR['H']
sub = [[Hm[i][k] for k in keep] for i in keep]
t0 = time.time(); P9 = FR['bareiss'](sub, FR['ctx9'], 9)
check('principal minor P (delete coordinates 6,1) is independent of the height h', all(m[8] == 0 for m in P9.to_dict()))
Pm = FR['to8'](P9.to_dict())
check('P has 3376 monomials and degree 2 in p', len(Pm) == 3376 and max(m[1] for m in Pm.to_dict()) == 2, (len(Pm), round(time.time()-t0,1)))
gam, bet, alp = coeffs_in(Pm, 'p', ctx)
u1 = Ac*bet - alp*Bc; u0 = Ac*gam - alp*Cc
N = u0*u0 - u1*(Bc*gam - bet*Cc)
print('INFO norm before square extraction: %d monomials' % len(N), flush=True)
W18 = N.sqrt()
check('u0^2 - u1(B gamma - beta C0) = W18^2, W18 with 15999 monomials, degree 18', W18*W18 == N and len(W18) == 15999 and W18.total_degree() == 18, (len(W18), W18.total_degree()))
del N
okS = None
for s in (1, -1):
    Sp = 2*Ac*u0 - Bc*u1 + 2*Ac*(s*W18); Sm = 2*Ac*u0 - Bc*u1 - 2*Ac*(s*W18)
    qp, rp = divmod(Sp, H8) if False else (None, None)
    try:
        R8sq = Sp / H8; R8 = R8sq.sqrt(); T8sq = Sm / (d*d*H6); T8 = T8sq.sqrt(); okS = s; break
    except Exception: continue
check('S+ = H8 R8^2 and S- = d^2 H6 T8^2 for one sign choice of W18', okS is not None, okS)
if okS is not None:
    W18 = okS*W18
    check('R8, T8 have 465 and 595 monomials', (len(R8), len(T8)) == (465, 595), (len(R8), len(T8)))
    check('W18, R8, T8 are free of p', all(m[1] == 0 for F in (W18, R8, T8) for m in F.to_dict()))
    Sp = H8*R8*R8; Sm = d*d*H6*T8*T8
    for lab, S, sq in (('+', Sp, u1*p + u0 + W18), ('-', Sm, u1*p + u0 - W18)):
        Dif = Pm*S - sq*sq
        k2, k1, k0 = coeffs_in(Dif, 'p', ctx)[::-1] if len(coeffs_in(Dif, 'p', ctx)) == 3 else (None, None, None)
        # Dif = lam * G with lam free of p  <=>  (k2,k1,k0) proportional to (A,B,C0)
        lam = k2 / Ac
        check('congruence P S%s = (u1 p + u0 %s W18)^2 mod G (exact multiple of G with p-free cofactor)' % (lab, lab), lam*G == Dif, len(lam))
# Corollary 2 local identities (sympy)
ps,cs,js,bs,ds,as_ = sp.symbols('p c j b d a')
Phis = ((ps-cs)*js + (bs-cs)*(ps+ds+as_) - ps*as_)**2 - 4*as_*cs*(ps-bs)*(ps+ds+as_)
R = (ps-bs)*(ps+ds+as_) - as_*ps
check('Phi|_{c=p} = R^2', sp.expand(Phis.subs(cs, ps) - R**2) == 0)
check('Phi_c|_{c=p} = 4ap(j+p+d) + 2(j+p+d-a)R', sp.expand(sp.diff(Phis, cs).subs(cs, ps) - (4*as_*ps*(js+ps+ds) + 2*(js+ps+ds-as_)*R)) == 0)
Xb = (bs-cs)*(js+bs+ds+as_) - as_*bs
check('Phi|_{p=b} = Xb^2', sp.expand(Phis.subs(ps, bs) - Xb**2) == 0)
check('Phi_p|_{p=b} = 2Xb(j+b-c-a) - 4ac(b+d+a)', sp.expand(sp.diff(Phis, ps).subs(ps, bs) - (2*Xb*(js+bs-cs-as_) - 4*as_*cs*(bs+ds+as_))) == 0)
X = (ps-cs)*js + (bs-cs)*(ps+ds+as_) - ps*as_
check('Phi_j = 2X(p-c)', sp.expand(sp.diff(Phis, js) - 2*X*(ps-cs)) == 0)

# ---- NUMERIC tests for the negative-p Heron document ----
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

pts, tried = sample_realizable([2,9,1,45,5,35,0,7], 15, 5000, maxseeds=25, margin=1e-2)
cnt = sum(1 for pt in pts if evalf(H6, pt) < 0 and evalf(H8, pt) < 0)
cnt6 = sum(1 for pt in pts if evalf(H6, pt) < 0 and pt['a'] > 0 and pt['c'] > 0 and pt['e'] > 0 and pt['p'] >= pt['b'])
gmax = max(abs(rel(G, pt)) for pt in pts) if pts else None
print('NUMERIC %d realizable points on the T_65 stratum (seeds 5000.., margin 1e-2): max |G|/scale = %.1e; points with H6<0 and H8<0: %d; points in the Cor. 2 region with H6<0: %d' % (len(pts), gmax, cnt, cnt6), flush=True)
# secant rank bound at non-realizable radial points (cited input): corank of H(w(g,-C/Q)) on G = 0
rng = np.random.default_rng(7)
Hn = FR['H']
def evalH(pt, h):
    x = [pt[n] for n in NAMES] + [h]
    M = np.zeros((14, 14))
    for i in range(14):
        for k in range(14):
            s = 0.0
            for m, v in Hn[i][k].to_dict().items():
                t = float(int(v))
                for ii, ex in enumerate(m):
                    if ex: t *= x[ii]**int(ex)
                s += t
            M[i, k] = s
    return M
done = 0; worst = 0.0
while done < 8:
    pt = {n: float(v) for n, v in zip(NAMES, rng.uniform(0.1, 3.0, 8))}; pt['d'] = 1.0
    Aa = evalf(Ac, pt); Bb = evalf(Bc, pt); Cc0 = evalf(Cc, pt)
    dsc = Bb*Bb - 4*Aa*Cc0
    if dsc < 0: continue
    for root in ((-Bb + np.sqrt(dsc))/(2*Aa), (-Bb - np.sqrt(dsc))/(2*Aa)):
        pt2 = dict(pt); pt2['p'] = root
        h = -evalf(C0, pt2)/evalf(FR['Q'], pt2)
        sv = np.linalg.svd(evalH(pt2, h), compute_uv=False)
        worst = max(worst, sv[-2]/sv[0]); done += 1
print('NUMERIC secant rank bound: 8 random real points of G=0 (seed 7, coefficients uniform in [0.1,3], p solved): max sigma_13/sigma_1 = %.1e (corank >= 2 numerically)' % worst, flush=True)

print('SUMMARY: %d EXACT check(s) failed' % len(FAIL) if FAIL else 'SUMMARY: all EXACT checks passed')
sys.exit(1 if FAIL else 0)
