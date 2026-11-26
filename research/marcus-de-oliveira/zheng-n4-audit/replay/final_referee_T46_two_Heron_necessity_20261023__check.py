#!/usr/bin/env python3
# Audit replay for supplement document final_referee_T46_two_Heron_necessity_20261023 (auditor C22-T47).
# Own code only: nothing from the supplement is executed, imported or read.
# Rebuilds the source vector c(g) from the paper's f_j (lib moc4lib), the cubic C of (A.2),
# the height equation, the 14x14 Hessian determinant and the radial factor G, then replays
# the document's exact claims. Prints PASS/FAIL per check; exits 1 if an EXACT check fails.
# Run: python -I -B final_referee_T46_two_Heron_necessity_20261023__check.py   (MOC4_LIB may point to the lib directory)
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

# ---- document-specific part: final_referee_T46_two_Heron_necessity_20261023 ----
GENS = [27,1,59,5,2,28,6,4]
T = frozenset(range(24))
for jg in GENS: T &= L.F_ZERO[jg]
check('the eight sources vanish exactly on the spectral tuple {1,3,5,6,7,8} = paper T_47 (literal)', sorted(T) == [1,3,5,6,7,8] and frozenset(T) == L.CONTACT_SETS[47])
gT = L.generators(T)
check('INFO: generators(T_47) has 10 elements; the document uses 8 of them (omits f0 = n1,4 and f3 = n14,34)', sorted(set(gT) - set(GENS)) == [0, 3], (gT))
check('the 8 sources are linearly independent', L.rank([L.F_VECS[jg] for jg in GENS]) == 8)
FR = get_frame(GENS)
ctx, Vv, G, C0, Qh = FR['ctx'], FR['V'], FR['G'], FR['C0'], FR['Q']
q,p,j,c,e,d,b,a = [Vv[n] for n in NAMES]
D_ = a+d; M_ = (c+d)*(a+d+e); L_ = M_ + b*D_
check('height constant C = -(q+d)L - jM', C0 == -(q+d)*L_ - j*M_, C0)
print('INFO Hessian prefactor (det H(w) = -4 * pref * G^2/Q^2): pref =', FR['pref'], '; G terms', len(G), 'degree', G.total_degree(), flush=True)
cj = coeffs_in(G, 'j', ctx); check('G is quadratic in j', len(cj) == 3)
Cc, Bc, Ac = cj
disc = Bc*Bc - 4*Ac*Cc
fl = disc.factor()
print('INFO disc_j G factorization:', fl[0], [(len(f), f.total_degree(), ex) for f, ex in fl[1]], flush=True)
qd = [f for f, ex in fl[1] if f == q + d]
others = [(f, ex) for f, ex in fl[1] if f != q + d]
check('disc_j G = (q+d)^2 * H6 * H8 (one sextic, one octic, content +-1)', len(qd) == 1 and dict((f.total_degree(), ex) for f, ex in others) == {6: 1, 8: 1} and abs(int(fl[0])) == 1)
H6 = [f for f, ex in others if f.total_degree() == 6][0]; H8 = [f for f, ex in others if f.total_degree() == 8][0]
T0 = L_*p + (b*D_ - M_)*q - b*D_*(c+e)
rhs = T0*T0 + 4*b*D_*M_*(q-c)*(q-e)
sg = 1 if H6 == rhs else (-1 if H6 == -rhs else 0)
check('(1): H6 = T0^2 + 4bDM(q-c)(q-e) (up to the sign normalisation of the factor)', sg != 0, sg)
H6 = sg*H6 if sg else H6; H8 = (sg*int(fl[0]))*H8 if sg else H8
check('disc_j G = (q+d)^2 H6 H8 with these signs', disc == (q+d)**2*H6*H8)
check('(H6)_p = 2 L T0', H6.derivative('p') == 2*L_*T0)
check('(H6)_q = 2(bD-M)T0 + 4bDM(2q-c-e)', H6.derivative('q') == 2*(b*D_-M_)*T0 + 4*b*D_*M_*(2*q-c-e))
# exceptional locus p = q = c = e = t
tctx = flint.fmpz_mpoly_ctx.get(('t','j','d','b','a'), 'deglex')
t_, j_, d_, b_, a_ = tctx.gens(); one = tctx.from_dict({(0,)*5: 1})
sub = {'q': (t_, one), 'p': (t_, one), 'c': (t_, one), 'e': (t_, one), 'j': (j_, one), 'd': (d_, one), 'b': (b_, one), 'a': (a_, one)}
H8t = chartsub(H8, sub, tctx)
Dt = a_+d_; Mt = (t_+d_)*(a_+d_+t_); Lt = Mt + b_*Dt
rest = H8t - 4*d_*Dt*b_*b_*Lt*t_*(t_+d_+b_)
try:
    W0 = rest.sqrt(); okW = True
except Exception:
    okW = False
check('on p=q=c=e=t: H8 - 4dDb^2 L t(t+d+b) is a perfect square W0^2', okW)
if okW:
    W0_0 = W0.subs({'t': 0}) if hasattr(W0, 'subs') else None
    check('W0 at t = 0 equals +-bdL', W0_0 is not None and (W0_0 == b_*d_*Lt.subs({'t': 0}) or W0_0 == -b_*d_*Lt.subs({'t': 0})), W0_0)
check('gcd(G, G_j) = 1', G.gcd(G.derivative('j')).is_one())
fG = G.factor(); check('G irreducible over Q', len(fG[1]) == 1 and fG[1][0][1] == 1)

print('SUMMARY: %d EXACT check(s) failed' % len(FAIL) if FAIL else 'SUMMARY: all EXACT checks passed')
sys.exit(1 if FAIL else 0)
