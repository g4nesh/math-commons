# ---------------------------------------------------------------------------------------
# Shared part (identical in the three K3 batch check scripts). Own code only; nothing from
# the supplement is executed or read. Library: ../lib/moc4lib.py.
# ---------------------------------------------------------------------------------------
import os, sys, time, itertools
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'lib'))
import moc4lib as L
import flint, sympy

FAILS = []
def check(name, ok, info=''):
    print(f"{'PASS' if ok else 'FAIL'}  {name}" + (f"  [{info}]" if info else ''), flush=True)
    if not ok: FAILS.append(name)

def kernel_basis_checks(T, basis, rows):
    """basis = list of f-indices (order of the document's coefficients); rows = permutation indices whose
    values of the basis profiles form the identity (None = not claimed)."""
    check('generators(T) as listed', True, f'generators={L.generators(T)}')
    check('contact rows have rank 7', L.rank([L.P[s] for s in T]) == 7)
    check('every basis profile vanishes on T', all(L.F_VALUES[j][s] == 0 for j in basis for s in T))
    check('basis profiles independent (rank 7)', L.rank([L.F_VECS[j] for j in basis]) == 7)
    if all(r is not None for r in rows):
        M = [[L.F_VALUES[j][r] for j in basis] for r in rows]
        check(f'rows {rows} of the basis form I_7', M == [[int(i == k) for k in range(7)] for i in range(7)])
    nonminor = [j for j in L.generators(T) if not L.F_IS_MINOR[j]]
    return nonminor

def lift_and_cubic(basis, names):
    """w(z) = sum_k x_k f_{basis[k]} + z*kappa (as polynomials); returns ctx, w, C(w)."""
    ctx = flint.fmpz_mpoly_ctx.get(tuple(names) + ('z',), 'deglex')
    g = ctx.gens()
    w = []
    for r in range(14):
        e = ctx.from_dict({})
        for k, j in enumerate(basis):
            cf = L.F_VECS[j][r]
            assert cf == int(cf)
            if cf: e = e + int(cf) * g[k]
        if L.KAPPA[r]: e = e + int(L.KAPPA[r]) * g[7]
        w.append(e)
    return ctx, w, L.cubic_C(w)

def hessian_det(ctx, w):
    ms, Csym = L.cubic_sympy()
    def lin(expr):
        p = sympy.Poly(expr, *ms); out = ctx.from_dict({})
        for mon, cf in p.terms():
            out = out + (int(cf) if sum(mon) == 0 else int(cf) * w[mon.index(1)])
        return out
    M = [[lin(sympy.diff(Csym, ms[i], ms[j])) for j in range(14)] for i in range(14)]
    n = 14; prev = ctx.from_dict({(0,) * 8: 1}); sign = 1; ndiv = 0
    for k in range(n - 1):
        if M[k][k] == 0:
            sw = next(i for i in range(k + 1, n) if M[i][k] != 0)
            M[k], M[sw] = M[sw], M[k]; sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                q, r = divmod(M[i][j] * M[k][k] - M[i][k] * M[k][j], prev)
                assert r == 0; ndiv += 1
                M[i][j] = q
        prev = M[k][k]
    return sign * M[n - 1][n - 1], ndiv

def radial_polynomial(ctx, Cw, det):
    """Returns (C0, Q, R) with C(w(z)) = C0 + Q z and Q^4 det H(w(-C0/Q)) = -4 Q^2 R^2 checked exactly."""
    Q = Cw.derivative(7); C0 = Cw - Q * ctx.gens()[7]
    assert Q.degrees()[7] == 0 and C0.degrees()[7] == 0
    dz = det.degrees()[7]
    parts = {}
    for mon, cf in det.terms():
        parts.setdefault(mon[7], {})[tuple(mon[:7]) + (0,)] = cf
    Ds = [ctx.from_dict(parts.get(i, {})) for i in range(dz + 1)]
    E4 = sum(((-C0) ** i) * (Q ** (4 - i)) * Ds[i] for i in range(dz + 1))   # Q^4 det H(w(-C0/Q))
    q, r = divmod(-E4, 4 * Q * Q)
    Rr = q.sqrt() if r == 0 else None
    ok = r == 0 and Rr is not None and Rr * Rr == q
    return C0, Q, Rr, dz, ok

def to7(P):
    return {tuple(int(e) for e in m[:7]): int(c) for m, c in P.to_dict().items()}

def cleared(R, coords, ctx):
    """prod_i den_i^{deg_i R} * R(num_i/den_i) for coords = [(num, den)] (7 coordinates)."""
    deg = [max(m[i] for m in R) for i in range(7)]
    one = ctx.from_dict({(0,) * ctx.nvars(): 1})
    pw = []
    for i, (n, e) in enumerate(coords):
        npw = [one]; epw = [one]
        for k in range(deg[i]): npw.append(npw[-1] * n); epw.append(epw[-1] * e)
        pw.append([npw[k] * epw[deg[i] - k] for k in range(deg[i] + 1)])
    tot = ctx.from_dict({})
    for m, cf in R.items():
        term = one * cf
        for i in range(7): term = term * pw[i][m[i]]
        tot = tot + term
    return tot

def poly_perm(ctx, P, idx):
    return ctx.from_dict({tuple(m[idx[i]] for i in range(len(idx))): c for m, c in P.to_dict().items()})

def cert(name, poly, marker, coeff, count):
    cs = poly.coeffs()
    allpos = all(c > 0 for c in cs)
    mk = poly.to_dict().get(marker, 0)
    check(f'{name}: all coefficients positive, marker present', allpos and mk == coeff,
          f'terms {len(cs)} (document {count}), marker coefficient {mk} (document {coeff})')
    check(f'{name}: term count equals document', len(cs) == count)

# f-permutation induced by the group (for the "literal actions" of cited atoms)
_FIDX = {L.F_VALUES[j]: j for j in range(82)}
def _fperm(g):
    t = L.g_index_map(g); out = []
    for j in range(82):
        v = [None] * 24
        for s in range(24): v[t[s]] = L.F_VALUES[j][s]
        out.append(_FIDX[tuple(v)])
    return out
_FP = None
def group_maps(ref, tgt, unordered_prefix=False):
    global _FP
    if _FP is None: _FP = [_fperm(g) for g in L.G_ELEMENTS]
    for g, p in zip(L.G_ELEMENTS, _FP):
        img = [p[j] for j in ref]
        if unordered_prefix:
            if sorted(img[:-1]) == sorted(tgt[:-1]) and img[-1] == tgt[-1]: return g
        elif img == list(tgt): return g
    return None

def numeric_stratum_test(T, seeds, restarts=10):
    import numpy as np
    Pn = L.P_numpy(); worst = None; found = 0
    for seed in seeds:
        pt = L.random_point_on_contact_stratum(T, seed=seed, tries=20)
        if pt is None: continue
        gaps = np.asarray(pt['gaps'])
        if gaps.min() < -1e-7 or max(abs(gaps[s]) for s in T) > 1e-7: continue
        found += 1
        cvec = np.linalg.lstsq(Pn, gaps, rcond=None)[0]
        v, _ = L.min_phi_over_SO4(list(cvec), restarts=restarts, seed=seed)
        v = v / max(1.0, float(np.abs(gaps).max()))
        worst = v if worst is None else min(worst, v)
    return found, worst

def numeric_atom(rng, gen, N=20, restarts=8):
    worst = None
    for _ in range(N):
        pairs = gen()
        c = [0.0] * 14
        for cf, j in pairs:
            for r in range(14): c[r] += cf * float(L.F_VECS[j][r])
        v, _ = L.min_phi_over_SO4(c, restarts=restarts, seed=int(rng.integers(1 << 30)))
        worst = v if worst is None else min(worst, v)
    return worst

def finish():
    print('SUMMARY:', 'all EXACT checks passed' if not FAILS else f'{len(FAILS)} EXACT check(s) FAILED: {FAILS}')
    sys.exit(1 if FAILS else 0)
# ---------------------------------------------------------------------------------------
# root_shared_coreH_uniform_20261023: seven-contact set H = {0,1,6,9,16,17,22}
# Usage: python -I -B this.py [--numeric]
# ---------------------------------------------------------------------------------------
T = {0, 1, 6, 9, 16, 17, 22}
check('T equals Table 5 auxiliary seven-set no. 21 of the paper', L.AUX7[21] == frozenset(T))
basis = [15, 3, 2, 5, 29, 33, 31]                  # coefficients (a,b,q,t,u,v,d)
names = ['a', 'b', 'q', 't', 'u', 'v', 'd']
nm = kernel_basis_checks(T, basis, [5, 14, 15, 11, 8, 3, 7])
check('generators(T) = {2,3,5,15,29,31,33} = basis, only f31 non-minor', sorted(L.generators(T)) == sorted(basis) and nm == [31])
check('SOC reference (f1,f9;f22) maps to (f15,f3;f31)', group_maps((1,9,22),(15,3,31),True) is not None)
check('old-H reference (f3,f15;f26) maps (aux,main,nonminor) to (f2,f15;f31),(f29,f15;f31),(f5,f3;f31),(f33,f3;f31)',
      all(group_maps((3,15,26), tg) is not None for tg in [(2,15,31),(29,15,31),(5,3,31),(33,3,31)]))
t0 = time.time()
ctx, w, Cw = lift_and_cubic(basis, names)
a, b, q, t, u, v, d, z = ctx.gens()
xi, psi, Pp, Rr_ = q+u, t+v, q*u, t*v
C_doc = Pp*(a-psi-d) + Rr_*(b-xi-d) - d*(xi+d)*(psi+d)
Q_doc = -a*b + a*xi + b*psi - xi*psi - d*(xi+psi) - d*d
check('(5): C(w(z)) = C + Q z with C, Q as printed', Cw == C_doc + Q_doc*z)
# (6) as a polynomial identity: C = P(a-d) + R(b-d) - P psi - R xi - d xi psi - d^2(xi+psi) - d^3
check('(6) decomposition: C = P(a-d) + R(b-d) - P psi - R xi - d xi psi - d^2 (xi+psi) - d^3',
      C_doc == Pp*(a-d) + Rr_*(b-d) - Pp*psi - Rr_*xi - d*xi*psi - d*d*(xi+psi) - d**3)
det, ndiv = hessian_det(ctx, w)
C0, Q, R, dz, ok = radial_polynomial(ctx, Cw, det)
G, r = divmod(R, d)
check('(7): Q^4 det H(w(-C/Q)) = -4 d^2 Q^2 G^2 with G a polynomial', ok and r == 0, f'Bareiss exact divisions {ndiv}, deg_z det H = {dz}, {time.time()-t0:.0f}s')
check('G: 499 terms, homogeneous of degree 8, coefficient 1 at d^8', len(G) == 499 and all(sum(m[:7]) == 8 for m in G.monoms()) and G.to_dict().get((0,0,0,0,0,0,8,0)) == 1)
check('(8): G symmetric under q<->u, t<->v, (a,b,q,t,u,v)<->(b,a,t,q,v,u)',
      poly_perm(ctx, G, [0,1,4,3,2,5,6,7]) == G and poly_perm(ctx, G, [0,1,2,5,4,3,6,7]) == G and poly_perm(ctx, G, [1,0,3,2,5,4,6,7]) == G)
Gd = to7(G)
pc = flint.fmpz_mpoly_ctx.get(('A','L','Y','K','U','V'), 'deglex')
A_, L_, Y_, K_, U_, V_ = pc.gens(); one = pc.from_dict({(0,)*6: 1}); zero = pc.from_dict({})
cert('(9) a,b<1', cleared(Gd, [(A_, one+A_), (L_, one+L_), (Y_, one), (K_, one), (U_, one), (V_, one), (one, one)], pc), (0,)*6, 1, 709)
cert('(9) wall a=1, b<1', cleared(Gd, [(one, one), (L_, one+L_), (Y_, one), (K_, one), (U_, one), (V_, one), (one, one)], pc), (0,)*6, 2, 235)
aa = one + A_; D = one + L_ + Y_; Vn = aa*D - L_; VD = (one+K_)*(aa*aa*D + K_*(aa*A_*D + L_))
check('(10): V_n and V_D have positive coefficients', all(c_ > 0 for c_ in Vn.coeffs()) and all(c_ > 0 for c_ in VD.coeffs()))
p = cleared(Gd, [(aa*aa*D, one), (L_+Y_, one), (K_*K_*Vn*Vn, VD*(one+U_)), (L_*L_, Y_*(one+V_)), (K_*K_*Vn*Vn*U_, VD*(one+U_)), (L_*L_*V_, Y_*(one+V_)), (aa*D, one)], pc)
cert('(11) dominant a>1, b>0', p, (0,0,4,0,0,0), 2, 476039)
D0 = (one+K_)*(aa + A_*K_)
cert('(13) a>1, b=0', cleared(Gd, [(aa, one), (zero, one), (K_*K_, D0*(one+U_)), (L_, one), (K_*K_*U_, D0*(one+U_)), (V_, one), (one, one)], pc), (0,)*6, 2, 2356)
sa, sb, lam, mu, X, m = sympy.symbols('a b lam mu X m', positive=True)
U, Vs = sympy.symbols('U V', nonnegative=True)
cp = (sympy.sqrt(X**2 + 4*X*m) - X) / 2
check('cap(X,m) solves X = lam^2/(m-lam)', sympy.simplify(cp**2/(m-cp) - X) == 0)
th = 1 - sa*sb; Ds = 1/th; L1 = sa*mu/th; Y1 = sa*(sb-mu)/th; K1 = lam/(1-lam-mu)
Vn1 = sa*Ds - L1; VD1 = (1+K1)*(sa**2*Ds + K1*(sa*(sa-1)*Ds + L1))
lst = [sa**2*Ds, L1+Y1, K1**2*Vn1**2/(VD1*(1+U)), L1**2/(Y1*(1+Vs)), K1**2*Vn1**2*U/(VD1*(1+U)), L1**2*Vs/(Y1*(1+Vs)), sa*Ds]
xi_, psi_ = lam**2/(sa-lam), mu**2/(sb-mu)
tgt = [sa, sb, xi_/(1+U), psi_/(1+Vs), xi_*U/(1+U), psi_*Vs/(1+Vs), 1]
check('(12): D = 1+L+Y and the inverse map recovers all seven normalized coefficients',
      sympy.simplify(Ds - 1 - L1 - Y1) == 0 and all(sympy.simplify(p_/(sa*Ds) - q_) == 0 for p_, q_ in zip(lst, tgt)))
K2 = lam/(1-lam)
check('(13): K^2/D0 = lam^2/(a-lam) with K = lam/(1-lam)', sympy.simplify(K2**2/((1+K2)*(sa+(sa-1)*K2)) - lam**2/(sa-lam)) == 0)
if '--numeric' in sys.argv:
    import numpy as np
    rng = np.random.default_rng(20261126)
    def gen_soc():
        x_ = float(rng.uniform(0.05, 5)); return [(x_, 15), (1/x_, 3), (1.0, 31)]
    print(f'NUMERIC atom SOC f15,f3;f31 at ab=d^2: N=20, restarts 8, worst min Phi = {numeric_atom(rng, gen_soc):.3e}')
    for aux, main in [(2,15),(29,15),(5,3),(33,3)]:
        def gen(aux=aux, main=main):
            m_ = float(rng.uniform(1.01, 5)); return [(1/(m_-1), aux), (m_, main), (1.0, 31)]
        print(f'NUMERIC atom old-H f{aux},f{main};f31 at X(m-d)=d^2: N=20, restarts 8, worst min Phi = {numeric_atom(rng, gen):.3e}')
    found, worst = numeric_stratum_test(T, range(12))
    print(f'NUMERIC theorem: {found} floating points with contacts containing T (seeds 0..11), worst min Phi/max gap = {worst:.3e}')
finish()
