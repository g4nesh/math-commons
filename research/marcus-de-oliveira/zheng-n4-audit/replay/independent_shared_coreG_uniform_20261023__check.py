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
# independent_shared_coreG_uniform_20261023: seven-contact set G = {0,1,2,7,9,10,23}
# Usage: python -I -B this.py [--numeric]
# ---------------------------------------------------------------------------------------
from fractions import Fraction as Fr
from math import factorial
T = {0, 1, 2, 7, 9, 10, 23}
check('T equals Table 5 auxiliary seven-set no. 10 of the paper', L.AUX7[10] == frozenset(T))
basis = [2, 3, 11, 12, 15, 21, 26]                 # coefficients (a,t,u,b,c,v,d)
names = ['a', 't', 'u', 'b', 'c', 'v', 'd']
nm = kernel_basis_checks(T, basis, [15, 14, 17, 8, 4, 16, None])
check('generators(T) = [2,3,4,11,12,15,21,26], only f26 non-minor', L.generators(T) == [2, 3, 4, 11, 12, 15, 21, 26] and nm == [26])
V = lambda j: L.F_VALUES[j]
x = sympy.symbols('a t u b c v d'); a_, t_, u_, b_, c_, v_, d_ = x
gv = [sum(cf * V(j)[s] for cf, j in zip(x, basis)) for s in range(24)]
check('coefficients: (a,t,u,b,c,v) = (g15,g14,g17,g8,g4,g16), d = g3-g16',
      all(sympy.expand(e) == 0 for e in [gv[15]-a_, gv[14]-t_, gv[17]-u_, gv[8]-b_, gv[4]-c_, gv[16]-v_, gv[3]-gv[16]-d_]))
check('(3): g21=a+d, g20=t+d, g11=u+d, g3=v+d', all(sympy.expand(e) == 0 for e in [gv[21]-a_-d_, gv[20]-t_-d_, gv[11]-u_-d_, gv[3]-v_-d_]))
check('(3): f4+f26 = f2+f3+f11+f21 in all 24 values', all(V(4)[s]+V(26)[s] == V(2)[s]+V(3)[s]+V(11)[s]+V(21)[s] for s in range(24)))
re_ = [(a_+d_)*V(2)[s] + (t_+d_)*V(3)[s] + (u_+d_)*V(11)[s] + b_*V(12)[s] + c_*V(15)[s] + (v_+d_)*V(21)[s] - d_*V(4)[s] for s in range(24)]
check('d<=0 rewrite equals g in all 24 values; profiles 2,3,4,11,12,15,21 are minor',
      all(sympy.expand(re_[s]-gv[s]) == 0 for s in range(24)) and all(L.F_IS_MINOR[j] for j in [2,3,4,11,12,15,21]))
# cited atoms: group maps from the reference families
check('SOC reference (f1,f9;f22) maps to each pair of {f2,f12,f15} with f26', all(group_maps((1,9,22),(i,j,26),True) is not None for i,j in [(2,12),(2,15),(12,15)]))
check('old-H reference (f3,f15;f26) maps (aux,main,nonminor) to (f3,f15;f26) and (f11,f12;f26)', group_maps((3,15,26),(3,15,26)) is not None and group_maps((3,15,26),(11,12,26)) is not None)
# cubic and radial polynomial
t0 = time.time()
ctx, w, Cw = lift_and_cubic(basis, names)
a, t, u, b, c, v, d, z = ctx.gens()
C_doc = (d+t)*(d+u)*(d+v) + a*b*c
Q_doc = (d+t)*(d+u) + a*b + a*c + b*c - c*t - b*u
check('C(w(z)) = C + Q z with C, Q as printed', Cw == C_doc + Q_doc*z)
det, ndiv = hessian_det(ctx, w)
C0, Q, R9, dz, ok = radial_polynomial(ctx, Cw, det)
check('(7): Q^4 det H(w(-C/Q)) = -4 Q^2 R9^2 with R9 a polynomial', ok, f'Bareiss exact divisions {ndiv}, deg_z det H = {dz}, {time.time()-t0:.0f}s')
check('R9: 536 terms, homogeneous of degree 9', len(R9) == 536 and all(sum(m[:7]) == 9 for m in R9.monoms()), f'terms {len(R9)}, degrees {R9.degrees()[:7]}')
check('(8): R9(a,t,u,b,c,v,d) = R9(a,u,t,c,b,v,d)', poly_perm(ctx, R9, [0,2,1,4,3,5,6,7]) == R9)
K = (u+b+d)*(d*d-a*b) + u*b*(b+d) + v*(b+d)*(u+d)
Rc0 = ctx.from_dict({m: cf for m, cf in R9.to_dict().items() if m[4] == 0})
check('R9|c=0 = (u+d)(a+d)(v+d)^2(t+d)^2 K', Rc0 == (u+d)*(a+d)*(v+d)**2*(t+d)**2*K)
R = to7(R9)
# bounded sector: polytope (10) in (a,b,c,l,n)
Hs = [((1,0,0,0,0),0), ((-1,0,0,0,0),1), ((0,-1,0,0,0),1), ((0,0,-1,0,0),1), ((-1,-1,-1,0,0),2),
      ((0,0,0,1,0),0), ((0,0,0,0,1),0), ((0,1,0,-1,0),0), ((0,0,1,0,-1),0), ((0,0,0,-1,-1),1)]
def ev(h, p): return sum(wi*xi for wi, xi in zip(h[0], p)) + h[1]
def rk(rows):
    M = [list(map(Fr, r)) for r in rows]; k = 0
    for col in range(5):
        piv = next((i for i in range(k, len(M)) if M[i][col] != 0), None)
        if piv is None: continue
        M[k], M[piv] = M[piv], M[k]
        for i in range(len(M)):
            if i != k and M[i][col] != 0:
                f = M[i][col] / M[k][col]; M[i] = [p - f*q for p, q in zip(M[i], M[k])]
        k += 1
    return k
def solve5(S):
    M = [list(map(Fr, Hs[i][0])) + [Fr(-Hs[i][1])] for i in S]
    if rk([r[:5] for r in M]) < 5: return None
    for col in range(5):
        piv = next(i for i in range(col, 5) if M[i][col] != 0)
        M[col], M[piv] = M[piv], M[col]
        for i in range(5):
            if i != col and M[i][col] != 0:
                f = M[i][col] / M[col][col]; M[i] = [p - f*q for p, q in zip(M[i], M[col])]
    return tuple(M[i][5] / M[i][i] for i in range(5))
cands = [solve5(S) for S in itertools.combinations(range(10), 5)]
verts = sorted({p for p in cands if p is not None and all(ev(h, p) >= 0 for h in Hs)})
check('polytope (10): 252 hyperplane 5-sets, exactly 13 vertices, all integer', len(cands) == 252 and len(verts) == 13 and all(q.denominator == 1 for p in verts for q in p))
def affrank(vs):
    vs = list(vs)
    return 0 if len(vs) <= 1 else rk([[verts[k][i]-verts[vs[0]][i] for i in range(5)] for k in vs[1:]])
def triangulate(face, dim, pick):
    face = sorted(face)
    if len(face) == dim + 1: return [tuple(face)]
    apex = pick(face)
    tight = {i for i in range(10) if all(ev(Hs[i], verts[k]) == 0 for k in face)}
    facets = set()
    for i in range(10):
        if i in tight: continue
        F = frozenset(k for k in face if ev(Hs[i], verts[k]) == 0)
        if apex not in F and affrank(F) == dim - 1: facets.add(F)
    return [(apex,) + s for F in facets for s in triangulate(F, dim - 1, pick)]
def vol(s):
    M = [[Fr(verts[k][i]-verts[s[0]][i]) for i in range(5)] for k in s[1:]]; dt = Fr(1)
    for col in range(5):
        piv = next((i for i in range(col, 5) if M[i][col] != 0), None)
        if piv is None: return Fr(0)
        if piv != col: M[col], M[piv] = M[piv], M[col]; dt = -dt
        dt *= M[col][col]
        for i in range(col+1, 5):
            f = M[i][col]/M[col][col]; M[i] = [p - f*q for p, q in zip(M[i], M[col])]
    return abs(dt) / factorial(5)
simp = triangulate(range(13), 5, min)
simp2 = triangulate(range(13), 5, max)
V1, V2 = sum(vol(s) for s in simp), sum(vol(s) for s in simp2)
check('pulling triangulation (apex = first vertex): 13 full-dimensional simplices; volume equals that of a second pulling triangulation (apex = last vertex)',
      len(simp) == 13 and all(vol(s) > 0 for s in simp) and V1 == V2, f'volume {V1}, second triangulation has {len(simp2)} simplices')
zc = flint.fmpz_mpoly_ctx.get(('z0','z1','z2','z3','z4','z5','V'), 'deglex')
Z = zc.gens(); one = zc.from_dict({(0,)*7: 1})
def form(cf):
    e = zc.from_dict({})
    for i, x_ in enumerate(cf):
        if x_: e = e + int(x_) * Z[i]
    return e
sizes = []; allok = True
for s in simp:
    vals = [verts[k] for k in s]
    Af, Bf, Cf, Lf, Nf = [form([vals[i][j] for i in range(6)]) for j in range(5)]
    Df = form([1]*6)
    p = cleared(R, [(Af, one), (Nf*Nf, Cf-Nf), (Lf*Lf, Bf-Lf), (Bf, one), (Cf, one), (Z[6]*Df, one), (Df, one)], zc)
    sizes.append(len(p))
    pos = all(q > 0 for q in p.coeffs())
    f4 = [[2 - sum(x_[:3]) for x_ in vals], [1 - x_[3] - x_[4] for x_ in vals], [x_[1]-x_[3] for x_ in vals], [x_[2]-x_[4] for x_ in vals]]
    supp = [frozenset(i for i in range(6) if f[i] > 0) for f in f4]
    v0 = [frozenset(i for i in range(6) if m[i] > 0) for m, q in p.terms() if m[6] == 0 and q > 0]
    strict = all(any(m <= frozenset(S) for m in v0) for r in range(1, 7) for S in itertools.combinations(range(6), r) if all(frozenset(S) & sp for sp in supp))
    allok &= pos and strict
check('(11): on each of the 13 simplices (B-L)^2(C-N)^2 R9(list) has positive coefficients and the strict-support condition holds', allok)
check('(11): table sizes equal the document (as a multiset) and total 352583', sorted(sizes) == sorted([24096, 31304, 24200, 28696, 28696, 24096, 31304, 24200, 32332, 29432, 29432, 28823, 15972]) and sum(sizes) == 352583, f'sizes {sizes}')
pc = flint.fmpz_mpoly_ctx.get(('A','L','N','V','X','Y'), 'deglex')
A_, L_, N_, V_, X_, Y_ = pc.gens(); one = pc.from_dict({(0,)*6: 1})
a0 = one + A_; D = one + L_ + X_ + N_ + Y_
p = cleared(R, [(a0*a0*D, one), (N_*N_, Y_), (L_*L_, X_), (L_+X_, one), (N_+Y_, one), (V_*a0*D, one), (a0*D, one)], pc)
cert('(12) dominant a > 1: X^2 Y^2 R9(list)', p, (0,0,0,0,2,2), 2, 56416)
pc2 = flint.fmpz_mpoly_ctx.get(('B','K','A','N','V','Y'), 'deglex')
B_, K_, A_, N_, V_, Y_ = pc2.gens(); one = pc2.from_dict({(0,)*6: 1})
b0 = one + B_; D = one + A_ + N_ + Y_
Vn = b0*D - N_; VD = (one+K_)*(b0*b0*D + K_*(b0*B_*D + N_))
check('(13): V_n and V_D have positive coefficients in (B,K,A,N,Y)', all(q > 0 for q in Vn.coeffs()) and all(q > 0 for q in VD.coeffs()))
p = cleared(R, [(A_, one), (N_*N_, Y_), (K_*K_*Vn*Vn, VD), (b0*b0*D, one), (N_+Y_, one), (V_*b0*D, one), (b0*D, one)], pc2)
cert('(13) dominant b > 1: Y^2 V_D^2 R9(list)', p, (0,0,0,0,0,2), 2, 141628)
# inverse maps (symbolic)
sa, sb, sc_, sl, sn, X, m = sympy.symbols('a b c l n X m', positive=True)
lam = (sympy.sqrt(X**2 + 4*X*m) - X) / 2
check('cap(X,m) solves X = lam^2/(m-lam)', sympy.simplify(lam**2/(m-lam) - X) == 0)
th = 1 - sa*(sb+sc_)
A1, L1, X1, N1, Y1 = sa-1, sa*sl/th, sa*(sb-sl)/th, sa*sn/th, sa*(sc_-sn)/th
a0s = 1 + A1; Ds = 1 + L1 + X1 + N1 + Y1
lst = [a0s**2*Ds, N1**2/Y1, L1**2/X1, L1+X1, N1+Y1, a0s*Ds]
tgt = [sa, sn**2/(sc_-sn), sl**2/(sb-sl), sb, sc_, 1]
check('(12) inverse map recovers (a, t, u, b, c, d=1)', all(sympy.simplify(p_/(a0s*Ds) - q_) == 0 for p_, q_ in zip(lst, tgt)))
th = 1 - sb*(sa+sc_)
B1, A1, N1, Y1, K1 = sb-1, sb*sa/th, sb*sn/th, sb*(sc_-sn)/th, sl/(1-sl-sn)
b0s = 1 + B1; Ds = 1 + A1 + N1 + Y1
Vn1 = b0s*Ds - N1; VD1 = (1+K1)*(b0s**2*Ds + K1*(b0s*(b0s-1)*Ds + N1))
lst = [A1, N1**2/Y1, K1**2*Vn1**2/VD1, b0s**2*Ds, N1+Y1, b0s*Ds]
check('(13) inverse map recovers (a, t, u, b, c, d=1)', all(sympy.simplify(p_/(b0s*Ds) - q_) == 0 for p_, q_ in zip(lst, tgt)))
if '--numeric' in sys.argv:
    import numpy as np
    rng = np.random.default_rng(20261125)
    def soc(i, j):
        def gen():
            x_ = float(rng.uniform(0.05, 5)); return [(x_, i), (1/x_, j), (1.0, 26)]
        return gen
    def oldh(aux, main):
        def gen():
            m_ = float(rng.uniform(1.01, 5)); return [(1/(m_-1), aux), (m_, main), (1.0, 26)]
        return gen
    for i, j in [(2,12),(2,15),(12,15)]:
        print(f'NUMERIC atom SOC f{i},f{j};f26 at xy=d^2: N=20, restarts 8, worst min Phi = {numeric_atom(rng, soc(i,j)):.3e}')
    for aux, main in [(3,15),(11,12)]:
        print(f'NUMERIC atom old-H f{aux},f{main};f26 at X(m-d)=d^2: N=20, restarts 8, worst min Phi = {numeric_atom(rng, oldh(aux,main)):.3e}')
    found, worst = numeric_stratum_test(T, range(12))
    print(f'NUMERIC theorem: {found} floating points with contacts containing T (lib random_point_on_contact_stratum, seeds 0..11), worst min Phi/max gap = {worst:.3e}')
finish()
