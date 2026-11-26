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
# root_shared_seven_A_uniform_20261023: seven-contact set {0,1,2,7,9,10,17}
# Usage: python -I -B this.py [--numeric]
# ---------------------------------------------------------------------------------------
T = {0, 1, 2, 7, 9, 10, 17}
check('T equals Table 5 auxiliary seven-set no. 9 of the paper', L.AUX7[9] == frozenset(T))
basis = [0, 2, 3, 12, 15, 21, 26]                  # coefficients (q,p,t,b,c,v,d)
names = ['q', 'p', 't', 'b', 'c', 'v', 'd']
nm = kernel_basis_checks(T, basis, [23, 15, 14, 8, 4, 16, 11])
check('generators(T) = basis, only f26 non-minor', sorted(L.generators(T)) == sorted(basis) and nm == [26])
check('SOC reference (f1,f9;f22) maps to each pair of {f2,f12,f15} with f26', all(group_maps((1,9,22),(i,j,26),True) is not None for i,j in [(2,12),(2,15),(12,15)]))
check('old-H reference (f11,f1;f34) and (f3,f15;f26) map (aux,main,nonminor) to (f3,f15;f26)', group_maps((11,1,34),(3,15,26)) is not None and group_maps((3,15,26),(3,15,26)) is not None)
gch = group_maps((2,0,12,70),(15,0,3,26))
check('chain reference (f2,f0,f12;f70) maps to (f15,f0,f3;f26) (main c, auxiliaries q,t)', gch is not None, f'g = {gch}')
t0 = time.time()
ctx, w, Cw = lift_and_cubic(basis, names)
q, p, t, b, c, v, d, z = ctx.gens()
C_doc = d*(d+v)*(d+q+t) + p*b*c
Q_doc = d*d + d*t + b*c + b*p + c*p - c*t
check('C(w(z)) = C + z Q with C, Q as printed', Cw == C_doc + Q_doc*z)
det, ndiv = hessian_det(ctx, w)
C0, Q, G, dz, ok = radial_polynomial(ctx, Cw, det)
check('(4): Q^4 det H(w(-C/Q)) = -4 Q^2 G^2 with G a polynomial', ok, f'Bareiss exact divisions {ndiv}, deg_z det H = {dz}, {time.time()-t0:.0f}s')
check('G: 486 terms, homogeneous of degree 9, coordinate degrees (3,2,2,2,2,3,9)', len(G) == 486 and all(sum(m[:7]) == 9 for m in G.monoms()) and list(G.degrees()[:7]) == [3,2,2,2,2,3,9])
Gc0 = ctx.from_dict({m: cf for m, cf in G.to_dict().items() if m[4] == 0})
check('G|c=0 = d(p+d)(v+d)^2(q+t+d)^2 {qd(b+v+d) + (b+d)(d^2-pb+vd)}', Gc0 == d*(p+d)*(v+d)**2*(q+t+d)**2*(q*d*(b+v+d) + (b+d)*(d*d-p*b+v*d)))
Gd = to7(G)
pc = flint.fmpz_mpoly_ctx.get(('P','Q','R','U','T','V'), 'deglex')
P_, Q_, R_, U_, T_, V_ = pc.gens(); one = pc.from_dict({(0,)*6: 1})
S = (one+P_)*(one+Q_+R_)
mains = {'R0': ((P_, one+P_), (Q_, one+Q_), (R_, one+R_)), 'Rp': ((one+P_, one), (Q_, S), (R_, S)),
         'Rb': ((Q_, S), (one+P_, one), (R_, S)), 'Rc': ((Q_, S), (R_, S), (one+P_, one))}
mono = lambda **kw: tuple(kw.get(x_, 0) for x_ in ['P','Q','R','U','T','V'])
expected = {('R0',0): (1784, mono(R=3), 1), ('R0',1): (3234, mono(R=2,U=2), 1), ('R0',2): (539, mono(), 1),
            ('Rp',0): (7969, mono(R=3), 2), ('Rp',1): (19440, mono(R=2,U=2), 2), ('Rp',2): (2780, mono(), 2),
            ('Rb',0): (7814, mono(R=3), 2), ('Rb',1): (18900, mono(R=2,U=2), 2), ('Rb',2): (2695, mono(), 2),
            ('Rc',0): (5025, mono(P=2), 2), ('Rc',1): (8400, mono(U=2), 16), ('Rc',2): (1260, mono(P=2), 16),
            ('c1',0): (718, mono(), 2), ('c1',1): (900, mono(U=2), 16), ('c1',2): (180, mono(), 16)}
for nmn, (pp, bb, cc) in mains.items():
    nc, ec = cc
    if nmn != 'Rc':
        rows = [((U_*ec, nc*(one+U_)), (T_, one)), (((one+U_)*ec, nc), ((one+U_)*T_*ec, nc*U_*(one+T_))), ((ec, nc), (T_, one))]
    else:
        rows = [((U_, one+U_), (T_, P_*(one+T_))), ((one+U_, one), ((one+U_)*T_, (P_+U_+P_*U_)*(one+T_))), ((one, one), (T_, P_*(one+T_)))]
    for k, (qq, tt) in enumerate(rows):
        n_, mk, mc = expected[(nmn, k)]
        cert(f'(9) chart {nmn} x row {k}', cleared(Gd, [qq, pp, tt, bb, cc, (V_, one), (one, one)], pc), mk, mc, n_)
Sp = one + P_ + Q_
for k, (qq, tt) in enumerate([((U_, one+U_), (T_, one)), ((one+U_, one), ((one+U_)*T_, U_*(one+T_))), ((one, one), (T_, one))]):
    n_, mk, mc = expected[('c1', k)]
    cert(f'(9) chart c=1 x row {k}', cleared(Gd, [qq, (P_, Sp), tt, (Q_, Sp), (one, one), (V_, one), (one, one)], pc), mk, mc, n_)
sp_, sb, sc_, sq, st = sympy.symbols('p b c q t', positive=True)
th = 1 - sp_*(sb+sc_); P1, Q1, R1 = sp_-1, sp_*sb/th, sp_*sc_/th
check('(5) R_p inverse map', all(sympy.simplify(e) == 0 for e in [Q1/((1+P1)*(1+Q1+R1)) - sb, R1/((1+P1)*(1+Q1+R1)) - sc_]))
U1 = sc_*sq - 1; T1 = U1*st/(sq - U1*st)
check('(6) middle row inverse map', all(sympy.simplify(e) == 0 for e in [(1+U1)/sc_ - sq, (1+U1)*T1/(sc_*U1*(1+T1)) - st]))
P2, U2 = sc_ - 1, sq - 1; W = P2+U2+P2*U2; T2 = W*st/(sq - W*st); T3 = P2*st/(1-P2*st)
check('(7) inverse maps (q>1 chain row; q<=1 old-H rows)', sympy.expand(W - (sc_*sq-1)) == 0 and sympy.simplify((1+U2)*T2/(W*(1+T2)) - st) == 0 and sympy.simplify(T3/(P2*(1+T3)) - st) == 0)
if '--numeric' in sys.argv:
    import numpy as np
    rng = np.random.default_rng(20261127)
    for i, j in [(2,12),(2,15),(12,15)]:
        def gen(i=i, j=j):
            x_ = float(rng.uniform(0.05, 5)); return [(x_, i), (1/x_, j), (1.0, 26)]
        print(f'NUMERIC atom SOC f{i},f{j};f26 at xy=d^2: N=20, restarts 8, worst min Phi = {numeric_atom(rng, gen):.3e}')
    def gen_h():
        m_ = float(rng.uniform(1.01, 5)); return [(1/(m_-1), 3), (m_, 15), (1.0, 26)]
    print(f'NUMERIC atom old-H f3,f15;f26 at t(c-d)=d^2: N=20, restarts 8, worst min Phi = {numeric_atom(rng, gen_h):.3e}')
    def gen_c():
        q_ = float(rng.uniform(0.1, 5)); t_ = float(rng.uniform(0.1, 5)); return [(q_, 0), (t_, 3), (1/q_ + 1/t_, 15), (1.0, 26)]
    print(f'NUMERIC atom chain f0,f3,f15;f26 at c=d^2/q+d^2/t: N=20, restarts 8, worst min Phi = {numeric_atom(rng, gen_c):.3e}')
    found, worst = numeric_stratum_test(T, range(12))
    print(f'NUMERIC theorem: {found} floating points with contacts containing T (seeds 0..11), worst min Phi/max gap = {worst:.3e}')
finish()
