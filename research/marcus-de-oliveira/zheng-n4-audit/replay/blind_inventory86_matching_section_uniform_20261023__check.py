"""Own exact replay for blind_inventory86_matching_section_uniform_20261023 (paper case T_87).

Nothing from the supplement is executed or read. Uses only the shared audit library
(../lib/moc4lib.py), built from the paper's definitions. Prints one line per check and
exits nonzero if an EXACT check fails. NUMERIC lines are floating-point evidence only.
"""
import os
import sys
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'lib'))
import moc4lib as L  # noqa: E402
import sympy as sp   # noqa: E402

FAIL = []


def check(name, ok, info=''):
    print(('PASS ' if ok else 'FAIL ') + name + (('  [' + info + ']') if info else ''))
    if not ok:
        FAIL.append(name)


T = L.CONTACT_SETS[87]
check('A1 paper T_87 = {4,5,7,9,12,14} (0-based lexicographic), as in the document', sorted(T) == [4, 5, 7, 9, 12, 14])

# ---------------------------------------------------------------- section 2: source
GENS = L.generators(T)
check('B1 generators(T_87) = {f0,f1,f2,f5,f12,f35,f45,f59,f70} as in (2)-(3)',
      GENS == [0, 1, 2, 5, 12, 35, 45, 59, 70], str(GENS))
check('B2 non-minor generators are exactly f35, f70', [j for j in GENS if not L.F_IS_MINOR[j]] == [35, 70])
BASIS = [0, 1, 2, 5, 12, 45, 59, 70]   # coefficients (a,b,c,e,f,g,k,d)
NAMES = ['a', 'b', 'c', 'e', 'f', 'g', 'k', 'd']
rowsT = [L.P[i] for i in sorted(T)]
check('B3 kernel of P[T,:] has dimension 8', 14 - L.rank(rowsT) == 8)
check('B4 the eight profiles of (2) are independent (so a basis of the kernel)',
      L.rank([L.F_VECS[j] for j in BASIS]) == 8)
V = L.F_VALUES
GAPIDX = [18, 16, 13, 10, 8, 3, 2]
ok = all(V[BASIS[k]][GAPIDX[l]] == (1 if k == l else 0) for k in range(8) for l in range(7))
check('B5 coefficients a,b,c,e,f,g,k equal the gaps at (18,16,13,10,8,3,2) exactly', ok)
lhs = [V[35][s] + V[70][s] for s in range(24)]
rhs = [V[1][s] + V[2][s] + V[5][s] + V[12][s] + V[45][s] + V[59][s] for s in range(24)]
check('B6 circuit (3): f35 + f70 = f1 + f2 + f5 + f12 + f45 + f59 on all 24 permutations', lhs == rhs)
NEG = [17, 15, 11, 6, 1, 0]
ok = True
for idx, name in zip(NEG, ['b', 'c', 'e', 'f', 'g', 'k']):
    row = {NAMES[k]: V[BASIS[k]][idx] for k in range(8) if V[BASIS[k]][idx] != 0}
    ok &= (row == {name: 1, 'd': 1})
check('B7 gaps at (17,15,11,6,1,0) are exactly b+d, c+d, e+d, f+d, g+d, k+d', ok)

# group action on generators
FIDX = {tuple(V[j]): j for j in range(82)}


def g_on_f(g, j):
    return FIDX[tuple(L.g_act_values(g, list(V[j])))]


STAB = [g for g in L.G_ELEMENTS if L.g_act_set(g, T) == T]
check('B8 stabilizer of T_87 in G has 12 elements', len(STAB) == 12, str(len(STAB)))
MINORS = [0, 1, 2, 5, 12, 45, 59]
swap = [g for g in STAB if g_on_f(g, 35) == 70 and g_on_f(g, 70) == 35
        and sorted(g_on_f(g, j) for j in MINORS) == MINORS]
check('B9 some stabilizer element exchanges f35 and f70 and permutes the seven minor generators',
      len(swap) > 0, f'{len(swap)} such elements, e.g. {swap[0] if swap else None}')
keep = [g for g in STAB if g_on_f(g, 70) == 70]
check('B10 six stabilizer elements fix f70, and each permutes the seven minor generators',
      len(keep) == 6 and all(sorted(g_on_f(g, j) for j in MINORS) == MINORS for g in keep))
edges = sorted(set(tuple(sorted(L.g_index_map(g)[j] for j in (4, 9))) for g in keep))
check('B11 images of the pair (4,9) under those six elements are the six listed edges',
      edges == [(4, 9), (4, 14), (5, 7), (5, 12), (7, 14), (9, 12)], str(edges))
# connected: cycle 4-9-12-5-7-14-4
adj = {v: set() for v in T}
for a_, b_ in edges:
    adj[a_].add(b_); adj[b_].add(a_)
seen, stack = {4}, [4]
while stack:
    v = stack.pop()
    for w in adj[v]:
        if w not in seen:
            seen.add(w); stack.append(w)
check('B12 the six edges connect all six contacts', seen == set(T))

# extra-contact table
TABLE = [((2, 10, 13), 68), ((3, 8, 16), 51), ((18,), 61)]
ok = True
info = []
for js, nu in TABLE:
    for j in js:
        Tj = set(T) | {j}
        g = next((g for g in L.G_ELEMENTS if L.g_act_set(g, L.CONTACT_SETS[nu]) <= Tj), None)
        ok &= g is not None
        info.append(f'{j}->T{nu}')
check('B13 extra contacts 2,10,13 / 3,8,16 / 18 contain an image of T_68 / T_51 / T_61', ok, ', '.join(info))
DZERO = [0, 1, 6, 11, 15, 17, 19, 20, 21, 22, 23]
check('B14 f70 > 0 exactly at the 11 remaining indices (d = 0 if such a gap vanishes)',
      sorted(s for s in range(24) if V[70][s] > 0) == DZERO and all(V[70][s] >= 0 for s in range(24)))
check('B15 the table rows partition the 18 noncontact indices',
      sorted([2, 10, 13, 3, 8, 16, 18] + DZERO) == sorted(set(range(24)) - set(T)))
STAGE = {68: 6, 51: 7, 61: 9, 87: 10}
check('B16 Appendix G: T_68, T_51, T_61 (stages 6, 7, 9) precede T_87 (stage 10)',
      all(STAGE[k] < STAGE[87] for k in (68, 51, 61)))

# ---------------------------------------------------------------- cubic identities
X = sp.symbols('a b c e f g k d')
z = sp.Symbol('z')
r, s, u, h = sp.symbols('r s u h')
ms = sp.symbols('m0:14')
Cm = sp.expand(L.cubic_C(ms))
GRAD = [sp.diff(Cm, mi) for mi in ms]
HESS_KAPPA = L.hess_kappa()


def lin_comb_vec(coeffs, offset=None):
    vec = [sp.Integer(0)] * 14
    for cf, j in zip(coeffs, BASIS):
        for t in range(14):
            vec[t] += cf * sp.Rational(L.F_VECS[j][t].numerator, L.F_VECS[j][t].denominator)
    vec[13] += z
    if offset is not None:
        vec = [vec[t] + offset[t] for t in range(14)]
    return [sp.expand(x) for x in vec]


def grads_at(vec):
    sub = dict(zip(ms, vec))
    return [sp.expand(gr.subs(sub, simultaneous=True)) for gr in GRAD]


def span_solve(polys, target, gens):
    """Find constant lambda (over Q(params)) with sum lambda_i polys_i == target."""
    lam = sp.symbols('l0:14')
    expr = sp.expand(sum(l_ * p_ for l_, p_ in zip(lam, polys)) - target)
    eqs = sp.Poly(expr, *gens).coeffs()
    sol = sp.solve(eqs, lam, dict=True)
    return sol


# (5): full-complex-equal locus m = L8 x + z kappa
G0 = grads_at(lin_comb_vec(X))
sol = span_solve(G0, -X[7] ** 2, list(X) + [z])
check('C5 eq. (5): -d^2 is a constant linear combination of the 14 partials on m = L8 x + z kappa',
      len(sol) > 0, str(sol[0] if sol else None))
# same with the doc's specific combination (d13 - d2), assuming the document's F is the paper's C
val = sp.expand(G0[13] - G0[2])
check('C5b with the paper cubic C, dC/dm13 - dC/dm2 = -d^2 exactly on that locus', val == -X[7] ** 2, str(val))

# leaf: (P v0)_i - (P v0)_4 = H_i, H4=0, H9=1, H5=r, H7=s, H12=u, H14=h
Tl = [4, 9, 5, 7, 12, 14]
targets = [0, 1, r, s, u, h]
vv = sp.symbols('v0:14')
eqs = [sum(L.P[i][t] * vv[t] for t in range(14)) - tg for i, tg in zip(Tl, targets)]
# sparse particular solution: set all but 6 pivot coordinates to zero
Mrows = sp.Matrix([[L.P[i][t] for t in range(14)] for i in Tl])
piv = Mrows.rref()[1]
solv = sp.solve([e.subs({vv[t]: 0 for t in range(14) if t not in piv}) for e in eqs], [vv[t] for t in piv], dict=True)[0]
V0 = [sp.expand(solv.get(vv[t], 0)) if t in piv else sp.Integer(0) for t in range(14)]
vals0 = [sp.expand(sum(L.P[i][t] * V0[t] for t in range(14))) for i in range(24)]
check('C6 own particular solution v0 gives the five fixed differences', all(sp.expand(vals0[i] - vals0[4] - tg) == 0 for i, tg in zip(Tl, targets)), str(V0))
check('C6b [L8 | kappa] has rank 9 and the five difference rows have rank 5',
      L.rank([L.F_VECS[j] for j in BASIS] + [L.KAPPA]) == 9 and
      L.rank([[L.P[i][t] - L.P[4][t] for t in range(14)] for i in (9, 5, 7, 12, 14)]) == 5)
GL = grads_at(lin_comb_vec(X, V0))
# z-coefficient structure: grad C(m + z kappa) = grad C(m) + z H_kappa m
mvec = lin_comb_vec(X, V0)
mvec0 = [sp.expand(x.subs(z, 0)) for x in mvec]
ok = all(sp.expand(GL[i] - grads_at(mvec0)[i] - z * sum(HESS_KAPPA[i][t] * mvec0[t] for t in range(14))) == 0 for i in range(14))
check('C7a grad C(m + z kappa) = grad C(m) + z H_kappa m on the leaf (only m0,m1,m5,m6 partials see z)', ok)
zc = [sp.expand(sum(HESS_KAPPA[i][t] * mvec0[t] for t in range(14))) for i in range(14)]
print('INFO z-coefficients of the partials on the leaf (own v0):', {i: zc[i] for i in range(14) if zc[i] != 0})

# (8): a constant combination of partials that depends on d only
dd = X[7]
lam = sp.symbols('l0:14')
expr = sp.expand(sum(l_ * p_ for l_, p_ in zip(lam, GL)))
P_ = sp.Poly(expr, *(list(X) + [z]))
eqs = []
for mon, cf in zip(P_.monoms(), P_.coeffs()):
    if any(mon[i] for i in range(7)) or mon[8]:   # any monomial containing a non-d variable or z
        eqs.append(cf)
solL = sp.solve(eqs, lam, dict=True)
found = None
if solL:
    q = sp.expand(expr.subs(solL[0]))
    # choose a normalization: make the d^2 coefficient -1 if possible
    free = sorted(q.free_symbols & set(lam), key=str)
    qd = sp.Poly(q, dd)
    found = (q, free)
check('C8 some constant combination of partials restricted to the leaf involves d only', solL != [] and found is not None and set(found[0].free_symbols) - set(lam) <= {dd, r, s, u, h})
q, free = found
print('INFO the d-only combination (free multipliers', free, '):', sp.factor(q))
# Extract a basis of the d-only space by setting free multipliers to unit vectors
basisq = []
for fs in free:
    qq = sp.expand(q.subs({f_: (1 if f_ == fs else 0) for f_ in free}))
    if qq != 0:
        basisq.append(qq)
print('INFO d-only combinations:', [sp.factor(b_) for b_ in basisq])
# compare with the document's monic quadratic up to a real shift d -> d + y(r,s,u,h)
docq = -dd ** 2 + (h + r + s - u - 1) * dd + r * (u - h)
yy = sp.Symbol('y')
match = None
for qq in basisq:
    pq = sp.Poly(qq, dd)
    if pq.degree() != 2:
        continue
    lead = pq.LC()
    shifted = sp.expand(qq.subs(dd, dd + yy) * (-1) / lead)
    sol_y = sp.solve(sp.Poly(sp.expand(shifted - docq), dd).coeffs(), yy, dict=True)
    if sol_y:
        match = (lead, sol_y)
check('C8b that combination is a quadratic in d with constant nonzero d^2 coefficient, equal to the document (8) after a real shift of d',
      match is not None, str(match))
DISC = sp.expand((h + r + s - u - 1) ** 2 + 4 * r * (u - h))
print('INFO discriminant of (8):', sp.factor(DISC))
# direct check with the document's indices assuming the paper cubic
val8 = sp.expand(GL[13] - GL[2])
print('INFO with own v0, dC/dm13 - dC/dm2 on the leaf =', sp.factor(val8))

# (7): some partial has z-coefficient linear with nonzero imaginary part at a strict point
print('INFO doc (7) z-coefficient (-d-e-g+h-u); own z-coefficients above differ by the choice of v0 and coordinates')

# ---------------------------------------------------------------- numeric: (8) on actual spectra
import random
import numpy as np


def numeric_leaf_test(trials=4000, seed=7):
    """Random complex spectra: compute m, decompose relative to T_87 with E = F9 - F4,
    check the quadratic (8)'s analogue (own v0) vanishes."""
    rng = np.random.default_rng(seed)
    Pn = np.array(L.P, dtype=float)
    Fv = np.array([[float(x) for x in L.F_VECS[j]] for j in BASIS]).T   # 14 x 8
    A = np.hstack([Fv, np.array([[float(x)] for x in L.KAPPA])])
    worst = 0.0
    fq = sp.lambdify((dd, r, s, u, h), basisq[0], 'numpy')
    v0f = sp.lambdify((r, s, u, h), V0, 'numpy')
    for _ in range(trials):
        a = rng.normal(size=4) + 1j * rng.normal(size=4)
        b = rng.normal(size=4) + 1j * rng.normal(size=4)
        m = np.array(L.matching_vector_np(a, b))
        F = Pn @ m
        E = F[9] - F[4]
        Hh = (F - F[4]) / E
        mm = m / E
        rr, ss, uu, hh = Hh[5], Hh[7], Hh[12], Hh[14]
        # m/E - v0(r,s,u,h) = L8 x + z kappa with complex (r,s,u,h): identity is polynomial so holds for complex parameters
        v0n = np.array(v0f(rr, ss, uu, hh), dtype=complex)
        sol, *_ = np.linalg.lstsq(A.astype(complex), mm - v0n, rcond=None)
        dval = sol[7]
        val = abs(fq(dval, rr, ss, uu, hh)) / (1 + abs(dval) ** 2)
        worst = max(worst, val)
    return worst


w = numeric_leaf_test()
check('NUMERIC D1 d-only quadratic vanishes at 4000 random complex spectra (seed 7), max normalized residual < 1e-6', w < 1e-6, f'{w:.2e}')


# ---------------------------------------------------------------- NUMERIC: brief main-statement test
def main_statement_test(T, seeds=(0, 1, 2, 3, 4, 5), extra=None):
    import numpy as np
    Pn = L.P_numpy()
    worst = None
    npts = 0
    for sd in seeds:
        pt = L.random_point_on_contact_stratum(T, seed=sd, tries=15, margin=1e-3)
        if pt is None:
            continue
        npts += 1
        cvec, *_ = np.linalg.lstsq(Pn, pt['gaps'], rcond=None)
        val, _ = L.min_phi_over_SO4(cvec, restarts=8, seed=sd)
        if extra is not None:
            extra(pt)
        worst = val if worst is None else min(worst, val)
    return npts, worst


_np_, _w_ = main_statement_test(T)
check('NUMERIC D9 main statement: stratum points found by least squares (seeds 0-5, margin 1e-3), min Phi over SO(4) by 8-start BFGS >= -1e-9',
      _w_ is None or _w_ >= -1e-9, f'{_np_} points, worst min Phi {_w_}')

print('SUMMARY', 'all EXACT checks pass' if not [f for f in FAIL if not f.startswith('NUMERIC')] else f'FAILED: {FAIL}')
sys.exit(1 if [f for f in FAIL if not f.startswith('NUMERIC')] else 0)
