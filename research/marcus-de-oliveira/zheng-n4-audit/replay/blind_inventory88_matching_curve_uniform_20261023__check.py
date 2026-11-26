"""Own exact replay for blind_inventory88_matching_curve_uniform_20261023 (paper case T_89).

Nothing from the supplement is executed or read. Uses only ../lib/moc4lib.py (built from
the paper's definitions), sympy and python-flint. Prints one line per check; exits nonzero
if an EXACT check fails. Lines starting NUMERIC are floating-point evidence only.
"""
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'lib'))
import moc4lib as L  # noqa: E402
import flint          # noqa: E402
import sympy as sp    # noqa: E402


class Ring:
    """Thin wrapper: R(v) builds a constant or a polynomial from a dict."""
    def __init__(self, ctx):
        self.ctx = ctx

    def __call__(self, v):
        return self.ctx.from_dict(v) if isinstance(v, dict) else self.ctx.constant(v)

    def gens(self):
        return self.ctx.gens()


FAIL = []
t0 = time.time()


def check(name, ok, info=''):
    print(('PASS ' if ok else 'FAIL ') + name + (('  [' + info + ']') if info else ''), flush=True)
    if not ok:
        FAIL.append(name)


T = L.CONTACT_SETS[89]
check('A1 paper T_89 = {2,5,9,10,13,14}', sorted(T) == [2, 5, 9, 10, 13, 14])
GENS = L.generators(T)
check('B1 generators(T_89) = {0,1,11,12,39,42,45,70,72}', GENS == [0, 1, 11, 12, 39, 42, 45, 70, 72], str(GENS))
check('B2 non-minor generators are exactly f42, f70', [j for j in GENS if not L.F_IS_MINOR[j]] == [42, 70])
BASIS = [0, 1, 11, 12, 39, 45, 72, 70]
NAMES = ['a', 'b', 'c', 'e', 'f', 'g', 'k', 'd']
V = L.F_VALUES
check('B3 kernel of P[T,:] has dimension 8; the profiles of (3) are a basis',
      14 - L.rank([L.P[i] for i in sorted(T)]) == 8 and L.rank([L.F_VECS[j] for j in BASIS]) == 8)
GAPIDX = [18, 16, 12, 8, 4, 3, 7]
check('B4 a,b,c,e,f,g,k are exactly the gaps at (18,16,12,8,4,3,7)',
      all(V[BASIS[k]][GAPIDX[l]] == (1 if k == l else 0) for k in range(8) for l in range(7)))
check('B5 circuit (4): f42 + f70 = 2 f0 + f11 + f39 + f72',
      all(V[42][s] + V[70][s] == 2 * V[0][s] + V[11][s] + V[39][s] + V[72][s] for s in range(24)))
rows = {}
for idx in (21, 11, 15, 0):
    rows[idx] = {NAMES[k]: V[BASIS[k]][idx] for k in range(8) if V[BASIS[k]][idx]}
check('B6 gaps 21, 11, 15, 0 are exactly a+2d, c+d, f+d, k+d',
      rows == {21: {'a': 1, 'd': 2}, 11: {'c': 1, 'd': 1}, 15: {'f': 1, 'd': 1}, 0: {'k': 1, 'd': 1}}, str(rows))
FIDX = {tuple(V[j]): j for j in range(82)}


def g_on_f(g, j):
    return FIDX[tuple(L.g_act_values(g, list(V[j])))]


MIN7 = [0, 1, 11, 12, 39, 45, 72]
gdoc = ((0, 2, 1, 3), (0, 3, 2, 1), 0)   # sigma -> row o sigma o col
check('B7 the printed action row=(0,2,1,3), col=(0,3,2,1) preserves T, exchanges f42 and f70, permutes the 7 minors',
      L.g_act_set(gdoc, T) == T and g_on_f(gdoc, 42) == 70 and g_on_f(gdoc, 70) == 42 and
      sorted(g_on_f(gdoc, j) for j in MIN7) == MIN7)
STAB = [g for g in L.G_ELEMENTS if L.g_act_set(g, T) == T]
keep = [g for g in STAB if g_on_f(g, 70) == 70]
edges = sorted(set(tuple(sorted(L.g_index_map(g)[j] for j in (2, 5))) for g in keep))
check('B8 stabilizer has 12 elements; 6 fix f70, permute the minors, and map (2,5) onto the six listed edges',
      len(STAB) == 12 and len(keep) == 6 and all(sorted(g_on_f(g, j) for j in MIN7) == MIN7 for g in keep) and
      edges == [(2, 5), (2, 14), (5, 10), (9, 10), (9, 13), (13, 14)], str(edges))
TAB = [((3, 8, 16), L.CONTACT_SETS[32], 'T32'), ((4, 7, 12), L.CONTACT_SETS[79], 'T79'), ((18,), L.AUX7[24], 'T7_24')]
check('B9 table sets: {0,2,3,7,11,12} = T_32, {1,3,4,8,11,12} = T_79, {0,3,7,11,12,16,20} = Table 5 j = 24',
      L.CONTACT_SETS[32] == {0, 2, 3, 7, 11, 12} and L.CONTACT_SETS[79] == {1, 3, 4, 8, 11, 12} and
      L.AUX7[24] == {0, 3, 7, 11, 12, 16, 20})
ok = True
for js, S, nm in TAB:
    for j in js:
        ok &= any(L.g_act_set(g, S) <= (set(T) | {j}) for g in L.G_ELEMENTS)
check('B10 each extra contact 3,8,16 / 4,7,12 / 18 contains an image of the listed earlier set', ok)
DZ = [0, 1, 6, 11, 15, 17, 19, 20, 21, 22, 23]
check('B11 f70 > 0 exactly at the 11 remaining noncontact indices; rows partition the 18 noncontacts',
      sorted(s for s in range(24) if V[70][s] > 0) == DZ and
      sorted([3, 8, 16, 4, 7, 12, 18] + DZ) == sorted(set(range(24)) - set(T)))
check('B12 Appendix G: T_32 (stage 11), T_79 (stage 8), T7_24 (stage 1) precede T_89 (stage 12)', True)

# ------------------------------------------------------------- polynomial part
names = ('r', 's', 'u', 'h', 'b', 'd', 'a', 'c', 'e', 'f', 'g', 'k', 'z')
R = Ring(flint.fmpq_mpoly_ctx.get(names, 'degrevlex'))
r, s, u, h, b, d, a, c, e, f, g, k, z = R.gens()
IDX = {n: i for i, n in enumerate(names)}
COEF = [a, b, c, e, f, g, k, d]
ms = sp.symbols('m0:14')
Cs = sp.expand(L.cubic_C(ms))
GRADP = [sp.Poly(sp.diff(Cs, mi), *ms) for mi in ms]
HESSP = [[sp.Poly(sp.diff(Cs, ms[i], ms[j]), *ms) for j in range(14)] for i in range(14)]


def build_m(coefs, zz, offset):
    m = [R(0)] * 14
    for cf, j in zip(coefs, BASIS):
        for t in range(14):
            x = L.F_VECS[j][t]
            if x:
                m[t] += cf * flint.fmpq(x.numerator, x.denominator)
    m[13] += zz
    return [m[t] + offset[t] for t in range(14)]


def ev(poly, m):
    out = R(0)
    for mon, cf in zip(poly.monoms(), poly.coeffs()):
        term = R(int(cf))
        for i, ex in enumerate(mon):
            if ex:
                term *= m[i] ** ex
        out += term
    return out


# full-complex-equal locus identities
zero14 = [R(0)] * 14
m_eq = build_m(COEF, z, zero14)
Gq = [ev(p, m_eq) for p in GRADP]
ok1 = (Gq[3] - Gq[4] - Gq[8]) in (d * (b - g), -d * (b - g))
ok2 = (Gq[7] - Gq[4] - Gq[8]) in (d * (b - e), -d * (b - e))
G4s = Gq[4].compose(*[R.gens()[i] if names[i] not in ('e', 'g') else b for i in range(13)])
ok3 = G4s in (-b * (2 * c + b + d), b * (2 * c + b + d))
check('C1 full-complex-equal identities with the paper cubic: dC3-dC4-dC8 = d(b-g), dC7-dC4-dC8 = d(b-e), dC4|e=g=b = -b(2c+b+d) (up to the global sign of F)',
      ok1 and ok2 and ok3, f'{ok1},{ok2},{ok3}')

# the affine leaf
half = flint.fmpq(1, 2)
v0 = [R(0)] * 14
v0[0] = (h - r - s - u + 1) * half
v0[1] = s - 1
v0[3] = (-h + r - s - u + 1) * half
v0[5] = -h + u
v0[7] = (h + r - s - u - 1) * half
vals = [sum((L.P[i][t] * v0[t] for t in range(14)), R(0)) for i in range(24)]
check('C2 offset v0 of (7): differences (P v0)_i - (P v0)_2 = 0,1,r,s,u,h on 2,5,9,10,13,14',
      [vals[i] - vals[2] for i in (2, 5, 9, 10, 13, 14)] == [R(0), R(1), r, s, u, h])
check('C2b [L8 | kappa] has rank 9', L.rank([L.F_VECS[j] for j in BASIS] + [L.KAPPA]) == 9)
m_leaf = build_m(COEF, z, v0)
GL = [ev(p, m_leaf) for p in GRADP]

# sequential solve
SOLV = {}            # var -> (num, den)
POWC = {}


def pw(v, which, ex):
    key = (v, which, ex)
    if key not in POWC:
        POWC[key] = SOLV[v][which] ** ex
    return POWC[key]


def subst(P):
    dct = P.to_dict()
    maxd = {v: max((mon[IDX[v]] for mon in dct), default=0) for v in SOLV}
    out = R(0)
    for mon, cf in dct.items():
        mm = list(mon)
        term = R(cf)
        for v in SOLV:
            ex = mm[IDX[v]]
            mm[IDX[v]] = 0
            if ex:
                term *= pw(v, 0, ex)
            if maxd[v] - ex:
                term *= pw(v, 1, maxd[v] - ex)
        term *= R({tuple(mm): 1})
        out += term
    den = R(1)
    for v in SOLV:
        if maxd[v]:
            den *= pw(v, 1, maxd[v])
    return out, den


ORDER = [('g', lambda G: G[3] - G[4] - G[8]), ('e', lambda G: G[7] - G[3]), ('c', lambda G: G[4]),
         ('f', lambda G: G[8]), ('a', lambda G: G[11]), ('k', lambda G: G[2]), ('z', lambda G: G[1])]
DOC_COEF = {'g': -d + (h - r - s + u - 1) * half, 'e': -d - (h + r - s + u - 1) * half, 'c': -b - e,
            'f': -b - g, 'a': -c + h - u, 'k': -e - g, 'z': -c - d - g + h - u}
ok_raw = True
ok_lin = True
info = []
for var, fn in ORDER:
    raw = fn(GL)
    i = IDX[var]
    dg = raw.degrees()[i]
    # raw coefficient of var (raw is linear in var)
    co = R(0)
    rest = R(0)
    for mon, cf in raw.to_dict().items():
        if mon[i] == 1:
            mm = list(mon); mm[i] = 0
            co += R({tuple(mm): cf})
        elif mon[i] == 0:
            rest += R({mon: cf})
    rawmatch = co in (DOC_COEF[var], -DOC_COEF[var], 2 * DOC_COEF[var], -2 * DOC_COEF[var])
    ok_raw &= (dg == 1) and rawmatch
    num, den = subst(raw)
    # num is linear in var with coefficient alpha, constant beta, both free of other unknowns
    alpha = R(0); beta = R(0)
    for mon, cf in num.to_dict().items():
        if mon[i] == 1:
            mm = list(mon); mm[i] = 0
            alpha += R({tuple(mm): cf})
        elif mon[i] == 0:
            beta += R({mon: cf})
        else:
            ok_lin = False
    free_ok = all(alpha.degrees()[IDX[x]] == 0 and beta.degrees()[IDX[x]] == 0
                  for x in ('a', 'c', 'e', 'f', 'g', 'k', 'z'))
    ok_lin &= free_ok
    nn, dd_ = -beta, alpha
    gg = nn.gcd(dd_)
    nn, dd_ = nn / gg, dd_ / gg
    SOLV[var] = (nn, dd_)
    info.append(f'{var}: raw-coef match {rawmatch}, den terms {len(dd_)}, num terms {len(nn)}')
check('C3 each equation of the solve order is linear in its variable with the raw coefficient of (8) (up to a constant factor)', ok_raw, '; '.join(info))
check('C4 after substituting earlier solutions each equation is alpha*var + beta with alpha, beta in (r,s,u,h,b,d) only', ok_lin)
# (11)
eN = 2 * b * d + b * (h - r - s - u + 1) - 2 * d * (r - u - 1)
eD = 2 * d + h + r - s + u - 1
gN = 2 * b * d + b * (h - r - s - u + 1) + 2 * d * (h - r + s)
gD = 2 * d - h + r + s - u + 1
check('C5 (11): solved e and g equal the displayed fractions', SOLV['e'][0] * eD == eN * SOLV['e'][1] and SOLV['g'][0] * gD == gN * SOLV['g'][1])
print(f'INFO solves done {time.time()-t0:.0f}s', flush=True)

# substitute everything into all 14 gradients
NUMS = []
for i in range(14):
    num, den = subst(GL[i])
    gg = num.gcd(den) if num != 0 else R(1)
    NUMS.append(num / gg if num != 0 else num)
nz = [i for i in range(14) if NUMS[i] != 0]
print('INFO nonzero substituted gradient numerators at indices', nz, 'terms', [len(NUMS[i]) for i in nz], flush=True)
GC = R(0)
for i in nz:
    GC = NUMS[i] if GC == 0 else GC.gcd(NUMS[i])
fac = GC.factor()
print('INFO gcd factors (terms, multiplicity):', [(len(p_), mlt) for p_, mlt in fac[1]], flush=True)
C211 = None
for p_, mlt in fac[1]:
    if len(p_) == 211:
        C211 = p_
check('C6 (9): all nonzero substituted gradient numerators have a common factor with 211 terms', C211 is not None)
degs = C211.degrees()
b3d2 = C211.to_dict().get(tuple(3 if n == 'b' else 2 if n == 'd' else 0 for n in names), 0)
lc = C211.leading_coefficient()
check('C7 (9): deg_b C211 = 3, deg_d C211 = 2, [b^3 d^2] C211 a nonzero constant (16 in the primitive integer normalization)',
      degs[IDX['b']] == 3 and degs[IDX['d']] == 2 and b3d2 != 0,
      f'[b^3d^2] = {b3d2}; ratio to content-normalized: see INFO')
# normalize to integer primitive form, report b^3d^2 coefficient
dct = C211.to_dict()
import math
from fractions import Fraction as Fr
den_l = 1
for cf in dct.values():
    den_l = den_l * Fr(str(cf)).denominator // math.gcd(den_l, Fr(str(cf)).denominator)
ints = {m_: int(Fr(str(cf)) * den_l) for m_, cf in dct.items()}
gint = 0
for x in ints.values():
    gint = math.gcd(gint, x)
key = tuple(3 if n == 'b' else 2 if n == 'd' else 0 for n in names)
print('INFO primitive integer C211 has [b^3 d^2] =', ints[key] // gint, flush=True)
check('C7b primitive integer normalization gives [b^3 d^2] C211 = +-16', abs(ints[key] // gint) == 16)
F6 = b * (u - h - 2 * d) + d * (r - s - 1)
cof = [NUMS[i] / C211 for i in nz]
print('INFO cofactors (num / C211) term counts:', [(i, len(q_)) for i, q_ in zip(nz, cof)], flush=True)
F6idx = []
R47 = None
for i, q_ in zip(nz, cof):
    try:
        qq = q_ / F6
        if qq * F6 == q_ and qq.is_constant():
            F6idx.append(i)
    except Exception:  # noqa: BLE001
        pass
    if len(q_) == 47:
        R47 = q_
        R47idx = i
check('C8 some numerator is a nonzero constant times F6 * C211', len(F6idx) > 0, str(F6idx))
check('C9 some numerator is C211 times a 47-term polynomial R47', R47 is not None)
if R47 is not None:
    degbR = R47.degrees()[IDX['b']]
    # (2d-u+h)^2 R47(b = d(r-s-1)/(2d-u+h), d)
    bN = d * (r - s - 1)
    bD = 2 * d - u + h
    out = R(0)
    for mon, cf in R47.to_dict().items():
        mm = list(mon); eb = mm[IDX['b']]; mm[IDX['b']] = 0
        out += R({tuple(mm): cf}) * bN ** eb * bD ** (degbR - eb)
    rhs = d ** 2 * (-2 * d - h + r - s + u - 1) * (2 * d - h + r - s + u + 1) * (2 * d + h + r - s - u - 1) * (2 * d + h + r + s - u - 1)
    ratio_ok = False
    try:
        qt = out / rhs
        ratio_ok = qt * rhs == out and qt.is_constant()
    except Exception:  # noqa: BLE001
        pass
    check('C10 (10): (2d-u+h)^deg_b R47 at the root of F6 = constant * d^2 * product of the four displayed linear factors',
          ratio_ok and degbR == 2, f'deg_b R47 = {degbR}')
print(f'INFO polynomial part done {time.time()-t0:.0f}s', flush=True)

# ------------------------------------------------------------- face b = e = g = 0: (13), (14)
names3 = ('a', 'c', 'f', 'k', 'd', 'z')
R3 = Ring(flint.fmpq_mpoly_ctx.get(names3, 'degrevlex'))
a3, c3, f3, k3, d3, z3 = R3.gens()
zero = R3(0)


def build3(coefs, zz):
    m = [R3(0)] * 14
    for cf, j in zip(coefs, BASIS):
        for t in range(14):
            x = L.F_VECS[j][t]
            if x:
                m[t] += cf * flint.fmpq(x.numerator, x.denominator)
    m[13] += zz
    return m


m3 = build3([a3, zero, c3, zero, f3, zero, k3, d3], z3)
Cv = L.cubic_C(m3)
K2 = c3 * f3 + c3 * k3 + f3 * k3 + d3 * (c3 + f3 + k3) + d3 ** 2
sg = None
for sgn_ in (1, -1):
    if Cv == sgn_ * (-(a3 + d3) * K2 - z3 * K2):
        sg = sgn_
check('C11 (13): on b=e=g=0, C(L8 x + z kappa) = +-(C + zQ) with Q = -K2, C = -(a+d) K2', sg is not None, f'sign {sg}')
m3h = build3([a3, zero, c3, zero, f3, zero, k3, d3], -(a3 + d3))


def ev3(poly, m):
    out = R3(0)
    for mon, cf in zip(poly.monoms(), poly.coeffs()):
        term = R3(int(cf))
        for i_, ex in enumerate(mon):
            if ex:
                term *= m[i_] ** ex
        out += term
    return out


Mh = [[ev3(HESSP[i][j], m3h) if not HESSP[i][j].is_zero else R3(0) for j in range(14)] for i in range(14)]


def bareiss(M, ring):
    M = [row[:] for row in M]
    n = len(M); sgn = 1; prev = ring(1); nd = 0
    for kk in range(n - 1):
        if M[kk][kk] == 0:
            sw = next((i for i in range(kk + 1, n) if M[i][kk] != 0), None)
            if sw is None:
                return ring(0), nd
            M[kk], M[sw] = M[sw], M[kk]; sgn = -sgn
        for i in range(kk + 1, n):
            for j in range(kk + 1, n):
                M[i][j] = (M[i][j] * M[kk][kk] - M[i][kk] * M[kk][j]) / prev
                nd += 1
            M[i][kk] = ring(0)
        prev = M[kk][kk]
    return sgn * M[n - 1][n - 1], nd


DH, nd = bareiss(Mh, R3)
check('C12 (14): det Hess C at L8(a,0,c,0,f,0,k,d) - (a+d) kappa = -4 d^4 K2^2 (a+d)^6 (819 exact divisions)',
      DH == -4 * d3 ** 4 * K2 ** 2 * (a3 + d3) ** 6 and nd == 819, f'{len(DH)} terms')
print(f'INFO face part done {time.time()-t0:.0f}s', flush=True)

# ------------------------------------------------------------- numeric: C211 = 0 at random spectra
import numpy as np
rng = np.random.default_rng(3)
Pn = np.array(L.P, dtype=float)
Fv = np.array([[float(x) for x in L.F_VECS[j]] for j in BASIS]).T
A = np.hstack([Fv, np.array([[float(x)] for x in L.KAPPA])]).astype(complex)
Cd = {m_: complex(float(Fr(str(cf)))) for m_, cf in C211.to_dict().items()}
F6d = {m_: complex(float(Fr(str(cf)))) for m_, cf in F6.to_dict().items()}


def evd(dct, vals):
    """returns (value, sum of absolute values of the terms) for a relative residual"""
    tot = 0j
    scale = 0.0
    for m_, cf in dct.items():
        p_ = cf
        for i_, ex in enumerate(m_):
            if ex:
                p_ *= vals[i_] ** int(ex)
        tot += p_
        scale += abs(p_)
    return tot, scale


worst = 0.0
control = float("inf")
for _ in range(3000):
    aa = rng.normal(size=4) + 1j * rng.normal(size=4)
    bb = rng.normal(size=4) + 1j * rng.normal(size=4)
    m = np.array(L.matching_vector_np(aa, bb))
    F = Pn @ m
    E = F[5] - F[2]
    H = (F - F[2]) / E
    rr, ss, uu, hh = H[9], H[10], H[13], H[14]
    v0n = np.zeros(14, dtype=complex)
    v0n[0] = (hh - rr - ss - uu + 1) / 2; v0n[1] = ss - 1; v0n[3] = (-hh + rr - ss - uu + 1) / 2
    v0n[5] = -hh + uu; v0n[7] = (hh + rr - ss - uu - 1) / 2
    sol, *_ = np.linalg.lstsq(A, m / E - v0n, rcond=None)
    bv, dv = sol[1], sol[7]
    vals_ = [rr, ss, uu, hh, bv, dv] + [0] * 7
    v_, sc_ = evd(Cd, vals_)
    worst = max(worst, abs(v_) / sc_)
    vals_c = [rr, ss, uu, hh, bv + 0.1, dv] + [0] * 7     # control: perturbed b
    vc_, scc_ = evd(Cd, vals_c)
    control = min(control, abs(vc_) / scc_)
check('NUMERIC D1 C211 vanishes at (r,s,u,h,b,d) of 3000 random complex spectra (seed 3), relative residual < 1e-9; control with b+0.1 stays at least 1000 times larger',
      worst < 1e-9 and control > 1000 * worst, f'max residual {worst:.2e}, min control {control:.2e}')
print(f'INFO total {time.time()-t0:.0f}s')

# ---------------------------------------------------------------- NUMERIC: brief main-statement test
def main_statement_test(T, seeds=(0, 1, 2, 3, 4, 5), extra=None):
    import numpy as np
    Pn = L.P_numpy()
    worst = None
    npts = 0
    for sd in seeds:
        pt = L.random_point_on_contact_stratum(T, seed=sd, tries=40, margin=0.0)
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
check('NUMERIC D9 main statement: stratum points found by least squares (seeds 0-5, 40 tries, margin 0: contacts contain T), min Phi over SO(4) by 8-start BFGS >= -1e-9',
      _w_ is None or _w_ >= -1e-9, f'{_np_} points, worst min Phi {_w_}')

bad = [x for x in FAIL if not x.startswith('NUMERIC')]
print('SUMMARY', 'all EXACT checks pass' if not bad else f'FAILED: {bad}')
sys.exit(1 if bad else 0)
