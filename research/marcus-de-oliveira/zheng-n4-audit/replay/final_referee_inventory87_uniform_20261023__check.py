"""Own exact replay for final_referee_inventory87_uniform_20261023 (paper case T_88).

Nothing from the supplement is executed or read. Uses only ../lib/moc4lib.py (built from
the paper's definitions), sympy and python-flint. Prints one line per check; exits nonzero
if an EXACT check fails. Lines starting NUMERIC are floating-point evidence only.
"""
import os
import sys
import time
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'lib'))
import moc4lib as L  # noqa: E402
import flint          # noqa: E402


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


T = L.CONTACT_SETS[88]
check('A1 paper T_88 = {1,5,6,10,13,14}', sorted(T) == [1, 5, 6, 10, 13, 14])
GENS = L.generators(T)
check('B1 generators(T_88) = {0,1,6,10,27,30,34,36,59} (the nine profiles of the document)',
      GENS == [0, 1, 6, 10, 27, 30, 34, 36, 59], str(GENS))
check('B2 non-minor generators are exactly f30, f34', [j for j in GENS if not L.F_IS_MINOR[j]] == [30, 34])
BASIS = [0, 1, 6, 10, 27, 36, 59, 30]
NAMES = ['q', 'a', 'v', 'b', 'w', 't', 'c', 'd']
V = L.F_VALUES
check('B3 kernel of P[T,:] has dimension 8 and the eight profiles of (1) are a basis',
      14 - L.rank([L.P[i] for i in sorted(T)]) == 8 and L.rank([L.F_VECS[j] for j in BASIS]) == 8)
GAPIDX = [20, 16, 11, 7, 12, 3, 0]
ok = all(V[BASIS[k]][GAPIDX[l]] == (1 if k == l else 0) for k in range(8) for l in range(7))
check('B4 (2): q,a,v,b,w,t,c are exactly the gaps g20,g16,g11,g7,g12,g3,g0', ok)
row19 = {NAMES[k]: V[BASIS[k]][19] for k in range(8) if V[BASIS[k]][19]}
check('B5 (2): g19 = q + d, hence d = g19 - g20', row19 == {'q': 1, 'd': 1}, str(row19))
check('B6 circuit (3): f30 + f34 = f0 + f6 + f27 + f36',
      all(V[30][s] + V[34][s] == V[0][s] + V[6][s] + V[27][s] + V[36][s] for s in range(24)))
ok = True
for idx, nm in zip([19, 15, 4, 8], ['q', 'v', 'w', 't']):
    row = {NAMES[k]: V[BASIS[k]][idx] for k in range(8) if V[BASIS[k]][idx]}
    ok &= row == {nm: 1, 'd': 1}
check('B7 gaps g19, g15, g4, g8 are exactly q+d, v+d, w+d, t+d', ok)
FIDX = {tuple(V[j]): j for j in range(82)}


def g_on_f(g, j):
    return FIDX[tuple(L.g_act_values(g, list(V[j])))]


gtr = (L.PERMS[0], L.PERMS[0], 1)
MIN7 = [0, 1, 6, 10, 27, 36, 59]
check('B8 the pure transpose (identity row/column, eps = 1) preserves T, swaps f30 and f34, permutes the seven minors',
      L.g_act_set(gtr, T) == T and g_on_f(gtr, 30) == 34 and g_on_f(gtr, 34) == 30 and
      sorted(g_on_f(gtr, j) for j in MIN7) == MIN7)
# the transpose action is the identity on spectra up to swapping a,b: check it preserves the class of matching vectors
# (tested in the library for all g); nothing further needed.

# SOC atoms (4): orbit of (f1, f9; f22) with pair unordered
REP = (1, 9, 22)
orb = set()
for g in L.G_ELEMENTS:
    i, j, k = (g_on_f(g, x) for x in REP)
    orb.add((frozenset((i, j)), k))
check('B9 SOC orbit of (f1,f9;f22) has 144 elements (pair unordered)', len(orb) == 144, str(len(orb)))
for pair in [(1, 10), (1, 59), (10, 59)]:
    check(f'B10 SOC tuple ({pair[0]},{pair[1]};30) lies in the proved SOC orbit with f30 as the compensated term',
          (frozenset(pair), 30) in orb)

# ------------------------------------------------- H3 allocation (5): exact algebra on samples + proof-read
import random
rng = random.Random(5)
ok = True
for _ in range(2000):
    a, b, c = (Fraction(rng.randint(0, 40), rng.randint(1, 9)) for _ in range(3))
    s = a + b + c
    m = max(a, b, c)
    if 2 * m <= s:
        dab, dac, dbc = (a + b - c) / 2, (a + c - b) / 2, (b + c - a) / 2
        # atom ab uses (dab at f_a, dab at f_b), capacity dab
        res_a, res_b, res_c = dab + dac, dab + dbc, dac + dbc
        ok &= (res_a, res_b, res_c) == (a, b, c) and min(dab, dac, dbc) >= 0 and dab + dac + dbc == s / 2
    else:
        M = m
        lv = sorted([a, b, c])[:2]
        S = sum(lv)
        if S > 0:
            # split M as M*l/S over the two atoms; capacity sqrt(M l^2 / S) each; sum sqrt(M S)
            xs = [M * l / S for l in lv]
            ok &= sum(xs) == M and all(x * l * S == M * l * l for x, l in zip(xs, lv))
check('C1 H3 allocation (5): resources sum exactly to a,b,c; capacities x*y = (alloc)^2 (2000 exact rational samples)', ok)

# ------------------------------------------------- cubic and Hessian identities (flint)
names = ('q', 'a', 'v', 'b', 'w', 't', 'c', 'd', 'H')
R = Ring(flint.fmpq_mpoly_ctx.get(names, 'degrevlex'))
q, a, v, b, w, t, c, d, H = R.gens()
COEF = [q, a, v, b, w, t, c, d]
mvec = [R(0)] * 14
for cf, j in zip(COEF, BASIS):
    for k in range(14):
        x = L.F_VECS[j][k]
        if x:
            mvec[k] += cf * flint.fmpq(x.numerator, x.denominator)
mvec[13] += H
Cval = L.cubic_C(mvec)
e1 = q + v + w + t
e2 = q * v + q * w + q * t + v * w + v * t + w * t
e3 = q * v * w + q * v * t + q * w * t + v * w * t
C0 = a * b * c + e3 + d * e2 + d ** 2 * e1 + d ** 3
Qf = a * b + a * c + b * c + d ** 2
sign = None
for sg in (1, -1):
    if Cval == sg * (C0 + H * Qf):
        sign = sg
check('C2 (7): C(L + H kappa) = +-(C + H Q) with the displayed C and Q (paper cubic C)', sign is not None, f'sign {sign}')

# For the Hessian identity (8) work at d = 1. Every polynomial involved is homogeneous
# (Hess entries are linear forms in (q,..,d,H); N has degree 22; Q^2 G^2 has degree 22), so
# the identity at d = 1 is equivalent to the homogeneous identity; G is re-homogenized below.
one_d = R(1)
mvec = [R(0)] * 14
for cf, j in zip([q, a, v, b, w, t, c, one_d], BASIS):
    for k in range(14):
        x = L.F_VECS[j][k]
        if x:
            mvec[k] += cf * flint.fmpq(x.numerator, x.denominator)
mvec[13] += H
C0 = a * b * c + e3 + e2 + e1 + 1
Qf = a * b + a * c + b * c + 1

# Hessian of the paper cubic as linear forms in m (exact, via sympy once)
import sympy as sp
ms = sp.symbols('m0:14')
Cs = sp.expand(L.cubic_C(ms))
HESS = [[sp.Poly(sp.diff(Cs, ms[i], ms[j]), *ms) for j in range(14)] for i in range(14)]


def lin_to_R(poly):
    out = R(0)
    for mon, cf in zip(poly.monoms(), poly.coeffs()):
        k = mon.index(1)
        out += int(cf) * mvec[k]
    return out


M = [[lin_to_R(HESS[i][j]) if not HESS[i][j].is_zero else R(0) for j in range(14)] for i in range(14)]


def bareiss_det(M):
    """Fraction-free Bareiss elimination with full pivoting on the sparsest nonzero entry
    of the remaining block (row and column swaps tracked in the sign); every division is
    exact and raises otherwise."""
    M = [row[:] for row in M]
    n = len(M)
    sgn = 1
    prev = R(1)
    ndiv = 0
    for k in range(n - 1):
        best = None
        for i in range(k, n):
            for j in range(k, n):
                if M[i][j] != 0 and (best is None or len(M[i][j]) < best[0]):
                    best = (len(M[i][j]), i, j)
        if best is None:
            return R(0), ndiv
        _, pi_, pj_ = best
        if pi_ != k:
            M[k], M[pi_] = M[pi_], M[k]
            sgn = -sgn
        if pj_ != k:
            for row in M:
                row[k], row[pj_] = row[pj_], row[k]
            sgn = -sgn
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                M[i][j] = (M[i][j] * M[k][k] - M[i][k] * M[k][j]) / prev
                ndiv += 1
            M[i][k] = R(0)
        prev = M[k][k]
    return sgn * M[n - 1][n - 1], ndiv


DET, ndiv = bareiss_det(M)
print(f'INFO det Hess C(L + H kappa): {len(DET)} terms, {ndiv} exact divisions, {time.time()-t0:.0f}s', flush=True)
check('C3 Bareiss: 819 exact divisions', ndiv == 819)
degH = DET.degrees()[8]
check('C4 det Hess C(L + H kappa) has degree <= 4 in H (H_kappa has rank 4)', degH <= 4, str(degH))
# N = Q^4 det at H = -C0/Q
coeffsH = [R(0)] * 5
for mon, cf in DET.to_dict().items():
    e = list(mon)
    kH = e[8]
    e[8] = 0
    coeffsH[kH] += R({tuple(e): cf})
N = sum((coeffsH[k] * (-C0) ** k * Qf ** (4 - k) for k in range(5)), R(0))
Q2 = Qf ** 2
G2 = None
try:
    G2 = N / (-4 * Q2)
    ok = G2 * (-4 * Q2) == N
except Exception as ex:   # noqa: BLE001
    ok = False
check('C5 (8): Q^4 det Hess(L - C kappa / Q) is divisible by -4 Q^2', ok)
G = None
try:
    G = G2.sqrt()
    ok = G * G == G2
except Exception:   # noqa: BLE001
    ok = False
check('C6 (8): the quotient is a perfect square G^2 (at d = 1)', ok)
# G (dehomogenized, d = 1) -> homogeneous of degree 9 in (q,a,v,b,w,t,c,d)
GD1 = G.to_dict()
dmon = tuple([0] * 9)
cd9 = GD1.get(dmon, 0)
if cd9 < 0:
    G = -G
    cd9 = -cd9
GD1 = G.to_dict()
check('C6b dehomogenized G has total degree <= 9', all(sum(m_) <= 9 for m_ in GD1))
Ghom = R({tuple(list(m_[:7]) + [9 - sum(m_), 0]): cf for m_, cf in GD1.items()})
check('C7 G has d^9 coefficient +1 after sign normalization', cd9 == 1, str(cd9))
check('C8 G has 2285 terms', len(Ghom) == 2285, str(len(Ghom)))
degs = tuple(Ghom.degrees()[:8])
check('C9 individual degrees of G in (q,a,v,b,w,t,c,d) are (3,2,3,2,3,3,2,9)', degs == (3, 2, 3, 2, 3, 3, 2, 9), str(degs))
negs = sum(1 for cf in G.coeffs() if cf < 0)
print(f'INFO G has {negs} negative coefficients in the original variables', flush=True)
check('C10 G is homogeneous of degree 9 (by construction from the homogeneous identity)', all(sum(m_[:8]) == 9 for m_ in Ghom.to_dict()))

# ------------------------------------------------- charts
names2 = ('q', 'v', 'w', 't', 'X', 'Y', 'Z')
S2 = Ring(flint.fmpq_mpoly_ctx.get(names2, 'degrevlex'))
Q0, Vv, Ww, Tt, Xx, Yy, Zz = S2.gens()
GD = G.to_dict()


def subst_chart(amap, den_list):
    """amap: dict var -> (numerator, list of (den_poly, power)) ; returns cleared polynomial.
    den_list: list of den polys; clearing exponent per den = max over terms."""
    # compute per-term denominator exponents
    terms = []
    maxe = [0] * len(den_list)
    for mon, cf in GD.items():
        qe, ae, ve, be, we, te, ce, de, _ = mon
        num = S2(cf) * Q0 ** qe * Vv ** ve * Ww ** we * Tt ** te
        dex = [0] * len(den_list)
        for var, ex in (('a', ae), ('b', be), ('c', ce)):
            nm, dens = amap[var]
            num *= nm ** ex
            for di, pw in dens:
                dex[di] += pw * ex
        terms.append((num, dex))
        maxe = [max(x, y) for x, y in zip(maxe, dex)]
    out = S2(0)
    for num, dex in terms:
        f = num
        for di, (dp, mx) in enumerate(zip(den_list, maxe)):
            f *= dp ** (mx - dex[di])
        out += f
    return out, maxe


def all_pos(poly):
    return all(cf > 0 for cf in poly.coeffs())


one = S2(1)
# balanced chart (10), d = 1
Dd = 1 + Xx + Yy + Zz
amap = {'a': (Yy + Zz, [(0, 1)]), 'b': (Xx + Zz, [(0, 1)]), 'c': (Xx + Yy, [(0, 1)])}
Gb, mx = subst_chart(amap, [Dd])
const = Gb.to_dict().get((0,) * 7, 0)
check('D1 balanced chart: cleared G has only positive coefficients and constant 1', all_pos(Gb) and const == 1,
      f'{len(Gb)} terms, D^{mx[0]}')
# dominant charts (12)-(13)
ell = 1 + Xx
Mnum = 1 + Yy * ell ** 2
ok_all = True
info = []
for main in ('a', 'b', 'c'):
    leaves = [x for x in ('a', 'b', 'c') if x != main]
    for order in (0, 1):
        l1, l2 = (leaves if order == 0 else leaves[::-1])
        amap = {main: (Mnum, [(0, 1), (1, 1)]), l1: (Zz, [(0, 1), (2, 1)]), l2: (one, [(0, 1), (2, 1)])}
        Gd, mx = subst_chart(amap, [ell, 1 + Yy, 1 + Zz])
        dct = Gd.to_dict()
        pureX = [(m_, cf) for m_, cf in dct.items() if all(m_[i] == 0 for i in (0, 1, 2, 3, 5, 6))]
        ok = all_pos(Gd) and len(pureX) > 0
        ok_all &= ok
        info.append(f'{main}:{l1},{l2}:{len(Gd)}terms,pureX={sorted(pureX)[:2]}')
check('D2 six dominant charts: every cleared G has only positive coefficients and a monomial in X alone',
      ok_all, '; '.join(info))
# solo charts
ok_all = True
info = []
for main in ('a', 'b', 'c'):
    zero = S2(0)
    amap = {x: ((Xx if x == main else zero), []) for x in ('a', 'b', 'c')}
    Gs, mx = subst_chart(amap, [])
    const = Gs.to_dict().get((0,) * 7, 0)
    ok = all_pos(Gs) and const == 1
    ok_all &= ok
    info.append(f'{main}:{len(Gs)}terms')
check('D3 three solo charts: positive coefficients and constant 1', ok_all, ', '.join(info))

# chart inverse checks (exact on rationals)
ok = True
for _ in range(500):
    a_, b_, c_ = (Fraction(rng.randint(0, 30), 31) for _ in range(3))
    # balanced points with a+b+c<2 and triangle inequality
    s_ = a_ + b_ + c_
    if 2 * max(a_, b_, c_) <= s_ and s_ < 2:
        x_, y_, z_ = (b_ + c_ - a_) / 2, (a_ + c_ - b_) / 2, (a_ + b_ - c_) / 2
        den = 1 - x_ - y_ - z_
        X_, Y_, Z_ = x_ / den, y_ / den, z_ / den
        D_ = 1 + X_ + Y_ + Z_
        ok &= ((Y_ + Z_) / D_, (X_ + Z_) / D_, (X_ + Y_) / D_) == (a_, b_, c_) and min(X_, Y_, Z_) >= 0
    elif 2 * max(a_, b_, c_) > s_:
        Mx = max(a_, b_, c_)
        Sx = s_ - Mx
        if Sx > 0 and Mx * Sx < 1:
            X_ = 1 / Sx - 1
            Y_ = (Mx * Sx - Sx ** 2) / (1 - Mx * Sx)
            l_ = 1 + X_
            ok &= X_ > 0 and Y_ >= 0 and (1 + Y_ * l_ ** 2) / (l_ * (1 + Y_)) == Mx
check('D4 chart inverses (10) and (12) recover a,b,c exactly (500 rational samples)', ok)

# numeric falsification of (8)-positivity: G > 0 at random sources with H3 < d
import numpy as np
rs = np.random.default_rng(11)
Gnum = {m_: float(cf) for m_, cf in GD.items()}


def Geval(xv):
    tot = 0.0
    for m_, cf in Gnum.items():
        p_ = cf
        for i_ in range(8):
            if m_[i_]:
                p_ *= xv[i_] ** int(m_[i_])
        tot += p_
    return tot


worst = 1e9
for _ in range(3000):
    xv = rs.exponential(size=8) * rs.choice([0, 1], size=8, p=[0.2, 0.8])
    xv[7] = 1.0
    aa, bb, cc = xv[1], xv[3], xv[6]
    s_ = aa + bb + cc
    m_ = max(aa, bb, cc)
    H3 = s_ / 2 if 2 * m_ <= s_ else np.sqrt(m_ * (s_ - m_))
    if H3 >= 1:
        xv[[1, 3, 6]] *= 0.99 / H3
    worst = min(worst, Geval(xv))
check('NUMERIC E1 G > 0 at 3000 random nonnegative sources with H3 < d = 1 (seed 11)', worst > 0, f'min {worst:.3e}')
print(f'INFO total time {time.time()-t0:.0f}s')

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


H3LOG = []


def h3_extra(pt):
    gp = pt['gaps']
    aa, bb, cc = gp[16], gp[7], gp[0]
    dd = gp[19] - gp[20]
    s_ = aa + bb + cc
    m_ = max(aa, bb, cc)
    h3 = s_ / 2 if 2 * m_ <= s_ else (m_ * (s_ - m_)) ** 0.5
    H3LOG.append((dd, h3))


_np_, _w_ = main_statement_test(T, extra=h3_extra)
check('NUMERIC D8 at stratum points with d > 0, d <= H3(a,b,c) (the document says realizable sources never have H3 < d)',
      all(dd <= h3 + 1e-7 for dd, h3 in H3LOG if dd > 0), str([(round(float(x), 4), round(float(y), 4)) for x, y in H3LOG]))
check('NUMERIC D9 main statement: stratum points found by least squares (seeds 0-5, 40 tries, margin 0: contacts contain T), min Phi over SO(4) by 8-start BFGS >= -1e-9',
      _w_ is None or _w_ >= -1e-9, f'{_np_} points, worst min Phi {_w_}')

bad = [f for f in FAIL if not f.startswith('NUMERIC')]
print('SUMMARY', 'all EXACT checks pass' if not bad else f'FAILED: {bad}')
sys.exit(1 if bad else 0)
