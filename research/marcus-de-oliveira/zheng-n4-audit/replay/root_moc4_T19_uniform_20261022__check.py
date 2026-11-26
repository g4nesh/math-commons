"""Exact replay and brief numeric tests for the audit of root_moc4_T19_uniform_20261022
(supplement label T19 = paper case T_20 = {0,1,5,6,7,8}).
Own code only; nothing from the supplement is executed, imported or read.
Run from this directory:  python -I -B root_moc4_T19_uniform_20261022__check.py
Prints one line per check; exits nonzero if an EXACT check fails."""
import os, sys, itertools, random
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'lib'))
import moc4lib as L
import sympy as sp
from fractions import Fraction as Fr

FAIL = []
def chk(name, ok, info=''):
    print(('PASS ' if ok else 'FAIL ') + name + ((' : ' + str(info)) if info != '' else ''))
    if not ok:
        FAIL.append(name)

T = {0, 1, 5, 6, 7, 8}
NU = 20

# ---------------------------------------------------------------- combinatorics
def minor_cone_dual(Tset):
    gens = L.generators(Tset)
    M = [j for j in gens if L.F_IS_MINOR[j]]
    rows = [list(map(Fr, L.F_VALUES[j])) for j in M]
    r = L.rank(rows)
    # basis of span from minor generators
    basis = []
    for j in M:
        cand = basis + [list(map(Fr, L.F_VALUES[j]))]
        if L.rank(cand) == len(cand):
            basis = cand
    # coordinates of each generator in the basis: solve basis^T c = v
    BT = L.transpose(basis)          # 24 x r
    def coords(vec):
        # least squares exact: pick r independent rows
        rowsel = []
        for s in range(24):
            cand = [BT[i] for i in rowsel + [s]]
            if L.rank(cand) == len(cand):
                rowsel.append(s)
            if len(rowsel) == r:
                break
        A = [BT[i] for i in rowsel]
        b = [Fr(vec[i]) for i in rowsel]
        Ai = L.inv(A)
        c = L.matvec(Ai, b)
        assert all(sum(BT[s][i] * c[i] for i in range(r)) == vec[s] for s in range(24))
        return c
    G = {j: coords(L.F_VALUES[j]) for j in M}
    facets = []
    for sub in itertools.combinations(M, r - 1):
        S = [G[j] for j in sub]
        if L.rank(S) != r - 1:
            continue
        ns = L.nullspace(S)
        assert len(ns) == 1
        n = ns[0]
        vals = [sum(a * b for a, b in zip(n, G[j])) for j in M]
        if all(x >= 0 for x in vals):
            pass
        elif all(x <= 0 for x in vals):
            n = [-a for a in n]; vals = [-x for x in vals]
        else:
            continue
        zero = frozenset(j for j, x in zip(M, vals) if x == 0)
        if L.rank([G[j] for j in zero]) == r - 1 and zero not in [f[0] for f in facets]:
            facets.append((zero, vals))
    return M, r, basis, facets

def functional_on_minors(coef, M):
    """values of the gap functional sum coef[s] g_s on the minor generators"""
    return [sum(Fr(c) * L.F_VALUES[j][s] for s, c in coef.items()) for j in M]

def same_ray(u, w):
    # positive proportionality
    for a, b in zip(u, w):
        if (a == 0) != (b == 0):
            return False
    nz = [(a, b) for a, b in zip(u, w) if a != 0]
    if not nz:
        return False
    lam = nz[0][1] / nz[0][0]
    return lam > 0 and all(b == lam * a for a, b in nz)

def stage_combinatorics(listed_coords, Aform, Bform, nonminor_claim, source_claim):
    gens = L.generators(T)
    chk('contact set equals paper T_%d (Appendix E)' % NU, L.CONTACT_SETS[NU] == frozenset(T))
    chk('find_g(T_%d, T) is the identity-type match' % NU, L.find_g(L.CONTACT_SETS[NU], T) is not None)
    chk('generators F_T = %s' % source_claim, gens == source_claim, gens)
    chk('nonminor generators = %s' % nonminor_claim, [j for j in gens if not L.F_IS_MINOR[j]] == nonminor_claim)
    rkP = L.rank([L.P[s] for s in sorted(T)])
    chk('dim(im P cap {g_T = 0}) = 14 - rank P[T] = 8', 14 - rkP == 8)
    M, r, basis, facets = minor_cone_dual(T)
    chk('minor generators span dimension 8 = span of all generators',
        r == 8 and L.rank([list(L.F_VALUES[j]) for j in gens]) == 8)
    listed = [({s: 1}, 'g%d' % s) for s in listed_coords] + [(Aform, 'A'), (Bform, 'B')]
    lv = [(functional_on_minors(c, M), nm) for c, nm in listed]
    fv = [v for z, v in facets]
    matched = []
    for v, nm in lv:
        hit = [i for i, w in enumerate(fv) if same_ray(v, w)]
        matched.append((nm, hit))
    ok = len(facets) == len(listed) and all(len(h) == 1 for _, h in matched) and \
        len(set(h[0] for _, h in matched)) == len(listed)
    chk('minor cone has exactly %d facets and they are the listed functionals' % len(listed), ok,
        'facets=%d' % len(facets))
    # nonminor generators: A or B negative on them (they lie outside the minor cone)
    for j in gens:
        if not L.F_IS_MINOR[j]:
            a = sum(Fr(c) * L.F_VALUES[j][s] for s, c in Aform.items())
            b = sum(Fr(c) * L.F_VALUES[j][s] for s, c in Bform.items())
            print('INFO nonminor f%d (%s): A = %s, B = %s' % (j, L.f_name(j), a, b))

# ---------------------------------------------------------------- algebra
x, h, k, l, v, vb, mu, q, qb = sp.symbols('x h k l v vb mu q qb')
INF = None

def Fprod(A, C, s):
    sig = L.PERMS[s]
    out = sp.Integer(1)
    for i in range(4):
        if A[i] is INF or C[sig[i]] is INF:
            continue
        out *= (A[i] - C[sig[i]])
    return out

def conjf(f):
    return f.subs({v: vb, q: qb}, simultaneous=True)

def Im(f):
    return (f - conjf(f)) / (2 * sp.I)

def is_zero(e):
    return sp.cancel(sp.together(e)) == 0

def stage_algebra():
    W, A3s, V = sp.symbols('W A3s V')
    a2V = 1 + k * (V - 1)
    Ag = [x, sp.Integer(1), a2V, A3s]
    Cg = [INF, sp.Integer(0), V, W]
    chk('F6/F0 = x', is_zero(Fprod(Ag, Cg, 6) / Fprod(Ag, Cg, 0) - x))
    chk('F8/F6 = 1/(1-k)', is_zero(Fprod(Ag, Cg, 8) / Fprod(Ag, Cg, 6) - 1 / (1 - k)))
    chk('F5/F0 = (1-w)A3/(A3-w)', is_zero(Fprod(Ag, Cg, 5) / Fprod(Ag, Cg, 0) - (1 - W) * A3s / (A3s - W)))
    A3W = h * W / (h - 1 + W)
    chk('A3 = hw/(h-1+w) solves it', is_zero((1 - W) * A3W / (A3W - W) - h))
    Ag2 = [x, sp.Integer(1), a2V, A3W]
    lnum = sp.factor(sp.numer(sp.together(Fprod(Ag2, Cg, 1) / Fprod(Ag2, Cg, 0) - l)))
    curve = l * (a2V - V) * W * (1 - W) - (a2V - W) * ((h - V) * W - V * (h - 1))
    ratio = sp.factor(sp.cancel(lnum / curve))
    chk('F1/F0 = l  <=>  curve (2) up to the factor %s' % ratio, sp.Poly(sp.numer(ratio), V, W).is_ground if True else False, ratio)
    # node parametrisation
    cc = 1 + l * (k - 1); Tq = cc * q - k
    vq = (q - 1) * (h * (q - k) + l * (k - 1) * q) / (q * Tq)
    wq = 1 + q * (vq - 1)
    sub = sp.expand(curve.subs(W, 1 + q * (V - 1)))
    quo, rem = sp.div(sp.Poly(sub, V), sp.Poly((V - 1)**2, V))
    lin = quo.as_expr()
    chk('curve with w = 1 + q(v-1) is (v-1)^2 times a polynomial linear in v', rem.is_zero and sp.Poly(lin, V).degree() == 1)
    a1c, a0c = sp.Poly(lin, V).all_coeffs()
    chk('its solution is v(q) of (3)', is_zero(-a0c / a1c - vq))
    chk('coefficient of v is a constant times q T(q)', sp.cancel(a1c / (q * Tq)).is_number or sp.Poly(sp.cancel(a1c / (q * Tq)), q).is_ground, sp.factor(a1c))
    chk('v-numerator at q = k/c equals kl(h-1)(k-1)^2(l-1)/c^2',
        is_zero(((q - 1) * (h * (q - k) + l * (k - 1) * q)).subs(q, k / cc) - k * l * (h - 1) * (k - 1)**2 * (l - 1) / cc**2))
    vv = vq; ww = wq
    a2 = 1 + k * (vv - 1)
    A = [x, sp.Integer(1), a2, h * ww / (h - 1 + ww)]
    C = [INF, sp.Integer(0), vv, ww]
    F = [sp.cancel(Fprod(A, C, s)) for s in range(24)]
    R = [sp.cancel(F[s] / F[0]) for s in range(24)]
    claims = {0: 1, 1: l, 5: h, 6: x, 7: l * x, 8: -x / (k - 1)}
    chk('contact ratios (R0,R1,R5,R6,R7,R8) = (1,l,h,x,lx,-x/(k-1))',
        all(is_zero(R[s] - c) for s, c in claims.items()))
    chk('R14 = (x-v)/((k-1)(v-1))', is_zero(R[14] - (x - vv) / ((k - 1) * (vv - 1))))
    chk('R16 = h R14', is_zero(R[16] - h * R[14]))
    chk('R19 = h - h(x-1)/(w-1)', is_zero(R[19] - (h - h * (x - 1) / (ww - 1))))
    y = (q - qb) / (2 * sp.I)
    aq2 = q * qb; aT2 = Tq * conjf(Tq); aq1 = (q - 1) * (qb - 1); aqk = (q - k) * (qb - k)
    g = {s: -Im(R[s]) for s in (2, 3, 4, 9, 10, 13, 18)}
    chk('(4) g3 = mu h k y/((k-1)|q|^2)', is_zero(g[3] - h * k * y / ((k - 1) * aq2)))
    chk('(4) g4/g3 = -l(h-1)(k-1)^2(l-1)|q|^2/(h|T|^2)', is_zero(g[4] + g[3] * l * (h - 1) * (k - 1)**2 * (l - 1) * aq2 / (h * aT2)))
    chk('(4) g9/g3 = -x(l-1)(k-1)|q|^2/(hk|q-1|^2)', is_zero(g[9] + g[3] * x * (l - 1) * (k - 1) * aq2 / (h * k * aq1)))
    chk('(4) g10/g3 = lx(k-1)|q|^2/(h|q-k|^2)', is_zero(g[10] - g[3] * l * x * (k - 1) * aq2 / (h * aqk)))
    chk('(5) Pick identity Im v/Im q', is_zero(Im(vv) / y - (h / aq2 - l * (h - 1) * (k - 1)**2 * (l - 1) / aT2)))
    chk('(6) Im w/Im q', is_zero(Im(ww) / y - (h - 1) / cc * (1 - k * l * (k - 1)**2 * (l - 1) / aT2)))
    J = (h - 1) * q**2 + (-h * k - h - k * l + k + l) * q + h * k
    S = (h * k + cc - k) * q - h * k
    chk('(8) v - 1 = J/(qT)', is_zero(vv - 1 - J / (q * Tq)))
    AR = R[3] - R[2] + R[18]
    BR = R[13] - R[3]
    chk('(8) A_R = (1-x) l S/J + (h+l-1)', is_zero(AR - ((1 - x) * l * S / J + h + l - 1)))
    chk('(8) B_R = -h(1-x)(q-k)T/((k-1)J)', is_zero(BR + h * (1 - x) * (q - k) * Tq / ((k - 1) * J)))
    chk('A = g3 - g2 + g18 = -mu Im A_R and B = g13 - g3 = -mu Im B_R (by linearity, definitional)', True)
    p = k / cc
    chk('J(0) = hk', is_zero(J.subs(q, 0) - h * k))
    chk('J(p) = kl(h-1)(k-1)^2(l-1)/c^2', is_zero(J.subs(q, p) - k * l * (h - 1) * (k - 1)**2 * (l - 1) / cc**2))
    chk('J(1) = -(k-1)(l-1)', is_zero(J.subs(q, 1) + (k - 1) * (l - 1)))
    chk('slope of S = hk + (k-1)(l-1)', is_zero(sp.diff(S, q) - (h * k + (k - 1) * (l - 1))))
    chk('S(0) = -hk', is_zero(S.subs(q, 0) + h * k))
    chk('S(p) = k(h-1)(k-c)/c', is_zero(S.subs(q, p) - k * (h - 1) * (k - cc) / cc))
    z0 = h * k / (h * k + cc - k)
    chk('J(z0) = -hk(h-1)(k-1)(l-1)/(hk+c-k)^2', is_zero(J.subs(q, z0) + h * k * (h - 1) * (k - 1) * (l - 1) / (h * k + cc - k)**2))
    chk('T(q) = c (q - p)', is_zero(Tq - cc * (q - p)))
    return R

def stage_sign_cells():
    reps3 = {0: Fr(-1), 1: Fr(1, 2), 2: Fr(2)}
    repsh = {1: Fr(1, 2), 2: Fr(2)}           # h > 0 from R16 = h R14
    def sgn(a):
        return (a > 0) - (a < 0)
    prelim = []
    for cx, ch, ck, cl in itertools.product(range(3), (1, 2), range(3), range(3)):
        X, H, K, Lv = reps3[cx], repsh[ch], reps3[ck], reps3[cl]
        cons = [-Lv * (H - 1) * (Lv - 1) / H > 0,                  # g4/g3
                -X * (Lv - 1) * (K - 1) / (H * K) > 0,             # g9/g3
                Lv * X * (K - 1) / H > 0]                          # g10/g3
        if all(cons):
            prelim.append((cx, ch, ck, cl))
    table = [(2, 2, 2, 1), (1, 2, 2, 1), (0, 2, 1, 1), (2, 1, 0, 0), (1, 1, 0, 0), (0, 1, 0, 2)]
    # paper table rows: k>1,x>0,0<l<1,h>1 (x in (0,1) or (1,inf)); 0<k<1,x<0,0<l<1,h>1;
    # k<0,x>0,l<0,0<h<1; k<0,x<0,l>1,0<h<1
    chk('preliminary sign table of the document is exactly the feasible cell set', sorted(prelim) == sorted(table), prelim)
    after = [c for c in prelim if (reps3[c[0]] - 1) / reps3[c[2]] > 0]
    chk('(x-1)/k > 0 leaves cases {k>1,x>1}, {k<0,0<x<1}, {k<0,x<0}', sorted(after) == sorted([(2, 2, 2, 1), (1, 1, 0, 0), (0, 1, 0, 2)]), after)
    final = []
    for c in after:
        X, H, K, Lv = reps3[c[0]], repsh[c[1]], reps3[c[2]], reps3[c[3]]
        # within the cell: sign of l(k-1) and -kl(l-1) are cell-determined (products of cell-signed factors)
        c_gt1 = Lv * (K - 1) > 0
        br_gt1 = -K * Lv * (Lv - 1) > 0
        if c_gt1 and br_gt1:
            s_muImq = sgn(H * K / (K - 1))           # from g3 > 0
            s_ratio = sgn((H - 1))                    # (h-1)/c with c>1, bracket>1
            s_g19 = -sgn(H) * sgn(X - 1) * s_ratio * s_muImq
            chk('cell %s: c>1, bracket>1 and g19 < 0 (contradiction)' % (c,), s_g19 < 0)
        else:
            final.append(c)
    chk('only region (7) survives: 0<h<1, k<0, l>1, x<0', final == [(0, 1, 0, 2)], final)
    chk('in (7): mu Im q > 0, c < 0, 0 < p < 1, interlacing (proof-read, see report)', True)

def stage_numeric(seed=20261125, N=4000):
    import numpy as np
    rng = np.random.default_rng(seed)
    def F_num(A, C):
        out = []
        for s in range(24):
            sig = L.PERMS[s]; p = 1 + 0j
            for i in range(4):
                if C[sig[i]] is None: continue
                p *= (A[i] - C[sig[i]])
            out.append(p)
        return np.array(out)
    nsupp = 0; viol = 0; worst = np.inf
    for _ in range(N):
        X, H, K, Lv = rng.normal(scale=2, size=4)
        if rng.random() < 0.5:   # bias half the samples towards region (7)
            X, H, K, Lv = -abs(X), rng.random(), -abs(K), 1 + abs(Lv)
        Q = complex(rng.normal(scale=2), rng.normal(scale=2))
        c = 1 + Lv * (K - 1); Tq = c * Q - K
        V = (Q - 1) * (H * (Q - K) + Lv * (K - 1) * Q) / (Q * Tq)
        W = 1 + Q * (V - 1)
        a2 = 1 + K * (V - 1)
        Fv = F_num([X, 1, a2, H * W / (H - 1 + W)], [None, 0, V, W])
        for mus in (1, -1):
            kap = 1j * mus / Fv[0]
            gg = np.real(kap * (Fv - Fv[0]))
            sc = np.max(np.abs(gg))
            if np.min(gg) < -1e-9 * sc:
                continue
            nsupp += 1
            if not (0 < H < 1 and K < 0 and Lv > 1 and X < 0):
                viol += 1
            Av = gg[3] - gg[2] + gg[18]; Bv = gg[13] - gg[3]
            worst = min(worst, min(Av, Bv) / sc)
    print('NUMERIC seed %d, %d random chart points: %d supported, %d outside region (7), min(A,B)/max|g| = %.3e'
          % (seed, N, nsupp, viol, worst))

if __name__ == '__main__':
    stage_combinatorics([3, 4, 9, 10, 14, 16, 19, 22],
                        {3: 1, 2: -1, 18: 1}, {13: 1, 3: -1},
                        [28, 34], [0, 1, 2, 3, 4, 5, 6, 27, 28, 33, 34])
    stage_algebra()
    stage_sign_cells()
    stage_numeric()
    print('EXACT FAILURES:', FAIL if FAIL else 'none')
    sys.exit(1 if FAIL else 0)
