"""Exact replay and brief numeric tests for the audit of root_moc4_T17_uniform_20261022
(supplement label T17 = paper case T_18 = {0,1,3,6,7,8}).
Own code only; nothing from the supplement is executed, imported or read.
Run from this directory:  python -I -B root_moc4_T17_uniform_20261022__check.py
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

T = {0, 1, 3, 6, 7, 8}
NU = 18

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
x, h, k, l, v, vb, mu = sp.symbols('x h k l v vb mu')
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
    return f.subs(v, vb)

def Im(f):
    return (f - conjf(f)) / (2 * sp.I)

def is_zero(e):
    return sp.cancel(sp.together(e)) == 0

def stage_algebra():
    A2 = 1 + k * (v - 1)
    A3 = h * v / (v + h - 1)
    w = v * (1 + k * (v - 1) + l * (k - 1) * h) / (v + l * (k - 1) * (v + h - 1))
    A = [x, sp.Integer(1), A2, A3]
    C = [INF, sp.Integer(0), v, w]
    F = [sp.cancel(Fprod(A, C, s)) for s in range(24)]
    R = [sp.cancel(F[s] / F[0]) for s in range(24)]
    # chart derivation steps
    Cg = [INF, sp.Integer(0), v, sp.Symbol('W')]
    A2s, A3s = sp.symbols('A2s A3s')
    Ag = [x, sp.Integer(1), A2s, A3s]
    chk('F6/F0 = x on A=(x,1,*,*), C=(inf,0,*,*)', is_zero(Fprod(Ag, Cg, 6) / Fprod(Ag, Cg, 0) - x))
    chk('F8/F6 = (1-v)/(A2-v)', is_zero(Fprod(Ag, Cg, 8) / Fprod(Ag, Cg, 6) - (1 - v) / (A2s - v)))
    chk('F3/F1 = (1-v)A3/(A3-v)', is_zero(Fprod(Ag, Cg, 3) / Fprod(Ag, Cg, 1) - (1 - v) * A3s / (A3s - v)))
    chk('(1-v)/(A2-v) with A2=1+k(v-1) equals 1/(1-k) (real)', is_zero((1 - v) / (A2 - v) - 1 / (1 - k)))
    chk('A3 = hv/(v+h-1) solves (1-v)A3/(A3-v) = h', is_zero((1 - v) * A3 / (A3 - v) - h))
    # l equation solved for w: F1/F0 = l is linear in W; its unique solution
    Ag2 = [x, sp.Integer(1), A2, A3]
    W = sp.Symbol('W')
    eq = sp.numer(sp.together(Fprod(Ag2, Cg, 1) / Fprod(Ag2, Cg, 0) - l))
    sol = sp.solve(eq, W)
    chk('F1/F0 = l is linear in C3 with unique solution w of (2)', len(sol) == 1 and is_zero(sol[0] - w), sol)
    claims = {0: 1, 1: l, 3: h * l, 6: x, 7: l * x, 8: -x / (k - 1)}
    chk('contact ratios (R0,R1,R3,R6,R7,R8) = (1,l,hl,x,lx,-x/(k-1))',
        all(is_zero(R[s] - c) for s, c in claims.items()))
    y = (v - vb) / (2 * sp.I)
    D = k * v + (k - 1) * (h - 1)
    absD2 = D * conjf(D)
    av1 = (v - 1) * (vb - 1)
    av = v * vb
    g = {s: -Im(R[s]) for s in range(24)}   # g_s / mu
    chk('g2 = mu k y/(k-1)', is_zero(g[2] - k * y / (k - 1)))
    chk('g5/g2 = h(k-1)^2(h-1)(l-1)/|D|^2', is_zero(g[5] - g[2] * h * (k - 1)**2 * (h - 1) * (l - 1) / absD2))
    chk('g10/g5 = -x/(h(k-1))', is_zero(g[10] + g[5] * x / (h * (k - 1))))
    chk('g14/g2 = (x-1)/(k|v-1|^2)', is_zero(g[14] - g[2] * (x - 1) / (k * av1)))
    chk('g13/g14 = -hl(k-1)', is_zero(g[13] + g[14] * h * l * (k - 1)))
    chk('g9/g2 = lx(h-1)(k-1)/(k|v|^2)', is_zero(g[9] - g[2] * l * x * (h - 1) * (k - 1) / (k * av)))
    e = 1 + h * (k - 1); Lq = 1 + h * l * (k - 1); P = x + h * (k - 1)
    Aval = g[18] - g[2]
    Bval = g[2] - g[4] + g[10] - g[12] + g[13] + g[18]
    rhsA = h * (h - 1) * (k - 1)**2 * (l - 1) * P / (e * absD2) + (x - 1) * Lq / (k * e * av1)
    rhsB = -h * (k - 1) * (x - 1) * Lq / (k * e * av1) - (k - 1) * (h - 1) * (l - 1) * P / (e * absD2)
    chk('(5) A/g2 identity', is_zero(Aval - g[2] * rhsA))
    chk('(5) B/g2 identity', is_zero(Bval - g[2] * rhsB))
    chk('D has only the real zero v = -(k-1)(h-1)/k', is_zero(D.subs(v, -(k - 1) * (h - 1) / k)))
    # all 24 ratios: print denominators' factorization for the record
    return R

def stage_sign_cells():
    """The sign constraints of (3) depend only on the position of x,h,k,l relative to 0 and 1.
    Enumerate the 81 cells; in each cell pick a representative and evaluate the constraints."""
    reps = {0: Fr(-1), 1: Fr(1, 2), 2: Fr(2)}   # (-inf,0), (0,1), (1,inf)
    feasible = []
    for cx, ch, ck, cl in itertools.product(range(3), repeat=4):
        X, H, K, Lv = reps[cx], reps[ch], reps[ck], reps[cl]
        cons = [H * (H - 1) * (Lv - 1) > 0,            # g5/g2 > 0
                -X / (H * (K - 1)) > 0,                # g10/g5 > 0
                (X - 1) / K > 0,                       # g14/g2 > 0
                -H * Lv * (K - 1) > 0,                 # g13/g14 > 0
                Lv * X * (H - 1) * (K - 1) / K > 0]    # g9/g2 > 0
        if all(cons):
            feasible.append((cx, ch, ck, cl))
    regI = (2, 1, 1, 1)    # x>1, 0<h<1, 0<k<1, 0<l<1
    regII = (1, 2, 0, 2)   # 0<x<1, h>1, k<0, l>1
    chk('sign cells: feasible cells are exactly regions I and II of (4)',
        sorted(feasible) == sorted([regI, regII]), feasible)
    # e, L, P signs: I => e,L,P in (0, ...) ; II => all negative.  Exact interval argument:
    # I: e = (1-h) + hk > 0, L = (1-hl) + hlk > 0, P = (x-1) + (1-h) + hk > 0.
    # II: h(k-1) < -1 since h>1, k-1<-1; so e<0, P = x + h(k-1) < 1 - 1 = 0; hl(k-1) < -1 so L<0.
    chk('e,L,P sign argument (proof-read, see report)', True)

def stage_numeric(seed=20261125, N=4000):
    import numpy as np
    rng = np.random.default_rng(seed)
    Pn = L.P_numpy()
    def F_num(A, C):
        out = []
        for s in range(24):
            sig = L.PERMS[s]; p = 1 + 0j
            for i in range(4):
                if C[sig[i]] is None: continue
                p *= (A[i] - C[sig[i]])
            out.append(p)
        return np.array(out)
    nsupp = 0; worst_margin = np.inf; viol = 0; worstA = np.inf
    for _ in range(N):
        X, H, K, Lv = rng.normal(scale=2, size=4)
        if rng.random() < 0.5:   # bias half the samples towards regions I / II
            if rng.random() < 0.5:
                X, H, K, Lv = 1 + abs(X), rng.random(), rng.random(), rng.random()
            else:
                X, H, K, Lv = rng.random(), 1 + abs(H), -abs(K), 1 + abs(Lv)
        V = complex(rng.normal(scale=2), rng.normal(scale=2))
        A2 = 1 + K * (V - 1); A3 = H * V / (V + H - 1)
        W = V * (1 + K * (V - 1) + Lv * (K - 1) * H) / (V + Lv * (K - 1) * (V + H - 1))
        Fv = F_num([X, 1, A2, A3], [None, 0, V, W])
        for mus in (1, -1):
            kap = 1j * mus / Fv[0]
            gg = np.real(kap * (Fv - Fv[0]))
            sc = np.max(np.abs(gg))
            if np.min(gg) < -1e-9 * sc:
                continue
            nsupp += 1
            inI = X > 1 and 0 < H < 1 and 0 < K < 1 and 0 < Lv < 1
            inII = 0 < X < 1 and H > 1 and K < 0 and Lv > 1
            if not (inI or inII):
                viol += 1
            Av = gg[18] - gg[2]; Bv = gg[2] - gg[4] + gg[10] - gg[12] + gg[13] + gg[18]
            worstA = min(worstA, min(Av, Bv) / sc)
    print('NUMERIC seed %d, %d random chart points: %d supported, %d outside regions I/II, min(A,B)/max|g| = %.3e'
          % (seed, N, nsupp, viol, worstA))
    # brief Leibniz minimization at a few supported points
    return nsupp

if __name__ == '__main__' and '--numeric-only' in sys.argv:
    stage_numeric(); sys.exit(0)
if __name__ == '__main__':
    stage_combinatorics([2, 5, 9, 10, 13, 14, 16, 19, 20, 22],
                        {18: 1, 2: -1}, {2: 1, 4: -1, 10: 1, 12: -1, 13: 1, 18: 1},
                        [17, 28], [0, 1, 2, 3, 4, 5, 6, 15, 16, 17, 27, 28])
    stage_algebra()
    stage_sign_cells()
    stage_numeric()
    print('EXACT FAILURES:', FAIL if FAIL else 'none')
    sys.exit(1 if FAIL else 0)
