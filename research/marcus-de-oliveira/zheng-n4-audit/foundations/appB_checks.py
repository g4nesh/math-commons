#!/usr/bin/env python3
"""Appendix B audit (Theorem B.1, Lemma B.2, Corollary B.3 of Zheng, "MOC for n = 4").

Own code only. Reads nothing from the supplement. Deterministic seeds.
Exit status: 0 all good; 1 an EXACT check failed; 3 a numeric anomaly was seen.

Levels: sections marked EXACT use Fractions / Gaussian rationals / sympy polynomial rings.
Sections marked NUMERIC use floats (numpy / scipy HiGHS) and are evidence only.

Indices 0..3 for rows, columns and eigenvalues, as in Appendix B.
"""
import itertools, random, sys, time
from fractions import Fraction as Fr

import numpy as np
from scipy.optimize import linprog
from sympy.polys.rings import ring
from sympy.polys.domains import QQ

T0 = time.time()
EXACT_FAIL = []
NUM_ANOM = []


def check(name, cond, info=""):
    if not cond:
        EXACT_FAIL.append(name)
        print("EXACT-FAIL", name, info, flush=True)


def anomaly(name, info=""):
    NUM_ANOM.append(name)
    print("NUMERIC-ANOMALY", name, info, flush=True)


PERMS = list(itertools.permutations(range(4)))
PIDX = {p: k for k, p in enumerate(PERMS)}
ID = (0, 1, 2, 3)


def sgn(p):
    s = 1
    for i in range(4):
        for j in range(i + 1, 4):
            if p[i] > p[j]:
                s = -s
    return s


def inv(p):
    q = [0] * 4
    for i, v in enumerate(p):
        q[v] = i
    return tuple(q)


def compose(a, b):  # (a o b)(i) = a(b(i))
    return tuple(a[b[i]] for i in range(4))


# ---------------------------------------------------------------- Gaussian rationals
class Q2:
    __slots__ = ("r", "i")

    def __init__(s, r, i=0):
        s.r = r if isinstance(r, Fr) else Fr(r)
        s.i = i if isinstance(i, Fr) else Fr(i)

    def __add__(s, o):
        o = _c(o); return Q2(s.r + o.r, s.i + o.i)
    __radd__ = __add__

    def __sub__(s, o):
        o = _c(o); return Q2(s.r - o.r, s.i - o.i)

    def __rsub__(s, o):
        return _c(o) - s

    def __mul__(s, o):
        o = _c(o); return Q2(s.r * o.r - s.i * o.i, s.r * o.i + s.i * o.r)
    __rmul__ = __mul__

    def __neg__(s):
        return Q2(-s.r, -s.i)

    def conj(s):
        return Q2(s.r, -s.i)

    def abs2(s):
        return s.r * s.r + s.i * s.i

    def inv(s):
        n = s.abs2(); return Q2(s.r / n, -s.i / n)

    def iszero(s):
        return s.r == 0 and s.i == 0

    def __eq__(s, o):
        o = _c(o); return s.r == o.r and s.i == o.i

    def __hash__(s):
        return hash((s.r, s.i))

    def __repr__(s):
        return f"({s.r}+{s.i}i)"


def _c(x):
    return x if isinstance(x, Q2) else Q2(x)


def mzeros(n, m=None):
    return [[Q2(0) for _ in range(m or n)] for _ in range(n)]


def eye(n):
    M = mzeros(n)
    for i in range(n):
        M[i][i] = Q2(1)
    return M


def mmul(A, B):
    n, k, m = len(A), len(B), len(B[0])
    return [[sum((A[i][t] * B[t][j] for t in range(k)), Q2(0)) for j in range(m)] for i in range(n)]


def madd(A, B, c=1):
    return [[A[i][j] + c * B[i][j] for j in range(len(A[0]))] for i in range(len(A))]


def ct(A):
    return [[A[j][i].conj() for j in range(len(A))] for i in range(len(A[0]))]


def det(M):
    n = len(M); A = [row[:] for row in M]; d = Q2(1)
    for c in range(n):
        p = next((r for r in range(c, n) if not A[r][c].iszero()), None)
        if p is None:
            return Q2(0)
        if p != c:
            A[c], A[p] = A[p], A[c]; d = -d
        d = d * A[c][c]; iv = A[c][c].inv()
        for r in range(c + 1, n):
            f = A[r][c] * iv
            if not f.iszero():
                for k in range(c, n):
                    A[r][k] = A[r][k] - f * A[c][k]
    return d


def minv(M):
    n = len(M); A = [row[:] + e for row, e in zip(M, eye(n))]
    for c in range(n):
        p = next(r for r in range(c, n) if not A[r][c].iszero())
        A[c], A[p] = A[p], A[c]
        iv = A[c][c].inv(); A[c] = [x * iv for x in A[c]]
        for r in range(n):
            if r != c and not A[r][c].iszero():
                f = A[r][c]; A[r] = [x - f * y for x, y in zip(A[r], A[c])]
    return [row[n:] for row in A]


def sub(M, rows, cols):
    return [[M[i][j] for j in cols] for i in rows]


def cayley_unitary(rng, n=4, real=False, amp=3):
    """U = (I-K)(I+K)^{-1}, K = iH (H Hermitian, Gaussian-integer entries) or K real skew."""
    K = mzeros(n)
    for j in range(n):
        if not real:
            K[j][j] = Q2(0, rng.randint(-amp, amp))
        for k in range(j + 1, n):
            x, y = rng.randint(-amp, amp), (0 if real else rng.randint(-amp, amp))
            # K = iH with H_jk = x + i y  ->  K_jk = -y + i x ; K_kj = -conj(K_jk)
            if real:
                K[j][k] = Q2(x); K[k][j] = Q2(-x)
            else:
                K[j][k] = Q2(-y, x); K[k][j] = Q2(y, x)
    I = eye(n)
    return mmul(madd(I, K, -1), minv(madd(I, K)))


def perm_matrix(p, phases=None):
    M = mzeros(4)
    for i in range(4):
        M[i][p[i]] = phases[i] if phases else Q2(1)
    return M


def is_unitary(U):
    P = mmul(U, ct(U))
    return all(P[i][j] == (Q2(1) if i == j else Q2(0)) for i in range(len(U)) for j in range(len(U)))


def blockdiag(A, B):
    n, m = len(A), len(B); M = mzeros(n + m)
    for i in range(n):
        for j in range(n):
            M[i][j] = A[i][j]
    for i in range(m):
        for j in range(m):
            M[n + i][n + j] = B[i][j]
    return M


# ---------------------------------------------------------------- linear algebra over Q
def _fr(x):
    return x if isinstance(x, (Fr, int)) else Fr(int(x.numerator), int(x.denominator))


def rank_Q(rows):
    A = [[Fr(_fr(x)) for x in r] for r in rows]; rk = 0
    ncol = len(A[0]) if A else 0
    for c in range(ncol):
        p = next((r for r in range(rk, len(A)) if A[r][c] != 0), None)
        if p is None:
            continue
        A[rk], A[p] = A[p], A[rk]
        for r in range(len(A)):
            if r != rk and A[r][c] != 0:
                f = A[r][c] / A[rk][c]; A[r] = [x - f * y for x, y in zip(A[r], A[rk])]
        rk += 1
    return rk


def nullvec_Q(cols):
    """Nonzero lambda with sum_k lambda_k cols[k] = 0 (cols are equal-length lists)."""
    m, k = len(cols[0]), len(cols)
    A = [[Fr(cols[j][i]) for j in range(k)] for i in range(m)]
    piv = []; r = 0
    for c in range(k):
        p = next((t for t in range(r, m) if A[t][c] != 0), None)
        if p is None:
            continue
        A[r], A[p] = A[p], A[r]
        f0 = A[r][c]; A[r] = [x / f0 for x in A[r]]
        for t in range(m):
            if t != r and A[t][c] != 0:
                f = A[t][c]; A[t] = [x - f * y for x, y in zip(A[t], A[r])]
        piv.append(c); r += 1
    free = next(c for c in range(k) if c not in piv)
    lam = [Fr(0)] * k; lam[free] = Fr(1)
    for row, c in enumerate(piv):
        lam[c] = -A[row][free]
    return lam


def pvec(p):
    return [1 if p[i] == j else 0 for i in range(4) for j in range(4)]


# ---------------------------------------------------------------- Theorem B.1 construction (exact)
def phi(p):
    r = p.index(0); k, l = [i for i in (1, 2, 3) if i != r]
    s = list(p); s[k], s[l] = p[l], p[k]
    return tuple(s)


def m_of_perm(sig):
    return [[1 if {sig[0], sig[i]} == ({1, 2, 3} - {j}) else 0 for j in (1, 2, 3)] for i in (1, 2, 3)]


def m_of_U(U):
    out = []
    for i in (1, 2, 3):
        row = []
        for j in (1, 2, 3):
            cols = [k for k in (1, 2, 3) if k != j]
            row.append(det(sub(U, [0, i], cols)).abs2())
        out.append(row)
    return out


def birkhoff_exact(W):
    W = [r[:] for r in W]; terms = []
    while any(W[i][j] != 0 for i in range(4) for j in range(4)):
        best = None
        for p in PERMS:
            if all(W[i][p[i]] > 0 for i in range(4)):
                mn = min(W[i][p[i]] for i in range(4))
                if best is None or mn > best[1]:
                    best = (p, mn)
        if best is None:
            return None
        p, mn = best; terms.append([p, mn])
        for i in range(4):
            W[i][p[i]] -= mn
    return terms


def caratheodory_exact(terms, maxk):
    while len(terms) > maxk:
        cols = [[1] + pvec(p) for p, _ in terms]
        lam = nullvec_Q(cols)
        if all(l <= 0 for l in lam):
            lam = [-l for l in lam]
        t = min(w / l for (p, w), l in zip(terms, lam) if l > 0)
        terms = [[p, w - t * l] for (p, w), l in zip(terms, lam)]
        terms = [[p, w] for p, w in terms if w != 0]
    return terms


def b1_weights_exact(U, stats=None):
    s = 1 - U[0][0].abs2()
    mu = {}
    if s == 0:
        mu[ID] = Fr(1); return mu, None
    m = m_of_U(U)
    W = [[Fr(0)] + [U[0][j].abs2() / s for j in (1, 2, 3)]]
    for i in (1, 2, 3):
        W.append([U[i][0].abs2() / s] + [m[i - 1][j] / s for j in range(3)])
    terms = birkhoff_exact(W)
    check("B.1 Birkhoff greedy found a permutation in the support", terms is not None)
    if stats is not None:
        stats.setdefault("greedy_terms", []).append(len(terms))
    terms = caratheodory_exact(terms, 9)
    if stats is not None:
        stats.setdefault("terms_after_caratheodory", []).append(len(terms))
    rec = [[Fr(0)] * 4 for _ in range(4)]
    for p, w in terms:
        check("B.1 Birkhoff weight positive", w > 0)
        check("B.1 W00=0 forces pi(0)!=0", p[0] != 0)
        for i in range(4):
            rec[i][p[i]] += w
        sg = phi(p); mu[sg] = mu.get(sg, Fr(0)) + s * w
    check("B.1 Birkhoff decomposition reproduces W", rec == W)
    mu[ID] = mu.get(ID, Fr(0)) + (1 - s)
    return mu, W


def lhs_B1(U, d, a, e):
    """det(diag(d,a1,a2,a3) - U diag(d,e1,e2,e3) U^*), exact."""
    Da = [d] + list(a); De = [d] + list(e)
    Uh = ct(U)
    M = [[sum((U[i][k] * De[k] * Uh[k][j] for k in range(4)), Q2(0)) for j in range(4)] for i in range(4)]
    A = [[(_c(Da[i]) if i == j else Q2(0)) - M[i][j] for j in range(4)] for i in range(4)]
    return det(A)


def rhs_mu(mu, a, e):
    tot = Q2(0)
    for sg, w in mu.items():
        pr = Q2(1)
        for i in range(4):
            pr = pr * (a[i] - e[sg[i]])
        tot = tot + w * pr
    return tot


def section_B1_exact(nsamp, rng):
    print("\n== Theorem B.1: exact replay on Gaussian-rational unitary matrices ==", flush=True)
    fam = []
    for _ in range(nsamp):
        fam.append(("cayley", cayley_unitary(rng)))
    # degenerate families
    for p in PERMS:
        ph = [Q2(1), Q2(0, 1), Q2(-1), Q2(0, -1)]
        fam.append(("perm", perm_matrix(p, ph)))
    for _ in range(4):
        fam.append(("1+3 block (u00 unimodular, s=0)", blockdiag([[Q2(0, 1)]], cayley_unitary(rng, 3))))
        fam.append(("3+1 block", blockdiag(cayley_unitary(rng, 3), [[Q2(-1)]])))
        fam.append(("2+2 block", blockdiag(cayley_unitary(rng, 2), cayley_unitary(rng, 2))))
        B = blockdiag(cayley_unitary(rng, 2), cayley_unitary(rng, 2))
        fam.append(("u00=0 (permuted 2+2 block)", mmul(perm_matrix((2, 1, 0, 3)), B)))
        fam.append(("u00=0 (permuted 3+1 block)", mmul(perm_matrix((3, 1, 2, 0)), blockdiag(cayley_unitary(rng, 3), [[Q2(1)]]))))
        fam.append(("real orthogonal", cayley_unitary(rng, real=True)))
        fam.append(("perm*2+2 block*perm", mmul(mmul(perm_matrix(rng.choice(PERMS)), blockdiag(cayley_unitary(rng, 2), cayley_unitary(rng, 2))), perm_matrix(rng.choice(PERMS)))))
    grid = list(itertools.product((0, 1), repeat=6))
    stats = {}
    maxsupp = 0
    for kind, U in fam:
        check("unitary " + kind, is_unitary(U))
        s = 1 - U[0][0].abs2()
        m = m_of_U(U)
        # (B.3) marginals
        for i in (1, 2, 3):
            check("(B.3) row marginal", sum(m[i - 1]) == s - U[i][0].abs2(), kind)
            check("(B.3) column marginal", sum(m[k][i - 1] for k in range(3)) == s - U[0][i].abs2(), kind)
        check("(B.3) total", sum(map(sum, m)) == 2 * s, kind)
        # Gram-matrix step for the first identity (rows 0,i restricted to [3])
        for i in (1, 2, 3):
            R = sub(U, [0, i], [1, 2, 3]); G = mmul(R, ct(R))
            G2 = [[Q2(s), -U[0][0] * U[i][0].conj()], [-U[i][0] * U[0][0].conj(), Q2(1 - U[i][0].abs2())]]
            check("p38 Gram matrix display", G == G2 or all(G[x][y] == G2[x][y] for x in range(2) for y in range(2)), kind)
            check("p38 det Gram = s-|ui0|^2", det(G) == Q2(s - U[i][0].abs2()), kind)
        # column identity sum_i |u0k uil - u0l uik|^2 = |u0k|^2+|u0l|^2
        for k, l in ((1, 2), (1, 3), (2, 3)):
            lhs = sum((U[0][k] * U[i][l] - U[0][l] * U[i][k]).abs2() for i in range(4))
            check("p38 column identity", lhs == U[0][k].abs2() + U[0][l].abs2(), kind)
        # complementary minor |det U[{0,k,l},[3]]|^2 = |u_i0|^2
        for i in (1, 2, 3):
            rows = [0] + [k for k in (1, 2, 3) if k != i]
            check("p39 complementary 3x3 minor", det(sub(U, rows, [1, 2, 3])).abs2() == U[i][0].abs2(), kind)
        mu, W = b1_weights_exact(U, stats)
        if W is not None:
            check("(B.4) W doubly stochastic", all(sum(r) == 1 for r in W) and all(sum(W[i][j] for i in range(4)) == 1 for j in range(4)) and all(x >= 0 for r in W for x in r), kind)
        else:
            check("s=0 case: row/col 0 trivial", all(U[0][j].iszero() and U[j][0].iszero() for j in (1, 2, 3)), kind)
        check("mu >= 0", all(w >= 0 for w in mu.values()), kind)
        check("sum mu = 1", sum(mu.values()) == 1, kind)
        check("(B.2) sum_{sigma(0)=0} mu = |u00|^2", sum(w for sg, w in mu.items() if sg[0] == 0) == U[0][0].abs2(), kind)
        nz = sum(1 for w in mu.values() if w != 0); maxsupp = max(maxsupp, nz)
        check("at most ten nonzero weights", nz <= 10, f"{kind} nz={nz}")
        mm = [[sum(w * m_of_perm(sg)[i][j] for sg, w in mu.items()) for j in range(3)] for i in range(3)]
        check("(B.5)", mm == m, kind)
        # (B.1) at d=0: both sides multi-affine in a1..a3,e1..e3 -> {0,1}^6 grid is a proof
        for g in grid:
            a = (0,) + g[:3]; e = (0,) + g[3:]
            if not (lhs_B1(U, Q2(0), [Q2(x) for x in g[:3]], [Q2(x) for x in g[3:]]) == rhs_mu(mu, [Q2(x) for x in a], [Q2(x) for x in e])):
                check("(B.1) at d=0 on {0,1}^6 grid", False, kind); break
        # (B.1) with d != 0 at random Gaussian-rational points (spot replay of the shift step)
        for _ in range(2):
            rq = lambda: Q2(Fr(rng.randint(-9, 9), rng.randint(1, 5)), Fr(rng.randint(-9, 9), rng.randint(1, 5)))
            d = rq(); a = [rq() for _ in range(3)]; e = [rq() for _ in range(3)]
            check("(B.1) with d != 0", lhs_B1(U, d, a, e) == rhs_mu(mu, [d] + a, [d] + e), kind)
    print(f"  samples: {nsamp} Cayley + {len(fam) - nsamp} degenerate; greedy Birkhoff terms max "
          f"{max(stats['greedy_terms'])}, after Caratheodory max {max(stats['terms_after_caratheodory'])}; "
          f"max nonzero weights {maxsupp}", flush=True)


# ---------------------------------------------------------------- generic symbolic checks for B.1
def section_B1_symbolic():
    print("\n== Theorem B.1: symbolic (generic) identities ==", flush=True)
    names = ["a1", "a2", "a3", "e1", "e2", "e3"] + [f"u{i}{j}" for i in range(4) for j in range(4)] + [f"v{i}{j}" for i in range(4) for j in range(4)]
    R, *g = ring(",".join(names), QQ)
    a = [R(0)] + g[0:3]; e = [R(0)] + g[3:6]
    u = [[g[6 + 4 * i + j] for j in range(4)] for i in range(4)]
    v = [[g[22 + 4 * i + j] for j in range(4)] for i in range(4)]   # v = entrywise conjugate, independent symbols
    M = [[sum(u[i][k] * e[k] * v[j][k] for k in range(4)) for j in range(4)] for i in range(4)]
    A = [[(a[i] if i == j else R(0)) - M[i][j] for j in range(4)] for i in range(4)]

    def rdet(X, n=None):
        n = len(X); tot = R(0)
        for p in itertools.permutations(range(n)):
            s = 1
            for i in range(n):
                for j in range(i + 1, n):
                    if p[i] > p[j]:
                        s = -s
            t = R(s)
            for i in range(n):
                t *= X[i][p[i]]
            tot += t
        return tot
    lhs = rdet(A)
    rhs = R(0)
    for k in (1, 2, 3):
        for Irest in itertools.combinations((1, 2, 3), k - 1):
            I = (0,) + Irest
            for J in itertools.combinations((1, 2, 3), k):
                t = R((-1) ** k)
                for i in range(4):
                    if i not in I:
                        t *= a[i]
                for j in J:
                    t *= e[j]
                t *= rdet([[u[i][j] for j in J] for i in I]) * rdet([[v[i][j] for j in J] for i in I])
                rhs += t
    check("generic Cauchy-Binet / principal-minor expansion behind (B.7) (a0=e0=0)", lhs == rhs)
    print("  (B.7) precursor: det(diag(a)-U diag(e) V^T) = sum_{I∋0,J⊂[3]} (-1)^|I| prod a prod e detU[I,J] detV[I,J]  (generic U,V): ",
          lhs == rhs, flush=True)

    # (B.7) + (B.3) => (B.8): linear algebra with symbols
    R2, a1, a2, a3, e1, e2, e3, s, x1, x2, y1, y2, m11, m12, m21, m22 = ring(
        "a1,a2,a3,e1,e2,e3,s,x1,x2,y1,y2,m11,m12,m21,m22", QQ)
    aa = [a1, a2, a3]; ee = [e1, e2, e3]
    x = [x1, x2, s - x1 - x2]; y = [y1, y2, s - y1 - y2]
    m = [[m11, m12, None], [m21, m22, None], [None, None, None]]
    m[0][2] = s - y[0] - m11 - m12
    m[1][2] = s - y[1] - m21 - m22
    m[2][0] = s - x[0] - m11 - m21
    m[2][1] = s - x[1] - m12 - m22
    m[2][2] = s - y[2] - m[2][0] - m[2][1]
    check("(B.3) column-3 marginal implied by the others", m[0][2] + m[1][2] + m[2][2] == s - x[2])
    astar = a1 * a2 * a3; estar = e1 * e2 * e3
    ahat = [aa[1] * aa[2], aa[0] * aa[2], aa[0] * aa[1]]
    ehat = [ee[1] * ee[2], ee[0] * ee[2], ee[0] * ee[1]]
    half = QQ(1, 2)
    H = [[ahat[i] * ehat[j] + astar * (ee[j] - half * (e1 + e2 + e3)) + estar * (aa[i] - half * (a1 + a2 + a3)) for j in range(3)] for i in range(3)]
    B7 = -astar * sum(ee[j] * x[j] for j in range(3)) + sum(ahat[i] * ehat[j] * m[i][j] for i in range(3) for j in range(3)) - estar * sum(aa[i] * y[i] for i in range(3))
    B8 = sum(H[i][j] * m[i][j] for i in range(3) for j in range(3))
    check("(B.8) = (B.7) modulo (B.3)", B7 == B8)
    print("  (B.7) and (B.3) imply (B.8) [symbolic, m_ij, |u0j|^2, |ui0|^2, s as free symbols]:", B7 == B8, flush=True)

    # permutation identity prod (a_i - e_sigma(i)) = sum H_ij m_ij(sigma), a0=e0=0
    R3, b1, b2, b3, f1, f2, f3 = ring("a1,a2,a3,e1,e2,e3", QQ)
    A3 = [R3(0), b1, b2, b3]; E3 = [R3(0), f1, f2, f3]
    astar = b1 * b2 * b3; estar = f1 * f2 * f3
    ah = [b2 * b3, b1 * b3, b1 * b2]; eh = [f2 * f3, f1 * f3, f1 * f2]
    H3 = [[ah[i] * eh[j] + astar * (E3[j + 1] - half * (f1 + f2 + f3)) + estar * (A3[i + 1] - half * (b1 + b2 + b3)) for j in range(3)] for i in range(3)]
    ok = True
    for sg in PERMS:
        pr = R3(1)
        for i in range(4):
            pr *= (A3[i] - E3[sg[i]])
        ms = m_of_perm(sg)
        ok &= (pr == sum(H3[i][j] * ms[i][j] for i in range(3) for j in range(3)))
    check("p39 permutation-matrix identity for all 24 sigma", ok)
    print("  prod(a_i-e_sigma(i)) = sum H_ij m_ij(sigma) for all 24 sigma:", ok, flush=True)

    # the nine H_ij are linearly independent (so (B.1) <=> (B.5), and (B.2) is forced)
    rows = []
    mons = sorted({mon for i in range(3) for j in range(3) for mon in H3[i][j].keys()})
    for i in range(3):
        for j in range(3):
            rows.append([H3[i][j].get(mon, 0) for mon in mons])
    rk = rank_Q(rows)
    check("H_ij linearly independent", rk == 9)
    print("  rank of the nine H_ij as polynomials:", rk, flush=True)

    # phi: involution on the 18 perms not fixing 0; m(phi(pi)) = lower-right block of P_pi
    ok = True
    nf = [p for p in PERMS if p[0] != 0]
    for p in nf:
        q = phi(p)
        ok &= (q[0] != 0 and phi(q) == p)
        ok &= m_of_perm(q) == [[1 if p[i] == j else 0 for j in (1, 2, 3)] for i in (1, 2, 3)]
    check("phi involution and m(phi(pi)) = lower-right block of P_pi", ok)
    print("  phi properties on the 18 permutations with pi(0)!=0:", ok, flush=True)
    # every sigma with sigma(0)!=0 has exactly two ones in m(sigma); sigma(0)=0 has none
    ok = all(sum(map(sum, m_of_perm(p))) == (2 if p[0] != 0 else 0) for p in PERMS)
    check("m(sigma) has two ones iff sigma(0)!=0", ok)
    # face {W00=0} of the Birkhoff polytope has affine dimension 8
    base = pvec(nf[0])
    rk = rank_Q([[x - y for x, y in zip(pvec(p), base)] for p in nf[1:]])
    check("face W00=0 has affine dimension 8", rk == 8)
    print("  affine dimension of face {W00=0}:", rk, flush=True)


# ---------------------------------------------------------------- Lemma B.2
def maxflow_exact(cap, s, t, n):
    flow = [[Fr(0)] * n for _ in range(n)]
    total = Fr(0)
    while True:
        par = [-1] * n; par[s] = s; q = [s]
        while q and par[t] == -1:
            x = q.pop(0)
            for y in range(n):
                if par[y] == -1 and cap[x][y] - flow[x][y] > 0:
                    par[y] = x; q.append(y)
        if par[t] == -1:
            return total, flow
        b = None; y = t
        while y != s:
            x = par[y]; r = cap[x][y] - flow[x][y]; b = r if b is None else min(b, r); y = x
        y = t
        while y != s:
            x = par[y]; flow[x][y] += b; flow[y][x] -= b; y = x
        total += b


EDGES = list(itertools.combinations(range(4), 2))


def b2_weights_exact(p, q, check_it=True):
    """Lemma B.2 law for e=(alpha,beta,gamma,gamma) with alpha-column p, beta-column q."""
    ap = [x.abs2() for x in p]; bp = [x.abs2() for x in q]
    w = {(i, j): (p[i] * q[j] - p[j] * q[i]).abs2() for i, j in EDGES}
    if check_it:
        check("B.2 sum w = 1", sum(w.values()) == 1)
        for i in range(4):
            check("B.2 incident mass", sum(w[tuple(sorted((i, j)))] for j in range(4) if j != i) == ap[i] + bp[i])
        for k in range(1, 5):
            for S in itertools.combinations(range(4), k):
                inside = sum(w[e] for e in EDGES if e[0] in S and e[1] in S)
                pS = sum(ap[i] for i in S); qS = sum(bp[i] for i in S)
                ip = sum((p[i] * q[i].conj() for i in S), Q2(0)).abs2()
                check("B.2 subset identity ||pS||^2||qS||^2-|<pS,qS>|^2", inside == pS * qS - ip)
                check("B.2 subset inequality", inside <= pS)
    # network: 0 source, 1..6 edge nodes, 7..10 vertices, 11 sink
    n = 12; cap = [[Fr(0)] * n for _ in range(n)]
    for k, e in enumerate(EDGES):
        cap[0][1 + k] = w[e]
        for v in e:
            cap[1 + k][7 + v] = Fr(2)
    for v in range(4):
        cap[7 + v][11] = ap[v]
    val, fl = maxflow_exact(cap, 0, 11, n)
    if check_it:
        check("B.2 max flow = 1", val == 1)
    theta = {}
    for k, (i, j) in enumerate(EDGES):
        theta[(i, j)] = fl[1 + k][7 + i]; theta[(j, i)] = fl[1 + k][7 + j]
    if check_it:
        check("B.2 theta >= 0", all(t >= 0 for t in theta.values()))
        for i, j in EDGES:
            check("B.2 theta_ij+theta_ji = w_ij", theta[(i, j)] + theta[(j, i)] == w[(i, j)])
        for i in range(4):
            check("B.2 out-mass a'_i", sum(theta[(i, j)] for j in range(4) if j != i) == ap[i])
            check("B.2 in-mass b'_i", sum(theta[(j, i)] for j in range(4) if j != i) == bp[i])
    mu = {}
    for (i, j), t in theta.items():
        if t == 0:
            continue
        rest = [k for k in range(4) if k not in (i, j)]
        for c in ((2, 3), (3, 2)):
            sg = [None] * 4; sg[i] = 0; sg[j] = 1; sg[rest[0]] = c[0]; sg[rest[1]] = c[1]
            sg = tuple(sg); mu[sg] = mu.get(sg, Fr(0)) + t / 2
    return mu


def lhs_general(U, Da, De):
    Uh = ct(U)
    M = [[sum((U[i][k] * De[k] * Uh[k][j] for k in range(4)), Q2(0)) for j in range(4)] for i in range(4)]
    return det([[(_c(Da[i]) if i == j else Q2(0)) - M[i][j] for j in range(4)] for i in range(4)])


def section_B2_exact(nsamp, rng):
    print("\n== Lemma B.2: exact replay ==", flush=True)
    fam = [("cayley", cayley_unitary(rng)) for _ in range(nsamp)]
    for _ in range(3):
        fam.append(("perm", perm_matrix(rng.choice(PERMS))))
        fam.append(("2+2 block", blockdiag(cayley_unitary(rng, 2), cayley_unitary(rng, 2))))
        fam.append(("real orthogonal", cayley_unitary(rng, real=True)))
        fam.append(("perm*block", mmul(perm_matrix(rng.choice(PERMS)), blockdiag([[Q2(1)]], cayley_unitary(rng, 3)))))
    # grid: a in {0,1}^4, alpha,beta in {0,1}, gamma in {0..4}; degrees (1,1,1,1,1,1,<=4) -> grid proves identity
    grid = [(g[:4], g[4], g[5], gm) for g in itertools.product((0, 1), repeat=6) for gm in range(5)]
    for kind, U in fam:
        p = [U[i][0] for i in range(4)]; q = [U[i][1] for i in range(4)]
        mu = b2_weights_exact(p, q)
        check("B.2 mu prob", all(w >= 0 for w in mu.values()) and sum(mu.values()) == 1, kind)
        for a, al, be, gm in grid:
            e = [al, be, gm, gm]
            l = lhs_general(U, [Q2(x) for x in a], [Q2(x) for x in e])
            r = rhs_mu(mu, [Q2(x) for x in a], [Q2(x) for x in e])
            if not l == r:
                check("Lemma B.2 identity on proving grid", False, kind); break
    # coincident a eigenvalues: det(diag(al,be,ga,ga) - U diag(e) U^*) via U^* and inverse perms
    for kind, U in fam[:max(3, nsamp // 3)]:
        Uh = ct(U)
        p = [Uh[i][0] for i in range(4)]; q = [Uh[i][1] for i in range(4)]
        nu = b2_weights_exact(p, q, check_it=False)
        mu = {}
        for sg, w in nu.items():
            mu[inv(sg)] = mu.get(inv(sg), Fr(0)) + w
        for e4, al, be, gm in grid:
            a = [al, be, gm, gm]
            if not lhs_general(U, [Q2(x) for x in a], [Q2(x) for x in e4]) == rhs_mu(mu, [Q2(x) for x in a], [Q2(x) for x in e4]):
                check("Lemma B.2 (coincident a, via U^* and inverse permutations)", False, kind); break
    print(f"  samples: {nsamp} Cayley + {len(fam) - nsamp} degenerate (e-coincidence); {max(3, nsamp // 3)} (a-coincidence); grid of 320 points each", flush=True)

    # symbolic: rank-two principal-minor formula and expectation identity
    names = "al,be,ga," + ",".join(f"p{i},pb{i},q{i},qb{i}" for i in range(4))
    R, al, be, ga, *rest = ring(names, QQ)
    P = rest[0::4]; PB = rest[1::4]; Qv = rest[2::4]; QB = rest[3::4]

    def rdet(X):
        n = len(X); tot = R(0)
        for pp in itertools.permutations(range(n)):
            t = R(sgn(tuple(pp) + tuple(range(n, 4))) if n < 4 else sgn(pp))
            for i in range(n):
                t *= X[i][pp[i]]
            tot += t
        return tot
    ok = True
    for k in range(1, 5):
        for I in itertools.combinations(range(4), k):
            Mx = [[(ga if i == j else R(0)) + (al - ga) * P[i] * PB[j] + (be - ga) * Qv[i] * QB[j] for j in I] for i in I]
            lhs = rdet(Mx)
            sa = sum(P[i] * PB[i] for i in I); sb = sum(Qv[i] * QB[i] for i in I)
            sw = sum((P[i] * Qv[j] - P[j] * Qv[i]) * (PB[i] * QB[j] - PB[j] * QB[i]) for i, j in itertools.combinations(I, 2))
            if k == 1:
                rhs = ga + (al - ga) * sa + (be - ga) * sb
            else:
                rhs = ga ** k + ga ** (k - 1) * ((al - ga) * sa + (be - ga) * sb) + ga ** (k - 2) * (al - ga) * (be - ga) * sw
            ok &= (lhs == rhs)
    check("p41 rank-two principal-minor formula (generic p,q, conj as independent symbols)", ok)
    print("  rank-two principal-minor formula, all 15 index sets, generic p,q:", ok, flush=True)
    R4, al, be, ga, PA, PBt, PBoth = ring("al,be,ga,PA,PB,PBoth", QQ)
    ok = True
    for k in range(2, 5):
        E = ga ** k * (1 - PA - PBt + PBoth) + al * ga ** (k - 1) * (PA - PBoth) + be * ga ** (k - 1) * (PBt - PBoth) + al * be * ga ** (k - 2) * PBoth
        F = ga ** k + ga ** (k - 1) * ((al - ga) * PA + (be - ga) * PBt) + ga ** (k - 2) * (al - ga) * (be - ga) * PBoth
        ok &= (E == F)
    check("p41 expectation of site product equals rank-two formula", ok)
    print("  expectation-of-placement identity (k=2..4):", ok, flush=True)


# ---------------------------------------------------------------- Corollary B.3
def relabel_B1_weights(U, i, j):
    """Weights for det(D_a - U D_e U^*) valid whenever a_i = e_j."""
    al = (i,) + tuple(k for k in range(4) if k != i)
    bt = (j,) + tuple(k for k in range(4) if k != j)
    U2 = [[U[al[k]][bt[l]] for l in range(4)] for k in range(4)]
    mu2, _ = b1_weights_exact(U2)
    out = {}
    for s2, w in mu2.items():
        sg = compose(bt, compose(s2, inv(al)))
        out[sg] = out.get(sg, Fr(0)) + w
    return out


def relabel_B2_C(U, i, j):
    """Weights valid whenever e_i = e_j."""
    others = [k for k in range(4) if k not in (i, j)]
    bt = (others[0], others[1], i, j)
    p = [U[r][bt[0]] for r in range(4)]; q = [U[r][bt[1]] for r in range(4)]
    mu2 = b2_weights_exact(p, q, check_it=False)
    out = {}
    for s2, w in mu2.items():
        sg = compose(bt, s2); out[sg] = out.get(sg, Fr(0)) + w
    return out


def relabel_B2_A(U, i, j):
    """Weights valid whenever a_i = a_j: det(D_a - U D_e U^*) = det(D_e - U^* D_a U)."""
    nu = relabel_B2_C(ct(U), i, j)
    out = {}
    for sg, w in nu.items():
        out[inv(sg)] = out.get(inv(sg), Fr(0)) + w
    return out


def section_B3_exact(nO, rng):
    print("\n== Corollary B.3: exact divisibility (B.9) for all 28 brackets ==", flush=True)
    names = ",".join(f"A{k}x,A{k}y" for k in range(4)) + "," + ",".join(f"C{k}x,C{k}y" for k in range(4))
    R, *g = ring(names, QQ)
    A = [(g[2 * k], g[2 * k + 1]) for k in range(4)]; C = [(g[8 + 2 * k], g[8 + 2 * k + 1]) for k in range(4)]
    br = lambda v, w: v[0] * w[1] - v[1] * w[0]
    Pi = {sg: br(A[0], C[sg[0]]) * br(A[1], C[sg[1]]) * br(A[2], C[sg[2]]) * br(A[3], C[sg[3]]) for sg in PERMS}
    Ra, *h = ring("a0,a1,a2,a3,e0,e1,e2,e3", QQ)
    aa = h[:4]; ee = h[4:]
    for t in range(nO):
        U = cayley_unitary(rng, real=True)
        O = [[QQ(U[i][j].r.numerator, U[i][j].r.denominator) for j in range(4)] for i in range(4)]
        check("O real", all(U[i][j].i == 0 for i in range(4) for j in range(4)))
        check("O orthogonal", is_unitary(U))
        check("det O = 1", det(U) == Q2(1))
        Delta = sum((R(sgn(sg)) * Pi[sg] * O[0][sg[0]] * O[1][sg[1]] * O[2][sg[2]] * O[3][sg[3]] for sg in PERMS), R(0))
        # chart identity: Leibniz sum = det(diag(a) - O diag(e) O^T)
        lhs = sum((Ra(sgn(sg)) * O[0][sg[0]] * O[1][sg[1]] * O[2][sg[2]] * O[3][sg[3]] * (aa[0] - ee[sg[0]]) * (aa[1] - ee[sg[1]]) * (aa[2] - ee[sg[2]]) * (aa[3] - ee[sg[3]]) for sg in PERMS), Ra(0))
        Mx = [[(aa[i] if i == j else Ra(0)) - sum(O[i][k] * ee[k] * O[j][k] for k in range(4)) for j in range(4)] for i in range(4)]
        dd = Ra(0)
        for sg in PERMS:
            tt = Ra(sgn(sg))
            for i in range(4):
                tt *= Mx[i][sg[i]]
            dd += tt
        check("B.3 chart: Delta_O = det(diag(a) - O diag(e) O^T)", lhs == dd)
        ndiv = 0; nnon = 0
        pairs = [("AC", i, j) for i in range(4) for j in range(4)] + [("AA", i, j) for i, j in EDGES] + [("CC", i, j) for i, j in EDGES]
        for typ, i, j in pairs:
            if typ == "AC":
                mu = relabel_B1_weights(U, i, j); d = br(A[i], C[j])
            elif typ == "CC":
                mu = relabel_B2_C(U, i, j); d = br(C[i], C[j])
            else:
                mu = relabel_B2_A(U, i, j); d = br(A[i], A[j])
            check("B.3 weights probability", all(w >= 0 for w in mu.values()) and sum(mu.values()) == 1, f"{typ}{i}{j}")
            F = Delta - sum((Pi[sg] * QQ(w.numerator, w.denominator) for sg, w in mu.items()), R(0))
            if F == 0:
                print(f"    note: Delta_O equals the weighted average identically for bracket {typ}{i}{j}", flush=True)
            qs, r = F.div([d])
            qq = qs[0] if qs else R(0)
            ok = (r == 0) and (qq * d == F)
            check("(B.9) divisibility", ok, f"{typ}{i}{j}")
            ndiv += ok
            # non-vacuity: without the common weights divisibility fails
            _, r2 = (Delta - Pi[ID]).div([d])
            nnon += (r2 != 0)
        print(f"  O #{t}: (B.9) divisible for {ndiv}/28 brackets; Delta_O - Pi_id NOT divisible for {nnon}/28 (control)", flush=True)


# ---------------------------------------------------------------- NUMERIC sections
def haar(rng, n=4):
    Z = (rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))) / np.sqrt(2)
    Q, Rm = np.linalg.qr(Z)
    return Q * (np.diag(Rm) / np.abs(np.diag(Rm)))


def haar_real(rng, n=4):
    Q, Rm = np.linalg.qr(rng.normal(size=(n, n)))
    Q = Q * np.sign(np.diag(Rm))
    if np.linalg.det(Q) < 0:
        Q[:, 0] *= -1
    return Q


GRID6 = np.array(list(itertools.product((0, 1), repeat=6)), dtype=float)
RHS6 = np.array([[np.prod([(np.r_[0, g[:3]][i] - np.r_[0, g[3:]][sg[i]]) for i in range(4)]) for sg in PERMS] for g in GRID6])


def lhs_grid6(U):
    Mats = np.zeros((64, 4, 4), dtype=complex)
    for t, g in enumerate(GRID6):
        Da = np.diag(np.r_[0, g[:3]]); De = np.diag(np.r_[0, g[3:]])
        Mats[t] = Da - U @ De @ U.conj().T
    return np.linalg.det(Mats)


def b1_weights_float(U, tol=1e-13):
    s = 1 - abs(U[0, 0]) ** 2
    mu = np.zeros(24)
    if s < 1e-12:
        mu[PIDX[ID]] = 1; return mu, 0
    W = np.zeros((4, 4))
    W[0, 1:] = abs(U[0, 1:]) ** 2; W[1:, 0] = abs(U[1:, 0]) ** 2
    for i in (1, 2, 3):
        for j in (1, 2, 3):
            cols = [k for k in (1, 2, 3) if k != j]
            W[i, j] = abs(np.linalg.det(U[np.ix_([0, i], cols)])) ** 2
    W /= s
    terms = []
    for _ in range(40):
        if W.max() < 1e-12:
            break
        best = None
        for p in PERMS:
            v = min(W[i, p[i]] for i in range(4))
            if v > tol and (best is None or v > best[1]):
                best = (p, v)
        if best is None:
            break
        p, v = best; terms.append([p, v])
        for i in range(4):
            W[i, p[i]] -= v
        W[W < tol] = 0
    while len(terms) > 9:  # Caratheodory reduction in floats
        Mx = np.array([[1] + pvec(p) for p, _ in terms], dtype=float).T
        _, _, Vt = np.linalg.svd(Mx); lam = Vt[-1]
        if lam.max() <= 0:
            lam = -lam
        tt = min(w / l for (p, w), l in zip(terms, lam) if l > 1e-12)
        terms = [[p, w - tt * l] for (p, w), l in zip(terms, lam)]
        terms = [[p, w] for p, w in terms if w > 1e-14]
    for p, w in terms:
        mu[PIDX[phi(p)]] += s * w
    mu[PIDX[ID]] += 1 - s
    return mu, len(terms)


def degenerate_family(rng):
    out = []
    for p in PERMS:
        P = np.zeros((4, 4), complex)
        for i in range(4):
            P[i, p[i]] = np.exp(1j * rng.uniform(0, 2 * np.pi))
        out.append(("perm", P))
    for _ in range(40):
        B = np.zeros((4, 4), complex); B[:2, :2] = haar(rng, 2); B[2:, 2:] = haar(rng, 2)
        out.append(("2+2", B))
        out.append(("perm*2+2*perm", np.eye(4)[list(rng.permutation(4))] @ B @ np.eye(4)[list(rng.permutation(4))]))
        C = np.zeros((4, 4), complex); C[0, 0] = np.exp(1j * rng.uniform(0, 6.3)); C[1:, 1:] = haar(rng, 3)
        out.append(("u00 unimodular", C))
        out.append(("u00=0", np.eye(4)[[1, 0, 2, 3]] @ C))
        eps = 10.0 ** rng.uniform(-9, -2)
        K = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4)); K = K - K.conj().T
        from scipy.linalg import expm
        out.append(("near s=0", C @ expm(eps * K)))
        out.append(("near permutation", out[rng.integers(0, 24)][1] @ expm(eps * K)))
        out.append(("real orthogonal", haar_real(rng).astype(complex)))
    F = np.array([[np.exp(2j * np.pi * j * k / 4) for k in range(4)] for j in range(4)]) / 2
    out.append(("DFT", F))
    return out


def lp_feasible(Amat, bvec):
    """min ||residual||_inf  s.t. mu>=0, sum mu = 1 ;  Amat mu ~ bvec."""
    m, n = Amat.shape
    # variables mu (n), t (1): minimize t, -t <= A mu - b <= t
    c = np.r_[np.zeros(n), 1.0]
    Aub = np.block([[Amat, -np.ones((m, 1))], [-Amat, -np.ones((m, 1))]])
    bub = np.r_[bvec, -bvec]
    Aeq = np.r_[np.ones(n), 0.0][None, :]
    res = linprog(c, A_ub=Aub, b_ub=bub, A_eq=Aeq, b_eq=[1.0], bounds=[(0, None)] * n + [(0, None)], method="highs")
    return res


def section_B1_numeric(n_constr, n_lp, rng):
    print("\n== Theorem B.1: NUMERIC (floats) ==", flush=True)
    fam = [("haar", haar(rng)) for _ in range(n_constr)] + degenerate_family(rng)
    worst = 0.0; worstB2 = 0.0; maxterms = 0; maxnz = 0
    for kind, U in fam:
        mu, nt = b1_weights_float(U)
        maxterms = max(maxterms, nt); maxnz = max(maxnz, int((mu > 1e-12).sum()))
        L = lhs_grid6(U); Rr = RHS6 @ mu
        scale = 1 + np.abs(L).max()
        err = np.abs(L - Rr).max() / scale
        worst = max(worst, err)
        worstB2 = max(worstB2, abs(sum(mu[PIDX[p]] for p in PERMS if p[0] == 0) - abs(U[0, 0]) ** 2))
        if err > 1e-8 or mu.min() < -1e-12:
            anomaly("B.1 float construction residual", f"{kind} err={err:.2e}")
    print(f"  paper's construction in floats: {n_constr} Haar + {len(fam) - n_constr} degenerate U; "
          f"max rel. residual on {{0,1}}^6 grid {worst:.2e}; max |(B.2) defect| {worstB2:.2e}; "
          f"max Birkhoff terms {maxterms}; max nonzero weights {maxnz}", flush=True)
    # independent LP feasibility (does not use the paper's reduction to the nine m_ij)
    fam2 = [("haar", haar(rng)) for _ in range(n_lp)] + degenerate_family(rng)
    worst = 0.0; maxnz = 0; worstB2 = 0.0; nfail = 0
    for kind, U in fam2:
        L = lhs_grid6(U)
        Amat = np.r_[RHS6, np.zeros_like(RHS6)]
        bvec = np.r_[L.real, L.imag]
        res = lp_feasible(Amat, bvec)
        if res.status != 0:
            nfail += 1; anomaly("B.1 LP status", f"{kind} status={res.status}"); continue
        mu = res.x[:24]; worst = max(worst, res.x[24] / (1 + np.abs(L).max()))
        maxnz = max(maxnz, int((mu > 1e-9).sum()))
        worstB2 = max(worstB2, abs(sum(mu[PIDX[p]] for p in PERMS if p[0] == 0) - abs(U[0, 0]) ** 2))
        if res.x[24] / (1 + np.abs(L).max()) > 1e-7:
            anomaly("B.1 LP infeasible", f"{kind} resid={res.x[24]:.2e}")
    print(f"  direct LP (min sup-residual of (B.1) on the 64-point proving grid): {n_lp} Haar + {len(fam2) - n_lp} degenerate U; "
          f"max min-residual {worst:.2e}; LP failures {nfail}; max LP-vertex support {maxnz}; "
          f"max |(B.2) defect| of LP solutions {worstB2:.2e} (B.2 is forced, see report)", flush=True)


def section_B2_numeric(n_lp, rng):
    print("\n== Lemma B.2: NUMERIC direct LP ==", flush=True)
    grid = [(g[:4], g[4], g[5], gm) for g in itertools.product((0, 1), repeat=6) for gm in range(5)]
    RHS = np.array([[np.prod([a[i] - [al, be, gm, gm][sg[i]] for i in range(4)]) for sg in PERMS] for a, al, be, gm in grid], dtype=float)
    fam = [("haar", haar(rng)) for _ in range(n_lp)] + degenerate_family(rng)[::3]
    worst = 0.0; nfail = 0
    for kind, U in fam:
        L = np.array([np.linalg.det(np.diag(a) - U @ np.diag([al, be, gm, gm]) @ U.conj().T) for a, al, be, gm in grid])
        res = lp_feasible(np.r_[RHS, np.zeros_like(RHS)], np.r_[L.real, L.imag])
        if res.status != 0:
            nfail += 1; anomaly("B.2 LP status", kind); continue
        r = res.x[24] / (1 + np.abs(L).max()); worst = max(worst, r)
        if r > 1e-7:
            anomaly("B.2 LP infeasible", f"{kind} resid={r:.2e}")
    print(f"  {n_lp} Haar + {len(fam) - n_lp} degenerate U, 320-point proving grid: max min-residual {worst:.2e}, LP failures {nfail}", flush=True)


def section_B3_numeric(ninst, rng):
    print("\n== Corollary B.3: NUMERIC curves t -> 0 (single vs. double coincidence) ==", flush=True)
    print("  real spectra (so the region {all Re g >= 0} is a nonempty sector), complex directions, tau = argmin Pi_sigma(0),", flush=True)
    print("  t = r exp(i theta), theta on a 6000-point grid plus 4001 points within 200 r of every boundary ray", flush=True)
    radii = [1e-2, 1e-3, 1e-4, 1e-5, 1e-6]
    th0 = np.linspace(0, 2 * np.pi, 6000, endpoint=False)
    u = np.linspace(-200, 200, 4001)
    sgns = np.array([sgn(p) for p in PERMS], float)
    summary = {}
    for mode in ("single", "double"):
        grow = []; nonempty = 0; examples = []; remainder = []
        for it in range(ninst):
            O = haar_real(rng)
            a = rng.normal(size=4).astype(complex); e = rng.normal(size=4).astype(complex)
            e[0] = a[0]
            if mode == "double":
                e[1] = a[1]
            da = rng.normal(size=4) + 1j * rng.normal(size=4); de = rng.normal(size=4) + 1j * rng.normal(size=4)
            Ofac = np.array([np.prod([O[k, p[k]] for k in range(4)]) for p in PERMS])
            mu, _ = b1_weights_float(O.astype(complex))  # Theorem B.1 law for the pair a0 = e0

            def Pi(t):
                A = a[:, None] + np.outer(da, t); C = e[:, None] + np.outer(de, t)
                return np.array([np.prod([A[k] - C[p[k]] for k in range(4)], axis=0) for p in PERMS])
            Pi0 = Pi(np.array([0.0]))[:, 0]
            tau = int(np.argmin(Pi0.real))
            eta = 1.0
            c = (da[0] - de[0]) * ((da[1] - de[1]) if mode == "double" else 1)
            k = 2 if mode == "double" else 1
            rays = [((np.angle(eta / c) + sgn_) / k) + 2 * np.pi * j / k for sgn_ in (np.pi / 2, -np.pi / 2) for j in range(k)]
            mins = []; rem = []
            for r in radii:
                th = np.concatenate([th0] + [ray + r * u for ray in rays])
                t = r * np.exp(1j * th)
                P = Pi(t)
                E = (a[0] + da[0] * t - e[0] - de[0] * t)
                if mode == "double":
                    E = E * (a[1] + da[1] * t - e[1] - de[1] * t)
                g = (eta * (P - P[tau]) / E).real
                Dl = (sgns * Ofac) @ P
                L = (eta * (Dl - P[tau]) / E).real
                mask = (g >= -1e-9 * (1 + np.abs(g))).all(axis=0)
                mins.append(L[mask].min() if mask.any() else np.nan)
                rem.append(np.abs(L - mu @ g).max())
            remainder.append(rem)
            if not np.isnan(mins[1]) and not np.isnan(mins[-1]):
                nonempty += 1
                neg_small = max(0.0, -mins[-1]); neg_big = max(1.0, -mins[1])
                grow.append(neg_small / neg_big)
                if len(examples) < 2 or (len(examples) < 4 and neg_small / neg_big > 100):
                    examples.append([f"{m:.3g}" for m in mins])
        grow = np.array(grow); remainder = np.array(remainder)
        summary[mode] = (nonempty, grow, remainder)
        print(f"  {mode}: {ninst} random instances; region nonempty at r=1e-3 and 1e-6 in {nonempty}; "
              f"growth ratio (neg. part of min L at r=1e-6)/max(1, neg. part at r=1e-3): max {grow.max() if len(grow) else float('nan'):.3g}, "
              f"instances with ratio > 100: {(grow > 100).sum() if len(grow) else 0}", flush=True)
        print(f"    max over instances of sup|L - sum mu Re g| (B.1 law for a0=e0) at r={radii}: "
              f"{[f'{x:.3g}' for x in remainder.max(axis=0)]}", flush=True)
        for ex in examples:
            print(f"    example: min L over region at r={radii}: {ex}", flush=True)
    ns, gs, rs = summary["single"]
    if len(gs) and gs.max() > 10:
        anomaly("B.3 single-coincidence negative part grows", f"{gs.max():.3g}")
    if rs[:, -1].max() > 10 * rs[:, 0].max() + 1:
        anomaly("B.3 single-coincidence remainder grows")


def main():
    rng = random.Random(20261125)
    nrng = np.random.default_rng(20261125)
    q = "--quick" in sys.argv  # smoke-test sizes
    section_B1_symbolic()
    section_B1_exact(4 if q else 60, rng)
    section_B2_exact(2 if q else 30, rng)
    section_B3_exact(1 if q else 2, rng)
    print(f"\n[exact sections done at {time.time() - T0:.1f}s]", flush=True)
    section_B1_numeric(200 if q else 20000, 50 if q else 3000, nrng)
    section_B2_numeric(20 if q else 600, nrng)
    section_B3_numeric(20 if q else 300, nrng)
    print(f"\nTotal time {time.time() - T0:.1f}s")
    print("EXACT failures:", EXACT_FAIL if EXACT_FAIL else "none")
    print("NUMERIC anomalies:", sorted(set(NUM_ANOM)) if NUM_ANOM else "none")
    sys.exit(1 if EXACT_FAIL else (3 if NUM_ANOM else 0))


if __name__ == "__main__":
    main()
