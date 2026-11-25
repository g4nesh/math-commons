#!/usr/bin/env python3
"""Exact Q(i) checks of the finite identities in reductions/README.md (commit 67528c1).

Review aid for mdo-reductions, written by vn-claude independently of the author's code.
Each trial uses small random Gaussian-integer data and exact Fractions, so a pass is an
exact identity check on that instance, not a proof of the general statement.

Checks per trial:
  Lemma 1   K skew-Hermitian => I+K and I-K invertible, C(K) unitary, entries in Q(i)
  Lemma 2   U unitary without eigenvalue -1 => K=(I-U)(I+U)^-1 skew-Hermitian, C(K)=U
  §3        reciprocal transform: det(A'+B') * det X * det Y = det(A+B), and for every
            permutation the transformed vertex times the same denominator is z_sigma
  §3        LUR invariance and swap symmetry, with exact rational unitaries

Bounds: n in 2..4, 60 trials, seed 20261125, entries in [-3, 3].
Usage: python -I research/marcus-de-oliveira/reviews/check_reductions_exact.py
"""
import random
from fractions import Fraction as F
from itertools import permutations

SEED, TRIALS = 20261125, 60
rng = random.Random(SEED)


class Q:
    __slots__ = ("r", "i")

    def __init__(self, r=0, i=0):
        self.r, self.i = F(r), F(i)

    def __add__(s, o): return Q(s.r + o.r, s.i + o.i)
    def __sub__(s, o): return Q(s.r - o.r, s.i - o.i)
    def __neg__(s): return Q(-s.r, -s.i)
    def __mul__(s, o): return Q(s.r * o.r - s.i * o.i, s.r * o.i + s.i * o.r)
    def conj(s): return Q(s.r, -s.i)
    def __eq__(s, o): return s.r == o.r and s.i == o.i
    def iszero(s): return s.r == 0 and s.i == 0

    def __truediv__(s, o):
        d = o.r * o.r + o.i * o.i
        n = s * o.conj()
        return Q(n.r / d, n.i / d)


Z, O = Q(0), Q(1)


def eye(n): return [[O if i == j else Z for j in range(n)] for i in range(n)]
def add(A, B): return [[a + b for a, b in zip(ra, rb)] for ra, rb in zip(A, B)]
def sub(A, B): return [[a - b for a, b in zip(ra, rb)] for ra, rb in zip(A, B)]
def adj(A): return [[A[j][i].conj() for j in range(len(A))] for i in range(len(A))]
def diag(v): return [[v[i] if i == j else Z for j in range(len(v))] for i in range(len(v))]
def scal(c, A): return [[c * a for a in r] for r in A]


def mul(A, B):
    n, m = len(A), len(B[0])
    out = []
    for i in range(n):
        row = []
        for j in range(m):
            s = Z
            for k in range(len(B)):
                s = s + A[i][k] * B[k][j]
            row.append(s)
        out.append(row)
    return out


def inv_det(A):
    """Exact Gauss-Jordan. Returns (inverse or None, determinant)."""
    n = len(A)
    M = [list(r) + e for r, e in zip(A, eye(n))]
    det = O
    for c in range(n):
        p = next((r for r in range(c, n) if not M[r][c].iszero()), None)
        if p is None:
            return None, Z
        if p != c:
            M[c], M[p] = M[p], M[c]
            det = -det
        piv = M[c][c]
        det = det * piv
        M[c] = [x / piv for x in M[c]]
        for r in range(n):
            if r != c and not M[r][c].iszero():
                f = M[r][c]
                M[r] = [x - f * y for x, y in zip(M[r], M[c])]
    return [row[n:] for row in M], det


def gi(): return Q(rng.randint(-3, 3), rng.randint(-3, 3))


def skew(n):
    K = [[Z] * n for _ in range(n)]
    for i in range(n):
        K[i][i] = Q(0, rng.randint(-3, 3))
        for j in range(i + 1, n):
            K[i][j] = gi()
            K[j][i] = -K[i][j].conj()
    return K


def cayley(K):
    n = len(K)
    inv, _ = inv_det(add(eye(n), K))
    return mul(sub(eye(n), K), inv)


def d(a, b, U):
    return inv_det(add(diag(a), mul(mul(U, diag(b)), adj(U))))[1]


fails, counts = [], {}
def check(name, ok):
    counts[name] = counts.get(name, 0) + 1
    if not ok:
        fails.append(name)


for _ in range(TRIALS):
    n = rng.randint(2, 4)
    I = eye(n)
    K = skew(n)
    check("Lemma1: I+K invertible", inv_det(add(I, K))[0] is not None)
    check("Lemma1: I-K invertible (used, not stated)", inv_det(sub(I, K))[0] is not None)
    U = cayley(K)
    check("Lemma1: C(K) unitary", mul(adj(U), U) == I)

    # Lemma 2 round trip on a rational unitary U without eigenvalue -1 (I+U invertible)
    invIU, _ = inv_det(add(I, U))
    if invIU is not None:
        K2 = mul(sub(I, U), invIU)
        check("Lemma2: K skew-Hermitian", add(K2, adj(K2)) == [[Z] * n for _ in range(n)])
        check("Lemma2: C(K)=U", cayley(K2) == U)

    a = [gi() for _ in range(n)]
    b = [gi() for _ in range(n)]
    base = d(a, b, U)

    # LUR with exact unit-modulus Gaussian-rational phases (3+4i)/5 etc.
    phases = [Q(1), Q(F(3, 5), F(4, 5)), Q(0, 1), Q(F(-5, 13), F(12, 13))]
    L = diag([rng.choice(phases) for _ in range(n)])
    R = diag([rng.choice(phases) for _ in range(n)])
    check("S3: LUR invariance", d(a, b, mul(mul(L, U), R)) == base)
    check("S3: swap a,b and U->U*", d(b, a, adj(U)) == base)

    # reciprocal transform with a pole-avoiding Gaussian-rational c
    while True:
        c = Q(F(rng.randint(-7, 7), 2), F(rng.randint(-7, 7), 3))
        if all(not (x - c).iszero() for x in a) and all(not (y + c).iszero() for y in b):
            break
    A = diag(a)
    B = mul(mul(U, diag(b)), adj(U))
    X, Y = sub(A, scal(c, I)), add(B, scal(c, I))
    Xi, dX = inv_det(X)
    Yi, dY = inv_det(Y)
    lhs = inv_det(add(Xi, Yi))[1] * dX * dY
    check("S3: reciprocal determinant identity", lhs == base)
    den = O
    for x in a:
        den = den * (x - c)
    for y in b:
        den = den * (y + c)
    ok = True
    for s in permutations(range(n)):
        zp, z = O, O
        for i in range(n):
            zp = zp * (O / (a[i] - c) + O / (b[s[i]] + c))
            z = z * (a[i] + b[s[i]])
        ok &= (zp * den == z)
    check("S3: reciprocal vertex scaling (all sigma)", ok)

for k, v in counts.items():
    print(f"{k:45s} {v:3d} instances")
print(f"seed={SEED} trials={TRIALS}: " + ("all exact checks passed" if not fails else f"FAILED: {sorted(set(fails))}"))
raise SystemExit(1 if fails else 0)
