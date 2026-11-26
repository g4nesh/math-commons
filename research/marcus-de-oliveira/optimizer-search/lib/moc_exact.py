"""Exact certification over Q(i) (fractions only).

Builder: float candidate (a, b, U) -> rational skew-Hermitian K with
U diag(R) ~ cayley(K), exact U = (I - K)(I + K)^{-1}, Gaussian-rational spectra,
exact d and all n! z_sigma, exact hull membership, and (if outside) a rational
direction (p, q) with  p Re d + q Im d  >  max_sigma (p Re z_sigma + q Im z_sigma).
The certificate JSON is checked by check_cert.py, which shares no code with this file.
"""
import itertools
import json
import math
from fractions import Fraction as Fr

import numpy as np


class G:
    """Gaussian rational x + i y."""
    __slots__ = ("x", "y")

    def __init__(self, x, y=0):
        self.x = Fr(x)
        self.y = Fr(y)

    def __add__(s, o):
        o = _g(o)
        return G(s.x + o.x, s.y + o.y)
    __radd__ = __add__

    def __sub__(s, o):
        o = _g(o)
        return G(s.x - o.x, s.y - o.y)

    def __rsub__(s, o):
        return _g(o) - s

    def __neg__(s):
        return G(-s.x, -s.y)

    def __mul__(s, o):
        o = _g(o)
        return G(s.x * o.x - s.y * o.y, s.x * o.y + s.y * o.x)
    __rmul__ = __mul__

    def __truediv__(s, o):
        o = _g(o)
        n = o.x * o.x + o.y * o.y
        return G((s.x * o.x + s.y * o.y) / n, (s.y * o.x - s.x * o.y) / n)

    def conj(s):
        return G(s.x, -s.y)

    def iszero(s):
        return s.x == 0 and s.y == 0

    def __eq__(s, o):
        o = _g(o)
        return s.x == o.x and s.y == o.y

    def enc(s):
        return [str(s.x), str(s.y)]

    def __complex__(s):
        return complex(float(s.x), float(s.y))


def _g(v):
    return v if isinstance(v, G) else G(v)


def rat(v, den):
    v = complex(v)
    return G(Fr(v.real).limit_denominator(den), Fr(v.imag).limit_denominator(den))


def mat_mul(A, B):
    n, k, m = len(A), len(B), len(B[0])
    return [[sum((A[i][t] * B[t][j] for t in range(k)), G(0)) for j in range(m)] for i in range(n)]


def adjoint(A):
    return [[A[j][i].conj() for j in range(len(A))] for i in range(len(A[0]))]


def solve(A, B):
    """A^{-1} B by Gauss-Jordan with nonzero pivoting."""
    n = len(A)
    M = [list(A[i]) + list(B[i]) for i in range(n)]
    w = len(M[0])
    for c in range(n):
        p = next(r for r in range(c, n) if not M[r][c].iszero())
        M[c], M[p] = M[p], M[c]
        piv = M[c][c]
        M[c] = [v / piv for v in M[c]]
        for r in range(n):
            if r != c and not M[r][c].iszero():
                f = M[r][c]
                M[r] = [M[r][j] - f * M[c][j] for j in range(w)]
    return [row[n:] for row in M]


def det(A):
    n = len(A)
    M = [list(r) for r in A]
    out = G(1)
    for c in range(n):
        p = next((r for r in range(c, n) if not M[r][c].iszero()), None)
        if p is None:
            return G(0)
        if p != c:
            M[c], M[p] = M[p], M[c]
            out = -out
        out = out * M[c][c]
        for r in range(c + 1, n):
            if not M[r][c].iszero():
                f = M[r][c] / M[c][c]
                M[r] = [M[r][j] - f * M[c][j] for j in range(n)]
    return out


def eye(n):
    return [[G(1) if i == j else G(0) for j in range(n)] for i in range(n)]


def exact_unitary_from_float(U, den=10 ** 12):
    """Rational skew-Hermitian K with cayley(K) ~ U diag(R); returns (K, Uex, R)."""
    from moc_param import cayley_coords
    Kf, R, gap = cayley_coords(np.asarray(U, complex))
    n = len(Kf)
    Xf = Kf.real + Kf.imag
    X = [[Fr(float(Xf[i, j])).limit_denominator(den) for j in range(n)] for i in range(n)]
    K = [[G((X[i][j] - X[j][i]) / 2, (X[i][j] + X[j][i]) / 2) for j in range(n)] for i in range(n)]
    I = eye(n)
    IpK = [[I[i][j] + K[i][j] for j in range(n)] for i in range(n)]
    ImK = [[I[i][j] - K[i][j] for j in range(n)] for i in range(n)]
    return K, solve(IpK, ImK), R


def is_unitary(U):
    n = len(U)
    P = mat_mul(U, adjoint(U))
    return all(P[i][j] == (G(1) if i == j else G(0)) for i in range(n) for j in range(n))


def d_exact(a, b, U):
    n = len(a)
    Ub = [[U[i][k] * b[k] for k in range(n)] for i in range(n)]
    B = mat_mul(Ub, adjoint(U))
    M = [[(a[i] if i == j else G(0)) + B[i][j] for j in range(n)] for i in range(n)]
    return det(M)


def z_exact(a, b):
    n = len(a)
    out = []
    for s in itertools.permutations(range(n)):
        z = G(1)
        for i in range(n):
            z = z * (a[i] + b[s[i]])
        out.append((s, z))
    return out


def membership(d, zs):
    """'inside' (interior), 'boundary', or 'outside' for the exact hull."""
    P = sorted(set((z.x, z.y) for _, z in zs))
    def cross(o, p, q):
        return (p[0] - o[0]) * (q[1] - o[1]) - (p[1] - o[1]) * (q[0] - o[0])
    lo, up = [], []
    for p in P:
        while len(lo) >= 2 and cross(lo[-2], lo[-1], p) <= 0:
            lo.pop()
        lo.append(p)
    for p in reversed(P):
        while len(up) >= 2 and cross(up[-2], up[-1], p) <= 0:
            up.pop()
        up.append(p)
    V = lo[:-1] + up[:-1] if len(P) > 1 else P
    q = (d.x, d.y)
    if len(V) == 1:
        return "boundary" if q == V[0] else "outside"
    if len(V) == 2:
        if cross(V[0], V[1], q) != 0:
            return "outside"
        lo_, hi_ = sorted([V[0], V[1]])
        return "boundary" if lo_ <= q <= hi_ else "outside"
    cs = [cross(V[k], V[(k + 1) % len(V)], q) for k in range(len(V))]
    if min(cs) > 0:
        return "inside"
    if min(cs) == 0:
        return "boundary"
    return "outside"


def build_certificate(a, b, U, theta, den=10 ** 12, dir_den=10 ** 9, out_json=None, U_exact=None):
    """Exact certificate attempt. U_exact (list of lists of G) bypasses the Cayley
    construction (used only to plant fakes). Returns dict with status
    'outside-certified' / 'not-certified' and the exact data."""
    n = len(a)
    if U_exact is None:
        K, Uex, R = exact_unitary_from_float(U, den)
        Kenc = [[v.enc() for v in row] for row in K]
    else:
        Uex, Kenc = U_exact, None
    ae = [rat(v, den) for v in a]
    be = [rat(v, den) for v in b]
    d = d_exact(ae, be, Uex)
    zs = z_exact(ae, be)
    memb = membership(d, zs)
    p = Fr(round(math.cos(theta) * dir_den), dir_den)
    q = Fr(round(math.sin(theta) * dir_den), dir_den)
    Ld = p * d.x + q * d.y
    Lz = [p * z.x + q * z.y for _, z in zs]
    mx = max(Lz)
    gap = Ld - mx
    cert = dict(kind="moc-outside-certificate", n=n,
                a=[v.enc() for v in ae], b=[v.enc() for v in be],
                U=[[v.enc() for v in row] for row in Uex], K=Kenc,
                direction=[str(p), str(q)], d=d.enc(), max_L_z=str(mx), gap=str(gap),
                n_perms=len(zs), membership=memb, unitary=is_unitary(Uex),
                status="outside-certified" if (gap > 0 and is_unitary(Uex)) else "not-certified")
    if out_json:
        with open(out_json, "w") as f:
            json.dump(cert, f, indent=1)
    return cert


def witness_json(a, b, U, name):
    """Input file for the repo verifier (certificates/verify.py schema)."""
    return {"schemaVersion": 1, "name": name, "a": [v.enc() for v in a],
            "b": [v.enc() for v in b], "U": [[v.enc() for v in row] for row in U]}
