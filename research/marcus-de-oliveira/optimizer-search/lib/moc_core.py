"""Floating-point evaluators for the Marcus-de Oliveira (MOC) search.

d(U) = det(D_a + U D_b U*), z_sigma = prod_i (a_i + b_sigma(i)), H = conv{z_sigma}.
The normalized signed margin is m = (signed distance from d to H) / diam(H),
positive outside H and negative inside. Permutations are 0-based tuples in
lexicographic order, sigma[i] = image of i.
"""
import itertools
import math
from functools import lru_cache

import numpy as np
from scipy.optimize import linear_sum_assignment

UNITARY_TOL = 1e-10


@lru_cache(maxsize=None)
def perms(n):
    """All n! permutations, lexicographic, shape (n!, n), int8. n <= 8."""
    if n > 8:
        raise ValueError("n > 8: use z_chunks / hull_lean (memory-lean path)")
    P = np.array(list(itertools.permutations(range(n))), dtype=np.int8)
    P.setflags(write=False)
    return P


def z_all(a, b):
    """All n! points z_sigma in perms(n) order (n <= 8)."""
    a = np.asarray(a, complex)
    b = np.asarray(b, complex)
    return (a[None, :] + b[perms(len(a))]).prod(axis=1)


def z_chunks(a, b):
    """Memory-lean enumeration for n = 8, 9: yields (sigmas, z) per value of
    sigma(0), in global lexicographic order; each chunk has (n-1)! rows."""
    a = np.asarray(a, complex)
    b = np.asarray(b, complex)
    n = len(a)
    Q = perms(n - 1)
    for j in range(n):
        rest = np.array([k for k in range(n) if k != j], dtype=np.int8)
        S = np.empty((len(Q), n), dtype=np.int8)
        S[:, 0] = j
        S[:, 1:] = rest[Q]
        z = (a[0] + b[j]) * (a[None, 1:] + b[S[:, 1:]]).prod(axis=1)
        yield S, z


def det_d(a, b, U):
    a = np.asarray(a, complex)
    b = np.asarray(b, complex)
    U = np.asarray(U, complex)
    return np.linalg.det(np.diag(a) + U @ np.diag(b) @ U.conj().T)


def unitarity_defect(U):
    U = np.asarray(U, complex)
    return float(np.abs(U @ U.conj().T - np.eye(len(U))).max())


# ---------------------------------------------------------------- hull geometry

def hull_indices(pts):
    """Strict convex hull (Andrew monotone chain) of complex points.
    Returns indices of hull vertices in counterclockwise order. Collinear
    points are dropped. 1 index if all points coincide, 2 if collinear."""
    pts = np.asarray(pts, complex)
    order = np.lexsort((pts.imag, pts.real))
    xs = pts.real[order].tolist()
    ys = pts.imag[order].tolist()
    idx = order.tolist()

    def cross(o, p, q):
        return (xs[p] - xs[o]) * (ys[q] - ys[o]) - (ys[p] - ys[o]) * (xs[q] - xs[o])

    lower, upper = [], []
    for k in range(len(idx)):
        while len(lower) >= 2 and cross(lower[-2], lower[-1], k) <= 0:
            lower.pop()
        lower.append(k)
    for k in range(len(idx) - 1, -1, -1):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], k) <= 0:
            upper.pop()
        upper.append(k)
    chain = lower[:-1] + upper[:-1]
    out = []
    for k in chain:
        if not out or (xs[k], ys[k]) != (xs[out[-1]], ys[out[-1]]):
            out.append(k)
    if len(out) > 1 and (xs[out[0]], ys[out[0]]) == (xs[out[-1]], ys[out[-1]]):
        out.pop()
    if not out:
        out = [0]
    return [idx[k] for k in out]


def _seg_nearest(p, u, v):
    e = v - u
    L2 = (e.real ** 2 + e.imag ** 2)
    if L2 == 0:
        return u
    t = ((p - u).real * e.real + (p - u).imag * e.imag) / L2
    t = min(1.0, max(0.0, t))
    return u + t * e


def signed_distance(d, V):
    """Signed distance from point d to the convex polygon with CCW vertices V
    (complex array). Returns (dist, theta): theta is a direction attaining
    max_theta [Re(e^{-i theta} d) - h_H(theta)] = dist."""
    V = np.asarray(V, complex)
    if len(V) == 1:
        r = d - V[0]
        return abs(r), (float(np.angle(r)) if r != 0 else 0.0)
    if len(V) == 2:
        q = _seg_nearest(d, V[0], V[1])
        r = d - q
        if abs(r) > 0:
            return abs(r), float(np.angle(r))
        return 0.0, float(np.angle((V[1] - V[0]) * -1j))
    E = np.roll(V, -1) - V
    nrm = -1j * E / np.abs(E)          # outward unit normals for CCW order
    s = (np.conj(nrm) * (d - V)).real  # n . (d - v_k)
    if (s <= 0).all():
        k = int(np.argmax(s))
        return float(s[k]), float(np.angle(nrm[k]))
    best = None
    for k in range(len(V)):
        q = _seg_nearest(d, V[k], V[(k + 1) % len(V)])
        r = abs(d - q)
        if best is None or r < best[0]:
            best = (r, q)
    return float(best[0]), float(np.angle(d - best[1]))


def diameter(V):
    V = np.asarray(V, complex)
    if len(V) < 2:
        return 0.0
    return float(np.abs(V[:, None] - V[None, :]).max())


def margin_points(d, z, n=None, tol_active=1e-9, allow_partial=False):
    """Signed margin of d against the point list z (all n! points!).
    Returns dict(m, dist, diam, theta, active, n_vertices, degenerate)."""
    z = np.asarray(z, complex)
    if n is not None and not allow_partial and len(z) != math.factorial(n):
        raise ValueError("incomplete permutation list: %d != %d!" % (len(z), n))
    hi = hull_indices(z)
    V = z[hi]
    diam = diameter(V)
    dist, theta = signed_distance(complex(d), V)
    proj = (np.exp(-1j * theta) * z).real
    h = proj.max()
    active = np.nonzero(proj >= h - tol_active * max(diam, 1e-300))[0].tolist()
    degenerate = len(V) < 3
    m = dist / diam if diam > 0 else float("nan")
    return dict(m=m, dist=dist, diam=diam, theta=theta, active=active,
                n_vertices=len(V), degenerate=degenerate,
                support_gap=float((np.exp(-1j * theta) * d).real - h))


def hull_lean(a, b):
    """Hull vertices of {z_sigma} for n = 8, 9 without holding all n! points.
    Returns (V, sigmas) with V the CCW hull vertices."""
    keepz, keeps = [], []
    for S, z in z_chunks(a, b):
        hi = hull_indices(z)
        keepz.append(z[hi])
        keeps.append(S[hi])
    Z = np.concatenate(keepz)
    S = np.concatenate(keeps)
    hi = hull_indices(Z)
    return Z[hi], S[hi]


def margin(a, b, U, tol_active=1e-9):
    """Full float evaluation of one configuration. Flags non-unitary U:
    ok is False when the unitarity defect exceeds UNITARY_TOL."""
    a = np.asarray(a, complex)
    n = len(a)
    d = det_d(a, b, U)
    if n <= 8:
        r = margin_points(d, z_all(a, b), n, tol_active)
    else:
        V, S = hull_lean(a, b)
        diam = diameter(V)
        dist, theta = signed_distance(d, V)
        r = dict(m=dist / diam if diam > 0 else float("nan"), dist=dist, diam=diam,
                 theta=theta, active=[tuple(int(x) for x in s) for s in S],
                 n_vertices=len(V), degenerate=len(V) < 3, support_gap=None,
                 note="n>8: active lists hull vertices only")
    r["d"] = complex(d)
    r["unitarity_defect"] = unitarity_defect(U)
    r["ok"] = r["unitarity_defect"] <= UNITARY_TOL and not r["degenerate"]
    return r


# ------------------------------------------------------- distances from trivial sets

def dist_genperm(U):
    """1 - max_sigma (1/n) sum_i |U_{i,sigma(i)}|^2; 0 iff U is a permutation
    times a diagonal unitary. Returns (value, sigma)."""
    W = np.abs(np.asarray(U)) ** 2
    r, c = linear_sum_assignment(-W)
    n = len(W)
    return float(1 - W[r, c].sum() / n), tuple(int(x) for x in c)


@lru_cache(maxsize=None)
def _subset_masks(n, k):
    return np.array([[1.0 if i in S else 0.0 for i in range(n)]
                     for S in itertools.combinations(range(n), k)])


def dist_blockdiag(U):
    """min over splits |S| = |T| = k (1 <= k < n) of sum_{i in S, j not in T}
    |U_ij|^2: 0 iff U is block diagonal after row and column permutations.
    Returns (value, (S, T))."""
    W = np.abs(np.asarray(U)) ** 2
    n = len(W)
    best = (np.inf, None)
    for k in range(1, n):
        R = _subset_masks(n, k)
        rows = R @ W.sum(axis=1)                    # sum_{i in S} rowsum_i
        L = rows[:, None] - R @ W @ R.T             # minus sum_{i in S, j in T}
        i, j = np.unravel_index(np.argmin(L), L.shape)
        if L[i, j] < best[0]:
            S = tuple(np.nonzero(R[i])[0].tolist())
            T = tuple(np.nonzero(R[j])[0].tolist())
            best = (float(L[i, j]), (S, T))
    return best


def min_separation(a, b):
    """(min_{i<j} |a_i - a_j|, min_{i<j} |b_i - b_j|) in canonical units:
    after the shift/scale normalization of moc_param.canonical_spectra
    (mean b = 0, joint rms of (a - mean a, b) equal to 1)."""
    a = np.asarray(a, complex)
    b = np.asarray(b, complex)
    t = b.mean()
    a2, b2 = a + t, b - t
    rms = math.sqrt(float((np.sum(np.abs(a2 - a2.mean()) ** 2) + np.sum(np.abs(b2) ** 2)) / (2 * len(a)))) or 1.0
    out = []
    for v in (a2, b2):
        D = np.abs(v[:, None] - v[None, :]) + np.diag(np.full(len(v), np.inf))
        out.append(float(D.min() / rms))
    return tuple(out)
