"""Core numerics for the Marcus-de Oliveira conjecture (MOC).

MOC: for normal A, B with eigenvalues a, b,
    det(A + B) in conv{ z_sigma = prod_i (a_i + b_sigma(i)) : sigma in S_n }.
WLOG A = diag(a), B = U diag(b) U^*.
"""
import itertools
import numpy as np
from scipy.optimize import linprog
from scipy.spatial import ConvexHull


def haar_unitary(n, rng):
    z = (rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))) / np.sqrt(2)
    q, r = np.linalg.qr(z)
    d = np.diag(r)
    return q * (d / np.abs(d))


def det_sum(a, b, U):
    return np.linalg.det(np.diag(a) + U @ np.diag(b) @ U.conj().T)


def perm_points(a, b):
    n = len(a)
    perms = list(itertools.permutations(range(n)))
    pts = np.array([np.prod([a[i] + b[s[i]] for i in range(n)]) for s in perms])
    return perms, pts


def subsets(n, k):
    return list(itertools.combinations(range(n), k))


def cauchy_binet(a, b, U):
    """sum_{|S|=|T|} |det U[S,T]|^2 a^{S^c} b^T."""
    n = len(a)
    tot = 0j
    for k in range(n + 1):
        for S in subsets(n, k):
            Sc = [i for i in range(n) if i not in S]
            aS = np.prod([a[i] for i in Sc]) if Sc else 1.0
            for T in subsets(n, k):
                m = abs(np.linalg.det(U[np.ix_(S, T)])) ** 2 if k else 1.0
                bT = np.prod([b[j] for j in T]) if T else 1.0
                tot += m * aS * bT
    return tot


def signed_hull_distance(z, pts):
    """>0 means z is OUTSIDE conv(pts) by that Euclidean distance; <=0 inside
    (negative = -distance to boundary).  Uses 2D hull facets."""
    P = np.column_stack([pts.real, pts.imag])
    try:
        h = ConvexHull(P)
    except Exception:
        # degenerate (collinear) hull: fall back to LP-ish projection
        return _degenerate_dist(z, P)
    eq = h.equations  # rows: n_x x + n_y y + c <= 0 inside, (n unit)
    v = eq[:, 0] * z.real + eq[:, 1] * z.imag + eq[:, 2]
    m = v.max()
    if m <= 0:
        return m
    # true Euclidean distance to polygon
    verts = P[h.vertices]
    # polygon edges in order
    order = h.vertices
    poly = P[order]
    best = np.inf
    q = np.array([z.real, z.imag])
    for i in range(len(poly)):
        p0, p1 = poly[i], poly[(i + 1) % len(poly)]
        d = p1 - p0
        t = np.clip(np.dot(q - p0, d) / np.dot(d, d), 0, 1)
        best = min(best, np.linalg.norm(q - (p0 + t * d)))
    return best


def _degenerate_dist(z, P):
    # points collinear: distance from z to segment spanned
    c = P.mean(0)
    X = P - c
    u, s, vt = np.linalg.svd(X)
    d = vt[0]
    t = X @ d
    p0, p1 = c + t.min() * d, c + t.max() * d
    q = np.array([z.real, z.imag])
    dd = p1 - p0
    if np.dot(dd, dd) == 0:
        return np.linalg.norm(q - p0)
    tt = np.clip(np.dot(q - p0, dd) / np.dot(dd, dd), 0, 1)
    return np.linalg.norm(q - (p0 + tt * dd))


# ---------- compound Birkhoff data ----------

def compound_vector(U, ks=None):
    """Vector of |det U[S,T]|^2 over all S,T with |S|=|T|=k, k in ks."""
    n = U.shape[0]
    if ks is None:
        ks = range(1, n)
    out = []
    for k in ks:
        for S in subsets(n, k):
            for T in subsets(n, k):
                out.append(abs(np.linalg.det(U[np.ix_(S, T)])) ** 2)
    return np.array(out)


def perm_matrix(s):
    n = len(s)
    P = np.zeros((n, n))
    for i in range(n):
        P[i, s[i]] = 1
    return P


def compound_perm_vectors(n, ks=None):
    perms = list(itertools.permutations(range(n)))
    V = np.array([compound_vector(perm_matrix(s), ks) for s in perms])
    return perms, V


def compound_birkhoff_lp(U, V=None, ks=None):
    """Is compound_vector(U) in conv(rows of V)?  Returns (feasible, slack)
    where slack = min over lambda of max violation (L_inf), via LP."""
    n = U.shape[0]
    if V is None:
        _, V = compound_perm_vectors(n, ks)
    p = compound_vector(U, ks)
    m, d = V.shape
    # variables: lambda (m), t (1). minimize t s.t. |V^T lambda - p| <= t, sum lambda = 1, lambda>=0
    c = np.zeros(m + 1); c[-1] = 1
    A_ub = np.vstack([np.hstack([V.T, -np.ones((d, 1))]),
                      np.hstack([-V.T, -np.ones((d, 1))])])
    b_ub = np.concatenate([p, -p])
    A_eq = np.hstack([np.ones((1, m)), np.zeros((1, 1))])
    res = linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=[1.0],
                  bounds=[(0, None)] * (m + 1), method="highs")
    return res.fun, res
