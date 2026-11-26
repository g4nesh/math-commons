"""Second variation of the MOC support inequalities at a lifted configuration.

Setting. (a, b) in C^n, U0 = U4 (+) 1 (block form, last index decoupled), and d(U0) on the
boundary of H = conv{z_sigma}. Off-block perturbation U = U0 expm(eps K_x) with K_x[i, n-1] = x_i,
K_x[n-1, i] = -conj(x_i), x in C^(n-1).

Facts used (checked numerically in check_lift.py):
  * d is even in x (conjugation by J = diag(1,..,1,-1) commutes with D_a, D_b and U0), so the
    first variation of d in x is 0 and the mixed second derivative with in-block / spectral
    directions is 0.
  * d(eps) = d0 + eps^2 w(x) + O(eps^4) with w(x) = x^T W x a complex quadratic form in the
    8 real coordinates (Re x, Im x).
  * The z_sigma do not depend on U. So the lifted contact is a strict local violation at second
    order iff w(x) leaves the tangent cone of H at d0 for some x, i.e. iff
    lambda_max(Re(e^{-i theta} W)) > 0 for an extreme outward normal theta of H at d0
    (the unique edge normal at an edge contact, the two adjacent edge normals at a vertex).
"""
import sys, os
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import numpy as np
import mpmath as mp
from scipy.linalg import expm
from moc_core import z_all, det_d, margin, perms


def hull_indices(z, merge=1e-12, coll=1e-12):
    """Robust CCW hull (replaces moc_core.hull_indices, which returned a non-convex polygon when
    many points nearly coincide, e.g. when a_i + b_j = 0 makes (n-1)! of the z vanish).
    Near-duplicates (within merge*scale) are merged first; near-collinear middle points
    (cross < coll * |p-o| * |q-o|) are dropped."""
    z = np.asarray(z, complex)
    scale = max(1e-300, float(np.max(np.abs(z - z.mean()))))
    order = np.lexsort((z.imag, z.real))
    keep = []
    for i in order:
        if not any(abs(z[i] - z[j]) <= merge * scale for j in keep[-64:]):
            keep.append(int(i))
    pts = z[keep]

    def cross(o, p, q):
        u, v = pts[p] - pts[o], pts[q] - pts[o]
        c = u.real * v.imag - u.imag * v.real
        return c if abs(c) > coll * abs(u) * abs(v) else 0.0
    lower, upper = [], []
    for k in range(len(pts)):
        while len(lower) >= 2 and cross(lower[-2], lower[-1], k) <= 0:
            lower.pop()
        lower.append(k)
    for k in range(len(pts) - 1, -1, -1):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], k) <= 0:
            upper.pop()
        upper.append(k)
    chain = lower[:-1] + upper[:-1]
    return [keep[k] for k in chain]


def check_convex(z, idx):
    V = np.asarray(z, complex)[idx]
    k = len(V)
    for i in range(k):
        e = V[(i + 1) % k] - V[i]
        nrm = -1j * e / abs(e)
        if np.max(((np.asarray(z) - V[i]) * np.conj(nrm)).real) > 1e-9 * max(1.0, np.max(np.abs(z))):
            return False
    return True


def robust_margin(d, z):
    """Sign-correct margin: max over hull edges of the outward signed distance of d to the
    edge line, divided by diam. Equals the true m for points inside, and is positive iff outside."""
    z = np.asarray(z, complex)
    idx = hull_indices(z)
    V = z[idx]
    k = len(V)
    diam = max(abs(V[i] - V[j]) for i in range(k) for j in range(k))
    best = -np.inf
    for i in range(k):
        e = V[(i + 1) % k] - V[i]
        nrm = -1j * e / abs(e)
        best = max(best, ((d - V[i]) * np.conj(nrm)).real)
    return best / diam


def K_off(x, n):
    K = np.zeros((n, n), complex)
    K[:n - 1, n - 1] = x
    K[n - 1, :n - 1] = -np.conj(x)
    return K


def w_quad(a, b, U0, x):
    """Second-order coefficient w(x) = (1/2) d''(0) of d(U0 expm(eps K_x)), analytic."""
    n = len(a)
    K = K_off(np.asarray(x, complex), n)
    Db = np.diag(b)
    C1 = K @ Db - Db @ K
    C2 = K @ C1 - C1 @ K
    M0 = np.diag(a) + U0 @ Db @ U0.conj().T
    N = U0 @ C1 @ U0.conj().T
    S = 0.5 * U0 @ C2 @ U0.conj().T
    X = np.linalg.solve(M0, N)
    Y = np.linalg.solve(M0, S)
    return np.linalg.det(M0) * (np.trace(Y) + 0.5 * np.trace(X) ** 2 - 0.5 * np.trace(X @ X))


def W_matrix(a, b, U0):
    """Complex symmetric 2(n-1) x 2(n-1) matrix W with w(x) = v^T W v, v = (Re x, Im x)."""
    n = len(a)
    m = 2 * (n - 1)

    def vec2x(v):
        return v[:n - 1] + 1j * v[n - 1:]
    E = np.eye(m)
    q = [w_quad(a, b, U0, vec2x(E[i])) for i in range(m)]
    W = np.zeros((m, m), complex)
    for i in range(m):
        W[i, i] = q[i]
        for j in range(i + 1, m):
            W[i, j] = W[j, i] = 0.5 * (w_quad(a, b, U0, vec2x(E[i] + E[j])) - q[i] - q[j])
    return W


def w_fd(a, b, U0, x, h=1e-4):
    n = len(a)
    K = K_off(np.asarray(x, complex), n)
    dp = det_d(a, b, U0 @ expm(h * K))
    dm = det_d(a, b, U0 @ expm(-h * K))
    d0 = det_d(a, b, U0)
    return (dp + dm - 2 * d0) / (2 * h * h)


def normals_at(d, z, tol=1e-9):
    """Extreme outward normal angles of H = conv z at the boundary point d.
    Returns (kind, [theta...], info). kind: 'vertex', 'edge' or 'interior'/'outside'."""
    z = np.asarray(z, complex)
    idx = hull_indices(z)
    V = z[idx]
    k = len(V)
    diam = max(abs(V[i] - V[j]) for i in range(k) for j in range(k))
    # vertex?
    dv = np.abs(V - d)
    iv = int(np.argmin(dv))
    if dv[iv] < tol * diam:
        e_in = V[iv] - V[iv - 1]
        e_out = V[(iv + 1) % k] - V[iv]
        th = [float(np.angle(-1j * e_in)), float(np.angle(-1j * e_out))]
        return "vertex", th, dict(vertex_perm=int(idx[iv]), diam=diam, dist=float(dv[iv]))
    # edge? pick the hull edge nearest to d as a segment
    best = None
    for i in range(k):
        u, v = V[i], V[(i + 1) % k]
        e = v - u
        nrm = -1j * e / abs(e)
        s = ((d - u) * np.conj(nrm)).real   # signed distance outward from the edge line
        t = ((d - u) * np.conj(e)).real / abs(e) ** 2
        tc = min(1.0, max(0.0, t))
        dseg = abs(d - (u + tc * e))
        if best is None or dseg < best[4]:
            best = (s, t, float(np.angle(nrm)), i, dseg)
    s, t, th, i, dseg = best
    if abs(s) < tol * diam and -1e-12 < t < 1 + 1e-12:
        return "edge", [th], dict(edge=(int(idx[i]), int(idx[(i + 1) % k])), t=float(t), diam=diam, dist=float(s))
    return "not-boundary", [], dict(diam=diam, dist=float(s))


def lift_curvature(a, b, U0, tol=1e-9):
    """For a boundary configuration: kind, normals, and lambda_max(Re(e^{-i th} W))/diam
    for each extreme normal th, plus the eigenvector."""
    z = z_all(a, b)
    d = det_d(a, b, U0)
    kind, ths, info = normals_at(d, z, tol)
    W = W_matrix(a, b, U0)
    out = []
    for th in ths:
        Hr = (np.exp(-1j * th) * W).real
        Hr = 0.5 * (Hr + Hr.T)
        ev, evec = np.linalg.eigh(Hr)
        out.append(dict(theta=th, lam_max=float(ev[-1] / info["diam"]),
                        lam_min=float(ev[0] / info["diam"]), eigs=(ev / info["diam"]).tolist(),
                        vmax=evec[:, -1].tolist()))
    return kind, info, out, W


# ----------------------------------------------------------------- high precision

def W_matrix_mp(a, b, U0, dps=40):
    """Same W in mpmath; U0 must be exactly unitary at this precision (pass mp matrices
    built e.g. from a Cayley transform)."""
    mp.mp.dps = dps
    n = len(a)
    m = 2 * (n - 1)
    a = [mp.mpc(x) for x in a]
    b = [mp.mpc(x) for x in b]
    Db = mp.diag(b)
    Uh = U0.transpose_conj()
    M0 = mp.diag(a) + U0 * Db * Uh

    def wq(x):
        K = mp.zeros(n, n)
        for i in range(n - 1):
            K[i, n - 1] = x[i]
            K[n - 1, i] = -mp.conj(x[i])
        C1 = K * Db - Db * K
        C2 = K * C1 - C1 * K
        N = U0 * C1 * Uh
        S = U0 * C2 * Uh * mp.mpf("0.5")
        Mi = M0 ** -1
        X = Mi * N
        Y = Mi * S
        tr = lambda A: sum(A[i, i] for i in range(n))
        return mp.det(M0) * (tr(Y) + tr(X) ** 2 / 2 - tr(X * X) / 2)

    def e(i):
        v = [mp.mpc(0)] * (n - 1)
        if i < n - 1:
            v[i] = mp.mpc(1)
        else:
            v[i - (n - 1)] = mp.mpc(0, 1)
        return v
    q = [wq(e(i)) for i in range(m)]
    W = mp.zeros(m, m)
    for i in range(m):
        W[i, i] = q[i]
        for j in range(i + 1, m):
            xi, xj = e(i), e(j)
            W[i, j] = W[j, i] = (wq([xi[k] + xj[k] for k in range(n - 1)]) - q[i] - q[j]) / 2
    return W
