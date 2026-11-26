"""High-precision (mpmath) re-evaluation of a saved configuration.

The float U is first written as U diag(R) = cayley(K) (moc_param.cayley_coords);
K is made exactly skew-Hermitian and U is rebuilt from K at `dps` digits, so it is
unitary to working precision. Then d, all n! points z_sigma, the hull and the
signed margin are recomputed at `dps` digits (n <= 8). Independent of moc_core's
hull code (own monotone chain on mpf).
"""
import itertools

import mpmath as mp
import numpy as np

from moc_param import cayley_coords


def _cayley_mp(K):
    n = K.rows
    I = mp.eye(n)
    return (I + K) ** -1 * (I - K)


def hp_config(a, b, U, dps=60):
    """Exactly-unitary (to dps digits) version of a float configuration."""
    mp.mp.dps = dps
    K, R, gap = cayley_coords(np.asarray(U, complex))
    n = len(a)
    Km = mp.matrix(n, n)
    for i in range(n):
        for j in range(n):
            Km[i, j] = mp.mpc(K[i, j].real, K[i, j].imag)
    Km = (Km - Km.H) / 2                 # exactly skew-Hermitian
    Um = _cayley_mp(Km)
    am = [mp.mpc(complex(v).real, complex(v).imag) for v in a]
    bm = [mp.mpc(complex(v).real, complex(v).imag) for v in b]
    return am, bm, Um, Km


def _hull(pts):
    P = sorted(set(pts))
    if len(P) <= 2:
        return P
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
    return lo[:-1] + up[:-1]


def _segdist(p, u, v):
    ex, ey = v[0] - u[0], v[1] - u[1]
    L2 = ex * ex + ey * ey
    t = 0 if L2 == 0 else ((p[0] - u[0]) * ex + (p[1] - u[1]) * ey) / L2
    t = min(1, max(0, t))
    qx, qy = u[0] + t * ex, u[1] + t * ey
    return mp.sqrt((p[0] - qx) ** 2 + (p[1] - qy) ** 2)


def hp_margin(a, b, U, dps=60):
    """Returns dict(m, dist, diam, d, unitarity_defect, dps) as mpf / strings."""
    am, bm, Um, Km = hp_config(a, b, U, dps)
    n = len(am)
    M = mp.diag(am) + Um * mp.diag(bm) * Um.H
    d = mp.det(M)
    udef = mp.mnorm(Um * Um.H - mp.eye(n), 1)
    pts = []
    for s in itertools.permutations(range(n)):
        z = mp.mpc(1)
        for i in range(n):
            z *= am[i] + bm[s[i]]
        pts.append((z.real, z.imag))
    V = _hull(pts)
    diam = max((mp.sqrt((p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2) for p in V for q in V), default=mp.mpf(0))
    p = (d.real, d.imag)
    if len(V) >= 3:
        s = []
        for k in range(len(V)):
            u, v = V[k], V[(k + 1) % len(V)]
            ex, ey = v[0] - u[0], v[1] - u[1]
            L = mp.sqrt(ex * ex + ey * ey)
            s.append(((p[0] - u[0]) * ey - (p[1] - u[1]) * ex) / L)   # outward normal (ey, -ex)
        if max(s) <= 0:
            dist = max(s)
        else:
            dist = min(_segdist(p, V[k], V[(k + 1) % len(V)]) for k in range(len(V)))
    elif len(V) == 2:
        dist = _segdist(p, V[0], V[1])
    else:
        dist = mp.sqrt((p[0] - V[0][0]) ** 2 + (p[1] - V[0][1]) ** 2)
    m = dist / diam if diam > 0 else mp.nan
    return dict(m=m, dist=dist, diam=diam, d=d, unitarity_defect=udef, dps=dps,
                n_vertices=len(V), m_str=mp.nstr(m, 30))
