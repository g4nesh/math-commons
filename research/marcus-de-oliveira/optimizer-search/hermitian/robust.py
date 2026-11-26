"""Hull-independent signed margin, used as a guard on moc_core.margin.

signed distance(d, H) = max_theta f(theta), f(theta) = Re(e^{-i theta} d) - max_sigma Re(e^{-i theta} z_sigma),
for H = conv{z_sigma} (all n! points). Every f(theta) is a lower bound. Candidate directions:
outward normals of all hull edges (inside case) and directions from the nearest points of the
hull boundary to d (outside case). Near-duplicate z (repeated eigenvalues) cannot flip the sign
here, because the support function is taken over all n! points and not from the polygon.
"""
import itertools
import os
import sys

import mpmath as mp
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
from moc_core import z_all, hull_indices, det_d, diameter, _seg_nearest  # noqa: E402


def robust_margin(a, b, U):
    a = np.asarray(a, complex)
    b = np.asarray(b, complex)
    d = complex(det_d(a, b, U))
    z = z_all(a, b)
    V = z[hull_indices(z)]
    diam = diameter(V)
    th = []
    for k in range(len(V)):
        E = V[(k + 1) % len(V)] - V[k]
        if abs(E) > 0:
            th.append(np.angle(-1j * E))
        q = _seg_nearest(d, V[k], V[(k + 1) % len(V)])
        if abs(d - q) > 0:
            th.append(np.angle(d - q))
    th = np.array(th)
    E = np.exp(-1j * th)
    f = (E * d).real - (E[:, None] * z[None, :]).real.max(axis=1)
    k = int(np.argmax(f))
    return dict(m=float(f[k] / diam), theta=float(th[k]), diam=diam, d=d)


def robust_margin_hp(a, b, U, theta_list, dps=60):
    """Lower bound max_theta f(theta)/diam at dps digits for the given directions
    (all n! points; U rebuilt exactly unitary via moc_hp.hp_config)."""
    from moc_hp import hp_config
    am, bm, Um, Km = hp_config(a, b, U, dps)
    n = len(am)
    d = mp.det(mp.diag(am) + Um * mp.diag(bm) * Um.H)
    zs = []
    for s in itertools.permutations(range(n)):
        p = mp.mpc(1)
        for i in range(n):
            p *= am[i] + bm[s[i]]
        zs.append(p)
    best = None
    for t in theta_list:
        e = mp.expj(-mp.mpf(t))
        f = (e * d).real - max((e * zz).real for zz in zs)
        best = f if best is None or f > best else best
    # diameter from all points (max over pairs of hull candidates is enough; use float hull indices)
    zf = np.array([complex(zz) for zz in zs])
    hi = hull_indices(zf)
    diam = max(abs(zs[i] - zs[j]) for i in hi for j in hi)
    return best / diam
