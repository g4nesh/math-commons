"""High-precision test of an edge contact found by joint5.py.

For a saved configuration (a, b, U) whose d lies (in floating point) on a hull
edge [z_sigma, z_tau] away from its endpoints, take the exact outward normal
theta of that edge and h = Re(e^{-i theta} z_sigma). With a, b fixed (their
float values taken as exact), maximize
    F(x) = [Re(e^{-i theta} d(U0 cayley(K(x)))) - h] / diam(H)
over the 25 real coordinates x by Newton's method on grad F = 0 at DPS digits
(finite-difference gradient and Hessian, eigenvalue-truncated pseudo-inverse,
since the critical set has positive dimension). A violation of MOC along this
edge would show up as F > 0; an exact contact as F -> 0 at the working
precision; a near miss as F converging to a negative value.

Usage: [SNAP=1] python refine_hp.py RUN.json INDEX [DPS] [STEPS]\nSNAP=1 first sets b_j = -a_i exactly for every pair with |a_i + b_j| < 1e-6.
"""
import itertools
import json
import os
import sys

import mpmath as mp
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
from moc_hp import hp_config, _hull  # noqa: E402


def main():
    f, idx = sys.argv[1], int(sys.argv[2])
    dps = int(sys.argv[3]) if len(sys.argv) > 3 else 100
    steps = int(sys.argv[4]) if len(sys.argv) > 4 else 6
    rec = json.load(open(f))["best"][idx]
    a = np.array([complex(*v) for v in rec["a"]])
    b = np.array([complex(*v) for v in rec["b"]])
    U = np.array(rec["U_re"]) + 1j * np.array(rec["U_im"])
    mp.mp.dps = dps
    am, bm, U0, _ = hp_config(a, b, U, dps)
    n = len(am)
    snapped = []
    if os.environ.get("SNAP") == "1":
        # zero-factor contacts: make a_i + b_j = 0 exact for the pairs that are ~0 in floating point
        for i in range(n):
            for j in range(n):
                if abs(am[i] + bm[j]) < 1e-6:
                    snapped.append([i, j, mp.nstr(abs(am[i] + bm[j]), 5)])
                    bm[j] = -am[i]
    I = mp.eye(n)
    Da, Db = mp.diag(am), mp.diag(bm)
    P = list(itertools.permutations(range(n)))
    z = []
    for s in P:
        v = mp.mpc(1)
        for i in range(n):
            v *= am[i] + bm[s[i]]
        z.append(v)
    V = _hull([(w.real, w.imag) for w in z])
    diam = max(mp.sqrt((p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2) for p in V for q in V)

    def dU(UU):
        return mp.det(Da + UU * Db * UU.H)

    d0 = dU(U0)
    # nearest hull edge to d0
    best = None
    for k in range(len(V)):
        u, w = V[k], V[(k + 1) % len(V)]
        ex, ey = w[0] - u[0], w[1] - u[1]
        L2 = ex * ex + ey * ey
        t = ((d0.real - u[0]) * ex + (d0.imag - u[1]) * ey) / L2
        tc = min(1, max(0, t))
        dist = mp.sqrt((d0.real - u[0] - tc * ex) ** 2 + (d0.imag - u[1] - tc * ey) ** 2)
        if best is None or dist < best[0]:
            best = (dist, k, t, u, w)
    _, k, t0, u, w = best
    normal = mp.mpc(w[1] - u[1], -(w[0] - u[0]))      # outward for CCW order
    e = mp.conj(normal) / abs(normal)                   # e^{-i theta}
    h = (e * mp.mpc(u[0], u[1])).real
    ends = [P[j] for j in range(len(z)) if abs(z[j] - mp.mpc(u[0], u[1])) < mp.mpf(10) ** (-dps // 2)
            or abs(z[j] - mp.mpc(w[0], w[1])) < mp.mpf(10) ** (-dps // 2)]
    on_line = [P[j] for j in range(len(z)) if abs((e * z[j]).real - h) < mp.mpf(10) ** (-dps // 2) * diam]

    def U_of(x):
        K = mp.matrix(n, n)
        for i in range(n):
            for j in range(n):
                X, XT = x[i * n + j], x[j * n + i]
                K[i, j] = mp.mpc((X - XT) / 2, (X + XT) / 2)
        return U0 * (mp.inverse(I + K) * (I - K))

    def F(x):
        return ((e * dU(U_of(x))).real - h) / diam

    m = n * n
    x = [mp.mpf(0)] * m
    hstep = mp.mpf(10) ** (-(dps // 4))
    hist = []
    for it in range(steps + 1):
        f0 = F(x)
        hist.append(mp.nstr(f0, 12))
        if it == steps:
            break
        g = mp.matrix(m, 1)
        H = mp.matrix(m, m)
        fp, fm = [None] * m, [None] * m
        for i in range(m):
            xp = list(x); xp[i] += hstep
            xm = list(x); xm[i] -= hstep
            fp[i], fm[i] = F(xp), F(xm)
            g[i] = (fp[i] - fm[i]) / (2 * hstep)
            H[i, i] = (fp[i] - 2 * f0 + fm[i]) / hstep ** 2
        for i in range(m):
            for j in range(i + 1, m):
                xpp = list(x); xpp[i] += hstep; xpp[j] += hstep
                xmm = list(x); xmm[i] -= hstep; xmm[j] -= hstep
                v = (F(xpp) - fp[i] - fp[j] + 2 * f0 - fm[i] - fm[j] + F(xmm)) / (2 * hstep ** 2)
                H[i, j] = H[j, i] = v
        E, Q = mp.eigsy(H)
        lam_max = max(abs(E[i]) for i in range(m))
        step = mp.matrix(m, 1)
        nkeep = 0
        for i in range(m):
            if abs(E[i]) > mp.mpf(10) ** (-12) * lam_max:
                nkeep += 1
                c = sum(Q[r, i] * g[r] for r in range(m)) / E[i]
                for r in range(m):
                    step[r] -= c * Q[r, i]
        x = [x[i] + step[i] for i in range(m)]
        print("it %d F %s |g| %s rank %d max|step| %s" % (it, mp.nstr(f0, 12), mp.nstr(mp.norm(g), 5), nkeep,
                                                       mp.nstr(max(abs(s) for s in step), 5)), flush=True)
    Uf = U_of(x)
    df = dU(Uf)
    tf = ((df.real - u[0]) * (w[0] - u[0]) + (df.imag - u[1]) * (w[1] - u[1])) / ((w[0] - u[0]) ** 2 + (w[1] - u[1]) ** 2)
    W2 = np.array([[float(abs(Uf[i, j]) ** 2) for j in range(n)] for i in range(n)])
    sys.path.insert(0, os.path.join(HERE, "..", "lib"))
    from moc_core import dist_genperm, dist_blockdiag
    out = dict(snapped_pairs=snapped, source=f, index=idx, restart=rec["restart"], seed=rec["seed"], dps=dps, F_history=hist,
               F_final=mp.nstr(F(x), 20), edge_position_t=mp.nstr(tf, 15), edge_endpoints=[list(s) for s in ends],
               perms_on_support_line=[list(s) for s in on_line],
               unitarity_defect=mp.nstr(mp.mnorm(Uf * Uf.H - I, 1), 5),
               dist_genperm=dist_genperm(W2 ** 0.5)[0], dist_blockdiag=dist_blockdiag(W2 ** 0.5)[0],
               absU2=W2.tolist(),
               a=[[mp.nstr(v.real, dps), mp.nstr(v.imag, dps)] for v in am],
               b=[[mp.nstr(v.real, dps), mp.nstr(v.imag, dps)] for v in bm],
               U=[[[mp.nstr(Uf[i, j].real, dps), mp.nstr(Uf[i, j].imag, dps)] for j in range(n)] for i in range(n)])
    tag = os.path.splitext(os.path.basename(f))[0]
    path = os.path.join(HERE, "contacts", "edge_%s_%d%s.json" % (tag, idx, "_snap" if snapped else ""))
    os.makedirs(os.path.dirname(path), exist_ok=True)
    json.dump(out, open(path, "w"), indent=1)
    print(json.dumps({k: out[k] for k in ("F_history", "F_final", "edge_position_t", "edge_endpoints",
                                          "perms_on_support_line", "unitarity_defect", "dist_genperm", "dist_blockdiag")}))


if __name__ == "__main__":
    main()
