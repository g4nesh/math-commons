"""Fixed-U search: hold a structured unitary U fixed and maximize the signed
margin m over the spectra (a, b) and the direction theta.

Smooth objective (same form as lib/moc_opt.py, with U frozen):
    F = [Re(e^{-i th} d) - LSE_tau(Re(e^{-i th} z_sigma))] / s,
where s is the rms spread of the z_sigma. For n <= 6 every permutation enters
the log-sum-exp; for n >= 7 a working set (the K permutations with the largest
projection, re-selected from the full list at each annealing stage) enters it,
and s is taken over a fixed random subsample. Every reported m comes from
moc_core.margin, which uses all n! points (hull_lean for n = 9).

Optional scale-invariant penalties keep the spectra away from the trivial
zero-margin sets: sep (minimum gap within a and within b, relative to the joint
rms) and fat (minimum sqrt(lambda2/lambda1) of the 2x2 covariance of a and of b,
which excludes collinear spectra, where the hull degenerates to a segment).
"""
import math
import os
import sys

import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib"))
from moc_core import perms, margin, z_chunks, min_separation  # noqa: E402
from moc_param import canonical_spectra  # noqa: E402


def adjugate(M):
    W, s, Vh = np.linalg.svd(M)
    n = len(s)
    ps = np.array([np.prod(np.delete(s, i)) for i in range(n)])
    dW = np.linalg.det(W)
    dV = np.linalg.det(Vh)
    return dV * dW * (Vh.conj().T @ np.diag(ps) @ W.conj().T), dW * dV * np.prod(s)


def all_perms(n):
    if n <= 8:
        return perms(n).astype(np.intp)
    out = []
    Q = perms(n - 1).astype(np.intp)
    for j in range(n):
        rest = np.array([k for k in range(n) if k != j], dtype=np.intp)
        S = np.empty((len(Q), n), dtype=np.intp)
        S[:, 0] = j
        S[:, 1:] = rest[Q]
        out.append(S)
    return np.concatenate(out)


def projections(a, b, th, n):
    """Re(e^{-i th} z_sigma) for all n! sigma, in all_perms order."""
    e = np.exp(-1j * th)
    if n <= 8:
        z = (a[None, :] + b[perms(n)]).prod(1)
        return (e * z).real
    return np.concatenate([(e * z).real for _, z in z_chunks(a, b)])


def zcloud_ratio(a, b, Ps):
    """sqrt(lambda_min / lambda_max) of the 2x2 covariance of the points z_sigma, sigma in Ps."""
    z = (a[None, :] + b[Ps]).prod(1)
    c = z - z.mean()
    X = np.stack([c.real, c.imag])
    ev = np.linalg.eigvalsh(X @ X.T)
    return math.sqrt(max(ev[0], 0) / ev[1]) if ev[1] > 0 else 0.0


def hull_width_ratio(V):
    """Minimum width / diameter of the convex polygon with CCW vertices V."""
    V = np.asarray(V, complex)
    if len(V) < 3:
        return 0.0
    diam = float(np.abs(V[:, None] - V[None, :]).max())
    best = np.inf
    for k in range(len(V)):
        e = V[(k + 1) % len(V)] - V[k]
        nrm = -1j * e / abs(e)
        best = min(best, float(np.max(np.abs((np.conj(nrm) * (V - V[k])).real))))
    return best / diam


def penalty(a, b, sep, fat, zfat=0.0, Ps=None):
    n = len(a)
    ca = a - a.mean()
    cb = b - b.mean()
    rms = math.sqrt((np.sum(np.abs(ca) ** 2) + np.sum(np.abs(cb) ** 2)) / (2 * n))
    if rms == 0:
        return 1e6
    p = 0.0
    iu = np.triu_indices(n, 1)
    for v in (a, b):
        if sep > 0:
            g = np.abs(v[:, None] - v[None, :])[iu] / rms
            h = np.clip(sep - g, 0, None)
            p += float(np.sum(h * h))
        if fat > 0:
            c = v - v.mean()
            X = np.stack([c.real, c.imag])
            C = X @ X.T
            ev = np.linalg.eigvalsh(C)
            f = math.sqrt(max(ev[0], 0) / ev[1]) if ev[1] > 0 else 0.0
            if f < fat:
                p += (fat - f) ** 2
    if zfat > 0:
        # offset guard: |mean a + mean b| <= 5 rms (a large common offset is the
        # constant-spectrum limit and destroys float accuracy, which grows like offset^n)
        off = abs(a.mean() + b.mean()) / rms
        if off > 5:
            p += (off - 5) ** 2
        f = zcloud_ratio(a, b, Ps)
        if f < zfat:
            p += (zfat - f) ** 2
    return p


class SpecProblem:
    def __init__(self, U, P_lse, P_s, tau, sep=0.0, fat=0.0, zfat=0.0, mu=100.0):
        self.U = np.asarray(U, complex)
        self.n = len(U)
        self.P = P_lse
        self.Ps = P_s
        self.tau, self.sep, self.fat, self.zfat, self.mu = tau, sep, fat, zfat, mu
        self.nev = 0

    def _zk(self, a, b, P):
        n = self.n
        Fm = a[None, :] + b[P]
        N = len(P)
        pre = np.ones((N, n), complex)
        suf = np.ones((N, n), complex)
        pre[:, 1:] = np.cumprod(Fm[:, :-1], axis=1)
        suf[:, :-1] = np.cumprod(Fm[:, ::-1][:, :-1], axis=1)[:, ::-1]
        return pre * suf, pre[:, -1] * Fm[:, -1]

    def fg(self, x):
        self.nev += 1
        n, U = self.n, self.U
        a = x[0:n] + 1j * x[n:2 * n]
        b = x[2 * n:3 * n] + 1j * x[3 * n:4 * n]
        th = x[-1]
        M = np.diag(a) + U @ np.diag(b) @ U.conj().T
        adj, d = adjugate(M)
        e = np.exp(-1j * th)
        Zk, z = self._zk(a, b, self.P)
        if self.Ps is self.P:
            Zks, zs = Zk, z
        else:
            Zks, zs = self._zk(a, b, self.Ps)
        Ns = len(zs)
        zbar = zs.mean()
        s2 = float(np.mean(np.abs(zs - zbar) ** 2))
        s = math.sqrt(s2) if s2 > 0 else 1e-300
        proj = (e * z).real
        T = self.tau * s
        mx = proj.max()
        w = np.exp((proj - mx) / T)
        Zs = w.sum()
        p = w / Zs
        lse = mx + T * math.log(Zs)
        num = (e * d).real - lse
        F = num / s
        P = self.P
        ga_d = np.diag(adj)
        gb_d = np.diag(U.conj().T @ adj @ U)
        ga_z = p @ Zk
        pw = p[:, None] * Zk
        gb_z = np.bincount(P.ravel(), weights=pw.real.ravel(), minlength=n) \
            + 1j * np.bincount(P.ravel(), weights=pw.imag.ravel(), minlength=n)
        dlse_dT = math.log(Zs) - float(p @ (proj - mx)) / T
        cz = np.conj(zs - zbar)
        Ha = cz @ Zks
        cw = cz[:, None] * Zks
        Hb = np.bincount(self.Ps.ravel(), weights=cw.real.ravel(), minlength=n) \
            + 1j * np.bincount(self.Ps.ravel(), weights=cw.imag.ravel(), minlength=n)
        coef = ((-self.tau * dlse_dT) / s - num / s ** 2) / (2 * s)
        g = np.zeros_like(x)
        ga = e * (ga_d - ga_z) / s + coef * (2.0 / Ns) * Ha
        gb = e * (gb_d - gb_z) / s + coef * (2.0 / Ns) * Hb
        g[0:n], g[n:2 * n] = ga.real, -ga.imag
        g[2 * n:3 * n], g[3 * n:4 * n] = gb.real, -gb.imag
        g[-1] = (float((e * d).imag) - float(p @ (e * z).imag)) / s
        obj, grad = -F, -g
        if self.sep > 0 or self.fat > 0 or self.zfat > 0:
            pk = (self.sep, self.fat, self.zfat, self.Ps)
            pv = penalty(a, b, *pk)
            if pv > 0:
                obj += self.mu * pv
                h = 1e-7
                for k in range(4 * n):
                    xp = x.copy()
                    xm = x.copy()
                    xp[k] += h
                    xm[k] -= h
                    ap = xp[0:n] + 1j * xp[n:2 * n]
                    bp = xp[2 * n:3 * n] + 1j * xp[3 * n:4 * n]
                    am = xm[0:n] + 1j * xm[n:2 * n]
                    bm = xm[2 * n:3 * n] + 1j * xm[3 * n:4 * n]
                    grad[k] += self.mu * (penalty(ap, bp, *pk) - penalty(am, bm, *pk)) / (2 * h)
        return obj, grad


def pack(a, b, th):
    return np.concatenate([a.real, a.imag, b.real, b.imag, [th]])


def unpack(x, n):
    return x[0:n] + 1j * x[n:2 * n], x[2 * n:3 * n] + 1j * x[3 * n:4 * n], x[-1]


TAUS = (0.1, 0.03, 0.01, 3e-3, 1e-3, 3e-4, 1e-4)


def maximize_spectra(U, a, b, rng, taus=TAUS, maxiter=250, sep=0.0, fat=0.0, zfat=0.0, K=4000, Ksub=5040):
    """Annealed L-BFGS over (a, b, theta) with U fixed. Returns (a, b, theta, nev)."""
    n = len(U)
    full = n <= 6
    Pall = all_perms(n) if not full else perms(n).astype(np.intp)
    Ps = Pall if full else Pall[rng.choice(len(Pall), size=min(Ksub, len(Pall)), replace=False)]
    nev = 0
    for tau in taus:
        a, b = canonical_spectra(a, b)
        th = margin(a, b, U)["theta"] if n <= 8 else _theta9(a, b, U)
        if full:
            Pl = Pall
            Psub = Pall
        else:
            pr = projections(a, b, th, n)
            idx = np.argpartition(-pr, K)[:K]
            Pl = Pall[idx]
            Psub = Ps
        prob = SpecProblem(U, Pl, Psub, tau, sep=sep, fat=fat, zfat=zfat)
        res = minimize(prob.fg, pack(a, b, th), jac=True, method="L-BFGS-B",
                       options=dict(maxiter=maxiter, gtol=1e-12, ftol=1e-15))
        a, b, th = unpack(res.x, n)
        nev += prob.nev
    a, b = canonical_spectra(a, b)
    return a, b, th, nev


def _theta9(a, b, U):
    return margin(a, b, U)["theta"]


def polish_fixed(U, a, b, th=None, rounds=6, K=3000, band=0.3, maxiter=500, sep=0.0, fat=0.0, zfat=0.0,
                 rng=None):
    """Active-set SLSQP with U fixed: maximize t subject to
    Re(e^{-i th}(d - z_sigma)) / s0 >= t on a working set, joint rms of the
    centred spectra = 1 and mean(b) = 0. Optional guard inequalities
    (relative gaps >= sep, fatness >= fat). Returns (a, b, th, info); the
    caller re-evaluates with moc_core.margin."""
    U = np.asarray(U, complex)
    n = len(U)
    a, b = canonical_spectra(a, b)
    r0 = margin(a, b, U)
    th = r0["theta"] if th is None else th
    Pall = all_perms(n)
    pr = projections(a, b, th, n)
    s0 = float(r0["diam"])
    order = np.argsort(-pr)
    work = set(order[:min(K, len(order))][pr[order[:min(K, len(order))]] >= pr.max() - band * s0].tolist())
    work |= set(order[:50].tolist())
    best = (r0["m"], a, b, th)
    info = dict(rounds=0)

    def eqf(y):
        aa, bb, _ = unpack(y[:4 * n + 1], n)
        rms2 = (np.sum(np.abs(aa - aa.mean()) ** 2) + np.sum(np.abs(bb) ** 2 - 0 * bb.real)) / (2 * n)
        rms2 = (np.sum(np.abs(aa - aa.mean()) ** 2) + np.sum(np.abs(bb - bb.mean()) ** 2)) / (2 * n)
        mb = bb.mean()
        return np.array([rms2 - 1.0, mb.real, mb.imag])

    def eqj(y):
        aa, bb, _ = unpack(y[:4 * n + 1], n)
        J = np.zeros((3, 4 * n + 2))
        ca, cb = aa - aa.mean(), bb - bb.mean()
        J[0, 0:n], J[0, n:2 * n] = ca.real / n, ca.imag / n
        J[0, 2 * n:3 * n], J[0, 3 * n:4 * n] = cb.real / n, cb.imag / n
        J[1, 2 * n:3 * n] = 1.0 / n
        J[2, 3 * n:4 * n] = 1.0 / n
        return J

    for rnd in range(rounds):
        idx = np.array(sorted(work))
        Pw = Pall[idx]
        N = len(Pw)
        rows = np.arange(N)[:, None]

        def parts(y):
            aa, bb, tt = unpack(y[:4 * n + 1], n)
            M = np.diag(aa) + U @ np.diag(bb) @ U.conj().T
            adj, d = adjugate(M)
            return aa, bb, tt, adj, d

        def cons(y):
            aa, bb, tt, adj, d = parts(y)
            zz = (aa[None, :] + bb[Pw]).prod(1)
            return (np.exp(-1j * tt) * (d - zz)).real / s0 - y[-1]

        def cons_jac(y):
            aa, bb, tt, adj, d = parts(y)
            e = np.exp(-1j * tt)
            Fm = aa[None, :] + bb[Pw]
            pre = np.ones((N, n), complex)
            suf = np.ones((N, n), complex)
            pre[:, 1:] = np.cumprod(Fm[:, :-1], axis=1)
            suf[:, :-1] = np.cumprod(Fm[:, ::-1][:, :-1], axis=1)[:, ::-1]
            Zk = pre * suf
            zz = pre[:, -1] * Fm[:, -1]
            Gb = np.zeros((N, n), complex)
            Gb[rows, Pw] = Zk
            ga = e * (np.diag(adj)[None, :] - Zk)
            gb = e * (np.diag(U.conj().T @ adj @ U)[None, :] - Gb)
            J = np.zeros((N, 4 * n + 2))
            J[:, 0:n], J[:, n:2 * n] = ga.real, -ga.imag
            J[:, 2 * n:3 * n], J[:, 3 * n:4 * n] = gb.real, -gb.imag
            J[:, 4 * n] = (e * (d - zz)).imag
            J[:, :4 * n + 1] /= s0
            J[:, -1] = -1.0
            return J

        cl = [dict(type="ineq", fun=cons, jac=cons_jac), dict(type="eq", fun=eqf, jac=eqj)]
        if sep > 0 or fat > 0 or zfat > 0:
            Ps = Pall if n <= 6 else Pall[(rng or np.random.default_rng(0)).choice(len(Pall), 5040, replace=False)]
            cl.append(dict(type="ineq", fun=lambda y: np.array([
                -penalty(*unpack(y[:4 * n + 1], n)[:2], sep, fat, zfat, Ps)]) + 1e-12))
        y0 = np.concatenate([pack(a, b, th), [float(np.min(cons(np.concatenate([pack(a, b, th), [0.0]]))))]])
        res = minimize(lambda y: -y[-1], y0, jac=lambda y: -np.eye(4 * n + 2)[-1],
                       method="SLSQP", constraints=cl, options=dict(maxiter=maxiter, ftol=1e-16))
        aa, bb, tt = unpack(res.x[:4 * n + 1], n)
        r = margin(aa, bb, U)
        info["rounds"] = rnd + 1
        if r["m"] > best[0]:
            best = (r["m"], aa, bb, tt)
        pr = projections(aa, bb, tt, n)
        M = np.diag(aa) + U @ np.diag(bb) @ U.conj().T
        dproj = (np.exp(-1j * tt) * np.linalg.det(M)).real
        viol = (dproj - pr) / s0 - res.x[-1]
        new = set(np.nonzero(viol < -1e-12)[0].tolist()) - work
        a, b, th = aa, bb, tt
        if not new:
            break
        work |= set(list(new)[:2000])
    info["working"] = len(work)
    m, a, b, th = best
    return a, b, th, info
