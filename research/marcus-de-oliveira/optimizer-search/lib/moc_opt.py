"""Optimizers that MAXIMIZE the signed margin m of d(U) against H = conv{z_sigma}.

Variables x = [Re a, Im a, Re b, Im b, X (n*n real), theta], U = U0 cayley(K(X)),
K(X) = (X - X^T)/2 + i(X + X^T)/2. U0 is a fixed base point that is moved to the
current U between annealing stages (so K stays small and Cayley stays well
conditioned). The direction theta is itself a variable, so max over theta is
exact; only the max over permutations is smoothed (log-sum-exp, temperature
tau * s, where s = rms spread of the z_sigma).

smooth objective  F = [Re(e^{-i th} d) - LSE_tau(Re(e^{-i th} z_sigma))] / s  <=  true signed distance / s.
"""
import json
import math
import time

import numpy as np
from scipy.optimize import minimize

from moc_core import (perms, margin, dist_genperm, dist_blockdiag, min_separation)
from moc_param import (skew_from_real, cayley, random_unitary, random_spectra,
                       canonical_spectra)


def pack(a, b, X, theta):
    return np.concatenate([np.real(a), np.imag(a), np.real(b), np.imag(b),
                           np.asarray(X, float).ravel(), [theta]])


def unpack(x, n):
    a = x[0:n] + 1j * x[n:2 * n]
    b = x[2 * n:3 * n] + 1j * x[3 * n:4 * n]
    X = x[4 * n:4 * n + n * n].reshape(n, n)
    return a, b, X, x[-1]


def _adjugate(M):
    W, s, Vh = np.linalg.svd(M)
    n = len(s)
    ps = np.array([np.prod(np.delete(s, i)) for i in range(n)])
    dW = np.linalg.det(W)
    dV = np.linalg.det(Vh)
    return dV * dW * (Vh.conj().T @ np.diag(ps) @ W.conj().T), dW * dV * np.prod(s)


def _gradX(C):
    """d Re tr(C dK) / dX for K = skew_from_real(X)."""
    return np.real(0.5 * (1 + 1j) * C.T + 0.5 * (-1 + 1j) * C)


def _hol_to_real(g):
    """Gradient of Re(h) w.r.t. (Re v, Im v) for holomorphic h with dh/dv = g."""
    return np.real(g), -np.imag(g)


class Problem:
    """Smooth objective with analytic gradient. Options (all optional):
    sep:   minimum separation |v_i - v_j| within a and within b, in canonical
           units (canonical_spectra: joint rms 1) (penalty if below)
    dperm: minimum dist_genperm(U)                 (penalty if below)
    dblock: minimum dist_blockdiag(U)               (penalty if below)
    mu:    penalty weight."""

    def __init__(self, n, U0=None, tau=0.03, sep=0.0, dperm=0.0, dblock=0.0, mu=100.0):
        self.n = n
        self.P = perms(n).astype(np.intp)
        self.U0 = np.eye(n, dtype=complex) if U0 is None else np.asarray(U0, complex)
        self.tau, self.sep, self.dperm, self.dblock, self.mu = tau, sep, dperm, dblock, mu
        self.nev = 0

    def config(self, x):
        a, b, X, th = unpack(x, self.n)
        return a, b, self.U0 @ cayley(skew_from_real(X)), th

    def fg(self, x):
        """Returns (-F + penalties, gradient) for minimization."""
        self.nev += 1
        n, P = self.n, self.P
        a, b, X, th = unpack(x, n)
        I = np.eye(n)
        K = skew_from_real(X)
        W = np.linalg.inv(I + K)
        V = (I - K) @ W
        U = self.U0 @ V
        Db = np.diag(b)
        M = np.diag(a) + U @ Db @ U.conj().T
        adj, d = _adjugate(M)
        e = np.exp(-1j * th)
        # --- z and its partial products
        Fm = a[None, :] + b[P]
        N = len(P)
        pre = np.ones((N, n), complex)
        suf = np.ones((N, n), complex)
        pre[:, 1:] = np.cumprod(Fm[:, :-1], axis=1)
        suf[:, :-1] = np.cumprod(Fm[:, ::-1][:, :-1], axis=1)[:, ::-1]
        Zk = pre * suf
        z = pre[:, -1] * Fm[:, -1]
        zbar = z.mean()
        s2 = float(np.mean(np.abs(z - zbar) ** 2))
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
        # --- gradients of num (holomorphic parts)
        ga_d = np.diag(adj).copy()                       # dd/da_k
        gb_d = np.diag(U.conj().T @ adj @ U).copy()      # dd/db_k
        ga_z = p @ Zk                                    # sum p dz/da_k
        gb_z = np.bincount(P.ravel(), weights=(p[:, None] * Zk).real.ravel(), minlength=n) \
            + 1j * np.bincount(P.ravel(), weights=(p[:, None] * Zk).imag.ravel(), minlength=n)
        # lse also depends on s through T: d lse / dT = log Zs - sum p (proj-mx)/T
        dlse_dT = math.log(Zs) - float(p @ (proj - mx)) / T
        # s gradient (holomorphic): ds2 = (2/N) Re(H dv)
        cz = np.conj(z - zbar)
        Ha = cz @ Zk
        Hb = np.bincount(P.ravel(), weights=(cz[:, None] * Zk).real.ravel(), minlength=n) \
            + 1j * np.bincount(P.ravel(), weights=(cz[:, None] * Zk).imag.ravel(), minlength=n)
        g = np.zeros_like(x)
        ra_num, ia_num = _hol_to_real(e * (ga_d - ga_z))
        rb_num, ib_num = _hol_to_real(e * (gb_d - gb_z))
        ra_s2, ia_s2 = _hol_to_real((2.0 / N) * Ha)
        rb_s2, ib_s2 = _hol_to_real((2.0 / N) * Hb)
        ds = 1.0 / (2 * s)
        # num depends on s via T = tau s:  dnum/ds = -tau * dlse_dT
        coef = (-self.tau * dlse_dT) / s - num / s ** 2
        g[0:n] = ra_num / s + coef * ds * ra_s2
        g[n:2 * n] = ia_num / s + coef * ds * ia_s2
        g[2 * n:3 * n] = rb_num / s + coef * ds * rb_s2
        g[3 * n:4 * n] = ib_num / s + coef * ds * ib_s2
        # U part: dd = tr(C dK) with adj' = U0^H adj U0 and V
        adjV = self.U0.conj().T @ adj @ self.U0
        Xm = Db @ V.conj().T @ adjV
        Ym = adjV @ V @ Db
        C = -W @ Xm @ (I + V) + (I + V.conj().T) @ Ym @ W.conj().T
        g[4 * n:4 * n + n * n] = (_gradX(e * C) / s).ravel()
        g[-1] = (float((e * d).imag) - float(p @ (e * z).imag)) / s
        obj = -F
        grad = -g
        # --- penalties
        if self.sep > 0:
            # separation in canonical units (stages start from canonical_spectra, joint rms 1)
            for v, off in ((a, 0), (b, 2 * n)):
                for i in range(n):
                    for j in range(i + 1, n):
                        r = abs(v[i] - v[j])
                        if r < self.sep:
                            h = self.sep - r
                            obj += self.mu * h * h
                            u = (v[i] - v[j]) / max(r, 1e-300)
                            gr = -2 * self.mu * h
                            grad[off + i] += gr * u.real
                            grad[off + n + i] += gr * u.imag
                            grad[off + j] -= gr * u.real
                            grad[off + n + j] -= gr * u.imag
        for kind, lim in (("perm", self.dperm), ("block", self.dblock)):
            if lim <= 0:
                continue
            Wabs = np.abs(U) ** 2
            if kind == "perm":
                val, sig = dist_genperm(U)
                G = np.zeros((n, n), complex)
                for i in range(n):
                    G[sig[i], i] = np.conj(U[i, sig[i]])
                scale = -2.0 / n
            else:
                val, (S, T2) = dist_blockdiag(U)
                mask = np.zeros((n, n))
                for i in S:
                    for j in range(n):
                        if j not in T2:
                            mask[i, j] = 1
                G = (mask * np.conj(U)).T
                scale = 2.0
            if val < lim:
                h = lim - val
                obj += self.mu * h * h
                Cp = -W @ (G @ self.U0) @ (I + V)
                grad[4 * n:4 * n + n * n] += (-2 * self.mu * h * scale * _gradX(Cp)).ravel()
        return obj, grad


def smooth_maximize(a, b, U, theta=None, taus=(0.1, 0.03, 0.01, 3e-3, 1e-3, 3e-4, 1e-4),
                    maxiter=400, **opts):
    """Annealed L-BFGS on the smooth objective. Returns (a, b, U, theta, info)."""
    n = len(a)
    if theta is None:
        theta = margin(a, b, U)["theta"]
    nev = 0
    for tau in taus:
        a, b = canonical_spectra(a, b)
        theta = margin(a, b, U)["theta"]
        prob = Problem(n, U0=U, tau=tau, **opts)
        x0 = pack(a, b, np.zeros((n, n)), theta)
        res = minimize(prob.fg, x0, jac=True, method="L-BFGS-B",
                       options=dict(maxiter=maxiter, gtol=1e-12, ftol=1e-15))
        a, b, U, theta = prob.config(res.x)
        nev += prob.nev
    return a, b, U, theta, dict(nev=nev)


def polish(a, b, U, theta=None, vary_spectra=True, band=0.05, rounds=6, maxiter=300):
    """Active-set polish: maximize t s.t. Re(e^{-i th}(d - z_sigma)) >= t s0 for
    sigma in a working set (SLSQP), where s0 = rms spread of the z_sigma at the
    start. With vary_spectra the spectra move under s(a,b) = s0 and fixed mean(b).
    The working set grows until no permutation outside it is violated.
    Returns (a, b, U, theta, info); the caller must re-evaluate m with margin()."""
    n = len(a)
    a, b = canonical_spectra(a, b)
    r0 = margin(a, b, U)
    theta = r0["theta"] if theta is None else theta
    P = perms(n).astype(np.intp)
    z = (a[None, :] + b[P]).prod(1)
    s0 = float(np.sqrt(np.mean(np.abs(z - z.mean()) ** 2)))
    proj = (np.exp(-1j * theta) * z).real
    work = set(np.nonzero(proj >= proj.max() - band * s0)[0].tolist())
    U0 = np.asarray(U, complex)
    nv = 4 * n + n * n + 1
    b_mean0 = b.mean()
    info = dict(rounds=0, working=len(work))
    best = (r0["m"], a, b, U, theta)
    for rnd in range(rounds):
        idx = np.array(sorted(work))
        y0 = np.concatenate([pack(a, b, np.zeros((n, n)), theta), [r0["support_gap"] / s0]])

        def parts(y):
            aa, bb, X, th = unpack(y[:nv], n)
            I = np.eye(n)
            K = skew_from_real(X)
            W = np.linalg.inv(I + K)
            V = (I - K) @ W
            UU = U0 @ V
            M = np.diag(aa) + UU @ np.diag(bb) @ UU.conj().T
            adj, d = _adjugate(M)
            return aa, bb, X, th, W, V, UU, adj, d

        def cons(y):
            aa, bb, X, th, W, V, UU, adj, d = parts(y)
            zz = (aa[None, :] + bb[P[idx]]).prod(1)
            return (np.exp(-1j * th) * (d - zz)).real / s0 - y[-1]

        def cons_jac(y):
            aa, bb, X, th, W, V, UU, adj, d = parts(y)
            e = np.exp(-1j * th)
            I = np.eye(n)
            Pw = P[idx]
            Fm = aa[None, :] + bb[Pw]
            J = np.zeros((len(idx), nv + 1))
            # d part (same for every row)
            gd = np.zeros(nv)
            ra, ia = _hol_to_real(e * np.diag(adj))
            rb, ib = _hol_to_real(e * np.diag(UU.conj().T @ adj @ UU))
            gd[0:n], gd[n:2 * n], gd[2 * n:3 * n], gd[3 * n:4 * n] = ra, ia, rb, ib
            adjV = U0.conj().T @ adj @ U0
            Db = np.diag(bb)
            C = -W @ (Db @ V.conj().T @ adjV) @ (I + V) + (I + V.conj().T) @ (adjV @ V @ Db) @ W.conj().T
            gd[4 * n:4 * n + n * n] = _gradX(e * C).ravel()
            gd[-1] = (e * d).imag
            J[:, :nv] = gd
            for r, sg in enumerate(Pw):
                F = Fm[r]
                zk = np.array([np.prod(np.delete(F, k)) for k in range(n)])
                ra, ia = _hol_to_real(e * zk)
                J[r, 0:n] -= ra
                J[r, n:2 * n] -= ia
                J[r, 2 * n + sg] -= ra
                J[r, 3 * n + sg] -= ia
                J[r, -2] -= (e * F.prod()).imag
            J[:, :nv] /= s0
            J[:, -1] = -1.0
            return J

        cl = [dict(type="ineq", fun=cons, jac=cons_jac)]
        if vary_spectra:
            def eq(y):
                aa, bb, _, _ = unpack(y[:nv], n)
                zz = (aa[None, :] + bb[P]).prod(1)
                s = np.sqrt(np.mean(np.abs(zz - zz.mean()) ** 2))
                m = bb.mean() - b_mean0
                return np.array([s / s0 - 1, m.real, m.imag])
            cl.append(dict(type="eq", fun=eq))
        else:
            def eqf(y):
                return y[:4 * n] - y0[:4 * n]
            cl.append(dict(type="eq", fun=eqf, jac=lambda y: np.eye(nv + 1)[:4 * n]))
        res = minimize(lambda y: -y[-1], y0, jac=lambda y: -np.eye(nv + 1)[-1],
                       method="SLSQP", constraints=cl,
                       options=dict(maxiter=maxiter, ftol=1e-16))
        aa, bb, X, th, W, V, UU, adj, d = parts(res.x)
        r = margin(aa, bb, UU)
        info["rounds"] = rnd + 1
        if r["m"] > best[0]:
            best = (r["m"], aa, bb, UU, th)
        zz = (aa[None, :] + bb[P]).prod(1)
        viol = (np.exp(-1j * th) * (d - zz)).real / s0 - res.x[-1]
        new = set(np.nonzero(viol < -1e-12)[0].tolist()) - work
        a, b, U0, theta, r0 = aa, bb, UU, th, r
        if not new:
            break
        work |= new
    info["working"] = len(work)
    _, a, b, U, theta = best
    return a, b, U, theta, info


def describe(a, b, U):
    """Margin plus distances from the trivial sets, JSON-ready."""
    r = margin(a, b, U)
    dp, sig = dist_genperm(U)
    db, st = dist_blockdiag(U)
    U = np.asarray(U, complex)
    return dict(m=r["m"], theta=r["theta"], active=r["active"][:50], n_active=len(r["active"]),
                d=[r["d"].real, r["d"].imag], diam=r["diam"], unitarity_defect=r["unitarity_defect"],
                ok=r["ok"], degenerate=r["degenerate"], dist_genperm=dp, dist_blockdiag=db,
                block_split=[list(st[0]), list(st[1])], sep=min_separation(a, b),
                a=[[float(v.real), float(v.imag)] for v in a], b=[[float(v.real), float(v.imag)] for v in b],
                U_re=U.real.tolist(), U_im=U.imag.tolist())


def run_restarts(n, restarts, seed, out_json=None, kind="mixed", do_polish=True,
                 time_limit=None, keep=20, taus=None, **opts):
    """Restart driver. Restart r uses rng = default_rng([seed, r]): random
    spectra (kind) and a Haar unitary, then smooth_maximize and (optionally)
    polish. Logs a summary and the `keep` best configurations as JSON."""
    t0 = time.time()
    recs, ms = [], []
    for r in range(restarts):
        if time_limit and time.time() - t0 > time_limit:
            break
        rng = np.random.default_rng([seed, r])
        a, b = random_spectra(n, rng, kind)
        U = random_unitary(n, rng)
        kw = dict(opts)
        if taus is not None:
            kw["taus"] = taus
        a, b, U, th, info = smooth_maximize(a, b, U, **kw)
        if do_polish:
            # polish is unconstrained: keep its result only if the constraint
            # options (sep, dperm, dblock) still hold at 0.9 of their limits
            try:
                pa, pb, pU, pth, pinfo = polish(a, b, U, th)
                okc = (dist_genperm(pU)[0] >= 0.9 * opts.get("dperm", 0)
                       and dist_blockdiag(pU)[0] >= 0.9 * opts.get("dblock", 0)
                       and min(min_separation(pa, pb)) >= 0.9 * opts.get("sep", 0))
                if okc:
                    a, b, U, th = pa, pb, pU, pth
            except (np.linalg.LinAlgError, ValueError):
                pass
        rec = describe(a, b, U)
        rec.update(seed=seed, restart=r)
        ms.append(rec["m"])
        recs.append(rec)
        recs.sort(key=lambda t: -t["m"] if t["m"] == t["m"] else 1e300)
        del recs[keep:]
    ms = np.array(ms, float)
    summary = dict(n=n, seed=seed, restarts_done=len(ms), kind=kind, polish=do_polish, opts=opts,
                   max_m=float(np.nanmax(ms)) if len(ms) else None,
                   n_m_gt_1e7=int(np.sum(ms > 1e-7)), n_m_gt_0=int(np.sum(ms > 0)),
                   quantiles=[float(q) for q in np.nanquantile(ms, [0, 0.5, 0.9, 0.99, 1])] if len(ms) else [],
                   n_nan=int(np.sum(np.isnan(ms))), seconds=time.time() - t0)
    out = dict(summary=summary, best=recs)
    if out_json:
        with open(out_json, "w") as f:
            json.dump(out, f, indent=1)
    return out
