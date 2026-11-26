"""Structured MOC search: one Hermitian matrix against a normal one, and related
structured spectral classes. Maximizes the signed margin m of d(U) against
H = conv{z_sigma} over U and over spectra restricted to a class.

Classes (y = real parameter vector, a and b built from y):
  herm      a real; b has k nonreal entries (first k), the other n-k real.
  linecirc  a real; b = c + exp(rho + i phi)            (b on a circle).
  circles   a = exp(i alpha); b = c + exp(rho + i beta)  (both on circles).
  unitary   a = exp(i alpha); b = exp(rho + i beta)      (A, B unitary up to scale).
  a3r       a real with 3 distinct values (multiplicities `mult`); b complex.
  a3c       a complex with 3 distinct values; b complex.
  outa      a real, a[n-1] = R fixed; b complex.
  outb      a real; b[n-1] = R exp(i psi), rest of b complex.
  outc      a complex, a[n-1] = R fixed; b complex.

The objective and its gradient come from lib/moc_opt.Problem (chain rule through
the class map). Polish: lib/moc_opt.polish with the spectra fixed (U, theta only),
which keeps the class. Usage: see main().
"""
import argparse
import json
import math
import os
import sys
import time

import numpy as np
from scipy.optimize import minimize

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
from moc_core import margin, dist_genperm, dist_blockdiag, min_separation  # noqa: E402
from moc_param import skew_from_real, cayley, random_unitary  # noqa: E402
from moc_opt import Problem, pack, polish, describe  # noqa: E402
from robust import robust_margin  # noqa: E402


class Cls:
    def __init__(self, name, n, k=None, mult=None, R=None):
        self.name, self.n, self.k, self.R = name, n, k, R
        self.mult = mult
        if name in ("a3r", "a3c"):
            g = []
            for j, mlt in enumerate(mult):
                g += [j] * mlt
            assert len(g) == n and len(mult) == 3
            self.g = np.array(g)
        self.p = len(self.init(np.random.default_rng(0)))

    # ------------------------------------------------------------ map y -> (a, b)
    def ab(self, y):
        n, nm = self.n, self.name
        p = len(y)
        Ja = np.zeros((n, p), complex)
        Jb = np.zeros((n, p), complex)
        I = np.eye(n)
        if nm == "herm":
            k = self.k
            a = y[:n].astype(complex)
            bi = np.zeros(n)
            bi[:k] = y[2 * n:2 * n + k]
            b = y[n:2 * n] + 1j * bi
            Ja[:, :n] = I
            Jb[:, n:2 * n] = I
            Jb[:k, 2 * n:2 * n + k] = 1j * np.eye(k)
        elif nm in ("linecirc", "circles"):
            if nm == "linecirc":
                a = y[:n].astype(complex)
                Ja[:, :n] = I
            else:
                ea = np.exp(1j * y[:n])
                a = ea
                Ja[:, :n] = np.diag(1j * ea)
            cr, ci, rho = y[n], y[n + 1], y[n + 2]
            eb = np.exp(rho + 1j * y[n + 3:2 * n + 3])
            b = cr + 1j * ci + eb
            Jb[:, n] = 1
            Jb[:, n + 1] = 1j
            Jb[:, n + 2] = eb
            Jb[:, n + 3:2 * n + 3] = np.diag(1j * eb)
        elif nm == "unitary":
            ea = np.exp(1j * y[:n])
            a = ea
            Ja[:, :n] = np.diag(1j * ea)
            rho = y[n]
            eb = np.exp(rho + 1j * y[n + 1:2 * n + 1])
            b = eb
            Jb[:, n] = eb
            Jb[:, n + 1:2 * n + 1] = np.diag(1j * eb)
        elif nm in ("a3r", "a3c"):
            G = np.zeros((n, 3))
            G[np.arange(n), self.g] = 1
            if nm == "a3r":
                a = (G @ y[:3]).astype(complex)
                Ja[:, :3] = G
                o = 3
            else:
                a = G @ y[:3] + 1j * (G @ y[3:6])
                Ja[:, :3] = G
                Ja[:, 3:6] = 1j * G
                o = 6
            b = y[o:o + n] + 1j * y[o + n:o + 2 * n]
            Jb[:, o:o + n] = I
            Jb[:, o + n:o + 2 * n] = 1j * I
        elif nm == "outa":
            a = np.concatenate([y[:n - 1], [self.R]]).astype(complex)
            Ja[:n - 1, :n - 1] = np.eye(n - 1)
            o = n - 1
            b = y[o:o + n] + 1j * y[o + n:o + 2 * n]
            Jb[:, o:o + n] = I
            Jb[:, o + n:o + 2 * n] = 1j * I
        elif nm == "outb":
            a = y[:n].astype(complex)
            Ja[:, :n] = I
            m1 = n - 1
            br, bi, psi = y[n:n + m1], y[n + m1:n + 2 * m1], y[n + 2 * m1]
            eo = self.R * np.exp(1j * psi)
            b = np.concatenate([br + 1j * bi, [eo]])
            Jb[:m1, n:n + m1] = np.eye(m1)
            Jb[:m1, n + m1:n + 2 * m1] = 1j * np.eye(m1)
            Jb[m1, n + 2 * m1] = 1j * eo
        elif nm == "outc":
            m1 = n - 1
            a = np.concatenate([y[:m1] + 1j * y[m1:2 * m1], [self.R]])
            Ja[:m1, :m1] = np.eye(m1)
            Ja[:m1, m1:2 * m1] = 1j * np.eye(m1)
            o = 2 * m1
            b = y[o:o + n] + 1j * y[o + n:o + 2 * n]
            Jb[:, o:o + n] = I
            Jb[:, o + n:o + 2 * n] = 1j * I
        else:
            raise ValueError(nm)
        return a, b, Ja, Jb

    def init(self, rng):
        n, nm = self.n, self.name
        g = rng.standard_normal
        if nm == "herm":
            return np.concatenate([g(n), g(n), g(self.k)])
        if nm == "linecirc":
            return np.concatenate([g(n), 0.5 * g(2), [0.0], 2 * np.pi * rng.random(n)])
        if nm == "circles":
            return np.concatenate([2 * np.pi * rng.random(n), 0.5 * g(2), [0.3 * g()],
                                   2 * np.pi * rng.random(n)])
        if nm == "unitary":
            return np.concatenate([2 * np.pi * rng.random(n), [0.3 * g()], 2 * np.pi * rng.random(n)])
        if nm == "a3r":
            return np.concatenate([g(3), g(n), g(n)])
        if nm == "a3c":
            return np.concatenate([g(6), g(n), g(n)])
        if nm == "outa":
            return np.concatenate([g(n - 1), g(n), g(n)])
        if nm == "outb":
            return np.concatenate([g(n), g(n - 1), g(n - 1), [2 * np.pi * rng.random()]])
        if nm == "outc":
            return np.concatenate([g(n - 1), g(n - 1), g(n), g(n)])
        raise ValueError(nm)

    def renorm(self, y):
        """Structure-preserving shift/scale (m unchanged); keeps values O(1)."""
        n, nm = self.n, self.name
        y = y.copy()
        if nm == "herm":
            t = y[:n].mean()
            y[:n] -= t
            y[n:2 * n] += t
            a, b, _, _ = self.ab(y)
            c = 1.0 / math.sqrt((np.sum(np.abs(a) ** 2) + np.sum(np.abs(b - b.mean()) ** 2)) / (2 * n))
            y *= c
        elif nm == "linecirc":
            t = y[:n].mean()
            y[:n] -= t
            y[n] += t
            a, b, _, _ = self.ab(y)
            c = 1.0 / math.sqrt((np.sum(np.abs(a) ** 2) + np.sum(np.abs(b - b.mean()) ** 2)) / (2 * n))
            y[:n] *= c
            y[n:n + 2] *= c
            y[n + 2] += math.log(c)
        elif nm in ("a3r", "a3c"):
            nv = 3 if nm == "a3r" else 6
            a, b, _, _ = self.ab(y)
            t = a.mean() if nm == "a3r" else -b.mean()
            t = t.real if nm == "a3r" else t
            # a -> a - t, b -> b + t
            y[:3] -= np.real(t)
            if nm == "a3c":
                y[3:6] -= np.imag(t)
            y[nv:nv + n] += np.real(t)
            y[nv + n:nv + 2 * n] += np.imag(t)
            a, b, _, _ = self.ab(y)
            c = 1.0 / math.sqrt((np.sum(np.abs(a - a.mean()) ** 2) + np.sum(np.abs(b) ** 2)) / (2 * n))
            y *= c
        return y

    def describe(self):
        return dict(cls=self.name, n=self.n, k=self.k, mult=self.mult, R=self.R)

    def bounds(self):
        """L-BFGS-B bounds on w = [y, X, theta]: the log-radius rho stays in [-2.5, 2.5]."""
        bd = [(None, None)] * (self.p + self.n * self.n + 1)
        if self.name in ("linecirc", "circles"):
            bd[self.n + 2] = (-2.5, 2.5)
        elif self.name == "unitary":
            bd[self.n] = (-2.5, 2.5)
        return bd


class SProblem:
    def __init__(self, C, U0, tau, **opts):
        self.C = C
        self.P = Problem(C.n, U0=U0, tau=tau, **opts)

    def split(self, w):
        p, n = self.C.p, self.C.n
        return w[:p], w[p:p + n * n].reshape(n, n), w[-1]

    def fg(self, w):
        n = self.C.n
        y, X, th = self.split(w)
        a, b, Ja, Jb = self.C.ab(y)
        obj, g = self.P.fg(pack(a, b, X, th))
        gy = (Ja.real.T @ g[0:n] + Ja.imag.T @ g[n:2 * n]
              + Jb.real.T @ g[2 * n:3 * n] + Jb.imag.T @ g[3 * n:4 * n])
        return obj, np.concatenate([gy, g[4 * n:]])


def s_smooth(C, y, U, taus=(0.1, 0.03, 0.01, 3e-3, 1e-3, 3e-4, 1e-4), maxiter=400, **opts):
    n = C.n
    nev = 0
    for tau in taus:
        y = C.renorm(y)
        a, b, _, _ = C.ab(y)
        th = margin(a, b, U)["theta"]
        prob = SProblem(C, U, tau, **opts)
        w0 = np.concatenate([y, np.zeros(n * n), [th]])
        res = minimize(prob.fg, w0, jac=True, method="L-BFGS-B", bounds=C.bounds(),
                       options=dict(maxiter=maxiter, gtol=1e-12, ftol=1e-15))
        y, X, th = prob.split(res.x)
        U = U @ cayley(skew_from_real(X))
        nev += prob.P.nev
    return y, U, th, nev


def record(C, y, U, extra):
    a, b, _, _ = C.ab(y)
    rec = describe(a, b, U)
    rec.update(C.describe())
    rec["y"] = [float(v) for v in y]
    rec.update(extra)
    # vertex vs edge contact: distance of d from the nearest active z, relative to diam
    from moc_core import z_all
    z = z_all(a, b)
    d = complex(*rec["d"])
    # min over ALL n! points: 0 means d coincides with some z_sigma (point contact)
    rec["d_to_nearest_z"] = float(np.abs(d - z).min() / rec["diam"])
    # hull-independent guard (moc_core.margin can flip the sign when z points nearly coincide)
    rec["m_robust"] = robust_margin(a, b, U)["m"]
    return rec


def run(C, restarts, seed, out_json, time_limit, keep=25, do_polish=True, **opts):
    t0 = time.time()
    recs, ms, nevs = [], [], 0
    done = 0
    for r in range(restarts):
        if time.time() - t0 > time_limit:
            break
        rng = np.random.default_rng([seed, r])
        y = C.init(rng)
        U = random_unitary(C.n, rng)
        try:
            y, U, th, nev = s_smooth(C, y, U, **opts)
        except (np.linalg.LinAlgError, ValueError, FloatingPointError):
            continue
        nevs += nev
        a, b, _, _ = C.ab(y)
        m0 = robust_margin(a, b, U)["m"]
        pol = None
        if do_polish:
            try:
                _, _, pU, pth, pinfo = polish(a, b, U, th, vary_spectra=False)
                okc = (dist_genperm(pU)[0] >= 0.9 * opts.get("dperm", 0)
                       and dist_blockdiag(pU)[0] >= 0.9 * opts.get("dblock", 0))
                mp_ = margin(a, b, pU)
                pol = robust_margin(a, b, pU)["m"]
                if okc and mp_["unitarity_defect"] <= 1e-10 and pol > m0:
                    U = pU
            except (np.linalg.LinAlgError, ValueError):
                pass
        rec = record(C, y, U, dict(seed=seed, restart=r, m_smooth=m0, m_polish=pol))
        ms.append(rec["m_robust"])
        done += 1
        recs.append(rec)
        recs.sort(key=lambda t: -t["m_robust"] if t["m_robust"] == t["m_robust"] else 1e300)
        del recs[keep:]
    ms = np.array(ms, float)
    summary = dict(C.describe(), seed=seed, restarts_requested=restarts, restarts_done=done,
                   opts={k: v for k, v in opts.items()}, polish=do_polish, evaluations=nevs,
                   max_m_robust=float(np.nanmax(ms)) if len(ms) else None,
                   n_m_gt_1e7=int(np.sum(ms > 1e-7)), n_m_gt_1e10=int(np.sum(ms > 1e-10)),
                   quantiles=[float(q) for q in np.nanquantile(ms, [0, 0.5, 0.9, 0.99, 1])] if len(ms) else [],
                   seconds=time.time() - t0)
    out = dict(summary=summary, best=recs)
    with open(out_json, "w") as f:
        json.dump(out, f, indent=1)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cls", required=True)
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--k", type=int, default=None)
    ap.add_argument("--mult", type=str, default=None)
    ap.add_argument("--R", type=float, default=None)
    ap.add_argument("--restarts", type=int, default=100)
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--time", type=float, default=540)
    ap.add_argument("--dperm", type=float, default=0.0)
    ap.add_argument("--dblock", type=float, default=0.0)
    ap.add_argument("--maxiter", type=int, default=400)
    ap.add_argument("--sep", type=float, default=0.0)
    ap.add_argument("--out", required=True)
    A = ap.parse_args()
    mult = tuple(int(v) for v in A.mult.split(",")) if A.mult else None
    C = Cls(A.cls, A.n, k=A.k, mult=mult, R=A.R)
    opts = {}
    if A.dperm > 0:
        opts["dperm"] = A.dperm
    if A.dblock > 0:
        opts["dblock"] = A.dblock
    if A.sep > 0:
        opts["sep"] = A.sep
    out = run(C, A.restarts, A.seed, A.out, A.time, maxiter=A.maxiter, **opts)
    s = out["summary"]
    print(json.dumps({k: s[k] for k in ("cls", "n", "k", "mult", "R", "seed", "restarts_done", "evaluations",
                                        "max_m_robust", "n_m_gt_1e7", "n_m_gt_1e10", "quantiles", "seconds")}))
    for rec in out["best"][:5]:
        print("  mrob=%.3e m=%.3e smooth=%.3e dgp=%.3f dbd=%.3f nact=%d dz=%.2e sep=%s r=%d" % (
            rec["m_robust"], rec["m"], rec["m_smooth"], rec["dist_genperm"], rec["dist_blockdiag"], rec["n_active"],
            rec["d_to_nearest_z"], ["%.3f" % v for v in rec["sep"]], rec["restart"]))


if __name__ == "__main__":
    main()
