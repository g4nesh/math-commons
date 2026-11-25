"""Batched adversarial search for MOC counterexamples with analytic gradients.
For each of N parallel starts we ascend
   J = (Re(e^{-i th} d) - smax_s Re(e^{-i th} z_s)) / W,   W = max v - min v,
over U in U(n) (via U <- U expm(H)), a, b in C^n and th.  J > 0 would be a counterexample."""
import sys, os, json, time, itertools
import numpy as np

def expm_skew(H):
    # H skew-Hermitian, batched: exp(H) = V exp(i w) V^H with iH Hermitian -> use eigh of (-1j*H)
    w, V = np.linalg.eigh(-1j * H)          # -iH Hermitian, H = i V diag(w) V^H
    return (V * np.exp(1j * w)[:, None, :]) @ np.conj(np.swapaxes(V, 1, 2))

def haar(N, n, rng):
    z = (rng.standard_normal((N, n, n)) + 1j * rng.standard_normal((N, n, n))) / np.sqrt(2)
    q, r = np.linalg.qr(z)
    d = np.diagonal(r, axis1=1, axis2=2)
    return q * (d / np.abs(d))[:, None, :]

class Problem:
    def __init__(self, n):
        self.n = n
        self.P = np.array(list(itertools.permutations(range(n))))
        self.Pinv = np.argsort(self.P, axis=1)

    def evaluate(self, U, a, b, th, tau, grad=True):
        n, P = self.n, self.P
        N = a.shape[0]
        M = np.einsum('kij,kj,klj->kil', U, b, np.conj(U))
        M[:, np.arange(n), np.arange(n)] += a
        dM = np.linalg.det(M)
        e = np.exp(-1j * th)
        F = a[:, None, :] + b[:, P]                 # (N, n!, n): factor for row i
        z = np.prod(F, axis=2)
        v = (e[:, None] * z).real
        dv = (e * dM).real
        vmax, vmin = v.max(1), v.min(1)
        W = vmax - vmin
        if tau is None:
            sm = vmax; pi = (v == vmax[:, None]).astype(float); pi /= pi.sum(1, keepdims=True)
            Lw = np.zeros(N)
        else:
            s = tau * W[:, None]
            ex = np.exp((v - vmax[:, None]) / s)
            sm = vmax + (s[:, 0]) * np.log(ex.sum(1)); pi = ex / ex.sum(1, keepdims=True)
            Lw = tau * (np.log(ex.sum(1)) - (pi * (v - vmax[:, None])).sum(1) / s[:, 0])
        J = (dv - sm) / W
        if not grad:
            return J, None
        Minv = np.linalg.inv(M)
        Y = np.conj(np.swapaxes(U, 1, 2)) @ Minv @ U
        ce = (e * dM)
        # d derivatives
        dd_a = ce[:, None] * np.diagonal(Minv, axis1=1, axis2=2)
        dd_b = ce[:, None] * np.diagonal(Y, axis1=1, axis2=2)
        G = b[:, :, None] * Y - Y * b[:, None, :]     # D_b Y - Y D_b
        X = ce[:, None, None] * G
        XH = np.conj(np.swapaxes(X, 1, 2))
        gH = (XH - X) / 2                              # skew-hermitian gradient of Re(e^{-ith} d)
        # z derivatives, weighted by pi : sum_s pi_s * e*z_s / F_si
        wz = pi * (e[:, None] * z)                     # (N, n!)
        R = wz[:, :, None] / F                         # (N, n!, n) d/d(factor i)
        dz_a = R.sum(1)
        dz_b = np.zeros_like(a)
        for j in range(n):
            # factor with b_j sits at row Pinv[s, j]
            dz_b[:, j] = R[np.arange(N)[:, None], np.arange(P.shape[0])[None, :], self.Pinv[None, :, j]].sum(1)
        # W derivative (subgradient): argmax - argmin
        imax, imin = v.argmax(1), v.argmin(1)
        def zgrad(idx):
            zz = e * z[np.arange(N), idx]
            Fi = F[np.arange(N), idx]                  # (N, n)
            ga = zz[:, None] / Fi
            gb = np.zeros_like(a)
            perm = P[idx]
            np.put_along_axis(gb, perm, ga, axis=1)
            return ga, gb
        ga1, gb1 = zgrad(imax); ga0, gb0 = zgrad(imin)
        dW_a, dW_b = ga1 - ga0, gb1 - gb0
        # theta derivatives
        dth_dv = (e * dM).imag
        dth_v = (e[:, None] * z).imag
        dth_sm = (pi * dth_v).sum(1)
        dth_W = dth_v[np.arange(N), imax] - dth_v[np.arange(N), imin]
        # assemble dJ (complex-valued derivatives c; ascent direction = conj(c))
        Wc = W[:, None]
        Lc = Lw[:, None]
        cA = (dd_a - dz_a - Lc * dW_a) / Wc - J[:, None] * dW_a / Wc
        cB = (dd_b - dz_b - Lc * dW_b) / Wc - J[:, None] * dW_b / Wc
        gth = (dth_dv - dth_sm - Lw * dth_W) / W - J * dth_W / W
        gU = gH / W[:, None, None]
        return J, (gU, np.conj(cA), np.conj(cB), gth)

def normalise(a, b):
    m = a.mean(1, keepdims=True)
    a = a - m; b = b + m
    s = np.maximum(np.abs(a).max(1, keepdims=True), np.abs(b).max(1, keepdims=True))
    return a / s, b / s

def run(n, N, steps, seed, out, lr=0.02):
    rng = np.random.default_rng(seed)
    pr = Problem(n)
    U = haar(N, n, rng)
    a = rng.standard_normal((N, n)) + 1j * rng.standard_normal((N, n))
    b = rng.standard_normal((N, n)) + 1j * rng.standard_normal((N, n))
    a, b = normalise(a, b)
    th = rng.uniform(0, 2 * np.pi, N)
    # Adam state
    mU = np.zeros_like(U); vU = np.zeros(N)
    ma = np.zeros_like(a); va = np.zeros(N); mb = np.zeros_like(b); vb = np.zeros(N)
    mt = np.zeros(N); vt = np.zeros(N)
    b1, b2 = 0.9, 0.999
    t0 = time.time(); hist = []
    for t in range(1, steps + 1):
        frac = t / steps
        tau = 0.05 * (1 - frac) + 0.0005 * frac
        J, (gU, gA, gB, gth) = pr.evaluate(U, a, b, th, tau)
        gnorm = np.sqrt((np.abs(gU)**2).sum((1, 2)) + (np.abs(gA)**2).sum(1) + (np.abs(gB)**2).sum(1) + gth**2) + 1e-30
        # simple normalised gradient ascent with momentum
        mU = b1 * mU + (1 - b1) * gU / gnorm[:, None, None]
        ma = b1 * ma + (1 - b1) * gA / gnorm[:, None]
        mb = b1 * mb + (1 - b1) * gB / gnorm[:, None]
        mt = b1 * mt + (1 - b1) * gth / gnorm
        step = lr * (1 - 0.9 * frac)
        H = step * mU
        H = (H - np.conj(np.swapaxes(H, 1, 2))) / 2
        U = U @ expm_skew(H)
        a = a + step * ma; b = b + step * mb; th = th + step * mt
        a, b = normalise(a, b)
        if t % 200 == 0 or t == steps:
            Je, _ = pr.evaluate(U, a, b, th, None, grad=False)
            k = int(np.nanargmax(Je))
            hist.append((t, float(np.nanmax(Je)), float(np.nanmedian(Je))))
            with open(out, 'w') as fh:
                json.dump({'n': n, 'N': N, 'seed': seed, 'step': t, 'hist': hist,
                           'best': float(Je[k]), 'U_re': U[k].real.tolist(), 'U_im': U[k].imag.tolist(),
                           'a': [a[k].real.tolist(), a[k].imag.tolist()],
                           'b': [b[k].real.tolist(), b[k].imag.tolist()], 'th': float(th[k]),
                           'top': np.sort(Je[np.isfinite(Je)])[-10:].tolist(),
                           'secs': time.time() - t0}, fh)
    return

if __name__ == '__main__':
    n, N, steps, seed = map(int, sys.argv[1:5])
    os.makedirs('moc/runs', exist_ok=True)
    run(n, N, steps, seed, 'moc/runs/b_n%d_s%d.json' % (n, seed))
