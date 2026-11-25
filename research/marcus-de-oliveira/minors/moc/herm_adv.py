"""Adversarial minimisation of the real-a certificate depth s*(U, a) over U in U(4) and real a.
Calibration: with complex a the same search should be able to find s* < 1 (since compound Birkhoff fails)."""
import sys, json, numpy as np, itertools
from scipy.linalg import expm
sys.path.insert(0, 'moc')
from herm_cert import HermCert, all_minor_sq, Phi
from moclib import haar_unitary
mode, seed, R, iters = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
hc = HermCert(4); rng = np.random.default_rng(seed)
starts = [np.load('moc/runs/deepest4.npy'), np.load('moc/runs/orth_best.npy')]
def depth_complex(U, a):
    # complex a: certificate must match Re and Im parts of Phi_a simultaneously
    import numpy as np
    from scipy.optimize import linprog
    def PhiC(xd): v = Phi(xd, a, 4); return np.concatenate([v.real, v.imag])
    V = np.array([PhiC(pd) for pd in hc.Pd]); p = PhiC(all_minor_sq(U)); c = V.mean(0); X = V - c
    u, sv, vt = np.linalg.svd(X, full_matrices=False); r = (sv > 1e-9 * sv[0]).sum(); B = vt[:r]
    Vr = X @ B.T; qr = (p - c) @ B.T; m = len(V)
    res = linprog(np.r_[np.zeros(m), -1.0], A_eq=np.vstack([np.hstack([Vr.T, -qr[:, None]]), np.r_[np.ones(m), 0][None]]),
                  b_eq=np.r_[np.zeros(r), 1.0], bounds=[(0, None)] * (m + 1), method='highs')
    return res.x[-1] if res.status == 0 else np.nan
def depth(U, a):
    return hc.depth(U, a)[0] if mode == 'real' else depth_complex(U, a)
def skew(n):
    H = rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n)); return (H - H.conj().T) / 2
best = (np.inf, None, None)
for rs in range(R):
    U = starts[rs % 2] if rs < 4 else haar_unitary(4, rng)
    a = rng.standard_normal(4) if mode == 'real' else rng.standard_normal(4) + 1j * rng.standard_normal(4)
    s = depth(U, a); stU, sta = 0.3, 0.5
    for it in range(iters):
        ok = False
        for _ in range(4):
            U2 = U @ expm(stU * skew(4) / 2); a2 = a + sta * (rng.standard_normal(4) if mode == 'real' else rng.standard_normal(4) + 1j * rng.standard_normal(4))
            s2 = depth(U2, a2)
            if s2 < s: U, a, s, ok = U2, a2, s2, True; break
        stU = stU * 1.25 if ok else stU * 0.7; sta = sta * 1.25 if ok else sta * 0.7
        if stU < 1e-6: break
    if s < best[0]: best = (s, U, a)
    print("  restart %d: s* = %.6f   (best %.6f)" % (rs, s, best[0]), flush=True)
print("RESULT mode=%s: min s* = %.6f" % (mode, best[0]))
np.save('moc/runs/herm_adv_%s_s%d_U.npy' % (mode, seed), best[1]); np.save('moc/runs/herm_adv_%s_s%d_a.npy' % (mode, seed), best[2])
