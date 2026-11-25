"""Real-a certificate depth for general n (random + adversarial), and complex-a calibration."""
import sys, numpy as np
from scipy.linalg import expm
sys.path.insert(0, 'moc')
from herm_cert import HermCert
from moclib import haar_unitary
n, seed, R, iters = map(int, sys.argv[1:5])
hc = HermCert(n); rng = np.random.default_rng(seed)
def skew():
    H = rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n)); return (H - H.conj().T) / 2
rand = [hc.depth(haar_unitary(n, rng), rng.standard_normal(n))[0] for _ in range(40)]
print("n=%d: rank of Phi_a on affine span = %d; 40 random (U, real a): min s* = %.4f" % (n, hc.depth(haar_unitary(n, rng), rng.standard_normal(n))[1], min(rand)), flush=True)
best = np.inf
for rs in range(R):
    U = haar_unitary(n, rng); a = rng.standard_normal(n); s = hc.depth(U, a)[0]; stU, sta = 0.3, 0.5
    for it in range(iters):
        ok = False
        for _ in range(4):
            U2 = U @ expm(stU * skew() / 2); a2 = a + sta * rng.standard_normal(n); s2 = hc.depth(U2, a2)[0]
            if s2 < s: U, a, s, ok = U2, a2, s2, True; break
        stU = stU * 1.25 if ok else stU * 0.7; sta = sta * 1.25 if ok else sta * 0.7
        if stU < 1e-6: break
    best = min(best, s); print("  restart %d: s* = %.6f (best %.6f)" % (rs, s, best), flush=True)
    np.save('moc/runs/herm_adv_n%d_s%d_U.npy' % (n, seed), U); np.save('moc/runs/herm_adv_n%d_s%d_a.npy' % (n, seed), a)
print("RESULT n=%d real a: min s* = %.6f" % (n, best))
