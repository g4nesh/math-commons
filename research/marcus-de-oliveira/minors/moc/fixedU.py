"""Fixed-U eigenvalue optimisation (no compound-LP filter), for large n."""
import sys, json, time, numpy as np
sys.path.insert(0, 'moc')
from bsearch import Problem, haar, normalise
n, seed, nU, N, steps = map(int, sys.argv[1:6])
rng = np.random.default_rng(seed); pr = Problem(n); res = []; t0 = time.time()
for u in range(nU):
    U = np.broadcast_to(haar(1, n, rng)[0], (N, n, n)).copy()
    a = rng.standard_normal((N, n)) + 1j * rng.standard_normal((N, n)); b = rng.standard_normal((N, n)) + 1j * rng.standard_normal((N, n))
    a, b = normalise(a, b); th = rng.uniform(0, 2 * np.pi, N)
    ma = np.zeros_like(a); mb = np.zeros_like(b); mt = np.zeros(N)
    for t in range(1, steps + 1):
        frac = t / steps; tau = 0.05 * (1 - frac) + 0.0005 * frac
        J, (gU, gA, gB, gth) = pr.evaluate(U, a, b, th, tau)
        gn = np.sqrt((np.abs(gA)**2).sum(1) + (np.abs(gB)**2).sum(1) + gth**2) + 1e-30
        ma = 0.9 * ma + 0.1 * gA / gn[:, None]; mb = 0.9 * mb + 0.1 * gB / gn[:, None]; mt = 0.9 * mt + 0.1 * gth / gn
        st = 0.03 * (1 - 0.9 * frac); a, b, th = a + st * ma, b + st * mb, th + st * mt; a, b = normalise(a, b)
    Je, _ = pr.evaluate(U, a, b, th, None, grad=False); res.append(float(np.nanmax(Je)))
    print("  n=%d U %d/%d: max J %.3e | overall %.3e  [%.0fs]" % (n, u + 1, nU, res[-1], max(res), time.time() - t0), flush=True)
json.dump(res, open('moc/runs/fixedU_n%d_s%d.json' % (n, seed), 'w'))
