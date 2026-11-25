"""Targeted search for general n: only U whose full compound vector lies OUTSIDE conv{p_sigma}
can host a counterexample.  Find such U by LP, then maximise J over (a,b,theta) with U fixed."""
import sys, json, time, numpy as np
sys.path.insert(0, 'moc')
from bsearch import Problem, haar, normalise
from moclib import compound_perm_vectors, compound_birkhoff_lp

if __name__ == '__main__':
    n, seed, nU, N, steps = map(int, sys.argv[1:6])
    rng = np.random.default_rng(seed)
    perms, V = compound_perm_vectors(n)
    bad = []; tried = 0; t0 = time.time()
    while len(bad) < nU:
        U = haar(1, n, rng)[0]; tried += 1
        gap, _ = compound_birkhoff_lp(U, V)
        if gap > 1e-7: bad.append((gap, U))
    print("n=%d: %d bad U found among %d Haar samples (%.1f%%), gaps max %.2e  [%.0fs]" %
          (n, nU, tried, 100.0 * nU / tried, max(g for g, _ in bad), time.time() - t0), flush=True)
    pr = Problem(n); results = []
    for idx, (gap, U0) in enumerate(bad):
        U = np.broadcast_to(U0, (N, n, n)).copy()
        a = rng.standard_normal((N, n)) + 1j * rng.standard_normal((N, n))
        b = rng.standard_normal((N, n)) + 1j * rng.standard_normal((N, n))
        a, b = normalise(a, b); th = rng.uniform(0, 2 * np.pi, N)
        ma = np.zeros_like(a); mb = np.zeros_like(b); mt = np.zeros(N)
        for t in range(1, steps + 1):
            frac = t / steps; tau = 0.05 * (1 - frac) + 0.0005 * frac
            J, (gU, gA, gB, gth) = pr.evaluate(U, a, b, th, tau)
            gn = np.sqrt((np.abs(gA)**2).sum(1) + (np.abs(gB)**2).sum(1) + gth**2) + 1e-30
            ma = 0.9 * ma + 0.1 * gA / gn[:, None]; mb = 0.9 * mb + 0.1 * gB / gn[:, None]; mt = 0.9 * mt + 0.1 * gth / gn
            st = 0.03 * (1 - 0.9 * frac)
            a, b, th = a + st * ma, b + st * mb, th + st * mt
            a, b = normalise(a, b)
        Je, _ = pr.evaluate(U, a, b, th, None, grad=False)
        k = int(np.nanargmax(Je))
        results.append({'gap': gap, 'J': float(Je[k]), 'U_re': U0.real.tolist(), 'U_im': U0.imag.tolist(),
                        'a': [a[k].real.tolist(), a[k].imag.tolist()], 'b': [b[k].real.tolist(), b[k].imag.tolist()], 'th': float(th[k])})
        Js = np.array([r['J'] for r in results])
        print("  %d/%d bad U: this max J %.3e | overall max %.3e, median %.3e  [%.0fs]" %
              (idx + 1, nU, Je[k], Js.max(), np.median(Js), time.time() - t0), flush=True)
        json.dump(results, open('moc/runs/badn%d_s%d.json' % (n, seed), 'w'))
