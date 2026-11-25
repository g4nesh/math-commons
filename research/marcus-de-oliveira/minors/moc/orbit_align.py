import sys, itertools, numpy as np
sys.path.insert(0, 'moc')
from reduced4 import PERMS, VERTS
from scipy.optimize import minimize
F = np.load('moc/facets4.npy')
onv = np.array([(np.abs(VERTS @ f[:-1] + f[-1]) < 1e-9) for f in F])
exo = [frozenset(PERMS[j] for j in range(24) if onv[i, j]) for i in range(len(F)) if onv[i].sum() == 13]
exo = set(exo)
comp = lambda p, q: tuple(p[q[i]] for i in range(4))
inv = lambda p: tuple(np.argsort(p))
orbit = set(); start = next(iter(exo)); stack = [start]
while stack:
    f = stack.pop()
    if f in orbit: continue
    orbit.add(f)
    for al in PERMS:
        stack.append(frozenset(comp(al, s) for s in f)); stack.append(frozenset(comp(s, al) for s in f))
    stack.append(frozenset(inv(s) for s in f))
print("exotic facets: %d; orbit of one facet under left/right S4 x transpose: %d; orbit subset of facet set: %s" % (len(exo), len(orbit), orbit <= exo))
# alignment re-check with varied scales / log-normal magnitudes
P = np.array(PERMS)
w = np.array([2,1,1,1,0,1,1,0,0,0,0,1,1,0,1,0,0,0,0,0,1,1,0,0], float)
wc = w - w.mean(); wc /= np.linalg.norm(wc)
def negcos(x):
    a = x[0:4] + 1j * x[4:8]; b = x[8:12] + 1j * x[12:16]
    v = (np.exp(1j * x[16]) * np.prod(a[None] + b[P], axis=1)).real; vc = v - v.mean(); nv = np.linalg.norm(vc)
    return 1.0 if nv < 1e-300 else float(vc @ wc / nv)
rng = np.random.default_rng(123); best = 0
for r in range(600):
    sc = np.exp(rng.uniform(-3, 3, 17)); x0 = rng.standard_normal(17) * sc
    res = minimize(negcos, x0, method='Nelder-Mead' if r % 3 == 0 else 'BFGS', options={'maxiter': 20000, 'xatol': 1e-12, 'fatol': 1e-14})
    best = min(best, res.fun)
print("alignment re-check (600 restarts, log-uniform scales, NM+BFGS): max cos = %.8f" % (-best))
