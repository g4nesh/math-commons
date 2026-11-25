"""How far beyond the exotic facet (along the deepest violator's direction) can a point go before some
MOC functional detects it?  x_t = foot + t * n_F ; 'd' := sum_sigma mu_sigma(x_t) z_sigma(a,b)."""
import sys, itertools, numpy as np
from scipy.optimize import minimize
sys.path.insert(0, 'moc')
from reduced4 import VERTS, PERMS
from bad4 import red_batch, EXn
P = np.array(PERMS)
U = np.load('moc/runs/deepest4.npy'); x = red_batch(U[None])[0]
vals = x @ EXn[:, :-1].T + EXn[:, -1]; k = vals.argmax(); nF = EXn[k, :-1]; depth = vals[k]
foot = x - depth * nF
A = np.vstack([VERTS.T, np.ones(24)])
def mu_of(y):  # minimum-norm affine representation
    return np.linalg.lstsq(A, np.append(y, 1.0), rcond=None)[0]
def Jfun(p, mu):
    a = p[0:4] + 1j * p[4:8]; b = p[8:12] + 1j * p[12:16]; e = np.exp(-1j * p[16])
    z = np.prod(a[None] + b[P], axis=1); v = (e * z).real; W = v.max() - v.min()
    if W < 1e-12: return 1.0
    return -(((e * (mu @ z)).real - v.max()) / W)
rng = np.random.default_rng(0)
def maxJ(mu, R=150):
    best = -np.inf
    for r in range(R):
        res = minimize(Jfun, rng.standard_normal(17), args=(mu,), method='BFGS', options={'maxiter': 3000})
        best = max(best, -res.fun)
    return best
if __name__ == "__main__":
  print("deepest violator depth (Euclidean, reduced coords) = %.5f" % depth)
  for t in [depth, 0.06, 0.1, 0.15, 0.2, 0.3, 0.5]:
      mu = mu_of(foot + t * nF)
      print("  t = %.4f (%.1fx actual depth): max J found = %+.4e   (min mu = %.3f)" % (t, t / depth, maxJ(mu), mu.min()), flush=True)
