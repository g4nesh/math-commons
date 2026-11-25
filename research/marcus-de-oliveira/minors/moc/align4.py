"""For n=4: can a MOC functional Re(e^{i th} f_{a,b}) align with an exotic facet normal?
On vertices: v_sigma = Re(e^{i th} z_sigma(a,b)).  Want v = C - t*w (t>0), w = facet weights."""
import sys; sys.path.insert(0,'moc')
import numpy as np, itertools
from scipy.optimize import minimize
from reduced4 import PERMS
P = np.array(PERMS)
w = np.array([2,1,1,1,0,1,1,0,0,0,0,1,1,0,1,0,0,0,0,0,1,1,0,0],float)
wc = w - w.mean(); wc /= np.linalg.norm(wc)
def zs(a,b):
    return np.prod(a[None,:]+b[P],axis=1)
def unpack(x):
    a = x[0:4]+1j*x[4:8]; b = x[8:12]+1j*x[12:16]; th=x[16]
    return a,b,th
def negcos(x):
    a,b,th = unpack(x)
    v = (np.exp(1j*th)*zs(a,b)).real
    vc = v - v.mean(); nv = np.linalg.norm(vc)
    if nv < 1e-300: return 1.0
    return float(vc@wc/nv)   # want this -> -1  (v ~ C - t w)
rng = np.random.default_rng(3)
best=[]
for r in range(400):
    x0 = rng.standard_normal(17)
    res = minimize(negcos, x0, method='BFGS', options={'maxiter':4000,'gtol':1e-12})
    best.append((res.fun, res.x))
best.sort(key=lambda t:t[0])
print("best cos similarities (want 1):", [round(-b[0],8) for b in best[:10]])
a,b,th = unpack(best[0][1])
# normalise: translate so that a has mean 0, scale
print("a =", np.round(a,5)); print("b =", np.round(b,5))
np.save('moc/align_best.npy', best[0][1])
