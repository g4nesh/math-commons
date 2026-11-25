"""High-accuracy zeros of f on U(4) (zero set {g0 = D[12,12] = D[13,14] = 0}).
Source: the sampled and walked zeros produced by ../star-exact (zeros_walk_t1.pkl), re-polished here by Newton on the exact
gradient, then filtered: keep a zero only if every cubic of ../star-exact/complement_t1.json vanishes there to < 1e-12."""
import sys, json, pickle, time, numpy as np
sys.path.insert(0, '../star-exact')
from fractions import Fraction as Fr
from polish import polish, grad
from zeros import mono as mono3
from exact_space import BLOCK_LIST as BL3
t0 = time.time()
Z = pickle.load(open('../star-exact/zeros_walk_t1.pkl', 'rb'))['Z']
comp3 = [np.array([[float(Fr(q)) for q in row] for row in c]) if c else None for c in json.load(open('../star-exact/complement_t1.json'))['complement']]
def resid3(U):
    m = mono3(U); return max(np.abs(c @ m[b]).max() for b, c in zip(BL3, comp3) if c is not None)
out = []; res = []
for U in Z:
    V = polish(polish(U, 1.0, iters=12), 1.0, iters=12); r = resid3(V); res.append(r)
    if r < 1e-12: out.append(V)
res = np.array(res)
print("re-polished %d zeros ; residual of the degree-3 vanishing conditions: median %.1e, max %.1e ; kept (< 1e-12): %d  [%.0fs]" % (len(Z), np.median(res), res.max(), len(out), time.time() - t0))
np.save('zeros_hq.npy', np.array(out))
