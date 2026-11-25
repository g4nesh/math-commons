"""Parameter-free certificate attempt for (*) in discriminant form.
(*)  <=>  [ g0 + a <= 0  ==>  4ab - (g0+a)^2 >= 0 ],   a = D[12,12], b = D[13,14].
Positivstellensatz ansatz:   4ab - (g0+a)^2  =  sigma0  -  (g0 + a) * sigma1 ,   sigma0, sigma1 sums of squares,
sigma0 of degree 4 and sigma1 of degree d1 in the entries of U (torus-symmetric, sampled identity)."""
import sys, time, itertools, numpy as np, cvxpy as cp, scipy.sparse as sps
from base import *
from sos2 import W, g0, a_, b_, bw
def build(d, U):
    mons = list(itertools.combinations_with_replacement(range(16), d)); blocks = {}
    for i, m in enumerate(mons):
        key = (tuple(np.bincount([v // 4 for v in m], minlength=4)), tuple(np.bincount([v % 4 for v in m], minlength=4))); blocks.setdefault(key, []).append(i)
    flat = np.array([u.reshape(-1) for u in U]); M = np.ones((len(U), len(mons)), complex)
    for i, m in enumerate(mons):
        for v in m: M[:, i] *= flat[:, v]
    slices = []; cols = []; pos = 0
    for b in blocks.values():
        n = len(b); iu = np.triu_indices(n)
        vals = np.real(np.conj(M[:, [b[p] for p in iu[0]]]) * M[:, [b[q] for q in iu[1]]]); vals[:, iu[0] != iu[1]] *= 2
        cols.append(vals); slices.append((pos, pos + len(iu[0]), n)); pos += len(iu[0])
    return np.hstack(cols), slices
def psd_cons(g, slices):
    cons = [g[[lo for (lo, hi, n) in slices if n == 1]] >= 0]
    for (lo, hi, n) in slices:
        if n == 1: continue
        iu = np.triu_indices(n); S = sps.lil_matrix((n * n, hi - lo))
        for k, (p, q) in enumerate(zip(*iu)): S[p * n + q, k] = 1; S[q * n + p, k] = 1
        cons.append(cp.reshape(sps.csr_matrix(S) @ g[lo:hi], (n, n), order='C') >> 0)
    return cons
if __name__ == '__main__':
    d1 = int(sys.argv[1]); K = int(sys.argv[2]); t0 = time.time()
    rng = np.random.default_rng(0); U = [haar(rng) for _ in range(K)]; X = np.array([x_U(u) for u in U])
    A4, s4 = build(4, U); A1, s1 = build(d1, U)
    w, V = np.linalg.eigh(A4 @ A4.T); keep = w > 1e-10 * w[-1]; P = V[:, keep].T
    print("degree-4 function space rank %d; sigma1 degree %d with %d unknowns  [%.0fs]" % (keep.sum(), d1, A1.shape[1], time.time() - t0), flush=True)
    for label, bv in [("D[13,14]", b_), ("D[13,23]", bw)]:
        h = X @ (g0 + a_); f = 4 * (X @ a_) * (X @ bv) - h**2
        g4 = cp.Variable(A4.shape[1]); g1 = cp.Variable(A1.shape[1]); gam = cp.Variable()
        cons = psd_cons(g4, s4) + psd_cons(g1, s1)
        cons.append((P @ A4) @ g4 - (P @ (h[:, None] * A1)) @ g1 == P @ f - gam * (P @ np.ones(K)))
        # residual of the part of h*sigma1 outside the degree-4 space must vanish as well: enforce on all samples in least squares sense
        prob = cp.Problem(cp.Maximize(gam), cons + [gam <= 1]); prob.solve(solver='CLARABEL')
        res = None
        if g4.value is not None: res = np.abs(A4 @ g4.value - h * (A1 @ g1.value) - (f - gam.value)).max()
        print("b = %s: best gamma = %s (%s); max residual of the identity on all samples %s  [%.0fs]" % (label, None if gam.value is None else "%+.3e" % gam.value, prob.status, None if res is None else "%.1e" % res, time.time() - t0), flush=True)
