"""Facial reduction for the t = 1 certificate.
1. From polished zeros: for each weight block, C_b = real coefficient vectors c with c . m_b(U) = 0 on the zero set
   (the cubics vanishing on the zero set).  Rationalise a row-reduced basis.
2. The reduced SDP  f = sum_b (C_b m_b)^* S_b (C_b m_b)  is solved in reduced_sdp.py and reduced_haar.py."""
import sys, time, json, pickle, numpy as np, cvxpy as cp
from fractions import Fraction as Fr
from numeric import *
from zeros import mono
def rref(A, tol=1e-6):
    A = A.copy(); r = 0; piv = []
    for c in range(A.shape[1]):
        if r == A.shape[0]: break
        k = r + np.argmax(np.abs(A[r:, c]))
        if abs(A[k, c]) < tol: continue
        A[[r, k]] = A[[k, r]]; A[r] /= A[r, c]
        for i in range(A.shape[0]):
            if i != r: A[i] -= A[i, c] * A[r]
        piv.append(c); r += 1
    return A[:r], piv
if __name__ == '__main__':
    t = float(sys.argv[1]); t0 = time.time()
    Z = pickle.load(open('zeros_polished_t%g.pkl' % t, 'rb'))['Z']; M = np.array([mono(U) for U in Z])
    comp = []; worst_round = 0; worst_res = 0; dens = set()
    for b in BLOCK_LIST:
        X = np.vstack([M[:, b].real, M[:, b].imag]); n = len(b)
        u, sv, vt = np.linalg.svd(X, full_matrices=True); sv = np.concatenate([sv, np.zeros(n)])[:n] / np.sqrt(len(Z))
        k = int((sv > 1e-5).sum()); Cn = vt[k:]
        if Cn.shape[0] == 0: comp.append([]); continue
        R, piv = rref(Cn)
        Rq = [[Fr(float(x)).limit_denominator(64) for x in row] for row in R]
        worst_round = max(worst_round, max(abs(float(q) - x) for rq, rr in zip(Rq, R) for q, x in zip(rq, rr)))
        Rf = np.array([[float(q) for q in row] for row in Rq])
        worst_res = max(worst_res, np.abs(Rf @ M[:, b].T).max()); dens |= {q.denominator for row in Rq for q in row}
        comp.append(Rq)
    rtot = sum(len(c) for c in comp)
    print("t = %g: complement dimension %d (kernel %d) ; max rounding move %.1e ; max |c.m(U)| over zeros %.1e ; denominators %s" % (t, rtot, 816 - rtot, worst_round, worst_res, sorted(dens)))
    hist = {}
    for b, c in zip(BLOCK_LIST, comp): hist[(len(b), len(c))] = hist.get((len(b), len(c)), 0) + 1
    print("   (block size, complement dim) counts:", dict(sorted(hist.items())))
    json.dump({'t': t, 'complement': [[[str(q) for q in row] for row in c] for c in comp]}, open('complement_t%g.json' % t, 'w'))
    # the reduced SDP itself is in reduced_sdp.py (rational sample) and reduced_haar.py (Haar sample)
