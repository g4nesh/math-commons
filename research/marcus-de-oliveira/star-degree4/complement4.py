"""Quartics vanishing on the zero set: for every weight block, the real coefficient vectors c with c . m_b(U) = 0 at all
high-accuracy zeros.  Row-reduced, rationalised, saved to complement4.json."""
import sys, json, time, numpy as np
from fractions import Fraction as Fr
from space4 import *
def rref(A, tol=1e-7):
    A = A.copy(); r = 0
    for c in range(A.shape[1]):
        if r == A.shape[0]: break
        k = r + np.argmax(np.abs(A[r:, c]))
        if abs(A[k, c]) < tol: continue
        A[[r, k]] = A[[k, r]]; A[r] /= A[r, c]
        for i in range(A.shape[0]):
            if i != r: A[i] -= A[i, c] * A[r]
        r += 1
    return A[:r]
if __name__ == '__main__':
    t0 = time.time(); Z = np.load('zeros_hq.npy'); M = mono_batch(Z); sizes = np.array([len(b) for b in BLOCK_LIST])
    allsv = []; svs = []
    for b in BLOCK_LIST:
        X = np.vstack([M[:, b].real, M[:, b].imag]); u, sv, vt = np.linalg.svd(X, full_matrices=True)
        sv = np.concatenate([sv, np.zeros(len(b))])[:len(b)] / np.sqrt(len(Z)); svs.append((sv, vt)); allsv.extend(sv.tolist())
    lg = np.log10(np.maximum(np.array(allsv), 1e-300))
    print("singular values (rms) of the monomial vectors on %d zeros, by decade:" % len(Z))
    for lo in range(-20, 1):
        c = ((lg >= lo) & (lg < lo + 1)).sum()
        if c: print("   1e%d .. 1e%d : %d" % (lo, lo + 1, c))
    print("   below 1e-20: %d" % (lg < -20).sum())
    comp = []; worst = 0; dens = set(); resid = 0
    for b, (sv, vt) in zip(BLOCK_LIST, svs):
        k = int((sv > 1e-6).sum()); Cn = vt[k:]
        if Cn.shape[0] == 0: comp.append([]); continue
        R = rref(Cn); Rq = [[Fr(float(x)).limit_denominator(64) for x in row] for row in R]
        worst = max(worst, max(abs(float(q) - x) for rq, rr in zip(Rq, R) for q, x in zip(rq, rr)))
        dens |= {q.denominator for row in Rq for q in row}
        resid = max(resid, np.abs(np.array([[float(q) for q in row] for row in Rq]) @ M[:, b].T).max()); comp.append(Rq)
    rtot = sum(len(c) for c in comp)
    print("complement dimension %d of %d (kernel %d) ; max rounding move %.1e ; denominators %s ; max |c.m(U)| over the zeros %.1e  [%.0fs]" % (rtot, len(MONS), len(MONS) - rtot, worst, sorted(dens), resid, time.time() - t0))
    hist = {}
    for b, c in zip(BLOCK_LIST, comp): hist[(len(b), len(c))] = hist.get((len(b), len(c)), 0) + 1
    print("(block size, complement dimension) counts:", dict(sorted(hist.items())))
    print("unknowns of the reduced problem: %d" % sum(len(c) * (len(c) + 1) // 2 for c in comp))
    json.dump({'complement': [[[str(q) for q in row] for row in c] for c in comp]}, open('complement4.json', 'w'))
