"""Rank of the relation space K4 at degree 4.

W4 : real space with basis  e_{P,Q} = Re(conj(m_P) m_Q),  P <= Q degree-4 monomials of equal torus weight.  dim W4 = 11208.
K4 : span of  Re( rho * conj(m_p') * m_q' ),  rho a homogeneous unitarity relation (rows/columns orthogonal, norms equal),
     m_p', m_q' degree-3 monomials.  Every element of K4 vanishes on U(4).
Claim checked here:  rank K4 >= R  (default R = 10254 = 11208 - 954).

Method.  Let Kmat be the (sparse, integer) generator matrix.  For random matrices A (R x #generators) and B (11208 x R)
over F_p, the R x R matrix  A * Kmat * B  is computed modulo p and shown to be invertible by LU factorisation modulo p.
Its rows are F_p-combinations of generators, so  rank_Q(K4) >= rank_p(Kmat) >= R.
Arithmetic: float64 holding exact integers; p < 2^19 keeps every intermediate below 2^53."""
import sys, time, itertools
from collections import defaultdict
import numpy as np, scipy.sparse as sps
from space4 import MONS, BLOCK_LIST, weight
P_MOD = 524287
def relations():
    rel = []
    for i in range(4):
        for j in range(4):
            if i != j:
                rel.append([(1, 4 * i + k, 4 * j + k) for k in range(4)]); rel.append([(1, 4 * k + i, 4 * k + j) for k in range(4)])
    for i in range(3):
        rel.append([(1, 4 * i + k, 4 * i + k) for k in range(4)] + [(-1, 4 * (i + 1) + k, 4 * (i + 1) + k) for k in range(4)])
        rel.append([(1, 4 * k + i, 4 * k + i) for k in range(4)] + [(-1, 4 * k + i + 1, 4 * k + i + 1) for k in range(4)])
    rel.append([(1, k, k) for k in range(4)] + [(-1, 4 * k, 4 * k) for k in range(4)])
    return rel
def generator_matrix():
    MIDX = {m: i for i, m in enumerate(MONS)}; PIDX = {}; col = 0
    for bl in BLOCK_LIST:
        for i in range(len(bl)):
            for j in range(i, len(bl)): PIDX[(bl[i], bl[j])] = col; col += 1
    M3 = list(itertools.combinations_with_replacement(range(16), 3)); W3 = [weight(m) for m in M3]; byw = defaultdict(list)
    for i, w in enumerate(W3): byw[w].append(i)
    rows = []; cols = []; vals = []; g = 0
    for rel in relations():
        c0, h0, a0 = rel[0]
        for qi, mq in enumerate(M3):
            wr = list(W3[qi][0]); wc = list(W3[qi][1])
            wr[h0 // 4] += 1; wc[h0 % 4] += 1; wr[a0 // 4] -= 1; wc[a0 % 4] -= 1
            if min(wr) < 0 or min(wc) < 0: continue
            for pi in byw.get((tuple(wr), tuple(wc)), []):
                mp = M3[pi]; row = {}
                for c, h, a in rel:
                    Q = MIDX[tuple(sorted(mq + (h,)))]; P = MIDX[tuple(sorted(mp + (a,)))]
                    k = PIDX[(P, Q) if P <= Q else (Q, P)]; row[k] = row.get(k, 0) + c
                for k, v in row.items():
                    if v: rows.append(g); cols.append(k); vals.append(v)
                g += 1
    return sps.csr_matrix((vals, (rows, cols)), shape=(g, col), dtype=np.float64)
def invertible_mod_p(B, p, nb=128, log=None):
    n = B.shape[0]
    for k0 in range(0, n, nb):
        k1 = min(k0 + nb, n)
        for k in range(k0, k1):
            nz = np.flatnonzero(B[k:, k])
            if len(nz) == 0: return False, k
            i = k + nz[0]
            if i != k: B[[k, i]] = B[[i, k]]
            inv = pow(int(B[k, k]), p - 2, p)
            B[k + 1:, k] = np.mod(B[k + 1:, k] * inv, p)
            if k + 1 < k1: B[k + 1:, k + 1:k1] = np.mod(B[k + 1:, k + 1:k1] - np.outer(B[k + 1:, k], B[k, k + 1:k1]), p)
        if k1 < n:
            for i in range(1, k1 - k0):
                B[k0 + i, k1:] = np.mod(B[k0 + i, k1:] - B[k0 + i, k0:k0 + i] @ B[k0:k0 + i, k1:], p)
            B[k1:, k1:] = np.mod(B[k1:, k1:] - B[k1:, k0:k1] @ B[k0:k1, k1:], p)
        if log and (k0 // nb) % 10 == 0: log(k1)
    return True, n
if __name__ == '__main__':
    R = int(sys.argv[1]) if len(sys.argv) > 1 else 10254; seed = int(sys.argv[2]) if len(sys.argv) > 2 else 1; t0 = time.time(); p = P_MOD
    Kmat = generator_matrix(); G, NWcols = Kmat.shape
    print("generators of K4: %d ; dim W4 = %d ; nonzeros %d ; max |entry| %d  [%.0fs]" % (G, NWcols, Kmat.nnz, int(np.abs(Kmat.data).max()), time.time() - t0), flush=True)
    rng = np.random.default_rng(seed); Om = rng.integers(0, p, size=(NWcols, R)).astype(np.float64); KT = Kmat.T.tocsr()
    B = np.empty((R, R)); maxabs = 0.0
    for lo in range(0, R, 1024):
        hi = min(lo + 1024, R); A = rng.integers(0, p, size=(G, hi - lo)).astype(np.float64)
        Y = KT @ A; maxabs = max(maxabs, np.abs(Y).max()); Y = np.mod(Y, p)                 # (NW x chunk): rows of A*Kmat, transposed
        Z = Y.T @ Om; maxabs = max(maxabs, np.abs(Z).max()); B[lo:hi] = np.mod(Z, p)
    print("B = A * Kmat * Omega assembled (%d x %d) ; largest intermediate %.3e (must stay below 2^53 = %.3e)  [%.0fs]" % (R, R, maxabs, 2.0**53, time.time() - t0), flush=True)
    assert maxabs < 2.0**53
    del Om
    ok, k = invertible_mod_p(B, p, log=lambda k1: print("   LU reached column %d  [%.0fs]" % (k1, time.time() - t0), flush=True))
    print("B invertible modulo %d: %s%s  [%.0fs]" % (p, ok, "" if ok else " (zero pivot column at %d)" % k, time.time() - t0))
    if ok: print("HENCE rank K4 >= %d ; with evaluation rank 954:  %d + 954 = %d = dim W4: %s" % (R, R, R + 954, R + 954 == NWcols))
