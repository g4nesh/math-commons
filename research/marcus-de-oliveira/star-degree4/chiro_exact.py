"""Exact recomputation of the orientation determinants behind the facet enumerations.

../pmsh/star/chiro.py computes, for every 12-subset of the 24 projected vertices, the determinant of the 12 x 12 matrix
of homogeneous coordinates as a quadratic c0 + c1 rho + c2 rho^2, from floating-point determinants at rho = 0, 1, -1
rounded to integers.  Here the same three determinants are computed in exact integer arithmetic (fraction-free Bareiss
elimination in int64) and compared with the stored coefficients.

Exactness.  Before every multiplication the entries involved are checked to be below 2^31 in absolute value, so no int64
operation overflows; every Bareiss division is checked to leave remainder zero.  A matrix whose pivot column is zero
below the diagonal has determinant zero and is excluded from the checks from then on.

Needs ../pmsh/star/chiro_coeffs.npy (run `python3 chiro.py` in ../pmsh/star first; 8 seconds)."""
import itertools, sys, time, numpy as np
sys.path.insert(0, '../pmsh/star')
from base import *
t0 = time.time()
Bc = np.array([[v[iB(s, j)] for s in range(3) for j in range(3)] for v in VX], dtype=np.int64)
Nn = np.zeros((24, 3, 3), dtype=np.int64)
for k, v in enumerate(VX):
    for a, (S, Sc) in enumerate(MATCH):
        for b, (T, Tc) in enumerate(MATCH): Nn[k, a, b] = v[iD(S, T)] + v[iD(S, Tc)]
def pts(rho): return np.hstack([np.ones((24, 1), dtype=np.int64), Bc, (Nn[:, 0, 0] + rho * Nn[:, 0, 1])[:, None], (Nn[:, 1, 0] + rho * Nn[:, 1, 1])[:, None]])
LIM = np.int64(2**31)
def bareiss(M):
    """exact determinants of the integer matrices M[i] (shape N x n x n), int64 with overflow and divisibility checks"""
    N, n, _ = M.shape; M = M.copy(); sign = np.ones(N, dtype=np.int64); prev = np.ones(N, dtype=np.int64); dead = np.zeros(N, dtype=bool); biggest = 0
    for k in range(n - 1):
        col = M[:, k:, k] != 0; has = col.any(1); r = col.argmax(1) + k
        dead |= ~has; idx = np.flatnonzero(has & (r != k) & ~dead)
        if len(idx): tmp = M[idx, k].copy(); M[idx, k] = M[idx, r[idx]]; M[idx, r[idx]] = tmp; sign[idx] = -sign[idx]
        M[dead] = 0; live = ~dead
        m = int(np.abs(M[live, k:, k:]).max()) if live.any() else 0; biggest = max(biggest, m)
        assert m < LIM, "an intermediate exceeds 2^31: int64 products are not guaranteed exact"
        piv = M[:, k, k].copy(); piv[dead] = 1
        sub = M[:, k + 1:, k + 1:] * piv[:, None, None] - M[:, k + 1:, k:k + 1] * M[:, k:k + 1, k + 1:]
        q, rem = np.divmod(sub, prev[:, None, None]); assert not rem[live].any(), "a Bareiss division is not exact"
        M[:, k + 1:, k + 1:] = q; prev = piv
    det = sign * M[:, n - 1, n - 1]; det[dead] = 0
    return det, biggest
if __name__ == '__main__':
    # self-test of the routine against Python big-integer determinants
    import random
    from fractions import Fraction as Fr
    rng = random.Random(5); T = np.array([[[rng.randint(-3, 3) for _ in range(12)] for _ in range(12)] for _ in range(300)], dtype=np.int64)
    T[:40, :, 3] = 0; T[40:80, 5] = T[40:80, 2]; T[80:120, 0, 0] = 0
    def det_exact(A):
        A = [[Fr(int(x)) for x in row] for row in A]; n = len(A); d = Fr(1)
        for k in range(n):
            pr = next((i for i in range(k, n) if A[i][k] != 0), None)
            if pr is None: return 0
            if pr != k: A[k], A[pr] = A[pr], A[k]; d = -d
            d *= A[k][k]
            for i in range(k + 1, n):
                f = A[i][k] / A[k][k]
                if f:
                    for j in range(k, n): A[i][j] -= f * A[k][j]
        return int(d)
    dt, _ = bareiss(T); assert all(int(dt[i]) == det_exact(T[i]) for i in range(300)); print("self-test: 300 random 12 x 12 integer matrices (120 of them degenerate or needing pivoting) agree with big-integer elimination", flush=True)
    C = np.load('../pmsh/star/chiro_coeffs.npy'); subs = np.array(list(itertools.combinations(range(24), 12)), dtype=np.int64)
    assert len(subs) == len(C) == 2704156
    mism = 0; biggest = 0; P = {r: pts(r) for r in (0, 1, -1)}
    for lo in range(0, len(subs), 120000):
        idx = subs[lo:lo + 120000]; c = C[lo:lo + 120000]
        for r in (0, 1, -1):
            d, b = bareiss(P[r][idx]); biggest = max(biggest, b)
            mism += int((d != c[:, 0] + r * c[:, 1] + r * r * c[:, 2]).sum())
    print("12-subsets: %d ; exact determinants at rho = 0, 1, -1 that differ from the stored quadratics: %d ; largest intermediate %d (< 2^31)  [%.0fs]" % (len(subs), mism, biggest, time.time() - t0))
    print("STORED ORIENTATION QUADRATICS ARE EXACT" if mism == 0 else "MISMATCH")
