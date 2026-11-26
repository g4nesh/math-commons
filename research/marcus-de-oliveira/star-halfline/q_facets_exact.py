"""The facets of the compound polytope Q, in exact integer arithmetic.

Q = conv{x(P_pi)} has 24 vertices and dimension 13.  Affine coordinates: B[s,j] (s, j <= 3) and the pairing-matrix entries
n_ab (a, b <= 2).  A 13-subset of the vertices spans a facet iff the other 11 lie weakly on one side of its hyperplane and
not all on it.  The orientation determinants of all 14-subsets (14 x 14 integer matrices) are computed by fraction-free
elimination in int64 with overflow and divisibility checks (routine of ../star-degree4/chiro_exact.py).
Expected: exactly 82 facets = 16 (B[s,j] = 0) + 18 (D[S,T] = 0) + 48 (the orbit of g0).
(../pmsh/star/verify_rest.py part (2) does the same with floating-point determinants rounded to integers.)"""
import itertools, sys, time, numpy as np
from math import comb
sys.path.insert(0, '../pmsh/star'); sys.path.insert(0, '../star-degree4')
from base import *
from chiro_exact import bareiss
t0 = time.time()
Y = np.array([[1] + [int(v[iB(s, j)]) for s in range(3) for j in range(3)] + [int(v[iD(MATCH[a][0], T)] + v[iD(MATCH[a][0], Tc)]) for a in range(2) for (T, Tc) in MATCH[:2]] for v in VX], dtype=np.int64)
assert Y.shape == (24, 14)
subs = np.array(list(itertools.combinations(range(24), 14)), dtype=np.int64); det = np.zeros(len(subs), dtype=np.int64); biggest = 0
for lo in range(0, len(subs), 100000):
    det[lo:lo + 100000], bg = bareiss(Y[subs[lo:lo + 100000]]); biggest = max(biggest, bg)
print("14-subsets: %d ; exact determinants computed ; max |det| = %d ; largest intermediate %d (< 2^31)  [%.0fs]" % (len(subs), np.abs(det).max(), biggest, time.time() - t0), flush=True)
s = np.sign(det); k = 13
binom = np.array([[comb(v, j) for j in range(15)] for v in range(24)], dtype=np.int64)
A = np.stack([binom[subs[:, j], j + 1] for j in range(k)], 1); Bm = np.stack([binom[subs[:, j + 1], j + 1] for j in range(k)], 1)
cA = np.hstack([np.zeros((len(subs), 1), dtype=np.int64), np.cumsum(A, 1)]); cB = np.hstack([np.zeros((len(subs), 1), dtype=np.int64), np.cumsum(Bm, 1)])
N13 = comb(24, 13); pos = np.zeros(N13, dtype=np.int64); neg = np.zeros(N13, dtype=np.int64); zm = np.zeros(N13, dtype=np.int64)
for i in range(14):
    rk = cA[:, i] + (cB[:, k] - cB[:, i]); sg = s * (-1) ** i
    pos += np.bincount(rk, weights=(sg > 0), minlength=N13).astype(np.int64); neg += np.bincount(rk, weights=(sg < 0), minlength=N13).astype(np.int64)
    np.bitwise_or.at(zm, rk[sg == 0], (np.int64(1) << subs[:, i])[sg == 0])
T_all = np.array(list(itertools.combinations(range(24), 13)), dtype=np.int64)
trk = sum(binom[T_all[:, j], j + 1] for j in range(13)); tmask = np.zeros(N13, dtype=np.int64); tmask[trk] = (np.int64(1) << T_all).sum(1)
one = (pos > 0) ^ (neg > 0); facets = set((tmask[one] | zm[one]).tolist())
zmask_of = lambda w: int(sum(1 << v for v in range(24) if VX[v] @ w == 0))
EX = list(exotic_facets().values()); ZF = {zmask_of(w) for w in EX}
ZB = {zmask_of(np.eye(NX, dtype=int)[iB(s_, j)]) for s_ in range(4) for j in range(4)}; ZD = {zmask_of(np.eye(NX, dtype=int)[iD(S, T)]) for S in P2 for T in P2}
ok = facets == (ZF | ZB | ZD) and len(ZF) == 48 and len(ZB) == 16 and len(ZD) == 18
print("facets of Q: %d ; equal to the union of %d (B = 0) + %d (D = 0) + %d (exotic): %s  [%.0fs]" % (len(facets), len(ZB), len(ZD), len(ZF), ok, time.time() - t0))
print("THE FACET LIST OF Q IS EXACT" if ok else "MISMATCH")
