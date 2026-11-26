"""How often is a Haar-random 5x5 unitary outside the orthogonal hull?  (floating point)

For U in U(5) the script evaluates Phi(P U Q) over all row and column permutations P, Q (only the 720 double cosets of the
symmetry group of Phi are needed) and records the largest value.  A value above 13/10 certifies, by Proposition B, that
r_5(U) is outside conv r_5(O(5)).  A value above 1 does the same if the numerically observed maximum 1 over O(5) is right.
usage: haar_fraction.py [seed] [number of unitaries]"""
import itertools, sys, numpy as np
from phi import TERMS
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 0); NU = int(sys.argv[2]) if len(sys.argv) > 2 else 300
P2 = list(itertools.combinations(range(5), 2)); IDX = {S: k for k, S in enumerate(P2)}
L1 = np.zeros((5, 5)); L2 = np.zeros((10, 10))
for S, T, c in TERMS:
    if len(S) == 1: L1[S[0], T[0]] += c
    else: L2[IDX[S], IDX[T]] += c
perms = [np.array(p) for p in itertools.permutations(range(5))]
def pair_perm(p): return np.array([IDX[tuple(sorted((p[a], p[b])))] for a, b in P2])
PP = [pair_perm(p) for p in perms]
# Phi(M[p][:, q]) = sum L1[i,j] B[p_i, q_j] + sum L2[S,T] D[p(S), q(T)]; use one representative q per coset of the affine group
aff = {tuple((a * i + b) % 5 for i in range(5)) for a in (1, 2, 3, 4) for b in range(5)}
reps = []; seen = set()
for q in itertools.permutations(range(5)):
    if q in seen: continue
    reps.append(q); seen.update(tuple(q[g[i]] for i in range(5)) for g in aff)
def haar():
    z = (rng.standard_normal((5, 5)) + 1j * rng.standard_normal((5, 5))) / np.sqrt(2); qq, r = np.linalg.qr(z); d = np.diag(r); return qq * (d / np.abs(d))
I = np.array(P2)
def best(U):
    B = np.abs(U)**2; a = U[I[:, 0]][:, I[:, 0]]; b = U[I[:, 1]][:, I[:, 1]]; c = U[I[:, 0]][:, I[:, 1]]; d = U[I[:, 1]][:, I[:, 0]]; D = np.abs(a * b - c * d)**2
    m = -np.inf
    for q in reps:
        qa = np.array(q); qp = pair_perm(qa); Bq = B[:, qa]; Dq = D[:, qp]
        for p, pp in zip(perms, PP): m = max(m, (L1 * Bq[p]).sum() + (L2 * Dq[pp]).sum())
    return m
vals = np.array([best(haar()) for _ in range(NU)])
print("%d Haar unitaries: largest relabelled Phi: median %.3f, maximum %.3f" % (NU, np.median(vals), vals.max()))
print("   above 13/10 (certified outside the orthogonal hull): %d ; above 1 (outside if the observed maximum over O(5) is right): %d" % ((vals > 1.3).sum(), (vals > 1 + 1e-9).sum()))
