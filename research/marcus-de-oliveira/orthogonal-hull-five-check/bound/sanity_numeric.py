"""Floating-point sanity checks around the bound (not part of the proof).

  1. Phi(V) = 37/25 and V unitary, for the matrix V of task mdo-orthogonal-hull-five
     (found by the agent ganesh-claude-minors), typed in from its README.
  2. The identity  Phi(O) = 1 - 6 det(O) E_S(O)  on O(5), where S is the set of
     permutations with Phi(P_sigma) = -5 and E_S(O) = sum_{sigma in S} sgn(sigma) prod_i O[i, sigma(i)].
  3. Phi on random orthogonal matrices and a few local maximisations.

Usage: python sanity_numeric.py [seed]   (numpy and scipy)
"""
import itertools, sys
import numpy as np
from scipy.optimize import minimize
from scipy.linalg import expm

rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 1)

def minor(M, I, J):
    return np.linalg.det(M[np.ix_(sorted(I), sorted(J))])

TERMS = []
for v in range(5):
    TERMS.append((-1, (v,), (v,)))
    S1, S2 = ((v + 1) % 5, (v - 1) % 5), ((v + 2) % 5, (v - 2) % 5)
    TERMS += [(1, S1, S2), (1, S2, S1)]
    for d in (1, 2):
        TERMS += [(-1, (v, (v + d) % 5), (v, (v - d) % 5)), (-1, (v, (v - d) % 5), (v, (v + d) % 5))]

def phi(M):
    return sum(s * abs(minor(M, I, J)) ** 2 for s, I, J in TERMS)

def perm_matrix(p):
    P = np.zeros((5, 5)); P[range(5), p] = 1; return P

def sign(p):
    s, p = 1, list(p)
    for i in range(5):
        while p[i] != i:
            j = p[i]; p[i], p[j] = p[j], p[i]; s = -s
    return s

# 1. V
s6, t = np.sqrt(6), 1j * np.sqrt(15)
V = np.array([[4, 4*s6, 4*s6, 4*s6, 4*s6],
              [4*s6, -1-t, -6-2*t, -6+2*t, 9+t],
              [4*s6, -6-2*t, -1+t, 9-t, -6+2*t],
              [4*s6, -6+2*t, 9-t, -1+t, -6-2*t],
              [4*s6, 9+t, -6+2*t, -6-2*t, -1-t]]) / 20
print("|V V* - I|_max =", np.abs(V @ V.conj().T - np.eye(5)).max(), "; Phi(V) =", phi(V), "; 37/25 =", 37 / 25)

# 2. S and the identity
perms = list(itertools.permutations(range(5)))
vals = {p: round(phi(perm_matrix(p))) for p in perms}
S = [p for p in perms if vals[p] == -5]
print("|S| =", len(S), "; identity in S:", (0, 1, 2, 3, 4) in S, "; (0 1) in S:", (1, 0, 2, 3, 4) in S,
      "; even perms in S:", sum(sign(p) == 1 for p in S))
def E_S(O):
    return sum(sign(p) * np.prod(O[range(5), p]) for p in S)
def haar_O():
    Q, R = np.linalg.qr(rng.standard_normal((5, 5)))
    return Q * np.sign(np.diag(R))
err = 0
for _ in range(200):
    O = haar_O()
    err = max(err, abs(phi(O) - (1 - 6 * np.linalg.det(O) * E_S(O))))
print("max |Phi(O) - (1 - 6 det(O) E_S(O))| over 200 Haar O in O(5):", err)

# 3. random values and local maxima
best = max(phi(haar_O()) for _ in range(20000))
print("max Phi over 20000 Haar orthogonal samples:", best)
def skew(a):
    A = np.zeros((5, 5)); A[np.triu_indices(5, 1)] = a; return A - A.T
tops = []
for _ in range(20):
    O0 = haar_O()
    r = minimize(lambda a: -phi(O0 @ expm(skew(a))), np.zeros(10), method="BFGS")
    tops.append(-r.fun)
print("20 local maximisations over O(5): largest", max(tops), "smallest", min(tops))
