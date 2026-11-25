"""SU(4) -> SO(6): the second compound in the real ('magic') basis.  Every coordinate of x(U) = (B, D) is then a
quadratic form in the 36 entries of R in SO(6).  This module builds those quadratic forms and checks them."""
import itertools, numpy as np
from base import *
def eps(S):
    p = list(S) + list(comp(S)); s = 1
    for i in range(4):
        for j in range(i + 1, 4):
            if p[i] > p[j]: s = -s
    return s
W = np.zeros((6, 6), complex)                 # columns: f+_M, f-_M  (M = 0,1,2)
for m, (S, Sc) in enumerate(MATCH):
    W[P2.index(S), 2 * m] = 1 / np.sqrt(2); W[P2.index(Sc), 2 * m] = eps(S) / np.sqrt(2)
    W[P2.index(S), 2 * m + 1] = 1j / np.sqrt(2); W[P2.index(Sc), 2 * m + 1] = -1j * eps(S) / np.sqrt(2)
def C2(U): return np.array([[np.linalg.det(U[np.ix_(S, T)]) for T in P2] for S in P2])
def R_of(U):
    U = U / np.linalg.det(U) ** 0.25
    R = W.conj().T @ C2(U) @ W
    assert np.abs(R.imag).max() < 1e-9, np.abs(R.imag).max()
    return R.real
def quad_D(S, T):
    """A (36x36 symmetric) with D[S,T] = vec(R)^T A vec(R),  vec(R)[6a+b] = R[a,b]"""
    v = np.array([W[P2.index(S), a] * np.conj(W[P2.index(T), b]) for a in range(6) for b in range(6)])
    return np.real(np.outer(v, v.conj()))
AD = {(S, T): quad_D(S, T) for S in P2 for T in P2}
def quad_B(s, j):
    """B[s,j] = (sum_{S contains s, T contains j} (D[S,T] - D[S,T^c]) + 1)/4  ->  (A, const)"""
    A = np.zeros((36, 36))
    for S in P2:
        if s not in S: continue
        for T in P2:
            if j not in T: continue
            A += (AD[(S, T)] - AD[(S, comp(T))]) / 4
    return A, 0.25
def quad_of(w, w0=0.0):
    """functional w.x + w0  ->  (Q, c0) with value vec(R)^T Q vec(R) + c0"""
    Q = np.zeros((36, 36)); c0 = w0
    for s in range(4):
        for j in range(4):
            if w[iB(s, j)]:
                A, c = quad_B(s, j); Q += w[iB(s, j)] * A; c0 += w[iB(s, j)] * c
    for S in P2:
        for T in P2:
            if w[iD(S, T)]: Q += w[iD(S, T)] * AD[(S, T)]
    return Q, c0
if __name__ == '__main__':
    rng = np.random.default_rng(0); err = 0; orth = 0
    for _ in range(50):
        U = haar(rng); R = R_of(U); v = R.reshape(-1); x = x_U(U)
        orth = max(orth, np.abs(R @ R.T - np.eye(6)).max(), abs(np.linalg.det(R) - 1))
        for S in P2:
            for T in P2: err = max(err, abs(v @ AD[(S, T)] @ v - x[iD(S, T)]))
        for s in range(4):
            for j in range(4):
                A, c = quad_B(s, j); err = max(err, abs(v @ A @ v + c - x[iB(s, j)]))
    print("R real orthogonal with det 1: max deviation %.1e ; quadratic forms reproduce all 52 coordinates: max error %.1e" % (orth, err))
