"""Shared exact data for the n=4 fiber analysis.
Coordinates x = (B, D): B[s,j] = |u_sj|^2 (16), D[S,T] = |det U[S,T]|^2 for 2-subsets S (rows), T (columns) (36).
Convention as in minors/pmsh.md: rows = positions, columns = eigen-indices."""
import itertools
import numpy as np
n = 4
P2 = list(itertools.combinations(range(n), 2))
MATCH = [((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))]
PERMS = list(itertools.permutations(range(n)))
NX = 52
def iB(s, j): return 4 * s + j
def iD(S, T): return 16 + 6 * P2.index(tuple(sorted(S))) + P2.index(tuple(sorted(T)))
def comp(S): return tuple(sorted(set(range(n)) - set(S)))
def mi(S):
    S = tuple(sorted(S))
    return next(k for k, (a, b) in enumerate(MATCH) if S in (a, b))
def x_perm(p):
    """integer compound vector of the permutation matrix with ones at (s, p[s])"""
    x = np.zeros(NX, dtype=int)
    for s in range(n): x[iB(s, p[s])] = 1
    for S in P2: x[iD(S, [p[s] for s in S])] = 1
    return x
def x_U(U):
    x = np.zeros(NX)
    for s in range(n):
        for j in range(n): x[iB(s, j)] = abs(U[s, j])**2
    for S in P2:
        for T in P2: x[iD(S, T)] = abs(np.linalg.det(U[np.ix_(S, T)]))**2
    return x
VX = np.array([x_perm(p) for p in PERMS])
def act(w, alpha, beta, transpose):
    """Functional w' with w'(x(U)) = w(x(P_alpha U P_beta)) (then transposed if requested)."""
    out = np.zeros(NX, dtype=w.dtype)
    for s in range(n):
        for j in range(n):
            s2, j2 = alpha[s], beta[j]
            out[iB(j2, s2) if transpose else iB(s2, j2)] += w[iB(s, j)]
    for S in P2:
        for T in P2:
            S2, T2 = [alpha[s] for s in S], [beta[t] for t in T]
            out[iD(T2, S2) if transpose else iD(S2, T2)] += w[iD(S, T)]
    return out
def exotic_facets():
    """The 48 exotic facets as integer functionals g_F(x) >= 0 on Q (vertex values 0,1,2), via the symmetry orbit of
    g0 = B11 + B33 + D[14,34] - D[12,14]  (1-based)."""
    g0 = np.zeros(NX, dtype=int)
    g0[iB(0, 0)] += 1; g0[iB(2, 2)] += 1; g0[iD((0, 3), (2, 3))] += 1; g0[iD((0, 1), (0, 3))] -= 1
    seen = {}
    for al in PERMS:
        for be in PERMS:
            for tr in (False, True):
                w = act(g0, al, be, tr)
                key = tuple(VX @ w)
                if key not in seen: seen[key] = w
    return seen
def n_matrix(x):
    """pairing matrix N[M,N'] = D[S,T] + D[S,T^c] and beta[M,N'] = D[S,T] - D[S,T^c] for the first-listed S, T."""
    Nm = np.zeros((3, 3)); be = np.zeros((3, 3))
    for a, (S, Sc) in enumerate(MATCH):
        for b, (T, Tc) in enumerate(MATCH):
            Nm[a, b] = x[iD(S, T)] + x[iD(S, Tc)]; be[a, b] = x[iD(S, T)] - x[iD(S, Tc)]
    return Nm, be
def shift(eta, q):
    """change of x under N -> N + eta q^T with B (hence all D[S,T]-D[S,T^c]) fixed"""
    dx = np.zeros(NX)
    for a, (S, Sc) in enumerate(MATCH):
        for b, (T, Tc) in enumerate(MATCH):
            for S_ in (S, Sc):
                for T_ in (T, Tc): dx[iD(S_, T_)] += eta[a] * q[b] / 2.0
    return dx
def haar(rng):
    z = (rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))) / np.sqrt(2)
    Q, R = np.linalg.qr(z); d = np.diag(R); return Q * (d / np.abs(d))
if __name__ == '__main__':
    ex = exotic_facets()
    print("exotic facets from the orbit of g0:", len(ex))
    vals = np.array(list(ex.keys()))
    print("vertex-value multisets:", {tuple(np.unique(v, return_counts=True)[1]) for v in vals})
