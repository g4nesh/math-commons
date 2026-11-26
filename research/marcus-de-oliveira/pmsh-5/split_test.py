"""Splitting test: can H be split along one row into Hermitian blocks of size n-1?

For Hermitian H = U D_c U* (n x n) and a row r, put rho_j = |U_rj|^2. A split is a family of Hermitian
H_j of size n-1 with spectrum c minus c_j such that, over subsets S of the other n-1 positions,
  (i)  sum_j rho_j det H_j[S,S]        = det H[S,S]
  (ii) sum_j rho_j c_j det H_j[S,S]    = det H[S+r, S+r].
A split plus PMSH(n-1) for every H_j gives a PMSH certificate for H (condition on the eigenvalue at r).
The diagonal case of (i), (ii) is the interlacing identity, so it always holds; the question is the
full identity. The script solves (i), (ii) by least squares over H_j = V_j D V_j*, V_j in U(n-1),
from random starts, trying rows until one succeeds; it reports the smallest residual (max norm).
usage: split_test.py n seed trials restarts
"""
import itertools
import sys

import numpy as np
from scipy.optimize import least_squares


def unitary(x, k):
    K = x[: k * k].reshape(k, k) + 1j * x[k * k:2 * k * k].reshape(k, k)
    Q, R = np.linalg.qr(K)
    return Q * (np.diag(R) / abs(np.diag(R)))

class Minors:
    def __init__(self, m):
        self.SUBS = [T for k in range(1, m + 1) for T in itertools.combinations(range(m), k)]
        byk = {}
        for i, T in enumerate(self.SUBS):
            byk.setdefault(len(T), []).append((i, T))
        self.IDX = [(np.array([i for i, _ in v]), np.array([T for _, T in v])) for v in byk.values()]
    def __call__(self, H):
        out = np.empty(len(self.SUBS))
        for ii, TT in self.IDX:
            out[ii] = np.linalg.det(H[TT[:, :, None], TT[:, None, :]]).real
        return out

def split_residual(H, c, U, row, restarts, rng, tol=1e-11):
    n = len(c); m1 = n - 1
    mm = Minors(m1)
    rest = [i for i in range(n) if i != row]
    rho = np.abs(U[row]) ** 2
    lhs1 = mm(H[np.ix_(rest, rest)])
    lhs2 = np.array([np.linalg.det(H[np.ix_([row] + [rest[s] for s in T], [row] + [rest[s] for s in T])]).real for T in mm.SUBS])
    specs = [np.delete(c, j) for j in range(n)]
    P = 2 * m1 * m1
    def res(z):
        M = [mm((V := unitary(z[j * P:(j + 1) * P], m1)) @ np.diag(specs[j]) @ V.conj().T) for j in range(n)]
        return np.r_[sum(rho[j] * M[j] for j in range(n)) - lhs1, sum(rho[j] * c[j] * M[j] for j in range(n)) - lhs2]
    best, bz = np.inf, None
    for r in range(restarts):
        sol = least_squares(res, rng.standard_normal(n * P), xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=20000)
        v = np.abs(sol.fun).max()
        if v < best: best, bz = v, sol.x
        if best < tol: break
    return best, bz


def main():
    n, seed, trials, restarts = map(int, sys.argv[1:5])
    gen = np.random.default_rng(seed)
    opt = np.random.default_rng(seed + 7777)
    worst = 0.0
    for t in range(trials):
        c = np.r_[0.0, np.sort(gen.random(n - 2)), 1.0]
        U = unitary(gen.standard_normal(2 * n * n), n)
        H = U @ np.diag(c) @ U.conj().T
        tried = []
        for row in range(n):
            b, _ = split_residual(H, c, U, row, restarts, opt)
            tried.append(b)
            if b < 1e-11:
                break
        worst = max(worst, min(tried))
        print(f"trial {t}: c={np.round(c, 4).tolist()} rows tried {len(tried)}, smallest residual {min(tried):.3e}", flush=True)
    print(f"worst over {trials} trials: {worst:.3e}")


if __name__ == "__main__":
    main()
