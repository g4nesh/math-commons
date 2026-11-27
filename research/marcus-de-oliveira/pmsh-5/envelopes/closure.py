"""Closure discs (README, Proposition M) and how much of U(5) they cover.

(1) Identity. Let U in U(5), R = C = {0,1,2}, and close output row 4 into input column 4 with feedback
zeta: V(zeta) = U[A,A] + U[A,4] zeta (1 - u zeta)^-1 U[4,A], A = {0,..,3}, u = U[4,4]. For |zeta| = 1,
V(zeta) is unitary and its block V(zeta)[R,C] = W + z x y^T with x = U[R,4], y = U[4,C],
z = zeta/(1 - u zeta). Every minor is alpha + z beta, so for each squared minor
  |alpha|^2 = mean over |zeta| = 1 of |alpha + z beta|^2  -  |beta|^2 / (1 - |u|^2).
Checked by quadrature.
(2) Coverage. Rotating the two outside rows and columns first, x runs over the range of I - WW* and y^T
over that of I - W^T conj(W). For the t-forms L_t of (E20), (E21), Proposition M proves L_t(W) >= 0
whenever some such (x, y) has L_t(|beta|^2) <= 0. Report the fraction of random (U, t) where the
minimum over (x, y) (BFGS, 8 starts) is <= 0, and the smallest L_t(W) among the rest.

usage: closure.py NUNITARIES SEED
"""
import itertools
import sys

import numpy as np
from scipy.optimize import minimize

SETS = [(S, T) for k in (1, 2, 3) for S in itertools.combinations(range(3), k) for T in itertools.combinations(range(3), k)]
IX = {st: i for i, st in enumerate(SETS)}


def minors(W):
    return np.array([np.linalg.det(W[np.ix_(S, T)]) for S, T in SETS])


def t_form(form, t):
    """Coefficients of L_t on the 19 squared minors: (E20) = min_t L_t, (E21) = min_t L_t."""
    c = np.zeros(len(SETS))

    def add(S, T, w):
        c[IX[(tuple(S), tuple(T))]] += w
    for S, T, w in (([0], [0], 1), ([0, 1], [0, 1], 1), ([0, 1], [1, 2], 1), ([1, 2], [0, 1], 1), ([0, 1, 2], [0, 1, 2], -1)):
        add(S, T, w)  # E-facet
    if form == 20:  # L_t = E + (1 + t) Y20 + A/t
        add([0, 1], [1, 2], 1 / t)
        for S, T, w in (([1], [1], 1), ([2], [1], 1), ([1, 2], [0, 1], -1), ([1, 2], [1, 2], -1)):
            add(S, T, (1 + t) * w)
    else:  # L_t = E + (1 + 1/t) A + t Y21
        add([0, 1], [0, 1], 1 + 1 / t)
        for S, T, w in (([1], [0], 1), ([1], [2], 1), ([0, 1], [0, 2], -1), ([1, 2], [0, 2], -1)):
            add(S, T, t * w)
    return c


def haar(n, rng):
    Z = rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))
    Q, R = np.linalg.qr(Z)
    return Q * (np.diag(R) / abs(np.diag(R)))


def check_identity(rng):
    err_u = err_m = 0.0
    for _ in range(50):
        U = haar(5, rng)
        A, u = list(range(4)), U[4, 4]
        W, x, y = U[:3, :3], U[:3, 4], U[4, :3]
        zetas = np.exp(2j * np.pi * np.arange(4096) / 4096)
        vals = []
        for zeta in zetas[::64]:
            V = U[np.ix_(A, A)] + np.outer(U[A, 4], U[4, A]) * zeta / (1 - u * zeta)
            err_u = max(err_u, np.abs(V @ V.conj().T - np.eye(4)).max())
            err_m = max(err_m, np.abs(V[:3, :3] - (W + np.outer(x, y) * zeta / (1 - u * zeta))).max())
        beta = minors(W + np.outer(x, y)) - minors(W)
        mean = np.mean([np.abs(minors(W + np.outer(x, y) * z / (1 - u * z))) ** 2 for z in zetas], axis=0)
        err_m = max(err_m, np.abs(np.abs(minors(W)) ** 2 - (mean - np.abs(beta) ** 2 / (1 - abs(u) ** 2))).max())
    print(f"(1) closure is unitary to {err_u:.1e}; block and mean-value identity to {err_m:.1e}")


def min_laplacian(W, c, rng):
    A, s, Bh = np.linalg.svd(W)
    Xp, Yp = A[:, 1:], Bh[1:, :].T
    m0 = minors(W)

    def f(p):
        x = Xp @ (p[0:2] + 1j * p[2:4])
        y = Yp @ (p[4:6] + 1j * p[6:8])
        beta = minors(W + np.outer(x, y)) - m0
        return c @ np.abs(beta) ** 2 / (np.linalg.norm(x) ** 2 * np.linalg.norm(y) ** 2)
    return min(minimize(f, rng.standard_normal(8), method="BFGS").fun for _ in range(8))


def main():
    nu, seed = int(sys.argv[1]), int(sys.argv[2])
    rng = np.random.default_rng(seed)
    check_identity(rng)
    for form in (20, 21):
        covered = total = 0
        worst = np.inf
        for _ in range(nu):
            W = haar(5, rng)[:3, :3]
            for t in (0.05, 0.3, 1.0, 3.0, 20.0):
                c = t_form(form, t)
                total += 1
                if min_laplacian(W, c, rng) <= 1e-12:
                    covered += 1
                else:
                    worst = min(worst, c @ np.abs(minors(W)) ** 2)
        print(f"(2) (E{form}): a superharmonic closure direction exists for {covered} of {total} (U, t); "
              f"smallest L_t(W) among the others {worst:.3f}", flush=True)


if __name__ == "__main__":
    main()
