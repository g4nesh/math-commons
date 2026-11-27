"""Two routes that do not work, with the data quoted in the README.

(1) Corners in a singular value. Fix V1, V2 in U(3) and s3, and let s2 run over [0, 1]: is the minimum of
(E20) or (E21) on W = V1 diag(1, s2, s3) V2 always at s2 in {0, 1}? (It would reduce blocks of U(5) to
blocks of U(4) and singular blocks.) Count lines with an interior minimum.
(2) Rank-one imaginary part. For a 3x3 contraction W, search diagonal phases with rank Im(D1 W D2) <= 1. Then
V = D1 W D2 = R + iN with R, N real and N of rank one, and for every tau > 0
  v(W) = v(V) = cos^2(th) v(R + tau N) + sin^2(th) v(R - N/tau),  tau = tan(th),
because each minor of R + zN is a + zb with a, b real. This reduces W to real matrices if both endpoints are
contractions for some tau. Report how often that holds (tau on a grid).

usage: failed_routes.py NLINES NMATRICES SEED
"""
import sys

import numpy as np
from scipy.optimize import minimize

from reduce_check import E20_block, E21_block, haar


def corners(nlines, rng):
    grid = np.linspace(0, 1, 101)
    for name, f in (("E20", E20_block), ("E21", E21_block)):
        interior = 0
        for _ in range(nlines):
            V1, V2, s3 = haar(3, rng), haar(3, rng), rng.uniform()
            vals = np.array([f(V1 @ np.diag([1, s, s3]) @ V2) for s in grid])
            interior += vals.min() < min(vals[0], vals[-1]) - 1e-9
        print(f"(1) {name}: interior minimum in s2 on {interior} of {nlines} lines", flush=True)


def rank_one_imaginary(nmat, rng):
    def V_of(W, ph):
        return np.exp(1j * np.r_[0.0, ph[:2]])[:, None] * W * np.exp(1j * ph[2:5])[None, :]

    def obj(ph, W):
        s = np.linalg.svd(V_of(W, ph).imag, compute_uv=False)
        return (s[1] ** 2 + s[2] ** 2) / (s[0] ** 2 + 1e-300)
    taus = np.exp(np.linspace(-6, 6, 2001))
    for units, label in ((0, "generic contractions"), (1, "blocks of U(5)"), (2, "blocks of U(4)")):
        found = good = 0
        worst = np.inf
        for _ in range(nmat):
            s = rng.uniform(size=3)
            s[:units] = 1.0
            W = haar(3, rng) @ np.diag(s) @ haar(3, rng)
            r = min((minimize(obj, rng.uniform(0, 2 * np.pi, 5), args=(W,), method="Nelder-Mead",
                              options=dict(xatol=1e-12, fatol=1e-20, maxiter=20000)) for _ in range(20)), key=lambda r: r.fun)
            if r.fun > 1e-14:
                continue
            found += 1
            V = V_of(W, r.x)
            Ul, sv, Vh = np.linalg.svd(V.imag)
            R, N = V.real, sv[0] * np.outer(Ul[:, 0], Vh[0])
            margin = max(min(1 - np.linalg.norm(R + t * N, 2), 1 - np.linalg.norm(R - N / t, 2)) for t in taus)
            good += margin >= -1e-9
            worst = min(worst, margin)
        print(f"(2) {label}: rank-one imaginary part found for {found} of {nmat}; both real endpoints contractions "
              f"for {good}; worst margin 1 - max norm {worst:+.1e}", flush=True)


if __name__ == "__main__":
    nlines, nmat, seed = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    rng = np.random.default_rng(seed)
    corners(nlines, rng)
    rank_one_imaginary(nmat, rng)
