"""Minimize the E-facet, (E20) and (E21), in their block forms (reduce_check.py), over 3x3 contractions W
with singular values (1, 1, s) (blocks of U(4)), (1, s, s') (blocks of U(5)) and (s, s', s'') (all
contractions, i.e. blocks of U(6)). At each of the best minimizers report the singular values, the E-facet
value, the two square-root partners A and Y, and how far W is from real up to diagonal phases (largest
normalized imaginary part of a 4-cycle product w_ij conj(w_kj) w_kl conj(w_il)).

usage: contraction_min.py NSTARTS SEED
"""
import itertools
import sys

import numpy as np
from scipy.optimize import minimize

from reduce_check import E20_block, E21_block, M


def E_block(W):
    return abs(W[0, 0]) ** 2 + M(W, [0, 1], [0, 1]) + M(W, [0, 1], [1, 2]) + M(W, [1, 2], [0, 1]) - abs(np.linalg.det(W)) ** 2


def partners(W, form):
    if form == "E20":
        return M(W, [0, 1], [1, 2]), abs(W[1, 1]) ** 2 + abs(W[2, 1]) ** 2 - M(W, [1, 2], [0, 1]) - M(W, [1, 2], [1, 2])
    return M(W, [0, 1], [0, 1]), abs(W[1, 0]) ** 2 + abs(W[1, 2]) ** 2 - M(W, [0, 1], [0, 2]) - M(W, [1, 2], [0, 2])


def cay(x):
    A = (x[:9] + 1j * x[9:18]).reshape(3, 3)
    A = A - A.conj().T
    return np.linalg.solve(np.eye(3) - A, np.eye(3) + A)


def contraction(x, units):
    s = 1 / (1 + np.exp(-np.clip(x[36:39], -50, 50)))
    s[:units] = 1.0
    return cay(x[:18]) @ np.diag(s) @ cay(x[18:36])


def realness(W):
    m = 0.0
    for i, k, j, l in itertools.product(range(3), repeat=4):
        p = W[i, j] * np.conj(W[k, j]) * W[k, l] * np.conj(W[i, l])
        if abs(p) > 1e-6:
            m = max(m, abs(p.imag) / abs(p))
    return m


def main():
    nstarts, seed = int(sys.argv[1]), int(sys.argv[2])
    rng = np.random.default_rng(seed)
    for name, f in (("E-facet", E_block), ("E20", E20_block), ("E21", E21_block)):
        for units, label in ((2, "blocks of U(4)"), (1, "blocks of U(5)"), (0, "all contractions")):
            res = sorted((minimize(lambda x: f(contraction(x, units)), rng.standard_normal(39), method="BFGS")
                          for _ in range(nstarts)), key=lambda r: r.fun)
            print(f"{name} over {label}: minimum {res[0].fun:+.3e}", flush=True)
            if name == "E-facet":
                continue
            for r in res[:3]:
                W = contraction(r.x, units)
                A, Y = partners(W, name)
                print(f"    value {r.fun:+.1e}  singular values {np.round(np.linalg.svd(W, compute_uv=False), 3)}  "
                      f"E {E_block(W):+.1e}  A {A:.1e}  Y {Y:.1e}  realness {realness(W):.0e}", flush=True)


if __name__ == "__main__":
    main()
