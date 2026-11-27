"""Zeros of (E20) and (E21) at blocks of U(5) that are not blocks of U(4) and not singular: minimize over
W = V1 diag(1, s2, s3) V2 with s2, s3 confined to [0.05, 0.95], and report the zeros found with their
E-facet value and square-root partners A, Y.

usage: interior_zeros.py NSTARTS SEED
"""
import sys

import numpy as np
from scipy.optimize import minimize

from contraction_min import E_block, cay, partners, realness
from reduce_check import E20_block, E21_block


def W_of(x):
    s = 0.05 + 0.9 / (1 + np.exp(-np.clip(x[36:38], -50, 50)))
    return cay(x[:18]) @ np.diag([1.0, s[0], s[1]]) @ cay(x[18:36])


def main():
    nstarts, seed = int(sys.argv[1]), int(sys.argv[2])
    rng = np.random.default_rng(seed)
    for name, f in (("E20", E20_block), ("E21", E21_block)):
        res = sorted((minimize(lambda x: f(W_of(x)), rng.standard_normal(38), method="BFGS") for _ in range(nstarts)),
                     key=lambda r: r.fun)
        zeros = [r for r in res if r.fun < 1e-8]
        print(f"{name}: {len(zeros)} of {nstarts} starts end at a zero with 0.05 <= s3 <= s2 <= 0.95", flush=True)
        for r in zeros[:4]:
            W = W_of(r.x)
            A, Y = partners(W, name)
            print(f"    value {r.fun:+.1e}  singular values {np.round(np.linalg.svd(W, compute_uv=False), 3)}  "
                  f"|det W| {abs(np.linalg.det(W)):.3f}  E {E_block(W):+.1e}  A {A:.1e}  Y {Y:.1e}  realness {realness(W):.0e}",
                  flush=True)


if __name__ == "__main__":
    main()
