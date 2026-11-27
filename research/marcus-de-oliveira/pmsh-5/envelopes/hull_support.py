"""Hull property (H): the squared-minor vector of every 3x3 block of a 5x5 unitary lies in the convex hull of
those of 3x3 blocks of 4x4 unitaries and of singular 3x3 contractions. Equivalently, the support functions
agree: h1(c) = max c.v(W) over F1 = {contractions with singular values (1, s2, s3)} equals
h2(c) = max over F2 = {singular values (1, 1, s)} u {det W = 0}. v(W) lists all 19 squared minors of W.

Directions: Gaussian, and Gaussian with the 3x3 coefficient forced negative (the sign it has in every
facet that matters here). Each maximum uses BFGS from several random starts, so a positive gap larger than
the optimizer tolerance would refute (H); agreement is evidence only.

usage: hull_support.py SEED NDIRECTIONS
"""
import itertools
import sys

import numpy as np
from scipy.optimize import minimize

SETS = [(S, T) for k in (1, 2, 3) for S in itertools.combinations(range(3), k) for T in itertools.combinations(range(3), k)]


def v(W):
    return np.array([abs(np.linalg.det(W[np.ix_(S, T)])) ** 2 for S, T in SETS])


def cay(x):
    A = (x[:9] + 1j * x[9:18]).reshape(3, 3)
    A = A - A.conj().T
    return np.linalg.solve(np.eye(3) - A, np.eye(3) + A)


def sig(t):
    return 1 / (1 + np.exp(-np.clip(t, -50, 50)))


FAMILIES = {
    "F1": lambda x: cay(x[:18]) @ np.diag([1.0, sig(x[36]), sig(x[37])]) @ cay(x[18:36]),
    "U(4) block (1,1,s)": lambda x: cay(x[:18]) @ np.diag([1.0, 1.0, sig(x[36])]) @ cay(x[18:36]),
    "singular (1,s,0)": lambda x: cay(x[:18]) @ np.diag([1.0, sig(x[36]), 0.0]) @ cay(x[18:36]),
    "singular (s,t,0)": lambda x: cay(x[:18]) @ np.diag([sig(x[37]), sig(x[36]), 0.0]) @ cay(x[18:36]),
}
STARTS = {"F1": 12, "U(4) block (1,1,s)": 12, "singular (1,s,0)": 5, "singular (s,t,0)": 5}


def main():
    seed, ndir = int(sys.argv[1]), int(sys.argv[2])
    rng = np.random.default_rng(seed)
    worst = -np.inf
    for i in range(ndir):
        c = rng.standard_normal(len(SETS))
        if i % 2:
            c[-1] = -abs(c[-1]) - 1.0
        h = {}
        for name, f in FAMILIES.items():
            h[name] = max(-minimize(lambda x: -c @ v(f(x)), rng.standard_normal(38), method="BFGS").fun
                          for _ in range(STARTS[name]))
        h2 = max(val for name, val in h.items() if not name.startswith("F1"))
        gap = h["F1"] - h2
        worst = max(worst, gap)
        arg = max((val, name) for name, val in h.items() if not name.startswith("F1"))[1]
        print(f"seed {seed} direction {i}: h1 {h['F1']:.7f} h2 {h2:.7f} (attained on {arg}) gap {gap:+.1e}", flush=True)
    print(f"seed {seed}: largest gap h1 - h2 = {worst:+.1e}")


if __name__ == "__main__":
    main()
