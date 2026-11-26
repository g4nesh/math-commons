"""Which non-minor facets of Q5 stay nonnegative on the squared minors of every 5x5 unitary?

Every facet of Q5 is an affine function g on aff(Q5) (gap vector on the 120 permutations). Squared-minor
vectors D(U) of unitaries lie in aff(Q5) (exact_checks.py), so g(D(U)) is defined. This minimizes g(D(U))
over U(5) for one representative of each non-minor orbit: the seven orbits A-G of the strata README and
the orbits N1-N9 of decompose.json, from NSTARTS random starts (BFGS in the QR parametrization).
A facet with minimum < 0 fails on unitaries; minimum 0 (to rounding) is numerical evidence that it holds.

usage: unitary_valid.py NSTARTS [seed] [LABEL,LABEL,...]
"""
import json
import os
import sys

import numpy as np
from scipy.optimize import minimize

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from decompose import BASIS, PAIRS, LABELS, STRATA, gap_from_zero_set  # noqa: E402

n = 5
B = np.array(BASIS.tolist(), dtype=float)


def affine_coeffs(gap):
    """Coefficients (on PAIRS, then constant) of an affine function with these values on the vertices."""
    Xc = np.c_[np.array([[float(v) for v in row] for row in BASIS.tolist()])]
    y, res, *_ = np.linalg.lstsq(Xc, gap.astype(float), rcond=None)
    assert np.abs(Xc @ y - gap).max() < 1e-9
    # BASIS columns are a subset of the columns of [X, 1]: map back
    import sympy as sp
    from decompose import X
    full = sp.Matrix(np.c_[X, np.ones(len(X), int)].tolist())
    piv = list(full.rref()[1])
    coef = np.zeros(len(PAIRS) + 1)
    coef[piv] = y
    return coef


def unitary(x):
    K = x[:n * n].reshape(n, n) + 1j * x[n * n:].reshape(n, n)
    Q, R = np.linalg.qr(K)
    return Q * (np.diag(R) / abs(np.diag(R)))


_I = np.array([[K[0], K[-1], J[0], J[-1]] for K, J in PAIRS])  # |K| = 1 rows repeat the index
_ONE = np.array([len(K) == 1 for K, J in PAIRS])


def dvec(U):
    a, b = U[_I[:, 0], _I[:, 2]], U[_I[:, 1], _I[:, 3]]
    c, d = U[_I[:, 0], _I[:, 3]], U[_I[:, 1], _I[:, 2]]
    return np.r_[np.abs(np.where(_ONE, a, a * b - c * d)) ** 2, 1.0]


def main():
    nstarts = int(sys.argv[1])
    rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 0)
    fac = json.load(open(os.path.join(STRATA, "facets_n5.json")))
    reps = [(lab, o["zero_set"]) for lab, o in zip(LABELS, [o for o in fac["orbits"] if not o["minor"]])]
    reps += [(o["label"], o["zero_set"]) for o in json.load(open(os.path.join(HERE, "decompose.json")))["new_orbits"]]
    if len(sys.argv) > 3:
        reps = [r for r in reps if r[0] in sys.argv[3].split(",")]
    out = {}
    for lab, Z in reps:
        coef = affine_coeffs(gap_from_zero_set(Z))
        f = lambda x: coef @ dvec(unitary(x))  # noqa: E731
        vals = [minimize(f, rng.standard_normal(2 * n * n), method="BFGS").fun for _ in range(nstarts)]
        out[lab] = min(vals)
        print(f"{lab}: min over U(5) {min(vals):+.3e}  (starts below -1e-6: {sum(v < -1e-6 for v in vals)}/{nstarts})", flush=True)
    if len(sys.argv) <= 3:
        json.dump(out, open(os.path.join(HERE, "unitary_valid.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
