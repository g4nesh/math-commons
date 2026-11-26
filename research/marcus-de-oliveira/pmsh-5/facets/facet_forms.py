"""Readable inequalities for the non-minor facets of Q5: a minimum-l1 representative of each facet's
affine function in the squared-minor coordinates x_{K->J} (|K| = |J| <= 2; x_{K->J} = |det U[K,J]|^2),
found by LP and checked exactly on all 120 vertices.

usage: facet_forms.py
"""
import json
import os
import sys
from fractions import Fraction

import numpy as np
from scipy.optimize import linprog

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from decompose import LABELS, PAIRS, STRATA, X, gap_from_zero_set  # noqa: E402


def name(K, J):
    return "x" + "".join(map(str, K)) + ">" + "".join(map(str, J))


def form(gap):
    A = np.c_[X, np.ones(len(X))]
    m = A.shape[1]
    cost = np.r_[np.ones(m - 1), 0.0, np.ones(m - 1), 0.0]
    r = linprog(cost, A_eq=np.c_[A, -A], b_eq=gap.astype(float), bounds=(0, None), method="highs")
    a = r.x[:m] - r.x[m:]
    q = [Fraction(float(x)).limit_denominator(100) for x in a]
    assert all(sum(q[j] * int(A[i, j]) for j in range(m)) == int(gap[i]) for i in range(len(gap))), "not exact"
    neg = [(-q[j], PAIRS[j]) for j in range(m - 1) if q[j] < 0]
    pos = [(q[j], PAIRS[j]) for j in range(m - 1) if q[j] > 0]
    lhs = " + ".join((f"{c}*" if c != 1 else "") + name(*p) for c, p in neg) or "0"
    rhs = " + ".join((f"{c}*" if c != 1 else "") + name(*p) for c, p in pos)
    if q[-1]:
        rhs += f" + {q[-1]}"
    return f"{lhs} <= {rhs}"


def main():
    fac = json.load(open(os.path.join(STRATA, "facets_n5.json")))
    reps = [(lab, o["zero_set"]) for lab, o in zip(LABELS, [o for o in fac["orbits"] if not o["minor"]])]
    reps += [(o["label"], o["zero_set"]) for o in json.load(open(os.path.join(HERE, "decompose.json")))["new_orbits"]]
    out = {}
    for lab, Z in reps:
        out[lab] = form(gap_from_zero_set(Z))
        print(f"{lab} ({len(Z)} vertices): {out[lab]}", flush=True)
    json.dump(out, open(os.path.join(HERE, "facet_forms.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
