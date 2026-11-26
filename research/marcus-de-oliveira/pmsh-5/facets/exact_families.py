"""Exact check of the closed-form spectral families quoted in the README.

For an orbit and a family (non-minor facet with weight 1, squared minors x_{K->J} with weights that are
rational functions of the Moebius coordinates (s, t) of the spectrum (0, 1, s, t, inf)), this checks at
random rational 1 < s < t, in exact arithmetic, that
    F = g + sum_a w_a(s, t) x_a   (values on the 120 permutations)
is the gap vector of a facet of P(c(s, t)) with vertex set Z of orbits.json: F is an affine function of the
vertices v_pi(c) (rank test), F >= 0, and F vanishes exactly on Z. Here c(s, t) = (0, 1/2, s/(s+1), t/(t+1), 1).

usage: exact_families.py [npoints]
"""
import itertools
import json
import os
import sys

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from decompose import PERMS, gap_from_zero_set  # noqa: E402

s, t = sp.symbols("s t")
FAMILIES = {  # orbit: {term: weight}; the non-minor facet is the one used at (s, t) = (3, 9)
    15: {"x01>13": (t - s) / s, "x01>34": (t - 1) * (t - s) / (s * (s - 1)), "x12>04": (t - s) / (s - 1),
         "x13>04": (t - s) / (s - 1), "x14>04": (t - s) / (s - 1), "x23>01": t - 1, "x24>12": (t - s) / (s - 1),
         "x34>12": (t - s) / (s - 1)},
    20: {"x01>12": (t - s) / (s - 1), "x04>01": t * (t - 1) / (t - s), "x04>14": (t - 1) / (t - s),
         "x14>01": t - 1, "x23>02": (t - s) / s},
    21: {"x01>03": (t - 1) / t, "x01>23": (t - 1) / (t - s), "x02>12": (t - s) / (s - 1),
         "x04>12": (t - s) / (s - 1), "x13>02": (t - s) / s, "x24>01": t - 1},
}
SUBS = [T for k in range(1, 6) for T in itertools.combinations(range(5), k)]


def xcol(name):
    K, J = (tuple(int(ch) for ch in part) for part in name[1:].split(">"))
    return sp.Matrix([int(tuple(sorted(p[i] for i in K)) == J) for p in PERMS])


def main():
    npts = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    import random
    rnd = random.Random(0)
    orbits = json.load(open(os.path.join(HERE, "orbits.json")))["orbits"]
    dec = {r["orbit"]: r for r in json.load(open(os.path.join(HERE, "decompose.json")))["decompositions"]}
    ok = True
    for k, fam in FAMILIES.items():
        Z = [i for i, x in enumerate(orbits[k]["gap"]) if x == 0]
        (ex,) = [t_ for t_ in dec[k]["terms"] if t_["kind"] != "minor"]
        g = sp.Matrix(gap_from_zero_set(ex["zero_set"]).tolist())
        cols = {name: xcol(name) for name in fam}
        for _ in range(npts):
            sv = sp.Rational(rnd.randint(11, 60), 10)
            tv = sv * sp.Rational(rnd.randint(11, 60), 10)
            c = [0, sp.Rational(1, 2), sv / (sv + 1), tv / (tv + 1), 1]
            V = sp.Matrix([[sp.prod([c[p[i]] for i in T]) for T in SUBS] + [1] for p in PERMS])
            F = g + sum((w.subs({s: sv, t: tv}) * cols[name] for name, w in fam.items()), sp.zeros(len(PERMS), 1))
            good = V.rank() == V.row_join(F).rank() and min(F) >= 0 and [i for i in range(len(PERMS)) if F[i] == 0] == Z
            ok &= good
            print(f"orbit {k}, (s, t) = ({sv}, {tv}): family member is the facet through Z: {good}", flush=True)
    print("ALL FAMILY CHECKS PASSED" if ok else "SOME FAMILY CHECK FAILED")


if __name__ == "__main__":
    main()
