"""PMSH(5) along the spectral family of one hard facet: maximize the facet over U(5) at many spectra.

The vertex set Z of a hard facet found at (s, t) = (3, 9) (orbits.json) spans a facet of P(c(s, t)) on a
whole region of Moebius classes (fit_families.py). For each (s, t) of a grid in that region this
recomputes the facet through Z at the shifted spectrum c(s, t) + 1 (normalized as in ../facet_classify.py:
l.(x - g) <= 1 with g the centroid of the vertices) and maximizes it over U(5) from NSTARTS random starts
with exact gradients (../facet_search.maximize_over_U). A value above 1 would be a violation of PMSH(5)
(to be confirmed independently); the README reports the largest value per orbit.

usage: family_min.py ORBIT[,ORBIT...] NSTARTS [seed]
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from facet_search import PM, maximize_over_U  # noqa: E402

P = PM(5)
GRID = [(s, s * r) for s in (1.1, 1.5, 3.0, 8.0, 30.0) for r in (1.1, 1.5, 3.0, 8.0, 30.0)]


def facet_at(Z, s, t):
    c = np.array([0, 0.5, s / (s + 1), t / (t + 1), 1.0]) + 1.0
    V = P.verts(c)
    g = V.mean(0)
    lh, *_ = np.linalg.lstsq(V[Z] - g, np.ones(len(Z)), rcond=None)
    vals = (V - g) @ lh
    if np.abs(vals[Z] - 1).max() > 1e-9 or np.delete(vals, Z).max() > 1 - 1e-9:
        return None
    return c, lh, g


def main():
    ks = [int(x) for x in sys.argv[1].split(",")]
    nstarts = int(sys.argv[2])
    rng = np.random.default_rng(int(sys.argv[3]) if len(sys.argv) > 3 else 0)
    orbits = json.load(open(os.path.join(HERE, "orbits.json")))["orbits"]
    for k in ks:
        Z = [i for i, x in enumerate(orbits[k]["gap"]) if x == 0]
        best, where, used = -np.inf, None, 0
        for s, t in GRID:
            fa = facet_at(Z, s, t)
            if fa is None:
                continue
            c, lh, g = fa
            used += 1
            val, _ = maximize_over_U(P, c, lh, g, nstarts, rng)
            if val > best:
                best, where = val, (s, t)
        print(f"orbit {k}: facet present at {used} of {len(GRID)} spectra; max over U(5) - 1 = {best - 1:+.2e} "
              f"at (s, t) = {where}", flush=True)


if __name__ == "__main__":
    main()
