"""Minimize the spectrum-free envelope inequalities of orbits 20 and 21 over U(5).

With w = t - 1 and lambda = (t - s)/(t - 1) every weight of the families of orbits 20 and 21
(exact_families.py) is nondecreasing in w, so the infimum over spectra is the limit w -> 0 followed by
AM-GM in lambda:
  (E20)  g + x04>01 + x04>14 + 2 sqrt(x01>12 (x04>01 + x04>14)) >= 0,
  (E21)  g + x01>23 + 2 sqrt(x01>23 (x02>12 + x04>12)) >= 0,
with g the non-minor facet (orbit E) of each decomposition. The family of facets of P5(c) holds on U(5)
for every spectrum iff its envelope does. This reports the minimum found over U(5) and the block
structure of the minimizer.

usage: envelope_min.py NSTARTS
"""
import json
import os
import sys

import numpy as np
from scipy.optimize import minimize

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))
from decompose import PAIRS, gap_from_zero_set  # noqa: E402
from facet_search import blocks  # noqa: E402
from unitary_valid import affine_coeffs, dvec, unitary  # noqa: E402

IDX = {"x" + "".join(map(str, K)) + ">" + "".join(map(str, J)): i for i, (K, J) in enumerate(PAIRS)}
ENVELOPES = {
    20: lambda d: d["x04>01"] + d["x04>14"] + 2 * np.sqrt(d["x01>12"] * (d["x04>01"] + d["x04>14"])),
    21: lambda d: d["x01>23"] + 2 * np.sqrt(d["x01>23"] * (d["x02>12"] + d["x04>12"])),
}


def main():
    nstarts = int(sys.argv[1])
    rng = np.random.default_rng(3)
    dec = {r["orbit"]: r for r in json.load(open(os.path.join(HERE, "decompose.json")))["decompositions"]}
    for k, env in ENVELOPES.items():
        (ex,) = [t for t in dec[k]["terms"] if t["kind"] != "minor"]
        coef = affine_coeffs(gap_from_zero_set(ex["zero_set"]))

        def f(x):
            D = dvec(unitary(x))
            return coef @ D + env({name: max(D[i], 0.0) for name, i in IDX.items()})

        res = [minimize(f, rng.standard_normal(50), method="BFGS") for _ in range(nstarts)]
        best = min(res, key=lambda r: r.fun)
        gmin = coef @ dvec(unitary(best.x))
        print(f"orbit {k}: min of envelope over U(5) {best.fun:+.3e} (facet E alone there {gmin:+.3e}); "
              f"minimizer blocks {blocks(unitary(best.x), 1e-3)}", flush=True)


if __name__ == "__main__":
    main()
