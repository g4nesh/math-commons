"""How the decomposition of one hard facet of P5(c) into facets of Q5 moves with the spectrum.

A Moebius map phi of the real line sends m(H) affinely onto m(phi(H)) (Jacobi's identity, and det(cH + d)
is constant on the isospectral set), and P(c) onto P(phi(c)) with the same vertex labels. So the vertex set
Z of a facet and its gap vector f (values on the 120 permutations, up to scale) depend only on the
Moebius class of c. Classes of 5 ordered points are parametrized by (0, 1, s, t, inf) with 1 < s < t;
a finite representative is c(s, t) = (0, 1/2, s/(s + 1), t/(t + 1), 1). The spectrum (0, 1/4, 1/2, 3/4, 1)
of orbits.json is (s, t) = (3, 9).

For an orbit of decompose.json this keeps Z fixed, recomputes the facet of P(c(s, t)) through Z (float),
checks that Z is still exactly its vertex set, and solves f = sum of the same Q5 facets with new weights
(least squares on the support found at (3, 9)). Prints the weights divided by the first non-minor weight,
and the residual.

usage: multipliers.py ORBIT s1,s2,... t1,t2,...
"""
import itertools
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from facet_search import PM  # noqa: E402

n = 5
P = PM(n)
PERMS = [tuple(p) for p in P.PA]


def xvec(K, J):
    return np.array([float(tuple(sorted(p[i] for i in K)) == tuple(J)) for p in PERMS])


def facet_gap(c, Z):
    """Gap vector (nonnegative, max 1) of the facet of P(c) through the vertices Z, or None."""
    V = np.c_[P.verts(c), -np.ones(len(PERMS))]
    _, s, vt = np.linalg.svd(V[Z])
    N = vt[int((s > 1e-10 * s[0]).sum()):].T
    W = V @ N
    u, s2, _ = np.linalg.svd(W, full_matrices=False)
    if len(s2) < 1 or (len(s2) > 1 and s2[1] > 1e-9 * s2[0]):
        return None  # Z does not span a hyperplane section of P(c)
    f = u[:, 0] * np.sign(u[:, 0][np.argmax(abs(u[:, 0]))])
    f = f / f.max()
    zero = set(np.flatnonzero(abs(f) < 1e-9).tolist())
    if f.min() < -1e-9 or zero != set(Z):
        return None  # Z is no longer the vertex set of a facet
    return f


def main():
    k = int(sys.argv[1])
    S = [float(x) for x in sys.argv[2].split(",")]
    T = [float(x) for x in sys.argv[3].split(",")]
    rec = json.load(open(os.path.join(HERE, "decompose.json")))["decompositions"][k]
    gap0 = np.array(json.load(open(os.path.join(HERE, "orbits.json")))["orbits"][k]["gap"])
    Z = np.flatnonzero(gap0 == 0)
    cols, names = [], []
    for t in rec["terms"]:
        if t["kind"] == "minor":
            cols.append(xvec(*t["pair"]))
            names.append("x" + "".join(map(str, t["pair"][0])) + ">" + "".join(map(str, t["pair"][1])))
        else:
            g = np.ones(len(PERMS))
            g[t["zero_set"]] = 0
            cols.append(None)  # filled below from the exact gap vector
            names.append(t["kind"])
    # exact gap vectors of the non-minor facets: recompute from their zero sets
    sys.path.insert(0, HERE)
    from decompose import gap_from_zero_set
    for i, t in enumerate(rec["terms"]):
        if t["kind"] != "minor":
            cols[i] = gap_from_zero_set(t["zero_set"]).astype(float)
    A = np.array(cols).T
    first = next(i for i, t in enumerate(rec["terms"]) if t["kind"] != "minor")
    print("orbit", k, "terms", names)
    for s in S:
        for t in T:
            if not 1 < s < t:
                continue
            c = np.array([0, 0.5, s / (s + 1), t / (t + 1), 1.0])
            f = facet_gap(c, Z)
            if f is None:
                print(f"s={s:g} t={t:g}: Z is not a facet")
                continue
            w, *_ = np.linalg.lstsq(A, f, rcond=None)
            res = np.abs(A @ w - f).max()
            print(f"s={s:g} t={t:g}: " + " ".join(f"{x / w[first]:.10g}" for x in w) + f"  residual {res:.1e}")


if __name__ == "__main__":
    main()
