#!/usr/bin/env python3
"""Find facets of the minor polytope Q_n by random linear programs and certify each one exactly.

A vertex of the slice {g in C : sum g = 1} of the gap cone C = L intersected with the nonnegative orthant
is an extreme ray of C, that is, a facet of Q_n. Each LP optimum (floating point) only proposes a
vertex set Z; the facet is then certified in exact rational arithmetic: the points v_s, s in Z, have
affine rank dim Q_n - 1, and the exact affine function vanishing on them has one strict sign on every
other vertex. Facets are reduced to orbits under the group of order 2 (n!)^2.

Usage:  python3 -I -B facet_search.py N SAMPLES SEED   (writes facets_nN.json next to this file)
This is a search: it proves that every listed facet is a facet, not that the list is complete.
"""
import itertools
import json
import os
import sys

import numpy as np
from flint import fmpq_mat
from scipy.optimize import linprog

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from minor_polytope import setup  # noqa: E402


def group(perms):
    n = len(perms[0])
    idx = {p: i for i, p in enumerate(perms)}
    inv = [idx[tuple(sorted(range(n), key=lambda i: s[i]))] for s in perms]
    rows = []
    for a in perms:
        for b in perms:
            row = [idx[tuple(a[s[b[i]]] for i in range(n))] for s in perms]   # s -> a o s o b
            rows.append(row)
            rows.append([row[inv[k]] for k in range(len(perms))])           # s -> a o s^-1 o b
    return np.array(rows, dtype=np.int16)


def canon(Z, G):
    ind = np.zeros((G.shape[0], G.shape[1]), dtype=np.uint8)
    np.put_along_axis(ind, G[:, sorted(Z)].astype(np.int64), 1, axis=1)
    packed = np.packbits(ind, axis=1)
    best = min(map(bytes, packed))
    return best.hex()


def certify(Z, pts, dimQ):
    """Exact: affine rank of Z is dimQ and the vanishing affine function is one-signed. Returns the
    values of that function on all vertices (integers, scaled), or None."""
    A = fmpq_mat([[1] + list(pts[s]) for s in sorted(Z)])
    if A.rank() != dimQ:
        return None
    # kernel of A: affine functions vanishing on Z (modulo those vanishing on all of Q)
    full = fmpq_mat([[1] + list(p) for p in pts])
    rref, rank = A.rref()
    ncols = A.ncols()
    piv, r = [], 0
    for c in range(ncols):
        if r < rank and rref[r, c] != 0:
            piv.append(c); r += 1
    free = [c for c in range(ncols) if c not in piv]
    for f in free:
        y = [0] * ncols
        y[f] = 1
        for i, c in enumerate(piv):
            y[c] = -rref[i, f]
        vals = [sum(y[c] * full[s, c] for c in range(ncols)) for s in range(len(pts))]
        if any(v != 0 for v in vals):
            pos = [v for v in vals if v != 0]
            if all(v > 0 for v in pos) or all(v < 0 for v in pos):
                zero = {s for s, v in enumerate(vals) if v == 0}
                return vals if zero == set(Z) else None
            return None
    return None


def main(n, samples, seed):
    perms, pairs, pts = setup(n)
    N = len(perms)
    dimQ = fmpq_mat([[1] + list(p) for p in pts]).rank() - 1
    minor_zero = {frozenset(s for s in range(N) if pts[s][c] == 0) for c in range(len(pairs))}
    G = group(perms)
    V = np.array(pts, dtype=float)
    B = np.hstack([V, np.ones((N, 1))])          # g = B y spans L
    rng = np.random.default_rng(seed)
    orbits = {}
    stats = {}
    for it in range(samples):
        w = rng.standard_normal(N)
        known = [o["zero_set"] for o in orbits.values()]
        if known and rng.random() < 0.7:
            # a random subset of a random known facet, moved by a random group element: explores neighbours
            Z0 = np.array(known[int(rng.integers(len(known)))])
            Z0 = G[int(rng.integers(G.shape[0]))][Z0]
            k = int(rng.integers(dimQ // 2, dimQ))
            T = rng.choice(Z0, size=min(k, len(Z0)), replace=False)
        else:
            T = rng.choice(N, size=int(rng.integers(0, dimQ)), replace=False)   # force g_T = 0
        Aeq = np.vstack([(np.ones(N) @ B)[None, :], B[T]])
        beq = np.concatenate([[1.0], np.zeros(len(T))])
        res = linprog(w @ B, A_ub=-B, b_ub=np.zeros(N), A_eq=Aeq, b_eq=beq,
                      bounds=[(None, None)] * B.shape[1], method="highs-ds")
        stats[res.status] = stats.get(res.status, 0) + 1
        if res.status != 0:
            continue
        g = B @ res.x
        Z = frozenset(int(s) for s in np.where(g < 1e-9)[0])
        if frozenset(Z) in minor_zero:
            key = "minor:" + canon(Z, G)
        else:
            key = canon(Z, G)
        if key in orbits:
            orbits[key]["hits"] += 1
            continue
        vals = certify(Z, pts, dimQ)
        orbits[key] = {"hits": 1, "vertices": len(Z), "minor": frozenset(Z) in minor_zero,
                       "certified": vals is not None, "zero_set": sorted(Z),
                       "first_sample": it}
        print(f"sample {it}: new orbit, {len(Z)} vertices, minor={frozenset(Z) in minor_zero}, certified={vals is not None}", flush=True)
    print("LP status counts (0 optimal, 2 infeasible):", stats)
    out = {"n": n, "samples": samples, "seed": seed, "dim_Q": dimQ, "group_order": int(G.shape[0]),
           "orbits": list(orbits.values())}
    with open(os.path.join(HERE, f"facets_n{n}.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print(f"dim Q = {dimQ}; {len(orbits)} facet orbits found in {samples} samples; "
          f"{sum(1 for o in orbits.values() if not o['minor'])} non-minor; "
          f"all certified: {all(o['certified'] for o in orbits.values())}")


if __name__ == "__main__":
    main(int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]))
