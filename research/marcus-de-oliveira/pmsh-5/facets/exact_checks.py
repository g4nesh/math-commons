"""Exact checks behind the README of this directory (integer and rational arithmetic, sympy).

1. R1 for n = 5: the squared-minor vector D(U) of every 5x5 unitary lies in aff(Q5). The affine relations
   of Q5 (in the coordinates x_{K,J}, |K| = |J| in {1, 2}, the others being complements) are checked to be
   spanned by identities valid for every unitary U:
     (i)  sum_T D[S,T] = 1 and sum_S D[S,T] = 1 (the compound matrices of U are unitary);
     (ii) sum_{T containing j} D[S,T] = sum_{s in S} |u_sj|^2 for |S| = 2, and the transposed identity
          (for the 2 x n matrix R = U[S,:] with orthonormal rows, sum_{T not containing j} |det R_T|^2 =
          det(I - r_j r_j^*) = 1 - |r_j|^2 by Cauchy-Binet; r_j is column j of R);
   with D[S,T] = D[S^c,T^c] (Jacobi's identity for unitary U) used to write 3x3 and 4x4 minors.
2. Every non-minor orbit N1-N9 of decompose.json is a facet of Q5: in the coordinates of ALL pairs
   (K, J) with 1 <= |K| = |J| <= 4 (complements kept as separate, equal coordinates), the vertex set has
   affine rank 41 = dim Q5 and the affine function vanishing on it is one-signed. (A code path independent
   of decompose.gap_from_zero_set, which works in a basis of 42 columns.)
3. Every decomposition of decompose.json is exact: sum of weight * (gap vector of the facet) equals the gap
   vector of orbits.json, with positive rational weights.

usage: exact_checks.py
"""
import itertools
import json
import os
import sys

import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from decompose import PAIRS, PERMS, X, gap_from_zero_set  # noqa: E402

n = 5


def check_r1():
    idx = {p: i for i, p in enumerate(PAIRS)}
    nc = len(PAIRS) + 1  # last column: constant
    rows = []

    def vec(terms, const):
        v = [0] * nc
        for (K, J), a in terms:
            v[idx[(K, J)]] += a
        v[-1] += const
        return v

    singles = [(i,) for i in range(n)]
    doubles = list(itertools.combinations(range(n), 2))
    for S in singles + doubles:
        rows.append(vec([((S, T), 1) for T in (singles if len(S) == 1 else doubles)], -1))
        rows.append(vec([((T, S), 1) for T in (singles if len(S) == 1 else doubles)], -1))
    for S in doubles:
        for j in range(n):
            rows.append(vec([((S, T), 1) for T in doubles if j in T] + [(((s,), (j,)), -1) for s in S], 0))
            rows.append(vec([((T, S), 1) for T in doubles if j in T] + [(((j,), (s,)), -1) for s in S], 0))
    I = sp.Matrix(rows)
    # all affine relations of Q5: left kernel of [X, 1]^T, i.e. vectors w with [X, 1] w = 0
    A = sp.Matrix(np.c_[X, np.ones(len(PERMS), int)].tolist())
    rel = A.nullspace()
    assert all(A * sp.Matrix(r) == sp.zeros(len(PERMS), 1) for r in rows), "an identity fails on a vertex"
    r_rel, r_id = len(rel), I.rank()
    r_both = sp.Matrix.vstack(I, sp.Matrix.hstack(*rel).T).rank()
    print(f"R1: aff(Q5) has {r_rel} independent affine relations (dim Q5 = {nc - 1 - r_rel}); the unitary "
          f"identities have rank {r_id}; together rank {r_both} -> {'spanned' if r_both == r_id == r_rel else 'NOT spanned'}")
    return r_both == r_id == r_rel


def check_new_facets():
    pairs = [(K, J) for k in range(1, n) for K in itertools.combinations(range(n), k)
             for J in itertools.combinations(range(n), k)]
    full = sp.Matrix([[1] + [int(tuple(sorted(p[i] for i in K)) == J) for K, J in pairs] for p in PERMS])
    dimQ = full.rank() - 1
    ok = True
    for o in json.load(open(os.path.join(HERE, "decompose.json")))["new_orbits"]:
        Z = o["zero_set"]
        A = full[Z, :]
        good = A.rank() == dimQ
        if good:
            # affine functions vanishing on Z, as value vectors on all vertices: must span one line
            V = sp.Matrix.hstack(*[full * v for v in A.nullspace()])
            vals = list(V[:, next(j for j in range(V.shape[1]) if any(V[:, j]))])
            good = V.rank() == 1 and (all(x >= 0 for x in vals) or all(x <= 0 for x in vals))
            good = good and sorted(i for i, x in enumerate(vals) if x == 0) == sorted(Z)
        ok &= good
        print(f"{o['label']}: {len(Z)} vertices, affine rank {A.rank() - 1}, facet of Q5: {good}")
    return ok


def check_decompositions():
    orbits = json.load(open(os.path.join(HERE, "orbits.json")))["orbits"]
    ok = True
    for rec in json.load(open(os.path.join(HERE, "decompose.json")))["decompositions"]:
        f = sp.Matrix(orbits[rec["orbit"]]["gap"])
        tot = sp.zeros(len(PERMS), 1)
        for t in rec["terms"]:
            w = sp.Rational(t["weight"])
            assert w > 0
            if t["kind"] == "minor":
                K, J = map(tuple, t["pair"])
                col = sp.Matrix([int(tuple(sorted(p[i] for i in K)) == J) for p in PERMS])
            else:
                col = sp.Matrix(gap_from_zero_set(t["zero_set"]).tolist())
            tot += w * col
        good = tot == f
        ok &= good
        print(f"orbit {rec['orbit']}: exact decomposition {good} "
              f"({', '.join(t['kind'] + ' ' + t['weight'] for t in rec['terms'] if t['kind'] != 'minor')}; "
              f"{sum(t['kind'] == 'minor' for t in rec['terms'])} squared minors)")
    return ok


if __name__ == "__main__":
    a = check_r1()
    b = check_new_facets()
    c = check_decompositions()
    print("ALL EXACT CHECKS PASSED" if a and b and c else "SOME CHECK FAILED")
