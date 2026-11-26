#!/usr/bin/env python3
"""Settle census candidates whose gap cone F_T lies in span N_T but whose squared-minor vectors lack
private permutations: enumerate the extreme rays of F_T = span(N_T) cap R^S_+ exactly (cddlib, GMP)
and report any that is not a squared-minor vector.  Usage: python3 -I -B resolve.py  (reads census_n5.json)"""
import json, os, sys
from fractions import Fraction
import cdd.gmp as cdd
from flint import fmpq_mat
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from minor_polytope import setup  # noqa: E402

perms, pairs, pts = setup(5); N = len(perms)
c = json.load(open(os.path.join(HERE, "census_n5.json")))
for u in c.get("unresolved", []):
    T = u["T"]
    NT = [k for k in range(len(pairs)) if all(pts[s][k] == 0 for s in T)]
    rref, rank = fmpq_mat([[pts[s][k] for s in range(N)] for k in NT]).rref()
    basis = [[Fraction(int(rref[i, s].p), int(rref[i, s].q)) for s in range(N)] for i in range(rank)]
    H = [[Fraction(0)] + [b[s] for b in basis] for s in range(N) if any(b[s] != 0 for b in basis)]
    poly = cdd.polyhedron_from_matrix(cdd.matrix_from_array(H, rep_type=cdd.RepType.INEQUALITY))
    gens = cdd.copy_generators(poly)
    minor_vecs = {tuple(pts[s][k] for s in range(N)) for k in NT}
    rays, bad = 0, []
    for row in gens.array:
        if row[0] != 0:
            continue
        g = [sum(Fraction(x) * b[s] for x, b in zip(row[1:], basis)) for s in range(N)]
        m = min(v for v in g if v > 0); gn = tuple(v / m for v in g)
        rays += 1
        if gn not in {tuple(Fraction(x) for x in mv) for mv in minor_vecs}:
            bad.append(sorted(s for s in range(N) if g[s] == 0))
    print(u["T"], "rank", rank, "|N_T|", len(NT), "rays", rays, "non-minor rays", len(bad), [len(z) for z in bad][:5])
    u["extreme_rays"] = rays; u["non_minor_rays"] = bad
json.dump(c, open(os.path.join(HERE, "census_n5.json"), "w"), indent=1)
