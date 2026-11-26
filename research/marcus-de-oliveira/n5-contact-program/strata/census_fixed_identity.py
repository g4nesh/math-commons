#!/usr/bin/env python3
"""Is there an 8-set of permutations meeting the support of every known non-minor facet of Q_5?
The group of order 28,800 is transitive on S_5, so such a set may be assumed to contain the identity;
supports containing the identity are then met.  Integer program (HiGHS).  Usage: python3 -I -B census_fixed_identity.py TIME_LIMIT"""
import json, os, sys
import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from facet_search import group  # noqa: E402
from minor_polytope import setup  # noqa: E402
perms, pairs, pts = setup(5); N = len(perms); G = group(perms)
ident = perms.index(tuple(range(5)))
f = json.load(open(os.path.join(HERE, "facets_n5.json")))
fam = set()
for o in f["orbits"]:
    if not o["minor"]:
        sup = sorted(set(range(N)) - set(o["zero_set"]))
        fam |= {tuple(sorted(G[g][sup])) for g in range(G.shape[0])}
fam = [s for s in fam if ident not in s]
A = np.zeros((len(fam), N))
for i, s in enumerate(fam):
    A[i, list(s)] = 1
lb = np.zeros(N); lb[ident] = 1
r = milp(np.ones(N), constraints=[LinearConstraint(A, lb=1), LinearConstraint(np.ones((1, N)), ub=8)],
         integrality=np.ones(N), bounds=Bounds(lb, 1), options={"time_limit": float(sys.argv[1])})
print(f"{len(fam)} supports avoid the identity; HiGHS: {r.message}; value {r.fun}; dual bound {getattr(r, 'mip_dual_bound', None)}")
if r.x is not None:
    print("8-set:", ["".join(map(str, perms[i])) for i in range(N) if r.x[i] > 0.5])
