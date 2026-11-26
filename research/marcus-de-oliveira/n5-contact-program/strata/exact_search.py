#!/usr/bin/env python3
"""Exact search for an 8-set of permutations containing the identity that meets the support of every known
non-minor facet of Q_5 (the only 8-sets that can be automatic; the group is transitive, so containing
the identity is no restriction).  Depth-first branching on an unmet support of fewest available
elements; a branch excludes the elements tried by its earlier siblings, so each candidate set is
reached at most once.  Integer arithmetic only.  Exhausting the tree proves that no such 8-set exists.
Usage: python3 -I -B exact_search.py SECONDS"""
import json, os, sys, time
import numpy as np
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
fam = sorted(s for s in fam if ident not in s)
M = np.zeros((len(fam), N), dtype=bool)
for i, s in enumerate(fam):
    M[i, list(s)] = True
limit = float(sys.argv[1]); t0 = time.time(); nodes = [0]; found = []
top_done = [0, 0]


def dfs(unmet, banned, chosen, left):
    nodes[0] += 1
    if not unmet.any():
        found.append(list(chosen)); return True
    if left == 0:
        return False
    if time.time() - t0 > limit:
        raise TimeoutError
    avail = M[unmet] & ~banned
    sizes = avail.sum(1)
    if sizes.min() == 0:
        return False
    # packing bound: pairwise disjoint unmet supports (on available elements) each need their own pick
    order = np.argsort(sizes)
    used = np.zeros(N, dtype=bool); packed = 0
    for i in order:
        if not (avail[i] & used).any():
            used |= avail[i]; packed += 1
            if packed > left:
                return False
    row = avail[int(order[0])]
    banned = banned.copy()
    for e in np.flatnonzero(row):
        if dfs(unmet & ~M[:, e], banned, chosen + [int(e)], left - 1):
            return True
        banned[e] = True
        if len(chosen) == 1:
            top_done[0] += 1
    return False


unmet = np.ones(len(fam), dtype=bool); banned = np.zeros(N, dtype=bool); banned[ident] = True
first = M[unmet] & ~banned
top_done[1] = int(first.sum(1).min())
try:
    res = dfs(unmet, banned, [ident], 7)
    status = "found" if res else "exhausted: no 8-set containing the identity meets every known non-minor support"
except TimeoutError:
    status = "time limit"
print(f"{len(fam)} supports; {status}; nodes {nodes[0]}; top-level branches done {top_done[0]} of {top_done[1]}; {time.time()-t0:.0f}s")
if found:
    print("8-set:", ["".join(map(str, perms[i])) for i in found[0]])
