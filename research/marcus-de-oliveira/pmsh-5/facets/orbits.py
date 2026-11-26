"""Orbits of the hard facets of P5(c), c = (0, 1/4, 1/2, 3/4, 1), under positions and reversal.

The spectrum c + 1 = (1, 5/4, 3/2, 7/4, 2) used in ../facet_classify.py is 4c + 4 = (4, 5, 6, 7, 8) up to
scaling, so every vertex v_pi of P(c) is an integer vector and facets can be checked exactly.
Symmetries of P(c) acting on the 120 permutations (itertools order):
  positions  pi -> pi o sigma   (H -> P H P^T permutes the subsets T),
  reversal   pi -> r o pi with r(j) = 4 - j   (H -> 12 I - H; m(12 - H) is an integer-linear image of m(H)).
They generate a group of order 240. A facet is identified by its vertex set Z (the permutations where it is
tight); the group acts on vertex sets.

For each hard facet of ../hard_n5_seed41_c0_0.25_0.5_0.75_1.npz (true facet, relaxation overshoot > 1e-6)
this computes Z, groups the 92 into orbits, and for one representative per orbit:
  * checks exactly that Z has affine rank dim P(c) - 1 = 25 and that the affine function vanishing on Z
    is strictly one-signed on the other vertices (integer arithmetic, sympy);
  * returns that function as an integer vector (l, h) with l.v_pi <= h on all vertices, and its gap
    vector on the 120 permutations: (h - l.v_pi) divided by the gcd of its entries;
  * reads whether the facet follows linearly from PMSH(4) (README N3: ../logs/n5_inherited_dense40_rows.txt
    and ../logs/n5_inherited_nan15_rows.txt, rows of the notinherited15 and nan15 files).
Output: orbits.json.

usage: orbits.py
"""
import itertools
import json
import os
import sys

import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from facet_search import PM  # noqa: E402

n = 5
P = PM(n)
PERMS = [tuple(p) for p in P.PA]
PIDX = {p: i for i, p in enumerate(PERMS)}
C4 = np.array([4, 5, 6, 7, 8])
VI = np.array([[int(np.prod([C4[p[t]] for t in T])) for T in P.SUBS] for p in PERMS], dtype=object)


def group():
    """All 240 maps on permutation indices: pi -> r^e o pi o sigma."""
    out = []
    for e in (0, 1):
        for sigma in PERMS:
            out.append([PIDX[tuple((4 - p[sigma[i]]) if e else p[sigma[i]] for i in range(n))] for p in PERMS])
    return out


def exact_facet(Z):
    """Integer (l, h, gap) with l.v = h on Z, l.v <= h on all vertices and gap = h - l.v_pi primitive,
    or None if Z is not the vertex set of a facet. (l, h) is determined only up to the equations of aff P."""
    rankP = sp.Matrix([[*VI[i], 1] for i in range(len(PERMS))]).rank()  # = dim P + 1
    if sp.Matrix([[*VI[i], 1] for i in Z]).rank() != rankP - 1:
        return None, rankP - 1
    for v in sp.Matrix([[*VI[i], -1] for i in Z]).nullspace():  # every v vanishes on Z
        vals = [sum(VI[i][j] * v[j] for j in range(len(P.SUBS))) - v[-1] for i in range(len(PERMS))]
        if any(x != 0 for x in vals):
            break
    sgn = 1 if max(vals) > 0 else -1
    if any(sgn * x < 0 for x in vals):
        return None, rankP - 1
    den = sp.ilcm(*[sp.Rational(x).q for x in v])
    w = [-sgn * int(sp.Rational(x) * den) for x in v]  # now l.v - h <= 0
    l, h = w[:-1], w[-1]
    gap = [h - sum(VI[i][j] * l[j] for j in range(len(P.SUBS))) for i in range(len(PERMS))]
    k = int(sp.igcd(*gap))
    gap = [x // k for x in gap]
    assert min(gap) == 0 and [i for i in range(len(PERMS)) if gap[i] == 0] == sorted(Z)
    return (l, h, gap), rankP - 1


def main():
    pre = os.path.join(os.path.dirname(HERE), "hard_n5_seed41_c0_0.25_0.5_0.75_1")
    d = np.load(pre + ".npz")
    cs = d["c"] + 1.0
    V = P.verts(cs)
    g = V.mean(0)
    hard = [i for i in range(len(d["relaxation"])) if d["facet"][i] and d["relaxation"][i] > 1 + 1e-6]
    zset = {i: frozenset(np.flatnonzero(np.abs((V - g) @ d["normals"][i] - 1) < 1e-7).tolist()) for i in hard}
    # inheritance status (README N3): survivors of 40 shifts
    survive = set()
    for suffix, rows in (("_notinherited15", "n5_inherited_dense40_rows.txt"), ("_nan15", "n5_inherited_nan15_rows.txt")):
        src = np.load(pre + suffix + ".npz")["source_index"]
        for line in open(os.path.join(os.path.dirname(HERE), "logs", rows)):
            k, _, _, vinh = line.split()
            if float(vinh) > 1 + 1e-6:
                survive.add(int(src[int(k)]))
    G = group()
    orbit_of, orbits = {}, []
    for i in hard:
        if zset[i] in orbit_of:
            continue
        orb = {frozenset(gm[j] for j in zset[i]) for gm in G}
        for z in orb:
            orbit_of[z] = len(orbits)
        orbits.append(dict(orbit=orb, members=[]))
    for i in hard:
        orbits[orbit_of[zset[i]]]["members"].append(i)
    out = []
    for k, o in enumerate(orbits):
        rep = min(sorted(z) for z in o["orbit"])
        fac, dimP = exact_facet(rep)
        status = {"inherited" if i not in survive else "not inherited" for i in o["members"]}
        rel = sorted(round(float(d["relaxation"][i]), 6) for i in o["members"])
        l, h, gap = fac if fac else (None, None, None)
        out.append(dict(orbit=k, size=len(o["orbit"]), harvested=len(o["members"]), status=sorted(status),
                        relaxation=sorted(set(rel)), vertices=len(rep), exact_facet=fac is not None,
                        rep=rep, l=l, h=h, gap=gap))
        print(f"orbit {k}: size {len(o['orbit'])}, harvested {len(o['members'])}, |Z| = {len(rep)}, "
              f"exact facet {fac is not None}, PMSH(4) {sorted(status)}, relaxation {sorted(set(rel))}", flush=True)
    print(f"dim P = {dimP}; {len(hard)} hard facets in {len(orbits)} orbits; "
          f"not inherited: {sum(1 for i in hard if i in survive)}")
    json.dump(dict(spectrum=[4, 5, 6, 7, 8], subsets=[list(T) for T in P.SUBS], orbits=out),
              open(os.path.join(HERE, "orbits.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
