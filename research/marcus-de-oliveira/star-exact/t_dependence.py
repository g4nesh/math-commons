"""How do the amplitude spaces depend on t?  For each t, list the distinct entries of the rationalised extra-kernel rows
and compare with t/(1+t).  Then test the hypothesis that the amplitude space of every block at t is spanned by
   (rows common to all t)  and  rows of the form  (1+t)*p - t*p'   with p, p' complement polynomials."""
import sys, json, pickle, numpy as np, sympy as sp
from fractions import Fraction as Fr
from exact_space import BLOCK_LIST
from reduced_sdp import load_comp
from second_reduction import rref_cp
comp = load_comp(1.0); ts = [Fr(x) for x in sys.argv[1:]]; kern = {}
for tq in ts:
    S = pickle.load(open('center_t%g.pkl' % float(tq), 'rb'))['S']; vals = {}
    for bi, (c, Sb) in enumerate(zip(comp, S)):
        if Sb is None: continue
        w, V = np.linalg.eigh(Sb); k = int((w < 1e-6).sum())
        if k == 0: continue
        R = rref_cp(V[:, :k].T)
        Rq = tuple(tuple(Fr(float(x)).limit_denominator(12) if abs(x) > 2e-3 else Fr(0) for x in row) for row in R)
        kern[(tq, bi)] = Rq
        for row in Rq:
            for q in row: vals[q] = vals.get(q, 0) + 1
    print("t = %-4s distinct entries of the extra-kernel rows: %s ; t/(1+t) = %s" % (tq, {str(k): v for k, v in sorted(vals.items())}, tq / (1 + tq)))
# hypothesis: replacing the entry t/(1+t) by a symbol gives the same pattern for every t
pat = {}
for (tq, bi), Rq in kern.items():
    r = tq / (1 + tq); key = tuple(tuple('r' if (q == r and q not in (0, 1)) else str(q) for q in row) for row in Rq)
    pat.setdefault(bi, set()).add(key)
same = sum(len(v) == 1 for v in pat.values())
print("blocks with extra kernel: %d ; blocks whose pattern is identical for all t after substituting r = t/(1+t): %d" % (len(pat), same))
shapes = {}
for v in pat.values():
    for key in v: shapes[key] = shapes.get(key, 0) + 1
for key, cnt in sorted(shapes.items(), key=lambda kv: -kv[1]): print("   pattern %s : %d" % (key, cnt))
