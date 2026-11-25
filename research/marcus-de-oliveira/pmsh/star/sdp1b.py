"""Level-1 relaxation over O(6) strengthened by the valid linear constraints B[s,j] >= 0."""
import numpy as np, cvxpy as cp
from base import *
from so6 import quad_of, quad_B
from sdp1 import vec_w, G0
def relax_min(Q, c0):
    Y = cp.Variable((36, 36), symmetric=True); cons = [Y >> 0]; idx = lambda a, b: 6 * a + b
    for b in range(6):
        for d in range(b, 6): cons.append(sum(Y[idx(a, b), idx(a, d)] for a in range(6)) == float(b == d))
    for a in range(6):
        for c in range(a, 6): cons.append(sum(Y[idx(a, b), idx(c, b)] for b in range(6)) == float(a == c))
    for s in range(4):
        for j in range(4):
            A, c = quad_B(s, j); cons.append(cp.trace(A @ Y) + c >= 0)
    prob = cp.Problem(cp.Minimize(cp.trace(Q @ Y) + c0), cons); prob.solve(solver='CLARABEL'); return prob.value
Q, c0 = quad_of(vec_w(G0)); print("with B >= 0: relaxation bound for g0: %.6f (true -0.067442)" % relax_min(Q, c0))
for t in [0.2, 1.0, 5.0]:
    w = vec_w([(t * c, k, a_, b_) for c, k, a_, b_ in G0] + [(t + 1, 'D', (0, 1), (0, 1)), (t * t, 'D', (0, 2), (0, 3))])
    Q, c0 = quad_of(w); print("t=%g: relaxation min of the hard form: %+.4e" % (t, relax_min(Q, c0)))
