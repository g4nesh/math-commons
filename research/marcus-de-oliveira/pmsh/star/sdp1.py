"""First-level SDP relaxation over the orthogonal group O(6):
   minimise <Q, Y> + c0   s.t.  Y PSD (36x36),  partial traces = I_6   (Y stands for vec(R) vec(R)^T).
A nonnegative optimum certifies  w.x(U) + w0 >= 0  for every U in U(4)  (dual: Q - I(x)Lc - Lr(x)I PSD)."""
import sys, numpy as np, cvxpy as cp
from base import *
from so6 import quad_of
def relax_min(Q, c0, solver='CLARABEL', verbose=False):
    Y = cp.Variable((36, 36), symmetric=True)
    cons = [Y >> 0]
    idx = lambda a, b: 6 * a + b
    for b in range(6):
        for d in range(b, 6):
            cons.append(sum(Y[idx(a, b), idx(a, d)] for a in range(6)) == (1.0 if b == d else 0.0))
    for a in range(6):
        for c in range(a, 6):
            cons.append(sum(Y[idx(a, b), idx(c, b)] for b in range(6)) == (1.0 if a == c else 0.0))
    prob = cp.Problem(cp.Minimize(cp.trace(Q @ Y) + c0), cons)
    prob.solve(solver=solver, verbose=verbose)
    return prob.value, Y.value
def vec_w(terms):
    w = np.zeros(NX)
    for coef, kind, a, b in terms:
        w[iB(a, b) if kind == 'B' else iD(a, b)] += coef
    return w
G0 = [(1, 'B', 0, 0), (1, 'B', 2, 2), (1, 'D', (0, 3), (2, 3)), (-1, 'D', (0, 1), (0, 3))]
if __name__ == '__main__':
    Q, c0 = quad_of(vec_w(G0)); v, _ = relax_min(Q, c0)
    print("relaxation lower bound for g0: %.6f   (true minimum over U(4): -0.067442)" % v)
    for t in [0.05, 0.2, 0.5, 1.0, 2.0, 5.0, 20.0]:
        for b in [((0, 2), (0, 3)), ((0, 2), (1, 2))]:
            w = vec_w([(t * c, k, a_, b_) for c, k, a_, b_ in G0] + [(t + 1, 'D', (0, 1), (0, 1)), (t * t, 'D', b[0], b[1])])
            Q, c0 = quad_of(w); v, _ = relax_min(Q, c0)
            print("t = %-5g  second minor D[%s,%s]: relaxation min of t*g0 + (t+1)*D[12,12] + t^2*D = %+.3e" % (t, "".join(str(i+1) for i in b[0]), "".join(str(i+1) for i in b[1]), v), flush=True)
