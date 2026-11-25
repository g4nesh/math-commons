"""Exact certificates showing PMSH is strictly weaker than compound Birkhoff (n=4).
H = O* D_c O*^T with Zheng's rational orthogonal O* (which violates an exotic compound facet exactly: psi = -158/2401).
For several rational spectra c we find lambda >= 0 on S_4, sum 1, with sum_pi lambda_pi c^{pi(T)} = det H[T,T] for ALL T,
verified in exact rational arithmetic (fractions.Fraction)."""
import itertools, json
from fractions import Fraction as Fr
import numpy as np
from scipy.optimize import linprog
import sympy as sp
O = [[Fr(v, 7) for v in row] for row in [[1, 4, 4, 4], [-4, 1, 4, -4], [-4, -4, 1, 4], [-4, 4, -4, 1]]]
n = 4
def matmul(A, B): return [[sum(A[i][k] * B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]
def det(M):
    M = [row[:] for row in M]; k = len(M); d = Fr(1)
    for i in range(k):
        p = next((r for r in range(i, k) if M[r][i] != 0), None)
        if p is None: return Fr(0)
        if p != i: M[i], M[p] = M[p], M[i]; d = -d
        d *= M[i][i]
        for r in range(i + 1, k):
            f = M[r][i] / M[i][i]
            for cc in range(i, k): M[r][cc] -= f * M[i][cc]
    return d
OT = [list(r) for r in zip(*O)]
assert matmul(OT, O) == [[Fr(int(i == j)) for j in range(n)] for i in range(n)]
Ts = [T for k in range(1, n) for T in itertools.combinations(range(n), k)]
perms = list(itertools.permutations(range(n)))
results = []
for c in [[Fr(1), Fr(2), Fr(3), Fr(5)], [Fr(1, 3), Fr(1), Fr(7, 2), Fr(4)], [Fr(-2), Fr(1, 2), Fr(1), Fr(9)], [Fr(1), Fr(11, 10), Fr(6, 5), Fr(10)]]:
    H = matmul(matmul(O, [[c[i] if i == j else Fr(0) for j in range(n)] for i in range(n)]), OT)
    m = [det([[H[i][j] for j in T] for i in T]) for T in Ts]
    V = [[np.prod([c[p[t]] for t in T]) for T in Ts] for p in perms]          # vertex vectors (exact Fractions)
    Vf = np.array([[float(x) for x in row] for row in V]); mf = np.array([float(x) for x in m])
    r = linprog(np.zeros(24), A_eq=np.vstack([Vf.T, np.ones(24)]), b_eq=np.r_[mf, 1.0], bounds=[(0, None)] * 24, method='highs-ds')
    assert r.status == 0, "float LP infeasible"
    supp = [i for i in range(24) if r.x[i] > 1e-9]
    # exact solve on the support (least-norm exact solution of the consistent system)
    A = sp.Matrix([[V[i][t] for i in supp] for t in range(len(Ts))] + [[1] * len(supp)])
    b = sp.Matrix([m[t] for t in range(len(Ts))] + [1])
    sol, params = A.gauss_jordan_solve(b)
    sol = sol.subs({p: 0 for p in params})
    lam = [Fr(int(sp.fraction(x)[0]), int(sp.fraction(x)[1])) for x in sol]
    full = [Fr(0)] * 24
    for i, v in zip(supp, lam): full[i] = v
    ok_eq = all(sum(full[i] * V[i][t] for i in range(24)) == m[t] for t in range(len(Ts))) and sum(full) == 1
    ok_nn = all(v >= 0 for v in full)
    results.append({'c': [str(x) for x in c], 'support': [''.join(map(str, perms[i])) for i in supp],
                    'lambda': {''.join(map(str, perms[i])): str(full[i]) for i in supp}, 'exact_equalities': ok_eq, 'nonnegative': ok_nn})
    print("c =", [str(x) for x in c], "| exact equalities:", ok_eq, "| nonnegative:", ok_nn, "| support size", len(supp))
json.dump({'O_star_times_7': [[1, 4, 4, 4], [-4, 1, 4, -4], [-4, -4, 1, 4], [-4, 4, -4, 1]], 'certificates': results}, open('moc/runs/pmsh_exact_zheng.json', 'w'), indent=1)
