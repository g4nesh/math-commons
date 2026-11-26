"""Symbolic closed form, in the ratio rho, of every hard facet attached to the base exotic facet
g0 = B11 + B33 + D[14,34] - D[12,14].   q = (rho, -1, 1-rho) sweeps all directions as rho ranges over R."""
import itertools, numpy as np, sympy as sp
from base import *
rho = sp.symbols('rho', real=True)
q = [rho, -1, 1 - rho]
g0 = np.zeros(NX, dtype=int)
g0[iB(0, 0)] += 1; g0[iB(2, 2)] += 1; g0[iD((0, 3), (2, 3))] += 1; g0[iD((0, 1), (0, 3))] -= 1
Z0 = set(np.where(VX @ g0 == 0)[0]); assert len(Z0) == 13
def shift_sym(k):
    eta = [0, 0, 0]; eta[k] = 1; eta[2] -= 1
    dx = [sp.Integer(0)] * NX
    for a, (S, Sc) in enumerate(MATCH):
        for b, (T, Tc) in enumerate(MATCH):
            for S_ in (S, Sc):
                for T_ in (T, Tc): dx[iD(S_, T_)] += sp.Rational(eta[a], 2) * q[b]
    return dx
sh = [shift_sym(0), shift_sym(1)]
gs = [sum(int(g0[k]) * sh[i][k] for k in range(NX)) for i in range(2)]
def name(k):
    k -= 16; S, T = P2[k // 6], P2[k % 6]
    return "D[%s,%s]" % ("".join(str(s + 1) for s in S), "".join(str(t + 1) for t in T))
def canon(k):
    """Jacobi: D[S,T] = D[S^c,T^c]; keep the representative whose row set contains row 1"""
    k0 = k - 16; S, T = P2[k0 // 6], P2[k0 % 6]
    return k if 0 in S else iD(comp(S), comp(T))
rows = {}
for a, b in itertools.combinations(range(16, NX), 2):
    if canon(a) != a or canon(b) != b: continue
    Z = [v for v in Z0 if VX[v, a] == 0 and VX[v, b] == 0]
    if len(Z) != 11: continue
    A = sp.Matrix([[sh[0][a], sh[0][b]], [sh[1][a], sh[1][b]]])
    if sp.simplify(A.det()) == 0: continue
    al = A.solve(-sp.Matrix(gs)); al = [sp.factor(sp.simplify(v)) for v in al]
    # ranges of rho where both coefficients are positive
    rng = sp.solve_univariate_inequality(al[0] > 0, rho, relational=False).intersect(sp.solve_univariate_inequality(al[1] > 0, rho, relational=False))
    if rng.is_empty: continue
    rows[(a, b)] = (al, rng, frozenset(Z))
print("hard-facet families attached to g0 (alpha, alpha' > 0 on the stated range of rho):")
for (a, b), (al, rng, Z) in sorted(rows.items(), key=lambda kv: str(kv[1][1])):
    print("  rho in %-22s  g0 + (%s)*%s + (%s)*%s >= 0" % (rng, al[0], name(a), al[1], name(b)))
import pickle; pickle.dump({k: (str(v[0][0]), str(v[0][1]), str(v[1])) for k, v in rows.items()}, open('families_g0.pkl', 'wb'))
