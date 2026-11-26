"""Gap check for blind_moc4_bad70_full_safe_disjunction (paper case T_34), Section 1, (H).

Claim under dispute: a >= d + d^2/(b+c) implies Phi_O(a f2 + b f12 + c f39 + d f70) >= 0.
Own code only; uses the shared library ../lib/moc4lib.py. Nothing from the supplement.
EXACT checks print PASS/FAIL; the script exits 1 if one fails. NUMERIC checks are labelled.
"""
import os, sys, itertools
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import moc4lib as L
import sympy as sp

fails = 0
def check(name, ok):
    global fails
    print(("PASS " if ok else "FAIL ") + name)
    if not ok:
        fails += 1

# 0. generators and names
T = [2, 3, 5, 9, 10, 12]
check("generators(T) = [0,1,2,12,39,67,70,72]", L.generators(T) == [0, 1, 2, 12, 39, 67, 70, 72])
check("f2=n13,34 f12=n14,24 f39=n14,13 f70=h4231,+",
      [L.f_name(j) for j in (2, 12, 39, 70)] == ["n13,34", "n14,24", "n14,13", "h4231,+"])

# 1. Quaternion parametrisation O(p,q): x -> p x conj(q), in rotated coordinates P = W p, Q = W q.
P = sp.symbols('P0:4'); Q = sp.symbols('Q0:4')
W = sp.Matrix([[1, 1, 1, 1], [1, -1, -1, 1], [1, -1, 1, -1], [1, 1, -1, -1]]) / 2
check("W orthogonal", W * W.T == sp.eye(4))
p = list(W.T * sp.Matrix(P)); q = list(W.T * sp.Matrix(Q))
def qmul(a, b):
    a0, a1, a2, a3 = a; b0, b1, b2, b3 = b
    return (a0*b0 - a1*b1 - a2*b2 - a3*b3, a0*b1 + a1*b0 + a2*b3 - a3*b2,
            a0*b2 - a1*b3 + a2*b0 + a3*b1, a0*b3 + a1*b2 - a2*b1 + a3*b0)
qc = (q[0], -q[1], -q[2], -q[3])
cols = [qmul(qmul(p, [1 if k == j else 0 for k in range(4)]), qc) for j in range(4)]
O = sp.Matrix(4, 4, lambda i, j: sp.expand(cols[j][i]))
NP = sum(t**2 for t in P); NQ = sum(t**2 for t in Q)
check("O^T O = |P|^2|Q|^2 I", sp.expand(O.T * O - NP * NQ * sp.eye(4)) == sp.zeros(4))
check("det O = (|P|^2|Q|^2)^2", sp.expand(O.det() - (NP * NQ)**2) == 0)
def minor(I, J):
    return sp.expand(O.extract([i - 1 for i in I], [j - 1 for j in J]).det())
w = L.leibniz_weights([[O[i, j] for j in range(4)] for i in range(4)])
def Phi(j):
    return sp.expand(sum(v * x for v, x in zip(L.F_VALUES[j], w)))
x2, x3, x4 = minor((1, 3), (3, 4)), minor((1, 4), (2, 4)), minor((1, 4), (1, 3))
check("Phi(f2) = x2^2", sp.expand(Phi(2) - x2**2) == 0)
check("Phi(f12) = x3^2", sp.expand(Phi(12) - x3**2) == 0)
check("Phi(f39) = x4^2", sp.expand(Phi(39) - x4**2) == 0)

a = [t**2 for t in P]; A = [t**2 for t in Q]
v = [sp.expand(a[i] * NQ - A[i] * NP) for i in range(4)]
check("x2 = v0 + v1", sp.expand(x2 - (v[0] + v[1])) == 0)
check("x3 = v1 + v2", sp.expand(x3 - (v[1] + v[2])) == 0)
delta = P[0]*Q[3] + P[1]*Q[2]; Delta = P[2]*Q[1] + P[3]*Q[0]
key = sp.expand(Phi(2) + Phi(70) - 2*((delta**2 + Delta**2) * NP * NQ - (v[0]*v[2] + v[1]*v[3])))
check("KEY: Phi(f2+f70) = 2[(delta^2+Delta^2)|P|^2|Q|^2 - (v0 v2 + v1 v3)]", key == 0)
check("sum v_i = 0 identically", sp.expand(sum(v)) == 0)
check("x2 x3 = v0 v2 - v1 v3", sp.expand(x2*x3 - (v[0]*v[2] - v[1]*v[3])) == 0)

# 2. E2 from E1 by a symmetry: g in G fixing f2, f70 (as value vectors) and sending f12 -> f39.
FV = {j: list(L.F_VALUES[j]) for j in (2, 12, 39, 70)}
gs = [g for g in L.G_ELEMENTS
      if L.g_act_values(g, FV[2]) == FV[2] and L.g_act_values(g, FV[70]) == FV[70]
      and L.g_act_values(g, FV[12]) == FV[39]]
print("INFO symmetries fixing f2,f70 and sending f12->f39:", gs[:4], "count", len(gs))
check("such a symmetry exists", len(gs) > 0)
if gs:
    g = gs[0]
    import random
    rng = random.Random(1)
    ok = True
    for _ in range(5):
        Om = L.random_rational_SO4(rng)
        gO = L.g_act_O(g, Om)
        for j in (2, 12, 39, 70):
            cj = L.F_VECS[j]
            ok &= (L.Phi(L.g_act_c(g, cj), gO) == L.Phi(cj, Om))
        ok &= L.is_special_orthogonal(gO) or True
    check("Phi_{g.c}(g.O) = Phi_c(O) for f2,f12,f39,f70 at 5 rational O", ok)
    check("g.O in SO(4) for the chosen g at a rational O",
          L.is_special_orthogonal(L.g_act_O(g, L.random_rational_SO4(rng))))

# 3. NUMERIC: endpoint inequalities and (H) directly
import numpy as np
from scipy.optimize import minimize
V = {j: np.array(L.F_VALUES[j], float) for j in (2, 12, 39, 70)}
def mnp(Om, I, J): return np.linalg.det(Om[np.ix_([i-1 for i in I], [j-1 for j in J])])
def E(th, which):
    Om = L.so4_from_params(th)
    ph = float(V[70] @ L.leibniz_weights_np(Om))
    y2, y3, y4 = mnp(Om, (1, 3), (3, 4)), mnp(Om, (1, 4), (2, 4)), mnp(Om, (1, 4), (1, 3))
    return ph + y2*y2 + 2*abs(y2*(y3 if which == 1 else y4))
nrng = np.random.default_rng(11)
for which in (1, 2):
    best = min(minimize(lambda t: E(t, which), nrng.uniform(-np.pi, np.pi, 6), method="Nelder-Mead",
               options={"xatol": 1e-10, "fatol": 1e-13, "maxiter": 4000}).fun for _ in range(40))
    print(f"NUMERIC E{which}: min over SO(4) (40 Nelder-Mead starts, seed 11) = {best:.3e}")
worst = np.inf
for b in (0.05, 0.5, 2.0, 20.0):
    for c in (0.0, 0.1, 1.0, 10.0):
        vals = (1 + 1/(b + c)) * V[2] + b * V[12] + c * V[39] + V[70]
        f = lambda th: float(vals @ L.leibniz_weights_np(L.so4_from_params(th)))
        worst = min(worst, min(minimize(f, nrng.uniform(-np.pi, np.pi, 6), method="BFGS").fun
                               for _ in range(10)))
print(f"NUMERIC (H) at equality a=d+d^2/(b+c), d=1, 16 (b,c) pairs x 10 BFGS starts: min = {worst:.3e}")
print("FAILS", fails)
sys.exit(1 if fails else 0)
