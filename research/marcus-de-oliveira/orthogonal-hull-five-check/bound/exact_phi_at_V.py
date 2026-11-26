"""Exact re-check (sympy) that the matrix V of task mdo-orthogonal-hull-five
(found by the agent ganesh-claude-minors; typed in here from its README) is
unitary and that Phi(V) = 37/25.  Written independently of that task's code.

Usage: python exact_phi_at_V.py
"""
import sympy as sp

s, t = sp.sqrt(6), sp.I * sp.sqrt(15)
V = sp.Matrix([[4, 4*s, 4*s, 4*s, 4*s],
               [4*s, -1-t, -6-2*t, -6+2*t, 9+t],
               [4*s, -6-2*t, -1+t, 9-t, -6+2*t],
               [4*s, -6+2*t, 9-t, -1+t, -6-2*t],
               [4*s, 9+t, -6+2*t, -6-2*t, -1-t]]) / 20

D = (V * V.H - sp.eye(5)).applyfunc(sp.expand)
print("V V^* = I exactly (every entry of V V^* - I expands to 0):", all(x == 0 for x in D))

def sqabs_minor(I, J):
    m = sp.expand(V.extract(sorted(I), sorted(J)).det())
    return sp.expand(m * sp.conjugate(m))

phi = 0
for v in range(5):
    phi += -sqabs_minor([v], [v])
    S1, S2 = [(v + 1) % 5, (v - 1) % 5], [(v + 2) % 5, (v - 2) % 5]
    phi += sqabs_minor(S1, S2) + sqabs_minor(S2, S1)
    for d in (1, 2):
        phi -= sqabs_minor([v, (v + d) % 5], [v, (v - d) % 5]) + sqabs_minor([v, (v - d) % 5], [v, (v + d) % 5])
phi = sp.expand(phi)
print("Phi(V) =", phi, "; Phi(V) - 37/25 expands to 0:", sp.expand(phi - sp.Rational(37, 25)) == 0)
