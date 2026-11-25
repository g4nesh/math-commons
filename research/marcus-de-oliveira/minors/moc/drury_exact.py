"""Exact form of the deepest exotic-facet violation on the skew-conference family O(y) = y I + x S."""
import sympy as sp
y = sp.symbols('y', real=True)
x = sp.sqrt((1 - y**2) / 3)
S = sp.Matrix([[0, 1, 1, -1], [-1, 0, -1, -1], [-1, 1, 0, 1], [1, 1, -1, 0]])
assert S + S.T == sp.zeros(4) and S * S == -3 * sp.eye(4)
O = y * sp.eye(4) + x * S
assert sp.simplify(O.T * O - sp.eye(4)) == sp.zeros(4)
m = lambda I, J: O.extract(I, J).det()**2
g = sp.expand(sp.simplify(O[0, 0]**2 + O[2, 2]**2 + m([0, 3], [2, 3]) - m([0, 1], [0, 3])))
print("g(y) on the family =", sp.factor(g))
crit = sp.solve(sp.diff(g, y), y)
vals = [(c, sp.nsimplify(c), sp.N(g.subs(y, c), 15)) for c in crit if c.is_real and -1 <= c <= 1]
for c, cs, v in vals: print("  critical y =", sp.simplify(c), "~", sp.N(c, 12), "  g =", sp.simplify(g.subs(y, c)), "~", v)
best = min(vals, key=lambda t: t[2])
print("minimum: y* =", sp.simplify(best[0]), " g_min =", sp.simplify(g.subs(y, best[0])), "=", sp.N(g.subs(y, best[0]), 15))
print("Drury functional 1 - g_min =", sp.simplify(1 - g.subs(y, best[0])), "=", sp.N(1 - g.subs(y, best[0]), 15))
print("check Zheng y=1/7: g =", sp.nsimplify(g.subs(y, sp.Rational(1, 7))))

# Minimal polynomials: eliminate y via u = y^2.  On the family, (2u - G)^2 = (16/27) u (1-u)^3, and y*^2 is a root of
# 16u^3 - 24u^2 + 36u - 1 (the critical-point equation).  Then substitute G = 1 - D and apply Cardano.
u, G, Dv = sp.symbols('u G D')
cubic_u = 16*u**3 - 24*u**2 + 36*u - 1
print("y*^2 is the unique real root of", cubic_u, ":", [sp.N(r, 15) for r in sp.Poly(cubic_u, u).nroots() if abs(sp.im(r)) < 1e-12])
res = sp.factor(sp.resultant(sp.expand((2*u - G)**2 - sp.Rational(16, 27) * u * (1 - u)**3), cubic_u, u))
print("resultant in G:", res)
gpoly = 4*G**3 - 12*G**2 + 14*G + 1
groot = [r for r in sp.Poly(gpoly, G).nroots(n=20) if abs(sp.im(r)) < 1e-15][0]
print("g_min = real root of", gpoly, "=", sp.N(groot, 16), "; matches family minimum:", abs(groot - sp.N(g.subs(y, best[0]), 20)) < 1e-14)
dpoly = sp.expand(gpoly.subs(G, 1 - Dv))
print("with G = 1 - D:", sp.factor(dpoly), " (i.e. 4D^3 + 2D - 7 = 0)")
cardano = sp.cbrt(sp.Rational(7, 8) + 11 * sp.sqrt(33) / 72) + sp.real_root(sp.Rational(7, 8) - 11 * sp.sqrt(33) / 72, 3)
print("Cardano form:", cardano, "=", sp.N(cardano, 16), "; residual", sp.N(4 * cardano**3 + 2 * cardano - 7, 20))
