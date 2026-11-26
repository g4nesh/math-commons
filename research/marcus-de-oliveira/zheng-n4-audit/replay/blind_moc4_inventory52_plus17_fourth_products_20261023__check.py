"""Own replay for blind_moc4_inventory52_plus17_fourth_products_20261023: chart coverage only.
The coefficient tables of this note are NOT replayed (see report)."""
import sys
import sympy as sp
fails = 0
def check(name, ok):
    global fails
    print(("PASS " if ok else "FAIL ") + name); fails += (not ok)
s, Ws, u, v, a, b, c = sp.symbols('s W_s u v a b c', positive=True)
a4 = s / (s + Ws); b4 = (s + Ws) * (1 + u) / s
check("(4): a b = 1 + u", sp.simplify(a4 * b4 - (1 + u)) == 0)
check("(4): 1 - a = Ws/(s+Ws) >= 0", sp.simplify(1 - a4 - Ws / (s + Ws)) == 0)
m = (b + 1) * (c + 1)
check("m - abc = (1-a)bc + b + c + 1", sp.expand(m - a * b * c - ((1 - a) * b * c + b + c + 1)) == 0)
na = 1 + 2 * s; da = 1 + s; nb = na + 2 * da * u; db = da * (1 + u)
A6 = na / da; B6 = nb / db
check("(6) row 2/3: a - 1 = s/(1+s), 2 - a = 1/(1+s)", sp.simplify(A6 - 1 - s / (1 + s)) == 0 and sp.simplify(2 - A6 - 1 / (1 + s)) == 0)
check("(6) row 3: b - a = u(2-a)/(1+u), 2 - b = (2-a)/(1+u)", sp.simplify(B6 - A6 - u * (2 - A6) / (1 + u)) == 0 and sp.simplify(2 - B6 - (2 - A6) / (1 + u)) == 0)
C7 = (nb + 3 * db * v) / (db * (1 + v))
check("(7): c - b = v(3-b)/(1+v), 3 - c = (3-b)/(1+v)", sp.simplify(C7 - B6 - v * (3 - B6) / (1 + v)) == 0 and sp.simplify(3 - C7 - (3 - B6) / (1 + v)) == 0)
print("UNCHECKED: all L_ij, controller, three-capacity and box coefficient tables of the note")
print("RESULT", "all EXACT checks passed" if fails == 0 else "%d failures" % fails)
sys.exit(1 if fails else 0)
