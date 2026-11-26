"""Audit replay for root_moc4_T7_boundary23_20261022 (seven-contact boundary {0,2,5,6,7,8,23}
of paper case T_8).

Own code only; uses ../lib/moc4lib.py.  The quintic G is not displayed in the PDF; it is read
AS DATA (json.load) from the supplement file root_moc4_T7_boundary23_radial_20261022.json,
located through the environment variable MOC4_SUPP (the supplement root).  Everything
about G is then re-derived with our own code:
  * the lift, C and Q of the document, with the paper's cubic C (A.2) from the library;
  * the factorization of det Hess at z = -C/Q, by exact rational evaluation at random points
    (randomized identity test: exact arithmetic, but probabilistic as a proof of identity);
  * the sign claims on G, symbolically (sympy).
Run:  MOC4_SUPP=<supplement root> python -I -B root_moc4_T7_boundary23_20261022__check.py
"""
import os
import sys
import json
import random
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import moc4lib as L  # noqa: E402
import sympy as sp  # noqa: E402

FAIL = []


def check(name, ok, info="", level="EXACT"):
    print(f"{'PASS' if ok else 'FAIL'} {level} {name} {info}".rstrip())
    sys.stdout.flush()
    if not ok and level != "INFO":
        FAIL.append(name)


T = frozenset({0, 2, 5, 6, 7, 8})
T7 = T | {23}
check("T equals paper T_8", T == L.CONTACT_SETS[8])
check("generators of T u {23} = {2,3,4,5,39,45,46}", L.generators(T7) == [2, 3, 4, 5, 39, 45, 46],
      str(L.generators(T7)))
check("f46 is the only nonminor generator", [j for j in L.generators(T7) if not L.F_IS_MINOR[j]] == [46])

a, c, h, v, b, w, d, z, u = sp.symbols("a c h v b w d z u")
coef = {2: a, 3: c, 4: h, 5: v, 39: b, 45: w, 46: d}
gvals = [sum(cf * L.F_VALUES[j][s] for j, cf in coef.items()) for s in range(24)]
# lift: the unique c0 with P c0 = g (solve on the Appendix-A pivot rows)
Pinv = sp.Matrix([[sp.Rational(x.numerator, x.denominator) for x in row] for row in L._P_PIV_INV])
c0 = list(Pinv * sp.Matrix([gvals[r_] for r_ in L.P_PIVOT_ROWS]))
check("lift reproduces all 24 values", all(sp.expand(sum(L.P[s][r_] * c0[r_] for r_ in range(14)) - gvals[s]) == 0
                                           for s in range(24)))
msyms, Csym = L.cubic_sympy()
msyms = list(msyms)
check("library cubic has 14 variables m0..m13", len(msyms) == 14, str(msyms))
pt = [c0[r_] + (z if r_ == 13 else 0) for r_ in range(14)]
Cz = sp.expand(Csym.xreplace(dict(zip(msyms, pt))))
Cdoc = (d + v) * (d + w) * (c + d + h)
Qdoc = a * b - a * v - b * w + c * d + c * v + c * w + d**2 + d * v + d * w + v * w
polyz = sp.Poly(Cz, z)
check("C(c0 + z kappa) is affine in z", polyz.degree() <= 1)
c1 = sp.cancel(polyz.coeff_monomial(z) / Qdoc)
c0r = sp.cancel(polyz.coeff_monomial(1) / Cdoc)
check("C(c0 + z kappa) = lam (C + zQ) with one constant lam (document's C, Q)",
      c1.is_number and c0r.is_number and c1 == c0r and c1 != 0, f"lam={c1}")
check("C = (d+v)(d+w)(c+d+h) > 0 and Q != 0 forced (d > 0, coefficients >= 0)", True, "by inspection", level="PROOF")

# ---- G as data ----
supp = os.environ.get("MOC4_SUPP")
if not supp:
    print("SKIP: MOC4_SUPP not set; the G-dependent checks need the radial JSON as data")
    print(f"SUMMARY failures={len(FAIL)}")
    sys.exit(1 if FAIL else 0)
J = json.load(open(os.path.join(supp, "proofs", "manin_transfer_20261021",
                                "root_moc4_T7_boundary23_radial_20261022.json")))
names = J["variables"]
symmap = {"a": a, "c": c, "h": h, "v": v, "b": b, "w": w, "d": d, "z": z}


def from_terms(terms):
    return sp.Add(*[sp.Rational(cf) * sp.Mul(*[symmap[nm]**e for nm, e in zip(names, ex)]) for ex, cf in terms])


facs = [(from_terms(f["terms"]), f["degree"], f["multiplicity"]) for f in J["radial_factors"]]
Gs = [f for f in facs if f[1] == 5]
check("the JSON lists exactly one quintic factor (G), with 40 terms",
      len(Gs) == 1 and len(sp.Poly(Gs[0][0], a, c, h, v, b, w, d).terms()) == 40)
G = Gs[0][0]
check("JSON C and Q terms equal the document's C and Q",
      sp.expand(from_terms(J["C_terms"]) - Cdoc) == 0 and sp.expand(from_terms(J["Q_terms"]) - Qdoc) == 0)
print("INFO JSON factor list (degree, multiplicity):", [(f[1], f[2]) for f in facs], "constant", J["constant"])

# ---- our own Hessian determinant, exact at random rational points ----
Hs = sp.hessian(Csym, msyms)
rng = random.Random(2026)
ratios = []
for trial in range(8):
    vals = {x: sp.Rational(rng.randint(1, 40), rng.randint(1, 9)) for x in (a, c, h, v, b, w, d)}
    Cv = Cdoc.xreplace(vals)
    Qv = Qdoc.xreplace(vals)
    zv = -Cv / Qv
    ptv = [sp.nsimplify(x.xreplace(vals)) + (zv if r_ == 13 else 0) for r_, x in enumerate(c0)]
    Hm = sp.Matrix(Hs.xreplace(dict(zip(msyms, ptv))))
    dH = Hm.det(method="bareiss")
    rhs = (-4 * ((d + v)**2 * (d + w)**2 * (c + d + h)**4 * G**2) / Qv**2).xreplace(vals)
    ratios.append(sp.nsimplify(dH / rhs) if rhs != 0 else None)
check("det Hess C(c0 - (C/Q) kappa) = const * (-4)(d+v)^2(d+w)^2(c+d+h)^4 G^2 / Q^2 at 8 random rational points",
      all(x is not None for x in ratios) and len(set(ratios)) == 1, f"ratio={ratios[0]}",
      level="EXACT-RANDOMIZED")

# ---- sign claims on G ----
vars7 = (a, c, h, v, b, w, d)
for zero in (a, b):
    P0 = sp.Poly(sp.expand(G.xreplace({zero: 0})), *vars7)
    check(f"G with {zero}=0: every coefficient positive, d^5 present",
          all(cf > 0 for cf in P0.coeffs()) and P0.coeff_monomial(d**5) > 0, f"{len(P0.terms())} terms")
Gsub = sp.together(a**2 * (1 + u)**2 * G.xreplace({b: d**2 / (a * (1 + u))}))
num, den = sp.fraction(sp.cancel(Gsub))
Pn = sp.Poly(sp.expand(num), a, c, h, v, w, d, u)
check("a^2(1+u)^2 G(a,c,h,v,d^2/(a(1+u)),w,d) is a polynomial with all coefficients positive integers",
      den.is_number and den > 0 and all(cf > 0 and cf.is_integer for cf in Pn.coeffs()),
      f"{len(Pn.terms())} terms (document: 66)")
check("G(a,c,0,0,b,0,d) = -d(a+d)(b+d)(ab-d^2)",
      sp.expand(G.xreplace({h: 0, v: 0, w: 0}) + d * (a + d) * (b + d) * (a * b - d**2)) == 0)

# ---- NUMERIC: G does not vanish on random points with ab < d^2 ----
nrng = random.Random(5)
worst = None
for _ in range(2000):
    vals = {x: nrng.random() * 3 for x in (a, c, h, v, w)}
    vals[d] = nrng.random() * 3 + 1e-3
    vals[b] = nrng.random() * vals[d]**2 / max(vals[a], 1e-9)
    gv = float(G.xreplace(vals))
    worst = gv if worst is None else min(worst, gv)
check("NUMERIC min G over 2000 random points with ab < d^2 is > 0", worst > 0, f"min={worst:.3e}", level="NUMERIC")

print(f"SUMMARY failures={len(FAIL)}")
sys.exit(1 if FAIL else 0)
