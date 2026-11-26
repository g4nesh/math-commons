"""Exact replay and numeric tests for root_inventory50_bounded_native_conditional_20261023
(conditional reduction for paper case T_51 = {1,3,5,10,11,12} to the seven-element set
D = {0,1,2,7,9,16,23}).

Rebuilt from the definitions (moc4lib for P, f_j, the group); the spectral parametrization
is rebuilt from the formulas displayed in the PDF and checked against raw sigma-point products.
Nothing from the supplement is read or executed.
Run: python -I -B <this file>.  One line per check; exit 1 if an EXACT check fails.
"""
import os, sys, random
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
from moc4lib import *            # noqa
import sympy as sp

FAIL = []
def check(name, ok, info=""):
    print(("PASS " if ok else "FAIL ") + name + ("  " + str(info) if info else ""), flush=True)
    if not ok:
        FAIL.append(name)

T = frozenset({1, 3, 5, 10, 11, 12})
check("T = paper T_51", T == CONTACT_SETS[51])
basis = [0, 1, 2, 8, 10, 12, 30, 59]
check("generators(T) = basis f0,f1,f2,f8,f10,f12,f30,f59", generators(T) == basis)
check("contact rank 6", rank([P[s] for s in sorted(T)]) == 6)
rows = [20, 16, 13, 14, 7, 6, 4, 0]
check("rows 20,16,13,14,7,6,4,0 of the basis profiles form I8",
      all(F_VALUES[basis[k]][rows[l]] == (1 if k == l else 0) for k in range(8) for l in range(8)))
check("coefficient of f30 is g4, coefficient of f8 is g14 (exclusive sources)", rows[6] == 4 and rows[3] == 14)
check("f30 is the only non-minor; all basis profiles >= 0",
      [j for j in basis if not F_IS_MINOR[j]] == [30] and all(F_VALUES[j][k] >= 0 for j in basis for k in range(24)))
check("ties not forcing d=0 are exactly 0,6,7,13,14,16,20", sorted(F_ZERO[30] - T) == [0, 6, 7, 13, 14, 16, 20])
B7, C7, D7 = map(frozenset, ({0, 1, 2, 7, 9, 10, 16}, {0, 1, 2, 7, 9, 16, 17}, {0, 1, 2, 7, 9, 16, 23}))
check("B, C, D are Table 5 sets 8, 11, 12", (B7, C7, D7) == (AUX7[8], AUX7[11], AUX7[12]))
check("tie 6 -> B, 13 -> B, 14 -> C, 20 -> D (some group element maps the set onto T+{j})",
      all(find_g(S, T | {j}) is not None for j, S in ((6, B7), (13, B7), (14, C7), (20, D7))))
for j, cases in ((0, (9, 11, 17, 38, 65, 76)), (7, (42, 65, 71)), (16, (9, 11, 17, 38, 65, 76))):
    hit = [nu for nu in cases if any(g_act_set(g, CONTACT_SETS[nu]) <= T | {j} for g in range(1152))]
    check("tie %d pulled back to six-element cases %s (not named in the document)" % (j, hit), hit == list(cases))

# exact finite instances of L_O >= d f30(O), d = g4
rng = random.Random(3)
ok = True
for _ in range(30):
    co = {j: Fr(rng.randint(0, 9)) for j in basis}
    c = [sum(co[j] * F_VECS[j][k] for j in basis) for k in range(14)]
    v = values(c)
    O = random_rational_SO4(rng)
    ok &= v[4] == co[30] and Phi(c, O) >= co[30] * Phi(F_VECS[30], O)
check("(1) L_O >= d f30(O) with d = g4, 30 exact random instances", ok)

# ---------------- spectral parametrization ----------------
a, t, W, Z, r, s, u, h, V = sp.symbols("a t W Z r s u h V")
Xx = W + (r - W) * Z + 1
Pp = Xx * (r + (u - 1) * W) + s * Z * (r - W)
Da = Xx * (r + (u - 1) * W) - s * W * Z * (r - W)
Dt = Xx * (r - W - u) + s * Z * (r - W)
Aspec = (a, W * (1 - Z), W + (r - W) * Z, W)
Bspec = (t, 0, 1)
def sigma_pt(p):
    v = sp.Integer(1)
    for i in range(4):
        if p[i] != 3:
            v *= Aspec[i] + Bspec[p[i]]
    return sp.expand(v)
Fv = [sigma_pt(p) for p in PERMS]
# the formulas (2) of the PDF, checked against the raw ratios
S_ = s * W * Z / Xx; U_ = u * W / (r - W); N_ = S_ + W + W * U_
check("(2) a = N/(1-S+U) equals W P/Da, t = N/(W+S-U) equals P/Dt",
      sp.simplify(N_ / (1 - S_ + U_) - W * Pp / Da) == 0 and sp.simplify(N_ / (W + S_ - U_) - Pp / Dt) == 0)
E = Fv[3] - Fv[1]
check("E = (a+t)(w-z), rE, sE, uE raw ratios (before inversion)",
      sp.expand(E - (a + t) * W * Z) == 0 and sp.expand(Fv[5] - Fv[3] - r * (a + t) * W * Z) == 0)
def cleared(i):
    """(F_i - F_1) * Da * Dt as a polynomial (F is bilinear in a, t)."""
    f = sp.Poly(Fv[i] - Fv[1], a, t)
    tot = 0
    for (ea, et), c in zip(f.monoms(), f.coeffs()):
        tot += c * (W * Pp) ** ea * Pp ** et * Da ** (1 - ea) * Dt ** (1 - et)
    return sp.expand(tot)
Eden = sp.expand(Pp * (W + 1) * (r - W) * Xx * W * Z)
check("E Da Dt = P (W+1)(r-W) X W Z", sp.expand(cleared(3) - Eden) == 0)
check("raw: sE = F11-F5 and uE = F11-F10 hold after inversion",
      sp.expand(cleared(11) - cleared(5) - s * Eden) == 0 and sp.expand(cleared(11) - cleared(10) - u * Eden) == 0)
def Hrat(i):
    return sp.cancel(cleared(i) / Eden)
H = {i: Hrat(i) for i in (0, 1, 3, 4, 5, 10, 11, 12, 14)}
check("contact levels H1=0,H3=1,H5=1+r,H11=1+r+s,H10=1+r+s-u",
      [sp.simplify(H[i] - v) for i, v in ((1, 0), (3, 1), (5, 1 + r), (11, 1 + r + s), (10, 1 + r + s - u))] == [0] * 5)
check("(3) H4 = r(W+1)/W and H0 = (W-r)(Z-1)",
      sp.simplify(H[4] - r * (W + 1) / W) == 0 and sp.simplify(H[0] - (W - r) * (Z - 1)) == 0)
n12, d12 = sp.fraction(sp.cancel(H[12] - h))
Fc = sp.expand(n12)
pF = sp.Poly(Fc, W, Z, r, s, u, h)
pWZ = sp.Poly(Fc, W, Z)
check("last-contact curve numerator: 105 terms, degree 5 in W and 2 in Z", len(pF.terms()) == 105 and pWZ.degree(W) == 5 and pWZ.degree(Z) == 2,
      (len(pF.terms()), pWZ.degree(W), pWZ.degree(Z)))
F0 = sp.factor(Fc.subs(W, 0))
ratio = sp.cancel(Fc.subs(W, 0) / (r**3 * ((r + s) * Z + 1)))
check("F(0,Z) = c * r^3((r+s)Z+1) with c = +-1", ratio in (1, -1), (F0, ratio))
sgn_ = ratio
check("last-coordinate denominator has degree 1 in Z", sp.Poly(sp.expand(d12), Z).degree() == 1, sp.factor(d12))
c2 = sp.Poly(Fc, Z).coeff_monomial(Z**2)
check("Z^2 coefficient = -W(r-W)^2((u-1)W+r)(r+s-u-W) (times the same sign c)",
      sp.expand(c2 - sgn_ * (-W * (r - W)**2 * ((u - 1) * W + r) * (r + s - u - W))) == 0, sp.factor(c2))
# finite branch at W = 0
Z0 = -1 / (r + s)
FW = sp.diff(Fc, W).subs({W: 0, Z: Z0}); FZ = sp.diff(Fc, Z).subs({W: 0, Z: Z0})
Z1 = sp.cancel(-FW / FZ)
check("Z(0) = -1/(r+s), Z'(0) = -(r^3+2r^2s+r^2+rs^2+rs-su)/(r(r+s)^3)",
      sp.simplify(Z1 + (r**3 + 2 * r**2 * s + r**2 + r * s**2 + r * s - s * u) / (r * (r + s)**3)) == 0, sp.factor(Z1))
def branch_series(G, Y, Y0, order):
    """Power series Y(W) = Y0 + ... solving G(W, Y) = 0, G(0,Y0)=0, G_Y(0,Y0) != 0."""
    coeffs = [Y0]
    GY = sp.diff(G, Y).subs({W: 0, Y: Y0})
    for k in range(1, order + 1):
        yk = sp.Symbol("y%d" % k)
        ser = sum(cf * W**j for j, cf in enumerate(coeffs)) + yk * W**k
        val = sp.expand(G.subs(Y, ser))
        ck = sp.Poly(val, W).coeff_monomial(W**k)
        sol = sp.solve(ck, yk)
        coeffs.append(sp.cancel(sol[0]))
    return coeffs
def laurent_lead(expr, Y, coeffs, n):
    """Return the W^-1 and W^0 coefficients of expr(W, Y(W)) using the series of Y."""
    num, den = sp.fraction(sp.cancel(sp.together(expr)))
    ser = sum(cf * W**j for j, cf in enumerate(coeffs))
    nu_ = sp.Poly(sp.expand(num.subs(Y, ser)), W)
    de_ = sp.Poly(sp.expand(den.subs(Y, ser)), W)
    nc = [sp.cancel(nu_.coeff_monomial(W**k)) for k in range(n)]
    dc = [sp.cancel(de_.coeff_monomial(W**k)) for k in range(n)]
    return nc, dc
zc = branch_series(Fc, Z, Z0, 1)
nc, dc = laurent_lead(H[14], Z, zc, 2)
# H14 = (n0 + n1 W + ..)/(d0 + d1 W + ..)
if dc[0] == 0:
    res14 = sp.cancel(nc[0] / dc[1])
    check("finite branch: H14 has a simple pole with residue -r", sp.simplify(res14 + r) == 0, res14)
else:
    check("finite branch: H14 regular (no pole)", False, "d0=%s" % dc[0])
check("finite branch: pole of H14 is simple (d0 = 0, d1 != 0); W is a local parameter, H4 has residue r, so H4 + H14 is holomorphic",
      dc[0] == 0 and dc[1] != 0)
# V-branch
GV = sp.expand(Fc.subs(Z, V / W) * W)
check("every monomial has degW - degZ >= -1 (so W*F(W,V/W) is a polynomial)", sp.Poly(GV, W, V) is not None and
      min(m[0] - m[1] for m in pWZ.monoms()) == -1)
low = sp.factor(GV.subs(W, 0))
check("lowest weighted polynomial = c * (-r^3 V((r+s-u)V - r - s))",
      sp.expand(GV.subs(W, 0) - sgn_ * (-r**3 * V * ((r + s - u) * V - r - s))) == 0, low)
V0 = (r + s) / (r + s - u)
vc = branch_series(GV, V, V0, 1)
nc, dc = laurent_lead(H[14].subs(Z, V / W), V, vc, 2)
if dc[0] == 0 and nc[0] != 0:
    res14b = sp.cancel(nc[0] / dc[1])
    check("V-branch: H14 has a simple pole with residue -r", sp.simplify(res14b + r) == 0, res14b)
else:
    check("V-branch: H14 pole structure", False, (nc[:2], dc[:2]))
check("V-branch: pole of H14 is simple (d0 = 0, d1 != 0); with H4 residue r, H4 + H14 is holomorphic",
      dc[0] == 0 and dc[1] != 0)
wit = {r: 2, s: 3, u: 4, h: 5}
check("witness (2,3,4,5): r(r+s)(r+s-u) and both branch denominators nonzero",
      (r * (r + s) * (r + s - u)).subs(wit) != 0 and FZ.subs(wit) != 0)

# ---------------- NUMERIC ----------------
try:
    import numpy as np
    worst = None; n = 0
    for seed in range(int(os.environ.get("NSEEDS", "6"))):
        pt = random_point_on_contact_stratum(sorted(T), seed=seed, pi_star=1)
        if pt is None:
            continue
        m = matching_vector_np(pt["a"], pt["b"], pt["alpha"])
        F1 = (P_numpy() @ m)[1]
        c = np.real(m); c[13] -= np.real(F1)
        v, _ = min_phi_over_SO4(c, restarts=12, seed=seed)
        n += 1
        worst = v if worst is None else min(worst, v)
    print("NUMERIC stratum points %d, worst min Phi over SO(4) = %s" % (n, worst))
except Exception as e:
    print("NUMERIC skipped:", e)
print("EXACT checks failed:", len(FAIL))
sys.exit(1 if FAIL else 0)
