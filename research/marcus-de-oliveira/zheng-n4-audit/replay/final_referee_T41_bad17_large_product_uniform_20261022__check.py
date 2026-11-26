"""Auditor replay for final_referee_T41_bad17_large_product_uniform_20261022 (T_42 bad17 branch, region ab >= d^2).
Own code only; uses ../lib/moc4lib.py. Prints one line per check; exits 1 if an EXACT check fails.
NUMERIC lines are floating-point evidence only."""
import os, sys, random
from fractions import Fraction as Fr
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib"))
import moc4lib as L
import sympy as sp
import flint

FAIL = []
def chk(name, ok, info=""):
    print(("PASS " if ok else "FAIL ") + name + (("  " + str(info)) if info else ""))
    if not ok: FAIL.append(name)

q, p, j, c, e, b, a, h = sp.symbols("q p j c e b a h")
SPEC = [(q, 3), (p, 1), (j, 53), (c, 6), (e, 15), (1, 17), (b, 16), (a, 4)]   # d = 1
g = [sum(co * sp.Rational(L.F_VECS[k][r].numerator, L.F_VECS[k][r].denominator) for co, k in SPEC) for r in range(14)]
w = list(g); w[13] = w[13] + h
Cw = sp.expand(L.cubic_C(w))
C0 = Cw.subs(h, 0); Q0 = sp.expand(sp.diff(Cw, h))
chk("C0 restriction is affine in h", sp.degree(Cw, h) == 1)
Cdoc = -(c + 1)*(a*b + a + a*e + b + b*q + 1 + e + j + q + e*j + e*q + j*q)
chk("C0 restricted cubic C equals the displayed formula", sp.expand(C0 - Cdoc) == 0)

X = a + q + 1; Y = b + e + 1; Tt = e + q + 1; V = a*b + (a + b + 1)*(p + 1)
x, y = sp.symbols("x y")
f = (p - a + b*x)**2 - 4*b*p/(b + 1)*((a + b + 1)*x - a)
gg = (p - b + a*y)**2 - 4*a*p/(a + 1)*((a + b + 1)*y - b)
r00 = -a*b - p*(a + b) + p**2 + 2*p + 2
rr = r00 + a*(a + p + 2)*y + b*(b + p + 2)*x - (2*V - a*b)*x*y
Fe = ((e + 1 + b)*(p - a) + e*b)**2 - 4*b*p*(e - a)*(e + 1 + b)
Fq = Fe.subs({a: b, b: a, q: e, e: q}, simultaneous=True)
xs = {x: e/Y, y: q/X}
chk("C3 Fe / Y^2 = f(e/Y)", sp.simplify(Fe/Y**2 - f.subs(xs)) == 0)
chk("C3 Fq / X^2 = g(q/X)", sp.simplify(Fq/X**2 - gg.subs(xs)) == 0)
R0 = sp.expand(r00*X*Y + a*(a + p + 2)*q*Y + b*(b + p + 2)*e*X - (2*V - a*b)*e*q)   # = r(e/Y, q/X) XY
chk("C3 R0 = r(e/Y, q/X) XY (polynomial, %d terms)" % len(sp.Poly(R0, q, p, e, b, a).terms()),
    sp.simplify(R0/(X*Y) - rr.subs(xs)) == 0)
numA = sp.Poly(sp.expand(R0**2 - Fe*Fq), q, p, e, b, a)
quoA, remA = sp.div(numA, sp.Poly(4*V, q, p, e, b, a))
A0 = quoA.as_expr()
chk("C3 A0 := (R0^2 - Fe Fq)/(4V) is a polynomial (zero remainder; %d terms); so R0^2 - 4VA0 = Fe Fq and A0/(X^2Y^2) = (r^2 - fg)/(4V)" % len(quoA.terms()),
    remA.is_zero)
D = X*Y*(c + 1) + Tt*c*j
Gc = sp.expand(V*D**2 + R0*j*D + A0*j**2)
Jn = sp.Symbol("Jn")
G0J = sp.expand(V*(X*Y)**2 + R0*Jn*X*Y + A0*Jn**2)
chk("C2 G(0,J) = V X^2Y^2 + R0 J XY + A0 J^2 (form (2) at c = 0)", sp.expand(Gc.subs(c, 0) - G0J.subs(Jn, j)) == 0)
# D^2 G(0, jXY/D) = V X^2Y^2 D^2 + R0 j X^2Y^2 D + A0 j^2 X^2Y^2 = (XY)^2 G(c,j): polynomial identity
chk("C2 Lemma 1 transfer D^2 G(0, jXY/D) = (XY)^2 G(c,j)",
    sp.expand(V*(X*Y)**2*D**2 + R0*j*X*Y*X*Y*D + A0*j**2*(X*Y)**2 - (X*Y)**2*Gc) == 0)
chk("C2 D - XY = c(XY + T j)", sp.expand(D - X*Y - c*(X*Y + Tt*j)) == 0)
chk("C4 swap (a,b),(x,y): r invariant, f <-> g",
    sp.expand(rr.subs({a: b, b: a, x: y, y: x}, simultaneous=True) - rr) == 0 and
    sp.simplify(f.subs({a: b, b: a, x: y, y: x}, simultaneous=True) - gg) == 0)

# --- pointwise exact test: full 14x14 Hessian determinant vs the compact G ---
import itertools
TENS = {}
ms, Cs = L.cubic_sympy()
for i in range(14):
    for jj in range(14):
        dij = sp.diff(Cs, ms[i], ms[jj])
        for k in range(14):
            vv = sp.diff(dij, ms[k])
            if vv != 0: TENS[(i, jj, k)] = int(vv)
def fdet(M):
    n = len(M)
    F = flint.fmpq_mat(n, n, [flint.fmpq(int(sp.numer(t)), int(sp.denom(t))) for row in M for t in row])
    d = F.det(); return sp.Rational(int(d.p), int(d.q))
rng = random.Random(11)
ratios = set(); npts = 0
for trial in range(25):
    pt = {s: sp.Rational(rng.randint(0, 30), rng.randint(1, 7)) for s in (q, p, j, c, e, b, a)}
    Qv = Q0.subs(pt)
    if Qv == 0: continue
    hv = -C0.subs(pt) / Qv
    wv = [sp.Rational(t.subs(pt)) for t in g]; wv[13] += hv
    H = [[0]*14 for _ in range(14)]
    for (i, jj, k), vv in TENS.items(): H[i][jj] += vv * wv[k]
    dt = fdet(H)
    Gv = Gc.subs(pt)
    if Gv == 0: continue
    ratios.add(-dt * Qv**2 / (4 * Gv**2 * (pt[c] + 1)**2)); npts += 1
chk("C1 at %d random rational points -det H(w) Q^2 / 4 = k (c+1)^2 G(c,j)^2 with one constant k (the factor (c+1)^2 > 0 is not displayed in the document)" % npts, len(ratios) == 1, "k = %s" % ratios)

# --- Lemma 2: exact coefficient tables ---
A_, B_, xi, eta, P_ = sp.symbols("A B xi eta P")
hL = 2*V*(p + (a*b - 1)/((a + 1)*(b + 1)))
Hx = (a + 1)*(b + 1)*(hL + rr)
fN = (b + 1)*f; gN = (a + 1)*gg
Ex = Hx**2 - (a + 1)*(b + 1)*fN*gN
for chart, sub in [("b>=1", {b: 1 + B_, a: 1 + B_ + A_}), ("b<=1", {b: 1/(1 + B_), a: 1 + B_ + A_})]:
    sub2 = dict(sub); sub2.update({x: xi/(1 + xi), y: eta/(1 + eta)})
    for nm, ex in [("H", Hx), ("E", Ex)]:
        num, den = sp.fraction(sp.cancel(sp.together(ex.subs(sub2, simultaneous=True))))
        P1 = sp.Poly(num, A_, B_, p, xi, eta); Pd = sp.Poly(den, A_, B_, p, xi, eta)
        cn = P1.coeffs(); cd = Pd.coeffs()
        sgn = 1 if all(t > 0 for t in cd) else (-1 if all(t < 0 for t in cd) else 0)
        okn = sgn != 0 and all(sgn * t > 0 for t in cn)
        chk("C5 Lemma 2 chart %s: %s numerator has %d coefficients, all of one strict sign matching the (single-signed) denominator" % (chart, nm, len(cn)), okn,
            "denominator = %s" % sp.factor(den))

# --- transport of the T24 coherent capacity ---
V24 = L.F_VALUES
def g_facet(gr, k):
    vals = L.g_act_values(gr, list(V24[k]))
    for kk in range(82):
        if list(V24[kk]) == vals: return kk
src = [59, 53, 55, 67, 54]; tgt = [1, 53, 4, 16, 17]
hits = [gr for gr in L.G_ELEMENTS if [g_facet(gr, k) for k in src] == tgt]
chk("B2 some group element maps [f59,f53,f55,f67,f54] to [f1,f53,f4,f16,f17] in order", bool(hits), "%d elements, e.g. %s" % (len(hits), hits[:2]))
Z5 = set(range(24))
for k in src: Z5 &= L.F_ZERO[k]
print("INFO common zero set of the T24 core rays", sorted(Z5), "T_25 =", sorted(L.CONTACT_SETS[25]), "contains T_25:", L.CONTACT_SETS[25] <= Z5)
stab = [gr for gr in L.G_ELEMENTS if L.g_act_set(gr, L.CONTACT_SETS[42]) == L.CONTACT_SETS[42] and gr != ((0,1,2,3),(0,1,2,3),0)]
chk("C4 the nontrivial stabilizer element of T_42 fixes f1,f53,f6,f17 and swaps f4<->f16, f3<->f15",
    len(stab) == 1 and [g_facet(stab[0], k) for k in [1, 53, 6, 17, 4, 16, 3, 15]] == [1, 53, 6, 17, 16, 4, 15, 3],
    [g_facet(stab[0], k) for k in [1, 53, 6, 17, 4, 16, 3, 15]] if stab else None)

# --- NUMERIC: Proposition 1 capacity on positive roots of the compact G ---
import math
import mpmath as mp
mp.mp.dps = 40
GP = sp.Poly(Gc, j)
coef_f = [sp.lambdify((q, p, c, e, b, a), GP.coeff_monomial(j**k), "mpmath") for k in (2, 1, 0)]
rng = random.Random(3); worst = None; nroot = 0
for trial in range(400):
    aa = rng.uniform(0.2, 8); bb = max(rng.uniform(1/aa, 8), 1/aa)
    pv = rng.choice([0.0, rng.uniform(0, 3)]); qv = rng.choice([0.0, rng.uniform(0, 5)]); ev = rng.choice([0.0, rng.uniform(0, 5)]); cv = rng.choice([0.0, rng.uniform(0, 5)])
    cf = [fn(qv, pv, cv, ev, bb, aa) for fn in coef_f]
    for rt in mp.polyroots(cf, maxsteps=200, extraprec=200):
        if abs(mp.im(rt)) < 1e-25 and mp.re(rt) > 0:
            nroot += 1
            jj = float(mp.re(rt)); Lc = (aa*bb - 1)/((aa + 1)*(bb + 1))
            slack = pv + Lc - 1/jj
            worst = slack if worst is None else min(worst, slack)
chk("NUMERIC capacity (1) 1/j <= p + L at %d positive roots j of the compact G (400 samples, seed 3, ab >= 1)" % nroot,
    worst is None or worst > -1e-9, "min slack = %s" % worst)

# --- NUMERIC: is the five-resource core safe at equality in (1)? (cited T24 coherent capacity) ---
import numpy as np
worstphi = None
rng = random.Random(4)
for trial in range(12):
    aa = rng.uniform(0.3, 6); bb = max(rng.uniform(1/aa, 6), 1/aa); pv = rng.uniform(0, 2)
    Lc = (aa*bb - 1)/((aa + 1)*(bb + 1)); jv = 1/(pv + Lc)
    cvec = [float(pv*L.F_VECS[1][r] + jv*L.F_VECS[53][r] + aa*L.F_VECS[4][r] + bb*L.F_VECS[16][r] + L.F_VECS[17][r]) for r in range(14)]
    m = L.min_phi_over_SO4(cvec, restarts=12, seed=trial)
    m = m[0] if isinstance(m, tuple) else m
    worstphi = m if worstphi is None else min(worstphi, m)
chk("NUMERIC min over SO(4) of Phi for p f1 + j f53 + a f4 + b f16 + f17 at equality in (1): 12 samples x 12 restarts",
    worstphi > -1e-7, "min = %s" % worstphi)
print("SUMMARY", "all EXACT checks passed" if not [f_ for f_ in FAIL if not f_.startswith("NUMERIC")] else "FAILED: %s" % FAIL)
sys.exit(1 if [f_ for f_ in FAIL if not f_.startswith("NUMERIC")] else 0)
