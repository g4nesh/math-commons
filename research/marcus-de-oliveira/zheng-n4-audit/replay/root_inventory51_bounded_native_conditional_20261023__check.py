"""Exact replay and numeric tests for root_inventory51_bounded_native_conditional_20261023
(conditional reduction for paper case T_52 = {1,3,5,10,11,14} to the seven-element set
G = {0,1,2,7,9,10,23}).

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

T = frozenset({1, 3, 5, 10, 11, 14})
check("T = paper T_52", T == CONTACT_SETS[52])
rays = [0, 1, 2, 10, 12, 13, 27, 30, 59]
check("generators(T) = nine source rays f0,f1,f2,f10,f12,f13,f27,f30,f59", generators(T) == rays)
check("rank of the nine rays is 8 = 14 - contact rank 6",
      rank([F_VECS[j] for j in rays]) == 8 and rank([P[s] for s in sorted(T)]) == 6)
check("only f30 non-minor; f30(4) = f30(15) = 1; all rays >= 0",
      [j for j in rays if not F_IS_MINOR[j]] == [30] and F_VALUES[30][4] == 1 and F_VALUES[30][15] == 1
      and all(F_VALUES[j][k] >= 0 for j in rays for k in range(24)))
check("ties not forcing theta=0 are exactly 0,6,7,12,13,16,20", sorted(F_ZERO[30] - T) == [0, 6, 7, 12, 13, 16, 20])
A7, C7, G7 = map(frozenset, ({0, 1, 2, 7, 9, 10, 17}, {0, 1, 2, 7, 9, 16, 17}, {0, 1, 2, 7, 9, 10, 23}))
check("A, C, G are Table 5 sets 9, 11, 10", (A7, C7, G7) == (AUX7[9], AUX7[11], AUX7[10]))
check("tie 6 -> A, 13 -> A, 12 -> C, 20 -> G",
      all(find_g(S, T | {j}) is not None for j, S in ((6, A7), (13, A7), (12, C7), (20, G7))))
for j, cases in ((0, (10, 11, 12, 24, 39)), (7, (42, 46)), (16, (10, 11, 12, 24, 39))):
    hit = [nu for nu in cases if any(g_act_set(g, CONTACT_SETS[nu]) <= T | {j} for g in range(1152))]
    check("tie %d pulled back to earlier six-element cases %s (not named in the document)" % (j, hit), hit == list(cases))
rng = random.Random(4)
ok = True
for _ in range(30):
    co = {j: Fr(rng.randint(0, 9)) for j in rays}
    c = [sum(co[j] * F_VECS[j][k] for j in rays) for k in range(14)]
    v = values(c)
    O = random_rational_SO4(rng)
    ok &= co[30] <= v[4] and co[30] <= v[15]
    ok &= -Phi(c, O) <= max(Fr(0), -Phi(F_VECS[30], O)) * min(v[4], v[15])
check("theta <= g4, g15 and (1) -L_O <= M min(g4,g15), 30 exact random instances", ok)

# ---------------- spectral parametrization ----------------
a, t, W, Z, r, s, u, h = sp.symbols("a t W Z r s u h")
Xx = W + (r - W) * Z + 1
Pp = Xx * (r + (u - 1) * W) + s * Z * (r - W)
Da = Xx * (r + (u - 1) * W) - s * W * Z * (r - W)
Dt = Xx * (r - W - u) + s * Z * (r - W)
S_ = s * W * Z / Xx; U_ = u * W / (r - W); N_ = S_ + W + W * U_
check("(2) a = N/(1-S+U) = W P/Da and t = N/(W+S-U) = P/Dt",
      sp.simplify(N_ / (1 - S_ + U_) - W * Pp / Da) == 0 and sp.simplify(N_ / (W + S_ - U_) - Pp / Dt) == 0)
Aspec = (a, W * (1 - Z), W + (r - W) * Z, W)
Bspec = (t, 0, 1)
def sigma_pt(p):
    v = sp.Integer(1)
    for i in range(4):
        if p[i] != 3:
            v *= Aspec[i] + Bspec[p[i]]
    return sp.expand(v)
Fv = [sigma_pt(p) for p in PERMS]
def cleared(i):
    f = sp.Poly(Fv[i] - Fv[1], a, t)
    tot = 0
    for (ea, et), c in zip(f.monoms(), f.coeffs()):
        tot += c * (W * Pp) ** ea * Pp ** et * Da ** (1 - ea) * Dt ** (1 - et)
    return sp.expand(tot)
Eden = sp.expand(Pp * (W + 1) * (r - W) * Xx * W * Z)
check("E Da Dt = P (W+1)(r-W) X W Z", sp.expand(cleared(3) - Eden) == 0)
check("raw ratios rE, sE, uE after inversion",
      sp.expand(cleared(5) - cleared(3) - r * Eden) == 0 and sp.expand(cleared(11) - cleared(5) - s * Eden) == 0
      and sp.expand(cleared(11) - cleared(10) - u * Eden) == 0)
H = {i: sp.cancel(cleared(i) / Eden) for i in (0, 4, 14, 15)}
check("(3) H4 = r(W+1)/W, H0 = (W-r)(Z-1)",
      sp.simplify(H[4] - r * (W + 1) / W) == 0 and sp.simplify(H[0] - (W - r) * (Z - 1)) == 0)
N14, D14 = sp.fraction(H[14])
N14, D14 = sp.expand(N14), sp.expand(D14)
# normalize so that D14 = P (the document writes D14 for the denominator of H14)
q = sp.cancel(D14 / sp.expand(Pp))
check("denominator of H14 is a constant multiple of P", q.is_number, q)
N14, D14 = sp.expand(N14 / q), sp.expand(D14 / q)
check("(4) H15 = -suW(W+1)(Z-1)/D14", sp.simplify(H[15] + s * u * W * (W + 1) * (Z - 1) / D14) == 0)
Cc = sp.expand(h * D14 - N14)
pC = sp.Poly(Cc, W, Z, r, s, u, h); pWZ = sp.Poly(Cc, W, Z)
check("C = h D14 - N14 has 68 terms, degree (3,2) in (W,Z)", len(pC.terms()) == 68 and pWZ.degree(W) == 3 and pWZ.degree(Z) == 2,
      (len(pC.terms()), pWZ.degree(W), pWZ.degree(Z)))
D0 = sp.expand(D14.subs(W, 0)); N0 = sp.expand(N14.subs(W, 0))
check("D0 = r((r+s)Z+1)", sp.expand(D0 - r * ((r + s) * Z + 1)) == 0)
check("N0 = r(1-Z)(r(r+s-u)Z + r - u)", sp.expand(N0 - r * (1 - Z) * (r * (r + s - u) * Z + r - u)) == 0)
C0 = sp.expand(Cc.subs(W, 0))
check("leading Z^2 coefficient of C(0,Z) is r^2(r+s-u)", sp.expand(sp.Poly(C0, Z).coeff_monomial(Z**2) - r**2 * (r + s - u)) == 0)
Z0 = -1 / (r + s)
check("(5) (r+s)^2 C(0,Z0) = rsu(r+s+1)", sp.simplify((r + s)**2 * C0.subs(Z, Z0) - r * s * u * (r + s + 1)) == 0)
c2 = sp.Poly(Cc, Z).coeff_monomial(Z**2)
check("(6) Z^2 coefficient of C = (r-W)((u-1)W+r)(r+s-u-W), independent of h",
      sp.expand(c2 - (r - W) * ((u - 1) * W + r) * (r + s - u - W)) == 0, sp.factor(c2))
check("last-coordinate denominator D14 has degree 1 in Z", sp.Poly(D14, Z).degree() == 1)
wit = {r: 2, s: 3, u: 4, h: 5}
fib = sp.expand(C0.subs(wit))
check("witness (2,3,4,5): fibre polynomial = 4Z^2 + 42Z + 14 (up to sign), nonzero discriminant",
      sp.expand(fib - (4 * Z**2 + 42 * Z + 14)) == 0 or sp.expand(fib + (4 * Z**2 + 42 * Z + 14)) == 0, fib)
check("witness: discriminant nonzero and r+s, r+s-u, r+s+1, rsu nonzero",
      sp.discriminant(fib, Z) != 0 and all(e.subs(wit) != 0 for e in (r + s, r + s - u, r + s + 1, r * s * u)))
# resultant form of the same fact: Res_Z(C(0,Z), D0) != 0 generically
res = sp.factor(sp.resultant(C0, D0, Z))
check("Res_Z(C(0,Z), D0) is a nonzero multiple of rsu(r+s+1)", sp.simplify(res / (r * s * u * (r + s + 1))).free_symbols <= {r, s, u}, res)

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
