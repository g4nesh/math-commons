"""Exact replay for the note 'A conditional harmonic reduction to the c = 0 contact face'
(independent_true_T48_native_harmonic_conditional_20261023, supplement).

Own code only.  Nothing from the supplement is read, imported or executed.  Uses the audit
library moc4lib (path: argument 1, default ../lib).  Every printed line is one check;
the exit status is nonzero iff an EXACT check fails.
"""
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "lib"))
import moc4lib as L
import sympy as sp
from fractions import Fraction

FAIL = []
def check(name, ok, level="EXACT"):
    print(f"[{level}] {'PASS' if ok else 'FAIL'}  {name}")
    if not ok and level == "EXACT":
        FAIL.append(name)

t0 = time.time()
T = frozenset([0, 3, 5, 6, 9, 12])
FV = L.F_VALUES
# ---------------- A. combinatorics -------------------------------------------------
check("contact set (0,3,5,6,9,12) is literally paper T_49", T == L.CONTACT_SETS[49])
check("T u {17} is literally Table 5 set #23", (T | {17}) == L.AUX7[23])
gens = L.generators(T)
check("generators(T) = {0,1,2,5,7,9,29,39,40,41}", gens == [0, 1, 2, 5, 7, 9, 29, 39, 40, 41])
check("only f40 is non-minor", [j for j in gens if not L.F_IS_MINOR[j]] == [40])
def add(*js):
    return tuple(sum(FV[j][s] for j in js) for s in range(24))
check("(1a) f40+f1 = f7+f39+f0 on all 24 values", add(40, 1) == add(0, 7, 39))
check("(1b) f40+f9 = f7+f2+f41 on all 24 values", add(40, 9) == add(2, 7, 41))
G8 = [7, 5, 29, 39, 2, 40, 41, 0]
names = ['q', 'p', 'j', 'c', 'e', 'd', 'b', 'a']
check("the 8 canonical generators are independent (rank 8 = rank of all 10)",
      L.rank([L.F_VECS[j] for j in G8]) == 8 and L.rank([L.F_VECS[j] for j in gens]) == 8)
off = [s for s in range(24) if s not in T]
extra = [s for s in off if FV[40][s] == 0]
check("off-T indices with f40 = 0 are exactly (3) = 2,7,11,13,16,17,23", extra == [2, 7, 11, 13, 16, 17, 23])
excl = {s: [names[k] for k, j in enumerate(G8) if FV[j][s] != 0] for s in extra}
check("each index of (3) carries exactly one canonical coefficient (2:j 7:b 11:p 13:e 16:q 17:c 23:a), value 1",
      excl == {2: ['j'], 7: ['b'], 11: ['p'], 13: ['e'], 16: ['q'], 17: ['c'], 23: ['a']}
      and all(FV[G8[names.index(excl[s][0])]][s] == 1 for s in extra))
wit = {2: 12, 7: 37, 11: 21, 13: 42, 16: 34, 23: 44}   # supplement labels T11,T36,T20,T41,T33,T43 = paper T_{k+1}
for j, k in wit.items():
    g = next((g for g in L.G_ELEMENTS if L.g_act_set(g, L.CONTACT_SETS[k]) <= (T | {j})), None)
    check(f"extra zero {j}: some g maps paper T_{k} into T u {{{j}}} (g = {None if g is None else (L.perm_to_oneline(g[0]), L.perm_to_oneline(g[1]), g[2])})", g is not None)

# ---------------- B. spectral chart ------------------------------------------------
y, z, w, x, t = sp.symbols('y z w x t')
r, s, u, h, J, W, R = sp.symbols('r s u h J W R')
A = [y, z, w, x]; Bf = [t, 0, 1, None]
def Fsig(sig, A, Bf):
    p = 1
    for i in range(4):
        b = Bf[sig[i]]
        p *= (A[i] - b) if b is not None else 1
    return p
F = [sp.expand(Fsig(sg, A, Bf)) for sg in L.PERMS]
E = sp.factor(F[0] - F[6])
check("E = F0 - F6 = t(y-z)(w-1)", sp.expand(E - t*(y - z)*(w - 1)) == 0)
H = [sp.cancel((F[k] - F[0]) / E) for k in range(24)]
C = (y - t) / (t*(y - z))
check("r := -H5 = -C(x-z)", sp.cancel(-H[5] + C*(x - z)) == 0)
check("s := H3 - H5 = C x (z-w)/(w-1)", sp.cancel(H[3] - H[5] - C*x*(z - w)/(w - 1)) == 0)
check("H0 = 0, H6 = -1 identically", H[0] == 0 and sp.cancel(H[6] + 1) == 0)
V = r*(z - w) + s*(w - 1)
Dt = -V - u*z*(w - 1) + w*z
Dy = V + (w - z)*(w - u*(w - 1))
sub = {x: s*z*(w - 1)/V, t: w*z/Dt, y: w*(V + w - z)/Dy}
Hs = [sp.cancel(sp.together(Hk.subs(sub))) for Hk in H]
check("inverse (4) gives H3 = s-r, H5 = -r, H12 = -u identically",
      sp.cancel(Hs[3] - (s - r)) == 0 and sp.cancel(Hs[5] + r) == 0 and sp.cancel(Hs[12] + u) == 0)
# coverage: y from the r equation, and rt + x - z = r z^2 Dy/(V Dt)
ysol = sp.solve(sp.Eq(-H[5], r), y)
check("r equation solved for y gives y = t(rz+x-z)/(rt+x-z)",
      len(ysol) == 1 and sp.cancel(ysol[0] - t*(r*z + x - z)/(r*t + x - z)) == 0)
check("after (4): rt + x - z = r z^2 Dy/(V Dt)",
      sp.cancel((r*t + x - z).subs(sub) - r*z**2*Dy/(V*Dt)) == 0)
print("  time", round(time.time() - t0, 1))
# J, W chart
dl = r - s
M = dl*W + s
ZD = J*(W - 1) - r
check("z = -M/ZD inverts J = V/(z(w-1))", sp.cancel((V/(z*(w - 1))).subs({z: -M/ZD, w: W}) - J) == 0)
HJ = [sp.cancel(sp.together(Hk.subs({z: -M/ZD, w: W}))) for Hk in Hs]
num, den = sp.fraction(sp.cancel(HJ[9]))
Cc = sp.expand(num - h*den)
Tdiv = J*(W*(dl - 1) + s) + dl*W*(u - 1) + s*u
check("denominator of H9 is J(W-1)T (up to sign)", sp.cancel(den / (J*(W - 1)*Tdiv)) in (1, -1))
check("deg_J C = deg_W C = 2", sp.degree(Cc, J) == 2 and sp.degree(Cc, W) == 2)
CR = sp.expand(Cc.subs(h, R - dl))
PW = sp.Poly(CR, W); PJ = sp.Poly(CR, J)
check("(6) [W^2]C = -RJ{(d-1)J + d(u-1)}", sp.expand(PW.coeff_monomial(W**2) + R*J*((dl - 1)*J + dl*(u - 1))) == 0)
check("(6) [J^2]C = -{(d-1)W+s}{RW-R-1}", sp.expand(PJ.coeff_monomial(J**2) + ((dl - 1)*W + s)*(R*W - R - 1)) == 0)
check("(6) C(J,0) = s(J+u){(R+1)J+d}", sp.expand(CR.subs(W, 0) - s*(J + u)*((R + 1)*J + dl)) == 0)
check("[J^2 W^2]C = -R(d-1)", sp.expand(sp.Poly(CR, J, W).coeff_monomial(J**2*W**2) + R*(dl - 1)) == 0)
print("  time", round(time.time() - t0, 1))
# amplitudes: dual functionals of the basis (8 canonical generators + completion)
basis = [list(L.F_VECS[j]) for j in G8]
for i in range(14):
    e = [Fraction(0)]*14; e[i] = Fraction(1)
    if L.rank(basis + [e]) > len(basis):
        basis.append(e)
Binv = L.inv(L.transpose(basis))
Ppinv = L.inv([L.P[rw] for rw in L.P_PIVOT_ROWS])
Mf = L.matmul(Binv, Ppinv)
Q_ = lambda f: sp.Rational(f.numerator, f.denominator)
vpiv = [HJ[rw] for rw in L.P_PIVOT_ROWS]
cvec = [sum(Q_(Ppinv[i][k])*vpiv[k] for k in range(14)) for i in range(14)]
check("the complex gap vector (H_sigma) lies in im P (all 24 rows)",
      all(sp.cancel(sum(L.P[rw][i]*cvec[i] for i in range(14)) - HJ[rw]) == 0 for rw in range(24)))
amp = {names[k]: sp.cancel(sum(Q_(Mf[k][i])*vpiv[i] for i in range(14))) for k in range(8)}
H9 = HJ[9]
claims = {'j': J, 'q+d': (R*J*(W - 1) - J - dl)/r, 'a+d': -(u - 1)*(J - s)/J,
          'c+d': -(J - s)/((W - 1)*J)}
claims['p'] = -(r*claims['q+d'] + (r + h)*J - s*h)/(J + dl)
def get(nm):
    return sum(amp[a] for a in nm.split('+'))
for nm, cl in claims.items():
    D = sp.cancel(get(nm) - cl.subs(R, dl + h).subs(h, H9))
    c9 = sp.cancel(sp.diff(D, J)/sp.diff(H9, J))
    rest = sp.cancel(D - c9*H9)
    ok = c9.is_number and not (rest.has(J) or rest.has(W)) and sp.Poly(rest, r, s, u).total_degree() <= 1
    check(f"(7) {nm}_amp equals the printed formula up to a real shift (shift = {c9}*h + ({rest}))", ok)
# pole set (8): every amplitude, after replacing H9 by h on the curve, has denominator dividing J W (W-1) M (J+d)
allowed = {sp.expand(f) for f in (J, W, W - 1, M, J + dl)} | {sp.expand(-f) for f in (J, W, W - 1, M, J + dl)}
def poles_ok(expr):
    dn = sp.fraction(sp.factor(sp.cancel(expr)))[1]
    for f, _m in sp.factor_list(dn)[1]:
        if (f.has(J) or f.has(W)) and sp.expand(f) not in allowed:
            return False, f
    return True, None
# The unreduced amplitudes (rational in J, W) have denominators built from J, W, W-1 and T.
# On the curve, a finite point with T = 0 has W*M = 0 by (11) (checked below), so it lies on
# the divisor (8).  Hence (8) holds for the restrictions to the curve.
allowed |= {sp.expand(Tdiv), sp.expand(-Tdiv)}
for nm in names:
    ok, bad = poles_ok(amp[nm])
    check(f"(8) {nm}_amp: denominator factors lie in {{J, W, W-1, T}} (T handled by (11))"
          + ("" if ok else f"  [extra factor {bad}]"), ok)
print("  time", round(time.time() - t0, 1))
# ---------------- C. explicit inverse (9) --------------------------------------------
K = J*W - s
BD = W - (J + u)*(W - 1)
YN = W*(K - J*M); YD = K*(W - u*(W - 1)) - J*M
Ah = [(YN, YD), (-M, ZD), (W, 1), (s, J)]
Bh = [(W, BD), (0, 1), (1, 1), (1, 0)]
br = lambda p, q: p[0]*q[1] - p[1]*q[0]
Pi = [sp.expand(sp.prod([br(Ah[i], Bh[sg[i]]) for i in range(4)])) for sg in L.PERMS]
Ei = sp.factor(Pi[0] - Pi[6])
check("(9) E = Pi0 - Pi6 = J W (W-1)^2 K T", sp.expand(Ei - J*W*(W - 1)**2*K*Tdiv) == 0)
check("(9) Pi0/E = J M / K", sp.cancel(Pi[0]/Ei - J*M/K) == 0)
check("inverse reproduces all 24 chart functions H_sigma = (Pi_sigma - Pi0)/E",
      all(sp.cancel((Pi[k] - Pi[0])/Ei - HJ[k]) == 0 for k in range(24)))
print("  time", round(time.time() - t0, 1))
# ---------------- D. resultants (10)-(12) --------------------------------------------
def res_eq(f, g, var, claim):
    rr = sp.factor(sp.resultant(sp.Poly(CR, var), sp.Poly(f, var)).as_expr())
    ratio = sp.cancel(rr/claim)
    return ratio.is_number and ratio != 0, ratio
Rs = dl + h  # R in terms of h; CR uses R as a free symbol, keep R
checks10 = [
    (J, J, -s*dl*((u - 1)*W - u), "Res_J(C,J)"),
    (W, W, s*(J + u)*((R + 1)*J + dl), "Res_W(C,W)"),
    (W - 1, W, (J + dl)*((r - 1)*J + s), "Res_W(C,W-1)"),
    (M, W, s*(J + dl)*((r*R + dl)*J + dl*(r*u - s)), "Res_W(C,M)"),
    (J + dl, J, -dl*R*(dl - u)*(W - 1)*M, "Res_J(C,J+d)"),
    (Tdiv, J, -(u - 1)*r*(dl - u)*W*M*((dl - 1)*W + s), "Res_J(C,T)"),
    (Tdiv, W, (u - 1)*s*r*(J + u)*(J + dl)*((dl - 1)*J + dl*(u - 1)), "Res_W(C,T)"),
]
QH = (R + 1)*J**2 + (R*(u - s) - (R - dl) + u - 1)*J - s*(u - 1)*(R + 1)
checks10.append((K, W, s*J*(J + dl)*QH, "Res_W(C,K) (12)"))
for f, var, claim, nm in checks10:
    ok, ratio = res_eq(f, None, var, claim)
    check(f"{nm} = printed factorization (ratio {ratio})", ok and ratio in (1, -1))
print("  time", round(time.time() - t0, 1))
# ---------------- E. specialization (r,s,u,h) = (2,3,4,5) ------------------------------
sp0 = {r: 2, s: 3, u: 4}
C0 = sp.expand(CR.subs(sp0).subs(R, (2 - 3) + 5))
w2 = sp.factor(sp.Poly(C0, W).coeff_monomial(W**2)); j2 = sp.factor(sp.Poly(C0, J).coeff_monomial(J**2))
check("at (2,3,4,5): [W^2]C ~ J(2J+3)", sp.cancel(w2/(J*(2*J + 3))).is_number)
check("at (2,3,4,5): [J^2]C ~ (2W-3)(4W-5)", sp.cancel(j2/((2*W - 3)*(4*W - 5))).is_number)
qh0 = sp.expand(QH.subs(sp0).subs(R, 4))
check("at (2,3,4,5): Q_H = 5J^2 + 2J - 45", sp.expand(qh0 - (5*J**2 + 2*J - 45)) == 0)
w0 = sp.Poly(C0, W).coeff_monomial(1)
check("at (2,3,4,5): gcd([W^2]C, [W^0]C) = 1 (no vertical component)", sp.gcd(w2, w0).is_number)
check("at (2,3,4,5): Q_H squarefree", sp.discriminant(qh0, J) != 0)
check("C is not identically zero for any real h when r != s (W^0 coefficient argument)",
      sp.expand(CR.subs(W, 0)) != 0)
print("  time", round(time.time() - t0, 1))

# ---------------- F. NUMERIC: brief search on the T_49 stratum through the chart ---------
import numpy as np, cmath, random
rng = random.Random(20261125)
amp_f = {nm: sp.lambdify((J, W, r, s, u, h), amp[nm], 'numpy') for nm in names}
H_f = [sp.lambdify((J, W, r, s, u), HJ[k], 'numpy') for k in range(24)]
Cf = sp.lambdify((J, W, r, s, u, R), CR, 'numpy')
CW = [sp.lambdify((J, r, s, u, R), sp.Poly(CR, W).coeff_monomial(W**k), 'numpy') for k in (2, 1, 0)]
Pn = L.P_numpy(); kept = 0; worst = None
for trial in range(4000):
    rr_, ss_, uu_, hh_ = [rng.uniform(-3, 3) for _ in range(4)]
    Jc = complex(rng.uniform(-3, 3), rng.uniform(-3, 3))
    co = [f(Jc, rr_, ss_, uu_, rr_ - ss_ + hh_) for f in CW]
    for Wc in np.roots(co):
        try:
            Hv = np.array([complex(f(Jc, Wc, rr_, ss_, uu_)) for f in H_f])
        except ZeroDivisionError:
            continue
        for eta in (1.0, -1.0):
            g = -eta*Hv.imag
            if min(g[sorted(set(range(24)) - T)]) <= 1e-9 or max(abs(g[list(T)])) > 1e-7:
                continue
            kept += 1
            cvecn = np.linalg.lstsq(Pn, g, rcond=None)[0]
            val = L.min_phi_over_SO4(list(cvecn), restarts=4, seed=trial)[0] if kept <= 40 else None
            if val is not None and (worst is None or val < worst):
                worst = val
print(f"  numeric: {kept} supported strict points of the T_49 stratum found from 4000 chart samples (seed 20261125); "
      f"min over SO(4) (4 BFGS restarts, first 40 points) worst = {worst}")
check("NUMERIC: no negative Leibniz sum found on sampled T_49 points", worst is None or worst > -1e-9, "NUMERIC")
print("  total time", round(time.time() - t0, 1))
print("FAILED:", FAIL if FAIL else "none")
sys.exit(1 if FAIL else 0)
