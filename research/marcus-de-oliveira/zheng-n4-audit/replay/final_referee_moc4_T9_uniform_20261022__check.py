# Replay for final_referee_moc4_T9_uniform_20261022 (primary proof of paper case T_10).
# Own code only; nothing from the supplement is read or run.  Exact checks use sympy's
# rational-function field over QQ (automatic gcd cancellation), so "==" is an exact identity.
# Numeric checks (marked NUMERIC) are floating point.  Exit code 1 iff an EXACT check fails.
import os, sys, itertools
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib"))
from moc4lib import (F_VALUES, F_IS_MINOR, P, rank, generators, CONTACT_SETS, find_g,
                     G_ELEMENTS, g_act_values, P_numpy, leibniz_weights_np, so4_from_params)
from sympy import field, QQ

FAILS = []
def chk(name, ok, level="EXACT"):
    print(("PASS" if ok else "FAIL"), level, name, flush=True)
    if level == "EXACT" and not ok:
        FAILS.append(name)

def lin(terms):
    out = [0] * 24
    for c, j in terms:
        for s in range(24):
            out[s] += c * F_VALUES[j][s]
    return out

# ---------------- A. contact set and generators ----------------
T9 = {0, 1, 3, 5, 7, 8}
chk("document set {0,1,3,5,7,8} equals paper T_10", CONTACT_SETS[10] == frozenset(T9))
G9 = generators(T9)
chk("generators(T) = {0,1,2,3,4,5,6,11,27,28}", G9 == [0, 1, 2, 3, 4, 5, 6, 11, 27, 28])
chk("only f28 is non-minor", [j for j in G9 if not F_IS_MINOR[j]] == [28])
dval = {j: F_VALUES[j][2] - F_VALUES[j][18] for j in G9}
chk("d = g2 - g18 is 1 on f28, -1 on f0,f3, 0 on the other generators",
    dval == {0: -1, 1: 0, 2: 0, 3: -1, 4: 0, 5: 0, 6: 0, 11: 0, 27: 0, 28: 1})
chk("f28 + f0 = f1 + f2 + f5 + f27", lin([(1, 28), (1, 0)]) == lin([(1, 1), (1, 2), (1, 5), (1, 27)]))
chk("f28 + f3 = f4 + f5 + f27", lin([(1, 28), (1, 3)]) == lin([(1, 4), (1, 5), (1, 27)]))
coef3 = [(lambda g: g[22], 1), (lambda g: g[19], 2), (lambda g: g[14], 4), (lambda g: g[20], 5),
         (lambda g: g[9], 6), (lambda g: g[6], 11), (lambda g: g[18], 27), (lambda g: g[2] - g[18], 28)]
chk("(3) g = g22 f1 + g19 f2 + g14 f4 + g20 f5 + g9 f6 + g6 f11 + g18 f27 + d f28 on every generator",
    all(lin([(cf(list(F_VALUES[j])), k) for cf, k in coef3]) == list(F_VALUES[j]) for j in G9))
chk("generators span the 8-dim space {Pc : (Pc)_T = 0}",
    14 - rank([P[s] for s in sorted(T9)]) == 8 and rank([list(F_VALUES[j]) for j in G9]) == 8)
chk("f28 is positive at index 2 (so d > 0 forces g2 > 0)", F_VALUES[28][2] == 1)
chk("real-v stratum: all generators of {0..5} u T are minor (auditor's argument)",
    all(F_IS_MINOR[j] for j in generators(set(range(6)) | T9)))

# ---------------- C. the spectral chart (exact) ----------------
PERMS = list(itertools.permutations(range(4)))
K, x, t, r, u, v, vb, z, w, q = field("x,t,r,u,v,vb,z,w,q", QQ)
def Fs(A, C):   # C[1] is the point at infinity: homogeneous factor -1
    out = []
    for s in PERMS:
        p = K(-1)
        for i in range(4):
            if s[i] != 1:
                p *= (A[i] - C[s[i]])
        out.append(p)
    return out
F = Fs([z, K(0), K(1), x], [w, None, v, q])
chk("F1-F0 = (z-w)(q-v)(x-1)", F[1] - F[0] == (z - w) * (q - v) * (x - 1))
chk("F3-F5 = -(z-w)(q-v)", F[3] - F[5] == -(z - w) * (q - v))
chk("F5-F0 = -x(z-w)(v-1)", F[5] - F[0] == -x * (z - w) * (v - 1))
qq = 1 + t * (v - 1); Y = x - qq
Fq = Fs([z, K(0), K(1), x], [w, None, v, qq])
R = [(f - Fq[0]) / ((z - w) * (v - 1)) for f in Fq]
rr = t * z * (v - x) / (z - w); uu = -Y * (v * z - v + w - z) / ((v - 1) * (z - w))
chk("r = R7 - R1 = t z (v-x)/(z-w)", R[7] - R[1] == rr)
chk("u = R8 = -Y(vz-v+w-z)/((v-1)(z-w))", R[8] == uu)
Qg = uu * t * (v - 1) * (v - x) + Y * (v * (rr - t) + t * x)
chk("inverse map: z Q = r v Y (so Q != 0 on the distinct chart)", z * Qg == rr * v * Y)
chk("inverse map: w = z[1 - t(v-x)/r]", w == z * (1 - t * (v - x) / rr))
def chart(vv):
    q_ = 1 + t * (vv - 1); Y_ = x - q_
    Q_ = u * t * (vv - 1) * (vv - x) + Y_ * (vv * (r - t) + t * x)
    z_ = r * vv * Y_ / Q_; w_ = z_ * (1 - t * (vv - x) / r)
    F_ = Fs([z_, K(0), K(1), x], [w_, None, vv, q_])
    return [(f - F_[0]) / ((z_ - w_) * (vv - 1)) for f in F_], Q_, Y_, q_
Ra, Qa, Ya, qa = chart(v); Rb, Qb, Yb, qb = chart(vb)
for j, c in {0: 0, 1: (t - 1) * (x - 1), 3: 1 - t - x, 5: -x, 7: r + (t - 1) * (x - 1), 8: u}.items():
    chk("contact response R%d" % j, Ra[j] == c)
gh = [(Ra[j] - Rb[j]) / (v - vb) for j in range(24)]      # Im R_j / Im v
y = x - 1; h = 1 - t
nrm = lambda a: (v - a) * (vb - a)
chk("(5) g2^ = -y/|v-1|^2", gh[2] == -y / nrm(1))
chk("(5) g6^ = r y h/(t|v-x|^2)", gh[6] == r * y * h / (t * nrm(x)))
chk("(5) g9^ = x y {r/|v-x|^2 - t u/|Y|^2}", gh[9] == x * y * (r / nrm(x) - t * u / (Ya * Yb)))
chk("(6) (x+t-1) g10^ = -x h {(x+t+u-1)/|v|^2 + t u y/|Y|^2}",
    (x + t - 1) * gh[10] == -x * h * ((x + t + u - 1) / (v * vb) + t * u * y / (Ya * Yb)))
d = gh[2] - gh[18]
chk("(7) g14 + d + t g6 + h g19 = 0", gh[14] + d + t * gh[6] + h * gh[19] == 0)
chk("(7) g9 + d + x g19 - y g22 = 0", gh[9] + d + x * gh[19] - y * gh[22] == 0)
chk("(7) g10 = g20 + d", gh[10] == gh[20] + d)
X = -x; T = -t
Cv = -r * (Ya + u) / Qa; Cvb = -r * (Yb + u) / Qb
m = (Cv - Cvb) / (v - vb)
chk("x<0: g14 = -(1+X)(1+T) m + T d (rho=-1 scale)", gh[14] == -(1 + X) * (1 + T) * m + T * d)
chk("x<0: g22 = -d + X m - X T u/|Y|^2", gh[22] == -d + X * m - X * T * u / (Ya * Yb))
Pv = qa * u * (v - x) + Ya * ((r + 1 - t) * v - x * (1 - t))
chk("R2 - R18 = -r P/Q + x(1-t) (real constant)", Ra[2] - Ra[18] == -r * Pv / Qa + x * (1 - t))
chk("R14 = u - r Y (Y+u)/Q", Ra[14] == u - r * Ya * (Ya + u) / Qa)
chk("P = Q + (Y+u)(v-x), hence P(theta) = (Y(theta)+u)(theta-x) at roots of Q", Pv == Qa + (Ya + u) * (v - x))
Qp = Qa.numer.as_expr(); import sympy as sp
Qpoly = sp.Poly(Qp, sp.Symbol("v"))
xs, ts, rs, us = sp.symbols("x t r u")
Xs, Ts, Us = sp.symbols("X T U")
sub = {xs: -Xs, ts: -Ts, us: -Us}
chk("Q leading coefficient = T(r+T+U) (x=-X,t=-T,u=-U)", sp.expand(Qpoly.LC().subs(sub) - Ts * (rs + Ts + Us)) == 0)
chk("Q constant coefficient = -T X (1+X+T+U)", sp.expand(Qpoly.coeff_monomial(1).subs(sub) + Ts * Xs * (1 + Xs + Ts + Us)) == 0)
chk("Q(x) = r x Y(x)", sp.expand(Qp.subs(sp.Symbol("v"), xs) - rs * xs * (xs - 1 - ts * (xs - 1))) == 0)
a0 = (1 + Xs + Ts) / Ts
Qa0 = sp.together(Qp.subs(sub).subs(sp.Symbol("v"), a0))
chk("Q(a0) = u t (a0-1)(a0-x) with Y(a0)=0", sp.simplify(Qa0 - (-Us) * (-Ts) * (a0 - 1) * (a0 + Xs)) == 0)

# ---------------- D. numeric tests ----------------
import numpy as np
from scipy.optimize import minimize
V = {j: np.array(F_VALUES[j], float) for j in range(82)}
def mins(fun, restarts, seed):
    rng = np.random.default_rng(seed); best = np.inf
    for _ in range(restarts):
        r_ = minimize(fun, rng.uniform(-np.pi, np.pi, 6), method="Nelder-Mead",
                      options={"maxiter": 4000, "xatol": 1e-10, "fatol": 1e-13})
        r_ = minimize(fun, r_.x, method="BFGS")
        best = min(best, r_.fun)
    return best
def soc(th):
    wv = leibniz_weights_np(so4_from_params(th))
    return float(V[28] @ wv) + 2 * np.sqrt(max(float(V[1] @ wv), 0) * max(float(V[11] @ wv), 0))
val = mins(soc, 40, 11)
chk("SOC capacity (1): min_O [Phi_f28 + 2 sqrt(Phi_f1 Phi_f11)] = %.2e >= -1e-9 (40 restarts, seed 11)" % val, val > -1e-9, "NUMERIC")
rng = np.random.default_rng(5); worst = np.inf
for trial in range(25):
    a_, B_ = np.exp(rng.uniform(-3, 3, 2))
    if a_ * B_ <= 1:
        continue
    c_ = (a_ * B_ + a_ + B_ + 1) / (a_ * B_ - 1)          # boundary aB(c-d) = d^2(a+B+c+d), d=1
    vals = a_ * V[4] + B_ * V[6] + c_ * V[2] + V[28]
    worst = min(worst, mins(lambda th: float(vals @ leibniz_weights_np(so4_from_params(th))), 6, trial))
chk("four-coefficient capacity (2) on its boundary: worst min_O Phi = %.2e >= -1e-9 (random (a,B), seed 5)" % worst, worst > -1e-9, "NUMERIC")
print("RESULT", "FAIL" if FAILS else "ALL EXACT CHECKS PASS", FAILS)
sys.exit(1 if FAILS else 0)
