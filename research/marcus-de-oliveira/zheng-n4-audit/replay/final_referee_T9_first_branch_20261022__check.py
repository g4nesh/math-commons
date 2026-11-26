# Replay for final_referee_T9_first_branch_20261022 (region I of paper case T_10).
# Own code only; nothing from the supplement is read or run.  EXACT checks are identities in
# sympy's rational-function field over QQ; NUMERIC checks are floating point.
import os, sys, itertools
from sympy import field, QQ

FAILS = []
def chk(name, ok, level="EXACT"):
    print(("PASS" if ok else "FAIL"), level, name, flush=True)
    if level == "EXACT" and not ok:
        FAILS.append(name)

PERMS = list(itertools.permutations(range(4)))
NEED = [0, 2, 6, 9, 14, 18, 19, 20, 22]
# ---------- spectral identities (1), (2), (3) in zeta = (v-x)/(v-1) ----------
K, X, T, r, U, ze, zb = field("X,T,r,U,ze,zb", QQ)
x = 1 + X; t = -T; u = -U; h = 1 + T; A0 = 1 + X
def Rlist(v):
    q = 1 + t * (v - 1); Y = x - q
    Qv = u * t * (v - 1) * (v - x) + Y * (v * (r - t) + t * x)
    z = r * v * Y / Qv; w = z * (1 - t * (v - x) / r)
    A = [z, K(0), K(1), x]; C = [w, None, v, q]
    F = {}
    for j in NEED:
        p = K(-1)
        for i in range(4):
            if PERMS[j][i] != 1:
                p *= (A[i] - C[PERMS[j][i]])
        F[j] = p
    return {j: (F[j] - F[0]) / ((z - w) * (v - 1)) for j in NEED}
Ra = Rlist((ze - x) / (ze - 1)); Rb = Rlist((zb - x) / (zb - 1))
dq = {j: (Ra[j] - Rb[j]) / (ze - zb) for j in NEED}
g = {j: dq[j] / dq[2] for j in NEED}          # gaps normalised by g2 = 1 (ratio is sign-free)
N = ze * zb
Qz = lambda s: r * (s - h) * (s - A0) - T * X * s * (s - h - U)
Kz = lambda s: r * ((X - U) * s - (h * X - U)) / Qz(s)
n = -(Kz(ze) - Kz(zb)) / (ze - zb)
d = -(ze * Kz(ze) - zb * Kz(zb)) / (ze - zb)
G = r / N; H = T * U / ((ze - h) * (zb - h))
chk("(1) g14 = h n - d", g[14] == h * n - d)
chk("(1) g19 = G - n", g[19] == G - n)
chk("(1) g9 = A0 (H - G)", g[9] == A0 * (H - G))
chk("(1) g18 = 1 - d", g[18] == 1 - d)
chk("(1) g6 = h G/T", g[6] == h * G / T)
chk("(1) g22 = [d + A0 (H - n)]/X", g[22] == (d + A0 * (H - n)) / X)
nz = lambda a: (ze - a) * (zb - a)
chk("(2) g20 = A0 h {T X (1+H) - H|zeta-1|^2}/(T|zeta-A0|^2) - d",
    g[20] == A0 * h * (T * X * (1 + H) - H * nz(1)) / (T * nz(A0)) - d)
xi = (ze + zb) / 2
S = X * H - A0 * G; L = G * (d + H - A0 * n) + X * n * H; M = d * X * (T + h * H)
B0 = T * n * S - G * H; Ec = T * G * X - G * H * h + T * d * S + T * n * X * (T + h * H); EL = Ec * L - G * M * (T * n + H)
chk("(3) N L = M", N * L == M)
chk("(3) 2 xi B0 = Ec - G N (T n + H)", 2 * xi * B0 == Ec - G * N * (T * n + H))
chk("Q/(zeta(zeta-h)) = G N - X(T+hH) + S conj(zeta)", Qz(ze) / (ze * (ze - h)) == G * N - X * (T + h * H) + S * zb)
chk("numerator/(zeta(zeta-h)) = G{X - H(N - 2xi + h)/T} conj(zeta) + G N H",
    r * ((X - U) * ze - (h * X - U)) / (ze * (ze - h)) == G * (X - H * (N - 2 * xi + h) / T) * zb + G * N * H)
chk("K = -d + n conj(zeta)", Kz(ze) == -d + n * zb)

# ---------- Lemma 1, Theorem 1, Theorem 2: polynomial identities ----------
K2, Gv, nv, dv, ev, Xv, Tv, av, cv, kv = field("G,n,d,e,X,T,a,c,k", QQ)
def mk(G, n, d, e, X, T):
    h = 1 + T; A0 = 1 + X; H = G + e
    S = X * H - A0 * G; L = G * (d + H - A0 * n) + X * n * H; M = d * X * (T + h * H)
    B0 = T * n * S - G * H; Ec = T * G * X - G * H * h + T * d * S + T * n * X * (T + h * H); EL = Ec * L - G * M * (T * n + H)
    return dict(H=H, L=L, M=M, B0=B0, EL=EL, F=h * EL - B0 * (M + h ** 2 * L))
h2 = 1 + Tv; A2 = 1 + Xv; k2 = Tv * Xv - 1
Vv = mk(Gv, nv, dv, ev, Xv, Tv); cexp = Gv - nv
chk("L = G d + G c + e(G + X n)", Vv["L"] == Gv * dv + Gv * cexp + ev * (Gv + Xv * nv))
chk("B0 = e(k n - c) - G(c + h n)", Vv["B0"] == ev * (k2 * nv - cexp) - Gv * (cexp + h2 * nv))
F0 = mk(Gv, nv, dv, K2(0), Xv, Tv)["F"]; F1 = mk(Gv, nv, dv, K2(1), Xv, Tv)["F"]; F2 = mk(Gv, nv, dv, K2(2), Xv, Tv)["F"]
chk("F is affine in e", F2 - 2 * F1 + F0 == 0 and Vv["F"] == F0 + ev * (F1 - F0))
nn = (av + dv) / h2; GG = cv + nn
Npoly = (A2**2*av**3 + h2*A2*(Xv+2)*cv*av**2 + av**2*((1-dv)*Xv*h2*A2 + dv*(Tv*Xv**2+3*Xv**2+4*Xv+2)) + h2**2*A2*cv**2*av
         + cv*av*((1-dv)*Xv*h2*(Tv+Xv+2) + 2*dv*h2*A2**2) + av*(dv*(1-dv)*Xv*(Tv+2)*A2 + dv**2*(Tv*Xv**2+3*Xv**2+3*Xv+1))
         + Xv*h2**2*cv**2 + Xv*dv*h2*(Xv+2)*cv + Xv*dv**2*A2)
dF = mk(GG, nn, dv, K2(1), Xv, Tv)["F"] - mk(GG, nn, dv, K2(0), Xv, Tv)["F"]
chk("dF/de = (T/h) N with the displayed N", dF == (Tv / h2) * Npoly)
q = k2 * nn - cv
estar = GG * (cv + h2 * nn) / q
chk("F(e*) = T h G^2 (A0 a + X d){A0 G a + X(a+d+c)}/q",
    mk(GG, nn, dv, estar, Xv, Tv)["F"] == Tv * h2 * GG**2 * (A2 * av + Xv * dv) * (A2 * GG * av + Xv * (av + dv + cv)) / q)
Rf = lambda e_: GG * (A2 * cv + dv) + e_ * (A2**2 * nn + A2 * cv - Xv * dv)
chk("R(e*) = A0 G^2 (A0 a + X d)/q", Rf(estar) == A2 * GG**2 * (A2 * av + Xv * dv) / q)
def P_of(c_, e_, n_, d_, T_, X_):
    h_ = 1 + T_; A_ = 1 + X_; k_ = T_ * X_ - 1
    p2 = 2*A_*h_*(e_+n_) + (k_+2*h_)*d_
    p1 = A_*h_*(e_**2 + (2-k_)*n_*e_ + n_**2) + (k_+2*h_)*d_*e_ + T_*A_*d_**2 - (k_+1)*d_ + (2*h_-T_*k_)*d_*n_
    p0 = -n_*(A_*e_+d_)*(h_*k_*(e_+n_) + T_*(k_+A_*(1-d_)))
    return A_*h_*c_**3 + p2*c_**2 + p1*c_ + p0
Hh = Vv["H"]; L = Vv["L"]; B0 = Vv["B0"]; M = Vv["M"]; EL = Vv["EL"]
LB0F20 = A2*h2*(Tv*Xv*(1+Hh)*L*B0 - Hh*(B0*M - EL + L*B0)) - dv*Tv*(B0*M - A2*EL + A2**2*L*B0)
Rr = Gv*(A2*cexp + dv) + ev*(A2**2*nv + A2*cexp - Xv*dv)
chk("(6) L B0 F20 = -T R P(c) (after N = M/L, 2 xi L B0 = EL)", LB0F20 == -Tv * Rr * P_of(cexp, ev, nv, dv, Tv, Xv))
TT = (kv + 1) / Xv; hh = 1 + TT; Wd = Xv + kv + 1
R0 = Xv**2 + 2*Xv*kv + 2*Xv + kv + 1
Q0 = (1-dv)*A2**2*(kv+1)*Wd + dv*Xv*kv*R0
Q1 = (1-dv)*A2*(kv+1)*Wd + 2*dv*kv*(Xv**2 + Xv*kv + Xv + kv + 1)
Q2 = (1-dv)*A2*Wd**2 + dv*kv*R0
Pend = (-kv*(A2*av + Xv*dv)*ev**2 - (dv*Q0 + A2*av*Q1 + Xv*kv*A2**2*av**2)/(A2*Wd)*ev
        - (dv**2*(1-dv)*(kv+1)*Wd*R0 + dv*(kv+1)*av*Q2 + Xv*dv*kv*A2*(kv+1)*av**2)/(A2*Wd**2))
chk("W = X h when T = (k+1)/X", Wd == Xv * hh)
chk("P(c*) equals the displayed quadratic in e", P_of(kv*dv/(A2*hh), ev, (av+dv)/hh, dv, TT, Xv) == Pend)
chk("Thm 2: (d+hc)(d+A0c) = T X d^2 + (c+d){A0 h c - (TX-1) d}",
    (dv + h2*cv)*(dv + A2*cv) == Tv*Xv*dv**2 + (cv + dv)*(A2*h2*cv - (Tv*Xv - 1)*dv))

# ---------- NUMERIC: sample region I and test the claimed bounds ----------
import numpy as np
def gaps(x, t, r, u, v):
    q = 1 + t*(v-1); Y = x - q
    Q = u*t*(v-1)*(v-x) + Y*(v*(r-t) + t*x)
    z = r*v*Y/Q; w = z*(1 - t*(v-x)/r)
    A = np.stack([z, 0*z, 0*z+1, 0*z+x], -1); C = np.stack([w, 0*w, v, q], -1)
    F = -np.ones(z.shape + (24,), complex)
    for kk, p in enumerate(PERMS):
        for i in range(4):
            if p[i] != 1:
                F[..., kk] *= (A[..., i] - C[..., p[i]])
    g = np.real((1j/((z-w)*(v-1)))[..., None]*(F - F[..., :1]))
    return g / g[..., 2:3]
rng = np.random.default_rng(2026); nS = 400000
Xn, Tn, rn, Un = np.exp(rng.uniform(-3, 3, (4, nS)))
vn = rng.normal(0, 2, nS) + 1j*rng.normal(0, 2, nS)*np.exp(rng.uniform(-3, 1, nS))
gn = gaps(1+Xn, -Tn, rn, -Un, vn); dn = gn[:, 2] - gn[:, 18]
ok = np.all(gn > -1e-12, axis=1) & (dn > 0)
lin5 = ((1+Xn)*(1+Tn)*gn[:, 19] - (Tn*Xn - 1)*dn)[ok]
socd = (gn[:, 6]*gn[:, 22] - dn**2)[ok]
chk("region I sample (seed 2026, %d draws, %d supporting with d>0): min (5) = %.3e, min g6 g22 - d^2 = %.3e"
    % (nS, ok.sum(), lin5.min(), socd.min()), ok.sum() > 0 and lin5.min() >= -1e-9 and socd.min() >= -1e-9, "NUMERIC")
print("RESULT", "FAIL" if FAILS else "ALL EXACT CHECKS PASS", FAILS)
sys.exit(1 if FAILS else 0)
