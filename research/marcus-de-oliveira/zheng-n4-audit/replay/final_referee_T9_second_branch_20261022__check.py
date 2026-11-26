# Replay for final_referee_T9_second_branch_20261022 (region x<0, t>1 of paper case T_10).
# Own code only; nothing from the supplement is read or run.  EXACT checks are identities in
# sympy's rational-function field over QQ; NUMERIC checks are floating point.
import sys, itertools
from sympy import field, QQ

FAILS = []
def chk(name, ok, level="EXACT"):
    print(("PASS" if ok else "FAIL"), level, name, flush=True)
    if level == "EXACT" and not ok:
        FAILS.append(name)

PERMS = list(itertools.permutations(range(4)))
NEED = [0, 2, 9, 10, 14, 18, 19, 20, 22]
def Rlist(K, v, x, t, r, u):
    q = 1 + t*(v-1); Y = x - q
    Qv = u*t*(v-1)*(v-x) + Y*(v*(r-t) + t*x)
    z = r*v*Y/Qv; w = z*(1 - t*(v-x)/r)
    A = [z, K(0), K(1), x]; C = [w, None, v, q]
    F = {}
    for j in NEED:
        p = K(-1)
        for i in range(4):
            if PERMS[j][i] != 1:
                p *= (A[i] - C[PERMS[j][i]])
        F[j] = p
    return {j: (F[j] - F[0])/((z - w)*(v - 1)) for j in NEED}, Qv
# ---------- Lemma 1 (rho = -1 scale: gaps = Im R / Im v) ----------
K, x, t, r, u, v, vb = field("x,t,r,u,v,vb", QQ)
Ra, Qa = Rlist(K, v, x, t, r, u); Rb, Qb = Rlist(K, vb, x, t, r, u)
g = {j: (Ra[j] - Rb[j])/(v - vb) for j in NEED}
d = g[2] - g[18]
h = 1 - t; k = x + t - 1; L = k + u
N0 = v*vb; Nx = (v-x)*(vb-x); NY = (k - t*v)*(k - t*vb)
M = -x*h*(k*L - t*L*(v+vb) + t*(t+u)*N0)
chk("L1: g10 N0 NY = M = -x h {kL - tL(v+vb) + t(t+u)N0}", g[10]*N0*NY == M)
D2 = -x*L*k + t*x*L*(v+vb) - t*(t*x+u)*N0
chk("L1: d |Qv|^2 = r^2 D2 + r t u L Nx", d*Qa*Qb == r**2*D2 + r*t*u*L*Nx)
chk("L1: D2 = M/h - t u (1-x) N0", D2 == M/h - t*u*(1-x)*N0)
chk("L1: D2 + L NY = M/x", D2 + L*NY == M/x)
chk("L1: g9 = x y {r/Nx - t u/NY}", g[9] == x*(x-1)*(r/Nx - t*u/NY))
chk("L1: g10 = g20 + d", g[10] == g[20] + d)
r0 = t*u*Nx/NY
chk("L1: d|Qv|^2/r = (r-r0) D2 + r0 (D2 + L NY)", d*Qa*Qb/r == (r - r0)*D2 + r0*(D2 + L*NY))
# ---------- Lemma 2 and (2),(3) in zeta coordinates ----------
K2, X, T, r2, U, ze, zb = field("X,T,r,U,ze,zb", QQ)
x2 = -X; t2 = T + 1; u2 = -U; A0 = X + 1
Ra, _ = Rlist(K2, (ze - x2)/(ze - 1), x2, t2, r2, u2); Rb, _ = Rlist(K2, (zb - x2)/(zb - 1), x2, t2, r2, u2)
dq = {j: (Ra[j] - Rb[j])/(ze - zb) for j in NEED}
gn = {j: dq[j]/dq[2] for j in NEED}          # normalised g2 = 1
N = ze*zb
Qz = lambda s: r2*(s+T)*(s+X) - t2*A0*s*(s+T-U)
Kz = lambda s: -r2*((A0+U)*s + T*A0 - U)/Qz(s)
G = r2/N; H = t2*U/((ze+T)*(zb+T))
m = (Kz(ze) - Kz(zb))/(ze - zb)
dd = -(ze*Kz(ze) - zb*Kz(zb))/(ze - zb)
chk("(2) g14 = T m - d", gn[14] == T*m - dd)
chk("(2) g19 = m + G", gn[19] == m + G)
chk("(2) g9 = X (G + H)", gn[9] == X*(G + H))
chk("(3) g18 = 1 - d", gn[18] == 1 - dd)
chk("(3) g22 = (X m - X H - d)/A0", gn[22] == (X*m - X*H - dd)/A0)
Vq = T*H - t2; C0 = G*N + A0*Vq; S0 = X*G + A0*H
chk("(4) d C0 = N (m S0 + G H)", dd*C0 == N*(m*S0 + G*H))
chk("Q/(zeta(zeta+T)) = C0 + S0 conj(zeta)", Qz(ze)/(ze*(ze+T)) == C0 + S0*zb)
chk("K = -d - m conj(zeta)", Kz(ze) == -dd - m*zb)
# ---------- Theorem 1 algebra (independent variables) ----------
K3, Gs, Hs, ms, ds, Xs, Ts, Ns, Vs = field("G,H,m,d,X,T,N,V", QQ)
A0s = Xs + 1; ts = Ts + 1; L0s = ts*ms + Gs; D0s = Ts*ms*(ms + Gs) - L0s*ds
C0s = Gs*Ns + A0s*Vs; S0s = Xs*Gs + A0s*Hs
chk("(5): (Xm-XH)C0 - N(mS0+GH) = A0[m(XV-NH) - H(GN+XV)]",
    (Xs*ms - Xs*Hs)*C0s - Ns*(ms*S0s + Gs*Hs) == A0s*(ms*(Xs*Vs - Ns*Hs) - Hs*(Gs*Ns + Xs*Vs)))
as_ = Ts*ms - ds; cs = ms + Gs; Bs = Xs*(Gs + Hs)
chk("capacity difference aB(c-d) - d^2(a+B+c+d) = X(G+H) D0 - L0 d^2",
    as_*Bs*(cs - ds) - ds**2*(as_ + Bs + cs + ds) == Xs*(Gs + Hs)*D0s - L0s*ds**2)
chk("T m (G+H) - L0 = G(Tm-1) + m(TH-t)", Ts*ms*(Gs + Hs) - L0s == Gs*(Ts*ms - 1) + ms*(Ts*Hs - ts))
chk("(G+H) D0 - L0 d (m-H) = (m+G){T m (G+H) - L0 d}", (Gs + Hs)*D0s - L0s*ds*(ms - Hs) == (ms + Gs)*(Ts*ms*(Gs + Hs) - L0s*ds))
# ---------- NUMERIC ----------
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
    return g/g[..., 2:3]
rng = np.random.default_rng(7); nS = 400000
Xn, Tn, rn, Un = np.exp(rng.uniform(-3, 3, (4, nS)))
vn = rng.normal(0, 2, nS) + 1j*rng.normal(0, 2, nS)*np.exp(rng.uniform(-3, 1, nS))
gq = gaps(-Xn, 1 + Tn, rn, -Un, vn); dn = gq[:, 2] - gq[:, 18]
ok = np.all(gq > -1e-12, axis=1) & (dn > 0)
a_, B_, c_ = gq[:, 14], gq[:, 9], gq[:, 19]
cap = (a_*B_*(c_ - dn) - dn**2*(a_ + B_ + c_ + dn))[ok]
chk("region II u<0 sample (seed 7, %d draws, %d supporting with d>0): min capacity difference %.3e"
    % (nS, ok.sum(), cap.min() if ok.any() else float('nan')), ok.any() and cap.min() > 0, "NUMERIC")
gq = gaps(-Xn, 1 + Tn, rn, Un, vn); dn = gq[:, 2] - gq[:, 18]
ok = np.all(gq > -1e-12, axis=1) & (dn > 1e-12)
chk("region x<0, t>1, u>0 sample (Lemma 1 says empty): %d supporting with d>0" % ok.sum(), ok.sum() == 0, "NUMERIC")
print("RESULT", "FAIL" if FAILS else "ALL EXACT CHECKS PASS", FAILS)
sys.exit(1 if FAILS else 0)
