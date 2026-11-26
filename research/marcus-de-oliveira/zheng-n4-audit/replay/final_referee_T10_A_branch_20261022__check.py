# Replay for final_referee_T10_A_branch_20261022 (branch A < 0 of paper case T_11, document set {0,1,2,4,7,9}).
# Own code only; nothing from the supplement is read or run.  EXACT checks are identities in
# sympy's rational-function field over QQ; NUMERIC checks are floating point.
import sys, itertools
from sympy import field, QQ, factor, Symbol

FAILS = []
def chk(name, ok, level="EXACT"):
    print(("PASS" if ok else "FAIL"), level, name, flush=True)
    if level == "EXACT" and not ok:
        FAILS.append(name)
PERMS = list(itertools.permutations(range(4)))

def chart(K, h, t, s, l, v, need=range(24), points=False):
    H = h - 1; k = 1 + t*H
    z = h*v/(1 + H*v)
    r = h*v*((k-1)*s*(l-1) - k)/((H*l*s*(k-1) - h*(k-1))*v + l*s*(k-1) - h)
    c = (z*(v-r) - t*r*(v-z))/((v-r) - t*(v-z))
    w = v + s*(z-c)*(v-r)/((1-c)*(z-r))
    if points:
        return z, r, c, w
    A = [None, K(0), K(1), z]; C = [w, v, r, c]          # A_0 = infinity: factor +1
    F = {}
    for j in set(need) | {0}:
        pr = K(1)
        for i in range(1, 4):
            pr *= (A[i] - C[PERMS[j][i]])
        F[j] = pr
    den = (z - c)*(v - r)
    return {j: (F[j] - F[0])/den for j in need}

K, h, t, s, l, v, vb = field("h,t,s,l,v,vb", QQ)
H = h - 1; k = 1 + t*H; D = h - k
# chart inverse ratios (the chart is the inverse of these real ratios; derivation is in the B-branch note)
z_, r_, c_, w_ = chart(K, h, t, s, l, v, points=True)
chk("chart inverse: h = z(v-1)/(v(z-1))", h == z_*(v-1)/(v*(z_-1)))
chk("chart inverse: t = (z-c)(v-r)/((r-c)(v-z))", t == (z_-c_)*(v-r_)/((r_-c_)*(v-z_)))
chk("chart inverse: s = (w-v)(1-c)(z-r)/((z-c)(v-r))", s == (w_-v)*(1-c_)*(z_-r_)/((z_-c_)*(v-r_)))
chk("chart inverse: l = (w-r) z/((w-v)(z-r))", l == (w_-r_)*z_/((w_-v)*(z_-r_)))
Ra = chart(K, h, t, s, l, v); Rb = chart(K, h, t, s, l, vb)
gh = [(Ra[j] - Rb[j])/(v - vb) for j in range(24)]      # Im R/Im v; g = -rho gh; rho = -1 below
chk("contacts of the chart are exactly T = {0,1,2,4,7,9} (generic parameters)", [j for j in range(24) if gh[j] == 0] == [0, 1, 2, 4, 7, 9])
g = gh; rho = -1
a = g[15]; b = g[8]; q = g[23]; G = g[3]; V = g[10]; e = g[3] - g[16]
chk("(2) g17 = g6 - b - e", g[17] == g[6] - b - e)
chk("(2) g14 = b - t(b+q+e)", g[14] == b - t*(b + q + e))
chk("(2) g18 = (1-t)(b+q+e)", g[18] == (1 - t)*(b + q + e))
dd = -(a + b + q + e)
chk("(3) q = D a - k d - h b with d = -(a+b+q+e)", q == D*a - k*dd - h*b)
J0 = l*D - h; N0 = k - s*(k-1)*(l-1); Cc = h*(k-1) + l*D
M1 = (D + s*H*(k-1))*v + s*(k-1) - D; M1b = (D + s*H*(k-1))*vb + s*(k-1) - D
M2 = Cc*v - J0; M2b = Cc*vb - J0
chk("(3) J0 g6 = -h b", J0*g[6] == -h*b)
chk("(4) a = -rho l s H/h", a == -rho*l*s*H/h)
chk("(4) G = rho h s D N0/|M1|^2", G == rho*h*s*D*N0/(M1*M1b))
chk("(4) b = -rho h J0 N0/|M2|^2", b == -rho*h*J0*N0/(M2*M2b))
chk("(5) k V = -rho h s/|v-1|^2 - D b", k*V == -rho*h*s/((v-1)*(vb-1)) - D*b)
# normalised a = 1
gN = [g[j]/g[15] for j in range(24)]
b = gN[8]; q = gN[23]; e = gN[3] - gN[16]; f = gN[14]; M = b + q + e; V = gN[10]; G = gN[3]
theta = l*D/h; beta = 1 - theta
chk("forced signs: h(1-f) + e + f = 0", h*(1 - f) + e + f == 0)
chk("t = (b-f)/M", t == (b - f)/M)
chk("g6 = b/beta", gN[6] == b/beta)
Delta = 1 + e - theta*M; alpha = t*H + theta; gamma = Delta + H; W = -k*Delta - D + theta*h
chk("alpha = C/h", alpha == Cc/h)
chk("(8) (1-t)/s = t Delta", (1 - t)/s == t*Delta)
chk("(8) W = -N0 Delta", W == -N0*Delta)
xi = alpha*v + beta; xib = alpha*vb + beta
chk("R = |xi|^2 = t beta W/(theta b)", xi*xib == t*beta*W/(theta*b))
eta = h/H
chk("(9) M1 = D(gamma xi + W)/(Delta alpha)", M1 == D*(gamma*xi + W)/(Delta*alpha))
chk("(9) G = eta Delta alpha^2 W/(theta |gamma xi + W|^2)", G == eta*Delta*alpha**2*W/(theta*(gamma*xi + W)*(gamma*xib + W)))
chk("(9) k V = D{eta alpha^2/(theta |xi-k|^2) - b}", k*V == D*(eta*alpha**2/(theta*(xi - k)*(xib - k)) - b))
# merged stratum C = 0
lm = -h*(k-1)/D
need = [0, 3, 8, 10, 14, 15, 16, 23]
Ra = chart(K, h, t, s, lm, v, need); Rb = chart(K, h, t, s, lm, vb, need)
gm = {j: (Ra[j] - Rb[j])/(v - vb) for j in need}; gm = {j: gm[j]/gm[15] for j in need}
b = gm[8]; q = gm[23]; e = gm[3] - gm[16]; V = gm[10]; G = gm[3]
theta = lm*D/h; beta = 1 - theta; p = -H
chk("C=0: theta = 1 - k", theta == 1 - k)
chk("(7) D = s H theta (1-b)", D == s*H*theta*(1 - b))
chk("(7) k e = b H theta (1 - 2 Re v) - h", k*e == b*H*theta*(1 - (v + vb)) - h)
chk("(7) k{e + (H theta/D) V |v-1|^2} = -b H theta |v|^2", k*(e + (H*theta/D)*V*(v-1)*(vb-1)) == -b*H*theta*v*vb)
wq = p*b*v + p*(1-b) - 1; wqb = p*b*vb + p*(1-b) - 1
chk("merged: G = p(p-1) beta b(1-b)/(theta Nw)", G == p*(p-1)*beta*b*(1-b)/(theta*wq*wqb))
chk("merged: 2 theta Re w = beta e + (2theta-1)(p-1) - p b theta", theta*(wq + wqb) == beta*e + (2*theta-1)*(p-1) - p*b*theta)
chk("merged: theta q = -(p-1)(1-b) + beta(b+e)", theta*q == -(p-1)*(1-b) + beta*(b+e))
chk("merged: s = (p-theta)/(p theta (1-b))", s == (p - theta)/(p*theta*(1 - b)))
# interval / SOS algebra in the gap variables (b, e, f, theta, M)
K3, b, e, f, th, M, dl = field("b,e,f,th,M,dl", QQ)
def derived(M):
    h = -(e+f)/(1-f); H = -(1+e)/(1-f); t = (b-f)/M; k = 1 + t*H; D = H*(1-t)
    al = t*H + th; De = 1 + e - th*M; ga = De + H; W = -k*De - D + th*h; eta = h/H; be = 1 - th
    R = t*be*W/(th*b); SV = eta*al**2/(th*b); SG = eta*De*al**2*W/(th*e)
    LG = (SG - ga**2*R - W**2)/(2*ga*W) if ga != 0 else None
    LV = (R + k**2 - SV)/(2*k)
    return dict(h=h, H=H, t=t, k=k, D=D, al=al, ga=ga, W=W, eta=eta, be=be, R=R, LG=LG, LV=LV)
dg = derived(M)
den = dg["LV"].denom.as_expr()
chk("L_V reduced denominator is 2 M b theta (1+e)(f-1): %s" % factor(den),
    factor(den/(2*Symbol("M")*Symbol("b")*Symbol("th")*(1+Symbol("e"))*(Symbol("f")-1))) == 1 or factor(den/(2*Symbol("M")*Symbol("b")*Symbol("th")*(1+Symbol("e"))*(Symbol("f")-1))) == -1)
Lc = th*(b+e)*M - b*(1+e)
chk("(11) L_G - L_V = -eta L alpha^2/(2 b e theta gamma)", dg["LG"] - dg["LV"] == -dg["eta"]*Lc*dg["al"]**2/(2*b*e*th*dg["ga"]))
chk("L - (e^2 - b) = {theta(b+e) - e} M + e q, q = M - b - e", Lc - (e**2 - b) == (th*(b+e) - e)*M + e*(M - b - e))
chk("f<1: gamma = -(f(1+e) + theta M (1-f))/(1-f)", dg["ga"] == -(f*(1+e) + th*M*(1-f))/(1-f))
Pp = ((1+e-e*th)*(f-b) - (b*(1+e) - th*(b+e)))**2 + 4*b*(1+e)**2*(1-th)*(f-b)
chk("SOS: L_V^2 - R = alpha^2 P/(4 b^2 theta^2 (1+e)^2)", dg["LV"]**2 - dg["R"] == dg["al"]**2*Pp/(4*b**2*th**2*(1+e)**2))
E = 1 + e + dl
dn = derived(E/th)
chk("t<0 branch: t = -theta(f-b)/E", dn["t"] == -th*(f-b)/E)
chk("W = delta - H + theta h - theta(f-b)H(1+delta)/E", dn["W"] == dl - dn["H"] + th*dn["h"] - th*(f-b)*dn["H"]*(1+dl)/E)
chk("W = delta - H + theta h + t H (1+delta)", dn["W"] == dl - dn["H"] + th*dn["h"] + dn["t"]*dn["H"]*(1+dl))
chk("W E = (1+delta)E - theta(f-b)H(1+delta) - beta h E", dn["W"]*E == (1+dl)*E - th*(f-b)*dn["H"]*(1+dl) - dn["be"]*dn["h"]*E)
Hh = dn["H"]; hh = dn["h"]
d0 = derived((1 + e + (1-b)*Hh)/th)
chk("E(delta0) = (f-b)H", 1 + e + (1-b)*Hh == (f-b)*Hh)
chk("W(delta0) = -beta b H", d0["W"] == -d0["be"]*b*Hh)
chk("alpha(delta0) = 0, R(delta0) = beta^2, L_V(delta0) = beta", d0["al"] == 0 and d0["R"] == d0["be"]**2 and d0["LV"] == d0["be"])
WH = th*hh - th*(f-b)*Hh*(1+Hh)/(1+e+Hh)
chk("W(delta=H) = theta h b/f and E(H) = f H", WH == th*hh*b/f and 1 + e + Hh == f*Hh)
# ---------- NUMERIC ----------
import numpy as np
def gaps(h, t, s, l, v):
    H = h - 1; k = 1 + t*H
    z = h*v/(1 + H*v)
    r = h*v*((k-1)*s*(l-1) - k)/((H*l*s*(k-1) - h*(k-1))*v + l*s*(k-1) - h)
    c = (z*(v-r) - t*r*(v-z))/((v-r) - t*(v-z))
    w = v + s*(z-c)*(v-r)/((1-c)*(z-r))
    A = np.stack([0*z, 0*z, 0*z+1, z], -1); C = np.stack([w, v, r, c], -1)
    F = np.ones(z.shape + (24,), complex)
    for kk, pp in enumerate(PERMS):
        for i in range(1, 4):
            F[..., kk] *= (A[..., i] - C[..., pp[i]])
    gg = np.real((1j/((z-c)*(v-r)))[..., None]*(F - F[..., :1]))
    return gg/np.abs(gg[..., 15:16])
rng = np.random.default_rng(2027); nS = 400000
hn, tn, sn, ln = rng.normal(0, 3, (4, nS)); vn = rng.normal(0, 2, nS) + 1j*rng.normal(0, 2, nS)
tot = 0; worst = np.inf
for sg in (1, -1):
    gg = sg*gaps(hn, tn, sn, ln, vn); Aval = gg[:, 16] - gg[:, 3]
    ok = np.all(gg > -1e-12, axis=1) & (Aval < -1e-9) & (gg[:, 8] > 1e-12)
    tot += int(ok.sum())
    if ok.any():
        worst = min(worst, float((gg[ok, 15]*gg[ok, 8] - Aval[ok]**2).min()))
chk("A<0 sample (seed 2027, 2 x %d draws, both signs of mu): %d supporting; min g15 g8 - e^2 = %.3e" % (nS, tot, worst), tot > 0 and worst >= -1e-9, "NUMERIC")
print("RESULT", "FAIL" if FAILS else "ALL EXACT CHECKS PASS", FAILS)
sys.exit(1 if FAILS else 0)
