"""Second-pass independent replay for the supplement document root_moc4_T4_uniform_20261022
(paper case T_4 = {0,1,4,5,7,8}).  Own code only; nothing from the supplement is executed,
imported or read.  Uses the shared audit library ../lib/moc4lib.py.

Run:  python -I -B root_moc4_T4_uniform_20261022__check.py
Prints one line per check; exits nonzero if any EXACT check fails.  NUMERIC lines are
floating-point evidence only.
"""
import os, sys, itertools, importlib.util
from fractions import Fraction
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("moc4lib", os.path.join(HERE, "..", "lib", "moc4lib.py"))
L = importlib.util.module_from_spec(spec); spec.loader.exec_module(L)

FAIL = []
def check(name, ok, info=""):
    print(("PASS " if ok else "FAIL ") + "EXACT " + name + (("  " + info) if info else ""))
    if not ok:
        FAIL.append(name)

V = L.F_VALUES
def lin(coeffs):
    return [sum(c * V[j][s] for j, c in coeffs.items()) for s in range(24)]

# ---------------------------------------------------------------- A. contact set, generators
T4 = frozenset({0, 1, 4, 5, 7, 8})
check("A1 T4 equals paper Appendix E T_4 literally", T4 == L.CONTACT_SETS[4])
gens = L.generators(T4)
check("A2 generators(T4) = {f0..f6,f11,f33,f34}", gens == [0, 1, 2, 3, 4, 5, 6, 11, 33, 34],
      str([L.f_name(j) for j in gens]))
# document table (0-based m_{I,J}) -> paper 1-based n_{I,J}
doc_tab = {0: ("1", "4"), 1: ("12", "34"), 2: ("13", "34"), 3: ("14", "34"), 4: ("1", "3"),
           5: ("12", "24"), 6: ("4", "1"), 11: ("13", "23"), 33: ("12", "13")}
check("A3 document minor table matches Appendix D", all(L.F_TABLE[j] == ("n",) + doc_tab[j] for j in doc_tab))
check("A4 f34 = f0 - f2 + f6 + f33 (values)", list(V[34]) == lin({0: 1, 2: -1, 6: 1, 33: 1}),
      L.f_name(34))
check("A5 relation f0-f1-f2-f3+f4 = 0", all(e == 0 for e in lin({0: 1, 1: -1, 2: -1, 3: -1, 4: 1})))
check("A6 relation -f0+f2-f6-f33+f34 = 0", all(e == 0 for e in lin({0: -1, 2: 1, 6: -1, 33: -1, 34: 1})))

# ---------------------------------------------------------------- C1. decomposition (2) on W
W = L.nullspace([L.P[s] for s in sorted(T4)])
ok2 = ok22 = ok18 = ok15 = True
for cvec in W:
    g = L.values(cvec)
    d = g[2] - g[13]
    dec = lin({0: g[19], 1: g[16], 3: g[14], 5: g[10], 6: g[15], 11: g[6], 33: g[13], 34: d})
    ok2 &= dec == list(g)
    ok22 &= g[22] == g[16] + g[19] + d
    ok18 &= g[18] == g[14] + g[19] + d
    ok15 &= g[15] == g[9] - d
check("C1 dim of {c : (Pc)_T4 = 0} is 8", len(W) == 8)
check("C2 decomposition (2) holds on all of W", ok2)
check("C3 g22 = g16 + g19 + d on W", ok22)
check("C4 g18 = g14 + g19 + d on W", ok18)
check("C5 g15 = g9 - d (a = g9 - d) is linear on W", ok15)

# d <= 0 branch: coefficients after adding A*(rel2) + tau*(rel1), d = -A
gs = {k: sp.Symbol(f"g{k}") for k in (19, 16, 14, 10, 15, 6, 13, 2, 9)}
Asym, tau = sp.symbols("A tau")
coef = {0: gs[19], 1: gs[16], 3: gs[14], 5: gs[10], 6: gs[15], 11: gs[6], 33: gs[13], 34: -Asym, 2: 0, 4: 0}
rel1 = {0: 1, 1: -1, 2: -1, 3: -1, 4: 1}
rel2 = {0: -1, 2: 1, 6: -1, 33: -1, 34: 1}
new = {j: coef[j] + Asym * rel2.get(j, 0) + tau * rel1.get(j, 0) for j in coef}
sub_lin = {gs[15]: gs[9] + Asym, gs[13]: gs[2] + Asym}   # g15 = g9 - d, g13 = g2 - d with d = -A
want = {0: gs[19] - Asym + tau, 1: gs[16] - tau, 2: Asym - tau, 3: gs[14] - tau, 4: tau,
        6: gs[9], 33: gs[2], 34: 0, 5: gs[10], 11: gs[6]}
check("C6 d<=0 replacement coefficients", all(sp.expand(new[j].subs(sub_lin) - want[j]) == 0 for j in want))

# ---------------------------------------------------------------- C2. four-coefficient bound (1)
# transport: g = (id, (0 2 1 3), 0), i.e. sigma -> sigma o (1 2) [0-based], the column swap
gT = ((0, 1, 2, 3), (0, 2, 1, 3), 0)
gv = lambda j: tuple(L.g_act_values(gT, list(V[j])))
check("C7 transport g maps T5 rays f0,f6,f10,f54 to T4 rays f0,f6,f5,f34",
      gv(0) == V[0] and gv(6) == V[6] and gv(10) == V[5] and gv(54) == V[34],
      f"{L.f_name(10)}->{L.f_name(5)}, {L.f_name(54)}->{L.f_name(34)}")
p = sp.symbols("p0:4"); q = sp.symbols("q0:4")
def qmul(a, b):
    a0, a1, a2, a3 = a; b0, b1, b2, b3 = b
    return (a0*b0 - a1*b1 - a2*b2 - a3*b3, a0*b1 + a1*b0 + a2*b3 - a3*b2,
            a0*b2 - a1*b3 + a2*b0 + a3*b1, a0*b3 + a1*b2 - a2*b1 + a3*b0)
E4 = [tuple(1 if k == j else 0 for k in range(4)) for j in range(4)]
cols = [qmul(qmul(p, E4[j]), q) for j in range(4)]
OQ = [[sp.expand(cols[j][i]) for j in range(4)] for i in range(4)]     # x -> p x q
npp = sum(e**2 for e in p); nqq = sum(e**2 for e in q)
OtO = sp.Matrix(OQ).T * sp.Matrix(OQ)
check("C8 O(p,q)^T O(p,q) = |p|^2|q|^2 I", all(sp.expand(OtO[i, j] - (npp * nqq if i == j else 0)) == 0
                                              for i in range(4) for j in range(4)))
def nf(e):
    e = sp.rem(sp.expand(e), p[0]**2 + p[1]**2 + p[2]**2 + p[3]**2 - 1, p[0])
    e = sp.rem(sp.expand(e), q[0]**2 + q[1]**2 + q[2]**2 + q[3]**2 - 1, q[0])
    return sp.expand(e)
def phi_sym(vals, O):
    tot = 0
    for s, sg, vv in zip(L.PERMS, L.SGN, vals):
        if vv:
            tot += sg * vv * O[0][s[0]] * O[1][s[1]] * O[2][s[2]] * O[3][s[3]]
    return tot
def quad(y):
    p0, p1, p2, p3, q0, q1, q2, q3 = y
    A = p0*q3 + p3*q0; B = p1*q2 - p2*q1
    D = p0**2 + p2**2 - q0**2 - q2**2; E = p1**2 - p2**2 - q0**2 + q3**2
    G = p0**2 - p3**2 - q1**2 + q2**2
    return A, B, D, E, G
y = tuple(p) + tuple(q)
A_, B_, D_, E_, G_ = quad(y)
H_ = 2 * (A_**2 + B_**2) - E_ * G_
tg = {0: (A_ - B_)**2, 6: (A_ + B_)**2, 10: D_**2, 54: H_}
check("C9 on unit p,q: Phi_f0=(A-B)^2, Phi_f6=(A+B)^2, Phi_f10=D^2, Phi_f54=H",
      all(nf(phi_sym(V[j], OQ) - tg[j]) == 0 for j in tg))
a, b, c, d = sp.symbols("a b c d")
Z = sp.Matrix([[a + b + 2*d, a - b, 0, 0], [a - b, a + b + 2*d, 0, 2*d], [0, 0, c, d], [0, 2*d, d, d]])
p0, p1, p2, p3, q0, q1, q2, q3 = y
S = [y, (-p1, p0, -p3, p2, q1, -q0, -q3, q2), (-p2, p3, p0, -p1, q2, q3, -q0, -q1),
     (-p3, -p2, p1, p0, q3, -q2, q1, -q0)]
rhs = 0
for yy in S:
    A1, B1, D1, E1, G1 = quad(yy)
    l = sp.Matrix([-yy[0] * B1, yy[0] * A1, -yy[7] * D1, -yy[7] * G1])
    m = sp.Matrix([yy[4] * B1, -yy[4] * A1, -yy[3] * D1, -yy[3] * E1])
    rhs += (l.T * Z * l)[0] + (m.T * Z * m)[0]
lhs = (npp + nqq) * (a * (A_ - B_)**2 + b * (A_ + B_)**2 + c * D_**2 + d * H_)
diff = sp.expand(lhs - rhs)
# mod N = |p|^2 - |q|^2: reduce p0^2 -> |q|^2 - p1^2 - p2^2 - p3^2
diffN = sp.expand(sp.rem(diff, p0**2 - (nqq - p1**2 - p2**2 - p3**2), p0))
check("C10 companion certificate (4): n[...] == sum l^T Z l + m^T Z m  (mod N)", diffN == 0)
# Z PSD under the capacity condition (Schur chain)
S1 = Z[:3, :3] - Z[:3, 3] * Z[3, :3] / d
S1w = sp.Matrix([[a + b + 2*d, a - b, 0], [a - b, a + b - 2*d, -2*d], [0, -2*d, c - d]])
S2 = S1[:2, :2] - S1[:2, 2] * S1[2, :2] / (c - d)
detS2 = sp.factor(sp.simplify(S2.det() - 4 * (a*b*(c - d) - d**2*(a + b + c + d)) / (c - d)))
check("C11 Schur complements of Z: S1 as stated, det S2 = 4[ab(c-d)-d^2(a+b+c+d)]/(c-d)",
      sp.simplify(S1 - S1w) == sp.zeros(3, 3) and detS2 == 0)
# capacity with d>0 forces c>d, a,b>0 (proof-read); T4-side Schur identity from the doc
lhsT4 = d - d**2/(a + d) - d**2/c - d**2/(b + d)
check("C12 doc's arrowhead Schur identity", sp.simplify(lhsT4 - d*(a*b*(c - d) - d**2*(a + b + c + d)) /
                                                        ((a + d)*c*(b + d))) == 0)

# ---------------------------------------------------------------- C3. spectral chart
z, w, v, vb, t, x, r, u = sp.symbols("z w v vb t x r u")
PERMS = L.PERMS
def Fsig(s, aa, cc):
    out = -1
    for i in range(4):
        if s[i] != 2:
            out *= (aa[i] - cc[s[i]])
    return sp.expand(out)
aa = [z, 0, 1, x]; cc = [w, v, None, t * v]
F = [Fsig(s, aa, cc) for s in PERMS]
R = [sp.cancel((F[k] - F[0]) / (v * (z - w))) for k in range(24)]
r_def = z * (w - v) * (1 - t * v) / (v * (z - w))
u_def = -(x - t * v) * (z * (1 + v - w) - v) / (v * (z - w))
claimR = {0: 0, 1: 1 - x, 4: t - x, 5: x * (t - 1), 7: r_def - x + 1, 8: u_def}
check("C13 contact ratios R_sigma on T4 = 0,1-x,t-x,x(t-1),r-x+1,u",
      all(sp.cancel(R[k] - e) == 0 for k, e in claimR.items()))
check("C14 R5 - R4 = -t (R1 - R0); R4 - R1 = t - 1 (basis of the reality argument)",
      sp.cancel(R[5] - R[4] + t * (R[1] - R[0])) == 0 and sp.cancel(R[4] - R[1] - (t - 1)) == 0)
zf = r*v*(t*v - x) / (r*t*v**2 - r*v*x + t**2*v**2 - t*u*v**2 - t*v*x - t*v + u*v + x)
wf = (-t*v + x)*(r - t*v + 1) / (-r*t*v + r*x + t*u*v - u)
sub = {r: r_def, u: u_def}
check("C15 z formula inverts (r,u): z_formula(r(z,w),u(z,w)) = z", sp.cancel(zf.subs(sub) - z) == 0)
check("C16 w formula inverts (r,u)", sp.cancel(wf.subs(sub) - w) == 0)
dz = sp.cancel(sp.denom(zf).subs(sub) + (v - w)*(x - t*v)*(t*v - 1)/(z - w))
dw = sp.cancel(sp.denom(wf).subs(sub) + (z - v)*(x - t*v)*(t*v - 1)/(v*(z - w)))
check("C17 denominators of z,w formulas = -(v-w)(x-tv)(tv-1)/(z-w), -(z-v)(x-tv)(tv-1)/(v(z-w))",
      dz == 0 and dw == 0)
Qs = lambda s: -t*(r + t - u)*s**2 + (r*t + r*x + t*x + t - u)*s - x*(r + 1)
check("C18 Q(v) = -(z-v)(v-w)(x-tv)(tv-1)/(v(z-w)) in original variables (so Q(v) != 0)",
      sp.cancel(Qs(v).subs(sub) + (z - v)*(v - w)*(x - t*v)*(t*v - 1)/(v*(z - w))) == 0)
# chart: R_sigma(v) with z,w substituted; scaled gap G = (R(v)-R(vb))/(v-vb) = Im R / Im v
need = [2, 6, 9, 10, 13, 14, 15, 16, 18, 19, 22]
Rc = {k: sp.cancel(R[k].subs({z: zf, w: wf})) for k in need}
G = {k: sp.cancel((Rc[k] - Rc[k].subs(v, vb)) / (v - vb)) for k in need}
dd = sp.cancel(G[2] - G[13]); bb = G[19]; cc_ = G[10]; aa_ = G[15]
yv, h = x - 1, 1 - t
Hh = h + u; Ee = Hh*(1 - t*x) - h*r*yv; Jj = Hh*(1 - t*x) + h*yv
check("C19 g2 = x/|v|^2 (scaled, rho=-1)", sp.cancel(G[2] - x/(v*vb)) == 0)
check("C20 g6 = r t (x-1)/|tv-1|^2", sp.cancel(G[6] - r*t*(x - 1)/((t*v - 1)*(t*vb - 1))) == 0)
check("C21 c = g10 = -t u x/|tv-x|^2", sp.cancel(cc_ + t*u*x/((t*v - x)*(t*vb - x))) == 0)
check("C22 g9 = (x-1) c and a = g9 - d", sp.cancel(G[9] - (x - 1)*cc_) == 0 and sp.cancel(aa_ - G[9] + dd) == 0)
check("C23 (3a) k = g16 = h c - t(b+d)", sp.cancel(G[16] - (h*cc_ - t*(bb + dd))) == 0)
check("C24 (3b) E e = t y H b - J d, e = g14", sp.cancel(Ee*G[14] - (t*yv*Hh*bb - Jj*dd)) == 0)
check("C25 (4) J = E + hy(r+1), J - xhH = -uy",
      sp.expand(Jj - Ee - h*yv*(r + 1)) == 0 and sp.expand(Jj - x*h*Hh + u*yv) == 0)
check("C26 g22, g18 identities hold in the chart too",
      sp.cancel(G[22] - G[16] - bb - dd) == 0 and sp.cancel(G[18] - G[14] - bb - dd) == 0)
K = t*yv*bb - x*h*dd
Dp = sp.cancel(dd*Qs(v)*Qs(vb)); Pp = sp.cancel(K*Qs(v)*Qs(vb))
def herm(e):
    num, den = sp.fraction(sp.cancel(e))
    if sp.Poly(den, v, vb).total_degree() != 0:
        return None
    P_ = sp.Poly(sp.expand(num / den), v, vb)
    if P_.degree(v) > 1 or P_.degree(vb) > 1:
        return None
    m00 = P_.coeff_monomial(1); m01 = P_.coeff_monomial(v); m10 = P_.coeff_monomial(vb)
    m11 = P_.coeff_monomial(v*vb)
    if sp.expand(m01 - m10) != 0:
        return None
    return sp.Matrix([[m00, m01], [m01, m11]])
MD = herm(Dp); MP = herm(Pp)
check("C27 D = d|Q(v)|^2 and P = K|Q(v)|^2 are Hermitian polynomials of bidegree <= (1,1)",
      MD is not None and MP is not None)
Lq = yv - r + u
if MD is not None and MP is not None:
    ev = lambda M: sp.expand((sp.Matrix([[1, x/t]]) * M * sp.Matrix([1, x/t]))[0])
    check("C28 (6) det M_D = -r^3 t^2 u x y L E", sp.expand(MD.det() + r**3*t**2*u*x*yv*Lq*Ee) == 0)
    check("C29 (6) D(x/t,x/t) = r u x y (hr - Hy)", sp.cancel(ev(MD) - r*u*x*yv*(h*r - Hh*yv)) == 0)
    check("C30 (7) det M_P = -r^2 t^2 u^2 x^2 h y L E", sp.expand(MP.det() + r**2*t**2*u**2*x**2*h*yv*Lq*Ee) == 0)
    check("C31 (7) P(x/t,x/t) = r u^2 x^2 y (y - h)", sp.cancel(ev(MP) - r*u**2*x**2*yv*(yv - h)) == 0)
check("C32 E = (h-y)(H+hy) + hyL", sp.expand(Ee - ((h - yv)*(Hh + h*yv) + h*yv*Lq)) == 0)
check("C33 hr-Hy<=0 bound: H(1-tx) - Hy^2 = H(h-y)(1+y); y - Hy/h + u = -u(y-h)/h",
      sp.expand(Hh*(1 - t*x) - Hh*yv**2 - Hh*(h - yv)*(1 + yv)) == 0 and
      sp.cancel(yv - Hh*yv/h + u + u*(yv - h)/h) == 0)
A0, B0, C0, D0, Y0, X0 = sp.symbols("A0 B0 C0 D0 Y0 X0")
check("C34 a=yc-d => ab(c-d)-d^2(a+b+c+d) = c{ybc - x d(b+d)} (x = y+1)",
      sp.expand(((Y0*C0 - D0)*B0*(C0 - D0) - D0**2*(Y0*C0 - D0 + B0 + C0 + D0)) -
                C0*(Y0*B0*C0 - (Y0 + 1)*D0*(B0 + D0))) == 0)
check("C35 ybc - xd(b+d) - (b+d)/h (tyb - xhd) = yb(c - t(b+d)/h)",
      sp.cancel((yv*B0*C0 - x*D0*(B0 + D0)) - (B0 + D0)/h*(t*yv*B0 - x*h*D0) - yv*B0*(C0 - t*(B0 + D0)/h)) == 0)

# ---------------------------------------------------------------- v real: minor-only entry face
ent = [k for k, s in enumerate(PERMS) if s[0] == 0]
check("C36 sigma(0)=0 permutations are 0..5", ent == [0, 1, 2, 3, 4, 5])
check("C37 F_sigma/(z-w) for sigma(0)=0 is a polynomial in v,t,x with integer coefficients",
      all(sp.cancel(F[k] / (z - w)).free_symbols <= {v, t, x} and sp.denom(sp.cancel(F[k] / (z - w))) == 1 for k in ent))
check("C38 generators of {0..5} (and of T4 u {0..5}) are all minor vectors",
      all(L.F_IS_MINOR[j] for j in L.generators(range(6))) and
      all(L.F_IS_MINOR[j] for j in L.generators(set(T4) | set(range(6)))))

# ---------------------------------------------------------------- D. numeric falsification
import numpy as np
from scipy.optimize import minimize
rng = np.random.default_rng(20261125)
Pn = L.P_numpy()
fv = {j: Pn @ np.array([float(e) for e in L.F_VECS[j]]) for j in (0, 5, 6, 34)}
def minphi(vals, restarts=12):
    f = lambda th: float(vals @ L.leibniz_weights_np(L.so4_from_params(th)))
    best = np.inf
    for _ in range(restarts):
        best = min(best, minimize(f, rng.uniform(-np.pi, np.pi, 6), method="BFGS").fun)
    return best
worst1 = np.inf
for _ in range(40):
    al, be, ga = rng.dirichlet([1, 1, 1])
    dcoef = 1.0
    acoef, bcoef, ccoef = 1/al - 1, 1/be - 1, 1/ga      # boundary of the capacity region
    vals = acoef*fv[6] + bcoef*fv[0] + ccoef*fv[5] + dcoef*fv[34]
    worst1 = min(worst1, minphi(vals) / (acoef + bcoef + ccoef + dcoef))
print(f"NUMERIC D1 inequality (1) on 40 random boundary tuples (12 BFGS restarts each, seed 20261125): "
      f"worst normalized min Phi = {worst1:.3e}")
# chart samples in the d>0 sign region; Mobius s -> 1/(s - s0) sends infinity to 0
zf_l = sp.lambdify((x, t, r, u, v), zf, "numpy"); wf_l = sp.lambdify((x, t, r, u, v), wf, "numpy")
PA = np.array(PERMS)
NS = 200000
half = NS // 2
xs = np.concatenate([1 + rng.exponential(1.0, half), rng.normal(0, 2, half)])
ts = np.concatenate([rng.uniform(0, 1, half), rng.normal(0, 1.5, half)])
rs = np.concatenate([rng.exponential(1.0, half), rng.normal(0, 2, half)])
us = np.concatenate([-rng.exponential(1.0, half), rng.normal(0, 2, half)])
vs = rng.normal(size=NS) + 1j * rng.normal(size=NS)
sg = np.concatenate([-np.ones(half), rng.choice([-1.0, 1.0], half)])
with np.errstate(all="ignore"):
    zs = zf_l(xs, ts, rs, us, vs); ws = wf_l(xs, ts, rs, us, vs)
    kaps = 1j * (sg / np.abs(vs.imag)) / (vs * (zs - ws))
    Ap = np.stack([zs, 0 * zs, 0 * zs + 1, xs + 0j], 1); Cp = np.stack([ws, vs, 0 * vs, ts * vs], 1)
    Fv = np.ones((NS, 24), complex)
    for k, s in enumerate(PERMS):
        for i in range(4):
            if s[i] != 2:
                Fv[:, k] *= Ap[:, i] - Cp[:, s[i]]
    Fv = -Fv
    gaps = (kaps[:, None] * (Fv - Fv[:, :1])).real
    scale = np.abs(gaps).max(1)
    adm = np.isfinite(scale) & (np.abs(gaps[:, list(T4)]).max(1) < 1e-9 * scale) & (gaps.min(1) >= -1e-12 * scale)
idx = np.flatnonzero(adm)
dpos = idx[(gaps[idx, 2] - gaps[idx, 13]) > 1e-9 * scale[idx]]
# sign claims of the d>0 branch at those points (rho normalized to -1 means sg = -1 here)
okc = True
for i in dpos:
    xx, tt, rr, uu = xs[i], ts[i], rs[i], us[i]
    gg = gaps[i] / (-(sg[i] / abs(vs[i].imag)) * vs[i].imag)     # scaled by (-rho)^{-1}
    bb_, dd_ = gg[19], gg[2] - gg[13]
    KK = tt * (xx - 1) * bb_ - xx * (1 - tt) * dd_
    okc &= (xx > 1) and (0 < tt < 1) and (rr > 0) and (uu < 0) and KK > -1e-9 * scale[i]
print(f"NUMERIC D2a chart sampling (seed 20261125, {NS} samples): {len(idx)} admissible T4 points, "
      f"{len(dpos)} with d>0; sign claims x>1, 0<t<1, r>0, u<0, K>0 hold at all d>0 points: {okc}")
worst2 = np.inf; maxrel = 0.0; tested = 0
for i in list(dpos[:25]) + [j for j in idx if j not in set(dpos)][:15]:
    s0 = complex(rng.normal(), rng.normal())
    A_pts = [zs[i], 0, 1, xs[i]]; C_pts = [ws[i], vs[i], None, ts[i] * vs[i]]
    mob = lambda s: 1 / (s - s0)
    a2 = [mob(s) for s in A_pts]; c2 = [mob(s) if s is not None else 0j for s in C_pts]
    F2 = np.array([np.prod([a2[k] - c2[s[k]] for k in range(4)]) for s in PERMS])
    lam = F2 / Fv[i]
    maxrel = max(maxrel, float(np.max(np.abs(lam / lam[0] - 1))))
    m = L.matching_vector_np(a2, [-e for e in c2], kaps[i] / lam[0])
    cg = np.real(m); cg[13] -= np.real(kaps[i] * Fv[i, 0])
    vals = Pn @ cg
    worst2 = min(worst2, minphi(vals, 8) / max(1e-12, vals.sum()))
    tested += 1
print(f"NUMERIC D2b min of Phi over SO(4) (8 BFGS restarts) at {tested} admissible points "
      f"({min(25, len(dpos))} with d>0): worst min Phi / sum(gaps) = {worst2:.3e}; "
      f"Mobius common-factor rel. error {maxrel:.1e}")
st = L.random_point_on_contact_stratum(T4, seed=7, margin=1e-2)
if st is not None:
    mv = L.matching_vector_np(st["a"], st["b"], st["alpha"])
    cg = np.real(mv); cg[13] -= np.real(L.F_values_np(st["a"], st["b"], st["alpha"])[st["pi_star"]])
    vals = Pn @ cg
    print(f"NUMERIC D3 lib stratum point (seed 7, margin 1e-2): contact {st['contact']}, "
          f"min Phi = {minphi(vals, 20):.3e}")
else:
    print("NUMERIC D3 lib stratum search returned None")

# D4: independent stratum points (lib least-squares search) with d > 0, mapped into the chart
n_d4 = 0; ok_d4 = True; worst4 = np.inf; ratio_spread = 0.0
for seed in range(40):
    st = L.random_point_on_contact_stratum(T4, seed=seed, tries=15, margin=1e-3)
    if st is None or st["contact"] != sorted(T4):
        continue
    gl = st["gaps"]
    if gl[2] - gl[13] <= 1e-6:
        continue
    A_ = list(st["a"]); C_ = [-e for e in st["b"]]
    phi_m = lambda s: (s - A_[1]) * (A_[2] - C_[2]) / ((s - C_[2]) * (A_[2] - A_[1]))
    xx = phi_m(A_[3]); vv = phi_m(C_[1]); tt = phi_m(C_[3]) / vv; zz = phi_m(A_[0]); ww = phi_m(C_[0])
    rr = zz * (ww - vv) * (1 - tt * vv) / (vv * (zz - ww))
    uu = -(xx - tt * vv) * (zz * (1 + vv - ww) - vv) / (vv * (zz - ww))
    reals = max(abs(xx.imag), abs(tt.imag), abs(rr.imag), abs(uu.imag))
    xx, tt, rr, uu = xx.real, tt.real, rr.real, uu.real
    kap = 1j * (-1 / vv.imag) / (vv * (zz - ww))
    Ap = [zz, 0, 1, xx]; Cp = [ww, vv, None, tt * vv]
    Fc = np.array([-np.prod([Ap[i] - Cp[s[i]] for i in range(4) if s[i] != 2]) for s in PERMS])
    gc = (kap * (Fc - Fc[0])).real
    off = [s for s in range(24) if s not in T4]
    rat = gc[off] / gl[off]
    ratio_spread = max(ratio_spread, float((rat.max() - rat.min()) / abs(rat.mean())))
    bb_, dd_ = gc[19], gc[2] - gc[13]
    KK = tt * (xx - 1) * bb_ - xx * (1 - tt) * dd_
    ok_d4 &= reals < 1e-6 and rat.min() > 0 and xx > 1 and 0 < tt < 1 and rr > 0 and uu < 0 and KK > 0
    aa4, cc4 = gc[15], gc[10]
    ok_d4 &= aa4 * bb_ * (cc4 - dd_) - dd_**2 * (aa4 + bb_ + cc4 + dd_) > 0
    mv = L.matching_vector_np(st["a"], st["b"], st["alpha"])
    cg = np.real(mv); cg[13] -= np.real(L.F_values_np(st["a"], st["b"], st["alpha"])[st["pi_star"]])
    vals = Pn @ cg
    worst4 = min(worst4, minphi(vals, 10) / vals.sum())
    n_d4 += 1
print(f"NUMERIC D4 lib stratum points with d>0 (seeds 0..39, margin 1e-3): {n_d4} points; chart gaps are a "
      f"positive multiple of the lib gaps (rel. spread {ratio_spread:.1e}) and x>1, 0<t<1, r>0, u<0, K>0, "
      f"capacity (1) strict all hold: {ok_d4}; worst min Phi / sum(gaps) = {worst4:.3e}")

print(f"SUMMARY exact checks failed: {len(FAIL)} {FAIL}")
sys.exit(1 if FAIL else 0)
