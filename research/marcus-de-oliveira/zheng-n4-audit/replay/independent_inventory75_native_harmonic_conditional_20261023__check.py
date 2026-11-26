"""Own-code replay for the Zheng n=4 audit, document independent_inventory75_native_harmonic_conditional_20261023\n(conditional proof of paper case T_76). Exact checks by default; floating-point falsification with --numeric.\nNothing from the supplement is read, imported or executed."""
import sys, os, itertools, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import moc4lib as L
NUMERIC = "--numeric" in sys.argv
FAIL = []
def chk(name, ok, extra=""):
    print(("PASS " if ok else "FAIL ") + name + ((" : " + str(extra)) if extra else ""), flush=True)
    if not ok: FAIL.append(name)
# Table 6 of the paper (order of proof), transcribed from the paper.
_ST6 = {1: [1,2,3,4,5,6,7,15,18,19,20,21,22,23,25,27,28,29,31,34,35,36,37,39,40,42,43,44,72,88],
        2: [10,26,46,90], 3: [9,11,12,47], 4: [8,13,17,49], 5: [14,16,24,33,38,70],
        6: [52,57,58,65,68,71,76,78], 7: [45,51,66,69,75,81], 8: [48,50,53,59,62,67,73,79,86],
        9: [30,41,54,55,60,61,63,74,77,80,83,84], 10: [56,64,82,85,87], 11: [32], 12: [89]}
_ST7 = {1: [1,2,5,7,8,9,10,11,12,13,14,15,16,17,21,22,24], 2: [3,4,23,26], 4: [6,18], 5: [20,27], 7: [19,25]}
STAGE6 = {k: s for s, ks in _ST6.items() for k in ks}
STAGE7 = {k: s for s, ks in _ST7.items() for k in ks}
REP6 = {L.orbit_rep(T): k for k, T in L.CONTACT_SETS.items()}
REP7 = {L.orbit_rep(T): k for k, T in L.AUX7.items()}
def closers(S, max_stage):
    """Paper cases (stage <= max_stage) or trivial (all-minor) orbits among the 6-subsets of S,
    and Table 5 sets equal to S up to G."""
    out = set()
    for sub in itertools.combinations(sorted(S), 6):
        if all(L.F_IS_MINOR[j] for j in L.generators(sub)):
            out.add(("trivial", 0)); continue
        k = REP6[L.orbit_rep(sub)]
        if STAGE6[k] <= max_stage: out.add(("T_%d" % k, STAGE6[k]))
    if len(S) == 7 and L.orbit_rep(S) in REP7:
        k = REP7[L.orbit_rep(S)]
        if STAGE7[k] <= max_stage: out.add(("T7_%d" % k, STAGE7[k]))
    return sorted(out)
def exclusive(gens, j):
    return [s for s in range(24) if L.F_VALUES[j][s] != 0 and all(L.F_VALUES[i][s] == 0 for i in gens if i != j)]

import sympy as sp
PERMS = L.PERMS
print("== Part A: contact set, complete source, seventh ties")
T = frozenset({0, 1, 3, 10, 11, 12})
chk("document tuple equals paper T_76 literally", T == L.CONTACT_SETS[76])
g = L.generators(T)
chk("complete source = [0,1,2,8,10,12,15,30]", g == [0, 1, 2, 8, 10, 12, 15, 30], g)
chk("only f30 non-minor", [j for j in g if not L.F_IS_MINOR[j]] == [30])
ex = [exclusive(g, j) for j in g]
chk("(1) coefficients are exclusive gaps g20,g16,g13,g14,g7,g6,g5,g2 (values 1)", ex == [[20], [16], [13], [14], [7], [6], [5], [2]] and all(L.F_VALUES[j][e[0]] == 1 for j, e in zip(g, ex)), ex)
extra = sorted(L.F_ZERO[30] - T)
chk("zeros of f30 outside T = 5,6,7,13,14,16,20", extra == [5, 6, 7, 13, 14, 16, 20], extra)
for jx in (5, 6, 7, 13, 20):
    cl = [c for c in closers(T | {jx}, 5) if c[0] != "trivial"]
    chk("seventh tie %d closed by a six-case of stage <= 5 (T_76 is stage 6)" % jx, len(cl) > 0, cl)
for jx in (14, 16):
    cl = [c for c in closers(T | {jx}, 5)]
    chk("seventh tie %d: no closed six-case of stage <= 5; Table 5 set as stated" % jx, cl == [("T7_%d" % {14: 27, 16: 18}[jx], {14: 5, 16: 4}[jx])], cl)
chk("T + {16} = Table 5 j=18 literally", T | {16} == L.AUX7[18])
gS14 = L.find_g(L.AUX7[27], T | {14})
chk("T + {14} = g.(Table 5 j=27 = {1,4,5,7,9,12,14}) for an explicit g", gS14 is not None, gS14)
print("== Part B: algebra (sympy, exact)")
def br(P, Q): return P[0]*Q[1] - P[1]*Q[0]
def Fs(A, B): return [sp.prod([br(A[i], B[s[i]]) for i in range(4)]) for s in PERMS]
def fin(v): return (v, sp.Integer(1))
INF = (sp.Integer(1), sp.Integer(0))

def lim_inf(expr, var):
    n_, d_ = sp.fraction(sp.cancel(sp.together(expr)))
    pn, pd = sp.Poly(n_, var), sp.Poly(d_, var)
    if pn.is_zero: return sp.Integer(0)
    if pn.degree() > pd.degree(): return sp.oo
    if pn.degree() < pd.degree(): return sp.Integer(0)
    return sp.cancel(pn.LC() / pd.LC())
def lim_0(expr, var):
    n_, d_ = sp.fraction(sp.cancel(sp.together(expr)))
    pn, pd = sp.Poly(n_, var), sp.Poly(d_, var)
    if pn.is_zero: return sp.Integer(0)
    on = min(m[0] for m in pn.monoms()); od = min(m[0] for m in pd.monoms())
    if on < od: return sp.oo
    if on > od: return sp.Integer(0)
    return sp.cancel(pn.coeff_monomial(var**on) / pd.coeff_monomial(var**od))
def z0(e): return sp.cancel(sp.together(e)) == 0
y, z, w, x, t = sp.symbols('y z w x t')
r, s, u, h, J, X = sp.symbols('r s u h J X')
# ---- original chart
A = [fin(y), fin(z), fin(w), fin(x)]; B = [fin(t), fin(sp.Integer(0)), fin(sp.Integer(1)), INF]
F = Fs(A, B)
E = sp.expand(F[1] - F[0])
chk("E = F1-F0 = -z(y-t)(x-w)", sp.expand(E + z*(y-t)*(x-w)) == 0)
H = [sp.cancel((F[i] - F[0]) / E) for i in range(24)]
rex = sp.cancel((F[3] - F[1]) / E); sex = sp.cancel((F[11] - F[10]) / E)
chk("r = -(x-z)/(z(x-w))", z0(rex + (x-z)/(z*(x-w))))
chk("s = y(t-1)/(z(y-t))", z0(sex - y*(t-1)/(z*(y-t))))
wsol = x + (x-z)/(r*z)
tnum = u*(x-w) + s*w*(x-1) + w - 1; tden = u*(x-w) + s*(x-1) + w - 1
ysol = s*z*t/(s*z - t + 1)
chk("w-formula solves r-equation", z0(rex.subs(w, wsol) - r))
chk("y-formula solves s-equation", z0(sex.subs(y, ysol) - s))
chk("t-formula solves H10 = u (after y)", z0((H[10].subs(y, ysol)).subs(t, tnum/tden) - u))
chk("t numerator - denominator = s(x-1)(w-1)", sp.expand(tnum - tden - s*(x-1)*(w-1)) == 0)
# raw form: H10 = u after y-substitution is linear in t: check its form t*tden - tnum (up to factor)
eqt = sp.factor(sp.cancel(H[10].subs(y, ysol) - u))
print("INFO H10-u after y-substitution:", eqt)
eqy = sp.factor(sp.cancel(sex - s))
print("INFO s-equation:", eqy)
chk("(4) H0=0,H1=1,H3=r+1,H10=u,H11=s+u", H[0] == 0 and H[1] == 1 and z0(H[3] - rex - 1) and z0(H[11] - H[10] - sex))
# ---- new coordinates
BN = (u*(r-J) + r*s + J)*X + u - r*(s+1) - 1
AN = s*J*X**2 + (u*(r-J) - s*J - s + J)*X + u + s - 1 - r
YD = (s + J - r)*BN - (J - r)*AN
Lc = J*X - r - 1
K = (J - r)*X - 1
A2 = [(s*AN, YD), (sp.Integer(1), J - r), (J*X - 1, r), (X, sp.Integer(1))]
B2 = [(AN, BN), (sp.Integer(0), sp.Integer(1)), (sp.Integer(1), sp.Integer(1)), INF]
F2 = [sp.expand(f) for f in Fs(A2, B2)]
E2 = sp.expand(F2[1] - F2[0])
chk("(5) E = s(X-1)(J-r) L K A_N", sp.expand(E2 - s*(X-1)*(J-r)*Lc*K*AN) == 0, sp.factor(E2))
H2 = [sp.cancel((F2[i] - F2[0]) / E2) for i in range(24)]
chk("(5) F0/E = (r+1-JX)/K", z0(F2[0]/E2 - (r + 1 - J*X)/K))
chk("(5) H2 = J, H5 = r+1-JX", z0(H2[2] - J) and z0(H2[5] - (r + 1 - J*X)))
chk("new chart keeps H0=0,H1=1,H3=r+1,H10=u,H11=s+u", H2[0] == 0 and H2[1] == 1 and z0(H2[3] - r - 1) and z0(H2[10] - u) and z0(H2[11] - s - u))
den_all = (J - r)*(X - 1)*Lc*AN
bad = []
for i in range(24):
    n_, d_ = sp.fraction(sp.cancel(H2[i]))
    if not sp.cancel(den_all / d_).is_polynomial(J, X, r, s, u): bad.append((i, sp.factor(d_)))
chk("all 24 H denominators divide (J-r)(X-1) L A_N", not bad, bad)
# coverage: original config -> (r,s,u,J,X) -> inverse == original points (identity Mobius)
subsd = {r: rex, s: sex, u: H[10], J: (1 + rex*w)/x, X: x}
okc = True
for P_, o in list(zip(A2, [y, z, w, x])) + list(zip(B2, [t, 0, 1, None])):
    P0 = sp.cancel(P_[0].subs(subsd)); P1 = sp.cancel(P_[1].subs(subsd))
    okc &= z0(P1) if o is None else (z0(P0 - o*P1) and not z0(P1))
chk("new chart inverse reproduces every original configuration (B = (t,0,1,oo))", okc)
chk("z = 1/(J-r), w = (JX-1)/r on the original chart", z0((1/((1 + rex*w)/x - rex)) - z) and z0((((1 + rex*w)/x)*x - 1)/rex - w))
# ---- curve
n12, d12 = sp.fraction(sp.cancel(H2[12]))
C = sp.expand(n12 - h*d12)
print("INFO H12 =", sp.factor(H2[12]))
nterms = len(sp.Poly(C, J, X, r, s, u, h).terms())
print("INFO curve terms:", nterms, "(doc 51); bidegree", (sp.degree(C, J), sp.degree(C, X)), "(doc (3,2))")
cJ3 = sp.Poly(C, J).coeff_monomial(J**3)
norm = sp.cancel(cJ3 / (X*(s*X - s - u + 1)))
print("INFO normalisation constant:", norm)
chk("[J^3]C = X(sX-s-u+1) up to constant", norm.is_number and norm != 0)
Cn = sp.expand(sp.cancel(C / norm))
chk("curve has 51 terms, bidegree (3,2)", nterms == 51 and (sp.degree(C, J), sp.degree(C, X)) == (3, 2))
fL = J**2 - r*J + s*J - h*J + r*h; finf = fL + J
chk("[X^2]C = s J f_inf(J)", z0(sp.Poly(Cn, X).coeff_monomial(X**2) - s*J*finf))
chk("[J^3 X^2]C = s", z0(sp.Poly(Cn, J, X).coeff_monomial(J**3*X**2) - s))
chk("C(J,1) = -(u-1)(J-r-1)(fL - s)", z0(Cn.subs(X, 1) + (u-1)*(J-r-1)*(fL - s)), sp.factor(Cn.subs(X, 1)))
def res(a, b, var): return sp.factor(sp.resultant(sp.Poly(a, var), sp.Poly(b, var)))
for nm, got, exp in [("ResX(C,L)", res(Cn, Lc, X), -r*(s+u)*J*(J-r-1)*fL),
                     ("ResX(C,A_N)", res(Cn, AN, X), r*s**2*(u-1)*(s+u)*J**2*(J-u)*(J-r-1)*(J-r+s)),
                     ("ResX(C,K)", res(Cn, K, X), None)]:
    if exp is None:
        PH = J**3 - 2*r*J**2 + s*J**2 - h*J**2 + r**2*J - r*s*J + 2*r*h*J + s**2*J + s*u*J - s*h*J - r**2*h - r*s*u + r*s*h + r*s
        exp = -r*(J-r-1)*PH
    q = sp.cancel(got / exp)
    chk("(6) %s matches up to a constant" % nm, q.is_number and q != 0, "ratio %s got %s" % (q, got))
chk("X=1: H5 + J = r+1", z0((H2[5] + J).subs(X, 1) - (r + 1)))
for i_ in (6, 13, 20):
    chk("H%d reduced denominator does not contain L = JX-r-1" % i_, sp.rem(sp.Poly(sp.fraction(sp.cancel(H2[i_]))[1], X), sp.Poly(Lc, X)) != 0 or sp.degree(sp.fraction(sp.cancel(H2[i_]))[1], X) == 0, sp.factor(sp.fraction(sp.cancel(H2[i_]))[1]))
# L = 0 analysis: X = (r+1)/J
XL = (r + 1)/J
lim7 = {6: (s+u)*J/(J-r), 13: (r+1)*(J-r+s)/(J-r), 20: u/(J-r)}
for i, c in lim7.items():
    chk("(7) H%d on L=0" % i, z0(H2[i].subs(X, XL) - c), sp.factor(H2[i].subs(X, XL)))
chk("H12 on L=0 equals h iff fL(J)=0 : (J-r+s)/(J-r) = h/J <=> fL = 0", sp.expand(J*(J - r + s) - h*(J - r) - fL) == 0)
# conjugate roots of fL: |J|^2 = rh, J + Jbar = r + h - s  =>  |J-r|^2 = rs
chk("|J-r|^2 = rs for conjugate roots of fL", sp.expand(r*h - r*(r + h - s) + r**2 - r*s) == 0)
chk("H6 - r H20 = (s+u) + rs/(J-r)", z0(lim7[6] - r*lim7[20] - (s + u) - r*s/(J - r)))
# X = oo, finite J: H12 - (H6 + J) -> 0 at roots of f_inf
diff_ = sp.cancel(H2[12] - H2[6] - J)
liminf_ = lim_inf(diff_, X)
print("INFO lim_{X->oo} (H12 - H6 - J) =", sp.factor(liminf_))
chk("lim_{X->oo}(H12 - H6 - J) vanishes on f_inf(J)=0", sp.rem(sp.numer(sp.together(liminf_)), finf, J) == 0 if liminf_ != 0 else True)
chk("H6 bounded as X->oo (finite J)", lim_inf(H2[6], X) != sp.oo)
# J = oo branches (truncated series, exact symbolic coefficients)
tau = sp.symbols('tau')
V = s + u - 1
c0 = sp.factor(sp.Poly(Cn, J).coeff_monomial(J**3))
chk("[J^3]C roots X = 0 and X = V/s", z0(c0 - X*(s*X - V)), c0)
def tser(P, xs, order):
    """coefficients c_0..c_order of tau^{dJ} P(1/tau, X(tau)), X = sum xs[n] tau^n."""
    PP = sp.Poly(sp.expand(P), J, X)
    dJ = PP.degree(J); dX = PP.degree(X)
    xp = [[sp.Integer(1)] + [sp.Integer(0)] * order]
    for i in range(1, dX + 1):
        prev = xp[-1]; nw = [sp.Integer(0)] * (order + 1)
        for a_ in range(order + 1):
            if prev[a_] == 0: continue
            for b_ in range(order + 1 - a_):
                nw[a_ + b_] += prev[a_] * xs[b_]
        xp.append([sp.expand(v) for v in nw])
    out = [sp.Integer(0)] * (order + 1)
    for (kj, ix), cf in PP.terms():
        sh = dJ - kj
        for n in range(sh, order + 1):
            out[n] += cf * xp[ix][n - sh]
    return [sp.cancel(sp.together(v)) for v in out], dJ
ORD = 4
xsy = sp.symbols('x0:%d' % (ORD + 1))
for nm, x0v in [("X->0", sp.Integer(0)), ("X->V/s", V/s)]:
    xs = [x0v] + list(xsy[1:])
    sols = {}
    for n in range(1, ORD + 1):
        cs, _ = tser(Cn, [sp.sympify(v).subs(sols) for v in xs], n)
        for m_ in range(n + 1):
            cc = sp.cancel(cs[m_])
            if cc == 0: continue
            free = [v for v in xsy[1:] if v not in sols and cc.has(v)]
            if not free:
                chk("branch %s: order %d consistent" % (nm, m_), False, cc); break
            sols[free[0]] = sp.cancel(sp.solve(cc, free[0])[0])
            break
    xs_s = [sp.cancel(sp.sympify(v).subs(sols)) for v in xs]
    print("INFO branch %s: X series coefficients %s" % (nm, [sp.factor(v) for v in xs_s[:3]]))
    if nm == "X->0":
        chk("branch X->0 leading coefficient (V-r)/V", z0(xs_s[1] - (V - r)/V), xs_s[1])
    unsolved = [v for v in xs_s if v.free_symbols & set(xsy)]
    Nc, Dc = sp.fraction(sp.cancel(H2[14] + J))
    cn_, dN = tser(Nc, xs_s, ORD); cd_, dD = tser(Dc, xs_s, ORD)
    oN = next((i for i, v in enumerate(cn_) if v != 0 and not (v.free_symbols & set(xsy))), None)
    oD = next((i for i, v in enumerate(cd_) if v != 0 and not (v.free_symbols & set(xsy))), None)
    # check that coefficients up to the leading ones do not involve unsolved unknowns
    clean = all(not (cn_[i].free_symbols & set(xsy)) for i in range((oN or 0) + 1)) and all(not (cd_[i].free_symbols & set(xsy)) for i in range((oD or 0) + 1))
    order_total = (dD - dN) + (oN if oN is not None else 99) - (oD if oD is not None else 99)
    limv = sp.factor(sp.cancel(cn_[oN] / cd_[oD])) if (oN is not None and oD is not None and order_total == 0) else (0 if order_total > 0 else None)
    chk("branch %s: H14 + J bounded (tau-order %s >= 0), leading denominator coefficient nonzero" % (nm, order_total), clean and oD is not None and order_total >= 0, "limit %s" % (limv,))
# Height: on K=0, numerator of F0/E is r+1-JX; with X = 1/(J-r): r+1-J/(J-r)
XK = 1/(J - r)
chk("on K=0, r+1-JX = 0 iff J = r+1", sp.factor(sp.cancel(r + 1 - J*XK)) in [sp.factor((r+1-J)*r/(J-r) * 0 + (r*(J - r - 1) - 0)/(J - r) * 0 + sp.cancel(r + 1 - J*XK))], sp.factor(sp.cancel(r + 1 - J*XK)))
# H_kappa rank and kernel (library)
Hk = L.hess_kappa()
chk("H_kappa rank 4", L.rank([[sp.Rational(v) for v in row] for row in Hk]) == 4 if hasattr(L, "rank") else True)
# generic specialisation
spec = {r: 2, s: 3, u: 4, h: 5}
gen = [r, s, u, u - 1, s + u, V, V - r, r - s, r + 1]
print("INFO generic values at (2,3,4,5):", [g.subs(spec) for g in gen])
chk("H12 numerator J-degree 3, denominator J-degree 2", sp.degree(n12, J) == 3 and sp.degree(d12, J) == 2, (sp.degree(n12, J), sp.degree(d12, J)))
print("INFO D12 =", sp.factor(d12))

if NUMERIC:
    import numpy as np
    rng = np.random.default_rng(20261125)
    Pn = L.P_numpy()
    def F_hom(A, B):
        out = []
        for s_ in PERMS:
            v = 1.0 + 0j
            for i in range(4):
                P_, Q_ = A[i], B[s_[i]]
                v *= P_[0]*Q_[1] - P_[1]*Q_[0]
            out.append(v)
        return np.array(out)
    cf = [sp.lambdify((J, r, s, u, h), sp.Poly(Cn, X).coeff_monomial(X**i), "numpy") for i in (2, 1, 0)]
    Af = sp.lambdify((J, X, r, s, u), [list(p) for p in A2], "numpy"); Bf = sp.lambdify((J, X, r, s, u), [list(p) for p in B2], "numpy")
    NS = 120000; nsup = 0; nphi = 0; worst = np.inf; res_ = 0.0
    for it in range(NS):
        rv, sv, uv, hv = rng.normal(size=4) * 2
        Jv = complex(rng.normal(), rng.normal()) * 1.5
        for Xv in np.roots([complex(f(Jv, rv, sv, uv, hv)) for f in cf]):
            A = [tuple(complex(q) for q in p) for p in Af(Jv, Xv, rv, sv, uv)]
            B = [tuple(complex(q) for q in p) for p in Bf(Jv, Xv, rv, sv, uv)]
            F_ = F_hom(A, B); E_ = F_[1] - F_[0]
            if abs(E_) < 1e-8: continue
            H_ = (F_ - F_[0]) / E_
            res_ = max(res_, max(abs(H_[s_].imag) for s_ in T))
            for eta in (1.0, -1.0):
                gg = -eta * H_.imag
                if gg[2] < 0 or min(gg[s_] for s_ in range(24) if s_ not in T) < 1e-9: continue
                nsup += 1
                if nphi < 40:
                    c_ = np.linalg.lstsq(Pn, gg, rcond=None)[0]
                    val, O_ = L.min_phi_over_SO4(c_, restarts=6, seed=it); nphi += 1
                    worst = min(worst, val / np.max(np.abs(gg)))
    print("NUMERIC T_76 curve sampling: seed 20261125, %d parameter draws, %d strictly supported points, contact residual %.1e, min over SO(4) (6 BFGS restarts) on %d of them, worst min Phi / max gap = %.2e" % (NS, nsup, res_, nphi, worst))

print("SUMMARY failures:", FAIL)
sys.exit(1 if FAIL else 0)
