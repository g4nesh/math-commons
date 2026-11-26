"""Audit replay for root_moc4_T10_B_branch_20261022 (branch B<0 of paper case T_11).

Own code only; nothing from the supplement is read or executed.  Uses ../lib/moc4lib.py.
Chart identities in the spectral parameters (h,t,s,l,mu,v) are tested by exact rational
evaluation at random points with v and conj(v) replaced by INDEPENDENT rationals v, vb
(a randomized identity test: exact arithmetic, not a symbolic proof).  The identities in
(h,t,Y,eta,delta) are proved symbolically.  Exits 1 if an EXACT check fails.
Run:  python -I -B root_moc4_T10_B_branch_20261022__check.py [--numeric]
"""
import os, sys, random
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib"))
from moc4lib import (PERMS, P, F_VALUES, F_IS_MINOR, CONTACT_SETS, generators, find_g, nullspace,
                     values, leibniz_weights_np, so4_from_params)
from fractions import Fraction as Fr
import sympy as sp

FAIL = []
def check(name, ok, level="EXACT"):
    print(("PASS" if ok else "FAIL"), level, name)
    if not ok and level.startswith("EXACT"):
        FAIL.append(name)

T = frozenset({0, 1, 2, 4, 7, 9})
FV = lambda j: F_VALUES[j]
g_ = find_g(CONTACT_SETS[11], T)
check("T = {0,1,2,4,7,9} is an image of paper T_11 = {0,1,3,5,6,8} (g = %s)" % (g_,), g_ is not None)
gens = generators(T)
check("generators(T) = {0,1,2,3,4,5,7,11,12,21,22,26}; non-minor f22, f26",
      gens == [0, 1, 2, 3, 4, 5, 7, 11, 12, 21, 22, 26] and [j for j in gens if not F_IS_MINOR[j]] == [22, 26])
WT = [values(c) for c in nullspace([P[k] for k in sorted(T)])]
on_WT = lambda f: all(f(y) == 0 for y in WT)
A = lambda y: y[16] - y[3]
B = lambda y: -y[3] + y[6] - y[12] + y[13] + y[18]
d = lambda y: -B(y)
terms = [(23, 1), (18, 3), (15, 4), (21, 5), (8, 7), (6, 11), (19, 21)]
check("canonical profile y = y23 f1 + y18 f3 + y15 f4 + y21 f5 + y8 f7 + y6 f11 + y19 f21 + d f22 on W_T",
      on_WT(lambda y: int(any(sum(y[k] * FV(j)[s] for k, j in terms) + d(y) * FV(22)[s] != y[s] for s in range(24)))))
AB = {j: (A(FV(j)), B(FV(j))) for j in gens}
print("INFO (A,B) on generators:", AB)
check("A+B >= 0 on every generator", all(a + b >= 0 for a, b in AB.values()))
check("g19 = G - d and g21 = V - b - d on W_T (G=g3, V=g10, b=g8)",
      on_WT(lambda y: y[19] - y[3] + d(y)) and on_WT(lambda y: y[21] - y[10] + y[8] + d(y)))
check("b = g8 = 0 boundary: {0,1,2,4,7,8} is paper T_7 and lies in T+{8}",
      frozenset({0, 1, 2, 4, 7, 8}) == CONTACT_SETS[7] and {0, 1, 2, 4, 7, 8} <= set(T) | {8})

# ---- spectral chart ---------------------------------------------------------------
def chart(h, t, s, l, v):
    z = h * v / (1 + (h - 1) * v); k = 1 + t * (h - 1)
    r = h * v * ((k - 1) * s * (l - 1) - k) / (((h - 1) * l * s * (k - 1) - h * (k - 1)) * v + l * s * (k - 1) - h)
    c0 = (z * (v - r) - t * r * (v - z)) / ((v - r) - t * (v - z))
    w = v + s * (z - c0) * (v - r) / ((1 - c0) * (z - r))
    return dict(z=z, k=k, r=r, c0=c0, w=w)
def Q(h, t, s, l, v):
    ch = chart(h, t, s, l, v)
    Asp = [None, 0, 1, ch["z"]]; Csp = [ch["w"], v, ch["r"], ch["c0"]]   # A_1 = infinity
    Pp = []
    for p in PERMS:
        e = 1
        for i in (1, 2, 3):
            e *= Asp[i] - Csp[p[i]]
        Pp.append(e)
    K = (ch["z"] - ch["c0"]) * (v - ch["r"])
    return ch, [(Pp[j] - Pp[0]) / K for j in range(24)], Pp
rnd = random.Random(7)
R = lambda: Fr(rnd.randint(-40, 40), rnd.randint(1, 17))
results = {}
n = 0
while n < 12:
    h, t, s, l, mu, v, vb = [R() for _ in range(7)]
    try:
        ch, Q1, Pp = Q(h, t, s, l, v); _, Q2, _ = Q(h, t, s, l, vb)
        g = [mu * (Q1[j] - Q2[j]) / 2 for j in range(24)]   # gap / i  (gap = Re(i mu Q))
        rho = -mu * (v - vb) / 2                            # rho / i  (rho = mu Im v)
        k = ch["k"]; D = h - k; C = h * (k - 1) + l * D; J0 = l * D - h; N0 = k - s * (k - 1) * (l - 1)
        M1 = lambda vv: (D + s * (h - 1) * (k - 1)) * vv + s * (k - 1) - D
        M2 = lambda vv: C * vv - J0
        a, b, c = g[15], g[8], g[18]; dd = -(-g[3] + g[6] - g[12] + g[13] + g[18])
        G, V, q = g[3], g[10], g[23]
        eta, delta, Y = b / a, dd / a, h * (h - 1) / C
        Pq = D * (h + s * (k - 1) * (h - l)) / (s * C * (k - 1)); NP = h + s * (k - 1) * (h - l)
        p = 1 - h; T0 = t - 1; Z = 1 - t * Y; E = p + k * Y; J = 1 + delta + Y * eta; Dn = J - Y
        W = T0 * p - k * delta - k * Y * eta; B0 = h * eta + k * delta
        xi = lambda vv: p * vv - E
        Rr = xi(v) * xi(vb); Rexi = (xi(v) + xi(vb)) / 2
        SG = h * Y * W * Dn / (Z * delta); SV = h * p * T0 * Y / (Z * B0)
        LG = (J**2 * Rr + Y**2 * W**2 - SG) / (2 * J * Y * W); LV = (SV - Rr - k**2 * Y**2) / (2 * k * Y)
        z = ch["z"]; c0 = ch["c0"]; r = ch["r"]; w = ch["w"]
        chk = {
            "D2 = (z-c0)(v-r), i.e. Q2 = 1 (contact 2 fixes kappa = i mu/((z-c0)(v-r)))": Q1[2] == 1,
            "contacts: g_j = 0 for j in T": all(g[j] == 0 for j in T),
            "h, t, s, l ratio formulas reproduce the parameters":
                h == z * (v - 1) / (v * (z - 1)) and t == (z - c0) * (v - r) / ((r - c0) * (v - z))
                and s == (w - v) * (1 - c0) * (z - r) / ((z - c0) * (v - r)) and l == (w - r) * z / ((w - v) * (z - r)),
            "Q1 = -1/(k-1) (contact 1)": Q1[1] == -1 / (k - 1),
            "(1) c = (t-1)(a+d)": c == (t - 1) * (a + dd),
            "(1) g19 = G - d, g21 = V - b - d": g[19] == G - dd and g[21] == V - b - dd,
            "(2) q = (h-k)a - kd - hb": q == (h - k) * a - k * dd - h * b,
            "(2) J0 g6 = -h b": J0 * g[6] == -h * b,
            "(3) a = -rho l s (h-1)/h": a == -rho * l * s * (h - 1) / h,
            "(3) G = rho h s D N0/|M1|^2": G == rho * h * s * D * N0 / (M1(v) * M1(vb)),
            "(3) b = -rho h J0 N0/|M2|^2": b == -rho * h * J0 * N0 / (M2(v) * M2(vb)),
            "(4) kV = -rho h s/|v-1|^2 - D b": k * V == -rho * h * s / ((v - 1) * (vb - 1)) - D * b,
            "r numerator = -h v N0": r * (((h - 1) * l * s * (k - 1) - h * (k - 1)) * v + l * s * (k - 1) - h) == -h * v * N0,
            "delta = P - Y eta": delta == Pq - Y * eta,
            "(l-1)NP + (h-l)N0 = C": (l - 1) * NP + (h - l) * N0 == C,
            "C = h k + J0": C == h * k + J0,
            "E = pZ + Y": E == p * Z + Y,
            "W = -N0 Dn": W == -N0 * Dn,
            "R = |xi|^2 = tYEW/(Z eta)": Rr == t * Y * E * W / (Z * eta),
            "G/a = hYW Dn/(Z|J xi - YW|^2)": G / a == h * Y * W * Dn / (Z * (J * xi(v) - Y * W) * (J * xi(vb) - Y * W)),
            "k g21/a = h^2 p/((-l)|xi+kY|^2) - B0": k * g[21] / a == h * h * p / ((-l) * (xi(v) + k * Y) * (xi(vb) + k * Y)) - B0,
            "g19/a = 2 delta J Y W/|J xi - YW|^2 (Re xi - LG)": g[19] / a == 2 * delta * J * Y * W / ((J * xi(v) - Y * W) * (J * xi(vb) - Y * W)) * (Rexi - LG),
            "g21/a = 2 Y B0/|xi + kY|^2 (LV - Re xi)": g[21] / a == 2 * Y * B0 / ((xi(v) + k * Y) * (xi(vb) + k * Y)) * (LV - Rexi),
        }
    except ZeroDivisionError:
        continue
    n += 1
    for kk, vv in chk.items():
        results.setdefault(kk, []).append(vv)
for kk, vv in results.items():
    check("%s  [exact at %d random rational points, v and vb independent]" % (kk, len(vv)), all(vv), "EXACT-pointwise")

# ---- symbolic identities in (h,t,Y,eta,delta) --------------------------------------
h, t, Y, eta, delta = sp.symbols("h t Y eta delta")
k = 1 + t * (h - 1); p = 1 - h; T0 = t - 1; Z = 1 - t * Y; E = p + k * Y; J = 1 + delta + Y * eta; Dn = J - Y
W = T0 * p - k * delta - k * Y * eta; B0 = h * eta + k * delta
R_ = t * Y * E * W / (Z * eta)
SG = h * Y * W * Dn / (Z * delta); SV = h * p * T0 * Y / (Z * B0)
LG = (J**2 * R_ + Y**2 * W**2 - SG) / (2 * J * Y * W); LV = (SV - R_ - k**2 * Y**2) / (2 * k * Y)
den = sp.factor(sp.denom(sp.together(sp.cancel(LV))))
check("LV has no k pole: reduced denominator = 2 eta (tY-1) B0", sp.simplify(den / (2 * eta * (t * Y - 1) * B0)) in (1, -1))
N = sp.cancel(-(LV - LG) * 2 * eta * delta * Z * J * B0 / h)
check("N(h,Y) := -(LV-LG) 2 eta delta Z J B0/h is a polynomial, of degree 1 in h and 1 in Y",
      sp.denom(sp.together(N)) == 1 and sp.Poly(N, h).degree() == 1 and sp.Poly(N, Y).degree() == 1)
Nh = sp.diff(N, h)
hs = (t - 1) * (1 + delta) / (t - 1 + eta + t * delta); K = eta * (t - 1 - delta) - t * delta**2
check("N_h(h,0) = -(eta+t delta)(eta+delta)(delta+1)", sp.expand(Nh.subs(Y, 0) + (eta + t * delta) * (eta + delta) * (delta + 1)) == 0)
check("N_h(h,1/t) = -eta(eta+t delta)(eta+t(delta+1)-1)/t", sp.simplify(Nh.subs(Y, 1 / t) + eta * (eta + t * delta) * (eta + t * (delta + 1) - 1) / t) == 0)
check("N(h*,0) = -(eta+delta)(delta+1)K/(eta+t delta+t-1)", sp.simplify(N.subs({Y: 0, h: hs}) + (eta + delta) * (delta + 1) * K / (eta + t * delta + t - 1)) == 0)
check("N(h*,1/t) = (eta+delta)(t delta - K)/t", sp.simplify(N.subs({Y: 1 / t, h: hs}) - (eta + delta) * (t * delta - K) / t) == 0)
check("q >= 0 <=> h <= h*: q/a = (t-1)(1+delta) - h(t-1+eta+t delta)",
      sp.expand((h - k) - k * delta - h * eta - ((t - 1) * (1 + delta) - h * (t - 1 + eta + t * delta))) == 0)
a, b, c, d_ = sp.symbols("a b c d")
check("ab(c-d) - d^2(a+b+c+d) = (a+d){b((t-1)a-d) - t d^2} when c = (t-1)(a+d)",
      sp.expand((a * b * (c - d_) - d_**2 * (a + b + c + d_)).subs(c, (t - 1) * (a + d_)) - (a + d_) * (b * ((t - 1) * a - d_) - t * d_**2)) == 0)

if "--numeric" in sys.argv:
    import numpy as np
    from scipy.optimize import minimize
    rng = np.random.default_rng(11)
    def minphi(y, restarts=6):
        y = np.asarray(y, float); best = np.inf
        for _ in range(restarts):
            rr = minimize(lambda th: float(y @ leibniz_weights_np(so4_from_params(th))), rng.uniform(-np.pi, np.pi, 6), method="BFGS")
            best = min(best, rr.fun)
        return best
    F = lambda j: np.array(FV(j), float)
    worst = np.inf
    for it in range(120):
        a_, b_ = rng.exponential(1, 2); dd = rng.uniform(0, 0.999) * np.sqrt(a_ * b_)
        cc = (a_ * b_ * dd + dd * dd * (a_ + b_ + dd)) / (a_ * b_ - dd * dd)
        if it % 3 == 1:
            cc *= 1 + rng.exponential(0.2)
        worst = min(worst, minphi(a_ * F(4) + b_ * F(7) + cc * F(3) + dd * F(22)) / (a_ + b_ + cc + dd))
    print("NUMERIC cited four-coefficient lemma (seed 11): 120 samples, 2/3 on the boundary ab(c-d)=d^2(a+b+c+d); worst normalized min Phi = %.3e" % worst)
    print("NUMERIC control outside the region (a=b=1, c=d=0.9): min Phi = %.4f" % minphi(F(4) + F(7) + 0.9 * F(3) + 0.9 * F(22), 10))
    def gaps_np(h, t, s, l, mu, v):
        ch = chart(h, t, s, l, v)
        Asp = [0, 1, ch["z"]]; Csp = [ch["w"], v, ch["r"], ch["c0"]]
        Pp = np.array([np.prod([Asp[i - 1] - Csp[p[i]] for i in (1, 2, 3)]) for p in PERMS])
        return ((1j * mu / ((ch["z"] - ch["c0"]) * (v - ch["r"]))) * (Pp - Pp[0])).real
    found = tried = 0; worstK = worstphi = np.inf
    while tried < 200000 and found < 40:
        tried += 1
        hh = rng.uniform(-3, 3); tt = rng.uniform(-3, 4); ss = rng.normal() * 2; ll = rng.normal() * 2
        vv = complex(rng.normal(), rng.normal()); mu = rng.choice([-1, 1])
        gg = gaps_np(hh, tt, ss, ll, mu, vv); gg = gg / np.abs(gg).max()
        if gg.min() < -1e-12:
            continue
        Bv = -gg[3] + gg[6] - gg[12] + gg[13] + gg[18]
        if Bv >= 0:
            continue
        found += 1
        a_, b_, c_, d2 = gg[15], gg[8], gg[18], -Bv
        worstK = min(worstK, a_ * b_ * (c_ - d2) - d2 * d2 * (a_ + b_ + c_ + d2)); worstphi = min(worstphi, minphi(gg))
    print("NUMERIC chart sampling (seed 11): %d draws, %d supported points with B<0; worst ab(c-d)-d^2(a+b+c+d) = %.3e; worst min Phi = %.3e" % (tried, found, worstK, worstphi))

print("SUMMARY: %d EXACT failures" % len(FAIL))
sys.exit(1 if FAIL else 0)
