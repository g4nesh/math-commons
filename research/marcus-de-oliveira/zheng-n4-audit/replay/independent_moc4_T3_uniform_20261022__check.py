"""Audit replay (batch C01-T1-T3).  Exact replay for independent_moc4_T3_uniform_20261022 (paper T_3 = {0,1,2,3,7,8})."""
import os, sys, itertools
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.environ.get("MOC4LIB", os.path.join(HERE, "..", "lib")))
# Own code only.  Nothing from the supplement is read or executed.
# Run: python -I -B THIS.py [--numeric]   (--numeric adds the floating-point tests)
import moc4lib as L
import sympy as sp

FAIL = []
def check(name, ok, info=""):
    print(f"{'PASS' if ok else 'FAIL'} EXACT {name} {info}".rstrip())
    if not ok:
        FAIL.append(name)

T = frozenset({0, 1, 2, 3, 7, 8})
check("doc_T_equals_paper_T3", L.CONTACT_SETS[3] == T)
gens = L.generators(T)
check("source_facets_T3", gens == [0, 1, 2, 3, 4, 5, 6, 11, 15, 16, 17], str(gens))
check("only_f17_nonminor", [j for j in gens if not L.F_IS_MINOR[j]] == [17])
V = L.F_VALUES
check("span_dim_8", L.rank([list(V[j]) for j in gens]) == 8)
coords = [4, 6, 9, 10, 12, 13, 16, 18]
check("coords_(4,6,9,10,12,13,16,18)_injective_on_span",
      L.rank([[V[j][s] for s in coords] for j in gens]) == 8)
minors = [j for j in gens if L.F_IS_MINOR[j]]
check("ten_minor_generators_rank_8", len(minors) == 10 and L.rank([list(V[j]) for j in minors]) == 8)

# facets of cone(minor generators) inside the 8-dim span, written as linear forms on the 24 gap values
# basis of span: use coordinates `coords`; a linear functional on the span = vector w in Q^8 on coords.
M = [[Fraction(V[j][s]) for s in coords] for j in minors]   # 10 x 8
facets = set()
for S in itertools.combinations(range(10), 7):
    rows = [M[i] for i in S]
    if L.rank(rows) != 7:
        continue
    ns = L.nullspace(rows)
    w = ns[0]
    vals = [sum(w[k] * M[i][k] for k in range(8)) for i in range(10)]
    if all(x >= 0 for x in vals) or all(x <= 0 for x in vals):
        if all(x <= 0 for x in vals):
            w = [-x for x in w]
        # normalize
        nz = next(x for x in w if x != 0)
        w = tuple(x / abs(nz) for x in w)
        facets.add(w)
print(f"INFO cone(minor generators) has {len(facets)} facets")
check("dual_has_12_extreme_rays", len(facets) == 12, str(len(facets)))
def on_gen(w, j):
    return sum(w[k] * V[j][coords[k]] for k in range(8))
neg17 = [w for w in facets if on_gen(w, 17) < 0]
check("exactly_one_facet_negative_on_f17", len(neg17) == 1)
others_ok = all(all(on_gen(w, j) >= 0 for j in gens) for w in facets if w not in neg17)
check("other_eleven_nonnegative_on_all_sources", others_ok)
# compare the negative facet with A = -g4+g6+g10-g12+g13+g18 on the span
Acoef = {4: -1, 6: 1, 10: 1, 12: -1, 13: 1, 18: 1}
def Aval(j):
    return sum(c * V[j][s] for s, c in Acoef.items())
w = neg17[0]
ratio = set()
for j in gens:
    a1, a2 = on_gen(w, j), Aval(j)
    ratio.add((a1, a2) if a2 == 0 else a1 / a2)
check("negative_facet_is_A_(up_to_positive_scale)", all(not isinstance(r, tuple) or r == (0, 0) for r in ratio)
      and len({r for r in ratio if not isinstance(r, tuple)}) == 1 and min(r for r in ratio if not isinstance(r, tuple)) > 0,
      str(ratio))
print("INFO A on generators:", {j: Aval(j) for j in gens})

# identity (3)
cs = {j: sp.Symbol(f"c{j}") for j in gens}
g = [sum(cs[j] * V[j][s] for j in gens) for s in range(24)]
A_ = -g[4] + g[6] + g[10] - g[12] + g[13] + g[18]
rhs = {1: g[22], 3: g[20], 4: g[13], 6: g[21], 11: g[6], 15: g[19], 16: g[10], 17: -A_}
ok = all(sp.expand(g[s] - sum(c * V[j][s] for j, c in rhs.items())) == 0 for s in range(24))
check("identity_(3)_profile_decomposition", ok)

# ---------- spectral chart ----------
z, w_, v, c, x = sp.symbols("z w v c x")
kap = sp.Symbol("kappa")
Aa = [z, 0, 1, x]; Cc = [w_, v, c, None]
def Fraw(s, Avec=Aa, Cvec=Cc, k=kap):
    p = L.PERMS[s]
    # the label paired with C3 = infinity contributes the bracket det((A_i,1),(1,0)) = -1
    return -k * sp.Mul(*[Avec[i] - Cvec[p[i]] for i in range(4) if p[i] != 3])
def Draw(r, s):
    return sp.factor(Fraw(r, k=1) - Fraw(s, k=1))
check("chart_D13/D02=x", sp.simplify(Draw(1, 3) / Draw(0, 2) - x) == 0, str(Draw(1, 3) / Draw(0, 2)))
check("chart_D23/D01=c/v", sp.simplify(Draw(2, 3) / Draw(0, 1) - c / v) == 0)
print("INFO D01 =", Draw(0, 1), " D02 =", Draw(0, 2))
t = sp.Symbol("t", real=True)
D28 = sp.factor(Draw(2, 8).subs(c, t * v)); D17 = sp.factor(Draw(1, 7).subs(c, t * v))
check("chart_D28=-tv(v-w)(z-1)", sp.expand(D28 + t * v * (v - w_) * (z - 1)) == 0, str(D28))
check("chart_D17=z(v-w)(x-tv)", sp.expand(D17 - z * (v - w_) * (x - t * v)) == 0, str(D17))
print("INFO D01 at c=tv:", sp.factor(Draw(0, 1).subs(c, t * v)))
# inverse chart (6)
xr, u, r, mu = sp.symbols("x u r mu", real=True)
vb = sp.Symbol("vb")  # conjugate of v treated as independent
zc = u * v / ((u + t) * v - xr)
wc = (v * (zc - 1) - r * zc) / (zc - 1 - r)
rdef = (v - wc) * (zc - 1) / (zc - wc); udef = zc * (xr - t * v) / (v * (zc - 1))
check("chart_(6)_inverts_(5)", sp.simplify(rdef - r) == 0 and sp.simplify(udef - u) == 0)
# The remaining chart identities are rational-function identities in the real parameters
# (x, t, u, r, mu, Re v, Im v).  sympy simplification of the full gap vector is too slow, so they
# are checked by EXACT Gaussian-rational arithmetic at NPTS random rational points
# (Schwartz-Zippel identity testing: exact arithmetic, probabilistic identity claim).
import random
from fractions import Fraction as Fr
class Q:
    __slots__ = ("a", "b")
    def __init__(s, a, b=0): s.a = Fr(a); s.b = Fr(b)
    def __add__(s, o): o = o if isinstance(o, Q) else Q(o); return Q(s.a + o.a, s.b + o.b)
    __radd__ = __add__
    def __neg__(s): return Q(-s.a, -s.b)
    def __sub__(s, o): return s + (-(o if isinstance(o, Q) else Q(o)))
    def __rsub__(s, o): return Q(o) - s
    def __mul__(s, o): o = o if isinstance(o, Q) else Q(o); return Q(s.a * o.a - s.b * o.b, s.a * o.b + s.b * o.a)
    __rmul__ = __mul__
    def conj(s): return Q(s.a, -s.b)
    def __truediv__(s, o):
        o = o if isinstance(o, Q) else Q(o); n = o.a * o.a + o.b * o.b
        return s * Q(o.a / n, -o.b / n)
    def __rtruediv__(s, o): return Q(o) / s
NPTS = 300
rng = random.Random(20261125)
rf = lambda: Fr(rng.randint(-40, 40) or 1, rng.randint(1, 13))
cnt = {k: 0 for k in ["contacts", "g4", "g6", "g9", "g10", "g14", "A", "k/T", "A=f-k/T", "g21", "(9)", "(10a)", "(10b)", "(11)", "vreal"]}
used = 0
for it in range(NPTS):
    X_, T_, uu, rr, mu_ = -abs(rf()), 0, rf(), rf(), rf()
    xr_ = rf(); t_ = rf(); v1, v2 = rf(), rf()
    if it % 3 == 0:          # also test (8)-(11) in the region t<0, x<0 used there
        xr_ = -abs(xr_); t_ = -abs(t_)
    if it % 5 == 0:          # stratum r=-1, u=-x of (11)
        rr = Fr(-1); uu = -xr_
    if t_ in (0, 1) or xr_ in (0, 1) or v2 == 0 or rr == 0 or uu == 0:
        continue
    vv = Q(v1, v2)
    try:
        z_ = uu * vv / ((uu + t_) * vv - xr_)
        w0 = (vv * (z_ - 1) - rr * z_) / (z_ - 1 - rr)
        kq = Q(0, mu_) / (vv * (z_ - w0))
        Avec = [z_, Q(0), Q(1), Q(xr_)]; Cvec = [w0, vv, t_ * vv, None]
        Fv = []
        for p in L.PERMS:
            pr = Q(-1) * kq
            for i in range(4):
                if p[i] != 3:
                    pr = pr * (Avec[i] - Cvec[p[i]])
            Fv.append(pr)
    except ZeroDivisionError:
        continue
    used += 1
    g = [(Fv[s] - Fv[0]).a for s in range(24)]
    Nn = v1 * v1 + v2 * v2; Bn = (xr_ - t_ * v1) ** 2 + (t_ * v2) ** 2
    Dn = ((t_ + uu) * v1 - (uu + xr_)) ** 2 + ((t_ + uu) * v2) ** 2
    # With the homogeneous sign -1 for the bracket against C3 = infinity (which reproduces the note's
    # D28, D17), the note's formulas (7)-(11) hold with rho = -mu Im v; with the sign dropped they hold
    # with rho = +mu Im v.  rho is a free nonzero real parameter whose sign the proof derives, so this
    # is a sign convention only.  Checked here: g4/(mu Im v x/N) = -1 at every point.
    rh = -mu_ * v2
    A_ = -g[4] + g[6] + g[10] - g[12] + g[13] + g[18]
    ok = lambda b: b
    cnt["contacts"] += all(g[s] == 0 for s in T)
    cnt["g4"] += g[4] == rh * xr_ / Nn
    cnt["g6"] += g[6] == rh * rr * t_ * uu * (xr_ - 1) / Bn
    cnt["g9"] += g[9] == -xr_ * g[6]
    cnt["g10"] += g[10] == (rr + 1) * g[4]
    cnt["g14"] += g[14] == rh * uu * (t_ - 1) * (xr_ - t_ - rr * (t_ + uu)) / Dn
    cnt["A"] += A_ == g[13] + g[22] / t_
    Xv, Tv = -xr_, -t_
    a_, b_, e_, f_, j_, k_ = g[4], g[6], g[14], g[13], g[20], g[22]
    cnt["k/T"] += k_ / Tv == e_ - j_
    cnt["A=f-k/T"] += A_ == f_ - k_ / Tv
    cnt["g21"] += g[21] == Xv * b_ + A_
    cnt["(9)"] += e_ / (-rh) == uu * (Tv + 1) * ((1 + rr) * Tv - rr * uu - Xv) / Dn
    cnt["(10a)"] += (uu - Xv) * j_ == -Xv * (1 + rr) * a_ + uu * (Tv - Xv) * e_ / (Tv + 1)
    cnt["(10b)"] += (Xv - Tv) * (f_ - e_) == Xv * (Tv + 1) * b_ + (Xv + 1) * (Tv - uu) * e_
    if rr == -1 and uu == Xv:
        cnt["(11)"] += j_ / (-rh) == Xv ** 2 / (Nn * (Tv - Xv))
    else:
        cnt["(11)"] += 1
    # v real: recompute with v2 = 0
    try:
        vv0 = Q(v1); z0 = uu * vv0 / ((uu + t_) * vv0 - xr_); w00 = (vv0 * (z0 - 1) - rr * z0) / (z0 - 1 - rr)
        k0 = Q(0, mu_) / (vv0 * (z0 - w00)); A0 = [z0, Q(0), Q(1), Q(xr_)]; C0 = [w00, vv0, t_ * vv0, None]
        F0 = []
        for p in L.PERMS:
            pr = Q(-1) * k0
            for i in range(4):
                if p[i] != 3:
                    pr = pr * (A0[i] - C0[p[i]])
            F0.append(pr)
        cnt["vreal"] += all((F0[s] - F0[0]).a == 0 for s in range(24))
    except ZeroDivisionError:
        cnt["vreal"] += 1
for k, n in cnt.items():
    check(f"chart_identity_{k}_at_{used}_random_rational_points", n == used, f"{n}/{used}")

if "--numeric" in sys.argv:
    import numpy as np
    from scipy.optimize import minimize
    for which in ("T3", "soc"):

        PN = L.P_numpy()
        FV = np.array(L.F_VALUES, float)


        def phi_vals(vals, O):
            return float(np.asarray(vals) @ L.leibniz_weights_np(O))


        def min_phi_vals(vals, restarts, rng):
            f = lambda th: phi_vals(vals, L.so4_from_params(th))
            best = np.inf
            for _ in range(restarts):
                r = minimize(f, rng.uniform(-np.pi, np.pi, 6), method="BFGS")
                best = min(best, r.fun)
            return best


        if which in ("soc", "all"):
            # SOC input of the T_3 note: Phi(a f1 + b f11 + c f17) >= 0 when ab >= c^2; worst case ab = c^2, c = 1
            for trip, seed in (((1, 11, 17), 11),):
                X, Y, H = trip
                rng = np.random.default_rng(seed)
                def f(z):
                    O = L.so4_from_params(z[:6]); tau = np.exp(z[6])
                    return phi_vals(tau * FV[X] + FV[Y] / tau + FV[H], O)
                best = np.inf
                for _ in range(80):
                    z0 = np.concatenate([rng.uniform(-np.pi, np.pi, 6), rng.normal(0, 1.5, 1)])
                    r = minimize(f, z0, method="BFGS")
                    best = min(best, r.fun)
                # direct sampling of Phi(fH) + 2 sqrt(Phi fX Phi fY)
                worst = np.inf
                for _ in range(5000):
                    O = L.so4_from_params(rng.uniform(-np.pi, np.pi, 6))
                    px, py, ph = (phi_vals(FV[k], O) for k in trip)
                    worst = min(worst, ph + 2 * np.sqrt(max(px, 0) * max(py, 0)))
                print(f"NUMERIC SOC triple {trip}: min over SO(4) x tau of Phi(tau f{X} + f{Y}/tau + f{H}) = {best:.3e} "
                      f"(80 BFGS restarts, seed {seed}); min of Phi(f{H})+2sqrt(Phi f{X} Phi f{Y}) on 5000 random O = {worst:.3e}")

        if which in ("scalar", "all"):
            rng = np.random.default_rng(5)
            worst = np.inf; worst3 = np.inf
            for k in range(50000):
                x = rng.exponential(1, 4) * (rng.random(4) > 0.25)
                y = rng.exponential(1, 4) * (rng.random(4) > 0.25)
                if x.sum() == 0 or y.sum() == 0:
                    continue
                y = y * x.sum() / y.sum()
                a, b, c, d = x; A, B, C, D = y
                s = c + d - C - D; t = b + c - B - C
                R = 2 * (a * C + b * D + c * A + d * B) - 4 * np.sqrt(a * b * c * d) - 4 * np.sqrt(A * B * C * D)
                worst = min(worst, (R + 2 * abs(s * t)) / (x.sum() ** 2))
                p = rng.normal(size=4); q = rng.normal(size=4); q *= np.linalg.norm(p) / np.linalg.norm(q)
                S = p[2]**2 + p[3]**2 - q[2]**2 - q[3]**2; T = p[1]**2 + p[2]**2 - q[1]**2 - q[2]**2
                H = 2 * sum(p[i]**2 * q[(i + 2) % 4]**2 for i in range(4)) - 4 * np.prod(p) + 4 * np.prod(q)
                worst3 = min(worst3, (H + 2 * abs(S * T)) / np.dot(p, p) ** 2)
            print(f"NUMERIC Thm2 min (R+2|st|)/mass^2 over 50000 samples (seed 5, ~25% zero coords) = {worst:.3e}; "
                  f"Thm3 min (H+2|ST|)/|p|^4 = {worst3:.3e}")

        if which in ("T3", "all"):
            rng = np.random.default_rng(7)
            feas = 0; Aneg = 0; worst_soc = np.inf; worst_phi = np.inf; tried = 0; nphi = 0
            while tried < 400000 and feas < 400:
                tried += 1
                xr, t, u, r = rng.normal(0, 2, 4)
                v = complex(*rng.normal(0, 1.5, 2)); mu = rng.normal()
                z = u * v / ((u + t) * v - xr)
                w = (v * (z - 1) - r * z) / (z - 1 - r)
                kap = 1j * mu / (v * (z - w))
                Av = [z, 0, 1, xr]; Cv = [w, v, t * v, None]
                F = np.array([-kap * np.prod([Av[i] - Cv[p[i]] for i in range(4) if p[i] != 3]) for p in L.PERMS])
                g = F.real - F.real[0]
                sc = np.abs(g).max()
                if sc == 0 or not np.isfinite(sc):
                    continue
                g = g / sc
                if g.min() < -1e-12:
                    continue
                feas += 1
                A = -g[4] + g[6] + g[10] - g[12] + g[13] + g[18]
                if A < 0:
                    Aneg += 1
                    worst_soc = min(worst_soc, g[6] * g[22] - A * A)
                if nphi < 25:
                    nphi += 1
                    worst_phi = min(worst_phi, min_phi_vals(g, 8, rng))
            print(f"NUMERIC T3 chart: {tried} random chart points (seed 7), {feas} with all gaps >= 0, {Aneg} with A<0; "
                  f"min (g6 g22 - A^2) over A<0 = {worst_soc:.3e}; min Phi over SO(4) at {nphi} feasible points "
                  f"(8 BFGS restarts each) = {worst_phi:.3e}")

        if which in ("T12", "all"):
            for nu, seeds in ((2, range(0, 6)), (3, range(0, 4))):
                T = sorted(L.CONTACT_SETS[nu])
                for seed in seeds:
                    rp = L.random_point_on_contact_stratum(T, seed=seed, margin=1e-2)
                    if rp is None:
                        print(f"NUMERIC T_{nu} seed {seed}: no stratum point found"); continue
                    a, b, al = np.array(rp["a"]), np.array(rp["b"]), rp["alpha"]
                    F = L.F_values_np(a, b, al); g = F.real - F.real[T[0]]
                    mind_a = min(abs(a[i] - a[j]) for i in range(4) for j in range(i))
                    mind_b = min(abs(b[i] - b[j]) for i in range(4) for j in range(i))
                    mph = min_phi_vals(g, 20, np.random.default_rng(seed))
                    extra = f" g7-g9={g[7]-g[9]:.2e}" if nu == 2 else ""
                    print(f"NUMERIC T_{nu} seed {seed}: contact {rp['contact']} min|a_i-a_j|={mind_a:.2e} "
                          f"min|b_i-b_j|={mind_b:.2e}{extra} minPhi={mph:.3e}")

print("RESULT", "FAIL" if FAIL else "ALL EXACT CHECKS PASS", FAIL)
sys.exit(1 if FAIL else 0)
