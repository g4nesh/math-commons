"""Independent replay of the proof of Case T_3 = {0,1,2,3,7,8}.

Document replayed: supplement, proofs/manin_transfer_20261021/
    independent_moc4_T3_uniform_20261022.pdf  (4 pages; read, nothing executed).

Run:  MOC4_SUPP=/path/to/supplement python replay_T3.py
(MOC4_SUPP is used only to read one JSON data file for the facet-numbering comparison;
without it that single comparison is skipped.)  Exit status 0 iff every check passes.

Exact parts use sympy rational-function fields / polynomial rings and Fractions.
The section marked NUMERIC uses numpy floating point and is a falsification attempt only.
"""
import random
import sys

sys.dont_write_bytecode = True

import replay_T3_common as C
from sympy.polys.domains import QQ
from sympy.polys.fields import field

T3 = (0, 1, 2, 3, 7, 8)
ck = C.Checker("replay_T3: Case T_3 = {0,1,2,3,7,8} (independent_moc4_T3_uniform_20261022.pdf)")

# ---------------------------------------------------------------------------------------
# A. objects
# ---------------------------------------------------------------------------------------
C.crosscheck_numbering(ck)

ck.section("A. contact set and generators (document p. 1)")
ck.check("document's T = {0,1,2,3,7,8} is literally the paper's T_3 (same lexicographic numbering; identity of G)",
         set(T3) == {0, 1, 2, 3, 7, 8})
gens = C.generators(T3)
ck.check("generators vanishing on T: {0,1,2,3,4,5,6,11,15,16,17}", gens == [0, 1, 2, 3, 4, 5, 6, 11, 15, 16, 17], gens)
ck.check("only f17 is not a minor vector", [j for j in gens if not C.IS_MINOR[j]] == [17])
C.verify_contact_cone(ck, T3)
COORD = (4, 6, 9, 10, 12, 13, 16, 18)
ck.check("span of the generators has dimension 8, and (g4,g6,g9,g10,g12,g13,g16,g18) are coordinates on it",
         C.rank([C.FVAL[j] for j in gens]) == 8
         and C.rank([[C.FVAL[j][s] for s in COORD] for j in gens]) == 8)


def coords(vals):
    return [vals[s] for s in COORD]


def A_of(vals):          # the functional (2) of the document
    return -vals[4] + vals[6] + vals[10] - vals[12] + vals[13] + vals[18]


minors = [j for j in gens if C.IS_MINOR[j]]
ck.check("the ten minor generators have rank 8", C.rank([C.FVAL[j] for j in minors]) == 8)
facets = C.cone_facets([coords(C.FVAL[j]) for j in minors])
ck.check("cone of the ten minors (in the 8-dim span) has exactly 12 facets", len(facets) == 12, len(facets))
A_coord = C.primitive([-1, 1, 0, 1, -1, 1, 0, 1])           # A in the coordinates COORD
neg = [l for l in facets if C.dot(l, coords(C.FVAL[17])) < 0]
ck.check("exactly one facet functional is negative on f17, and it is A = -g4+g6+g10-g12+g13+g18",
         len(neg) == 1 and neg[0] == A_coord, neg)
ck.check("each of the other eleven facet functionals is >= 0 on all eleven generators",
         all(C.dot(l, coords(C.FVAL[j])) >= 0 for l in facets if l != A_coord for j in gens))
ck.info("=> for g in N_T (cone of the 11 generators, Lemma 4.2): g in cone(minors) <=> A(g) >= 0.")

ck.section("A. identity (3) on the span (document p. 2)")
coef = {1: lambda g: g[22], 3: lambda g: g[20], 4: lambda g: g[13], 6: lambda g: g[21],
        11: lambda g: g[6], 15: lambda g: g[19], 16: lambda g: g[10], 17: lambda g: -A_of(g)}
ok = True
for k in gens:
    g = C.FVAL[k]
    rec = [sum(coef[j](g) * C.FVAL[j][s] for j in coef) for s in range(24)]
    ok = ok and tuple(rec) == g
ck.check("(3): g = g22 f1 + g20 f3 + g13 f4 + g21 f6 + g6 f11 + g19 f15 + g10 f16 - A f17 on every generator (hence on the span)", ok)
ck.info("coefficients of f3, f4, f6, f15, f16 are gaps (>= 0); the SOC triple is (f1: g22, f11: g6; f17: -A).")

# ---------------------------------------------------------------------------------------
# B. the universal inequality used
# ---------------------------------------------------------------------------------------
ck.section("B. SOC member used: a f1 + b f11 + c f17, ab >= c^2 (document p. 1, eq. (1))")
rng = random.Random(20261125)
g = C.orbit_check(ck, "SOC", (1, 9, 22), (1, 11, 17), expected_orbit=144)
C.covariance_check(ck, "SOC (1,9;22) -> (1,11;17)", g, [(1, 1), (9, 11), (22, 17)], rng)
C.replay_SOC_representative(ck)

# ---------------------------------------------------------------------------------------
# C. the chart
# ---------------------------------------------------------------------------------------
ck.section("C. normalisation and contact factorisations (document p. 2)")
C.infinity_convention_check(ck, "C", 3)


def products(A, Cc):
    one = A[1] * 0 + 1
    F = []
    for s in C.PERMS:
        p = one
        for i in range(4):
            if Cc[s[i]] is None:
                continue
            p = p * (A[i] - Cc[s[i]])
        F.append(p)
    return F


K0, z, w, v, c, x, t = field("z,w,v,c,x,t", QQ)
F = products([z, K0(0), K0(1), x], [w, v, c, None])
D = lambda i, j: F[i] - F[j]
ck.check("D13 = x D02 and v D23 = c D01 identically (complex z, w, v, c, x)",
         D(1, 3) == x * D(0, 2) and v * D(2, 3) == c * D(0, 1))
ck.check("D02 = (z-w)(c-v), D01 = -v(z-w)(1-x): nonzero for distinct points",
         D(0, 2) == (z - w) * (c - v) and D(0, 1) == -v * (z - w) * (1 - x))
Ft = products([z, K0(0), K0(1), x], [w, v, t * v, None])
Dt = lambda i, j: Ft[i] - Ft[j]
ck.check("with c = t v: F_8 - F_2 = -t v (v-w)(z-1) and F_7 - F_1 = z (v-w)(x - t v) "
         "[the document's D28, D17, with the orientation D_ij = F_j - F_i]",
         Dt(2, 8) == t * v * (v - w) * (z - 1) and Dt(1, 7) == -z * (v - w) * (x - t * v))
r_def = (v - w) * (z - 1) / (z - w)
u_def = z * (x - t * v) / (v * (z - 1))
ck.check("(F_2 - F_8)/(v(z-w)) = t r and (F_1 - F_7)/(v(z-w)) = -u r with r, u as in (5)",
         Dt(2, 8) / (v * (z - w)) == t * r_def and Dt(1, 7) / (v * (z - w)) == -u_def * r_def)
ck.check("contact ratios R_sigma = (F_sigma - F_0)/(v(z-w)) on T = (0,1,2,3,7,8): 0, 1-x, 1-t, 1-tx, (1-x) + u r, (1-t) - t r",
         [(Ft[s] - Ft[0]) / (v * (z - w)) for s in T3]
         == [K0(0), 1 - x, 1 - t, 1 - t * x, (1 - x) + u_def * r_def, (1 - t) - t * r_def])
ck.info("PROOF-READ: kappa = i m/(v(z-w)), m in C*: contacts <=> m(1-x), m(1-t), m(1-tx) real; differences give m x(1-t), m t(1-x) real,")
ck.info("  so x, t real (t != 1, x != 1), then m real; contacts 7, 8 <=> u r and t r real <=> r, u real (t != 0, r != 0).")
ck.check("chart (6) is exhaustive: z((u+t)v - x) = u v and w(z-1-r) = v(z-1) - r z hold identically for r, u defined by (5)",
         z * ((u_def + t) * v - x) == u_def * v and w * (z - 1 - r_def) == v * (z - 1) - r_def * z)
ck.check("if z - 1 - r = 0 then (z-1)(z-v) = 0:  (z-1-r)(z-w) = (z-1)(z-v) identically",
         (z - 1 - r_def) * (z - w) == (z - 1) * (z - v))

K, x, t, u, r, v, vb = field("x,t,u,r,v,vb", QQ)


def chart(x, t, u, r, v):
    z = u * v / ((u + t) * v - x)
    w = (v * (z - 1) - r * z) / (z - 1 - r)
    return z, w


def R_of(x, t, u, r, v):
    z, w = chart(x, t, u, r, v)
    F = products([z, v * 0, v * 0 + 1, x], [w, v, t * v, None])
    return [(f - F[0]) / (v * (z - w)) for f in F]


def gaps_over_rho(x, t, u, r, v, vb):
    """gamma_sigma with g_sigma = Re(kappa (F_sigma - F_0)) = rho * gamma_sigma, where
    kappa = i mu/(v(z-w)), rho = mu Im v, and v, vb stand for v and its conjugate:
    Re(i mu R(v)) = -mu Im R(v) = -rho (R(v) - R(vb))/(v - vb) for R with real coefficients."""
    R1, R2 = R_of(x, t, u, r, v), R_of(x, t, u, r, vb)
    return [-(a - b) / (v - vb) for a, b in zip(R1, R2)]


zc, wc = chart(x, t, u, r, v)
ck.check("(6) inverts (5): r = (v-w)(z-1)/(z-w) and u = z(x-tv)/(v(z-1)) on the chart",
         (v - wc) * (zc - 1) / (zc - wc) == r and zc * (x - t * v) / (v * (zc - 1)) == u)
gam = gaps_over_rho(x, t, u, r, v, vb)
ck.check("on the chart the gaps vanish identically exactly on T", [s for s in range(24) if gam[s] == 0] == list(T3))
N = v * vb
B = (x - t * v) * (x - t * vb)
Dn = ((t + u) * v - (u + x)) * ((t + u) * vb - (u + x))
z2, w2 = chart(t, t, -t, r, v)
ck.check("D = 0 with nonreal v needs t+u = u+x = 0, and then z = w (= v): collision", z2 == w2 and z2 == v)

ck.section("C. identities (7) (document p. 3); g_sigma = rho * gamma_sigma")
ck.check("g4 = rho x/N", gam[4] == x / N)
ck.check("g6 = rho r t u (x-1)/B", gam[6] == r * t * u * (x - 1) / B)
ck.check("g9 = -x g6", gam[9] == -x * gam[6])
ck.check("g10 = (r+1) g4", gam[10] == (r + 1) * gam[4])
ck.check("g14 = rho u (t-1)[x - t - r(t+u)]/D", gam[14] == u * (t - 1) * (x - t - r * (t + u)) / Dn)
ck.check("A = g13 + g22/t", A_of(gam) == gam[13] + gam[22] / t)
rec = [sum(coef[j](gam) * C.FVAL[j][s] for j in coef) for s in range(24)]
ck.check("the chart's gap vector satisfies (3) in all 24 entries (so it lies in the span)", rec == gam)

ck.section("C. identities (8)-(11) after x = -X, t = -T, rho = -1 (document pp. 3-4)")
X, T = -x, -t
G = [-y for y in gam]                       # gaps at rho = -1
a_, b_, e_, f_, j_, k_ = G[4], G[6], G[14], G[13], G[20], G[22]
A_ = A_of(G)
ck.check("(8a) k/T = e - j", k_ / T == e_ - j_)
ck.check("(8b) A = f - k/T", A_ == f_ - k_ / T)
ck.check("(8c) g21 = X b + A", G[21] == X * b_ + A_)
ck.check("(9)  e = u(T+1)[(1+r)T - r u - X]/D", e_ == u * (T + 1) * ((1 + r) * T - r * u - X) / Dn)
ck.check("(10a) (u - X) j = -X(1+r) a + u(T-X) e/(T+1)", (u - X) * j_ == -X * (1 + r) * a_ + u * (T - X) * e_ / (T + 1))
ck.check("(10b) (X - T)(f - e) = X(T+1) b + (X+1)(T-u) e", (X - T) * (f_ - e_) == X * (T + 1) * b_ + (X + 1) * (T - u) * e_)
ck.check("a = g4 = X/N and b = g6 = -r T u (X+1)/B at rho = -1",
         a_ == X / N and b_ == -r * T * u * (X + 1) / B)
K2, x2, t2, v2, vb2 = field("x,t,v,vb", QQ)
gam2 = gaps_over_rho(x2, t2, -x2, K2(-1), v2, vb2)
ck.check("(11) on the stratum r = -1, u = X: j = g20 = X^2/(|v|^2 (T - X)) at rho = -1",
         -gam2[20] == (-x2) ** 2 / (v2 * vb2 * ((-t2) - (-x2))))
ck.info("sign bookkeeping at rho = -1, T, X > 0: b = -r u T (X+1)/B > 0 <=> r u < 0, as the document states.")

ck.section("C. the closing inequality (document p. 4)")
ck.info("PROOF-READ: T >= X, alpha = -A > 0: alpha = k/T - f <= k/T and alpha = X b - g21 <= X b <= T b, so A^2 <= b k = g6 g22.")
ck.info("Then (3) has coefficients g22 (f1), g6 (f11), -A (f17) with g22 g6 >= A^2: the SOC member (1,11;17) applies.")

# ---------------------------------------------------------------------------------------
# C'. EXACT instances: supporting sources OUTSIDE the minor cone, from finite spectra
# ---------------------------------------------------------------------------------------
ck.section("C'. EXACT instances built from the paper's definitions (finite Gaussian-rational spectra)")
from fractions import Fraction as Fr

GQ = C.GQ


def exact_T3(label, xq, tq, uq, rq, vq):
    zq = uq * vq / ((uq + tq) * vq - xq)
    wq = (vq * (zq - 1) - rq * zq) / (zq - 1 - rq)
    a_s, b_s, alpha, gp = C.exact_paper_gaps([zq, GQ(0), GQ(1), GQ(xq)], [wq, vq, tq * vq, None], 0, 1)
    contact = C.exact_instance_report(ck, label, gp, T3)
    Aq = A_of(gp)
    cf = {j: fn(gp) for j, fn in coef.items()}
    ck.check("%s: A = %s < 0 (not in the minor cone); (3) reproduces all 24 gaps with coefficients >= 0; "
             "SOC condition g22 g6 - A^2 = %s >= 0" % (label, C.fl(Aq), C.fl(cf[1] * cf[11] - Aq * Aq)),
             Aq < 0 and all(val >= 0 for val in cf.values())
             and [sum(cf[j] * C.FVAL[j][s] for j in cf) for s in range(24)] == gp
             and cf[1] * cf[11] >= Aq * Aq)
    Nq, Bq = vq.norm2(), (xq - tq * vq).norm2()
    Dq = ((tq + uq) * vq - (uq + xq)).norm2()
    rho = gp[4] * Nq / xq
    ck.check("%s: with rho := g4 N/x = %s < 0, the formulas (7) for g6, g9, g10, g14 hold at this point" % (label, C.fl(rho)),
             rho < 0 and gp[6] == rho * rq * tq * uq * (xq - 1) / Bq and gp[9] == -xq * gp[6]
             and gp[10] == (rq + 1) * gp[4] and gp[14] == rho * uq * (tq - 1) * (xq - tq - rq * (tq + uq)) / Dq)
    return contact


exact_T3("instance (x,t,u,r,v) = (-1,-3,4,-1/2,i/2)", Fr(-1), Fr(-3), Fr(4), Fr(-1, 2), GQ(0, Fr(1, 2)))
cb = exact_T3("boundary instance (x,t,u,r,v) = (-1,-1,5/2,-1,(-1+i)/2) with T = X and r = -1",
              Fr(-1), Fr(-1), Fr(5, 2), Fr(-1), GQ(Fr(-1, 2), Fr(1, 2)))
ck.check("the boundary instance has extra contacts beyond T (non-strict contacts are exercised)", len(cb) > 6, cb)

# ---------------------------------------------------------------------------------------
# D. NUMERIC falsification attempts (floating point; evidence only)
# ---------------------------------------------------------------------------------------
ck.section("D. NUMERIC falsification attempts on the intermediate claims (numpy, seeds printed)")
import numpy as np

PA = np.array(C.PERMS)


def sample_T3(xs, ts, us, rs, vs, rngn):
    """Gaps of the configuration (A, C) = ((z,0,1,x), (w,v,tv,inf)) computed from FINITE
    spectra: a random Moebius map makes all eight points finite, a_i = m(A_i),
    b_j = -m(C_j), F_sigma = prod_i (a_i + b_sigma(i)) (the paper's z_sigma), and the
    direction kappa is chosen so that Re kappa(F_1 - F_0) = 0."""
    zs = us * vs / ((us + ts) * vs - xs)
    ws = (vs * (zs - 1) - rs * zs) / (zs - 1 - rs)
    n = len(xs)
    al, be, ga, de = [rngn.normal(size=n) + 1j * rngn.normal(size=n) for _ in range(4)]
    mob = lambda s: (al * s + be) / (ga * s + de)
    a = np.stack([mob(zs), mob(0 * zs), mob(0 * zs + 1), mob(xs + 0j)])
    b = -np.stack([mob(ws), mob(vs), mob(ts * vs), al / ga])
    Fs = np.ones((24, n), dtype=complex)
    for k, s in enumerate(C.PERMS):
        for i in range(4):
            Fs[k] *= a[i] + b[s[i]]
    d1 = Fs[1] - Fs[0]
    kap = 1j * np.conj(d1) / np.abs(d1)
    g = (kap * (Fs - Fs[0])).real
    g = g / np.max(np.abs(g), axis=0)
    return g


def run_family(name, seed, n, gen):
    rngn = np.random.default_rng(seed)
    xs, ts, us, rs, vs = gen(rngn, n)
    g = sample_T3(xs, ts, us, rs, vs, rngn)
    tol = 1e-9
    contact_err = np.max(np.abs(g[list(T3)]), axis=0)
    good = np.isfinite(contact_err)
    worst_contact = float(np.max(contact_err[good]))
    # orient: the supporting sign, if any
    stats = dict(n=int(good.sum()), supp=0, t_pos=0, r_pos=0, r_neg=0, A_neg=0)
    viol = []
    worst = dict(tpos_A=np.inf, rpos_A=np.inf, TX=np.inf, soc=np.inf)
    for sign in (1.0, -1.0):
        gg = sign * g
        sup = good & (np.min(gg, axis=0) >= -tol)
        idx = np.nonzero(sup)[0]
        stats["supp"] += len(idx)
        for i in idx:
            h = gg[:, i]
            xv, tv, uv, rv = xs[i], ts[i], us[i], rs[i]
            Av = -h[4] + h[6] + h[10] - h[12] + h[13] + h[18]
            if tv > 0:
                stats["t_pos"] += 1
                worst["tpos_A"] = min(worst["tpos_A"], Av)
                if Av < -1e-8:
                    viol.append(("t>0 but A<0", i))
            else:
                if not (xv < 0 and rv >= -1 - 1e-9 and rv * uv < 0):
                    viol.append(("sign claims x<0, r>=-1, ru<0 fail", i))
                if rv > 0:
                    stats["r_pos"] += 1
                    worst["rpos_A"] = min(worst["rpos_A"], Av)
                    if Av < -1e-8:
                        viol.append(("r>0 but A<0", i))
                else:
                    stats["r_neg"] += 1
                    worst["TX"] = min(worst["TX"], (-tv) - (-xv))
                    if (-tv) < (-xv) - 1e-9:
                        viol.append(("T<X in branch -1<=r<0", i))
            if Av < 0:
                stats["A_neg"] += 1
                m = h[6] * h[22] - Av * Av
                worst["soc"] = min(worst["soc"], m)
                if m < -1e-8:
                    viol.append(("A<0 and A^2 > g6 g22", i))
    ck.check("NUMERIC %s: seed %d, %d samples, max |contact gap| %.1e, %d supporting "
             "(t>0: %d, r>0: %d, -1<=r<0: %d, A<0: %d); no violation of the intermediate claims"
             % (name, seed, stats["n"], worst_contact, stats["supp"], stats["t_pos"], stats["r_pos"],
                stats["r_neg"], stats["A_neg"]),
             not viol and worst_contact < 1e-7, viol[:3])
    ck.info("   min A (t>0) %.3e | min A (r>0) %.3e | min T-X (-1<=r<0) %.3e | min g6 g22 - A^2 (A<0) %.3e"
            % (worst["tpos_A"], worst["rpos_A"], worst["TX"], worst["soc"]))
    return stats


def gen_broad(rngn, n):
    sg = lambda: rngn.choice([-1.0, 1.0], size=n)
    mag = lambda: np.exp(rngn.normal(scale=1.2, size=n))
    return (sg() * mag(), sg() * mag(), sg() * mag(), sg() * mag(),
            rngn.normal(size=n) + 1j * rngn.normal(size=n))


def gen_target(rngn, n):
    """Region of the document's last branch: x, t < 0, -1 <= r < 0, u > 0, T >= X."""
    X = np.exp(rngn.normal(scale=1.0, size=n))
    T = X * (1 + np.exp(rngn.normal(scale=1.5, size=n)) * rngn.integers(0, 2, size=n))   # T = X half the time
    r = -rngn.uniform(0, 1, size=n)
    r[: n // 10] = -1.0                                                                 # boundary r = -1
    u = np.exp(rngn.normal(scale=1.2, size=n))
    u[n // 10: n // 5] = X[n // 10: n // 5]                                             # boundary u = X
    return (-X, -T, u, r, rngn.normal(size=n) + 1j * rngn.normal(size=n))


def gen_rpos(rngn, n):
    X = np.exp(rngn.normal(scale=1.0, size=n))
    T = np.exp(rngn.normal(scale=1.0, size=n))
    r = np.exp(rngn.normal(scale=1.2, size=n))
    u = -np.exp(rngn.normal(scale=1.2, size=n))
    return (-X, -T, u, r, rngn.normal(size=n) + 1j * rngn.normal(size=n))


s1 = run_family("broad chart sample", 31001, 200000, gen_broad)
s2 = run_family("targeted branch t<0, -1<=r<0 (incl. r=-1, u=X, T=X)", 31002, 200000, gen_target)
s3 = run_family("targeted branch t<0, r>0", 31003, 100000, gen_rpos)
ck.check("NUMERIC: the test of (4) is not vacuous (supporting samples with A < 0 were found)",
         s1["A_neg"] + s2["A_neg"] > 0, s1["A_neg"] + s2["A_neg"])

ck.finish()
