"""Audit replay for independent_moc4_T23_uniform_20261022 (paper case T_24 = {2,3,4,6,7,8}).

Own code only; uses ../lib/moc4lib.py.  Nothing from the supplement is read or executed.
Prints one line per check; exits nonzero if an EXACT check fails.  --numeric adds
floating-point falsification tests.
"""
import os
import sys
import itertools

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import moc4lib as L  # noqa: E402
import sympy as sp  # noqa: E402

FAIL = []


def check(name, ok, info=""):
    print(f"{'PASS' if ok else 'FAIL'} EXACT {name} {info}".rstrip())
    if not ok:
        FAIL.append(name)


T = frozenset({2, 3, 4, 6, 7, 8})
check("contact_set_equals_paper_T24", L.CONTACT_SETS[24] == T)
gens = L.generators(T)
check("generators_T24", gens == [0, 1, 2, 3, 4, 5, 6, 53, 54, 67], str(gens))
check("only_f54_nonminor", [j for j in gens if not L.F_IS_MINOR[j]] == [54])

# canonical decomposition on the subspace {Pc vanishes on T}: test on every generator
COEF = {0: 22, 2: 13, 3: 12, 5: 10, 6: 17, 53: 16, 67: 1}   # f_k coefficient = g_index


def canon(vals):
    out = {k: vals[i] for k, i in COEF.items()}
    out[54] = vals[5] - vals[16]
    return out


ok = True
ok2 = True
for j in gens:
    v = L.F_VALUES[j]
    cf = canon(v)
    rec = [sum(cf[k] * L.F_VALUES[k][s] for k in cf) for s in range(24)]
    ok = ok and rec == list(v)
    ok2 = ok2 and (v[5] - v[16] == v[9] - v[17])
check("canonical_decomposition_identity", ok)
check("delta=g5-g16=g9-g17_on_span", ok2)
print(f"INFO rank of generators = {L.rank([L.F_VECS[j] for j in gens])} (= 14 - 6, so the span is the whole contact subspace)")
dlt = {j: L.F_VALUES[j][5] - L.F_VALUES[j][16] for j in gens}
check("delta_on_generators", dlt == {j: (1 if j == 54 else -1 if j in (1, 4) else 0) for j in gens}, str(dlt))


def eqv(l, r):
    from collections import Counter
    def vec(cnt):
        out = [0] * 14
        for j, c in Counter(cnt).items():
            for k in range(14):
                out[k] += c * L.F_VECS[j][k]
        return out
    return vec(l) == vec(r)


check("f54+f1=f0+f6+f53", eqv([54, 1], [0, 6, 53]))
check("f54+f4=f2+f3+f6+f53", eqv([54, 4], [2, 3, 6, 53]))

# extra ties j with f54(j) = 0
free = [j for j in range(24) if j not in T and L.F_VALUES[54][j] == 0]
check("extra_ties_with_f54_zero", free == [1, 10, 12, 13, 16, 17, 22], str(free))
# extra tie j with f54(j) > 0 forces delta <= 0: delta = -(sum of coefficients that are positive at j)?
# check: for such j, g_j = delta*f54(j) + (nonneg terms) = 0 forces delta <= 0 (all other canonical
# coefficients are >= 0 when delta>0 is assumed).
okj = all(L.F_VALUES[54][j] > 0 for j in range(24) if j not in T and j not in free)
check("all_other_extra_ties_have_f54_positive", okj)

# boundary inputs through the group
claimed = {1: (15, {1, 2, 3, 4, 6, 7}), 10: (20, {2, 3, 4, 6, 8, 10}), 12: (13, {2, 3, 4, 6, 8, 12}),
           22: (12, {2, 3, 4, 6, 8, 22}), 16: (25, {2, 3, 4, 7, 8, 16})}
for j, (k, img) in claimed.items():
    g = L.find_g(L.CONTACT_SETS[k], img)
    check(f"tie_{j}_image_{sorted(img)}_is_g.T_{k}_inside_T+{j}", g is not None and img <= (T | {j}), str(g))
for j in (13, 17):
    S7 = T | {j}
    hits = [k for k, v in L.AUX7.items() if L.find_g(v, S7) is not None]
    sub = sorted({k for k, v in L.CONTACT_SETS.items()
                  for S in itertools.combinations(sorted(S7), 6) if L.find_g(v, S) is not None})
    print(f"INFO T+{j}: equals Table-5 seven-set(s) {hits}; contains six-subsets in orbits of paper cases {sub}")

# ---------------- the rational chart --------------------------------------------------
z, x, c, v, w = sp.symbols("z x c v w")
A = [z, 0, 1, x]
Cc = [c, v, w, None]   # C_3 = infinity: the factor (A_i - inf) is normalized to 1 (a global sign, absorbed in kappa)


def F(k):
    p = L.PERMS[k]
    return sp.Mul(*[(1 if Cc[p[i]] is None else A[i] - Cc[p[i]]) for i in range(4)])


X = 1 - x
D0 = sp.factor(F(3) - F(2))
check("F3-F2=(z-c)wX", sp.simplify(D0 - (z - c) * w * X) == 0, str(D0))
s_ = x * (1 - v) / (w * X)
t_ = c * (z - v) / (w * (z - c))
r_ = (z - 1) * (c - v) / ((z - c) * X)
u_ = (c - w) * (z - v) / (w * (z - c) * X)
q_ = r_ - u_
Hraw = {j: sp.cancel((F(j) - F(2)) / ((z - c) * w * X)) for j in range(24)}
want = {2: 0, 3: 1, 4: s_, 6: r_ - u_, 7: r_ + t_ - u_, 8: r_, 9: r_ + t_ - u_ * X}
for j, e in want.items():
    check(f"H{j}_formula", sp.cancel(Hraw[j] - e) == 0)

# inverse parametrization
S, Tt, R, U, XX = sp.symbols("s t r u X")
Q = R - U
d = 1 - Tt + U * XX
ell = Tt * (1 - XX) + S * XX
h = 1 - U * XX / Tt
n = Tt + Q * XX
zz = 1 - R * XX / d
ww = (1 - XX) * n / ell
cc = ww / h
vv = 1 - S * XX * n / ell
sub = {z: zz, x: 1 - XX, c: cc, v: vv, w: ww}
for nm, e, tgt in [("s", s_, S), ("t", t_, Tt), ("r", r_, R), ("u", u_, U)]:
    check(f"inverse_chart_returns_{nm}", sp.cancel(e.subs(sub, simultaneous=True) - tgt) == 0)
# forward identities used for exhaustiveness
check("w=ch", sp.cancel((c * (1 - u_ * X / t_)) - w) == 0)
check("z-v=th(z-c)", sp.cancel(t_ * (1 - u_ * X / t_) * (z - c) - (z - v)) == 0)
check("r=-d(z-1)/X", sp.cancel(-(1 - t_ + u_ * X) * (z - 1) / X - r_) == 0)

Jp = Q * S * U * XX**2 + (-Q * S * Tt + Q * Tt + S * Tt * U - S * U + Tt * U) * XX - Tt * (R + S * (Tt - 1))
check("z-c=XJ/((t-uX)ell d)", sp.cancel((zz - cc) - XX * Jp / ((Tt - U * XX) * ell * d)) == 0)

print("INFO computing chart forms of H1,H5,H10,H13,H16", flush=True)
Hc = {j: sp.cancel(sp.together(Hraw[j].subs(sub, simultaneous=True))) for j in (1, 5, 10, 13, 16)}
Hc[9] = R + Tt - U * XX
A1 = S - Q - Tt
L1 = Q * ((1 - XX) - 1) - Tt
check("H1=s+A1/L1", sp.cancel(Hc[1] - (S + A1 / L1)) == 0)
check("H10=(r(q+t)-su)/q-(tr/q)A1/L1", sp.cancel(Hc[10] - ((R * (Q + Tt) - S * U) / Q - Tt * R / Q * A1 / L1)) == 0)
H5lim = sp.limit(Hc[5].subs(XX, 1 - sp.Symbol("xx")), sp.Symbol("xx"), sp.oo)
check("H5->1_as_x->inf", sp.simplify(H5lim - 1) == 0, str(H5lim))
# which H_j have J in the denominator
Jpoly = sp.Poly(Jp, XX)
withJ = []
for j in (1, 5, 9, 10, 13, 16):
    den = sp.denom(sp.together(Hc[j]))
    if sp.Poly(den, XX).degree() >= 2 and sp.rem(sp.Poly(den, XX), Jpoly).is_zero:
        withJ.append(j)
print(f"INFO among H1,H5,H9,H10,H13,H16 the ones with J dividing the denominator: {withJ}")
check("H1_H5_H9_H10_have_no_J_pole", not ({1, 5, 9, 10} & set(withJ)))
check("H13_H16_have_J_pole", {13, 16} <= set(withJ))
# Structural pole statement (all 24 H_j): labels have denominators d (z), ell (w, v), ell*h (c);
# 1/((z-c) w X) = (t-uX) ell d / (X J) * ell / ((1-X) n).  Every factor except J is a polynomial
# of degree <= 1 in X with real coefficients, so every nonreal finite pole is a zero of J.
check("label_denominators_real_linear", all(sp.Poly(f, XX).degree() <= 1 for f in (d, ell, h, n, Tt - U * XX, 1 - XX, XX)))

# residue ratio at a root of J: res H13 / res H16 = z/(z-1)
def num_den_wrt_J(H):
    H = sp.together(H)
    N, D = sp.fraction(H)
    K = sp.cancel(D / Jp)
    return N, K


N13, K13 = num_den_wrt_J(Hc[13])
N16, K16 = num_den_wrt_J(Hc[16])
expr = sp.together(N13 * K16 * (zz - 1) - N16 * K13 * zz)
num = sp.numer(expr)
rem = sp.rem(sp.Poly(sp.expand(num), XX), Jpoly)
check("residue_ratio_H13/H16=z/(z-1)_mod_J", rem.is_zero)

# raw numerator of H16 at c=v=z
F16 = sp.factor(F(16) - F(2))
print(f"INFO F16-F2 raw = {F16}")
check("F16-F2_at_c=v=z", sp.factor((F(16) - F(2)).subs({c: z, v: z})) in
      (sp.factor((1 - z) * (x - z) * (z - w)), sp.factor(-(1 - z) * (x - z) * (z - w))),
      str(sp.factor((F(16) - F(2)).subs({c: z, v: z}))))

print(f"SUMMARY exact checks failed: {len(FAIL)} {FAIL}")

if "--numeric" in sys.argv:
    import numpy as np
    # (a) SOC lemma: Phi_{f54}(O) + 2 sqrt(Phi_{f5}(O) Phi_{f2}(O)) >= 0 on SO(4)
    from scipy.optimize import minimize
    Pn = L.P_numpy()
    v54, v5, v2 = (Pn @ np.array([float(t) for t in L.F_VECS[j]]) for j in (54, 5, 2))

    def soc(th):
        wts = L.leibniz_weights_np(L.so4_from_params(th))
        return float(v54 @ wts + 2 * np.sqrt(max(0.0, v5 @ wts) * max(0.0, v2 @ wts)))
    rng = np.random.default_rng(7)
    best = min(minimize(soc, rng.uniform(-np.pi, np.pi, 6), method="Nelder-Mead",
                        options={"maxiter": 4000, "xatol": 1e-10, "fatol": 1e-12}).fun for _ in range(60))
    print(f"NUMERIC SOC_lemma min over SO(4) of Phi54+2sqrt(Phi5 Phi2): {best:.3e} (60 Nelder-Mead starts, seed 7)")
    # (b) chart sampling: real s,t,r,u, complex X; supported points; min Phi
    labf = sp.lambdify((S, Tt, R, U, XX), [zz, cc, vv, ww], "numpy")

    def lam(s0, t0, r0, u0, X0):
        z0, c0, v0, w0 = labf(s0, t0, r0, u0, X0)
        Aa = [z0, 0, 1, 1 - X0]
        Cv = [c0, v0, w0, None]
        Fv = [np.prod([(1 if Cv[pp[i]] is None else Aa[i] - Cv[pp[i]]) for i in range(4)]) for pp in L.PERMS]
        den = (z0 - c0) * w0 * X0
        return [(f - Fv[2]) / den for f in Fv]
    rng = np.random.default_rng(11)
    nsup, worst, worst_fail = 0, np.inf, 0
    for it in range(40000):
        s0, t0, r0, u0 = rng.normal(size=4) * rng.choice([0.3, 1, 3], size=4)
        X0 = complex(rng.normal() * 2, rng.normal() * 2)
        try:
            Hv = np.array(lam(s0, t0, r0, u0, X0), dtype=complex)
        except ZeroDivisionError:
            continue
        for mu in (1.0, -1.0):
            g = -mu * Hv.imag
            if np.all(g >= -1e-12) and np.all(np.isfinite(g)):
                nsup += 1
                dl = g[5] - g[16]
                if dl > 0 and dl * dl > g[10] * g[13]:
                    worst_fail += 1
                    if worst_fail <= 30:
                        cvec = np.linalg.lstsq(Pn, g / np.max(g), rcond=None)[0]
                        val, _ = L.min_phi_over_SO4(cvec, restarts=6, seed=it)
                        worst = min(worst, val)
    print(f"NUMERIC chart sampling: 40000 draws (seed 11), supported={nsup}, "
          f"in SOC-failure region={worst_fail}, worst minPhi over first 30 of them={worst:.3e}")
    # (c) library stratum points
    wl, fl = np.inf, 0
    for seed in range(10):
        pt = L.random_point_on_contact_stratum(T, seed=seed, tries=10, pi_star=2)
        if pt is None:
            continue
        fl += 1
        cvec = np.linalg.lstsq(Pn, pt["gaps"], rcond=None)[0]
        val, _ = L.min_phi_over_SO4(cvec, restarts=8, seed=seed)
        wl = min(wl, val)
    print(f"NUMERIC lib stratum points for T_24: {fl}/10 found, worst minPhi={wl:.3e}")

sys.exit(1 if FAIL else 0)
