"""Audit replay for independent_moc4_T8_polar_reduction_20261022
(supplement label T8 = paper case T_9 = {0,1,2,3,11,12}; also Section 1 of
independent_moc4_T8_uniform_20261022, which repeats this text).

Own code only; uses the shared library ../lib/moc4lib.py.  Nothing from the supplement is
read or executed.  Prints one line per check; exits nonzero if an EXACT check fails.
Run:  python -I -B independent_moc4_T8_polar_reduction_20261022__check.py
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import moc4lib as L  # noqa: E402
import sympy as sp  # noqa: E402

FAIL = []


def check(name, ok, info=""):
    print(f"{'PASS' if ok else 'FAIL'} EXACT {name} {info}".rstrip())
    if not ok:
        FAIL.append(name)


T = frozenset({0, 1, 2, 3, 11, 12})
z, w, v, x, t, r, u, eta, d, tau, zeta = sp.symbols("z w v x t r u eta d tau zeta")
Pp = 1 - t

# ---------------------------------------------------------------------------------------
# 0. Contact set, generators, coordinates
# ---------------------------------------------------------------------------------------
check("contact set is paper T_9", L.CONTACT_SETS[9] == T)
GENS = L.generators(T)
check("generators(T_9) = f0,f1,f2,f7,f8,f9,f10,f12,f15,f18 with f18 = h(3412,-) the only non-minor",
      GENS == [0, 1, 2, 7, 8, 9, 10, 12, 15, 18]
      and [j for j in GENS if not L.F_IS_MINOR[j]] == [18] and L.f_name(18) == "h3412,-")
BASIS = [1, 2, 7, 8, 9, 12, 15, 18]           # doc's coefficient names K,L,J,Q,C,B,h,A
NAMES = ["K", "L", "J", "Q", "C", "B", "h", "A"]
SY = sp.symbols(NAMES)
FORMS = [sum(sp.Integer(L.F_VALUES[k][s]) * SY[i] for i, k in enumerate(BASIS)) for s in range(24)]
# the 8 basis vectors are independent and span N_T's linear hull: f0 and f10 expressed in them
Kc, Lc, Jc, Qc, Cc, Bc, hc, Ac = SY
def expr_in_basis(j):
    sol = sp.solve([FORMS[s] - L.F_VALUES[j][s] for s in range(24)], SY, dict=True)
    return sol
s0 = expr_in_basis(0); s10 = expr_in_basis(10)
check("f0 = f1+f2+f7+f15-f18 and f10 = f1+f8+f9+f15-f18 (unique)",
      s0 == [{Kc: 1, Lc: 1, Jc: 1, Qc: 0, Cc: 0, Bc: 0, hc: 1, Ac: -1}]
      and s10 == [{Kc: 1, Lc: 0, Jc: 0, Qc: 1, Cc: 1, Bc: 0, hc: 1, Ac: -1}])
check("gap forms: g4=A+h, g6=B, g7=C, g8=A+B+J+Q, g9=Q, g13=A+C+L, g16=2A+C+J+K, g18=h, "
      "g10=A+C+J, g14=A+J+Q, g15=A+L+Q, g17=A+K, g20=J, g21=L, g23=K",
      [FORMS[s] for s in (4, 6, 7, 8, 9, 13, 16, 18, 10, 14, 15, 17, 20, 21, 23)]
      == [Ac + hc, Bc, Cc, Ac + Bc + Jc + Qc, Qc, Ac + Cc + Lc, 2 * Ac + Cc + Jc + Kc, hc,
          Ac + Cc + Jc, Ac + Jc + Qc, Ac + Lc + Qc, Ac + Kc, Jc, Lc, Kc])


def coords_from_pivots(R):
    """(K,L,J,Q,C,B,h,A) from pivotal values R = {4,6,7,8,9,13,16,18}."""
    A_ = R[4] - R[18]; B_ = R[6]; C_ = R[7]; Q_ = R[9]
    J_ = R[8] - A_ - B_ - Q_; L_ = R[13] - A_ - C_; K_ = R[16] - 2 * A_ - C_ - J_; h_ = R[18]
    return dict(zip(SY, [K_, L_, J_, Q_, C_, B_, h_, A_]))


def full_vector(R):
    sub = coords_from_pivots(R)
    return [sp.simplify(f.subs(sub)) for f in FORMS]


# Source-cone theorem: g in N_T and g18 >= g4 (coefficient A <= 0)  =>  minor combination.
# Auditor's proof: with A' = -A >= 0 and any lam in [max(0, A'-min(Q,C)), min(A',L,J)],
# g = lam f0 + (A'-lam) f10 + (K-A') f1 + (L-lam) f2 + (J-lam) f7 + (Q-A'+lam) f8
#     + (C-A'+lam) f9 + B f12 + (h-A') f15, all coefficients >= 0 because
# g17=A+K, g4=A+h, g10, g13, g14, g15 >= 0.  Check the identity exactly:
lam = sp.Symbol("lam"); Ap = -Ac
vec = {0: lam, 10: Ap - lam, 1: Kc - Ap, 2: Lc - lam, 7: Jc - lam, 8: Qc - Ap + lam,
       9: Cc - Ap + lam, 12: Bc, 15: hc - Ap}
ok = all(sp.expand(sum(c * L.F_VALUES[j][s] for j, c in vec.items()) - FORMS[s]) == 0
         for s in range(24))
check("source-cone identity (A<=0 branch is a minor combination; auditor's decomposition)", ok)

# ---------------------------------------------------------------------------------------
# 1. Chart, contacts, curve (1), inverse, pivotal formulas
# ---------------------------------------------------------------------------------------
Avec = [z, 0, 1, x]
Cvec = [w, v, t * v, None]      # None = point at infinity (its factor is a common constant)


def Fsig(s, Av=Avec, Cv=Cvec):
    p = sp.Integer(1)
    for i in range(4):
        j = L.PERMS[s][i]
        if Cv[j] is None:
            continue
        p *= (Av[i] - Cv[j])
    return sp.expand(p)


H = [sp.cancel((Fsig(s) - Fsig(0)) / (v * (z - w))) for s in range(24)]
check("H0..H3 = 0, 1-x, P, 1-tx identically on the chart",
      [sp.expand(H[j]) for j in range(4)] == [0, 1 - x, 1 - t, 1 - t * x])
check("r := P - H12 = z(v-1)(tv-w)/(v(z-w))",
      sp.cancel(Pp - H[12] - z * (v - 1) * (t * v - w) / (v * (z - w))) == 0)

# Auditor's lemma: contacts {0,1,2,3} with 8 distinct points force the two-circle chart.
# Normalise A2=0, A3=1, C4=oo; free A4=al, C2=p, C3=q.  D_s = F_s - F_0 (common factor (Z-W)):
al, p_, q_ = sp.symbols("al p q")
Av2 = [sp.Symbol("Z"), 0, 1, al]; Cv2 = [sp.Symbol("W"), p_, q_, None]
D = [sp.factor((Fsig(s, Av2, Cv2) - Fsig(0, Av2, Cv2)) / (Av2[0] - Cv2[0])) for s in (1, 2, 3)]
check("D1,D2,D3 = p(1-al), p-q, p-q*al", [sp.expand(D[0] - p_ * (1 - al)), sp.expand(D[1] - (p_ - q_)),
                                          sp.expand(D[2] - (p_ - q_ * al))] == [0, 0, 0])
lm, mu_ = sp.symbols("lm mu")
# D2/D1 = lm, D3/D1 = mu real; eliminate s=q/p: quadratic in al has root al=1
quad = sp.expand((1 - (1 - lm * (1 - al)) * al) - mu_ * (1 - al))
check("eliminant -lm al^2 + (lm-1+mu) al + (1-mu) has root al=1, other root -(1-mu)/lm (real)",
      sp.expand(quad.subs(al, 1)) == 0 and sp.simplify(quad.subs(al, -(1 - mu_) / lm)) == 0)

a_ = x - r + u + r * (x - 1) / (v - 1) + u / (t * v - 1)
b_ = r * (t * v - x) / (v - 1)
c_ = x + u * v / (t * v - 1)
winv = t * v + r * (t * v - x) / d
zinv = (t * v * d + r * (t * v - x)) / (d + (v - 1) * (t * v - x) / v)
check("inverse: H12 = P - r identically", sp.cancel(H[12].subs({w: winv, z: zinv}) - (Pp - r)) == 0)
check("inverse: H18 - H4 = d identically", sp.cancel((H[18] - H[4]).subs({w: winv, z: zinv}) - d) == 0)
eta_d = sp.cancel((H[7] - (1 - x)).subs({w: winv, z: zinv}))
check("transfer d = (v-1)(eta+b)/(P v)", sp.cancel(d - (v - 1) * (eta_d + b_) / (Pp * v)) == 0)
rr = Pp - H[12]; dd = H[18] - H[4]
check("converse: on any chart point, d := H18-H4 gives w and z by the inverse formulas",
      sp.cancel(dd * (w - t * v) - rr * (t * v - x)) == 0
      and sp.cancel(z * (dd + (v - 1) * (t * v - x) / v) - (t * v * dd + rr * (t * v - x))) == 0)
num11 = sp.numer(sp.together(sp.cancel((H[11] - u).subs({w: winv, z: zinv}).subs(d, (v - 1) * (eta + b_) / (Pp * v)))))
quadE = sp.numer(sp.together(eta ** 2 - a_ * eta - b_ * c_))
ratio = sp.factor(sp.cancel(num11 / quadE))
check("H11 = u  <=>  eta^2 - a eta - b c = 0 (numerator is quadratic times eta-free factor)",
      not ratio.has(eta), f"factor={ratio}")

QPOLY = sp.Poly(sp.expand(quadE), eta)


def reduce_on_curve(e):
    """value of a rational function of (z,w,v) on the curve, as al + be*eta."""
    e = sp.cancel(sp.cancel(e.subs({w: winv, z: zinv})).subs(d, (v - 1) * (eta + b_) / (Pp * v)))
    n, dn = sp.fraction(sp.together(e))
    nr = sp.Poly(n, eta).rem(QPOLY); dr = sp.Poly(dn, eta).rem(QPOLY)
    n1, n0 = [nr.coeff_monomial(eta ** k) for k in (1, 0)]
    d1, d0 = [dr.coeff_monomial(eta ** k) for k in (1, 0)]
    pr = sp.Poly(sp.expand((n1 * eta + n0) * (d1 * (a_ - eta) + d0)), eta).rem(
        sp.Poly(eta ** 2 - a_ * eta - b_ * c_, eta))
    nrm = sp.cancel(-d1 ** 2 * b_ * c_ + d0 * d1 * a_ + d0 ** 2)
    return sp.factor(sp.cancel(pr.coeff_monomial(1) / nrm)), sp.factor(sp.cancel(pr.coeff_monomial(eta) / nrm))


dexpr = (v - 1) * (eta + b_) / (Pp * v)
hh = u / (t * v - 1)
H6c = (t * v - 1) * eta / (t * v - x)
CLAIM = {4: x / v, 6: H6c, 7: eta, 8: H6c + (hh - eta) / x, 9: hh, 13: r * (x - 1) / (v - 1),
         16: -Pp * dexpr + (1 - 1 / x) * eta + (-1 + 1 / x + t * (v - 1)) * hh, 18: x / v + dexpr}
CONST = {}
for j, cl in CLAIM.items():
    al_, be_ = reduce_on_curve(H[j] - cl)
    CONST[j] = al_
    check(f"pivotal H{j} = displayed + real constant on the curve", be_ == 0 and not al_.has(v),
          f"const={al_}")
for j, val in zip([0, 1, 2, 3, 11, 12], [0, 1 - x, 1 - t, 1 - t * x, u, 1 - t - r]):
    al_, be_ = reduce_on_curve(H[j])
    check(f"contact value H{j} constant on the curve", be_ == 0 and sp.simplify(al_ - val) == 0)
check("H4 = 1-t-x at v = oo (constant term)", sp.simplify(CONST[4] - (1 - t - x)) == 0)

# ---------------------------------------------------------------------------------------
# 2. Punctures and residue table
# ---------------------------------------------------------------------------------------
# Poles of the pivotal formulas: denominators v, v-1, tv-1, tv-x, x (constant) and eta,
# whose poles come from a, b c (poles at v=1, 1/t only).
dens = set()
for cl in list(CLAIM.values()) + [a_, b_ * c_]:
    for fac in sp.factor_list(sp.denom(sp.together(cl)))[1]:
        dens.add(sp.factor(fac[0]))
check("denominators of pivotal formulas and of a, bc only vanish at v in {0,1,1/t,x/t}",
      all(sp.solve(f, v) in ([0], [1], [1 / t], [x / t], []) for f in dens), str(dens))
# behaviour at v = oo: a, b, c finite
check("a, b, c finite at v=oo", all(sp.limit(e, v, sp.oo) not in (sp.oo, -sp.oo, sp.zoo)
                                     for e in (a_, b_, c_)))


def residues(v0, branch):
    """Residues of the 8 pivotal H at v=v0 on a branch: ('finite', eta0) or ('pole',)."""
    out = {}
    if branch[0] == "pole":
        lead = sp.factor(sp.limit((v - v0) * a_, v, v0))   # eta ~ lead/(v-v0)
        # leading-order equation zeta^2 - (tau a)|0 zeta - (tau^2 bc)|0 = 0 with (tau^2 bc)|0 = 0
        bc0 = sp.limit((v - v0) ** 2 * b_ * c_, v, v0)
        assert bc0 == 0 and lead != 0
        for j, cl in CLAIM.items():
            e = sp.cancel((tau * cl).subs(v, v0 + tau).subs(eta, zeta / tau)).subs(zeta, lead)
            out[j] = sp.factor(sp.limit(e, tau, 0))
        return out, lead
    eta0 = branch[1]
    for j, cl in CLAIM.items():
        e = sp.cancel((tau * cl).subs(v, v0 + tau).subs(eta, eta0))
        out[j] = sp.factor(sp.limit(e, tau, 0))
    return out, eta0


def finite_roots(v0):
    a0 = sp.limit((v - v0) * a_, v, v0); bc0 = sp.limit((v - v0) * b_ * c_, v, v0)
    if a0 == 0 and bc0 == 0:
        return sp.solve(sp.cancel(eta ** 2 - a_ * eta - b_ * c_).subs(v, v0), eta)
    # one pole branch: finite root = -(bc)/a leading
    return [sp.cancel(-bc0 / a0)]


R1 = r * (x - 1); U = u / t
TABLE = [
    ("0, eta=-rx", 0, ("finite", -r * x), {4: x, 6: 0, 7: 0, 8: 0, 9: 0, 13: 0, 16: 0, 18: x}, "0"),
    ("0, eta=x", 0, ("finite", x), {4: x, 6: 0, 7: 0, 8: 0, 9: 0, 13: 0, 16: x * (1 + r),
                                     18: -x * (r + t) / Pp}, "R21"),
    ("1, eta finite", 1, ("finite", None), {4: 0, 6: 0, 7: 0, 8: 0, 9: 0, 13: R1, 16: 0, 18: 0}, "0"),
    ("1, eta pole", 1, ("pole",), {4: 0, 6: -Pp * R1 / (t - x), 7: R1, 8: -Pp * R1 / (t - x) - R1 / x,
                                   9: 0, 13: R1, 16: (1 - 1 / x) * R1, 18: 0}, "0"),
    ("1/t, eta finite", 1 / t, ("finite", None), {4: 0, 6: 0, 7: 0, 8: U / x, 9: U, 13: 0,
                                                  16: (1 / x - t) * U, 18: 0}, "0"),
    ("1/t, eta pole", 1 / t, ("pole",), {4: 0, 6: 0, 7: U, 8: 0, 9: U, 13: 0, 16: 0, 18: U}, "R18"),
]
check("v=0 roots of the quadratic are x and -rx",
      set(sp.solve(sp.cancel(eta ** 2 - a_ * eta - b_ * c_).subs(v, 0), eta)) == {x, -r * x})
for name, v0, br, want, Aclaim in TABLE:
    if br[0] == "finite" and br[1] is None:
        br = ("finite", finite_roots(v0)[0])
    got, _ = residues(v0, br)
    same = all(sp.simplify(got[j] - want[j]) == 0 for j in want)
    Rfull = full_vector(got)
    AR = sp.simplify(got[18] - got[4])
    target = {"0": 0, "R21": Rfull[21], "R18": Rfull[18]}[Aclaim]
    check(f"residue row '{name}' and A(R) = {Aclaim}", same and sp.simplify(AR - target) == 0,
          f"R={[got[j] for j in (4,6,7,8,9,13,16,18)]}")
# v = x/t: roots 0 and a(x/t); pole only in H6, H8 (equal residues), on the branch eta -> a(x/t)
ax = sp.cancel(a_.subs(v, x / t))
check("v=x/t: b(x/t)=0 so the roots are 0 and a(x/t)", sp.cancel(b_.subs(v, x / t)) == 0)
got, _ = residues(x / t, ("finite", ax))
Kx = sp.factor((x - 1) * ax / t)
check("residue row 'x/t' = (0,K,0,K,0,0,0,0), K=(x-1)a(x/t)/t, A(R)=0",
      all(sp.simplify(got[j]) == 0 for j in (4, 7, 9, 13, 16, 18))
      and sp.simplify(got[6] - Kx) == 0 and sp.simplify(got[8] - Kx) == 0)
got0, _ = residues(x / t, ("finite", 0))
check("v=x/t, branch eta->0: tau*H_j -> 0 for j != 6,8 (H6 handled by simple zero of eta)",
      all(sp.simplify(got0[j]) == 0 for j in (4, 7, 9, 13, 16, 18)))
# ramified case a(x/t) = c(x/t) = 0
sol = sp.solve([sp.numer(sp.together(ax)), sp.numer(sp.together(c_.subs(v, x / t)))], [r, u], dict=True)
check("a(x/t)=c(x/t)=0  <=>  r = x-t, u = -t(x-1)",
      len(sol) == 1 and sp.expand(sol[0][r] - (x - t)) == 0 and sp.expand(sol[0][u] + t * (x - 1)) == 0, str(sol))
ram = {r: x - t, u: -t * (x - 1)}
check("ramified case: a, b, c each divisible by (tv - x)",
      all(sp.cancel(sp.cancel(e.subs(ram)) / (t * v - x)).is_rational_function(v)
          and sp.factor(sp.numer(sp.together(e.subs(ram)))).as_poly(v).eval(x / t) == 0
          for e in (a_, b_, c_)))
# local parameter at coalescence: F_v = -b'(x/t) c(x/t)
Fq = eta ** 2 - a_ * eta - b_ * c_
Fv = sp.diff(Fq, v).subs(eta, 0)
check("coalescent root (a(x/t)=0): F_v = -b'(x/t) c(x/t)",
      sp.simplify((Fv - (-sp.diff(b_, v) * c_)).subs(v, x / t).subs(r, sp.solve(ax, r)[0])) == 0)

# ---------------------------------------------------------------------------------------
# 3. Nonvanishing of d, w, z-w for nonreal v
# ---------------------------------------------------------------------------------------
check("d=0 => eta=-b; quadratic gives b(b+a-c)=0; b+a-c = v(t-1)((rt+u)v-(r+u))/((v-1)(tv-1))",
      sp.cancel((b_ ** 2 + a_ * b_ - b_ * c_) - b_ * (b_ + a_ - c_)) == 0
      and sp.cancel(b_ + a_ - c_ - v * (t - 1) * ((r * t + u) * v - (r + u)) / ((v - 1) * (t * v - 1))) == 0)
check("w=0 => eta=-b/t (via transfer); quadratic = (b/t^2)(b+ta-t^2 c); b+ta-t^2c = x(t-1)(r-tv+t)/(v-1)",
      sp.cancel(dexpr.subs(eta, -b_ / t) - (-r * (t * v - x) / (t * v))) == 0
      and sp.cancel((b_ / t) ** 2 + a_ * b_ / t - b_ * c_ - (b_ / t ** 2) * (b_ + t * a_ - t ** 2 * c_)) == 0
      and sp.cancel(b_ + t * a_ - t ** 2 * c_ - x * (t - 1) * (r - t * v + t) / (v - 1)) == 0)
zmw = sp.factor(sp.cancel(zinv - winv))
check("z - w = (v-1)(x-tv) w/(v d + (v-1)(tv-x)): z=w with w != 0 finite forces (v-1)(tv-x)=0",
      sp.cancel(zinv - winv - (v - 1) * (x - t * v) * winv / (v * d + (v - 1) * (t * v - x))) == 0,
      f"z-w={zmw}")

# ---------------------------------------------------------------------------------------
# 4. Endpoint inventory
# ---------------------------------------------------------------------------------------
nonminor_ok = []
for j in range(24):
    if j in T:
        continue
    if any(not L.F_IS_MINOR[k] for k in L.generators(T | {j})):
        nonminor_ok.append(j)
check("seventh contacts keeping a non-minor generator = {6,7,9,18,20,21,23}",
      nonminor_ok == [6, 7, 9, 18, 20, 21, 23], str(nonminor_ok))
MAPS = [L.g_index_map(g) for g in L.G_ELEMENTS]


def contains_case(S, Tk):
    for gi, tm in enumerate(MAPS):
        if frozenset(tm[i] for i in Tk) <= S:
            return L.G_ELEMENTS[gi]
    return None


for j, k in [(6, 5), (6, 6), (7, 4), (7, 7), (21, 5), (23, 4)]:
    g = contains_case(T | {j}, L.CONTACT_SETS[k])
    check(f"T_9 u {{{j}}} contains a G-image of paper T_{k} (Stage 1)", g is not None, f"g={g}")
for j, k in [(9, 1), (18, 3), (20, 4)]:
    g = L.find_g(L.AUX7[k], T | {j})
    check(f"T_9 u {{{j}}} is a G-image of Table 5 set T^(7)_{k}", g is not None, f"g={g}")
g = contains_case(T | {9}, L.CONTACT_SETS[8])
print(f"INFO T_9 u {{9}} also contains paper T_8 (Stage 4, i.e. later than T_9): g={g}")

print("SUMMARY", "FAILED: " + ", ".join(FAIL) if FAIL else "all EXACT checks passed")
sys.exit(1 if FAIL else 0)
