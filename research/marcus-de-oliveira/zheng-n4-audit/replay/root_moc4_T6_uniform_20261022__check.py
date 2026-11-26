"""Audit replay for root_moc4_T6_uniform_20261022 (primary proof of paper case
T_7 = {0,1,2,4,7,8}).

Own code only; uses the shared library ../lib/moc4lib.py.  Nothing from the supplement is
read or executed.  Prints one line per check; exits nonzero if an EXACT check fails.
Symbolic identities are rational-function identities in (v, vb = conj v) and the real
parameters, checked with sympy.cancel.
Run:  python -I -B root_moc4_T6_uniform_20261022__check.py
"""
import os
import sys
import random
import itertools
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import moc4lib as L  # noqa: E402
import sympy as sp  # noqa: E402

FAIL = []


def check(name, ok, info="", level="EXACT"):
    print(f"{'PASS' if ok else 'FAIL'} {level} {name} {info}".rstrip())
    sys.stdout.flush()
    if not ok and level == "EXACT":
        FAIL.append(name)


T = frozenset({0, 1, 2, 4, 7, 8})
check("T equals paper T_7", T == L.CONTACT_SETS[7])
gens = L.generators(T)
check("generators of T = {0,1,2,3,4,5,6,11,21,22,23,26}",
      gens == [0, 1, 2, 3, 4, 5, 6, 11, 21, 22, 23, 26], str(gens))
check("nonminor generators exactly 22, 26", [j for j in gens if not L.F_IS_MINOR[j]] == [22, 26])
check("f22 and f26 are positive at permutation 3", L.F_VALUES[22][3] > 0 and L.F_VALUES[26][3] > 0)

NS = L.nullspace([L.P[s] for s in sorted(T)])
Wvals = [L.values(c) for c in NS]
COORD = [3, 6, 9, 10, 12, 13, 16, 18]
check("dim W = 8 and (g3,g6,g9,g10,g12,g13,g16,g18) are coordinates",
      len(NS) == 8 and L.rank([[w_[k] for k in COORD] for w_ in Wvals]) == 8)
q = sp.symbols("q0:8")
Mq = sp.Matrix([[Wvals[a][k] for a in range(8)] for k in COORD])
lam = Mq.inv() * sp.Matrix(q)
GQ = [sp.expand(sum(lam[a] * Wvals[a][s] for a in range(8))) for s in range(24)]
gq = dict(zip(COORD, q))
check("(2) g19 = g6 - g12 + g13 + g18 on W", sp.expand(GQ[19] - (gq[6] - gq[12] + gq[13] + gq[18])) == 0)

minor = [j for j in gens if L.F_IS_MINOR[j]]
Gm = [[Fraction(L.F_VALUES[j][k]) for k in COORD] for j in minor]
check("10 minor generators, rank 8", len(minor) == 10 and L.rank(Gm) == 8)
rays = set()
for sub in itertools.combinations(range(len(minor)), 7):
    rows = [Gm[i] for i in sub]
    if L.rank(rows) != 7:
        continue
    n = L.nullspace(rows)[0]
    vals = [sum(n[k] * Gm[i][k] for k in range(8)) for i in range(len(minor))]
    if all(x <= 0 for x in vals):
        n = [-x for x in n]
    elif not all(x >= 0 for x in vals):
        continue
    m = next(x for x in n if x != 0)
    rays.add(tuple(x / abs(m) for x in n))
check("minor cone has exactly 14 extreme dual rays", len(rays) == 14, f"got {len(rays)}")
Gall = {j: [Fraction(L.F_VALUES[j][k]) for k in COORD] for j in gens}
bad = [r_ for r_ in rays if any(sum(r_[k] * Gall[j][k] for k in range(8)) < 0 for j in gens)]
A_lin = -gq[3] + gq[9] + gq[16]
B_lin = -gq[3] + gq[6] + gq[9] - gq[12] + gq[13] + gq[18]


def as_lin(r_):
    return sp.expand(sum(sp.Rational(r_[k].numerator, r_[k].denominator) * q[k] for k in range(8)))


badlin = [as_lin(r_) for r_ in bad]


def prop(e1, e2):
    rr = sp.simplify(e1 / e2)
    return rr.is_number and rr > 0


check("exactly two dual rays negative on a source generator; they are A, B of (1)",
      len(bad) == 2 and any(prop(b, A_lin) for b in badlin) and any(prop(b, B_lin) for b in badlin),
      str(badlin))

# ---------------- the spectral chart (exact rational-function arithmetic, python-flint) ----
# Every product F_sigma has real (integer) coefficients in v and the real parameters, so with
# kappa = i mu / ((z-c)(v-r)) and X_sigma = (F_sigma - F_0)/((z-c)(v-r)):
#   g_sigma = Re(i mu X_sigma(v)) = (i mu / 2) (X_sigma(v) - X_sigma(w)),  w = conj(v),
#   rho = mu Im v = (i mu / 2) (-(v - w)).
# All identities checked below are homogeneous in the common factor i mu / 2, so we work
# with ghat = X(v) - X(w) and rhohat = -(v - w): exact identities in Q(v, w, r, s, h, t, eta).
import flint  # noqa: E402

CTX = flint.fmpq_mpoly_ctx.get(("v", "w", "r", "s", "h", "t", "eta"), "lex")
_v, _w, _r, _s, _h, _t, _eta = CTX.gens()
ONE = CTX.from_dict({(0,) * 7: 1})


class RF:
    __slots__ = ("n", "d")

    def __init__(self, n, d=None):
        if isinstance(n, RF):
            self.n, self.d = n.n, n.d
            return
        if not isinstance(n, flint.fmpq_mpoly):
            n = ONE * flint.fmpq(n)
        if d is None:
            d = ONE
        g = n.gcd(d)
        if g != ONE and not g.is_zero():
            n = n.exact_division(g) if hasattr(n, "exact_division") else n // g
            d = d // g
        self.n, self.d = n, d

    @staticmethod
    def c(x):
        return x if isinstance(x, RF) else RF(x)

    def __add__(s1, o):
        o = RF.c(o)
        return RF(s1.n * o.d + o.n * s1.d, s1.d * o.d)

    __radd__ = __add__

    def __neg__(s1):
        return RF(-s1.n, s1.d)

    def __sub__(s1, o):
        return s1 + (-RF.c(o))

    def __rsub__(s1, o):
        return RF.c(o) - s1

    def __mul__(s1, o):
        o = RF.c(o)
        return RF(s1.n * o.n, s1.d * o.d)

    __rmul__ = __mul__

    def __truediv__(s1, o):
        o = RF.c(o)
        if o.n.is_zero():
            raise ZeroDivisionError
        return RF(s1.n * o.d, s1.d * o.n)

    def __rtruediv__(s1, o):
        return RF.c(o) / s1

    def __pow__(s1, k):
        out = RF(1)
        for _ in range(k):
            out = out * s1
        return out

    def iszero(s1):
        return s1.n.is_zero()


v, w, r, s, h, t, eta = (RF(x) for x in (_v, _w, _r, _s, _h, _t, _eta))


def eqz(a, b=0):
    return (RF.c(a) - RF.c(b)).iszero()



# Stage 1 (small, sympy): forward direction of the chart, generic C0 = w0
zz, cc, vv_, rr_, w0, hS, tS, sS = sp.symbols("zz cc vv rr w0 hS tS sS")


def products(Al, Cl):
    out = []
    for p in L.PERMS:
        f = 1
        for i in (1, 2, 3):
            f = f * (Al[i] - Cl[p[i]])
        out.append(f)
    return out


Fg = products([None, 0, 1, zz], [w0, vv_, rr_, cc])
hdef = zz * (vv_ - 1) / (vv_ * (zz - 1))
tdef = (cc - zz) * (vv_ - rr_) / ((cc - rr_) * (vv_ - zz))
check("(F0-F1)/(F2-F4) = v(z-1)/[z(v-1)] = 1/h", sp.cancel((Fg[0] - Fg[1]) / (Fg[2] - Fg[4]) - 1 / hdef) == 0)
R021 = sp.cancel((Fg[0] - Fg[2]) / (Fg[0] - Fg[1]))
okt = any(sp.cancel(R021 - cand) == 0 for cand in (tdef * (hdef - 1), -tdef * (hdef - 1),
          1 / (tdef * (hdef - 1)), -1 / (tdef * (hdef - 1)), tdef * (1 / hdef - 1), -tdef * (1 / hdef - 1),
          tdef * (hdef - 1) / hdef, -tdef * (hdef - 1) / hdef))
print("INFO (F0-F2)/(F0-F1) =", sp.factor(R021))
check("(F0-F2)/(F0-F1) is +-t(h-1), +-t(h-1)/h or a reciprocal (so t real once h real)", okt)
check("F0 - F2 = (z-c)(r-v)", sp.expand(Fg[0] - Fg[2] - (zz - cc) * (rr_ - vv_)) == 0)
check("(F0-F8)/(F0-F2) = 1 + r(w-v)/(v-r)  (contact 8 <=> s real)",
      sp.cancel((Fg[0] - Fg[8]) / (Fg[0] - Fg[2]) - 1 - rr_ * (w0 - vv_) / (vv_ - rr_)) == 0)
# contact 7 with w = v + s(v-r)/r, z, c from (h, t): (F0-F7)/(F0-F2) = s X + 1/(k-1),
# X = (1-c)(z-r)/(r(z-c)) = (k/r - 1/v)/(k-1) - 1, k = 1 + t(h-1)
# contact 7, exact rational functions in (v, r, s, h, t) (python-flint; CTX names reused)
rF, sF, hF, tF = RF(_r), RF(_s), RF(_h), RF(_t)
vF = RF(_v)
zH = hF * vF / (1 + (hF - 1) * vF)
cH = (zH * (vF - rF) - tF * rF * (vF - zH)) / ((vF - rF) - tF * (vF - zH))
wS = vF + sF * (vF - rF) / rF
kH = 1 + tF * (hF - 1)
Xr = (1 - cH) * (zH - rF) / (rF * (zH - cH))
check("X = (1-c)(z-r)/[r(z-c)] = (k/r - 1/v)/(k-1) - 1", eqz(Xr, (kH / rF - 1 / vF) / (kH - 1) - 1))
F0 = (0 - vF) * (1 - rF) * (zH - cH)
F2 = (0 - rF) * (1 - vF) * (zH - cH)
F7 = (0 - wS) * (1 - cH) * (zH - rF)
check("(F0-F7)/(F0-F2) = -(s X + 1/(k-1))  (contact 7 <=> X real <=> l real, given s, k real)",
      eqz((F0 - F7) / (F0 - F2), -(sF * Xr + 1 / (kH - 1))))
check("PERMS[7] = (1,0,3,2), PERMS[2] = (0,2,1,3)", L.PERMS[7] == (1, 0, 3, 2) and L.PERMS[2] == (0, 2, 1, 3))

# variables of CTX are reused with new meanings: (v, w=conj v, r->k, s->s, h->h, t->l, eta->spare)
k = r
l = t
spare = eta

# Stage 2: the full chart (3) in (h, k, l, s, v), exact rational functions (python-flint)
def chart_X(vv):
    z = h * vv / (1 + (h - 1) * vv)
    rr = k * vv / (1 + l * vv)
    tt = (k - 1) / (h - 1)
    c = (z * (vv - rr) - tt * rr * (vv - z)) / ((vv - rr) - tt * (vv - z))
    ww = vv + s * (vv - rr) / rr
    Cl = [ww, vv, rr, c]
    Al = [None, RF(0), RF(1), z]
    F = []
    for p in L.PERMS:
        f = RF(1)
        for i in (1, 2, 3):
            f = f * (Al[i] - Cl[p[i]])
        F.append(f)
    den = (z - c) * (vv - rr)
    return [(F[j] - F[0]) / den for j in range(24)], rr, ww, tt


Xv, rv, wv, tt = chart_X(v)
Xw, _, _, _ = chart_X(w)
g = [Xv[j] - Xw[j] for j in range(24)]
check("chart (3): all six contact gaps vanish identically", all(g[j].iszero() for j in T))
check("r(w-v)/(v-r) = s in the chart", eqz(rv * (wv - v) / (v - rv), s))
check("k/r - 1/v = l in the chart", eqz(k / rv - 1 / v, l))
rho = -(v - w)
Dd = h - k
m = h - l - 1
n = k - l - 1
al = h * k - h * l - k
Mv = al * v - Dd
Mw = al * w - Dd
e = g[14]
forms = {3: rho * h * Dd * n * k / ((k - 1) * Mv * Mw), 6: -rho * s / (k * v * w),
         10: rho * h * s * n / ((k - 1) * k * (v - 1) * (w - 1)), 14: -rho * (k + l * s) / k,
         18: -Dd / (k - 1) * e, 9: (1 + s * m / Dd) * g[3]}
for j, f in forms.items():
    check(f"(4)-(5) g{j}", eqz(g[j], f))


def lin_eval(expr, vals):
    out = RF(0)
    for qi, val in zip(q, vals):
        cf = expr.coeff(qi)
        if cf != 0:
            out = out + RF(flint.fmpq(int(sp.numer(cf)), int(sp.denom(cf)))) * val
    return out


qs = [g[j] for j in COORD]
Aval = lin_eval(A_lin, qs)
Bval = lin_eval(B_lin, qs)
check("(2) holds for chart gaps", eqz(g[19], g[6] - g[12] + g[13] + g[18]))
check("chart gaps lie in W (linear reconstruction, all 24)", all(eqz(g[j], lin_eval(GQ[j], qs)) for j in range(24)))
check("(6) A = g16 + (sm/D) g3", eqz(Aval, g[16] + s * m / Dd * g[3]))
check("(6) A = g23 + e/t", eqz(Aval, g[23] + e / tt))
check("(7) B = g19 + (sm/D) g3", eqz(Bval, g[19] + s * m / Dd * g[3]))
check("(7) B = g15 - e/t", eqz(Bval, g[15] - e / tt))
check("G/g10 = D k^2 |v-1|^2/(s |M|^2)", eqz(g[3] / g[10], Dd * k * k * (v - 1) * (w - 1) / (s * Mv * Mw)))
check("g6/e = s/(|v|^2 (k+ls))", eqz(g[6] / e, s / (v * w * (k + l * s))))
check("g10/g6 = -h(k-l-1)|v|^2/((k-1)|v-1|^2)", eqz(g[10] / g[6], -h * n * v * w / ((k - 1) * (v - 1) * (w - 1))))
check("g18/e = -D/(k-1)", eqz(g[18] / e, -Dd / (k - 1)))
# interpolation identity, W a fresh variable (spare), l = h - 1 + W
Wv = spare
lW = h - 1 + Wv
qq = s / Dd
alW = h * k - h * lW - k
check("interpolation (k+ls)/D = (1-zeta)k/D - zeta alpha/(D W), zeta = qW, l = h-1+W",
      eqz((k + lW * s) / Dd, (1 - qq * Wv) * k / Dd - qq * Wv * alW / (Dd * Wv)))
check("alpha/D = 1 - h(1 + W/D) when l = h-1+W", eqz(alW / Dd, 1 - h * (1 + Wv / Dd)))
check("k - l - 1 = -D - W when l = h-1+W", eqz(k - lW - 1, -Dd - Wv))

print(f"SUMMARY failures={len(FAIL)}")
sys.exit(1 if FAIL else 0)
