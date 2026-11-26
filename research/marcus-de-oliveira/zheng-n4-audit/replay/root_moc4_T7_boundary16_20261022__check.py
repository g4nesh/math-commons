"""Audit replay for root_moc4_T7_boundary16_20261022 (seven-contact face
{0,2,5,6,7,8,16} of paper case T_8).

Own code only; uses the shared library ../lib/moc4lib.py.  Nothing from the supplement is
read or executed.  Prints one line per check; exits nonzero if an EXACT check fails.
Run:  python -I -B root_moc4_T7_boundary16_20261022__check.py
"""
import os
import sys
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


T = frozenset({0, 2, 5, 6, 7, 8})
check("T equals paper T_8", T == L.CONTACT_SETS[8])
check("generators of T u {16} = {0,2,3,5,6,39,45,46}", L.generators(T | {16}) == [0, 2, 3, 5, 6, 39, 45, 46])
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



# CTX names reused: v -> E, w -> Eb = conj E, r -> r, s -> s, h -> x, t -> t, eta -> spare.
# Gaps (mu = -1): g_j = Im R_j(E) = (R_j(E) - R_j(Eb))/(2i); we use ghat_j = R_j(E) - R_j(Eb).
E, Eb, x = v, w, h


def source_R(EE):
    P0 = EE * EE - ((1 + s) * x - 1 - t) * EE + s * x * (x - 1)
    vv = r - t * ((r - 1) * EE + s * (x - 1)) / P0
    ww = vv * (EE - s * x) / (vv * EE - s)
    zz = (ww * (s + vv) - vv) / (s - 1 + ww)
    Al = [None, vv, x * vv, zz]
    Cl = [r, RF(1), RF(0), ww]
    F = []
    for p in L.PERMS:
        f = RF(-1)          # bracket of A0 = infinity with C_j is -1 (document convention)
        for i in (1, 2, 3):
            f = f * (Al[i] - Cl[p[i]])
        F.append(f)
    Dn = vv * (zz - ww)
    return [(F[j] - F[0]) / Dn for j in range(24)], vv, ww, zz, P0, F, Dn


RE, vE, wE, zE, P0E, FE, DnE = source_R(E)
RB, _, _, _, _, _, _ = source_R(Eb)
gh = [RE[j] - RB[j] for j in range(24)]
bh = gh[1] - gh[3]
dh = gh[12] - gh[16] - gh[18]

# ---------------- boundary j = 16 (s = 1), chart identities ----------------
Rr = -r          # r = -R
Xx = -x          # x = -X
check("j = 16 forces s = 1: g16 = (1-s) g12 identically (exact, general s)", eqz(gh[16], (1 - s) * gh[12]))
Lp = (r - 1) * (r * x - t - 1) + s * (r * x - 1)
Jp = r * (1 - s) * x + r + 2 * s - t - 2
c0 = r * x + s - t - 1
# specialise s = 1 by rebuilding the source curve with s = 1
s_saved = s
s = RF(1)
RE1, vE1, wE1, zE1, _, _, _ = source_R(E)
RB1, _, _, _, _, _, _ = source_R(Eb)
g1h = [RE1[j] - RB1[j] for j in range(24)]
b1 = g1h[1] - g1h[3]
d1 = g1h[12] - g1h[16] - g1h[18]
L1 = Rr * (1 - Rr * Xx) + (Rr + 1) * t
J1 = -Rr - t
Z1E = r * E - (r * x - t)
Z1B = r * Eb - (r * x - t)
check("at s = 1: L = R(1-RX) + (R+1)t, J = -R-t, c0 = rx - t",
      eqz((r - 1) * (r * x - t - 1) + (r * x - 1), L1) and
      eqz(r + 2 - t - 2, J1))
wE_ = Z1E + J1
wB_ = Z1B + J1
check("d/b = t R^2/|w|^2 (w = Z + J)", eqz(d1 / b1, t * Rr * Rr / (wE_ * wB_)))
check("g3/b = (1+X) L R^2/|(R+1)w + R^2(1+X)|^2",
      eqz(g1h[3] / b1, (1 + Xx) * L1 * Rr * Rr / (((Rr + 1) * wE_ + Rr * Rr * (1 + Xx)) * ((Rr + 1) * wB_ + Rr * Rr * (1 + Xx)))))
check("g23 = g9 - d and g9 = R b at s = 1, so g23/b = R - d/b", eqz(g1h[9], Rr * b1))
s = s_saved

# ---------------- the two coefficient intervals (sympy) ----------------
R_, X_, t_, d_, p_ = sp.symbols("R X t d p")      # p = Re w, b = 1
L_ = R_ * (1 - R_ * X_) + (R_ + 1) * t_
w2 = t_ * R_**2 / d_                               # |w|^2 from d = t R^2/|w|^2
# g3 >= d  <=>  (1+X) L R^2 >= d |(R+1) w + R^2 (1+X)|^2
lhs3 = (1 + X_) * L_ * R_**2 - d_ * ((R_ + 1)**2 * w2 + 2 * (R_ + 1) * R_**2 * (1 + X_) * p_ + R_**4 * (1 + X_)**2)
LG = R_ * (1 + X_) * (1 - R_ * X_) / ((R_ + 1) * d_) + t_ * (X_ - R_) / d_ - R_**2 * (1 + X_)**2 / (R_ + 1)
check("g3 >= d  <=>  2(1+X) Re w <= LG  (lhs3 = d R^2 (R+1) [LG - 2(1+X)p])",
      sp.simplify(lhs3 - d_ * R_**2 * (R_ + 1) * (LG - 2 * (1 + X_) * p_)) == 0)
# g10 >= 0 via (6) at s = 1: r(x-1) = R(1+X) > 0 times bracket
r_, x_ = -R_, -X_
Zabs2 = w2 + 2 * (R_ + t_) * p_ + (R_ + t_)**2       # Z = w + R + t
ReZ = p_ + R_ + t_
br = (r_ * x_ - 1) * Zabs2 + 2 * L_ * ReZ + L_ * ((r_ - 1) * (1 + t_) + 1)
LV = R_ * (1 - R_ * X_) / d_ - R_ * X_ * (R_ + 1) - t_ * (X_ - R_)
check("g10 >= 0  <=>  2(1+X) Re w >= LV  (bracket = t R [2(1+X)p - LV])",
      sp.simplify(br - t_ * R_ * (2 * (1 + X_) * p_ - LV)) == 0)
check("LG - LV = (X-R)(1+d)L/((R+1)d)", sp.simplify(LG - LV - (X_ - R_) * (1 + d_) * L_ / ((R_ + 1) * d_)) == 0)

print(f"SUMMARY failures={len(FAIL)}")
sys.exit(1 if FAIL else 0)
