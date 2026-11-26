"""Audit replay for independent_moc4_T2_uniform_20261022 (primary proof of paper case
T_6 = {0,1,2,4,6,8}).

Own code only; uses the shared library ../lib/moc4lib.py.  Nothing from the supplement is
read or executed.  Prints one line per check; exits nonzero if an EXACT check fails.
Symbolic identities are checked as rational-function identities in (v, w = conj v) and
the real parameters, with sympy.cancel.
Run:  python -I -B independent_moc4_T2_uniform_20261022__check.py
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


T = frozenset({0, 1, 2, 4, 6, 8})
check("T equals paper T_6", T == L.CONTACT_SETS[6])
gens = L.generators(T)
check("generators of T = {0,1,2,3,4,5,6,9,10,21,22,23,25}",
      gens == [0, 1, 2, 3, 4, 5, 6, 9, 10, 21, 22, 23, 25], str(gens))
check("nonminor generators are exactly 22, 25",
      [j for j in gens if not L.F_IS_MINOR[j]] == [22, 25])
check("f1, f9 minor; f22 nonminor", L.F_IS_MINOR[1] and L.F_IS_MINOR[9] and not L.F_IS_MINOR[22])
check("rank of generator span = 8", L.rank([list(L.F_VALUES[j]) for j in gens]) == 8)

# ---------------- linear reconstruction from q ----------------
NS = L.nullspace([L.P[s] for s in sorted(T)])
Wvals = [L.values(c) for c in NS]
COORD = [3, 7, 9, 10, 12, 13, 16, 18]
q = sp.symbols("q0:8")
Mq = sp.Matrix([[Wvals[a][k] for a in range(8)] for k in COORD])
lam = Mq.inv() * sp.Matrix(q)
GQ = [sp.expand(sum(lam[a] * Wvals[a][s] for a in range(8))) for s in range(24)]
q0, q1, q2, q3, q4, q5, q6, q7 = q
claimed = {5: q0, 11: -q1 + q2 + q3, 14: q4, 15: -q0 - q1 + q2 + q5, 17: -q0 - q1 + q2 + q6,
           19: -q1 - q4 + q5 + q7, 20: -q1 + q3 + q7,
           21: -q0 - 2 * q1 + q2 + q3 - q4 + q5 + q7, 22: -q1 - q4 + q6 + q7,
           23: -q0 - q1 + q2 - q4 + q6 + q7}
check("reconstruction formulas g5..g23 from q", all(sp.expand(GQ[k] - v) == 0 for k, v in claimed.items()))
check("contact coordinates vanish", all(GQ[k] == 0 for k in T))
# generators lie in W
check("all 13 generators lie in {g_T = 0}", all(all(L.F_VALUES[j][s] == 0 for s in T) for j in gens))

# ---------------- dual of the minor cone ----------------
minor = [j for j in gens if L.F_IS_MINOR[j]]
Gm = [[Fraction(L.F_VALUES[j][k]) for k in COORD] for j in minor]   # generators in q-coords
check("11 minor generators have rank 8", L.rank(Gm) == 8)
rays = set()
for sub in itertools.combinations(range(len(minor)), 7):
    rows = [Gm[i] for i in sub]
    if L.rank(rows) != 7:
        continue
    n = L.nullspace(rows)[0]
    vals = [sum(n[k] * Gm[i][k] for k in range(8)) for i in range(len(minor))]
    if all(x >= 0 for x in vals):
        pass
    elif all(x <= 0 for x in vals):
        n = [-x for x in n]
    else:
        continue
    # normalize
    m = next(x for x in n if x != 0)
    n = tuple(x / abs(m) for x in n)
    rays.add(n)
check("minor cone has exactly 18 extreme dual rays", len(rays) == 18, f"got {len(rays)}")
Gall = {j: [Fraction(L.F_VALUES[j][k]) for k in COORD] for j in gens}
bad = [r for r in rays if any(sum(r[k] * Gall[j][k] for k in range(8)) < 0 for j in gens)]
A_expr = -q0 - 2 * q1 + q2 + q3 - q4 + q6 + q7
B_expr = -q0 - q1 + q2 - q4 + q5 + q7


def as_lin(r):
    return sp.expand(sum(sp.Rational(r[k].numerator, r[k].denominator) * q[k] for k in range(8)))


badlin = [as_lin(r) for r in bad]


def prop(e1, e2):
    return sp.simplify(e1 / e2).is_number and sp.simplify(e1 / e2) > 0


check("exactly two dual rays are negative on a source generator, and they are A and B of (2)",
      len(bad) == 2 and any(prop(b, A_expr) for b in badlin) and any(prop(b, B_expr) for b in badlin),
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


def chart_X(vv, rr, ss, hh, tt):
    """X_sigma = (F_sigma - F_0)/((z-c)(v-r)) for the chart (3), as RF in the given inputs."""
    z = hh * vv / (1 + (hh - 1) * vv)
    c = (z * (vv - rr) - tt * rr * (vv - z)) / ((vv - rr) - tt * (vv - z))
    Cl = [rr + ss * (vv - rr), vv, rr, c]
    Al = [None, RF(0), RF(1), z]
    F = []
    for p in L.PERMS:
        f = RF(1)
        for i in (1, 2, 3):
            f = f * (Al[i] - Cl[p[i]])
        F.append(f)
    den = (z - c) * (vv - rr)
    return [(F[j] - F[0]) / den for j in range(24)], z, c


def chart_gaps(vv, ww, rr, ss, hh, tt):
    Xv, z, c = chart_X(vv, rr, ss, hh, tt)
    Xw, zw, cw = chart_X(ww, rr, ss, hh, tt)
    return [Xv[j] - Xw[j] for j in range(24)], z, c


G, z, c = chart_gaps(v, w, r, s, h, t)
check("chart (3) with kappa = i mu/((z-c)(v-r)): gaps at T vanish identically",
      all(G[j].iszero() for j in T))
check("inverse h = z(v-1)/(v(z-1))", eqz(z * (v - 1) / (v * (z - 1)), h))
check("inverse t = (c-z)(v-r)/((c-r)(v-z))", eqz((c - z) * (v - r) / ((c - r) * (v - z)), t))
check("1+(h-1)v = (v-1)/(z-1)", eqz(1 + (h - 1) * v, (v - 1) / (z - 1)))
check("(v-r)-t(v-z) = (v-r)(z-r)/(c-r)", eqz((v - r) - t * (v - z), (v - r) * (z - r) / (c - r)))

# four circles from the four real ratios, with generic symbolic spectra (sympy, small)
a0, a1, a2, a3, c0, c1, c2, c3 = sp.symbols("a0 a1 a2 a3 c0 c1 c2 c3")
AA = [a0, a1, a2, a3]
CC = [c0, c1, c2, c3]


def Fg(kk):
    p = L.PERMS[kk]
    return sp.Mul(*[AA[i] - CC[p[i]] for i in range(4)])


def Dg(i, j):
    return sp.expand(Fg(i) - Fg(j))


def cr(x1, x2, x3, x4):
    return (x1 - x3) * (x2 - x4) / ((x1 - x4) * (x2 - x3))


ratios = [((0, 6), (2, 8)), ((0, 2), (6, 8)), ((2, 4), (0, 1)), ((1, 4), (0, 2))]
circ = [(a0, a1, a2, c2), (a0, c0, c1, c2), (a1, a2, a3, c1), (a3, c1, c2, c3)]
okc = True
for (pa, pb), quad in zip(ratios, circ):
    R = sp.cancel(Dg(*pa) / Dg(*pb))
    print("INFO ratio D%d%d/D%d%d =" % (pa + pb), sp.factor(R))
    found = False
    for perm in itertools.permutations(quad):
        X = cr(*perm)
        for sgn_ in (1, -1):
            if sp.cancel(R - sgn_ * X) == 0:
                found = True
                break
        if found:
            break
    okc &= found
check("each of D06/D28, D02/D68, D24/D01, D14/D02 is +- a cross ratio of the stated four points", okc)

# ---------------- gap formulas ----------------
rho = -(v - w)
tau = 1 / t
N = v * w
O1 = (v - 1) * (w - 1)
aa = h * r - h - r
D = (aa * v + r) * (aa * w + r)
H = h * r * (r - 1)
K = H / D
J = (aa * N + r * (v + w) - r) / D
Lq = r / ((h - 1) * N)
Mq_ = h * (r - 1) / ((h - 1) * O1)
alpha = 1 + J
beta = 1 + h * J
E = rho * s * (tau - 1)
g = G
formulas = {3: rho * K * (tau - 1), 7: rho * Lq * tau * (s - 1), 10: rho * Mq_ * tau * (s - 1),
            12: -rho * s, 9: g[7] + s * g[3], 11: g[10] + s * g[3], 18: g[7] + E, 20: g[10] + E,
            15: g[12] + E * J, 17: g[12] + E * h * J, 19: g[3] + E * (alpha - K),
            22: g[3] + E * (beta - K), 21: g[10] + E * alpha, 23: g[7] + E * beta}
for kk, f in formulas.items():
    check(f"gap formula g{kk}", eqz(g[kk], f))


def lin_eval(expr, vals):
    """evaluate a sympy linear form in q0..q7 on RF values"""
    out = RF(0)
    for qi, val in zip(q, vals):
        cf = expr.coeff(qi)
        if cf != 0:
            out = out + RF(flint.fmpq(int(sp.numer(cf)), int(sp.denom(cf)))) * val
    return out


qs = [g[kk] for kk in COORD]
Aval = lin_eval(A_expr, qs)
Bval = lin_eval(B_expr, qs)
check("A = g10 + E beta", eqz(Aval, g[10] + E * beta))
check("B = g7 + E alpha", eqz(Bval, g[7] + E * alpha))
check("chart gaps satisfy the linear reconstruction (all 24)",
      all(eqz(G[kk], lin_eval(GQ[kk], qs)) for kk in range(24)))
check("(4) K + aJ = 1", eqz(K + aa * J, 1))
check("(4) M - L = J D/((h-1) N O1)", eqz(Mq_ - Lq, J * D / ((h - 1) * N * O1)))
check("A + B = g21 + g23", eqz(Aval + Bval, g[21] + g[23]))
check("AB - g21 g23 = -E g7 J^2 D/(L N O1)",
      eqz(Aval * Bval - g[21] * g[23], -E * g[7] * J * J * D / (Lq * N * O1)))
check("L/M = r O1/(h (r-1) N)", eqz(Lq / Mq_, r * O1 / (h * (r - 1) * N)))

# ---------------- Case III: p = -tau, q = 1/(s-1), scale q/(-rho) ----------------
p3 = -tau
q3 = 1 / (s - 1)
sc3 = q3 / (-rho)
check("Case III g19 = (p+1)((q+1)alpha - K)", eqz(g[19] * sc3, (p3 + 1) * ((q3 + 1) * alpha - K)))
check("Case III g22 = (p+1)((q+1)beta - K)", eqz(g[22] * sc3, (p3 + 1) * ((q3 + 1) * beta - K)))
check("Case III: scaled E = (p+1)(q+1)", eqz(E * sc3, (p3 + 1) * (q3 + 1)))

# ---------------- Case I: p = 1/(tau-1), q = -1/s, scale pq/rho ----------------
p1 = 1 / (tau - 1)
q1_ = -1 / s
sc1 = p1 * q1_ / rho
R1 = (1 + p1) * (1 + q1_)
bI = -Lq * R1
cI = -Mq_ * R1
caseI = {3: K * q1_, 12: p1, 7: bI, 10: cI, 9: bI - K, 11: cI - K, 18: bI - 1, 20: cI - 1,
         15: p1 - J, 17: p1 - h * J, 19: K * q1_ - (alpha - K), 22: K * q1_ - (beta - K),
         21: cI - alpha, 23: bI - beta}
gI = [g[kk] * sc1 for kk in range(24)]
for kk, f in caseI.items():
    check(f"Case I scaled g{kk}", eqz(gI[kk], f))
check("Case I A = c - beta, B = b - alpha", eqz(Aval * sc1, cI - beta) and eqz(Bval * sc1, bI - alpha))
check("outer: K + L alpha = r^2 (r-1)|1+(h-1)v|^2/((h-1) N D)",
      eqz(K + Lq * alpha, r * r * (r - 1) * (1 + (h - 1) * v) * (1 + (h - 1) * w) / ((h - 1) * N * D)))
check("outer: K + M beta = h r (r-1)^2 |1+(h-1)v|^2/((h-1) O1 D)",
      eqz(K + Mq_ * beta, h * r * (r - 1) ** 2 * (1 + (h - 1) * v) * (1 + (h - 1) * w) / ((h - 1) * O1 * D)))


# swap A1 <-> A2 : the plane map x -> 1 - x with labels A1, A2 exchanged
def quantities(vv, ww, rr, hh):
    aa_ = hh * rr - hh - rr
    D_ = (aa_ * vv + rr) * (aa_ * ww + rr)
    K_ = hh * rr * (rr - 1) / D_
    J_ = (aa_ * vv * ww + rr * (vv + ww) - rr) / D_
    L_ = rr / ((hh - 1) * vv * ww)
    M_ = hh * (rr - 1) / ((hh - 1) * (vv - 1) * (ww - 1))
    return J_, L_, M_, K_


Js, Ls, Ms, Ks = quantities(1 - v, 1 - w, 1 - r, 1 / h)
check("swap: (J, L, M, K) -> (hJ, M, L, K)", eqz(Js, h * J) and eqz(Ls, Mq_) and eqz(Ms, Lq) and eqz(Ks, K))
tr12 = (0, 2, 1, 3)
Gs, _, _ = chart_gaps(1 - v, 1 - w, 1 - r, s, 1 / h, t)
# with mu -> -mu (and rhohat -> -rhohat) the swapped configuration reproduces the gaps
okp = all(eqz(-Gs[L.perm_index(L.compose(L.PERMS[kk], tr12))], G[kk]) for kk in range(24))
check("swap (r,h,v,mu) -> (1-r,1/h,1-v,-mu), s,t fixed: g'_{sigma o (12)} = g_sigma, so rho is preserved", okp)
check("swap fixes T", L.g_act_set((L.PERMS[0], tr12, 0), T) == T)
T1 = 1 + (2 - h) * J
ll = -Lq
mm = -Mq_
check("(5) K - l T1 = (l/m)(K - m alpha)", eqz(K - ll * T1, (ll / mm) * (K - mm * alpha)))
check("l T1 - m alpha = (m - l)(K/m - alpha)", eqz(ll * T1 - mm * alpha, (mm - ll) * (K / mm - alpha)))


def comb(d):
    return [sum((cf * L.F_VALUES[j][kk] for j, cf in d.items()), RF(0)) for kk in range(24)]


rem6 = comb({1: eta - h * J, 3: eta, 5: cI - alpha, 9: eta + 1, 10: bI - eta - 1, 22: J - eta,
             23: eta + 1 - K})
sub6 = comb({4: p1 - J, 21: K * q1_ - (alpha - K)})
ok6 = all(eqz(gI[kk], sub6[kk] + rem6[kk]) for kk in range(24))
check("(6) g = (p-J)E02 + [Kq-(alpha-K)]E31 + remainder (6), all 24 coordinates, eta free", ok6)
check("E02 = f4 = n_{1,3}, E31 = f21 = n_{4,2}",
      L.F_VALUES[4] == tuple(1 if p_[0] == 2 else 0 for p_ in L.PERMS) and
      L.F_VALUES[21] == tuple(1 if p_[3] == 1 else 0 for p_ in L.PERMS))
check("(eta-hJ)(eta+1)-(J-eta)^2 = T1 eta - J(J+h)",
      eqz((eta - h * J) * (eta + 1) - (J - eta) * (J - eta), T1 * eta - J * (J + h)))
check("K - 1 = -aJ and 1 + a = (r-1)(h-1) (so K-1 <= J when 0<r<1, h<0, J>0)",
      eqz(K - 1, -aa * J) and eqz(1 + aa, (r - 1) * (h - 1)))
check("(1+J)^2 - T1 = J(J+h)", eqz(alpha * alpha - T1, J * (J + h)))

# ---------------- Case II: p = tau/(1-tau), q = s/(1-s), scale (1+p)(1+q)/(-rho) -------
p2 = tau / (1 - tau)
q2 = s / (1 - s)
sc2 = (1 + p2) * (1 + q2) / (-rho)
bII = Lq * p2
cII = Mq_ * p2
caseII = {3: K * (1 + q2), 12: q2 * (1 + p2), 7: bII, 10: cII, 15: q2 * (p2 + alpha),
          17: q2 * (p2 + beta), 19: K + q2 * alpha, 22: K + q2 * beta, 21: cII + q2 * alpha,
          23: bII + q2 * beta}
gII = [g[kk] * sc2 for kk in range(24)]
for kk, f in caseII.items():
    check(f"Case II scaled g{kk}", eqz(gII[kk], f))
check("Case II A = c + q beta, B = b + q alpha",
      eqz(Aval * sc2, cII + q2 * beta) and eqz(Bval * sc2, bII + q2 * alpha))
x_ = -J
Pp = h * x_ - 1
T2 = 2 * Pp + alpha
check("T2 = (2h-1)x - 1 = P + (h-1)x", eqz(T2, (2 * h - 1) * x_ - 1) and eqz(T2, Pp + (h - 1) * x_))
check("(7) K - M T2 = (M/L)(K - L P)", eqz(K - Mq_ * T2, (Mq_ / Lq) * (K - Lq * Pp)))
check("M T2 - L P = (L - M)(K/L - P)", eqz(Mq_ * T2 - Lq * Pp, (Lq - Mq_) * (K / Lq - Pp)))
Cst = Lq * p2 / q2
Bst = Mq_ * p2 / q2
sub8 = comb({4: q2 * (p2 - Pp), 21: K - q2 * Pp})
remII = [(gII[kk] - sub8[kk]) / q2 for kk in range(24)]
remII_sw = [remII[L.perm_index(L.compose(L.PERMS[kk], tr12))] for kk in range(24)]
rhs8 = comb({1: eta + alpha, 3: eta + 1, 5: Cst - Pp, 9: eta, 10: Bst - eta, 22: Pp - eta, 23: eta + K})
ok8 = all(eqz(remII_sw[kk], rhs8[kk]) for kk in range(24))
check("(8) after subtracting q(p-P)E02 + (K-qP)E31, dividing by q, swapping A1<->A2 coordinates; eta free", ok8)
if not ok8:
    print("INFO (8) without the coordinate swap holds:", all(eqz(remII[kk], rhs8[kk]) for kk in range(24)))
    diff = [kk for kk in range(24) if not eqz(remII_sw[kk], rhs8[kk])]
    print("INFO (8) differing coordinates:", diff)
etaII = Pp * Pp / T2
check("eta(eta+alpha) = (P-eta)^2 at eta = P^2/T2", eqz(etaII * (etaII + alpha), (Pp - etaII) ** 2))
check("(P + alpha)^2 = P^2 + alpha T2 (so eta >= -alpha)", eqz((Pp + alpha) ** 2, Pp * Pp + alpha * T2))

print(f"SUMMARY failures={len(FAIL)}")
sys.exit(1 if FAIL else 0)
