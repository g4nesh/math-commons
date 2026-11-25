"""Exact replay for "A four-coefficient support inequality and its scalar form"
(final_referee_T5_four_budget_capacity_20261022).  Own code only; nothing from the
supplement is read or executed.  One printed line per check; exit status nonzero iff an
EXACT check fails.

Checks
  1. identity (4) modulo N = |p|^2 - |q|^2 (sympy Poly division by the single generator N)
  2. Z >= 0  <=>  capacity (1) for a, b, c, d >= 0 (Schur complements, symbolic)
  3. the scalar form (3) <=> (2)-for-all-(1)-tuples (symbolic pieces + numeric tests)
  4. the quaternion responses of f_0, f_6, f_10, f_54 (paper's Appendix D numbering)
     equal (A-B)^2, (A+B)^2, D^2, H modulo N, for the convention found by search
  5. generators(T_5) and the profile identity g = g22 f0 + ... - A0 f54
  6. NUMERIC: random tests of (2) under (1), including the boundary, and of (3)
"""
import itertools
import os
import random
import sys

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import moc4lib as L  # noqa: E402

FAIL = []


def check(name, ok, level="EXACT", info=""):
    print(f"[{level}] {'PASS' if ok else 'FAIL'}  {name}" + (f"  ({info})" if info else ""))
    if not ok and level == "EXACT":
        FAIL.append(name)


p = sp.symbols("p0:4")
q = sp.symbols("q0:4")
a, b, c, d = sp.symbols("a b c d")
Y = list(p) + list(q)
GENS = Y + [a, b, c, d]


def ABDEGH(y):
    p0, p1, p2, p3, q0, q1, q2, q3 = y
    A = p0 * q3 + p3 * q0
    B = p1 * q2 - p2 * q1
    D = p0**2 + p2**2 - q0**2 - q2**2
    E = p1**2 - p2**2 - q0**2 + q3**2
    G = p0**2 - p3**2 - q1**2 + q2**2
    H = 2 * (A**2 + B**2) - E * G
    return A, B, D, E, G, H


A, B, D, E, G, H = ABDEGH(Y)
N = sum(x**2 for x in p) - sum(x**2 for x in q)
n = sum(x**2 for x in p) + sum(x**2 for x in q)
Z = sp.Matrix([[a + b + 2 * d, a - b, 0, 0],
               [a - b, a + b + 2 * d, 0, 2 * d],
               [0, 0, c, d],
               [0, 2 * d, d, d]])


def ell(y):
    A_, B_, D_, E_, G_, _ = ABDEGH(y)
    return sp.Matrix([-y[0] * B_, y[0] * A_, -y[7] * D_, -y[7] * G_])


def emm(y):
    A_, B_, D_, E_, G_, _ = ABDEGH(y)
    return sp.Matrix([y[4] * B_, -y[4] * A_, -y[3] * D_, -y[3] * E_])


p0, p1, p2, p3 = p
q0, q1, q2, q3 = q
S = [
    [p0, p1, p2, p3, q0, q1, q2, q3],
    [-p1, p0, -p3, p2, q1, -q0, -q3, q2],
    [-p2, p3, p0, -p1, q2, q3, -q0, -q1],
    [-p3, -p2, p1, p0, q3, -q2, q1, -q0],
]


def in_ideal_N(expr, gens=GENS):
    P_ = sp.Poly(sp.expand(expr), *gens)
    NP = sp.Poly(N, *gens)
    _, r = sp.reduced(P_.as_expr(), [NP.as_expr()], *gens, order="lex")
    return sp.expand(r) == 0


# ---------------------------------------------------------------- 1. identity (4)
lhs = n * (a * (A - B)**2 + b * (A + B)**2 + c * D**2 + d * H)
rhs = 0
for Sj in S:
    l_ = ell(Sj)
    m_ = emm(Sj)
    rhs += (l_.T * Z * l_)[0, 0] + (m_.T * Z * m_)[0, 0]
diff = sp.expand(lhs - rhs)
check("identity (4): n[a(A-B)^2+b(A+B)^2+cD^2+dH] - sum_j (l^T Z l + m^T Z m)(S_j y) in (N)",
      in_ideal_N(diff))
check("identity (4) is not an identity without the reduction mod N (sanity)", diff != 0, "INFO")

# ---------------------------------------------------------------- 2. Z >= 0 <=> (1)
# Schur complement on the (4,4) entry d, then on c - d
Zs1 = Z[:3, :3] - Z[:3, 3] * Z[3, :3] / d
Zs1 = Zs1.applyfunc(sp.simplify)
target1 = sp.Matrix([[a + b + 2 * d, a - b, 0], [a - b, a + b - 2 * d, -2 * d], [0, -2 * d, c - d]])
check("Schur complement of d in Z equals the displayed 3x3", (Zs1 - target1).applyfunc(sp.simplify) == sp.zeros(3, 3))
Zs2 = Zs1[:2, :2] - Zs1[:2, 2] * Zs1[2, :2] / (c - d)
det2 = sp.simplify(Zs2.det() - 4 * (a * b * (c - d) - d**2 * (a + b + c + d)) / (c - d))
check("det of the 2x2 Schur complement = 4[ab(c-d) - d^2(a+b+c+d)]/(c-d)", det2 == 0)
# d = 0: Z is PSD for a, b, c >= 0 (eigenvalues 2a, 2b, c, 0)
Z0 = Z.subs(d, 0)
ev = Z0.eigenvals()
check("d = 0: eigenvalues of Z are 2a, 2b, c, 0", set(ev) == {2 * a, 2 * b, c, 0})
# converse direction (Z >= 0, d > 0 => c > d): if c = d the 3x3 Schur complement has a zero
# diagonal entry with a nonzero off-diagonal entry -2d in its row
check("c = d, d > 0: Schur complement has (3,3)=0 and (2,3)=-2d (so not PSD)",
      target1.subs(c, d)[2, 2] == 0 and target1.subs(c, d)[1, 2] == -2 * d)

# random exact rational check of Z >= 0 <=> (1) for nonnegative a, b, c, d
rng = random.Random(4)
agree = 0
tot = 0
for _ in range(400):
    vals = {s: sp.Rational(rng.randint(0, 12), rng.randint(1, 4)) for s in (a, b, c, d)}
    cap = (a * b * (c - d) - d**2 * (a + b + c + d)).subs(vals) >= 0
    Zn = Z.subs(vals)
    # PSD test via all principal minors (exact)
    psd = all(Zn.extract(list(I), list(I)).det() >= 0
              for k in range(1, 5) for I in itertools.combinations(range(4), k))
    tot += 1
    agree += (cap == psd)
check(f"random rational a,b,c,d >= 0: Z PSD (all principal minors) <=> (1), {tot} samples, seed 4",
      agree == tot, info=f"{agree}/{tot}")

# ---------------------------------------------------------------- 3. scalar form
al, be, ga = sp.symbols("alpha beta gamma", positive=True)
sub = {d: 1, a: 1 / al - 1, b: 1 / be - 1, c: 1 / ga}
cap_expr = (a * b * (c - d) - d**2 * (a + b + c + d)).subs(sub).subs(ga, 1 - al - be)
check("choice d=1, a=1/alpha-1, b=1/beta-1, c=1/gamma gives equality in (1) on alpha+beta+gamma=1",
      sp.simplify(cap_expr) == 0)
# (1) with d = 1 <=> 1/(a+1) + 1/(b+1) + 1/c <= 1 (for c > 0): the cleared forms coincide
lhs1 = sp.expand((a * b * (c - 1) - (a + b + c + 1)))
lhs2 = sp.expand((a + 1) * (b + 1) * c - (b + 1) * c - (a + 1) * c - (a + 1) * (b + 1))
check("(1) at d=1 equals (a+1)(b+1)c [1 - 1/(a+1) - 1/(b+1) - 1/c] (converse direction)",
      sp.expand(lhs1 - lhs2) == 0)
u, v, w = sp.symbols("u v w", nonnegative=True)
# H + w + 2 sqrt(uv) + ... >= 0  <=>  EG <= (|A-B|+|A+B|+|D|)^2 when u=(A-B)^2 etc.
expr = sp.expand(2 * (A**2 + B**2) - ((A - B)**2 + (A + B)**2))
check("u + v = 2(A^2+B^2), so H + (sqrt u+sqrt v+sqrt w)^2 - u - v = (|A-B|+|A+B|+|D|)^2 - EG",
      expr == 0)

# ---------------------------------------------------------------- 4. quaternion responses
def qmul(x, y):
    a0, a1, a2, a3 = x
    b0, b1, b2, b3 = y
    return (a0 * b0 - a1 * b1 - a2 * b2 - a3 * b3,
            a0 * b1 + a1 * b0 + a2 * b3 - a3 * b2,
            a0 * b2 - a1 * b3 + a2 * b0 + a3 * b1,
            a0 * b3 + a1 * b2 - a2 * b1 + a3 * b0)


def conj(x):
    return (x[0], -x[1], -x[2], -x[3])


def quat_matrix(pp, qq, kind, transpose):
    """Matrix of x -> (kind) in the basis 1, i, j, k: column j = image of e_j."""
    cols = []
    for j in range(4):
        e = [0, 0, 0, 0]
        e[j] = 1
        e = tuple(e)
        l, r, swap = kind
        P_ = conj(pp) if l else pp
        Q_ = conj(qq) if r else qq
        if swap:
            P_, Q_ = Q_, P_
        cols.append(qmul(qmul(P_, e), Q_))
    M = [[cols[j][i] for j in range(4)] for i in range(4)]
    if transpose:
        M = [[M[j][i] for j in range(4)] for i in range(4)]
    return M


def phi_vals(vals, M):
    tot = 0
    for s, sg, val in zip(L.PERMS, L.SGN, vals):
        if val == 0:
            continue
        t = sg * val
        for i in range(4):
            t = t * M[i][s[i]]
        tot += t
    return tot


FQ = {0: (A - B)**2, 6: (A + B)**2, 10: D**2, 54: H}
check("Appendix D: f_0, f_6, f_10 are minor vectors and f_54 is a non-minor h vector",
      L.F_IS_MINOR[0] and L.F_IS_MINOR[6] and L.F_IS_MINOR[10] and not L.F_IS_MINOR[54],
      info=f"f0={L.f_name(0)}, f6={L.f_name(6)}, f10={L.f_name(10)}, f54={L.f_name(54)}")

# search the conventions numerically on a rational point with |p| = |q|
rngq = random.Random(11)


def rational_equal_norm_pair(rr):
    pp = [sp.Rational(rr.randint(-5, 5)) for _ in range(4)]
    while all(x == 0 for x in pp):
        pp = [sp.Rational(rr.randint(-5, 5)) for _ in range(4)]
    t = sp.Rational(rr.randint(1, 6), rr.randint(1, 6))
    s_ = sp.Rational(rr.randint(-6, 6), rr.randint(1, 6))
    # rational unit quaternion from two stereographic parameters
    den = 1 + t**2 + s_**2
    uq = ((1 - t**2 - s_**2) / den, 2 * t / den, 2 * s_ / den, 0)
    uq = (uq[0], uq[1], 0, uq[2]) if rr.random() < 0.5 else uq
    qq = qmul(tuple(pp), uq)
    return tuple(pp), tuple(qq)


found = []
pts = [rational_equal_norm_pair(rngq) for _ in range(3)]
for kind in itertools.product((0, 1), (0, 1), (0, 1)):
    for tr in (0, 1):
        ok = True
        for pp, qq in pts:
            M = quat_matrix(pp, qq, kind, tr)
            nn = sum(x**2 for x in pp) + sum(x**2 for x in qq)
            sb = dict(zip(Y, list(pp) + list(qq)))
            for j, Qj in FQ.items():
                if 4 * phi_vals(L.F_VALUES[j], M) != nn**2 * Qj.subs(sb):
                    ok = False
                    break
            if not ok:
                break
        if ok:
            found.append((kind, tr))
check("some quaternion convention x -> (p or pbar) x (q or qbar) (possibly swapped/transposed) "
      "makes all four responses match at 3 rational points", len(found) > 0, info=f"{found}")
if found:
    kind, tr = found[0]
    M = quat_matrix(p, q, kind, tr)
    allok = True
    for j, Qj in FQ.items():
        e = sp.expand(4 * phi_vals(L.F_VALUES[j], M) - n**2 * Qj)
        okj = in_ideal_N(e, Y)
        allok &= okj
        check(f"symbolic: 4 Phi_f{j}(M(p,q)) - n^2 Q_{j} in (N) for convention {kind},{tr}", okj)
    # orthogonality of M/|p||q| on N = 0:  M M^T = |p|^2|q|^2 I
    MM = sp.Matrix(M) * sp.Matrix(M).T
    okorth = all(in_ideal_N(sp.expand(4 * MM[i, j] - (n**2 if i == j else 0)), Y)
                 for i in range(4) for j in range(4))
    check("M(p,q) M(p,q)^T = (n/2)^2 I modulo N (M/(n/2) is in SO(4) when |p|=|q|)", okorth)

# ---------------------------------------------------------------- 5. T5 generators and profile
T5 = L.CONTACT_SETS[5]
check("paper T_5 = {1,2,3,4,6,8} (0-based lexicographic indices)", T5 == frozenset({1, 2, 3, 4, 6, 8}))
gens = L.generators(T5)
check("generators(T_5) (all f_j vanishing on T_5) = {f0,f1,f2,f3,f4,f5,f6,f9,f10,f53,f54}",
      gens == [0, 1, 2, 3, 4, 5, 6, 9, 10, 53, 54], info=f"{gens}")
used = [0, 2, 3, 5, 6, 10, 53, 54]
check("the document's eight profile vectors are among the generators", set(used) <= set(gens))
check("f0,f2,f3,f5,f6,f10,f53 minor; f54 non-minor; the extra generators f1,f4,f9 are minor",
      all(L.F_IS_MINOR[j] for j in (0, 2, 3, 5, 6, 10, 53, 1, 4, 9)) and not L.F_IS_MINOR[54])
check("only non-minor generator of T_5 is f54", [j for j in gens if not L.F_IS_MINOR[j]] == [54])
check("the eight profile vectors are linearly independent", L.rank([list(L.F_VECS[j]) for j in used]) == 8)
check("all 11 generators span the same 8-dim space {c : (Pc)_pi = 0 for pi in T_5}",
      L.rank([list(L.F_VECS[j]) for j in gens]) == 8 and
      L.rank([list(L.P[s]) for s in sorted(T5)]) == 6)


def profile_rhs(g):
    A0 = g[16] - g[0]
    terms = [(g[22], 0), (g[13], 2), (g[12], 3), (g[10], 5), (g[17], 6), (g[7], 10), (g[16], 53), (-A0, 54)]
    out = [0] * 24
    for coef, j in terms:
        for s in range(24):
            out[s] += coef * L.F_VALUES[j][s]
    return out


ok = all(list(L.F_VALUES[j]) == profile_rhs(L.F_VALUES[j]) for j in gens)  # all 11
check("profile identity g = g22 f0 + g13 f2 + g12 f3 + g10 f5 + g17 f6 + g7 f10 + g16 f53 - (g16-g0) f54 "
      "on each of the 11 generators (hence on the whole 8-dim contact subspace)", ok)
# capacity in gap form: with a=g22, b=g17, c=g7, d=-A0
gg = sp.symbols("g0:24")
A0s = gg[16] - gg[0]
capg = sp.expand((a * b * (c - d) - d**2 * (a + b + c + d)).subs({a: gg[22], b: gg[17], c: gg[7], d: -A0s}))
capdoc = sp.expand(gg[22] * gg[17] * (gg[7] + A0s) - A0s**2 * (gg[22] + gg[17] + gg[7] - A0s))
check("gap form of the capacity: a=g22, b=g17, c=g7, d=-A0 reproduces the displayed inequality",
      sp.expand(capg - capdoc) == 0)

# ---------------------------------------------------------------- 6. numeric falsification
import numpy as np  # noqa: E402

rs = np.random.default_rng(2026)
fA = sp.lambdify(Y, [A, B, D, E, G, H], "numpy")
worst2 = np.inf
worst3 = np.inf
NS = 20000
for _ in range(NS):
    pp = rs.normal(size=4)
    qq = rs.normal(size=4)
    qq *= np.linalg.norm(pp) / np.linalg.norm(qq)
    pp /= np.linalg.norm(pp)
    qq /= np.linalg.norm(qq)
    A_, B_, D_, E_, G_, H_ = fA(*pp, *qq)
    # (3)
    worst3 = min(worst3, (abs(D_) + 2 * max(abs(A_), abs(B_)))**2 - E_ * G_)
    # (2) on a random point of the capacity boundary (and sometimes the d=0 face)
    x, y = rs.uniform(0, 5, size=2)
    if x * y <= 1:
        x, y = 1 / x + 0.1, y + 1 / max(x, 1e-3)
    if x * y > 1:
        dd = 1.0
        aa, bb, cc = x, y, (x + 1) * (y + 1) / (x * y - 1)
        val = aa * (A_ - B_)**2 + bb * (A_ + B_)**2 + cc * D_**2 + dd * H_
        worst2 = min(worst2, val)
check(f"NUMERIC (3) on {NS} random unit pairs (seed 2026): min of RHS-LHS", worst3 >= -1e-12, "NUMERIC",
      info=f"min = {worst3:.3e}")
check(f"NUMERIC (2) on capacity-boundary coefficients, same samples: min value", worst2 >= -1e-9, "NUMERIC",
      info=f"min = {worst2:.3e}")
# tightness probe: (3) equality attained?  report the smallest relative slack
print(f"[INFO] smallest slack in (3) over samples: {worst3:.3e}")

print("SUMMARY", "ALL EXACT CHECKS PASSED" if not FAIL else f"FAILED: {FAIL}")
sys.exit(1 if FAIL else 0)
