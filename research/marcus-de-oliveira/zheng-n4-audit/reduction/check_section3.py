#!/usr/bin/env python3
"""Exact replays for the Section 3 / Section 4.2-4.3 / Section 6 reduction audit
(Zheng, "The Marcus-de Oliveira conjecture for n = 4", 2026).

Everything here is derived from the paper's printed definitions only:
  * the 14 noncrossing matchings B_0..B_13 of eq. (A.1),
  * the cubic C of eq. (A.2),
  * the entry-subspace coordinates of Appendix A.2.
No supplement file is read or executed; the shared audit library is not used.

Run:  python -I -B check_section3.py      (exit code 0 iff every check passes)
Writes check_section3_P.json (the 24 x 14 matrix P) next to this script.
"""
import itertools
import json
import os
import random
import sys
from functools import lru_cache

import sympy as sp

FAIL = []


def check(name, cond, detail=""):
    tag = "PASS" if cond else "FAIL"
    print(f"[{tag}] {name}" + (f"  ({detail})" if detail else ""))
    if not cond:
        FAIL.append(name)


# ---------------------------------------------------------------- permutations
PERMS = list(itertools.permutations([1, 2, 3, 4]))  # lexicographic, pi_0 = id
PIDX = {p: k for k, p in enumerate(PERMS)}


def sgn(p):
    s = 1
    for i in range(len(p)):
        for j in range(i + 1, len(p)):
            if p[i] > p[j]:
                s = -s
    return s


def compose(p, q):  # (p o q)(i) = p(q(i)), one-line tuples on 1..4
    return tuple(p[q[i] - 1] for i in range(4))


def inverse(p):
    r = [0] * 4
    for i in range(4):
        r[p[i] - 1] = i + 1
    return tuple(r)


# ---------------------------------------------------------------- straightening
BASIS = [
    ((1, 2), (3, 4), (5, 6), (7, 8)), ((1, 2), (3, 4), (5, 8), (6, 7)),
    ((1, 2), (3, 6), (4, 5), (7, 8)), ((1, 2), (3, 8), (4, 5), (6, 7)),
    ((1, 2), (3, 8), (4, 7), (5, 6)), ((1, 4), (2, 3), (5, 6), (7, 8)),
    ((1, 4), (2, 3), (5, 8), (6, 7)), ((1, 6), (2, 3), (4, 5), (7, 8)),
    ((1, 6), (2, 5), (3, 4), (7, 8)), ((1, 8), (2, 3), (4, 5), (6, 7)),
    ((1, 8), (2, 3), (4, 7), (5, 6)), ((1, 8), (2, 5), (3, 4), (6, 7)),
    ((1, 8), (2, 7), (3, 4), (5, 6)), ((1, 8), (2, 7), (3, 6), (4, 5)),
]  # eq. (A.1), transcribed from the page image (p. 35)
BIDX = {b: k for k, b in enumerate(BASIS)}


def crosses(p, q):
    (a, c), (b, d) = sorted([p, q])
    return a < b < c < d


@lru_cache(maxsize=None)
def straighten(key):
    """key: tuple of sorted pairs (i<j). Returns tuple of (noncrossing key, coeff)."""
    pairs = list(key)
    for s in range(4):
        for t in range(s + 1, 4):
            if crosses(pairs[s], pairs[t]):
                (a, c), (b, d) = sorted([pairs[s], pairs[t]])
                rest = [pairs[u] for u in range(4) if u not in (s, t)]
                out = {}
                # Pluecker: [ac][bd] = [ab][cd] + [ad][bc]   (a<b<c<d)
                for new in ([(a, b), (c, d)], [(a, d), (b, c)]):
                    for k2, v in straighten(tuple(sorted(rest + new))):
                        out[k2] = out.get(k2, 0) + v
                return tuple((k2, v) for k2, v in out.items() if v)
    return ((key, 1),)


def expand_matching(pairs):
    """Signed product of brackets [i j] (arbitrary order) -> coefficient vector in B basis."""
    s = 1
    norm = []
    for i, j in pairs:
        if i == j:
            return [0] * 14
        if i > j:
            i, j = j, i
            s = -s
        norm.append((i, j))
    vec = [0] * 14
    for k2, v in straighten(tuple(sorted(norm))):
        vec[BIDX[k2]] += s * v
    return vec


# noncrossing matchings of 1..8 (should be exactly the 14 of (A.1))
allm = set()
def gen(rem, cur):
    if not rem:
        allm.add(tuple(sorted(cur)))
        return
    a = rem[0]
    for b in rem[1:]:
        gen([r for r in rem if r not in (a, b)], cur + [(a, b)])
gen(list(range(1, 9)), [])
nc = {m for m in allm if not any(crosses(p, q) for p, q in itertools.combinations(m, 2))}
check("A.1 lists exactly the 14 noncrossing perfect matchings of 1..8",
      nc == set(BASIS) and len(allm) == 105, f"{len(nc)} noncrossing of {len(allm)}")

P = sp.Matrix([expand_matching([(i + 1, 4 + s[i]) for i in range(4)]) for s in PERMS])

# ---------------------------------------------------------------- polynomial cross-check of P
X = sp.symbols("x1:9")
Y = sp.symbols("y1:9")
GENS = X + Y


def br(i, j):
    return sp.Poly(X[i - 1] * Y[j - 1] - X[j - 1] * Y[i - 1], *GENS, domain="ZZ")


def mpoly(pairs):
    r = sp.Poly(1, *GENS, domain="ZZ")
    for i, j in pairs:
        r = r * br(i, j)
    return r


Bpoly = [mpoly(b) for b in BASIS]
ok = True
for k, s in enumerate(PERMS):
    lhs = mpoly([(i + 1, 4 + s[i]) for i in range(4)])
    rhs = sp.Poly(0, *GENS, domain="ZZ")
    for a in range(14):
        if P[k, a]:
            rhs = rhs + int(P[k, a]) * Bpoly[a]
    ok &= (lhs - rhs).is_zero
check("P (straightening) satisfies p_sigma = sum_a P[sigma,a] B_a as polynomials in 16 homogeneous variables", ok)
check("P is an integer matrix", all(v == int(v) for v in P))
check("rank P = 14", P.rank() == 14)
rows = [0, 1, 2, 3, 4, 6, 7, 8, 9, 10, 12, 13, 16, 18]
dsub = P.extract(rows, list(range(14))).det()
check("det of the 14x14 submatrix with rows 0,1,2,3,4,6,7,8,9,10,12,13,16,18 is 1 (App. A)", dsub == 1, f"det = {dsub}")

kappa = sp.zeros(14, 1)
kappa[13] = 1
check("P kappa = (1,...,1) with kappa = e_13 (3.2)", P * kappa == sp.ones(24, 1))


def Psi(v):
    """matching vector (B_0..B_13 values) of a configuration v = list of 8 (x,y) pairs"""
    def b(i, j):
        return v[i - 1][0] * v[j - 1][1] - v[j - 1][0] * v[i - 1][1]
    out = []
    for m in BASIS:
        r = 1
        for i, j in m:
            r *= b(i, j)
        out.append(sp.expand(r))
    return sp.Matrix(out)


kconf = [(0, 1)] * 4 + [(1, 1)] * 4
check("kappa = Psi(first four columns (0,1), last four (1,1)), so kappa in X", Psi(kconf) == kappa)

# ---------------------------------------------------------------- the cubic (A.2)
m = sp.symbols("m0:14")
m0, m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13 = m
C = (m13 * (m0 * m6 - m1 * m5) - m0 * m9 * (m0 + m1 + m2 + m3 + m5 + m7 + m9) - m1 * m5 * m9
     + m10 * (-m0 * m3 - m0 * m9 + m1 * m2 + m1 * m7)
     + m11 * (-m0 * m7 - m0 * m9 + m2 * m5 + m3 * m5)
     + m12 * (-m0 * m9 + m1 * m7 - m2 * m6 - m2 * m9 + m3 * m5 + m3 * m7)
     + m2 * m10 * m11 - m3 * m8 * m10 - m4 * m7 * m11 + m4 * m8 * m9)
grad = [sp.diff(C, mi) for mi in m]


def eval_on_B(expr):
    pe = sp.Poly(expr, *m)
    r = sp.Poly(0, *GENS, domain="ZZ")
    for mon, cf in pe.terms():
        t = sp.Poly(int(cf), *GENS, domain="ZZ")
        for a, e in enumerate(mon):
            for _ in range(e):
                t = t * Bpoly[a]
        r = r + t
    return r


ok = all(eval_on_B(g).is_zero for g in grad)
check("grad C vanishes identically on symbolic matching vectors (16 homogeneous vars) => X subset Sing(C)", ok)
check("C itself vanishes on matching vectors (Euler)", eval_on_B(C).is_zero)

# independence of the 14 partials via the listed quadratic monomials (informational; only X c Sing(C) is used)
mons = [m10 * m11, m9 * m12, m8 * m10, m8 * m9, m7 * m12, m4 * m11, m4 * m9, m4 * m8,
        m4 * m7, m3 * m12, m3 * m8, m3 * m7, m2 * m12, m1 * m5]
CM = sp.Matrix([[sp.Poly(g, *m).coeff_monomial(mo) for mo in mons] for g in grad])
check("App. A.1: 14x14 coefficient matrix of partials on the listed monomials has det +-1",
      abs(CM.det()) == 1, f"det = {CM.det()}")

# Hessian at kappa, (3.10)-(3.12)
kap_sub = {mi: (1 if i == 13 else 0) for i, mi in enumerate(m)}
H = sp.hessian(C, m).subs(kap_sub)
Hexp = sp.zeros(14, 14)
Hexp[0, 6] = Hexp[6, 0] = 1
Hexp[1, 5] = Hexp[5, 1] = -1
check("Hess(C)(kappa) = e0 e6^t + e6 e0^t - e1 e5^t - e5 e1^t (3.10)", H == Hexp)
check("rank Hess(C)(kappa) = 4", H.rank() == 4)
ns = H.nullspace()
kerdim = len(ns)
kerok = all(v[0] == 0 and v[1] == 0 and v[5] == 0 and v[6] == 0 for v in ns)
check("ker Hess(C)(kappa) = {m0=m1=m5=m6=0} (dim 10)", kerdim == 10 and kerok)
check("grad C(kappa) = 0", all(g.subs(kap_sub) == 0 for g in grad))
t = sp.Symbol("t")
shift = {mi: mi + (t if i == 13 else 0) for i, mi in enumerate(m)}
pol_ok = all(sp.expand(grad[i].subs(shift, simultaneous=True) - grad[i] - t * (H * sp.Matrix(m))[i]) == 0
             for i in range(14))
check("polarization grad C(x + t kappa) = grad C(x) + t H_kappa x, exactly (3.12)", pol_ok)

# ---------------------------------------------------------------- entry subspace E(Q), (3.7), (3.8)
q = sp.Matrix(4, 4, lambda i, j: sp.Symbol(f"q{i+1}{j+1}"))
rhs = sp.Matrix([sum(q[i, s[i] - 1] for i in range(4)) for s in PERMS])
Psub = P.extract(rows, list(range(14)))
EQ = Psub.LUsolve(rhs.extract(rows, [0]))
EQ = EQ.applyfunc(sp.expand)
check("P E(Q) = (sum_i q_{i,sigma(i)})_sigma is consistent on all 24 rows (3.7)",
      (P * EQ - rhs).applyfunc(sp.expand) == sp.zeros(24, 1))
Qc = lambda i, j: q[i - 1, j - 1]
paperE = [0, 0, Qc(1, 3) - Qc(1, 4) - Qc(2, 3) + Qc(2, 4), Qc(1, 2) - Qc(1, 3) - Qc(2, 2) + Qc(2, 3),
          Qc(1, 1) - Qc(1, 2) - Qc(2, 1) + Qc(2, 2), 0, 0, Qc(2, 3) - Qc(2, 4) - Qc(3, 3) + Qc(3, 4),
          Qc(3, 3) - Qc(3, 4) - Qc(4, 3) + Qc(4, 4), Qc(2, 2) - Qc(2, 3) - Qc(3, 2) + Qc(3, 3),
          Qc(2, 1) - Qc(2, 2) - Qc(3, 1) + Qc(3, 2), Qc(3, 2) - Qc(3, 3) - Qc(4, 2) + Qc(4, 3),
          Qc(3, 1) - Qc(3, 2) - Qc(4, 1) + Qc(4, 2), Qc(1, 4) + Qc(2, 3) + Qc(3, 2) + Qc(4, 1)]
check("E(Q) coordinates equal those printed in App. A.2", all(sp.expand(EQ[i] - paperE[i]) == 0 for i in range(14)))
qs = list(q)
Jac = sp.Matrix([[sp.diff(EQ[i], s) for s in qs] for i in range(14)])
check("dim E = 10 (rank of Q -> E(Q))", Jac.rank() == 10)


def cayley(rng):
    K = sp.zeros(4, 4)
    for i in range(4):
        for j in range(i + 1, 4):
            v = sp.Rational(rng.randint(-9, 9), rng.randint(1, 7))
            K[i, j], K[j, i] = v, -v
    I4 = sp.eye(4)
    O = (I4 - K) * (I4 + K).inv()
    assert O * O.T == I4 and O.det() == 1
    return O


def Phi(c, O, extended=False):
    """Leibniz sum (3.4); extended=True divides by det O (paper's extension to O(4), Sec. 4.2)"""
    Pc = P * c
    s = 0
    for k, p in enumerate(PERMS):
        prod = 1
        for i in range(4):
            prod *= O[i, p[i] - 1]
        s += sgn(p) * Pc[k] * prod
    return sp.expand(s / O.det()) if extended else sp.expand(s)


rng = random.Random(20261125)
ok38 = True
okneg = True
for _ in range(3):
    O = cayley(rng)
    lhs = Phi(EQ, O)
    rhs38 = sum(q[i, j] * O[i, j] ** 2 for i in range(4) for j in range(4))
    ok38 &= sp.expand(lhs - rhs38) == 0
    Om = O * sp.diag(-1, 1, 1, 1)
    okneg &= sp.expand(Phi(EQ, Om) + sum(q[i, j] * Om[i, j] ** 2 for i in range(4) for j in range(4))) == 0
check("(3.8) Phi_{E(Q)}(O) = sum q_ij O_ij^2 for symbolic Q, 3 rational O in SO(4)", ok38)
check("(3.8) fails with sign -1 when det O = -1 (so SO(4) is essential)", okneg)

# ---------------------------------------------------------------- minor vectors and (4.2)
Pinv_rows = Psub.inv()


def vec_from_values(vals):
    c = Pinv_rows * sp.Matrix([vals[r] for r in rows])
    return c if P * c == sp.Matrix(vals) else None


minors = {}
for size in (1, 2):
    for I in itertools.combinations(range(1, 5), size):
        for J in itertools.combinations(range(1, 5), size):
            vals = [1 if set(p[i - 1] for i in I) == set(J) else 0 for p in PERMS]
            minors[(I, J)] = vec_from_values(vals)
check("all 16+36 indicator vectors 1{sigma(I)=J} lie in P(V_R) (n_{I,J} exist)",
      all(v is not None for v in minors.values()))
distinct = {tuple(v) for v in minors.values()}
check("34 distinct minor vectors after complementary identification", len(distinct) == 34, f"{len(distinct)}")
ok42 = True
for _ in range(2):
    O = cayley(rng)
    for (I, J), c in minors.items():
        sub = O.extract([i - 1 for i in I], [j - 1 for j in J])
        ok42 &= sp.simplify(Phi(c, O) - sub.det() ** 2) == 0
check("(4.2) Phi_{n_{I,J}}(O) = det(O[I,J])^2 for all 52 (I,J), 2 rational O in SO(4)", ok42)

# ---------------------------------------------------------------- identity (3.5) on exact Gaussian-rational data
def grat(r):
    return sp.Rational(r.randint(-12, 12), r.randint(1, 5)) + sp.I * sp.Rational(r.randint(-12, 12), r.randint(1, 5))


ok35 = True
okF = True
for trial in range(3):
    O = cayley(rng)
    a = [grat(rng) for _ in range(4)]
    b = [grat(rng) for _ in range(4)]
    alpha = grat(rng)
    v = [(a[i], 1) for i in range(4)] + [(-b[j], 1) for j in range(4)]
    x = (alpha * Psi(v)).applyfunc(sp.expand)
    Fx = (P * x).applyfunc(sp.expand)
    for k, s in enumerate(PERMS):
        okF &= sp.expand(Fx[k] - alpha * sp.prod([a[i] + b[s[i] - 1] for i in range(4)])) == 0
    Da, Db = sp.diag(*a), sp.diag(*b)
    D = sp.expand((Da + O * Db * O.T).det())
    for k, pstar in enumerate(PERMS):
        c = x.applyfunc(sp.re) - sp.re(Fx[k]) * kappa
        lhs = Phi(c, O)
        rhs35 = sp.re(sp.expand(alpha * D - Fx[k]))
        ok35 &= sp.simplify(lhs - rhs35) == 0
        ok35 &= sp.simplify(P * c - Fx.applyfunc(lambda z: sp.re(z - Fx[k]))) == sp.zeros(24, 1)
check("(3.1) F_sigma(alpha Psi(v(a,b))) = alpha prod (a_i + b_sigma(i)) (3 exact trials)", okF)
check("(3.5) Phi_{c(x)}(O) = Re(alpha det(D_a + O D_b O^t) - F_pi*(x)) and (P c)_sigma = Re g_sigma, 3 trials x 24 references", ok35)

# fourth-root rescaling used in the proof of Prop. 4.3
tt = sp.Symbol("tt")
a = [grat(rng) for _ in range(4)]
b = [grat(rng) for _ in range(4)]
v1 = [(a[i], 1) for i in range(4)] + [(-b[j], 1) for j in range(4)]
v2 = [(tt * a[i], 1) for i in range(4)] + [(-tt * b[j], 1) for j in range(4)]
check("Psi(v(t a, t b)) = t^4 Psi(v(a,b)) (justifies t^4 = alpha in Prop. 4.3)",
      (Psi(v2) - tt ** 4 * Psi(v1)).applyfunc(sp.expand) == sp.zeros(14, 1))

# ---------------------------------------------------------------- dim X = 6 by exact Jacobian rank
xs, ys = sp.symbols("X1:9"), sp.symbols("Y1:9")
vsym = [(xs[i], ys[i]) for i in range(8)]
Bsym = Psi(vsym)
pt = [(3, 1), (-1, 2), (5, 7), (2, -3), (4, 1), (-2, 5), (1, 1), (7, -2)]
sub = {}
for i in range(8):
    sub[xs[i]], sub[ys[i]] = pt[i]
J = sp.Matrix([[sp.diff(Bsym[a], w) for w in list(xs) + list(ys)] for a in range(14)]).subs(sub)
rJ = J.rank()
check("rank dPsi = 6 at a rational configuration => dim X >= 6", rJ == 6, f"rank {rJ}")
fields = []
for Xg in (sp.Matrix([[1, 0], [0, -1]]), sp.Matrix([[0, 1], [0, 0]]), sp.Matrix([[0, 0], [1, 0]])):
    vx = [0] * 16
    for i in range(8):
        w = Xg * sp.Matrix(pt[i])
        vx[i], vx[8 + i] = w[0], w[1]
    fields.append(vx)
for i in range(7):  # torus {prod t = 1}: scale column i up, column 8 down
    vx = [0] * 16
    vx[i], vx[8 + i] = pt[i]
    vx[7], vx[15] = -pt[7][0], -pt[7][1]
    fields.append(vx)
Fm = sp.Matrix(fields).T
check("10 group directions (sl2 + 7-dim torus {prod t_i = 1}) lie in ker dPsi and are independent => dim X <= 6",
      (J * Fm) == sp.zeros(14, 10) and Fm.rank() == 10)

# ---------------------------------------------------------------- symmetry group action on V_C
def action_matrix(tau):
    """M with M Psi(v) = Psi(v o tau), (v o tau)_i = v_{tau(i)}; tau dict on 1..8"""
    # row a = coordinates of B_a(v o tau) = prod [tau(i) tau(j)] in the basis B
    rows_ = [expand_matching([(tau[i], tau[j]) for i, j in Bm]) for Bm in BASIS]
    return sp.Matrix(rows_)


def gens():
    out = []
    for rho in [(2, 1, 3, 4), (2, 3, 4, 1)]:
        out.append(("first block " + str(rho), {**{i + 1: rho[i] for i in range(4)}, **{4 + j: 4 + j for j in range(1, 5)}}))
        out.append(("second block " + str(rho), {**{i: i for i in range(1, 5)}, **{4 + j + 1: 4 + rho[j] for j in range(4)}}))
    out.append(("swap blocks", {**{i: 4 + i for i in range(1, 5)}, **{4 + i: i for i in range(1, 5)}}))
    return out


okM = True
okX = True
okP = True
okPhi = True
okFlip = True
vr = [(rng.randint(-5, 5), rng.randint(-5, 5)) for _ in range(8)]
for name, tau in gens():
    M = action_matrix(tau)
    okM &= (M * kappa == kappa) and all(e == int(e) for e in M)
    vt = [vr[tau[i] - 1] for i in range(1, 9)]
    okX &= (M * Psi(vr) - Psi(vt)).applyfunc(sp.expand) == sp.zeros(14, 1)
    # induced action on permutations: F_sigma(M x) = F_{h(sigma)}(x); find h by matching rows of P M and P
    PM = P * M
    hmap = {}
    for k in range(24):
        hits = [l for l in range(24) if PM.row(k) == P.row(l)]
        hmap[k] = hits[0] if len(hits) == 1 else None
    okP &= all(v is not None for v in hmap.values()) and len(set(hmap.values())) == 24
    if name == "swap blocks":
        okP &= all(PERMS[hmap[k]] == inverse(PERMS[k]) for k in range(24))
    # Leibniz covariance: Phi_{M c}(O') = Phi_c(O) for an explicit O' in O(4), and the det -1 column flip
    O = cayley(rng)
    cvec = sp.Matrix([sp.Rational(rng.randint(-9, 9), rng.randint(1, 4)) for _ in range(14)])
    # Phi_{Mc}(O') = sum_k (P c)_{h(k)} w_k(O'), so we need w_k(O') = w_{h(k)}(O)
    if name == "swap blocks":
        Op = O.T
    elif name.startswith("first"):
        r = tuple(tau[i] for i in range(1, 5))
        Op = sp.Matrix(4, 4, lambda i, j: O[r[i] - 1, j])
    else:
        r = tuple(tau[4 + j] - 4 for j in range(1, 5))
        Op = sp.Matrix(4, 4, lambda i, j: O[i, r[j] - 1])
    lhs = Phi(M * cvec, Op, extended=True)
    okPhi &= sp.simplify(lhs - Phi(cvec, O, extended=True)) == 0
    if Op.det() == -1:
        Of = Op * sp.diag(-1, 1, 1, 1)
        okFlip &= sp.simplify(Phi(M * cvec, Of) - lhs) == 0
check("S4 x S4 x| C2 generators act on V_C by integer matrices fixing kappa", okM)
check("the action is induced by relabelling columns, so it maps X onto X (M Psi(v) = Psi(v o tau))", okX)
check("each generator permutes the 24 rows of P (P M = Pi P); block swap sends F_sigma to F_{sigma^-1}", okP)
check("(4.7) Leibniz covariance on O(4) with the det-normalised extension", okPhi)
check("det -1 images: flipping one column returns to SO(4) with the same extended Leibniz sum", okFlip)

# C transforms by sign under adjacent transpositions (informational, App. A.1)
okC = True
for i in range(1, 8):
    tau = {j: j for j in range(1, 9)}
    tau[i], tau[i + 1] = i + 1, i
    M = action_matrix(tau)
    okC &= (M * M == sp.eye(14))
    msub = dict(zip(m, list(M * sp.Matrix(m))))
    okC &= sp.expand(C.subs(msub, simultaneous=True) + C) == 0
check("T_i^2 = I and C(T_i m) = -C(m) for the 7 adjacent transpositions (informational)", okC)

# ---------------------------------------------------------------- dump P
here = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(here, "check_section3_P.json"), "w") as fh:
    json.dump({"description": "P[sigma][a]: coefficient of B_a (eq. A.1) in prod_i [i,4+sigma(i)]; "
               "rows = permutations of 1..4 in lexicographic order (one-line), produced by "
               "Pluecker straightening and checked as a polynomial identity",
               "perms": ["".join(map(str, p)) for p in PERMS],
               "P": [[int(P[k, a]) for a in range(14)] for k in range(24)],
               "kappa": [int(kappa[a]) for a in range(14)]}, fh, indent=1)
print("wrote check_section3_P.json")

print()
print("SUMMARY:", "ALL PASS" if not FAIL else f"{len(FAIL)} FAILED: {FAIL}")
sys.exit(1 if FAIL else 0)
