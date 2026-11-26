"""Independent exact checks of the dimension-five functional Phi of task mdo-orthogonal-hull-five.

Phi, the unitary V and the claim r(V) not in conv r(O(5)) are due to agent ganesh-claude-minors
(GitHub owner g4nesh). This file re-implements everything from the formulas; it imports nothing
from that task. Standard library only, plus an optional sympy cross-check of Part C.

Usage:  python3 check_exact.py          (reads V_exact.json from this directory; about 10 seconds)

Part A  Phi on the 120 permutation matrices; the set S where it equals -5 and its structure.
Part B  The Laplace-Jacobi identity  sum_{sigma(I)=J} sgn(sigma) prod_i O[i,sigma(i)] = det(O) det(O[I,J])^2
        and its consequence Phi(O) = 1 - 6 det(O) E_S(O), checked exactly on rational orthogonal
        matrices of both determinants (Cayley transforms).
Part C  V is unitary and Phi(V) = 37/25, in exact arithmetic in Q(sqrt6, i sqrt15).
Part D  V is a critical point of Phi on U(5); the Hessian there is negative semidefinite with kernel
        exactly the 9-dimensional tangent space of the orbit {D1 V D2 : D1, D2 diagonal unitary},
        so V is a strict local maximum of Phi on U(5) modulo that orbit.
"""
import itertools
import json
import os
import random
from fractions import Fraction as Fr

N = 5
PERMS = list(itertools.permutations(range(N)))
ok_all = True


def report(name, cond):
    global ok_all
    ok_all = ok_all and bool(cond)
    print(("PASS  " if cond else "FAIL  ") + name)


# ---------------------------------------------------------------- the functional
def phi_terms():
    """Phi = sum of coef * |det M[I,J]|^2 over the returned (coef, I, J)."""
    t = []
    for v in range(N):
        t.append((-1, (v,), (v,)))
    for v in range(N):
        A = tuple(sorted({(v + 1) % N, (v - 1) % N}))
        B = tuple(sorted({(v + 2) % N, (v - 2) % N}))
        t += [(1, A, B), (1, B, A)]
    for v in range(N):
        for d in (1, 2):
            P = tuple(sorted({v, (v + d) % N}))
            Q = tuple(sorted({v, (v - d) % N}))
            t += [(-1, P, Q), (-1, Q, P)]
    return t


TERMS = phi_terms()


def sign(p):
    p, s = list(p), 1
    for i in range(len(p)):
        while p[i] != i:
            j = p[i]
            p[i], p[j] = p[j], p[i]
            s = -s
    return s


def cycle_type(p):
    seen, ct = set(), []
    for i in range(N):
        if i not in seen:
            L, j = 0, i
            while j not in seen:
                seen.add(j)
                j = p[j]
                L += 1
            ct.append(L)
    return tuple(sorted(ct, reverse=True))


def comp(a, b):
    return tuple(a[b[i]] for i in range(N))


def inv(a):
    r = [0] * N
    for i, x in enumerate(a):
        r[x] = i
    return tuple(r)


def maps(p, I, J):
    return tuple(sorted(p[i] for i in I)) == J


# ---------------------------------------------------------------- Part A
print("Part A: Phi on permutation matrices")
report("Phi has 5 + 10 + 20 = 35 terms with 35 distinct (I,J)", len(TERMS) == 35 and len({(I, J) for _, I, J in TERMS}) == 35)
VAL = {p: sum(c for c, I, J in TERMS if maps(p, I, J)) for p in PERMS}
S = sorted(p for p in PERMS if VAL[p] == -5)
SS = set(S)
report("Phi(P_sigma) takes only the values 1 and -5, each on 60 permutations",
       sorted(set(VAL.values())) == [-5, 1] and len(S) == 60)
AGL = {tuple((a * x + b) % N for x in range(N)) for a in range(1, N) for b in range(N)}
D5 = {tuple((a * x + b) % N for x in range(N)) for a in (1, N - 1) for b in range(N)}
relabel_ok = all(sorted((c, tuple(sorted(g[i] for i in I)), tuple(sorted(g[j] for j in J))) for c, I, J in TERMS) == sorted(TERMS) for g in AGL)
transp_ok = sorted((c, J, I) for c, I, J in TERMS) == sorted(TERMS)
report("Phi is invariant under M -> P_g M P_g^T for g in AGL(1,5), and under transposition", relabel_ok and transp_ok)
by_type = {}
for p in S:
    by_type[cycle_type(p)] = by_type.get(cycle_type(p), 0) + 1
print("      cycle types in S:", by_type)
transpositions = {p for p in PERMS if cycle_type(p) == (2, 1, 1, 1)}
three_cycles = {p for p in PERMS if cycle_type(p) == (3, 1, 1)}
four_cycles_nonaffine = {p for p in PERMS if cycle_type(p) == (4, 1) and p not in AGL}
report("S = D5 (x -> +-x+b) u {10 transpositions} u {20 three-cycles} u {20 four-cycles not in AGL(1,5)}",
       SS == D5 | transpositions | three_cycles | four_cycles_nonaffine)
report("identity and the transposition (0 1) lie in S", tuple(range(N)) in SS and (1, 0, 2, 3, 4) in SS)
report("S has 30 even and 30 odd permutations, so S is not A5", sum(sign(p) for p in S) == 0)
report("S is closed under inversion (Phi invariant under transposition)", all(inv(p) in SS for p in S))
report("S is not a subgroup (not closed under composition)", not all(comp(a, b) in SS for a in S for b in S))
normalizer = {g for g in PERMS if all(comp(comp(g, p), inv(g)) in SS for p in S)}
report("{g : g S g^-1 = S} is exactly AGL(1,5)", normalizer == AGL)
report("no g != e has gS = S or Sg = S (S is not a union of left or right cosets of a nontrivial subgroup)",
       [g for g in PERMS if all(comp(g, p) in SS for p in S)] == [tuple(range(N))]
       and [g for g in PERMS if all(comp(p, g) in SS for p in S)] == [tuple(range(N))])
dcs, seen = [], set()
for p in PERMS:
    if p not in seen:
        dc = {comp(comp(g, p), h) for g in AGL for h in AGL}
        seen |= dc
        dcs.append(dc)
report("AGL(1,5) has two double cosets (sizes 20, 100); S meets both and contains neither",
       sorted(len(d) for d in dcs) == [20, 100] and all(d & SS and not d <= SS for d in dcs))

# ---------------------------------------------------------------- exact linear algebra over Q
def matmul(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]


def det_q(A):
    A = [row[:] for row in A]
    n, d = len(A), Fr(1)
    for c in range(n):
        piv = next((r for r in range(c, n) if A[r][c] != 0), None)
        if piv is None:
            return Fr(0)
        if piv != c:
            A[c], A[piv] = A[piv], A[c]
            d = -d
        d *= A[c][c]
        for r in range(c + 1, n):
            f = A[r][c] / A[c][c]
            for k in range(c, n):
                A[r][k] -= f * A[c][k]
    return d


def inverse_q(A):
    n = len(A)
    M = [list(A[i]) + [Fr(int(i == j)) for j in range(n)] for i in range(n)]
    for c in range(n):
        piv = next(r for r in range(c, n) if M[r][c] != 0)
        M[c], M[piv] = M[piv], M[c]
        pv = M[c][c]
        M[c] = [x / pv for x in M[c]]
        for r in range(n):
            if r != c and M[r][c] != 0:
                f = M[r][c]
                M[r] = [x - f * y for x, y in zip(M[r], M[c])]
    return [row[n:] for row in M]


def minor_q(A, I, J):
    return det_q([[A[i][j] for j in J] for i in I]) if I else Fr(1)


# ---------------------------------------------------------------- Part B
print("Part B: the Laplace-Jacobi identity on rational orthogonal matrices")
rng = random.Random(20261125)
allIJ = [(I, J) for r in range(N + 1) for I in itertools.combinations(range(N), r) for J in itertools.combinations(range(N), r)]
ok_minor, ok_phi, ok_orth = True, True, True
for trial in range(3):
    K = [[Fr(0)] * N for _ in range(N)]
    for i in range(N):
        for j in range(i + 1, N):
            K[i][j] = Fr(rng.randint(-4, 4), rng.randint(1, 3))
            K[j][i] = -K[i][j]
    Id = [[Fr(int(i == j)) for j in range(N)] for i in range(N)]
    O1 = matmul([[Id[i][j] - K[i][j] for j in range(N)] for i in range(N)], inverse_q([[Id[i][j] + K[i][j] for j in range(N)] for i in range(N)]))
    O2 = [[-O1[i][j] if j == 0 else O1[i][j] for j in range(N)] for i in range(N)]   # determinant -1
    for O in (O1, O2):
        ok_orth &= matmul(O, [list(r) for r in zip(*O)]) == Id
        dO = det_q(O)
        prod = {}
        for p in PERMS:
            x = Fr(sign(p))
            for i in range(N):
                x *= O[i][p[i]]
            prod[p] = x
        for I, J in allIJ:
            lhs = sum(prod[p] for p in PERMS if maps(p, I, J))
            ok_minor &= lhs == dO * minor_q(O, I, J) ** 2
        E_S = sum(prod[p] for p in S)
        phi_O = sum(c * minor_q(O, I, J) ** 2 for c, I, J in TERMS)
        ok_phi &= phi_O == 1 - 6 * dO * E_S
report("6 rational orthogonal matrices (3 with det +1, 3 with det -1) are exactly orthogonal", ok_orth)
report("sum_{sigma(I)=J} sgn(sigma) prod_i O[i,sigma(i)] = det(O) det(O[I,J])^2 for all 252 (I,J), both components", ok_minor)
report("Phi(O) = 1 - 6 det(O) E_S(O), with E_S(O) = sum_{sigma in S} sgn(sigma) prod_i O[i,sigma(i)]", ok_phi)
print("      (general proof: Laplace expansion along the rows I gives the left side as eps*det O[I,J]*det O[I^c,J^c];")
print("       Jacobi's complementary-minor identity with O^-1 = O^T gives det O[I^c,J^c] = eps*det(O)*det O[I,J].)")


# ---------------------------------------------------------------- the field Q(s, t), s^2 = 6, t^2 = -15
class F:
    """a + b s + c t + d u with s = sqrt6, t = i sqrt15, u = s t;  st = u, su = 6t, tu = -15s, u^2 = -90."""
    __slots__ = ("a", "b", "c", "d")

    def __init__(self, a=0, b=0, c=0, d=0):
        self.a, self.b, self.c, self.d = Fr(a), Fr(b), Fr(c), Fr(d)

    def __add__(x, y):
        return F(x.a + y.a, x.b + y.b, x.c + y.c, x.d + y.d)

    def __sub__(x, y):
        return F(x.a - y.a, x.b - y.b, x.c - y.c, x.d - y.d)

    def __neg__(x):
        return F(-x.a, -x.b, -x.c, -x.d)

    def __mul__(x, y):
        if not isinstance(y, F):
            y = F(y)
        a, b, c, d, e, f, g, h = x.a, x.b, x.c, x.d, y.a, y.b, y.c, y.d
        return F(a * e + 6 * b * f - 15 * c * g - 90 * d * h,
                 a * f + b * e - 15 * c * h - 15 * d * g,
                 a * g + c * e + 6 * b * h + 6 * d * f,
                 a * h + d * e + b * g + c * f)

    __rmul__ = __mul__

    def conj(x):
        return F(x.a, x.b, -x.c, -x.d)

    def re(x):          # real part, an element of Q(sqrt6)
        return F(x.a, x.b)

    def im_over_t(x):   # imaginary part divided by sqrt15, as an element of Q(sqrt6): x = re + t * im_over_t
        return F(x.c, x.d)

    def is_zero(x):
        return x.a == 0 and x.b == 0 and x.c == 0 and x.d == 0

    def inv_real(x):    # inverse of a nonzero element of Q(sqrt6)
        assert x.c == 0 and x.d == 0
        n = x.a * x.a - 6 * x.b * x.b
        return F(x.a / n, -x.b / n)

    def sign_real(x):   # sign of a + b sqrt6
        assert x.c == 0 and x.d == 0
        a, b = x.a, x.b
        if a >= 0 and b >= 0:
            return 0 if a == 0 and b == 0 else 1
        if a <= 0 and b <= 0:
            return -1
        # opposite signs: compare a^2 with 6 b^2
        return (1 if a > 0 else -1) if a * a > 6 * b * b else (1 if b > 0 else -1)

    def __repr__(x):
        return "(%s + %s s + %s t + %s st)" % (x.a, x.b, x.c, x.d)


ZERO, ONE, T = F(), F(1), F(0, 0, 1)

here = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(here, "V_exact.json")) as fh:
    data = json.load(fh)
V = [[F(*[Fr(x) for x in e]) for e in row] for row in data["V"]]


def mm(A, B):
    return [[sum((A[i][k] * B[k][j] for k in range(N)), ZERO) for j in range(N)] for i in range(N)]


def adj(A):
    return [[A[j][i].conj() for j in range(N)] for i in range(N)]


def fminor(M, I, J):
    if len(I) == 1:
        return M[I[0]][J[0]]
    return M[I[0]][J[0]] * M[I[1]][J[1]] - M[I[0]][J[1]] * M[I[1]][J[0]]


def phi_F(M):
    tot = ZERO
    for c, I, J in TERMS:
        m = fminor(M, I, J)
        tot = tot + c * (m * m.conj())
    return tot


# ---------------------------------------------------------------- Part C
print("Part C: the unitary V (due to ganesh-claude-minors)")
VV = mm(V, adj(V))
report("V V^* = I exactly in Q(sqrt6, i sqrt15)",
       all((VV[i][j] - (ONE if i == j else ZERO)).is_zero() for i in range(N) for j in range(N)))
pv = phi_F(V)
report("Phi(V) = 37/25 exactly", pv.a == Fr(37, 25) and pv.b == 0 and pv.c == 0 and pv.d == 0)
sq = {}
for i in range(N):
    for j in range(N):
        m = V[i][j] * V[i][j].conj()
        sq[(m.a, m.b)] = sq.get((m.a, m.b), 0) + 1
report("|V_ij|^2 = 1/25 on the diagonal (5 entries) and 6/25 off it (20 entries)",
       sq == {(Fr(1, 25), Fr(0)): 5, (Fr(6, 25), Fr(0)): 20} and all((V[i][i] * V[i][i].conj()).a == Fr(1, 25) for i in range(N)))
print("      Phi(V) - 1 = 12/25 is the excess over the value 1 that Phi takes at the permutation matrices outside S.")

try:  # independent second implementation with sympy
    import sympy as sp
    s6, t15 = sp.sqrt(6), sp.I * sp.sqrt(15)
    Vs = sp.Matrix(N, N, lambda i, j: sum(sp.Rational(x) * b for x, b in zip(data["V"][i][j], (1, s6, t15, s6 * t15))))
    uni = sp.simplify(sp.expand(Vs * Vs.H - sp.eye(N))) == sp.zeros(N, N)
    ph = sum(c * sp.expand(sp.Abs(sp.expand(Vs[list(I), list(J)].det())) ** 2) for c, I, J in TERMS)
    ph = sp.nsimplify(sp.simplify(ph))
    report("sympy cross-check: V V^* = I and Phi(V) = 37/25", uni and ph == sp.Rational(37, 25))
except ImportError:
    print("      (sympy not installed; skipped the sympy cross-check)")

# ---------------------------------------------------------------- Part D
print("Part D: V is a critical point and a strict local maximum modulo phases")
# real basis of u(5) with entries in the field: E_jk - E_kj, t(E_jk + E_kj) (j<k), t E_jj
basis = []
for j in range(N):
    for k in range(j + 1, N):
        X = [[ZERO] * N for _ in range(N)]
        X[j][k], X[k][j] = ONE, -ONE
        basis.append(X)
for j in range(N):
    for k in range(j + 1, N):
        X = [[ZERO] * N for _ in range(N)]
        X[j][k], X[k][j] = T, T
        basis.append(X)
for j in range(N):
    X = [[ZERO] * N for _ in range(N)]
    X[j][j] = T
    basis.append(X)


def taylor2(X):
    """g(e) = Phi(exp(eX) V) = Phi(V) + e*g1 + e^2*g2 + O(e^3); returns (g1, g2) in Q(sqrt6)."""
    A = mm(X, V)
    B = mm(X, A)                       # X^2 V; exp(eX)V = V + e A + e^2 B/2 + ...
    g1, g2 = ZERO, ZERO
    for c, I, J in TERMS:
        p0 = fminor(V, I, J)
        if len(I) == 1:
            p1 = A[I[0]][J[0]]
            p2 = B[I[0]][J[0]] * Fr(1, 2)
        else:
            (i, k), (j, l) = I, J
            p1 = A[i][j] * V[k][l] + V[i][j] * A[k][l] - A[i][l] * V[k][j] - V[i][l] * A[k][j]
            p2 = (B[i][j] * V[k][l] + V[i][j] * B[k][l] - B[i][l] * V[k][j] - V[i][l] * B[k][j]) * Fr(1, 2) \
                + A[i][j] * A[k][l] - A[i][l] * A[k][j]
        g1 = g1 + c * 2 * (p0.conj() * p1).re()
        g2 = g2 + c * ((p1.conj() * p1).re() + 2 * (p0.conj() * p2).re())
    return g1, g2


report("first derivative of Phi(exp(eX)V) at e=0 vanishes for all 25 basis directions X of u(5)",
       all(taylor2(X)[0].is_zero() for X in basis))


def addm(X, Y):
    return [[X[i][j] + Y[i][j] for j in range(N)] for i in range(N)]


q = [taylor2(X)[1] for X in basis]
H = [[None] * 25 for _ in range(25)]          # H[a][b] = bilinear form with H[a][a] = g2(X_a)
for a in range(25):
    H[a][a] = q[a]
    for b in range(a + 1, 25):
        H[a][b] = H[b][a] = (taylor2(addm(basis[a], basis[b]))[1] - q[a] - q[b]) * Fr(1, 2)


def coords(X):
    """coordinates of a skew-Hermitian X in the basis above, as elements of Q(sqrt6)."""
    out = [X[j][k].re() for j in range(N) for k in range(j + 1, N)]
    out += [X[j][k].im_over_t() for j in range(N) for k in range(j + 1, N)]
    out += [X[j][j].im_over_t() for j in range(N)]
    return out


# tangent vectors of the orbit: left phases t E_jj (basis 20..24) and right phases V (t E_jj) V^*
kern = [coords(basis[20 + j]) for j in range(N)]
for j in range(N):
    D = [[ZERO] * N for _ in range(N)]
    D[j][j] = T
    kern.append(coords(mm(mm(V, D), adj(V))))


def Hv(vec):
    return [sum((H[a][b] * vec[b] for b in range(25)), ZERO) for a in range(25)]


report("the Hessian kills all 10 orbit tangent vectors (5 left and 5 right phase directions)",
       all(all(x.is_zero() for x in Hv(k)) for k in kern))


def rank_real(rows):
    rows = [r[:] for r in rows]
    rk, col, n = 0, 0, len(rows[0])
    while rk < len(rows) and col < n:
        piv = next((r for r in range(rk, len(rows)) if not rows[r][col].is_zero()), None)
        if piv is None:
            col += 1
            continue
        rows[rk], rows[piv] = rows[piv], rows[rk]
        iv = rows[rk][col].inv_real()
        for r in range(len(rows)):
            if r != rk and not rows[r][col].is_zero():
                f = rows[r][col] * iv
                rows[r] = [x - f * y for x, y in zip(rows[r], rows[rk])]
        rk += 1
        col += 1
    return rk


report("the orbit tangent space has dimension 9", rank_real(kern) == 9)
# complement: 16 of the 20 off-diagonal basis vectors, chosen so that with the orbit tangents they span u(5)
comp_idx = None
for drop in itertools.combinations(range(20), 4):
    W = [a for a in range(20) if a not in drop]
    rows = kern + [[ONE if b == a else ZERO for b in range(25)] for a in W]
    if rank_real(rows) == 25:
        comp_idx = W
        break
report("found 16 basis directions complementing the orbit tangent space", comp_idx is not None)
Hc = [[H[a][b] for b in comp_idx] for a in comp_idx]
piv_signs = []
M = [r[:] for r in Hc]
for c in range(16):
    p = M[c][c]
    piv_signs.append(p.sign_real())
    if p.is_zero():
        break
    ip = p.inv_real()
    for r in range(c + 1, 16):
        f = M[r][c] * ip
        M[r] = [x - f * y for x, y in zip(M[r], M[c])]
report("Hessian on that complement is negative definite (16 negative pivots in exact LDL^T)",
       len(piv_signs) == 16 and all(s == -1 for s in piv_signs))
print("      Hence Phi is constant on the 9-dimensional orbit of V and strictly smaller near V off it:")
print("      V is a strict local maximum of Phi on U(5) modulo phases.  (Global maximality is NOT checked here.)")

print()
print("ALL EXACT CHECKS PASSED" if ok_all else "SOME CHECK FAILED")
