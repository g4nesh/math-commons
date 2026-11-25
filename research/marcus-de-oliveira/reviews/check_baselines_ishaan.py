#!/usr/bin/env python3
"""Exact Q(i) re-check of baselines/README.md (commit a0537fd): forest theorem and triangle.

Review aid for mdo-baselines, written by ishaan-air. It shares no code with
baselines/forest.py or certificates/verify.py. Standard library only; deterministic.

Checks:
  (1)   d(P_p) = z_{p^-1}, using the objective's z_sigma = prod_i (a_i + b_sigma(i))
  (2)   single arbitrary U(2) gate: block determinant formula
  (7)   determinant = sum |det U[S,T]|^2 a^(S^c) b^T, on non-forest unitaries
  (5)   every squared minor of a forest product equals the independent-transposition mass,
        with each gate's magnitude q drawn independently from seven exact values
  (6,8) determinant identity, also with monomial end factors and pushed-forward weights
  zeros every forest product with n >= 3 has a zero entry (scope note in the review)
  §4    triangle numbers; a Birkhoff decomposition of |U_ij|^2 gives other weights that
        match every squared minor of the same U
  ctrl  the product weights fail for a repeated edge and for the triangle
  scope the 4x4 matrix O_* quoted in literature/README.md has psi(O_*) = -158/2401 < 0 while
        psi(P) >= 0 for all 24 permutation matrices, so no spectrum-independent weights
        match all of its squared minors

Bounds: n in 1..6, at most 5 gates, 150 forest cases, seed 20261125.
Usage: python3 -I research/marcus-de-oliveira/reviews/check_baselines_ishaan.py
"""
from fractions import Fraction as F
from itertools import combinations, permutations, product
import random

SEED = 20261125
rng = random.Random(SEED)
Z, ONE = (F(0), F(0)), (F(1), F(0))


def add(x, y): return (x[0] + y[0], x[1] + y[1])
def sub(x, y): return (x[0] - y[0], x[1] - y[1])
def mul(x, y): return (x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0])
def conj(x): return (x[0], -x[1])
def abs2(x): return x[0] * x[0] + x[1] * x[1]
def inv(x): n = abs2(x); return (x[0] / n, -x[1] / n)
def scal(f): return (F(f), F(0))


def prod(xs):
    out = ONE
    for x in xs:
        out = mul(out, x)
    return out


def matmul(A, B):
    return [[prod_sum(A[i], [B[k][j] for k in range(len(B))]) for j in range(len(B[0]))]
            for i in range(len(A))]


def prod_sum(row, col):
    out = Z
    for x, y in zip(row, col):
        out = add(out, mul(x, y))
    return out


def adj(A): return [[conj(A[j][i]) for j in range(len(A))] for i in range(len(A[0]))]
def eye(n): return [[ONE if i == j else Z for j in range(n)] for i in range(n)]
def diag(v): return [[v[i] if i == j else Z for j in range(len(v))] for i in range(len(v))]


def det(M):
    M, n, d = [r[:] for r in M], len(M), ONE
    for c in range(n):
        p = next((r for r in range(c, n) if M[r][c] != Z), None)
        if p is None:
            return Z
        if p != c:
            M[c], M[p], d = M[p], M[c], (-d[0], -d[1])
        d, pinv = mul(d, M[c][c]), inv(M[c][c])
        for r in range(c + 1, n):
            f = mul(M[r][c], pinv)
            for k in range(c, n):
                M[r][k] = sub(M[r][k], mul(f, M[c][k]))
    return d


def det_leibniz(M):
    n, out = len(M), Z
    for p in permutations(range(n)):
        sign = (-1) ** sum(p[i] > p[j] for i in range(n) for j in range(i + 1, n))
        t = prod(M[i][p[i]] for i in range(n))
        out = add(out, (sign * t[0], sign * t[1]))
    return out


def d_ab(a, b, U):
    n = len(a)
    M = matmul(matmul(U, diag(b)), adj(U))
    return [[add(M[i][j], a[i]) if i == j else M[i][j] for j in range(n)] for i in range(n)]


def z_sigma(a, b, sigma):  # objective convention: row i receives b_sigma(i)
    return prod(add(a[i], b[sigma[i]]) for i in range(len(a)))


def inverse(p):
    q = [0] * len(p)
    for j, i in enumerate(p):
        q[i] = j
    return tuple(q)


def perm_matrix(p, phases=None):  # column-to-row: P e_j = phase_j e_p(j)
    n = len(p)
    P = [[Z] * n for _ in range(n)]
    for j in range(n):
        P[p[j]][j] = phases[j] if phases else ONE
    return P


PYTH = [(F(3, 5), F(4, 5)), (F(5, 13), F(12, 13)), (F(8, 17), F(15, 17)), (F(7, 25), F(24, 25)),
        (F(20, 29), F(21, 29)), (F(1), F(0)), (F(0), F(1))]
PHASES = [ONE, (F(0), F(1)), (F(3, 5), F(4, 5)), (F(5, 13), F(-12, 13)), (F(-8, 17), F(15, 17)),
          (F(-1), F(0))]


def gate(c, s):
    ph = [rng.choice(PHASES) for _ in range(4)]
    g = matmul(matmul(diag(ph[:2]), [[scal(c), scal(s)], [scal(-s), scal(c)]]), diag(ph[2:]))
    assert matmul(g, adj(g)) == eye(2) and matmul(adj(g), g) == eye(2)
    return g


def embed(n, i, j, g):
    G = eye(n)
    for r, R in enumerate((i, j)):
        for c, C in enumerate((i, j)):
            G[R][C] = g[r][c]
    return G


def product_of(n, gates):  # gates in application order: V = G_m ... G_1
    V = eye(n)
    for i, j, g in gates:
        V = matmul(embed(n, i, j, g), V)
    return V


def product_weights(n, gates):  # enumerate eps in {0,1}^m, compose transpositions explicitly
    out = {}
    for eps in product((0, 1), repeat=len(gates)):
        p, w = list(range(n)), F(1)
        for e, (i, j, g) in zip(eps, gates):
            q = abs2(g[0][1])
            w *= q if e else 1 - q
            if e:
                p = [j if x == i else i if x == j else x for x in p]
        out[tuple(p)] = out.get(tuple(p), F(0)) + w
    return out


def check_minors(U, weights):
    n, count = len(U), 0
    for k in range(n + 1):
        for S in combinations(range(n), k):
            for T in combinations(range(n), k):
                lhs = abs2(det([[U[i][j] for j in T] for i in S]))
                rhs = sum((w for p, w in weights.items() if set(p[j] for j in T) == set(S)), F(0))
                assert lhs == rhs, (S, T, lhs, rhs)
                count += 1
    return count


def expansion(a, b, U):  # right side of (7)
    n, out = len(a), Z
    for k in range(n + 1):
        for S in combinations(range(n), k):
            for T in combinations(range(n), k):
                m = abs2(det([[U[i][j] for j in T] for i in S]))
                if m:
                    t = mul(prod(a[i] for i in range(n) if i not in S), prod(b[j] for j in T))
                    out = add(out, (m * t[0], m * t[1]))
    return out


def spectrum(n):
    pool = [Z, Z, ONE] + [(F(rng.randint(-3, 3)), F(rng.randint(-3, 3))) for _ in range(4)]
    return [rng.choice(pool) for _ in range(n)]


def random_forest(n):
    order = list(range(n))
    rng.shuffle(order)
    edges = [(order[v], order[rng.randrange(v)]) for v in range(1, n)]
    edges = [e for e in edges if rng.random() < 0.8][:5]
    rng.shuffle(edges)
    return edges


def tree_distance(n, edges, i, j):
    nbr = {v: set() for v in range(n)}
    for x, y in edges:
        nbr[x].add(y), nbr[y].add(x)
    seen, frontier, dist = {i}, {i}, 0
    while frontier:
        if j in frontier:
            return dist
        frontier = {y for x in frontier for y in nbr[x]} - seen
        seen |= frontier
        dist += 1
    return None


def birkhoff(B):  # exact greedy Birkhoff decomposition of a doubly stochastic matrix
    n, B, out = len(B), [r[:] for r in B], {}
    while any(any(r) for r in B):
        p = next(p for p in permutations(range(n)) if all(B[p[j]][j] > 0 for j in range(n)))
        w = min(B[p[j]][j] for j in range(n))
        for j in range(n):
            B[p[j]][j] -= w
        out[p] = out.get(p, F(0)) + w
    return out


def main():
    # (1), (2), (7) on small instances
    for n in range(1, 5):
        for _ in range(10):
            a, b = spectrum(n), spectrum(n)
            p = tuple(rng.sample(range(n), n))
            assert det(d_ab(a, b, perm_matrix(p))) == z_sigma(a, b, inverse(p))
    for _ in range(40):
        n = rng.randint(2, 5)
        i, j = rng.sample(range(n), 2)
        c, s = rng.choice(PYTH)
        g, a, b = gate(c, s), spectrum(n), spectrum(n)
        t = abs2(g[0][0])
        two = add(mul(scal(t), mul(add(a[i], b[i]), add(a[j], b[j]))),
                  mul(scal(1 - t), mul(add(a[i], b[j]), add(a[j], b[i]))))
        rest = prod(add(a[k], b[k]) for k in range(n) if k not in (i, j))
        assert det(d_ab(a, b, embed(n, i, j, g))) == mul(two, rest)
    tri = [(0, 1), (1, 2), (0, 2)]
    for _ in range(10):  # (7) on triangle and repeated-edge products (not forests)
        gates = [(i, j, gate(*rng.choice(PYTH))) for i, j in tri + [(0, 1)]]
        U, a, b = product_of(3, gates), spectrum(3), spectrum(3)
        assert det(d_ab(a, b, U)) == expansion(a, b, U) == det_leibniz(d_ab(a, b, U))
    # (5), (6), (8) and forced zeros on forest products
    cases = minors = 0
    for n, reps in ((1, 5), (2, 20), (3, 35), (4, 40), (5, 35), (6, 15)):
        for _ in range(reps):
            edges = random_forest(n)
            gates = []
            for x, y in edges:
                i, j = (x, y) if rng.random() < 0.5 else (y, x)
                gates.append((i, j, gate(*rng.choice(PYTH))))
            V, W = product_of(n, gates), product_weights(n, gates)
            assert matmul(V, adj(V)) == eye(n) and sum(W.values()) == 1
            assert all(w >= 0 for w in W.values())
            minors += check_minors(V, W)
            if n >= 3:
                assert any(V[i][j] == Z for i in range(n) for j in range(n))
                for i in range(n):
                    for j in range(i + 1, n):
                        dist = tree_distance(n, edges, i, j)
                        if dist is not None and dist >= 2:
                            assert V[i][j] == Z or V[j][i] == Z
            l, r = tuple(rng.sample(range(n), n)), tuple(rng.sample(range(n), n))
            U = matmul(matmul(perm_matrix(l, [rng.choice(PHASES) for _ in range(n)]), V),
                       perm_matrix(r, [rng.choice(PHASES) for _ in range(n)]))
            Wp = {}
            for p, w in W.items():
                key = tuple(l[p[r[j]]] for j in range(n))
                Wp[key] = Wp.get(key, F(0)) + w
            minors += check_minors(U, Wp)
            a, b = spectrum(n), spectrum(n)
            for M, weights in ((V, W), (U, Wp)):
                rhs = Z
                for p, w in weights.items():
                    zv = z_sigma(a, b, inverse(p))
                    rhs = add(rhs, (w * zv[0], w * zv[1]))
                lhs = d_ab(a, b, M)
                assert det(lhs) == rhs
                if n <= 5:
                    assert det_leibniz(lhs) == rhs
            cases += 1
    # §4 triangle obstruction and a different set of weights for the same U
    c, s = F(3, 5), F(4, 5)
    R = [[scal(c), scal(s)], [scal(-s), scal(c)]]
    gates = [(0, 1, R), (1, 2, R), (0, 2, R)]  # U = R13 R23 R12
    U, W = product_of(3, gates), product_weights(3, gates)
    ret = sum((w for p, w in W.items() if p[0] == 0), F(0))
    assert U[0][0] == scal(F(109, 125)) and abs2(U[0][0]) == F(11881, 15625)
    assert ret == c ** 4 + s ** 6 == F(6121, 15625)
    assert abs2(U[0][0]) - ret == 2 * c * c * s ** 3 == F(1152, 3125)
    a, b = [Z, ONE, ONE], [ONE, Z, Z]
    assert det(d_ab(a, b, U)) == scal(F(11881, 15625))
    assert {z_sigma(a, b, sg) for sg in permutations(range(3))} == {Z, ONE}
    try:
        check_minors(U, W)
        raise SystemExit("product weights unexpectedly matched the triangle")
    except AssertionError:
        pass
    other = birkhoff([[abs2(U[i][j]) for j in range(3)] for i in range(3)])
    tri_minors = check_minors(U, other)
    for cs in PYTH[:5]:  # U_11 = c^2 + s^3 for every Pythagorean pair tried
        Rc = [[scal(cs[0]), scal(cs[1])], [scal(-cs[1]), scal(cs[0])]]
        Uc = product_of(3, [(0, 1, Rc), (1, 2, Rc), (0, 2, Rc)])
        assert Uc[0][0] == scal(cs[0] ** 2 + cs[1] ** 3)
    # control: a repeated edge also breaks the product weights
    g = gate(F(3, 5), F(4, 5))
    try:
        check_minors(product_of(2, [(0, 1, g), (0, 1, g)]), product_weights(2, [(0, 1, g), (0, 1, g)]))
        raise SystemExit("product weights unexpectedly matched the repeated edge")
    except AssertionError:
        pass
    # scope: literature/README.md's O_* (Zheng, Ex. 1.3) has no spectrum-independent weights
    O = [[scal(F(x, 7)) for x in row] for row in
         ((1, 4, 4, 4), (-4, 1, 4, -4), (-4, -4, 1, 4), (-4, 4, -4, 1))]
    assert matmul(adj(O), O) == eye(4) and det(O) == ONE

    def psi2(V):  # 2*psi; index sets are 1-based in the source, 0-based here
        m = lambda S, T: abs2(det([[V[i][j] for j in T] for i in S]))
        return (sum(abs2(V[i][i]) for i in range(4)) - m((0, 1), (0, 2)) - m((0, 1), (1, 3))
                + m((0, 2), (0, 1)) + m((0, 2), (2, 3)))
    assert psi2(O) / 2 == F(-158, 2401)
    assert all(psi2(perm_matrix(p)) >= 0 for p in permutations(range(4)))
    print(f"seed {SEED}: {cases} forest cases, {minors} squared minors, all exact checks passed; "
          f"triangle: U11=109/125, |U11|^2=11881/15625, product weights 6121/15625, "
          f"Birkhoff weights on {len(other)} permutations match all {tri_minors} squared minors; "
          f"O_*: psi=-158/2401 < 0 <= psi(P) for all 24 permutations")


if __name__ == "__main__":
    main()
