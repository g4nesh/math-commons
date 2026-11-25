#!/usr/bin/env python3
"""Independent exact recheck of the mdo-adversarial witnesses (reviewer: ishaan-air).

Standard library only; exact Gaussian-rational arithmetic with fractions.Fraction.
Inputs are typed from Kovacec, DMUC p2611 v2, section 4 (not read from the fixtures);
the fixtures are only compared against them afterwards.

Checks, per witness (Goyal; Zhang-Ma):
  1. U^T U = I (U is real), D_a + U D_b U^* equals the matrix printed in the source,
     and det equals the printed value.
  2. sigma-points in lexicographic order equal the printed coordinates.
  3. p = det lies in conv{z_sigma}: explicit Caratheodory triangle with exact
     nonnegative barycentric coordinates (no facet enumeration).
  4. p lies on no transposition edge, and closed rays from p in several lattice
     directions meet no closed edge (so p is an exterior point in the sense of
     the source's Definition 6b, hence outside cl J(A,B)).
  5. the printed inflex points are crossings of two transposition edges.
Source formulas (p.3-4): the Lemma's (i,j) entry, the Corollary determinant
formula, and the Corollary 3 outside-product example (343 vs 441).
Run from the repository root:
  python3 -I research/marcus-de-oliveira/reviews/check_adversarial_ishaan.py
"""
import json
import sys
from fractions import Fraction as Q
from itertools import combinations, permutations
from pathlib import Path

FIX = Path("research/marcus-de-oliveira/certificates/fixtures")
fails = []


def ok(cond, msg):
    print(("PASS " if cond else "FAIL ") + msg)
    if not cond:
        fails.append(msg)


def C(x, y=0):
    return (Q(x), Q(y))


def add(x, y): return (x[0] + y[0], x[1] + y[1])
def sub(x, y): return (x[0] - y[0], x[1] - y[1])
def mul(x, y): return (x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0])
def conj(x): return (x[0], -x[1])
def scal(q, x): return (q * x[0], q * x[1])


def det(m):  # Laplace expansion; n <= 4 here
    if len(m) == 1:
        return m[0][0]
    tot = C(0)
    for j in range(len(m)):
        minor = [row[:j] + row[j + 1:] for row in m[1:]]
        term = mul(m[0][j], det(minor))
        tot = add(tot, term) if j % 2 == 0 else sub(tot, term)
    return tot


def matmul(A, B):
    n = len(A)
    return [[sum_c(mul(A[i][k], B[k][j]) for k in range(n)) for j in range(n)] for i in range(n)]


def sum_c(it):
    s = C(0)
    for x in it:
        s = add(s, x)
    return s


def adj(A):
    return [[conj(A[j][i]) for j in range(len(A))] for i in range(len(A))]


def diag(v):
    n = len(v)
    return [[v[i] if i == j else C(0) for j in range(n)] for i in range(n)]


def cross(o, a, b):
    return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])


def on_seg(p, a, b):
    return cross(a, b, p) == 0 and min(a[0], b[0]) <= p[0] <= max(a[0], b[0]) \
        and min(a[1], b[1]) <= p[1] <= max(a[1], b[1])


def ray_meets(p, d, a, b):
    """closed ray p + t d (t >= 0) vs closed segment [a, b], exact."""
    e = (b[0] - a[0], b[1] - a[1])
    w = (a[0] - p[0], a[1] - p[1])
    den = d[0] * e[1] - d[1] * e[0]
    if den == 0:  # parallel
        if w[0] * d[1] - w[1] * d[0] != 0:
            return False
        return max((q[0] - p[0]) * d[0] + (q[1] - p[1]) * d[1] for q in (a, b)) >= 0
    t = (w[0] * e[1] - w[1] * e[0]) / den
    s = (w[0] * d[1] - w[1] * d[0]) / den
    return t >= 0 and 0 <= s <= 1


def crossing(a, b, c, d):
    r, s = (b[0] - a[0], b[1] - a[1]), (d[0] - c[0], d[1] - c[1])
    den = r[0] * s[1] - r[1] * s[0]
    if den == 0:
        return None
    w = (c[0] - a[0], c[1] - a[1])
    t, u = (w[0] * s[1] - w[1] * s[0]) / den, (w[0] * r[1] - w[1] * r[0]) / den
    return (a[0] + t * r[0], a[1] + t * r[1]) if 0 <= t <= 1 and 0 <= u <= 1 else None


def witness(name, fixture, a, b, scale, Uint, printed_mat, printed_det, printed_z, rays, inflex):
    print(f"== {name}")
    n = len(a)
    U = [[C(Q(x, scale)) for x in row] for row in Uint]
    fx = json.loads((FIX / fixture).read_text(encoding="utf-8"))
    same = [tuple(map(Q, v)) for v in fx["a"]] == a and [tuple(map(Q, v)) for v in fx["b"]] == b \
        and [[tuple(map(Q, v)) for v in row] for row in fx["U"]] == U
    ok(same, "fixture a, b, U equal the inputs typed from p2611 v2 p.8")
    I = diag([C(1)] * n)
    ok(matmul(adj(U), U) == I, "U^* U = I exactly")
    M = [[add(diag(a)[i][j], matmul(matmul(U, diag(b)), adj(U))[i][j]) for j in range(n)] for i in range(n)]
    pm, pd = printed_mat
    ok(M == [[scal(Q(1, pd), C(*e)) for e in row] for row in pm], f"D_a + U D_b U^* = (1/{pd}) x printed matrix")
    p = det(M)
    ok(p == printed_det, f"det = {p[0]} + {p[1]} i, as printed")
    perms = list(permutations(range(n)))  # lexicographic, as in the source
    z = []
    for s in perms:
        v = C(1)
        for i in range(n):
            v = mul(v, add(a[i], b[s[i]]))
        z.append(v)
    ok(z == [C(*q) for q in printed_z], "sigma-points z1..z6 equal the printed coordinates")
    # 3. Caratheodory triangle with exact barycentric coordinates
    tri = None
    for i, j, k in combinations(range(len(z)), 3):
        D = cross(z[i], z[j], z[k])
        if D == 0:
            continue
        l1, l2, l3 = cross(p, z[j], z[k]) / D, cross(z[i], p, z[k]) / D, cross(z[i], z[j], p) / D
        if min(l1, l2, l3) >= 0:
            tri = (i + 1, j + 1, k + 1, l1, l2, l3)
            break
    ok(tri is not None, f"p in conv(z): triangle z{tri[0]},z{tri[1]},z{tri[2]} weights {tri[3]}, {tri[4]}, {tri[5]}"
       if tri else "p in conv(z)")
    # 4. edges = pairs of permutations differing in exactly two positions
    edges = [(s, t) for s, t in combinations(range(len(perms)), 2)
             if sum(x != y for x, y in zip(perms[s], perms[t])) == 2]
    ok(len(edges) == n * (n - 1) // 2 * len(perms) // 2, f"{len(edges)} transposition edges (Lemma 5 count)")
    ok(not any(on_seg(p, z[s], z[t]) for s, t in edges), "p lies on no edge")
    for d in rays:
        dd = (Q(d[0]), Q(d[1]))
        ok(not any(ray_meets(p, dd, z[s], z[t]) for s, t in edges), f"closed ray p + t{d} meets no closed edge")
    xs = {}
    for (s1, t1), (s2, t2) in combinations(edges, 2):
        q = crossing(z[s1], z[t1], z[s2], z[t2])
        if q is not None and q not in z:
            xs.setdefault(q, []).append(f"[z{s1+1},z{t1+1}]x[z{s2+1},z{t2+1}]")
    for q in inflex:
        ok(q in xs, f"inflex point {q[0]} + {q[1]} i is the crossing {xs.get(q)}")
    return z, p, edges


# Zhang-Ma (p2611 v2 p.8)
zm = witness("Zhang-Ma", "zhang-ma-2026.json",
             [C(-3, -1), C(-5, -2), C(4, -3)], [C(-2, 1), C(4, 3), C(-2, -5)], 7,
             [[-3, -2, 6], [-2, -6, -3], [-6, 3, -2]],
             ([[(-221, -208), (72, 132), (-36, 60)], [(72, 132), (-127, -31), (-108, -72)],
               [(-36, 60), (-108, -72), (152, -104)]], 49),
             C(Q(2864, 49), Q(1840, 49)),
             [(-30, -50), (280, 280), (-130, 10), (-28, -56), (232, 376), (24, -20)],
             [(0, -1), (1, 0), (1, -1), (1, 1)], [(Q(7672, 155), Q(4424, 155))])
z, p, _ = zm
cr = cross(z[1], z[3], p)
ok(cr != 0, f"Zhang-Ma: p is off the line z2 z4 (cross product {cr}); source says it 'seems to lie on [z2,z4]'")

# Goyal (p2611 v2 p.8)
witness("Goyal", "goyal-2026.json",
        [C(0), C(1), C(0, 1)], [C(-1, -1), C(1), C(-1, 1)], 3,
        [[-1, 2, 2], [2, -1, 2], [2, 2, -1]],
        ([[(-1, 3), (-4, 6), (8, 0)], [(-4, 6), (2, 0), (-4, -6)], [(8, 0), (-4, -6), (-1, 6)]], 9),
        C(Q(2, 3), Q(2, 9)),
        [(6, -2), (2, 0), (2, 1), (0, -1), (0, 2), (2, -2)],
        [(-1, -1), (-1, 0), (-1, 1)], [(Q(1), Q(0))])

# Source formulas, p.3-4. Rational point of S^3 with all coordinates distinct and nonzero.
print("== p2611 v2 Lemma / Corollary / Corollary 3")
c, bb, s, d = Q(2, 5), Q(1, 5), Q(2, 5), Q(4, 5)
ok(c * c + bb * bb + s * s + d * d == 1, "(c,b,s,d) = (2,1,2,4)/5 lies on S^3")
n, i, j = 4, 0, 2
G = diag([C(1)] * n)
G[i][i], G[i][j], G[j][i], G[j][j] = C(c, bb), C(s, d), C(-s, d), C(c, -bb)
bv = [C(1, 2), C(-3, 1), C(2, -5), C(7, 3)]
av = [C(2, -1), C(0, 4), C(-1, -1), C(5, 2)]
H = matmul(matmul(G, diag(bv)), adj(G))
true_ij = mul(mul(C(c, bb), C(s, d)), sub(bv[j], bv[i]))
proof_ij = add(mul(mul(C(c, bb), C(-s, -d)), bv[i]), mul(mul(C(s, d), C(c, -bb)), bv[j]))
display_ij = mul(mul(C(bb, -c), C(d, -s)), sub(bv[i], bv[j]))
ok(H[i][j] == true_ij, "(G D_b G^*)_{ij} = (c+ib)(s+id)(b_j-b_i)")
ok(H[i][j] == display_ij, "the displayed Lemma matrix entry (b-ic)(d-is)(b_i-b_j) is correct")
ok(H[i][j] != proof_ij, "proof subcase k=i,l=j as printed, (c+ib)(-s-id)b_i+(s+id)(c-ib)b_j, is NOT the entry (expected)")
lhs = det([[add(diag(av)[r][k], H[r][k]) for k in range(n)] for r in range(n)])
Mo = mul(add(av[1], bv[1]), add(av[3], bv[3]))
rhs = mul(Mo, add(mul(add(av[i], bv[i]), add(av[j], bv[j])),
                  scal(d * d + s * s, mul(sub(av[i], av[j]), sub(bv[i], bv[j])))))
ok(lhs == rhs, "Corollary: det(D_a + G D_b G^*) = M((a_i+b_i)(a_j+b_j)+(a_i-a_j)(b_i-b_j)(d^2+s^2)) exactly")
a3, b3, pi = [1, 2, 3], [4, 5, 6], [2, 1, 0]  # pi = (3,2,1) in one-line form, i=1, j=2, G = I
true3 = (a3[0] + b3[pi[0]]) * (a3[1] + b3[pi[1]]) * (a3[2] + b3[pi[2]])
printed3 = (a3[2] + b3[2]) * (a3[0] + b3[pi[0]]) * (a3[1] + b3[pi[1]])
ok((true3, printed3) == (343, 441), f"Corollary 3 example: true {true3}, printed outside product gives {printed3}")

print(f"{'ALL CHECKS PASS' if not fails else str(len(fails)) + ' FAILURES'}")
sys.exit(1 if fails else 0)
