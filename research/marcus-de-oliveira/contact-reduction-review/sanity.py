"""Independent finite boundary checks for the contact-reduction review.

Python standard library only; no imports from the author's checker.
These fixtures do not prove the all-dimension or analytic statements.
"""

from fractions import Fraction as F
from itertools import permutations
from math import factorial, prod


def bracket(v, w):
    return v[0] * w[1] - v[1] * w[0]


def phi(v, w, p):
    return prod(bracket(v[i], w[p[i]]) for i in range(len(p)))


def parity(p):
    return (-1) ** sum(p[i] > p[j] for i in range(len(p)) for j in range(i + 1, len(p)))


def slice_recovery(n):
    s = F(7, 3)
    ps = {i: F(-i - 3) for i in range(2, n)}
    qs = {j: F(j + 3) for j in range(1, n)}
    v = [(0, s), (1, 0)] + [(ps[i], 1) for i in range(2, n)]
    w = [(1, 1)] + [(qs[j], 1) for j in range(1, n)]
    for j in range(1, n):
        p = [0, j] + [k for k in range(1, n) if k != j]
        t = p.copy()
        t[0], t[1] = t[1], t[0]
        assert phi(v, w, t) / phi(v, w, p) == qs[j]
    for i in range(2, n):
        p = [None] * n
        p[0], p[i] = 0, 1
        remaining = iter(range(2, n))
        for k in range(n):
            if p[k] is None:
                p[k] = next(remaining)
        t = p.copy()
        t[0], t[i] = t[i], t[0]
        rho = phi(v, w, p) / phi(v, w, t)
        assert qs[1] * (rho - 1) / (rho * qs[1] - 1) == ps[i]
    assert -phi(v, w, range(n)) / prod(ps[i] - qs[i] for i in range(2, n)) == s


def singular_line_point(n):
    # The author's affine chart formula divides by Q+t. Check its excluded
    # point t=-Q directly in homogeneous coordinates, without taking a limit.
    c, delta = F(2), F(5, 3)
    b = [F(j + 1) for j in range(n)]
    q = prod(c + bj for bj in b)
    a = [c] * (n - 1) + [c + delta]
    v = [(c, 1)] * (n - 1) + [(delta, 0)]
    w = [(-bj, 1) for bj in b]
    for p in permutations(range(n)):
        z = prod(a[i] + b[p[i]] for i in range(n))
        assert phi(v, w, p) == z - q


def balanced_two_line_boundary(n, k):
    # n vectors on each of two perpendicular lines, equal norms, moment n*I.
    # Includes k=0,n as well as the semistable collision boundary.
    e, f = (F(1), F(0)), (F(0), F(1))
    v = [e] * k + [f] * (n - k)
    w = [e] * (n - k) + [f] * k
    vectors = v + w
    assert sum(u[0] ** 2 for u in vectors) == n
    assert sum(u[1] ** 2 for u in vectors) == n
    assert sum(u[0] * u[1] for u in vectors) == 0
    values = [phi(v, w, p) for p in permutations(range(n))]
    nonzero = [x for x in values if x]
    assert len(nonzero) == factorial(k) * factorial(n - k)
    assert set(nonzero) == {(-1) ** (n - k)}


def signed_admissible_line(n):
    pp = list(permutations(range(n)))
    mu = [F(1, factorial(n)) + 2 * parity(p) for p in pp]
    assert min(mu) < 0 and sum(mu) == 1
    for i in range(n):
        for j in range(n):
            assert sum(m for p, m in zip(pp, mu) if p[i] == j) == F(1, n)
    # Arbitrary nonconstant assignment-cost vector lies in im(A^T).
    costs = [[F((i + 2) ** (j + 1), i + j + 1) for j in range(n)] for i in range(n)]
    x = [sum(costs[i][p[i]] for i in range(n)) for p in pp]
    d = sum(m * y for m, y in zip(mu, x))
    assert d == sum(x) / factorial(n)
    assert min(x) <= d <= max(x)


def punctured_curve(t):
    # z=(1-t^2+2it)/(1+t^2) on C*, H1=(z+1/z)/2, H2=1-H1.
    # For 0<=t<=1/2 the exact arc reaches the finite contact at z=1.
    re, im = (1 - t * t) / (1 + t * t), 2 * t / (1 + t * t)
    assert re * re + im * im == 1
    h1, h2 = re, 1 - re
    gap = 1 - 3 * h1
    assert h1 > 0 and h2 >= 0 and h1 + h2 == 1
    assert F(-2) <= gap <= F(-4, 5)
    if t == 0:
        assert h2 == 0 and gap == -2
    else:
        assert h2 > 0


def main():
    groups = {}
    for n in range(2, 8):
        slice_recovery(n)
    groups['rational_slice_recovery_n2_to_n7'] = 6
    for n in range(2, 7):
        singular_line_point(n)
    groups['excluded_affine_chart_point_n2_to_n6'] = 5
    for n in range(2, 7):
        for k in range(n + 1):
            balanced_two_line_boundary(n, k)
    groups['balanced_two_line_configurations_n2_to_n6'] = 25
    for n in range(3, 7):
        signed_admissible_line(n)
    groups['negative_admissible_weights_n3_to_n6'] = 4
    for k in range(21):
        punctured_curve(F(k, 40))
    groups['meromorphic_two_puncture_arc_samples'] = 21
    for label, count in groups.items():
        print(f'PASS {label}: {count} cases')
    print(f'PASS all {sum(groups.values())} finite exact sanity cases')


if __name__ == '__main__':
    main()
