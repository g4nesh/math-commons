"""Exact checks for the reduction of (E20) and (E21) to one 3x3 block, and their identification with (star)
on blocks of 4x4 unitaries (README, Propositions B and S).

Every linear identity below is between affine functions of squared minors. It is checked on all
permutation matrices (vertices of Q5, resp. Q4) in exact integer arithmetic; since D(U) lies in the affine
hull of the vertices for every unitary U (R1 of ../facets/README.md for n = 5, R1 of
../../pmsh/equivalence.md for n = 4), the identity then holds on all of U(5), resp. U(4). The square-root
terms are compared as minors, symbol by symbol. A floating-point check on random unitaries follows.

usage: reduce_check.py
"""
import itertools

import numpy as np


def sq(P, K, J):
    """Squared minor |det P[K,J]|^2 (empty minor = 1) for a permutation matrix given as a tuple."""
    return int(tuple(sorted(P[k] for k in K)) == tuple(sorted(J)))


def x5(P, s):
    """Q5 coordinate written as in ../facets: 'x01>12' means |det U[{0,1},{1,2}]|^2."""
    K, J = s[1:].split(">")
    return sq(P, [int(c) for c in K], [int(c) for c in J])


def block_terms(P, R, C, n):
    """Squared minors of W = U[R, C] in the labels of the README: rows (a, p, q), cols (a', g, d).
    |det W|^2 is written as the complementary minor (Jacobi), so every term is a Q_n coordinate."""
    a, p, q = R
    a_, g, d = C
    Rc = [i for i in range(n) if i not in R]
    Cc = [j for j in range(n) if j not in C]
    return {
        "w_aa'": sq(P, [a], [a_]), "w_pg": sq(P, [p], [g]), "w_qg": sq(P, [q], [g]),
        "w_pa'": sq(P, [p], [a_]), "w_pd": sq(P, [p], [d]),
        "M(ap;a'g)": sq(P, [a, p], [a_, g]), "M(ap;gd)": sq(P, [a, p], [g, d]), "M(ap;a'd)": sq(P, [a, p], [a_, d]),
        "M(pq;a'g)": sq(P, [p, q], [a_, g]), "M(pq;gd)": sq(P, [p, q], [g, d]), "M(pq;a'd)": sq(P, [p, q], [a_, d]),
        "det": sq(P, Rc, Cc) if Rc else sq(P, R, C),
    }


def E_W(b):  # the E-facet as a function of the block
    return b["w_aa'"] + b["M(ap;a'g)"] + b["M(ap;gd)"] + b["M(pq;a'g)"] - b["det"]


def Y20_W(b):  # defect of column g in rows p, q
    return b["w_pg"] + b["w_qg"] - b["M(pq;a'g)"] - b["M(pq;gd)"]


def Y21_W(b):  # defect of row p in columns a', d
    return b["w_pa'"] + b["w_pd"] - b["M(ap;a'd)"] - b["M(pq;a'd)"]


def check_n5():
    R20, C20 = (1, 0, 4), (3, 1, 2)
    R21, C21 = (1, 0, 3), (2, 3, 1)
    bad = 0
    for P in itertools.permutations(range(5)):
        g20 = x5(P, "x1>3") + x5(P, "x01>12") + x5(P, "x01>13") + x5(P, "x04>13") - x5(P, "x23>04")
        g21 = x5(P, "x1>2") + x5(P, "x01>13") + x5(P, "x01>23") + x5(P, "x03>23") - x5(P, "x24>04")
        b20, b21 = block_terms(P, R20, C20, 5), block_terms(P, R21, C21, 5)
        bad += g20 != E_W(b20)
        bad += (x5(P, "x04>01") + x5(P, "x04>14")) != Y20_W(b20)
        bad += x5(P, "x01>12") != b20["M(ap;gd)"]
        bad += g21 != E_W(b21)
        bad += (x5(P, "x02>12") + x5(P, "x04>12")) != Y21_W(b21)
        bad += x5(P, "x01>23") != b21["M(ap;a'g)"]
    print(f"n = 5, 120 vertices: E-facets, Y-terms and square-root partners agree with the block forms; failures {bad}")
    return bad


def star_terms(P, r, c):
    """(star) with rows 1..4 -> r and columns 1..4 -> c: (linear part, partner a = D[12,12], partner b = D[13,14]).
    (star): D[13,13] + D[14,14] + D[13,23] + D[12,12] + D[13,14] + 2 sqrt(D[12,12] D[13,14]) - D[12,34] >= 0."""
    D = lambda i, j, k, l: sq(P, [r[i - 1], r[j - 1]], [c[k - 1], c[l - 1]])
    lin = D(1, 3, 1, 3) + D(1, 4, 1, 4) + D(1, 3, 2, 3) + D(1, 2, 1, 2) + D(1, 3, 1, 4) - D(1, 2, 3, 4)
    return lin, (r[0], r[1], c[0], c[1]), (r[0], r[2], c[0], c[3])


def check_n4():
    a, p, q, o = 0, 1, 2, 3
    a_, g, d, o_ = 0, 1, 2, 3
    bad = 0
    for P in itertools.permutations(range(4)):
        b = block_terms(P, (a, p, q), (a_, g, d), 4)
        # (E20) on U(4): linear part E + Y, partners M(ap;gd) and Y = D[pq;g o'].
        bad += Y20_W(b) != sq(P, [p, q], [g, o_])
        lin, pa, pb = star_terms(P, (p, q, a, o), (g, o_, a_, d))
        bad += (E_W(b) + Y20_W(b)) != lin
        # (E21) on U(4): linear part E + A, partners A = M(ap;a'g) and Y = D[po;a'd].
        bad += Y21_W(b) != sq(P, [p, o], [a_, d])
        lin, pa2, pb2 = star_terms(P, (q, o, a, p), (o_, d, a_, g))
        bad += (E_W(b) + b["M(ap;a'g)"]) != lin
        # Jacobi on U(4) identifies the (star) partners with the (E21) partners.
        bad += sq(P, [q, o], [o_, d]) != b["M(ap;a'g)"]
        bad += sq(P, [q, a], [o_, g]) != sq(P, [p, o], [a_, d])
    # The (E20) partners are literally the (star) partners.
    assert (pa, pb) == ((p, q, g, o_), (p, a, g, d))
    assert (pa2, pb2) == ((q, o, o_, d), (q, a, o_, g))
    print(f"n = 4, 24 vertices: (E20) and (E21) on 3x3 blocks of U(4) are (star) images; failures {bad}")
    return bad


def M(W, r, c):
    return abs(np.linalg.det(W[np.ix_(r, c)])) ** 2


def E20_block(W):
    Y = abs(W[1, 1]) ** 2 + abs(W[2, 1]) ** 2 - M(W, [1, 2], [0, 1]) - M(W, [1, 2], [1, 2])
    E = abs(W[0, 0]) ** 2 + M(W, [0, 1], [0, 1]) + M(W, [0, 1], [1, 2]) + M(W, [1, 2], [0, 1]) - abs(np.linalg.det(W)) ** 2
    return E + Y + 2 * np.sqrt(M(W, [0, 1], [1, 2]) * max(Y, 0.0))


def E21_block(W):
    Y = abs(W[1, 0]) ** 2 + abs(W[1, 2]) ** 2 - M(W, [0, 1], [0, 2]) - M(W, [1, 2], [0, 2])
    E = abs(W[0, 0]) ** 2 + M(W, [0, 1], [0, 1]) + M(W, [0, 1], [1, 2]) + M(W, [1, 2], [0, 1]) - abs(np.linalg.det(W)) ** 2
    A = M(W, [0, 1], [0, 1])
    return E + A + 2 * np.sqrt(A * max(Y, 0.0))


def haar(n, rng):
    Z = rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))
    Q, R = np.linalg.qr(Z)
    return Q * (np.diag(R) / abs(np.diag(R)))


def check_float():
    rng = np.random.default_rng(0)
    m = lambda U, K, J: abs(np.linalg.det(U[np.ix_(K, J)])) ** 2
    err5 = err4 = 0.0
    for _ in range(2000):
        U = haar(5, rng)
        Y = m(U, [0, 4], [0, 1]) + m(U, [0, 4], [1, 4])
        e20 = m(U, [1], [3]) + m(U, [0, 1], [1, 2]) + m(U, [0, 1], [1, 3]) + m(U, [0, 4], [1, 3]) - m(U, [2, 3], [0, 4]) \
            + Y + 2 * np.sqrt(m(U, [0, 1], [1, 2]) * Y)
        Y = m(U, [0, 2], [1, 2]) + m(U, [0, 4], [1, 2])
        A = m(U, [0, 1], [2, 3])
        e21 = m(U, [1], [2]) + m(U, [0, 1], [1, 3]) + A + m(U, [0, 3], [2, 3]) - m(U, [2, 4], [0, 4]) + A + 2 * np.sqrt(A * Y)
        err5 = max(err5, abs(e20 - E20_block(U[np.ix_([1, 0, 4], [3, 1, 2])])),
                   abs(e21 - E21_block(U[np.ix_([1, 0, 3], [2, 3, 1])])))
        V = haar(4, rng)
        D = lambda r, c, i, j, k, l: m(V, [r[i - 1], r[j - 1]], [c[k - 1], c[l - 1]])
        star = lambda r, c: (D(r, c, 1, 3, 1, 3) + D(r, c, 1, 4, 1, 4) + D(r, c, 1, 3, 2, 3)
                             + (np.sqrt(D(r, c, 1, 2, 1, 2)) + np.sqrt(D(r, c, 1, 3, 1, 4))) ** 2 - D(r, c, 1, 2, 3, 4))
        err4 = max(err4, abs(E20_block(V[:3, :3]) - star((1, 2, 0, 3), (1, 3, 0, 2))),
                   abs(E21_block(V[:3, :3]) - star((2, 3, 0, 1), (3, 2, 0, 1))))
    print(f"floating point, 2000 Haar U(5) and U(4): max |original - block form| {err5:.1e}, "
          f"max |block form - (star) image| {err4:.1e}")
    return err5 < 1e-12 and err4 < 1e-12


def main():
    ok = check_n5() == 0 and check_n4() == 0 and check_float()
    print("ALL CHECKS PASSED" if ok else "CHECK FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
