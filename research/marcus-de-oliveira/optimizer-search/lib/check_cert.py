#!/usr/bin/env python3
"""Independent checker for an MOC outside-certificate (standard library only).

Shares no code with the builder. Usage: python3 check_cert.py cert.json
Exit 0 iff: U U* = I exactly; with d = det(D_a + U D_b U*) and the real
functional L(w) = p Re w + q Im w, L(d) > L(z_sigma) for every one of the n!
permutations sigma (enumerated here), and the claimed gap equals the recomputed one.
Method: everything is scaled to Gaussian integers; d by fraction-free Bareiss
elimination over Z[i]; z_sigma by a recursive permutation generator.
"""
import json
import sys
from fractions import Fraction
from math import factorial, lcm


def parse(c):
    return Fraction(c[0]), Fraction(c[1])


def cmul(u, v):
    return (u[0] * v[0] - u[1] * v[1], u[0] * v[1] + u[1] * v[0])


def cadd(u, v):
    return (u[0] + v[0], u[1] + v[1])


def cexact_div(u, v):
    n = v[0] * v[0] + v[1] * v[1]
    x = u[0] * v[0] + u[1] * v[1]
    y = u[1] * v[0] - u[0] * v[1]
    if x % n or y % n:
        raise ArithmeticError("Bareiss division not exact")
    return (x // n, y // n)


def bareiss_det(M):
    n = len(M)
    A = [list(r) for r in M]
    sign, prev = 1, (1, 0)
    for k in range(n - 1):
        if A[k][k] == (0, 0):
            sw = next((r for r in range(k + 1, n) if A[r][k] != (0, 0)), None)
            if sw is None:
                return (0, 0)
            A[k], A[sw] = A[sw], A[k]
            sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                t = cadd(cmul(A[i][j], A[k][k]), cmul((-A[i][k][0], -A[i][k][1]), A[k][j]))
                A[i][j] = cexact_div(t, prev)
        prev = A[k][k]
    d = A[n - 1][n - 1]
    return (sign * d[0], sign * d[1])


def all_perms(n):
    if n == 0:
        yield ()
        return
    for p in all_perms(n - 1):
        for pos in range(n):
            yield p[:pos] + (n - 1,) + p[pos:]


def check(cert):
    n = cert["n"]
    a = [parse(v) for v in cert["a"]]
    b = [parse(v) for v in cert["b"]]
    U = [[parse(v) for v in row] for row in cert["U"]]
    p, q = Fraction(cert["direction"][0]), Fraction(cert["direction"][1])
    if len(a) != n or len(b) != n or len(U) != n or any(len(r) != n for r in U):
        return False, "shape mismatch"
    DU = lcm(*[x.denominator for r in U for e in r for x in e])
    DS = lcm(*[x.denominator for v in a + b for x in v])
    Ui = [[(int(e[0] * DU), int(e[1] * DU)) for e in r] for r in U]
    ai = [(int(v[0] * DS), int(v[1] * DS)) for v in a]
    bi = [(int(v[0] * DS), int(v[1] * DS)) for v in b]
    # exact unitarity: Ui Ui* = DU^2 I
    for i in range(n):
        for j in range(n):
            s = (0, 0)
            for k in range(n):
                s = cadd(s, cmul(Ui[i][k], (Ui[j][k][0], -Ui[j][k][1])))
            if s != ((DU * DU if i == j else 0), 0):
                return False, "U is not unitary (entry %d,%d)" % (i, j)
    # Mt = DU^2 D_a~ + U~ D_b~ U~*  so that  D_a + U D_b U* = Mt / (DU^2 DS)
    Mt = []
    for i in range(n):
        row = []
        for j in range(n):
            s = (DU * DU * ai[i][0], DU * DU * ai[i][1]) if i == j else (0, 0)
            for k in range(n):
                s = cadd(s, cmul(cmul(Ui[i][k], bi[k]), (Ui[j][k][0], -Ui[j][k][1])))
            row.append(s)
        Mt.append(row)
    dt = bareiss_det(Mt)                  # d = dt / (DU^2 DS)^n
    scale_z = DU ** (2 * n)               # z = zt / DS^n  ->  zt * DU^(2n) / (DU^2 DS)^n
    Ld = p * dt[0] + q * dt[1]
    best, count = None, 0
    for s in all_perms(n):
        count += 1
        z = (1, 0)
        for i in range(n):
            z = cmul(z, cadd(ai[i], bi[s[i]]))
        Lz = (p * z[0] + q * z[1]) * scale_z
        if best is None or Lz > best:
            best = Lz
    if count != factorial(n):
        return False, "permutation enumeration incomplete"
    common = Fraction(1, (DU * DU * DS) ** n)
    gap = (Ld - best) * common
    if gap <= 0:
        return False, "no strict separation: gap = %s" % gap
    if "gap" in cert and Fraction(cert["gap"]) != gap:
        return False, "claimed gap %s != recomputed %s" % (cert["gap"], gap)
    return True, "certified outside: gap = %s over %d permutations" % (gap, count)


if __name__ == "__main__":
    with open(sys.argv[1]) as f:
        ok, msg = check(json.load(f))
    print(("PASS: " if ok else "FAIL: ") + msg)
    sys.exit(0 if ok else 1)
