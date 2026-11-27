#!/usr/bin/env python3
"""Independent check of pmscan on small bases.

For each pair (M, p) on the command line this script
  1. recomputes the cyclotomic polynomial Phi_M from the Moebius product formula,
  2. checks that the factors used by frontier.py are monic of degree ord_M(p) and
     multiply to Phi_M mod p,
  3. computes every principal minor of F_M over F_p[x]/(h) by plain Gaussian
     elimination, one subset at a time, with its own field arithmetic,
  4. compares the set of vanishing minors with the list printed by pmscan.

It shares no arithmetic code with pmscan.c and uses a different algorithm
(one elimination per subset instead of the Schur-complement recursion).

Usage: python3 verify_small.py BIN_DIR M:p [M:p ...]
Standard library only.
"""
import sys
from math import gcd

sys.path.insert(0, ".")
import frontier as fr  # only for the factor lists and for running pmscan


def moebius(n):
    r, d = 1, 2
    while d * d <= n:
        if n % d == 0:
            n //= d
            if n % d == 0:
                return 0
            r = -r
        d += 1
    return -r if n > 1 else r


def mul_int(a, b):
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            r[i + j] += x * y
    return r


def div_exact_int(a, b):
    a = a[:]
    q = [0] * (len(a) - len(b) + 1)
    for s in range(len(q) - 1, -1, -1):
        c = a[s + len(b) - 1] // b[-1]
        q[s] = c
        for i, y in enumerate(b):
            a[s + i] -= c * y
    assert not any(a)
    return q


def phi_poly(M):
    num, den = [1], [1]
    for d in range(1, M + 1):
        if M % d == 0:
            mu = moebius(M // d)
            fac = [-1] + [0] * (d - 1) + [1]
            if mu == 1:
                num = mul_int(num, fac)
            elif mu == -1:
                den = mul_int(den, fac)
    return div_exact_int(num, den)


class Field:
    """F_p[x]/(h), elements are tuples of f digits; multiplication through a
    discrete-log table built here from schoolbook multiplication."""

    def __init__(self, p, h):
        self.p, self.f = p, len(h) - 1
        self.h = h
        self.zero = tuple([0] * self.f)
        self.one = tuple([1] + [0] * (self.f - 1))
        q = p ** self.f
        gen = None
        for cand in range(1, q):
            g = self._from_int(cand)
            x, k = g, 1
            while x != self.one:
                x, k = self._slow_mul(x, g), k + 1
                if k > q:
                    raise AssertionError("not a field")
            if k == q - 1:
                gen = g
                break
        assert gen is not None or q == 2
        self.exp, self.log = [], {}
        x = self.one
        for k in range(q - 1):
            self.exp.append(x)
            assert x not in self.log
            self.log[x] = k
            x = self._slow_mul(x, gen) if gen else x
        assert len(self.log) == q - 1
        self.q1 = q - 1

    def _from_int(self, v):
        d = []
        for _ in range(self.f):
            d.append(v % self.p)
            v //= self.p
        return tuple(d)

    def _slow_mul(self, a, b):
        p, f, h = self.p, self.f, self.h
        t = [0] * (2 * f - 1)
        for i, x in enumerate(a):
            for j, y in enumerate(b):
                t[i + j] += x * y
        for k in range(2 * f - 2, f - 1, -1):
            c = t[k] % p
            t[k] = 0
            for j in range(f):
                t[k - f + j] -= c * h[j]
        return tuple(x % p for x in t[:f])

    def mul(self, a, b):
        if a == self.zero or b == self.zero:
            return self.zero
        return self.exp[(self.log[a] + self.log[b]) % self.q1]

    def inv(self, a):
        return self.exp[(-self.log[a]) % self.q1]

    def sub(self, a, b):
        p = self.p
        return tuple((x - y) % p for x, y in zip(a, b))

    def zeta(self):
        if self.f == 1:
            return ((-self.h[0]) % self.p,)
        return tuple([0, 1] + [0] * (self.f - 2))


def is_singular(K, rows):
    n = len(rows)
    A = [r[:] for r in rows]
    for c in range(n):
        piv = next((r for r in range(c, n) if A[r][c] != K.zero), None)
        if piv is None:
            return True
        A[c], A[piv] = A[piv], A[c]
        inv = K.inv(A[c][c])
        for r in range(c + 1, n):
            if A[r][c] != K.zero:
                fac = K.mul(A[r][c], inv)
                A[r] = [K.sub(x, K.mul(fac, y)) for x, y in zip(A[r], A[c])]
    return False


def check(M, p, bin_dir):
    f, facs = fr.cyclotomic_factors(M, p)
    target = [c % p for c in phi_poly(M)]
    prod = [1]
    for h in facs:
        assert len(h) == f + 1 and h[-1] == 1
        prod = [c % p for c in mul_int(prod, h)]
    assert prod == target, "factors do not multiply to Phi_M mod p"
    assert fr.mult_order(p, M) == f and gcd(M, p) == 1
    total = 0
    for h in facs:
        K = Field(p, h)
        z = K.zeta()
        zp = [K.one]
        for _ in range(M):
            zp.append(K.mul(zp[-1], z))
        assert zp[M] == K.one and all(zp[d] != K.one for d in range(1, M))
        mine = set()
        for mask in range(1, 1 << M):
            idx = [i for i in range(M) if mask >> i & 1]
            if is_singular(K, [[zp[(i * j) % M] for j in idx] for i in idx]):
                mine.add(mask)
        theirs = set(fr.run_chart(M, p, f, h, bin_dir, ["list"])["zero_masks"])
        assert mine == theirs, "zero sets differ for h=%s" % h
        total += len(mine)
        print("  (M,p)=(%d,%d) h=%s: %d vanishing minors, identical to pmscan" % (M, p, h, len(mine)))
    return total


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(2)
    for arg in sys.argv[2:]:
        M, p = map(int, arg.split(":"))
        check(M, p, sys.argv[1])
    print("all pairs agree")
