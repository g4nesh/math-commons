#!/usr/bin/env python3
"""Second-order lifting: shared code.

Notation (see README.md).  M is a positive integer, q is an odd prime not dividing M.
A chart is an irreducible factor h of Phi_M mod q; the residue field is
F = F_q[x]/(h) and eta = x mod h is a primitive M-th root of unity in F.
An index set K in Z_M x Z_q is given by its fibres: fib[a] = sorted list of the
b in {0..q-1} with (a, b) in K.  c_a = len(fib[a]).  Level sets
T_m = {a : c_a > m} for m = 0, 1, ..., and e = sum_a C(c_a, 2).

second_order(...) returns (e, pre, A0, Lam) with the following meaning: in the
ring Z[zeta_(qM)] / (chart) = F[t]/(t^(q-1)), where zeta_q = 1 + t, the
principal minor on K equals
        pre * t^e * (A0 + O(t))                 in general, and
        pre * t^(e+1) * (Lam + O(t))            when A0 = 0.
A0 is the first-order product of principal minors of F_M over the level sets;
Lam is the second-order term of README.md, Theorem 3.

Standard library only.  The cyclotomic factorisation is imported from the
frontier folder of the same contributor.
"""
import os
import sys
from math import comb, factorial

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "frontier", "ganesh-claude-minors"))
import frontier as fr  # noqa: E402  (cyclotomic_factors only)


class Field:
    """F_p[x]/(h) through discrete logarithms: the element g^k is the int k
    (0 <= k < q-1) and zero is -1."""

    def __init__(self, p, h):
        self.p, self.f, self.h = p, len(h) - 1, h
        f = self.f
        q = p ** f
        self.q1 = q - 1

        def smul(a, b):
            t = [0] * (2 * f - 1)
            for i, x in enumerate(a):
                if x:
                    for j, y in enumerate(b):
                        t[i + j] += x * y
            for k in range(2 * f - 2, f - 1, -1):
                c = t[k] % p
                t[k] = 0
                if c:
                    for j in range(f):
                        t[k - f + j] -= c * h[j]
            return tuple(x % p for x in t[:f])

        one = tuple([1] + [0] * (f - 1))
        primes, r, d = [], q - 1, 2
        while d * d <= r:
            if r % d == 0:
                primes.append(d)
                while r % d == 0:
                    r //= d
            d += 1
        if r > 1:
            primes.append(r)

        def spow(a, e):
            res, base = one, a
            while e:
                if e & 1:
                    res = smul(res, base)
                base = smul(base, base)
                e >>= 1
            return res

        gen = None
        for cand in range(1, q):
            g, v = [], cand
            for _ in range(f):
                g.append(v % p)
                v //= p
            g = tuple(g)
            if all(spow(g, (q - 1) // r) != one for r in primes) or q == 2:
                gen = g
                break
        exp, log = [], {}
        x = one
        for k in range(q - 1):
            exp.append(x)
            log[x] = k
            x = smul(x, gen)
        assert len(log) == q - 1
        self.zech = []
        for k in range(q - 1):
            e = list(exp[k])
            e[0] = (e[0] + 1) % p
            self.zech.append(log.get(tuple(e), -1))
        self.negl = 0 if p == 2 else (q - 1) // 2
        self.zero, self.one = -1, 0
        self.ints = [-1] + [log[tuple([n] + [0] * (f - 1))] for n in range(1, p)]
        self.eta = log[((-h[0]) % p,)] if f == 1 else log[tuple([0, 1] + [0] * (f - 2))]

    def from_int(self, n):
        return self.ints[n % self.p]

    def mul(self, a, b):
        return -1 if a < 0 or b < 0 else (a + b) % self.q1

    def inv(self, a):
        return (-a) % self.q1

    def neg(self, a):
        return a if a < 0 else (a + self.negl) % self.q1

    def add(self, a, b):
        if a < 0:
            return b
        if b < 0:
            return a
        z = self.zech[(b - a) % self.q1]
        return -1 if z < 0 else (a + z) % self.q1

    def sub(self, a, b):
        return self.add(a, self.neg(b))

    def det(self, rows):
        n = len(rows)
        A = [r[:] for r in rows]
        d = 0
        for c in range(n):
            piv = next((r for r in range(c, n) if A[r][c] >= 0), None)
            if piv is None:
                return -1
            if piv != c:
                A[c], A[piv] = A[piv], A[c]
                d = self.neg(d)
            d = self.mul(d, A[c][c])
            inv = self.inv(A[c][c])
            for r in range(c + 1, n):
                if A[r][c] >= 0:
                    fac = self.mul(A[r][c], inv)
                    Ar, Ac = A[r], A[c]
                    for j in range(c + 1, n):
                        if Ac[j] >= 0:
                            Ar[j] = self.sub(Ar[j], self.mul(fac, Ac[j]))
        return d


class Chart:
    """A chart of (M, q): the field, the powers of eta, and cached minors of F_M."""

    def __init__(self, M, q, h):
        self.M, self.q, self.h = M, q, h
        self.K = Field(q, h)
        K = self.K
        self.ep = [0]
        for _ in range(M):
            self.ep.append(K.mul(self.ep[-1], K.eta))
        assert self.ep[M] == 0 and all(self.ep[d] != 0 for d in range(1, M))
        self._cache = {}

    def minor(self, R, C):
        """det F_M[R, C] in the chart; R and C are sorted lists."""
        key = (tuple(R), tuple(C))
        if key not in self._cache:
            M, ep = self.M, self.ep
            self._cache[key] = self.K.det([[ep[(a * b) % M] for b in C] for a in R])
        return self._cache[key]

    def cofactor(self, T, a, a2):
        """Signed cofactor of F_M[T, T] at (a, a2); T sorted."""
        v = self.minor([x for x in T if x != a], [x for x in T if x != a2])
        return self.K.neg(v) if (T.index(a) + T.index(a2)) % 2 else v

    def bordered(self, T, a, a2):
        """(-1)^(pos(a)+pos(a2)) det F_M[T u a, T u a2] for a, a2 outside the sorted list T."""
        R, C = sorted(T + [a]), sorted(T + [a2])
        v = self.minor(R, C)
        return self.K.neg(v) if (R.index(a) + C.index(a2)) % 2 else v


def charts(M, q):
    f, facs = fr.cyclotomic_factors(M, q)
    return [Chart(M, q, h) for h in facs]


def level_sets(fib, q):
    return [sorted(a for a in fib if len(fib[a]) > m) for m in range(q + 1)]


def second_order(ch, fib):
    """(e, pre, A0, Lam) for the chart ch and the index set with fibres fib."""
    K, q = ch.K, ch.q
    T = level_sets(fib, q)
    e = sum(comb(len(v), 2) for v in fib.values())
    pre = 1
    for v in fib.values():
        for i in range(len(v)):
            for j in range(i + 1, len(v)):
                pre = pre * (v[j] - v[i]) ** 2 % q
    for m in range(q):
        pre = pre * pow(pow(factorial(m), q - 2, q), len(T[m]), q) % q
    d = [ch.minor(T[m], T[m]) for m in range(q + 1)]
    A0 = 0
    for m in range(q):
        A0 = K.mul(A0, d[m])
    sig = {a: K.from_int(sum(v)) for a, v in fib.items()}
    Lam = -1
    for c in range(1, q):
        shell = [a for a in T[c - 1] if a not in T[c]]
        if not shell:
            continue
        rest = 0
        for m in range(q):
            if m not in (c - 1, c):
                rest = K.mul(rest, d[m])
        if rest < 0:
            continue
        Q = -1
        for a in shell:
            if sig[a] < 0:
                continue
            for a2 in shell:
                if sig[a2] < 0:
                    continue
                term = K.mul(K.mul(sig[a], sig[a2]), K.mul(ch.cofactor(T[c - 1], a, a2), ch.bordered(T[c], a, a2)))
                Q = K.add(Q, term)
        Lam = K.add(Lam, K.mul(K.mul(K.from_int(pow(c, q - 2, q)), rest), Q))
    return e, K.from_int(pre), A0, Lam


def direct_series(ch, fib):
    """The principal minor on K in F[t]/(t^(q-1)) with zeta_q = 1 + t, as the list of
    its q-1 coefficients.  Expansion along rows with a memo over column sets, so
    it is meant for |K| up to about 11."""
    K, q, M, ep = ch.K, ch.q, ch.M, ch.ep
    L = q - 1
    idx = [(a, b) for a in sorted(fib) for b in fib[a]]
    n = len(idx)
    pw = {}

    def one_plus_t(nexp):
        nexp %= q
        if nexp not in pw:
            pw[nexp] = [K.from_int(comb(nexp, k)) for k in range(L)]
        return pw[nexp]

    rows = [[[K.mul(ep[(a * a2) % M], x) for x in one_plus_t(b * b2)] for (a2, b2) in idx] for (a, b) in idx]

    def smul(u, v):
        r = [-1] * L
        for i, x in enumerate(u):
            if x < 0:
                continue
            for j in range(L - i):
                if v[j] >= 0:
                    r[i + j] = K.add(r[i + j], K.mul(x, v[j]))
        return r

    memo = {}

    def go(i, cols):
        if i == n:
            return [0] + [-1] * (L - 1)
        key = cols
        if key in memo:
            return memo[key]
        acc = [-1] * L
        for pos, j in enumerate(cols):
            term = smul(rows[i][j], go(i + 1, cols[:pos] + cols[pos + 1:]))
            if pos % 2:
                acc = [K.sub(x, y) for x, y in zip(acc, term)]
            else:
                acc = [K.add(x, y) for x, y in zip(acc, term)]
        memo[key] = acc
        return acc

    return go(0, tuple(range(n)))
