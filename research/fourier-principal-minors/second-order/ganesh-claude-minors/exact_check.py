#!/usr/bin/env python3
"""Exact check of README.md in characteristic zero, to finite q-adic precision.

Written by a separate referee pass (another agent of the same session, same
GitHub owner) that was given only README.md; it imports nothing from so_lib.py.
Unlike validate.py it is not limited to e + 1 <= q - 2: fibres go up to size q,
so the levels c up to q - 1 and the coefficient Gamma_(q-1) are exercised.


Works in O/q^r = (Z/q^r)[t]/(Phi_q(1+t)) for pairs (M, q) with q = 1 mod M, so
that W = Z_q and eta is a Teichmuller lift.  Checks
  L1   Lemma 1: v(m! Gamma_m - pi^m) >= m+1 for all 0 <= m <= q-1, Gamma_0 = 1,
       Gamma_(q-1) = -q/(q-1), and zeta^x = sum Gamma_m [x]^m.
  P2   Proposition 2 by brute force over all pairs of fillings (tiny cases).
  T3   Theorem 3(i): v(D_K - U (A0 + pi Lam)) >= e+2, v(U) = e, with
       U = prod V(x_a)^2 prod Gamma_m^|T_m|, signs exactly as in the README.
  C4   Corollary 4 row sums and the displayed formula for Lam (q odd).
"""
import itertools
import random
import sys
from math import comb, factorial


class Ring:
    def __init__(self, q, r):
        self.q, self.r, self.Q, self.L = q, r, q ** r, q - 1
        self.red = [comb(q, i + 1) for i in range(q - 1)]  # t^(q-1) = -sum red[i] t^i
        self.zero = tuple([0] * self.L)
        self.one = tuple([1] + [0] * (self.L - 1))
        self.pi = self.mul_raw([0, 1]) if self.L > 1 else (-(2) % self.Q,)

    def const(self, n):
        return tuple([n % self.Q] + [0] * (self.L - 1))

    def mul_raw(self, coeffs):
        """reduce an arbitrary-length coefficient list"""
        c = list(coeffs) + [0] * max(0, self.L - len(coeffs))
        L, Q = self.L, self.Q
        for k in range(len(c) - 1, L - 1, -1):
            x = c[k] % Q
            c[k] = 0
            if x:
                for i in range(L):
                    c[k - L + i] -= x * self.red[i]
        return tuple(x % Q for x in c[:L])

    def add(self, u, v):
        return tuple((x + y) % self.Q for x, y in zip(u, v))

    def sub(self, u, v):
        return tuple((x - y) % self.Q for x, y in zip(u, v))

    def mul(self, u, v):
        L = self.L
        c = [0] * (2 * L - 1)
        for i, x in enumerate(u):
            if x:
                for j, y in enumerate(v):
                    c[i + j] += x * y
        return self.mul_raw(c)

    def smul(self, n, u):
        return tuple((n * x) % self.Q for x in u)

    def pw(self, u, n):
        res = self.one
        for _ in range(n):
            res = self.mul(res, u)
        return res

    def val(self, u):
        """pi-adic valuation, capped at r(q-1) (returns that cap for 0)."""
        q, best = self.q, self.r * (self.q - 1)
        for i, x in enumerate(u):
            x %= self.Q
            if x == 0:
                continue
            k = 0
            while x % q == 0:
                x //= q
                k += 1
            best = min(best, i + (q - 1) * k)
        return best


def det_generic(rows, add, sub, mul, zero, one):
    n = len(rows)
    memo = {}

    def go(i, cols):
        if i == n:
            return one
        if cols in memo:
            return memo[cols]
        acc = zero
        for pos, j in enumerate(cols):
            term = mul(rows[i][j], go(i + 1, cols[:pos] + cols[pos + 1:]))
            acc = sub(acc, term) if pos % 2 else add(acc, term)
        memo[cols] = acc
        return acc

    return go(0, tuple(range(n)))


def idet(rows, Q):
    if not rows:
        return 1
    return det_generic(rows, lambda a, b: (a + b) % Q, lambda a, b: (a - b) % Q,
                       lambda a, b: (a * b) % Q, 0, 1)


def teich(x, q, r):
    Q = q ** r
    v = x % Q
    for _ in range(r):
        v = pow(v, q, Q)
    return v


def gammas(R):
    """Solve zeta^x = sum_m Gamma_m [x]^m by Lagrange interpolation over Z/q^r."""
    q, Q = R.q, R.Q
    tau = [teich(x, q, R.r) for x in range(q)]
    zeta = R.add(R.one, R.pi)
    G = [R.zero] * q
    for x in range(q):
        # L_x(T) = prod_{y != x} (T - tau_y) / (tau_x - tau_y)
        poly = [1]
        den = 1
        for y in range(q):
            if y == x:
                continue
            new = [0] * (len(poly) + 1)
            for i, c in enumerate(poly):
                new[i + 1] = (new[i + 1] + c) % Q
                new[i] = (new[i] - c * tau[y]) % Q
            poly = new
            den = den * (tau[x] - tau[y]) % Q
        inv = pow(den, -1, Q)
        zx = R.pw(zeta, x)
        for m in range(q):
            G[m] = R.add(G[m], R.smul(poly[m] * inv % Q, zx))
    return G, tau


def check_lemma1(q, r):
    R = Ring(q, r)
    G, tau = gammas(R)
    zeta = R.add(R.one, R.pi)
    assert G[0] == R.one
    for x in range(q):
        acc = R.zero
        for m in range(q):
            acc = R.add(acc, R.smul(1 if m == 0 else pow(tau[x], m, R.Q), G[m]))
        assert acc == R.pw(zeta, x)
    for m in range(q):
        diff = R.sub(R.smul(factorial(m), G[m]), R.pw(R.pi, m))
        assert R.val(G[m]) == m, (q, m, R.val(G[m]))
        assert R.val(diff) >= m + 1, (q, m, R.val(diff))
    assert G[q - 1] == R.const((-q) * pow(q - 1, -1, R.Q)), "Gamma_(q-1)"
    return R, G, tau


def elem_of_order(M, q, r):
    for g in range(1, q):
        if all(pow(g, d, q) != 1 for d in range(1, M)) and pow(g, M, q) == 1:
            return teich(g, q, r)
    raise ValueError


def setup(M, q, r):
    R, G, tau = check_lemma1(q, r)
    eta = elem_of_order(M, q, r)
    ep = [pow(eta, k, R.Q) for k in range(M)]
    return R, G, tau, ep


def FM_minor(ep, M, Rw, Cl, Q):
    return idet([[ep[(a * b) % M] for b in Cl] for a in Rw], Q)


def direct(R, ep, M, fib):
    q = R.q
    zeta = R.add(R.one, R.pi)
    zp = [R.pw(zeta, n) for n in range(q)]
    idx = [(a, b) for a in sorted(fib) for b in fib[a]]
    rows = [[R.smul(ep[(a * a2) % M], zp[(b * b2) % q]) for (a2, b2) in idx] for (a, b) in idx]
    return det_generic(rows, R.add, R.sub, R.mul, R.zero, R.one)


def theorem3(R, G, tau, ep, M, fib):
    q, Q = R.q, R.Q
    T = [sorted(a for a in fib if len(fib[a]) > m) for m in range(q + 1)]
    e = sum(comb(len(v), 2) for v in fib.values())
    d = [FM_minor(ep, M, T[m], T[m], Q) for m in range(q)]
    A0 = 1
    for m in range(q):
        A0 = A0 * d[m] % Q
    sig = {a: sum(tau[b] for b in v) % Q for a, v in fib.items()}
    Lam = 0
    Lam_diff = 0  # the "differences" part of Corollary 4
    sq = 0
    for c in range(1, q):
        shell = [a for a in T[c - 1] if a not in T[c]]
        rest = 1
        for m in range(q):
            if m not in (c - 1, c):
                rest = rest * d[m] % Q
        cinv = pow(c, -1, Q)
        Qm = {}
        for a in shell:
            for a2 in shell:
                Tc1 = T[c - 1]
                cof = FM_minor(ep, M, [x for x in Tc1 if x != a], [x for x in Tc1 if x != a2], Q)
                if (Tc1.index(a) + Tc1.index(a2)) % 2:
                    cof = -cof
                Rw, Cl = sorted(T[c] + [a]), sorted(T[c] + [a2])
                E = FM_minor(ep, M, Rw, Cl, Q)
                if (Rw.index(a) + Cl.index(a2)) % 2:
                    E = -E
                Qm[a, a2] = cof * E % Q
        for a in shell:  # Corollary 4, row sums
            assert sum(Qm[a, a2] for a2 in shell) % Q == d[c - 1] * d[c] % Q, "row sum"
            assert all(Qm[a, a2] == Qm[a2, a] for a2 in shell), "symmetry"
            sq = (sq + sig[a] * sig[a] * cinv) % Q
        s = sum(sig[a] * sig[a2] * Qm[a, a2] for a in shell for a2 in shell)
        Lam = (Lam + cinv * rest * s) % Q
        s2 = sum((sig[a] - sig[a2]) ** 2 * Qm[a, a2] for a in shell for a2 in shell)
        Lam_diff = (Lam_diff + cinv * rest * s2) % Q
    if q != 2:
        assert Lam == (A0 * sq - pow(2, -1, Q) * Lam_diff) % Q, "Corollary 4 formula"
    U = R.one
    for v in fib.values():
        x = [tau[b] for b in v]
        V = 1
        for i in range(len(x)):
            for j in range(i + 1, len(x)):
                V = V * (x[j] - x[i]) % Q
        U = R.smul(V * V % Q, U)
    for m in range(q):
        U = R.mul(U, R.pw(G[m], len(T[m])))
    return e, U, A0, Lam, T


def check_T3(R, G, tau, ep, M, fib):
    e, U, A0, Lam, T = theorem3(R, G, tau, ep, M, fib)
    cap = R.r * (R.q - 1)
    assert e + 2 < cap, "precision too low"
    assert R.val(U) == e
    D = direct(R, ep, M, fib)
    pred = R.mul(U, R.add(R.const(A0), R.smul(Lam, R.pi)))
    v = R.val(R.sub(D, pred))
    assert v >= e + 2, ("Theorem 3(i) fails", M, R.q, fib, e, v)
    return e, A0 % R.q == 0, Lam % R.q == 0, R.val(D)


def check_P2(R, G, tau, ep, M, fib):
    """Brute force over pairs of fillings, signs epsilon as defined in the README."""
    q, Q = R.q, R.Q
    A = sorted(fib)
    D = direct(R, ep, M, fib)
    fills = list(itertools.product(*[list(itertools.combinations(range(q), len(fib[a]))) for a in A]))

    def info(I):
        elems = [(a, m) for a, Ia in zip(A, I) for m in Ia]
        inv = sum(1 for (a, m) in elems for (a2, m2) in elems if a < a2 and m > m2)
        lev = [sorted(a for a, Ia in zip(A, I) if m in Ia) for m in range(q)]
        vd = 1
        for a, Ia in zip(A, I):
            x = [tau[b] for b in fib[a]]
            vd = vd * idet([[1 if m == 0 else pow(xi, m, Q) for m in Ia] for xi in x], Q) % Q
        return (-1) ** inv, lev, vd

    infos = [info(I) for I in fills]
    total = R.zero
    for (s1, l1, v1) in infos:
        for (s2, l2, v2) in infos:
            if any(len(x) != len(y) for x, y in zip(l1, l2)):
                continue
            coef = s1 * s2 * v1 * v2 % Q
            g = R.one
            for m in range(q):
                coef = coef * FM_minor(ep, M, l1[m], l2[m], Q) % Q
                g = R.mul(g, R.pw(G[m], len(l1[m])))
            total = R.add(total, R.smul(coef, g))
    assert total == D, ("Proposition 2 fails", M, q, fib)


def random_fib(M, q, rng, maxpts, maxc):
    while True:
        supp = rng.sample(range(M), rng.randint(1, min(M, 5)))
        fib = {a: sorted(rng.sample(range(q), rng.randint(1, min(q, maxc)))) for a in supp}
        if sum(len(v) for v in fib.values()) <= maxpts:
            return fib


def main():
    rng = random.Random(7)
    for q in (2, 3, 5, 7, 11, 13):
        check_lemma1(q, 3)
    print("L1 ok for q in 2,3,5,7,11,13 (all m up to q-1, incl. Gamma_(q-1) and q=2)")

    # Proposition 2, brute force
    n = 0
    for (M, q) in ((1, 2), (2, 3), (2, 5), (4, 5), (3, 7)):
        R, G, tau, ep = setup(M, q, 3)
        for _ in range(12):
            fib = random_fib(M, q, rng, 5, 3)
            check_P2(R, G, tau, ep, M, fib)
            n += 1
    print("P2 ok on %d configurations" % n)

    # Theorem 3(i) and Corollary 4
    tot = obs = 0
    hist = {}
    for (M, q, r, maxpts, maxc, cnt) in (
        (1, 2, 6, 2, 2, 4), (2, 3, 8, 6, 3, 40), (2, 5, 6, 8, 5, 40), (4, 5, 6, 9, 5, 80),
        (3, 7, 5, 9, 7, 40), (6, 7, 5, 9, 7, 80), (5, 11, 4, 9, 6, 60), (10, 11, 4, 9, 6, 120),
        (12, 13, 3, 9, 5, 60),
    ):
        R, G, tau, ep = setup(M, q, r)
        k = 0
        while k < cnt:
            fib = random_fib(M, q, rng, maxpts, maxc)
            e = sum(comb(len(v), 2) for v in fib.values())
            if e + 2 >= r * (q - 1):
                continue
            e, a0z, lz, vD = check_T3(R, G, tau, ep, M, fib)
            k += 1
            tot += 1
            obs += a0z
            mc = max(len(v) for v in fib.values())
            hist[mc] = hist.get(mc, 0) + 1
            if q != 2:
                if not a0z:
                    assert vD == e
                elif not lz:
                    assert vD == e + 1
        print("T3/C4 ok (M,q)=(%d,%d): %d configs" % (M, q, cnt), flush=True)
    print("T3 total %d, with A0 = 0 mod q: %d; max fibre size histogram %s" % (tot, obs, sorted(hist.items())))


if __name__ == "__main__":
    main()
