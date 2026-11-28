#!/usr/bin/env python3
"""The valuation statements of Theorems A, C and Proposition F, tested q-adically.

validate.py works modulo the prime, where only the residues of the window
coefficients are visible.  Here the principal minor D_K itself is computed in
(Z/q^r)[zeta_q], for base 10 above 11 (q = 1 mod M, so eta is an integer root
of unity modulo q^r and each chart is a choice of eta), by a division-free
determinant algorithm, and its exact valuation at pi = zeta - 1 is compared
with the prediction:

  A. Theorem A(iv): v(D_K) = e + j0 when the least order j0 of a nonzero window
     coefficient is at most q - 2, and v(D_K) >= e + q - 1 when there is none.
  C. Theorem C: for a block of multiplicity delta, v(D_K) >= e + w, and
     v(D_K) = e + j0 when j0 <= w + q - 2 - m_max.  Bottom blocks with
     delta = 4 have w = 10 > q - 2 = 9, so they test the range beyond the window.
  F. Proposition F: for Ferrers sets, v(D_K) = e + (q - 1) sum_m v_q(d_(T_m)).
  G. A measurement for the outlook of README.md: the valuation of D_K for sets
     that contain a product block S x B, |B| = s + 1, and have a single fibre
     outside S.  No theorem of the report covers the window-blind ones.

The window coefficients come from window.c through ho_lib; the q-adic
arithmetic below shares no code with it.

Usage: python3 padic_check.py PATH_TO_window_BINARY
Standard library only; the random seed is fixed.
"""
import itertools
import random
import sys
from math import comb

import ho_lib as ho

so = ho.so
M, Q, R = 10, 11, 12
MOD = Q ** R
ZERO = [0] * Q


def cmul(x, y):
    """Product in (Z/q^r)[zeta]/(zeta^q - 1)."""
    out = [0] * Q
    for i, a in enumerate(x):
        if a:
            for j, b in enumerate(y):
                if b:
                    k = i + j
                    out[k - Q if k >= Q else k] += a * b
    return [v % MOD for v in out]


def bird_det(A):
    """Determinant over a commutative ring without division (Bird's algorithm):
    F_1 = A, F_(k+1) = mu(F_k) A, det A = (-1)^(n-1) F_n[0][0], where mu keeps the
    strict upper triangle and replaces the diagonal entry i by minus the sum of the
    later diagonal entries."""
    n = len(A)
    F = [row[:] for row in A]
    for _ in range(n - 1):
        suffix = [ZERO] * (n + 1)
        for i in range(n - 1, -1, -1):
            suffix[i] = [(u + v) % MOD for u, v in zip(suffix[i + 1], F[i][i])]
        G = [[ZERO] * n for _ in range(n)]
        for i in range(n):
            mu_row = [None] * n
            mu_row[i] = [(-v) % MOD for v in suffix[i + 1]]
            for j in range(i + 1, n):
                mu_row[j] = F[i][j]
            for j in range(n):
                acc = [0] * Q
                for k in range(i, n):
                    m = mu_row[k]
                    if any(m):
                        acc = [u + v for u, v in zip(acc, cmul(m, A[k][j]))]
                G[i][j] = [v % MOD for v in acc]
        F = G
    d = F[0][0]
    return d if n % 2 == 1 else [(-v) % MOD for v in d]


def valuation(x):
    """Valuation at pi of the class of x in Z_q[zeta], or None if x = 0 to the precision."""
    y = [(c - x[Q - 1]) % MOD for c in x[:Q - 1]]          # modulo Phi_q
    c = [0] * (Q - 1)                                       # zeta = 1 + pi
    for k, a in enumerate(y):
        for i in range(k + 1):
            c[i] = (c[i] + a * comb(k, i)) % MOD
    best = None
    for i, a in enumerate(c):
        if a:
            v = 0
            while a % Q == 0:
                a //= Q
                v += 1
            best = i + (Q - 1) * v if best is None else min(best, i + (Q - 1) * v)
    return best


def minor_valuation(eta, fib):
    K = [(a, b) for a in sorted(fib) for b in fib[a]]
    A = []
    for a, b in K:
        row = []
        for a2, b2 in K:
            x = [0] * Q
            x[(b * b2) % Q] = pow(eta, a * a2, MOD)
            row.append(x)
        A.append(row)
    return valuation(bird_det(A))


def level_minor_valuation(eta, T):
    """v_q(d_T) for T in Z_M (an integer determinant modulo q^r)."""
    if not T:
        return 0
    A = [[[pow(eta, a * a2, MOD)] + [0] * (Q - 1) for a2 in T] for a in T]
    v = valuation(bird_det(A))
    assert v is not None and v % (Q - 1) == 0
    return v // (Q - 1)


def block_fib(rng, S, s, delta, slack=1):
    outside = [a for a in range(M) if a not in S]
    while True:
        fib = {}
        for a in S:
            fib[a] = sorted(rng.sample(range(Q), s + delta if rng.random() < 0.75 else min(s + delta + slack, Q)))
        for a in outside:
            c = s if rng.random() < 0.5 else rng.randint(0, s)
            if c and rng.random() < 0.6:
                fib[a] = sorted(rng.sample(range(Q), c))
        T = so.level_sets(fib, Q)
        if [m for m in range(Q) if T[m] == sorted(S)] == list(range(s, s + delta)):
            return fib


def main():
    binary = sys.argv[1]
    rng = random.Random(1111)
    charts = so.charts(M, Q)
    data = []
    for ch in charts:
        res = (-ch.h[0]) % Q
        eta = pow(res, Q ** (R - 1), MOD)                    # Teichmueller lift
        assert pow(eta, M, MOD) == 1 and eta % Q == res
        sing = [list(S) for size in (3, 4, 5) for S in itertools.combinations(range(M), size) if ch.minor(list(S), list(S)) < 0]
        data.append((ch, eta, sing))
    limit = (Q - 1) * R

    def first_order(ch, fib, L):
        lam = ho.window_c(binary, ch, [fib], L)[0]
        return next((j for j, x in enumerate(lam) if x >= 0), None)

    # A
    n = by = 0
    hist = {}
    for _ in range(60):
        ch, eta, sing = rng.choice(data)
        if rng.random() < 0.5:
            fib = {a: sorted(rng.sample(range(Q), rng.randint(1, 3))) for a in range(M) if rng.random() < 0.5}
        else:
            fib = block_fib(rng, rng.choice(sing), rng.randint(0, 2), 1)
        e = ho.excess(fib)
        if not fib or e + 2 * Q > limit or sum(len(v) for v in fib.values()) > 26:
            continue
        j0 = first_order(ch, fib, Q - 1)
        v = minor_valuation(eta, fib)
        if j0 is None:
            assert v is None or v >= e + Q - 1
            by += 1
        else:
            assert v == e + j0
        hist[j0] = hist.get(j0, 0) + 1
        n += 1
    print("A. Theorem A(iv), base 10 above 11: %d index sets; v(D_K) = e + j0 in every case with j0 <= 9 (least orders: %s); v(D_K) >= e + 10 in the %d cases without one" % (
        n, ", ".join("%s: %d" % (k, hist[k]) for k in sorted(k for k in hist if k is not None)), by))

    # C
    print("C. Theorem C: blocks of multiplicity delta")
    total = beyond = 0
    for pos, s, delta in [("interior", 1, 2), ("interior", 2, 2), ("interior", 1, 3), ("bottom", 0, 2), ("bottom", 0, 3), ("bottom", 0, 4)]:
        ms = [delta - i for i in range(delta)] if pos == "bottom" else [min(i + 1, delta - i) for i in range(delta)]
        w, mmax = sum(ms), max(ms)
        jmax = w + Q - 2 - mmax
        n = exact = out = 0
        orders = {}
        tries = 0
        while n < 12 and tries < 200:
            tries += 1
            ch, eta, sing = rng.choice(data)
            S = rng.choice(sing)
            fib = block_fib(rng, S, s, delta)
            e = ho.excess(fib)
            if e + jmax + 2 * Q > limit or sum(len(v) for v in fib.values()) > 30:
                continue
            j0 = first_order(ch, fib, jmax + 2)
            v = minor_valuation(eta, fib)
            assert v is None or v >= e + w
            assert j0 is None or j0 >= w
            if j0 is not None and j0 <= jmax:
                assert v == e + j0
                exact += 1
                orders[j0] = orders.get(j0, 0) + 1
                beyond += j0 > Q - 2
            else:
                out += 1
            n += 1
        total += n
        print("   %s block, s = %d, delta = %d: w = %d, covered orders %d..%d; %d sets; v(D_K) >= e + w in all; v(D_K) = e + j0 in the %d with j0 in the covered range (j0: %s); %d outside it" % (
            pos, s, delta, w, w, jmax, n, exact, ", ".join("%d: %d" % kv for kv in sorted(orders.items())), out))
    print("   total %d sets; %d of them with j0 beyond the window (j0 > 9)" % (total, beyond))

    # F
    n = singular = 0
    while n < 30:
        ch, eta, sing = rng.choice(data)
        order = rng.sample(range(Q), Q)
        S = rng.choice(sing)
        s = rng.randint(0, 2)
        size = {a: rng.randint(s + 1, s + 3) if a in S else rng.randint(0, s) for a in range(M)}
        fib = {a: sorted(order[:c]) for a, c in size.items() if c}
        e = ho.excess(fib)
        T = so.level_sets(fib, Q)
        pred = e + (Q - 1) * sum(level_minor_valuation(eta, T[m]) for m in range(Q))
        if pred + Q > limit or sum(len(v) for v in fib.values()) > 30:
            continue
        assert minor_valuation(eta, fib) == pred
        singular += pred > e
        n += 1
    print("F. Proposition F: %d Ferrers sets (%d with a singular level set): v(D_K) = e + 10 * sum of v_11(d_T) in every case" % (n, singular))
    # G: a measurement, not a test of a theorem
    dist = {}
    n = 0
    while n < 48:
        ch, eta, sing = rng.choice(data)
        S = rng.choice(sing)
        s = rng.choice([1, 1, 2])
        B = rng.sample(range(Q), s + 1)
        rest = [b for b in range(Q) if b not in B]
        sizes = [s + 1] + [rng.randint(s + 2, s + 3) for _ in S[1:]]
        rng.shuffle(sizes)
        fib = {a: sorted(B + rng.sample(rest, c - s - 1)) for a, c in zip(S, sizes)}
        fib[rng.choice([a for a in range(M) if a not in S])] = sorted(rng.sample(range(Q), s))
        T = so.level_sets(fib, Q)
        if [m for m in range(Q) if T[m] == sorted(S)] != [s] or level_minor_valuation(eta, S) != 1:
            continue
        e = ho.excess(fib)
        j0 = first_order(ch, fib, 14)
        v = minor_valuation(eta, fib)
        assert v is not None
        key = ("no nonzero window coefficient to order 13" if j0 is None else "least order %d" % j0, v - e)
        assert j0 is None or (j0 <= Q - 2 and v - e == j0)
        dist[key] = dist.get(key, 0) + 1
        n += 1
    print("G. Measurement (no theorem of the report applies to the blind cases): %d sets with S = T_s singular of multiplicity one, S x B inside K for |B| = s + 1, and one fibre outside S: %s" % (
        n, "; ".join("%s, v(D_K) - e = %d: %d sets" % (k[0], k[1], c) for k, c in sorted(dist.items()))))
    print("all checks passed")


if __name__ == "__main__":
    main()
