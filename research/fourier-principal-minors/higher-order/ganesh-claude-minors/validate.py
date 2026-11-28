#!/usr/bin/env python3
"""Exact checks of the statements of README.md.

  1. Theorem A against the second-order report: the coefficients of tau^e and
     tau^(e+1) of the window polynomial are A_0 and Lambda, also when A_0 != 0.
  2. Theorem A against determinants computed directly with zeta = 1 + t in
     FF[t]/(t^(q-1)): the two agree after the substitution tau = log(1 + t).
  3. The C program against the Python reference.
  4. Theorem B: for a level set of multiplicity two, Lambda_0 = Lambda_1 = 0 and
     Lambda_2 = rest * Q_s * Q_(s+2) / (s (s + 2)).
  5. Theorem C: for a block of multiplicity delta, Lambda_j = 0 for j < w.
  6. Proposition F (Ferrers sets) and Theorem D (peeling), as identities modulo
     a prime l = 1 (mod N).
  7. Proposition E (product block), as an identity modulo a prime l = 1 (mod N).

Usage: python3 validate.py PATH_TO_window_BINARY
Standard library only; random seeds are fixed.
"""
import itertools
import random
import sys

import ho_lib as ho

so = ho.so


def rand_fib(rng, M, q, maxc, dens):
    fib = {}
    for a in range(M):
        if rng.random() < dens:
            fib[a] = sorted(rng.sample(range(q), rng.randint(1, maxc)))
    return fib


def log_series(q):
    """tau = log(1 + t) modulo (q, t^(q-1)): coefficients of t^0 .. t^(q-2)."""
    return [0] + [(-1) ** (k - 1) * pow(k, q - 2, q) % q for k in range(1, q - 1)]


def singular_sets(ch, sizes, limit):
    out = []
    for r in sizes:
        for S in itertools.combinations(range(ch.M), r):
            if ch.minor(list(S), list(S)) < 0:
                out.append(list(S))
                if len(out) >= limit:
                    return out
    return out


def block_fib(rng, M, q, S, s, delta, pin=0.6):
    """Fibres with T_m = S exactly for s <= m < s + delta."""
    outside = [a for a in range(M) if a not in S]
    while True:
        fib = {}
        for a in S:
            c = s + delta if rng.random() < pin else rng.randint(s + delta, q)
            fib[a] = sorted(rng.sample(range(q), c))
        for a in outside:
            c = s if rng.random() < pin else rng.randint(0, s)
            if c:
                fib[a] = sorted(rng.sample(range(q), c))
        T = so.level_sets(fib, q)
        if [m for m in range(q) if T[m] == sorted(S)] == list(range(s, s + delta)):
            return fib


def block_weight(q, s, delta):
    """(w, m_max) of Theorem C for a block at levels s .. s + delta - 1."""
    if s == 0 and s + delta == q:
        return None, None
    if s == 0:
        ms = [delta - i for i in range(delta)]
    elif s + delta == q:
        ms = [i + 1 for i in range(delta)]
    else:
        ms = [min(i + 1, delta - i) for i in range(delta)]
    return sum(ms), max(ms)


def main():
    binary = sys.argv[1]
    rng = random.Random(20261128)

    # 1 and 3
    print("1. Window polynomial against the second-order report (coefficients of tau^e and tau^(e+1)); 3. C against Python")
    tot = [0, 0, 0, 0]
    for M, q in [(6, 7), (10, 3), (10, 11), (14, 13), (14, 29), (15, 7), (21, 5)]:
        charts = so.charts(M, q)
        L = min(q - 1, 8)
        n = a0 = a1 = nz = cp = 0
        fibs, chosen = [], []
        while len(fibs) < 60:
            fib = rand_fib(rng, M, q, min(3, q), 0.5)
            if fib:
                fibs.append(fib)
                chosen.append(rng.choice(charts))
        for ch in charts:
            idx = [i for i in range(len(fibs)) if chosen[i] is ch]
            dets = ho.window_c(binary, ch, [fibs[i] for i in idx], L)
            for i, d in zip(idx, dets):
                lam = ho.lam_coeffs(ch, fibs[i], L, d)
                e, pre, A0, Lam = so.second_order(ch, fibs[i])
                n += 1
                a0 += lam[0] == A0
                a1 += L < 2 or lam[1] == Lam
                nz += A0 >= 0
                if n % 4 == 0:
                    cp += d == ho.series_det(ch.K, ho.build_H(ch, fibs[i], L), L)
                else:
                    cp += 1
        print("   (%d, %d): %d configurations (%d with A_0 != 0); Lambda_0 = A_0 in %d, Lambda_1 = Lambda in %d; C = Python in %d of %d compared" % (M, q, n, nz, a0, a1, cp - (n - n // 4), n // 4))
        assert a0 == n and a1 == n and cp == n
        tot[0] += n
        tot[1] += nz
        tot[2] += n // 4
    print("   total %d configurations, %d with A_0 != 0, %d C/Python comparisons: all agree" % (tot[0], tot[1], tot[2]))

    # 2
    print("2. Window polynomial against direct determinants with zeta = 1 + t (all coefficients below t^(q-1))")
    total = informative = 0
    for M, q in [(6, 7), (10, 11), (14, 13), (15, 7), (21, 5), (10, 3)]:
        charts = so.charts(M, q)
        logt = log_series(q)
        n = inf = 0
        while n < 30:
            fib = rand_fib(rng, M, q, 2, 0.35)
            size = sum(len(v) for v in fib.values())
            e = ho.excess(fib)
            if not fib or size > 9 or e > q - 2:
                continue
            ch = rng.choice(charts)
            K = ch.K
            L = q - 1
            detH = ho.series_det(K, ho.build_H(ch, fib, L), L)
            v = 1
            for b in fib.values():
                for i in range(len(b)):
                    for j in range(i + 1, len(b)):
                        v = v * (b[j] - b[i]) % q
            vn = K.from_int(v)
            P = [-1] * L                      # P(tau) modulo tau^(q-1)
            for j in range(L - e):
                P[e + j] = K.mul(vn, detH[j])
            # substitute tau = log(1 + t)
            acc = [-1] * L
            power = [1] + [0] * (L - 1)       # tau^w as integers mod q
            for w in range(L):
                if P[w] >= 0:
                    for i in range(L):
                        if power[i]:
                            acc[i] = K.add(acc[i], K.mul(P[w], K.from_int(power[i])))
                nxt = [0] * L
                for i in range(L):
                    if power[i]:
                        for k in range(L - i):
                            nxt[i + k] = (nxt[i + k] + power[i] * logt[k]) % q
                power = nxt
            direct = so.direct_series(ch, fib)
            assert acc == direct
            n += 1
            inf += sum(1 for x in direct if x >= 0) > 1
        print("   (%d, %d): %d configurations agree in all %d coefficients (%d with more than one nonzero coefficient)" % (M, q, n, q - 1, inf))
        total += n
        informative += inf
    print("   total %d configurations: all agree" % total)

    # 4
    print("4. Theorem B on level sets of multiplicity two")
    total = nonzero = 0
    for M, q, sizes in [(15, 7, [6]), (10, 11, [3, 4]), (14, 13, [3, 4, 5]), (14, 29, [3, 4])]:
        charts = so.charts(M, q)
        n = nz = 0
        for ch in charts:
            carriers = singular_sets(ch, sizes, 3)
            fibs, ss = [], []
            for S in carriers:
                for _ in range(12):
                    s = rng.randint(1, min(q - 3, 4))
                    fibs.append(block_fib(rng, M, q, S, s, 2))
                    ss.append(s)
            if not fibs:
                continue
            L = 4
            for fib, s, d in zip(fibs, ss, ho.window_c(binary, ch, fibs, L)):
                lam = ho.lam_coeffs(ch, fib, L, d)
                pred = ho.separated_raises(ch, fib, [s, s + 2])
                assert lam[0] < 0 and lam[1] < 0 and lam[2] == pred
                n += 1
                nz += pred >= 0
        print("   (%d, %d): %d (configuration, chart) pairs: Lambda_0 = Lambda_1 = 0 and Lambda_2 equals the product formula in all; nonzero in %d" % (M, q, n, nz))
        total += n
        nonzero += nz
    print("   total %d pairs, product nonzero in %d" % (total, nonzero))

    # 5
    print("5. Theorem C: order of vanishing for a block of multiplicity delta on base 15 above 7 (both charts), and base 10 above 11")
    C = [3, 5, 6, 9, 10, 12]
    total = 0
    for M, q, S_of in [(15, 7, lambda ch: C), (10, 11, lambda ch: singular_sets(ch, [3, 4], 1)[0])]:
        charts = so.charts(M, q)
        for s in range(0, q):
            for delta in range(1, q - s + 1):
                w, mmax = block_weight(q, s, delta)
                if w is None or w > 16 or (M == 10 and (delta > 3 or s > 4)):
                    continue
                L = w + 1
                n = lead = 0
                for ch in charts:
                    S = S_of(ch)
                    fibs = [block_fib(rng, M, q, S, s, delta) for _ in range(8)]
                    for fib, d in zip(fibs, ho.window_c(binary, ch, fibs, L)):
                        assert all(x < 0 for x in d[:w])
                        n += 1
                        lead += d[w] >= 0
                total += n
                print("   (%d, %d) block at levels %d..%d: w = %d; %d (configuration, chart) pairs with Lambda_j = 0 for all j < w; Lambda_w != 0 in %d" % (M, q, s, s + delta - 1, w, n, lead))
    print("   total %d pairs" % total)

    # 6
    print("6. Proposition F and Theorem D modulo a prime l = 1 (mod N)")
    for M, Mp in [(15, 7), (21, 5), (35, 3), (10, 11), (7, 15)]:
        mod = ho.ModN(M, Mp)
        l = mod.l
        ok = 0
        for _ in range(30):
            order = rng.sample(range(Mp), Mp)
            size = {a: rng.randint(0, Mp) for a in range(M)}
            K = [(a, order[i]) for a in range(M) for i in range(size[a])]
            T = [[a for a in range(M) if size[a] > m] for m in range(Mp)]
            Delta = [mod.delta(order[:m]) for m in range(Mp + 1)]
            pred = 1
            for m in range(Mp):
                pred = pred * mod.d(T[m]) % l * pow(Delta[m + 1] * pow(Delta[m], l - 2, l) % l, len(T[m]), l) % l
            D = mod.minor(K)
            assert D == pred and D != 0
            ok += 1
        print("   Ferrers sets in Z_%d x Z_%d: product formula holds and is nonzero modulo %d in %d of %d" % (M, Mp, l, ok, 30))
    for M, q in [(15, 7), (10, 11), (21, 5)]:
        mod = ho.ModN(M, q)
        l = mod.l
        ok = 0
        for _ in range(30):
            r = rng.randint(1, q - 1)
            order = rng.sample(range(q), q)
            X, rest = order[:r], order[r:]
            size = {a: rng.randint(0, q) for a in range(M)}
            fib = {a: sorted(X[:min(c, r)] + rng.sample(rest, max(c - r, 0))) for a, c in size.items() if c}
            T = [[a for a in range(M) if size[a] > m] for m in range(q)]
            for j in range(r):
                assert sorted(a for a in fib if X[j] in fib[a]) == T[j]
            K = [(a, b) for a in sorted(fib) for b in fib[a]]
            Kt = [(a, b) for a in T[r] for b in sorted(set(fib[a]) | set(X))] if r < q else []
            lhs, rhs = mod.minor(K), mod.minor(Kt)
            for j in range(r):
                dj = mod.delta(X[:j + 1]) * pow(mod.delta(X[:j]), l - 2, l) % l
                lhs = lhs * pow(dj, len(T[r]), l) % l * mod.d(T[r]) % l
                rhs = rhs * pow(dj, len(T[j]), l) % l * mod.d(T[j]) % l
            assert lhs == rhs
            ok += 1
        print("   Peeling identity of Theorem D in Z_%d x Z_%d holds modulo %d in %d of %d random sets with r Ferrers levels" % (M, q, l, ok, 30))
    mod = ho.ModN(15, 7)
    ok = 0
    for _ in range(40):
        s = rng.randint(0, 5)
        r = s + 1
        order = rng.sample(range(7), 7)
        X, rest = order[:r], order[r:]
        size = {a: rng.randint(s + 1, 7) if a in C else rng.randint(0, s) for a in range(15)}
        size[rng.choice(C)] = s + 1
        fib = {a: sorted(X[:min(c, r)] + rng.sample(rest, max(c - r, 0))) for a, c in size.items() if c}
        assert so.level_sets(fib, 7)[s] == C
        K = [(a, b) for a in sorted(fib) for b in fib[a]]
        assert mod.minor(K) != 0
        ok += 1
    print("   Order 105: %d of %d random members with C = T_s and s + 1 Ferrers levels have a nonzero minor modulo %d (as Theorem D predicts)" % (ok, 40, mod.l))
    # 7
    print("7. Proposition E modulo a prime l = 1 (mod N)")
    for M, q in [(15, 7), (10, 11), (21, 5)]:
        mod = ho.ModN(M, q)
        l = mod.l
        ok = 0
        for _ in range(25):
            S = sorted(rng.sample(range(M), rng.randint(1, M - 1)))
            B = sorted(rng.sample(range(q), rng.randint(1, q - 1)))
            fib = {}
            for a in range(M):
                v = sorted({b for b in range(q) if rng.random() < 0.4} | (set(B) if a in S else set()))
                if v:
                    fib[a] = v
            K = [(a, b) for a in sorted(fib) for b in fib[a]]
            rest = [(a, b) for (a, b) in K if not (a in S and b in B)]
            FS = ho.inv_mod([[mod.ep[a * a2 % M] for a2 in S] for a in S], l)
            ZB = ho.inv_mod([[mod.zp[b * b2 % q] for b2 in B] for b in B], l)
            PS, PB = {}, {}
            mat = []
            for a, b in rest:
                row = []
                for a2, b2 in rest:
                    if (a, a2) not in PS:
                        PS[(a, a2)] = sum(mod.ep[a * S[i] % M] * FS[i][j] % l * mod.ep[S[j] * a2 % M] for i in range(len(S)) for j in range(len(S))) % l
                    if (b, b2) not in PB:
                        PB[(b, b2)] = sum(mod.zp[b * B[i] % q] * ZB[i][j] % l * mod.zp[B[j] * b2 % q] for i in range(len(B)) for j in range(len(B))) % l
                    f, z = mod.ep[a * a2 % M], mod.zp[b * b2 % q]
                    row.append((f * (z - PB[(b, b2)]) + (f - PS[(a, a2)]) * PB[(b, b2)]) % l)
                mat.append(row)
            lhs = mod.minor(K)
            rhs = pow(mod.d(S), len(B), l) * pow(mod.delta(B), len(S), l) % l * ho.det_mod(mat, l) % l
            assert lhs == rhs
            ok += 1
        print("   Z_%d x Z_%d: product-block identity holds modulo %d in %d of 25 random sets containing S x B" % (M, q, l, ok))
    print("all checks passed")


if __name__ == "__main__":
    main()
