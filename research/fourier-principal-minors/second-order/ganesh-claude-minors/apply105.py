#!/usr/bin/env python3
"""The second-order criterion on the family that P(105) reduces to.

Z_105 = Z_15 x Z_7, C = {3, 5, 6, 9, 10, 12} (the punctured cross).  The family
consists of the K whose level set {a : c_a >= t} equals C for some t in 1..7.
For such K put
    U = {a in C : c_a = t}          (upper shell),
    L = {a not in C : c_a = t - 1}  (lower shell, for t >= 2),
so that T_(t-1) = C, T_t = C minus U and T_(t-2) = C u L.  With sigma_a the sum
of the fibre over a (mod 7), Theorem 1 of README.md gives, in each of the two
charts above 7,
    Lambda = D_rest * ( t^-1 d(C u L) Q_U(sigma_U) + (t-1)^-1 d(C minus U) Q_L(sigma_L) ),
where D_rest is the product of d(T_m) over the levels other than t-2, t-1, t,
the first term is absent for t = 7 and has d(C u L) replaced by 1 for t = 1, and
the second is absent for t = 1 and has d(C minus U) replaced by 1 for t = 7.
K is "resolved" when Lambda != 0 in at least one chart; then its minor is nonzero.

The script
  0. recomputes, by brute force in Python, that the singular subsets of Z_15 above
     7 are the 30 translates of C and of its complement, in both charts;
  1. checks the bracket formula against so_lib.second_order on random members;
  2. treats the single-layer lifts (t = 1, all fibres of size one) exhaustively;
  3. counts exactly the resolved members of the slice t = 1 (K inside C x Z_7);
  4. computes exactly the share of the family in which C is a level set of
     multiplicity at least two (there Lambda vanishes identically);
  5. samples the family uniformly and classifies the samples.

Usage: python3 apply105.py [SAMPLES]      (default 300000; about six minutes)
Standard library only.
"""
import itertools
import random
import sys
from math import comb

import so_lib as so

M, Q7 = 15, 7
C = [3, 5, 6, 9, 10, 12]
OUT = [a for a in range(M) if a not in C]


def W(lo, hi):
    """Number of subsets of Z_7 with size in [lo, hi]."""
    return sum(comb(7, c) for c in range(max(lo, 0), min(hi, 7) + 1))


class Tables:
    def __init__(self, ch):
        self.ch, K = ch, ch.K
        self.K = K
        assert ch.minor(C, C) < 0
        self.QU, self.dU, self.QL, self.dL = {}, {}, {}, {}
        for r in range(0, 7):
            for U in itertools.combinations(C, r):
                rest = [a for a in C if a not in U]
                self.dU[U] = ch.minor(rest, rest)
                self.QU[U] = {(a, a2): K.mul(ch.cofactor(C, a, a2), ch.bordered(rest, a, a2)) for a in U for a2 in U}
        for r in range(0, 10):
            for L in itertools.combinations(OUT, r):
                big = sorted(C + list(L))
                self.dL[L] = ch.minor(big, big)
                self.QL[L] = {(a, a2): K.mul(ch.cofactor(big, a, a2), ch.bordered(C, a, a2)) for a in L for a2 in L}

    def form(self, Qm, sig):
        K, acc = self.K, -1
        for (a, a2), coef in Qm.items():
            s = sig[a] * sig[a2] % Q7
            if s and coef >= 0:
                acc = K.add(acc, K.mul(K.from_int(s), coef))
        return acc

    def bracket(self, t, U, L, sig):
        K = self.K
        up = low = -1
        if t <= 6 and U:
            fac = 0 if t == 1 else self.dL[L]
            up = K.mul(K.mul(K.from_int(pow(t, 5, Q7)), fac), self.form(self.QU[U], sig))
        if t >= 2 and L:
            fac = 0 if t == 7 else self.dU[U]
            low = K.mul(K.mul(K.from_int(pow(t - 1, 5, Q7)), fac), self.form(self.QL[L], sig))
        return K.add(up, low)


def sample_sizes(t, rng):
    hi = [c for c in range(t, 8) for _ in range(comb(7, c))]
    lo = [c for c in range(0, t) for _ in range(comb(7, c))]
    return {a: rng.choice(hi) if a in C else rng.choice(lo) for a in range(M)}


def main():
    n_samples = int(sys.argv[1]) if len(sys.argv) > 1 else 300000
    charts = so.charts(M, Q7)
    assert len(charts) == 2

    # 0. singular subsets of Z_15 above 7
    full = set(range(M))
    expected = {frozenset((a + s) % M for a in C) for s in range(M)}
    expected |= {frozenset(full - S) for S in expected}
    for ch in charts:
        found = set()
        for mask in range(1, 1 << M):
            S = [a for a in range(M) if mask >> a & 1]
            if ch.minor(S, S) < 0:
                found.add(frozenset(S))
        assert found == expected
    singular = expected
    print("0. In both charts above 7 the singular subsets of Z_15 are the 30 translates of C and of its complement (brute force over 2^15 subsets).")

    tabs = [Tables(ch) for ch in charts]

    # 1. bracket formula against the general second-order term
    rng = random.Random(105)
    checked = 0
    for _ in range(60):
        t = rng.randint(1, 3)
        size = sample_sizes(t, rng)
        fib = {a: sorted(rng.sample(range(Q7), c)) for a, c in size.items() if c}
        U = tuple(a for a in C if size[a] == t)
        L = tuple(a for a in OUT if size[a] == t - 1) if t >= 2 else ()
        sig = {a: sum(v) % Q7 for a, v in fib.items()}
        T = so.level_sets(fib, Q7)
        for tab in tabs:
            ch, K = tab.ch, tab.K
            e, pre, A0, Lam = so.second_order(ch, fib)
            assert A0 < 0
            rest = 0
            for m in range(Q7):
                if m not in (t - 2, t - 1, t):
                    rest = K.mul(rest, ch.minor(T[m], T[m]))
            assert Lam == K.mul(rest, tab.bracket(t, U, L, sig))
            checked += 1
    print("1. Bracket formula equals the general second-order term on %d random (member, chart) pairs." % checked)

    # 2. single-layer lifts
    unresolved = [b for b in itertools.product(range(Q7), repeat=6)
                  if all(tab.bracket(1, tuple(C), (), dict(zip(C, b))) < 0 for tab in tabs)]
    per_chart = [sum(1 for b in itertools.product(range(Q7), repeat=6) if tab.bracket(1, tuple(C), (), dict(zip(C, b))) < 0) for tab in tabs]
    assert all(len(set(b)) == 1 for b in unresolved)
    print("2. Single-layer lifts: Lambda vanishes for %d and %d of the 7^6 vectors b in the two charts, and in both charts only for the %d constant b." % (per_chart[0], per_chart[1], len(unresolved)))

    # 3. the slice t = 1, exactly
    total = resolved = 0
    by_size = {}
    for r in range(0, 7):
        for U in itertools.combinations(C, r):
            weight = 120 ** (6 - r)
            n_res = 0
            for s in itertools.product(range(Q7), repeat=r):
                if r and any(tab.bracket(1, U, (), dict(zip(U, s))) >= 0 for tab in tabs):
                    n_res += 1
            total += weight * Q7 ** r
            resolved += weight * n_res
            a, b = by_size.get(r, (0, 0))
            by_size[r] = (a + n_res, b + Q7 ** r)
    assert total == 127 ** 6
    print("3. Slice t = 1 (K inside C x Z_7, all six fibres nonempty): %d sets, %d resolved (%.4f%%), %d not." % (total, resolved, 100.0 * resolved / total, total - resolved))
    print("   resolved (U, sigma_U) classes by |U|: " + ", ".join("|U|=%d: %d of %d" % (r, by_size[r][0], by_size[r][1]) for r in sorted(by_size)))

    # 4. multiplicity of C, exactly
    weights = {t: W(t, 7) ** 6 * W(0, t - 1) ** 9 for t in range(1, 8)}
    fam = sum(weights.values())
    mult2 = 0
    for t in range(1, 8):
        p_u = (W(t + 1, 7) / W(t, 7)) ** 6 if t <= 6 else 0.0   # U empty: C is also T_t
        p_l = (W(0, t - 2) / W(0, t - 1)) ** 9 if t >= 2 else 0.0  # L empty: C is also T_(t-2)
        mult2 += weights[t] * (1 - (1 - p_u) * (1 - p_l))
    print("4. Family size %d (2^%.2f).  C has multiplicity >= 2 as a level set in a share %.4f%% of it; there Lambda = 0 identically." % (fam, __import__("math").log2(fam), 100.0 * mult2 / fam))

    # 5. uniform sample of the family
    rng = random.Random(7105)
    ts = [t for t in range(1, 8)]
    cum = [weights[t] for t in ts]
    counts = {"resolved": 0, "multiplicity >= 2": 0, "another singular level": 0, "Lambda = 0 in both charts": 0}
    for _ in range(n_samples):
        t = rng.choices(ts, weights=cum)[0]
        size = sample_sizes(t, rng)
        U = tuple(a for a in C if size[a] == t)
        L = tuple(a for a in OUT if size[a] == t - 1) if t >= 2 else ()
        if (t <= 6 and not U) or (t >= 2 and not L):
            counts["multiplicity >= 2"] += 1
            continue
        levels = [frozenset(a for a in range(M) if size[a] > m) for m in range(Q7)]
        if any(levels[m] in singular for m in range(Q7) if m not in (t - 2, t - 1, t) and levels[m]):
            counts["another singular level"] += 1
            continue
        sig = {a: rng.randrange(Q7) for a in U + L}   # fibre sums of sizes 1..6 are uniform
        if any(tab.bracket(t, U, L, sig) >= 0 for tab in tabs):
            counts["resolved"] += 1
        else:
            counts["Lambda = 0 in both charts"] += 1
    print("5. Uniform sample of %d members: " % n_samples + "; ".join("%s: %d (%.4f%%)" % (k, v, 100.0 * v / n_samples) for k, v in counts.items()))


if __name__ == "__main__":
    main()
