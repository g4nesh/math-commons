#!/usr/bin/env python3
"""Theorems A-D of README.md on the family that P(105) reduces to.

Z_105 = Z_15 x Z_7, C = {3, 5, 6, 9, 10, 12}.  The family consists of the index
sets K with T_m = C for some level m.  For such K the levels with T_m = C form
a block s, ..., s + delta - 1; fibres over C have size >= s + delta and fibres
outside C have size <= s.

  0. Census of the family by (s, delta), counting every index set once.
  1. Multiplicity one: uniform sample; the second-order term first (through the
     tables of the second-order report), then the window polynomial for the
     members it leaves.
  2. Multiplicity two and more: uniform samples of every class; least order of
     a nonzero window coefficient over the two charts, against the range in
     which Theorems A and C prove that the minor is nonzero.
  3. Totals.
  4. Structured members: Ferrers shapes, members with one fibre in each shell,
     and one explicit member that no criterion of this report decides.

Usage: python3 apply105.py PATH_TO_window_BINARY [N1] [N2]
       N1 = sample size for step 1 (default 3000000), N2 = sample size per large
       class in step 2 (default 20000).
Standard library only; random seeds are fixed.
"""
import importlib.util
import itertools
import os
import random
import sys
from math import comb, log2

import ho_lib as ho

so = ho.so
spec = importlib.util.spec_from_file_location(
    "so_apply105", os.path.join(ho.HERE, "..", "..", "second-order", "ganesh-claude-minors", "apply105.py"))
so105 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(so105)

M, Q = 15, 7
C = so105.C
OUT = so105.OUT
W = so105.W
SUBSETS = {c: list(itertools.combinations(range(Q), c)) for c in range(Q + 1)}


def class_size(s, d):
    return (W(s + d, 7) ** 6 - W(s + d + 1, 7) ** 6) * (W(0, s) ** 9 - W(0, s - 1) ** 9)


def block_weight(s, d):
    if s == 0 and s + d == Q:
        return None, None
    if s == 0:
        ms = [d - i for i in range(d)]
    elif s + d == Q:
        ms = [i + 1 for i in range(d)]
    else:
        ms = [min(i + 1, d - i) for i in range(d)]
    return sum(ms), max(ms)


def sample_sizes(rng, s, d):
    hi = [c for c in range(s + d, 8) for _ in range(comb(7, c))]
    lo = [c for c in range(0, s + 1) for _ in range(comb(7, c))]
    while True:
        size = {a: rng.choice(hi) for a in C}
        if s + d > 7 or any(size[a] == s + d for a in C):
            break
    while True:
        low = {a: rng.choice(lo) for a in OUT}
        if s == 0 or any(low[a] == s for a in OUT):
            break
    size.update(low)
    return size


def sample_member(rng, s, d):
    size = sample_sizes(rng, s, d)
    return {a: list(rng.choice(SUBSETS[c])) for a, c in size.items() if c}


def first_orders(binary, charts, fibs, L):
    """For each set: the least j < L with Lambda_j != 0 in some chart, or None."""
    per = [[next((j for j, x in enumerate(d) if x >= 0), None) for d in ho.window_c(binary, ch, fibs, L)] for ch in charts]
    return [min((v for v in t if v is not None), default=None) for t in zip(*per)]


def peel_depths(fib):
    """(bottom, top): the number of consecutive level sets T_0, T_1, ... (respectively
    T_6, T_5, ...) that occur as slices {a : x in K_a} for distinct x."""
    T = so.level_sets(fib, Q)
    sl = [sorted(a for a in fib if x in fib[a]) for x in range(Q)]
    out = []
    for levels in (range(Q), range(Q - 1, -1, -1)):
        used, r = set(), 0
        for m in levels:
            x = next((x for x in range(Q) if x not in used and sl[x] == T[m]), None)
            if x is None:
                break
            used.add(x)
            r += 1
        out.append(r)
    return tuple(out)


def peelable(fib, s, d):
    """Theorem D applies to K (bottom) or to its complement (top)."""
    bottom, top = peel_depths(fib)
    return bottom >= s + d or top >= Q - s


def show(dist):
    return ", ".join("%s: %d" % ("none" if k is None else k, v) for k, v in sorted(dist.items(), key=lambda kv: (kv[0] is None, kv[0])))


def main():
    binary = sys.argv[1]
    n1 = int(sys.argv[2]) if len(sys.argv) > 2 else 3000000
    n2 = int(sys.argv[3]) if len(sys.argv) > 3 else 20000
    charts = so.charts(M, Q)
    assert len(charts) == 2 and all(ch.minor(C, C) < 0 for ch in charts)
    cs = set(C)
    assert all(cs & {(a + t) % M for a in C} for t in range(M))   # no translate of C is disjoint from C

    # 0. census
    classes = {(s, d): class_size(s, d) for s in range(Q) for d in range(1, Q - s + 1)}
    fam = sum(classes.values())
    pairs = sum(W(t, 7) ** 6 * W(0, t - 1) ** 9 for t in range(1, 8))
    print("0. Census.  Index sets K with some level set equal to C: %d (2^%.3f)." % (fam, log2(fam)))
    print("   Counted with the multiplicity of C as a level set: %d (2^%.3f); this is the figure the second-order report called the family size." % (pairs, log2(pairs)))
    by_mult = {}
    for (s, d), n in classes.items():
        by_mult[d] = by_mult.get(d, 0) + n
    print("   Share by multiplicity delta of C: " + "; ".join("%d: %.3g%%" % (d, 100.0 * by_mult[d] / fam) for d in sorted(by_mult)))
    print("   Largest classes (s, delta): " + "; ".join("(%d, %d): %.4g%%" % (s, d, 100.0 * n / fam) for (s, d), n in sorted(classes.items(), key=lambda kv: -kv[1])[:8]))

    # 1. multiplicity one
    rng = random.Random(1051)
    tabs = [so105.Tables(ch) for ch in charts]
    ss = list(range(Q))
    wts = [classes[(s, 1)] for s in ss]
    left, n_second = [], 0
    for _ in range(n1):
        s = rng.choices(ss, weights=wts)[0]
        size = sample_sizes(rng, s, 1)
        fib = {a: rng.choice(SUBSETS[c]) for a, c in size.items() if c}
        t = s + 1
        U = tuple(a for a in C if size[a] == t) if t <= 6 else ()
        Lo = tuple(a for a in OUT if size[a] == t - 1) if t >= 2 else ()
        sig = {a: sum(fib[a]) % Q for a in U + Lo}
        if any(tab.bracket(t, U, Lo, sig) >= 0 for tab in tabs):
            n_second += 1
        else:
            left.append((s, {a: list(v) for a, v in fib.items()}))
    orders = first_orders(binary, charts, [f for _, f in left], Q - 1)
    dist = {}
    for o in orders:
        dist[o] = dist.get(o, 0) + 1
    unresolved1 = [(s, f) for (s, f), o in zip(left, orders) if o is None]
    assert all(o is None or o >= 2 for o in orders)
    print("1. Multiplicity one (%.4f%% of the family), uniform sample of %d members:" % (100.0 * by_mult[1] / fam, n1))
    print("   second-order term nonzero in some chart: %d (%.4f%%); not: %d (%.4f%%)." % (n_second, 100.0 * n_second / n1, len(left), 100.0 * len(left) / n1))
    print("   For those %d, least order j <= 5 of a nonzero window coefficient over the two charts: %s." % (len(left), show(dist)))
    print("   Left by Theorem A: %d of %d; of these, %d are peelable (Theorem D)." % (len(unresolved1), n1, sum(peelable(f, s, 1) for s, f in unresolved1)))

    # 2. multiplicity two and more
    print("2. Multiplicity two and more, by class (s, delta); w and the last order j_max covered by Theorems A and C:")
    rng = random.Random(1052)
    rows = []
    for (s, d), n in sorted(classes.items(), key=lambda kv: -kv[1]):
        if d == 1:
            continue
        w, mmax = block_weight(s, d)
        share = n / fam
        if w is None:
            print("   (%d, %d): share %.3g%%: the product set C x Z_7; nonzero by Proposition F." % (s, d, 100.0 * share))
            rows.append((share, 1, 1, 0))
            continue
        if mmax >= Q - 1:
            print("   (%d, %d): share %.3g%%: w = %d, m_max = %d = q - 1; outside Theorem C, not sampled." % (s, d, 100.0 * share, w, mmax))
            rows.append((share, 1, 0, 1))
            continue
        jmax = max(Q - 2, w + Q - 2 - mmax)
        size = n2 if share > 1e-5 else max(n2 // 20, 200)
        fibs = [sample_member(rng, s, d) for _ in range(size)]
        orders = first_orders(binary, charts, fibs, jmax + 1)
        dist = {}
        for o in orders:
            dist[o] = dist.get(o, 0) + 1
        bad = [f for f, o in zip(fibs, orders) if o is None]
        assert all(o is None or o >= w for o in orders)
        rows.append((share, size, size - len(bad), len(bad)))
        print("   (%d, %d): share %.3g%%; w = %d, j_max = %d; sample %d; least order: %s; left: %d (peelable: %d)" % (
            s, d, 100.0 * share, w, jmax, size, show(dist), len(bad), sum(peelable(f, s, d) for f in bad)))

    # 3. totals
    und = by_mult[1] / fam * len(unresolved1) / n1 + sum(sh * bad / size for sh, size, ok, bad in rows if size > 1)
    not_sampled = sum(sh for sh, size, ok, bad in rows if size == 1 and bad == 1)
    total_samples = n1 + sum(size for sh, size, ok, bad in rows if size > 1)
    total_left = len(unresolved1) + sum(bad for sh, size, ok, bad in rows if size > 1)
    main_n = n1 + sum(size for sh, size, ok, bad in rows if size > 1 and sh > 1e-5)
    main_left = len(unresolved1) + sum(bad for sh, size, ok, bad in rows if size > 1 and sh > 1e-5)
    print("3. Totals.  Sampled members: %d; left undecided by Theorems A and C: %d, of which %d lie in the classes with share above 0.001%% (%d samples)." % (
        total_samples, total_left, main_left, main_n))
    print("   Weighted by class size, the undecided share of the family in these samples is %.3g%%; the classes not sampled have total share %.3g%%." % (100.0 * und, 100.0 * not_sampled))
    print("   A uniform sample cannot see thin structured subfamilies; see step 4.")

    # 4. structured members
    print("4. Structured members.")
    rng = random.Random(1054)
    fibs = []
    for _ in range(200):
        s = rng.randint(0, 6)
        order = rng.sample(range(Q), Q)
        size = {a: rng.randint(s + 1, 7) if a in C else rng.randint(0, s) for a in range(M)}
        size[rng.choice(C)] = s + 1
        fibs.append({a: sorted(order[:c]) for a, c in size.items() if c})
    orders = first_orders(binary, charts, fibs, 14)
    mod = ho.ModN(M, Q)
    nonzero = sum(mod.minor([(a, b) for a in sorted(f) for b in f[a]]) != 0 for f in fibs)
    print("   Ferrers members (nested fibres), 200 random: window coefficients of order < 14 all zero in both charts for %d; minor nonzero modulo %d for %d (Proposition F)." % (
        sum(o is None for o in orders), mod.l, nonzero))
    print("   Members with a single fibre in each shell next to C (C of multiplicity one at level s), 1500 random per s, window to order 13:")
    example = None
    for s in range(Q):
        fibs = []
        for _ in range(1500):
            size = {a: rng.randint(min(s + 2, Q), Q) for a in C}
            size[rng.choice(C)] = s + 1
            size.update({a: rng.randint(0, max(s - 1, 0)) for a in OUT})
            if s >= 1:
                size[rng.choice(OUT)] = s
            fibs.append({a: sorted(rng.sample(range(Q), c)) for a, c in size.items() if c})
        orders = first_orders(binary, charts, fibs, 14)
        tab = {}
        for f, o in zip(fibs, orders):
            key = ("order <= 5" if o is not None and o <= 5 else "order 6..13" if o is not None else "all zero", "peelable" if peelable(f, s, 1) else "not peelable")
            tab[key] = tab.get(key, 0) + 1
            if o is None and not peelable(f, s, 1) and example is None:
                example = (s, f)
        print("     s = %d: %s" % (s, "; ".join("%s, %s: %d" % (k[0], k[1], v) for k, v in sorted(tab.items()))))
    if example:
        s, f = example
        D = mod.minor([(a, b) for a in sorted(f) for b in f[a]])
        print("   An explicit member decided by none of Theorems A-D: C = T_%d, fibres %s." % (s, {a: f[a] for a in sorted(f)}))
        print("   Its window coefficients of order < 14 vanish in both charts and it is not peelable; its minor is %s modulo %d, so it is not a zero of the determinant." % ("nonzero" if D else "ZERO", mod.l))


if __name__ == "__main__":
    main()
