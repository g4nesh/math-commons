#!/usr/bin/env python3
"""Where the window polynomial is blind at order 105: structured members.

Z_105 = Z_15 x Z_7, C = {3, 5, 6, 9, 10, 12}; members with C = T_s of
multiplicity one.  For four structured families the script tabulates, for
every s, the least order of a nonzero window coefficient (over the two charts,
orders 0..13) against two structural properties:

  double block  K contains C x B with |B| = s + 1, AND the fibres outside C lie
                in a set B' of size s.
  nested        a double block with B' inside B: every fibre outside C is
                contained in every fibre over C.  Then the minor is nonzero by
                Corollary G1 of README.md.  Every member to which Theorem D
                applies (K or its complement peelable) is nested.

Families (C at level s, everything else random):
  one    one fibre of size s + 1 over C and one fibre of size s outside C;
         all other fibres over C have size >= s + 2, outside C size <= s - 1;
  const  all fibres of size s + 1 over C are equal and all fibres of size s
         outside C are equal;
  block  K contains C x B for a random B of size s + 1, otherwise uniform;
  double K contains C x B and its fibres outside C lie in B', for random B, B'
         of sizes s + 1 and s.

The script also counts exactly the members of multiplicity one that have a
double block.

Every member whose window coefficients of order <= 5 all vanish and which is
not nested is undecided by this report.  For each of them the minor is
evaluated modulo a prime l = 1 (mod 105); a nonzero value shows that it is not
a zero of the determinant.

Usage: python3 blind105.py PATH_TO_window_BINARY [SAMPLES_PER_CELL]   (default 1500)
Standard library only; the random seed is fixed.
"""
import random
import sys
from math import comb, log2

import ho_lib as ho

so = ho.so
M, Q = 15, 7
C = [3, 5, 6, 9, 10, 12]
OUT = [a for a in range(M) if a not in C]


def first_orders(binary, charts, fibs, L):
    per = [[next((j for j, x in enumerate(d) if x >= 0), None) for d in ho.window_c(binary, ch, fibs, L)] for ch in charts]
    return [min((v for v in t if v is not None), default=None) for t in zip(*per)]


def peelable(fib, s):
    T = so.level_sets(fib, Q)
    sl = [sorted(a for a in fib if x in fib[a]) for x in range(Q)]
    for levels, need in ((range(Q), s + 1), (range(Q - 1, -1, -1), Q - s)):
        used, r = set(), 0
        for m in levels:
            x = next((x for x in range(Q) if x not in used and sl[x] == T[m]), None)
            if x is None:
                break
            used.add(x)
            r += 1
        if r >= need:
            return True
    return False


def double_block(fib, s):
    common = set(range(Q))
    for a in C:
        common &= set(fib.get(a, []))
    union = set()
    for a in OUT:
        union |= set(fib.get(a, []))
    return len(common) >= s + 1 and len(union) <= s


def nested(fib, s):
    common = set(range(Q))
    for a in C:
        common &= set(fib.get(a, []))
    union = set()
    for a in OUT:
        union |= set(fib.get(a, []))
    return len(common) >= s + 1 and len(union) <= s and union <= common


def fam_one(rng, s):
    size = {a: rng.randint(min(s + 2, Q), Q) for a in C}
    size[rng.choice(C)] = s + 1
    size.update({a: rng.randint(0, max(s - 1, 0)) for a in OUT})
    if s >= 1:
        size[rng.choice(OUT)] = s
    return {a: sorted(rng.sample(range(Q), c)) for a, c in size.items() if c}


def fam_const(rng, s):
    size = {a: rng.randint(s + 1, Q) for a in C}
    for a in rng.sample(C, rng.randint(1, 6)):
        size[a] = s + 1
    size.update({a: rng.randint(0, s) for a in OUT})
    if s >= 1:
        for a in rng.sample(OUT, rng.randint(1, 9)):
            size[a] = s
    BU = sorted(rng.sample(range(Q), s + 1))
    BL = sorted(rng.sample(range(Q), s))
    fib = {}
    for a, c in size.items():
        if c:
            fib[a] = BU if (a in C and c == s + 1) else BL if (a not in C and c == s) else sorted(rng.sample(range(Q), c))
    return fib


def fam_block(rng, s):
    B = rng.sample(range(Q), s + 1)
    rest = [b for b in range(Q) if b not in B]
    size = {a: rng.randint(s + 1, Q) for a in C}
    size[rng.choice(C)] = s + 1
    size.update({a: rng.randint(0, s) for a in OUT})
    if s >= 1:
        size[rng.choice(OUT)] = s
    return {a: sorted(B + rng.sample(rest, c - s - 1)) if a in C else sorted(rng.sample(range(Q), c)) for a, c in size.items() if c}


def fam_double(rng, s):
    B = rng.sample(range(Q), s + 1)
    rest = [b for b in range(Q) if b not in B]
    Bp = rng.sample(range(Q), s)
    size = {a: rng.randint(s + 1, Q) for a in C}
    size[rng.choice(C)] = s + 1
    size.update({a: rng.randint(0, s) for a in OUT})
    if s >= 1:
        size[rng.choice(OUT)] = s
    return {a: sorted(B + rng.sample(rest, c - s - 1)) if a in C else sorted(rng.sample(Bp, c)) for a, c in size.items() if c}


def count_double(s):
    """Number of K with C = T_s of multiplicity one and a double block."""
    over = comb(Q, s + 1) * ((2 ** (6 - s)) ** 6 - (2 ** (6 - s) - 1) ** 6)   # fibres over C contain B, one equals B
    outside = comb(Q, s) * ((2 ** s) ** 9 - (2 ** s - 1) ** 9) if s else 1     # fibres outside C inside B', one equals B'
    return over * outside


def main():
    binary = sys.argv[1]
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 1500
    rng = random.Random(1055)
    charts = so.charts(M, Q)
    mod = ho.ModN(M, Q)
    undecided = zeros = 0
    blind_without_block = block_not_blind = 0
    example = None
    for name, fam in (("one", fam_one), ("const", fam_const), ("block", fam_block), ("double", fam_double)):
        print("family %s, %d members per level s; columns: least order <= 5 | least order 6..13 | all zero to order 13" % (name, n))
        for s in range(Q):
            fibs = [fam(rng, s) for _ in range(n)]
            for f in fibs[:50]:
                assert so.level_sets(f, Q)[s] == C
            orders = first_orders(binary, charts, fibs, 14)
            tab = {}
            for f, o in zip(fibs, orders):
                col = 0 if (o is not None and o <= 5) else 1 if o is not None else 2
                pe, pb = nested(f, s), double_block(f, s)
                assert pe or not peelable(f, s)
                row = "nested" if pe else "double block, not nested" if pb else "no double block"
                tab.setdefault(row, [0, 0, 0])[col] += 1
                if col == 2 and not pb:
                    blind_without_block += 1
                if col != 2 and pb:
                    block_not_blind += 1
                if col and not pe:
                    undecided += 1
                    if mod.minor([(a, b) for a in sorted(f) for b in f[a]]) == 0:
                        zeros += 1
                    if example is None and col == 2 and s == 1 and name == "one":
                        example = f
            print("   s = %d: " % s + "; ".join("%s: %d | %d | %d" % (row, *tab[row]) for row in ("nested", "double block, not nested", "no double block") if row in tab))
    print("members undecided by Theorem A and Corollary G1 (no nonzero coefficient of order <= 5, not nested): %d; of these, minor zero modulo %d: %d" % (undecided, mod.l, zeros))
    print("members with all coefficients to order 13 zero and no double block: %d" % blind_without_block)
    print("members with a double block and some nonzero coefficient of order <= 13: %d" % block_not_blind)
    counts = [count_double(s) for s in range(Q)]
    print("exact number of members with C = T_s of multiplicity one and a double block, s = 0..6: %s; total %d (2^%.2f)" % (counts, sum(counts), log2(sum(counts))))
    nest = [count_double(s) * (s + 1) // comb(Q, s) for s in range(Q)]
    for s in range(Q):
        assert count_double(s) * (s + 1) % comb(Q, s) == 0
    print("of these, nested (B' one of the s + 1 subsets of B of size s), s = 0..6: %s; total %d (2^%.2f)" % (nest, sum(nest), log2(sum(nest))))
    if example:
        print("example of an undecided member (family one, s = 1): fibres %s" % {a: example[a] for a in sorted(example)})


if __name__ == "__main__":
    main()
