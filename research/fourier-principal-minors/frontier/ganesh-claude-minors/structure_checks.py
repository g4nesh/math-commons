#!/usr/bin/env python3
"""Checks for the two structural statements in README.md.

1. Punctured star.  For pairwise coprime divisors d_1, ..., d_k >= 2 of a
   square-free M, let H_j be the subgroup of Z_M of order d_j and
   S = (H_1 u ... u H_k) minus {0}.  Then
       |det F_M[S, S]| = |c| * prod_j d_j^((d_j - 2)/2),
       c = (k - 1) * prod_j d_j - sum_j prod_{i != j} d_i.
   Checked here in floating point for a list of cases, and exactly in finite
   characteristic: for every prime r dividing c and not M, S is singular in
   every chart above r.

2. Order 105.  In both charts above 7, the vanishing principal minors of F_15
   are exactly the 15 translates of the punctured cross (H_3 u H_5) minus {0}
   and the 15 translates of its complement.

Usage: python3 structure_checks.py BIN_DIR
Standard library only.
"""
import cmath
import sys

sys.path.insert(0, ".")
import frontier as fr


def det_float(A):
    n = len(A)
    A = [r[:] for r in A]
    d = 1
    for c in range(n):
        piv = max(range(c, n), key=lambda r: abs(A[r][c]))
        if abs(A[piv][c]) < 1e-9:
            return 0
        if piv != c:
            A[c], A[piv] = A[piv], A[c]
            d = -d
        d *= A[c][c]
        for r in range(c + 1, n):
            f = A[r][c] / A[c][c]
            A[r] = [x - f * y for x, y in zip(A[r], A[c])]
    return d


def star(M, ds):
    S = set()
    for d in ds:
        S |= {(M // d) * i % M for i in range(1, d)}
    assert len(S) == sum(d - 1 for d in ds)
    return sorted(S)


def star_constant(ds):
    prod = 1
    for d in ds:
        prod *= d
    return (len(ds) - 1) * prod - sum(prod // d for d in ds)


def mask_of(S):
    return sum(1 << a for a in S)


def main(bin_dir):
    print("1a. punctured-star formula, floating point")
    for M, ds in [(6, [2, 3]), (10, [2, 5]), (14, [2, 7]), (15, [3, 5]), (21, [3, 7]), (35, [5, 7]),
                  (30, [2, 3, 5]), (30, [5, 6]), (30, [3, 10]), (30, [2, 15]), (42, [2, 3, 7]),
                  (42, [6, 7]), (42, [3, 14]), (70, [2, 5, 7]), (105, [3, 5, 7])]:
        S = star(M, ds)
        w = cmath.exp(2j * cmath.pi / M)
        val = abs(det_float([[w ** ((a * b) % M) for b in S] for a in S]))
        c = star_constant(ds)
        pred = abs(c)
        for d in ds:
            pred *= d ** ((d - 2) / 2)
        assert abs(val - pred) <= 1e-8 * max(1.0, pred), (M, ds, val, pred)
        print("   M=%d, orders %s: c = %d, |minor| = %.6g" % (M, ds, c, val))

    print("1b. punctured stars are singular in every chart above each prime r | c, r not dividing M")
    for M, ds in [(10, [2, 5]), (14, [2, 7]), (15, [3, 5]), (21, [3, 7]), (22, [2, 11])]:
        c = star_constant(ds)
        S = mask_of(star(M, ds))
        for r in sorted(fr.factorint(abs(c))):
            if M % r == 0:
                continue
            f, facs = fr.cyclotomic_factors(M, r)
            for h in facs:
                zeros = set(fr.run_chart(M, r, f, h, bin_dir, ["list"])["zero_masks"])
                assert S in zeros, (M, ds, r, h)
            print("   M=%d, orders %s, c = %d: singular in all %d chart(s) above %d" % (M, ds, c, len(facs), r))

    print("2. vanishing principal minors of F_15 above 7")
    M, r = 15, 7
    cross = star(M, [3, 5])
    full = (1 << M) - 1
    translates = {mask_of([(a + s) % M for a in cross]) for s in range(M)}
    expected = translates | {full ^ m for m in translates}
    assert len(expected) == 30
    f, facs = fr.cyclotomic_factors(M, r)
    for h in facs:
        zeros = set(fr.run_chart(M, r, f, h, bin_dir, ["list"])["zero_masks"])
        assert zeros == expected, h
        print("   chart h=%s: the 30 vanishing minors are the translates of %s and of its complement" % (h, cross))
    print("all checks passed")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(2)
    main(sys.argv[1])
