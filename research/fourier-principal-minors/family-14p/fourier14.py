"""Check the norm and flag criteria for the order-14 Fourier matrix.

Uses only Python 3's standard library. See the accompanying research note for
why affine symmetries and complementary minors cover every subset.
"""

from collections import defaultdict
from itertools import combinations
from math import gcd


ORDER = 14
DEGREE = 6


def determinant_polynomial(indices):
    """Leibniz determinant by column subsets in Z[x]/(x^14 - 1)."""
    count = len(indices)
    values = [[0] * ORDER for _ in range(1 << count)]
    values[0][0] = 1
    for mask, polynomial in enumerate(values[:-1]):
        row = indices[bin(mask).count('1')]
        for col, index in enumerate(indices):
            bit = 1 << col
            if mask & bit:
                continue
            inversions = bin(mask >> (col + 1)).count('1')
            sign = (-1) ** inversions
            shift = row * index
            result = values[mask | bit]
            for power, coefficient in enumerate(polynomial):
                result[(power + shift) % ORDER] += sign * coefficient
    return values[-1]


def cyclotomic_coordinates(polynomial):
    """Reduce with x^7=-1 and Phi_14(x)=x^6-x^5+x^4-x^3+x^2-x+1."""
    low = [polynomial[i] - polynomial[i + 7] for i in range(7)]
    return [low[i] - (-1) ** i * low[6] for i in range(DEGREE)]


def integer_determinant(matrix):
    """Fraction-free elimination, including exact division and row exchanges."""
    rows = [row.copy() for row in matrix]
    last = len(rows) - 1
    previous, sign = 1, 1
    for col in range(last):
        if rows[col][col] == 0:
            pivot = next((k for k in range(col + 1, last + 1) if rows[k][col]), None)
            if pivot is None:
                return 0
            rows[col], rows[pivot] = rows[pivot], rows[col]
            sign *= -1
        for i in range(col + 1, last + 1):
            for j in range(col + 1, last + 1):
                value = rows[i][j] * rows[col][col] - rows[i][col] * rows[col][j]
                rows[i][j], remainder = divmod(value, previous)
                if remainder:
                    raise ArithmeticError('Bareiss division was not exact')
        previous = rows[col][col]
    return sign * rows[last][last]


def norm(polynomial):
    # Column j describes multiplication of x^j by the determinant polynomial.
    columns = []
    for shift in range(DEGREE):
        shifted = polynomial[-shift:] + polynomial[:-shift] if shift else polynomial
        columns.append(cyclotomic_coordinates(shifted))
    return abs(integer_determinant([list(row) for row in zip(*columns)]))


def prime_factors(value):
    if value == 0:
        raise ValueError('Zero norm has no finite list of prime divisors')
    factors = set()
    divisor = 2
    while divisor * divisor <= value:
        if value % divisor == 0:
            factors.add(divisor)
            while value % divisor == 0:
                value //= divisor
        divisor = 3 if divisor == 2 else divisor + 2
    if value > 1:
        factors.add(value)
    return factors


def remainder(polynomial, divisor, prime):
    """Remainder modulo a monic polynomial over F_prime, low coefficients first."""
    work = list(polynomial)
    degree = len(divisor) - 1
    for k in range(len(work) - 1, degree - 1, -1):
        leading = work[k] % prime
        for j in range(degree):
            work[k - degree + j] -= leading * divisor[j]
    return tuple(v % prime for v in work[:degree])


def prime_charts(prime):
    """The only primes needing the flag test here split in degrees one or two."""
    cyclotomic = [1, -1, 1, -1, 1, -1, 1]
    linear = [(-r, 1) for r in range(prime) if remainder(cyclotomic, (-r, 1), prime) == (0,)]
    if len(linear) == DEGREE:
        return linear
    # If p=-1 mod 14, Frobenius pairs each primitive root with its inverse.
    if prime % ORDER != ORDER - 1 or linear:
        raise ValueError(f'Unsupported splitting type at {prime}')
    quadratics = [(1, -t, 1) for t in range(prime)
                  if remainder(cyclotomic, (1, -t, 1), prime) == (0, 0)]
    assert len(quadratics) * 2 == DEGREE
    return quadratics


def singular_charts(polynomial, charts, prime):
    return sum(1 << j for j, divisor in enumerate(charts)
               if not any(remainder(polynomial, divisor, prime)))


def chain_covers_charts(prime, candidate_subsets):
    charts = prime_charts(prime)
    full_set = (1 << ORDER) - 1
    full_cover = (1 << len(charts)) - 1
    bad_charts = [0] * (full_set + 1)
    for indices in candidate_subsets:
        subset = sum(1 << i for i in indices)
        polynomial = determinant_polynomial(indices)
        conjugate = polynomial[:1] + polynomial[:0:-1]
        for mask, poly in ((subset, polynomial), (full_set ^ subset, conjugate)):
            bad = singular_charts(poly, charts, prime)
            assert bad and bad_charts[mask] in (0, bad)
            bad_charts[mask] = bad

    # Track all unions of bad charts attained by nested subsets ending at S.
    reachable = [set() for _ in bad_charts]
    reachable[0].add(0)
    for subset in range(1, full_set + 1):
        bits = subset
        while bits:
            bit = bits & -bits
            bits ^= bit
            reachable[subset].update(prev | bad_charts[subset]
                                     for prev in reachable[subset ^ bit])
        if full_cover in reachable[subset]:
            return True
    return False


def main():
    units = [a for a in range(ORDER) if gcd(a, ORDER) == 1]
    factors = set()
    valuations = defaultdict(lambda: [0] * (ORDER + 1))
    divisibility = defaultdict(set)
    for size in range(ORDER // 2 + 1):
        seen = set()
        representatives = 0
        by_size = set()
        for indices in combinations(range(ORDER), size):
            if indices in seen:
                continue
            orbit = {tuple(sorted((a * i + b) % ORDER for i in indices))
                     for a in units for b in range(ORDER)}
            seen.update(orbit)
            representatives += 1
            absolute_norm = norm(determinant_polynomial(indices))
            primes = prime_factors(absolute_norm)
            by_size.update(primes)
            for prime in primes - {2, 7}:
                divisibility[prime].update(orbit)
                value, exponent = absolute_norm, 0
                while value % prime == 0:
                    value //= prime
                    exponent += 1
                for rank in (size, ORDER - size):
                    valuations[prime][rank] = max(valuations[prime][rank], exponent)
        factors.update(by_size)
        print(f'Size {size}: {len(seen)} subsets, {representatives} orbits, primes {sorted(by_size)}')
    # Jacobi complementary minors can change divisibility only at primes of 14.
    # Both primes actually occur, for instance in the full Fourier determinant.
    print('Prime support for all sizes:', sorted(factors | prime_factors(ORDER)))
    unresolved = []
    for prime, bounds in sorted(valuations.items()):
        if prime in (3, 5):  # P(42) and P(70) have prior proofs.
            continue
        budget = sum(sorted(bounds, reverse=True)[:prime])
        if budget >= DEGREE:
            covered = chain_covers_charts(prime, divisibility[prime])
            print(f'p={prime}, budget={budget}, chain covering every chart: {covered}')
            if covered:
                unresolved.append(prime)
    print('Square-free multipliers remaining after lifting and prior results:', unresolved)


if __name__ == '__main__':
    main()
