"""Exact radical-field replay of the small-ratio zero-entry family in README.

No floating point: a number is a rational linear combination of square roots
of squarefree positive integers. This checker is a finite regression; the
one-parameter proof is the explicit orthonormal-frame calculation in README.
"""
import argparse
from fractions import Fraction as F
from functools import lru_cache
from math import gcd, isqrt
from pathlib import Path
import json
import platform
import time


@lru_cache(None)
def squarefree(n):
    if n == 0:
        return 0, 1
    r = isqrt(n)
    if r * r == n:
        return r, 1
    outside, inside, prime = 1, 1, 2
    while prime * prime <= n:
        power = 0
        while n % prime == 0:
            n //= prime
            power += 1
        outside *= prime ** (power // 2)
        if power % 2:
            inside *= prime
        prime += 1
    if n > 1:
        inside *= n
    return outside, inside


class Radical:
    def __init__(self, terms=0):
        if not isinstance(terms, dict):
            terms = {1: F(terms)}
        self.terms = {k: F(v) for k, v in terms.items() if v}

    def __add__(self, other):
        other = other if isinstance(other, Radical) else Radical(other)
        terms = self.terms.copy()
        for k, v in other.terms.items():
            terms[k] = terms.get(k, F(0)) + v
        return Radical(terms)

    __radd__ = __add__

    def __neg__(self):
        return Radical({k: -v for k, v in self.terms.items()})

    def __sub__(self, other):
        return self + (-other if isinstance(other, Radical) else -F(other))

    def __mul__(self, other):
        other = other if isinstance(other, Radical) else Radical(other)
        terms = {}
        for r, a in self.terms.items():
            for s, b in other.terms.items():
                common = gcd(r, s)
                key = (r // common) * (s // common)
                terms[key] = terms.get(key, F(0)) + a * b * common
        return Radical(terms)

    __rmul__ = __mul__

    def __eq__(self, other):
        other = other if isinstance(other, Radical) else Radical(other)
        return self.terms == other.terms

    def encode(self):
        return {str(k): str(v) for k, v in sorted(self.terms.items())}


def root(value):
    value = F(value)
    assert value >= 0
    a, r = squarefree(value.numerator)
    b, s = squarefree(value.denominator)
    return Radical({r: F(a)}) * Radical({s: F(b, value.denominator)})


def inner(a, b):
    return sum((x * y for x, y in zip(a, b)), Radical())


def combine(coefficients, basis):
    return [sum((coefficient * row[j] for coefficient, row in zip(coefficients, basis)), Radical())
            for j in range(4)]


def squared_minor(U, rows, columns):
    i, j = rows
    k, ell = columns
    determinant = U[i][k] * U[j][ell] - U[i][ell] * U[j][k]
    return determinant * determinant


def family(mu, nu):
    assert mu > 0 and mu < 1 and mu * mu + nu * nu == 1
    q = root(F(8, 9))
    invq = q * F(9, 8)
    d = root(9 + mu * mu)
    invd = d * (1 / (9 + mu * mu))
    basis = [[Radical(1), Radical(), Radical(), Radical()],
             [Radical(), Radical(), invq * F(2, 3), invq * F(-2, 3)],
             [Radical(), q, invq * F(-2, 9), invq * F(-2, 9)]]
    first = [Radical(0), Radical(F(1, 3)), Radical(F(2, 3)), Radical(F(2, 3))]
    a = combine([invd * (-mu * nu), invd * -3, invd * (mu * mu)], basis)
    b = combine([invd * (-3 * nu), invd * mu, invd * (3 * mu)], basis)
    c = combine([Radical(mu), Radical(), Radical(nu)], basis)
    return [first, a, b, c]


def run():
    start = time.process_time()
    records = []
    for denominator in [2, 5, 10, 100]:
        k = F(1, denominator)
        mu, nu = 2 * k / (1 + k * k), (1 - k * k) / (1 + k * k)
        U = family(mu, nu)
        for i in range(4):
            for j in range(4):
                assert inner(U[i], U[j]) == int(i == j)
        A = squared_minor(U, (0, 1), (0, 1))
        C = squared_minor(U, (0, 2), (0, 3))
        g = U[0][0] * U[0][0] + U[2][2] * U[2][2] + squared_minor(U, (0, 3), (2, 3)) - squared_minor(U, (0, 1), (0, 3))
        expected_A = mu * mu * nu * nu / (9 * (9 + mu * mu))
        expected_C = 4 * nu * nu / (9 + mu * mu)
        assert A == expected_A and C == expected_C and g == -4 * expected_A
        assert expected_A / expected_C == mu * mu / 36
        sqrtAC = 2 * mu * nu * nu / (3 * (9 + mu * mu))
        assert sqrtAC * sqrtAC == expected_A * expected_C
        margin = -3 * expected_A + 2 * sqrtAC
        assert margin == expected_A * (12 / mu - 3) and margin > 0
        records.append({'mu': str(mu), 'nu': str(nu), 'A': str(expected_A),
                        'C': str(expected_C), 'g0': str(-4 * expected_A),
                        't_star': str(mu / 6), 'star_margin': str(margin)})
    return {'cases': len(records), 'unitarity_and_identity_failures': 0,
            'arithmetic': 'exact rational combinations of squarefree radicals',
            'python': platform.python_version(), 'cpu_seconds': time.process_time() - start,
            'records': records}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-report', action='store_true')
    args = parser.parse_args()
    result = run()
    if args.write_report:
        Path(__file__).with_name('zero-stratum-report.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))
