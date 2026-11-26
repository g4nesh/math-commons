#!/usr/bin/env python3
"""Exact all-complex-fourth slice and singleton-marginal obstruction.

Python standard library only. Discovery used floating LP; this checker does not.
"""
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, permutations
from pathlib import Path
import argparse
import json
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'four-cycle-spectral'))
import spectral as sp
import fiber_certificate as archived
Q, Z, O = sp.Q, sp.Z, sp.O
PERMS = tuple(permutations(range(4)))
FREE = (frozenset(), frozenset((0,)), frozenset((2,)), frozenset((0, 2)))


def exact(value):
    if isinstance(value, bool) or not isinstance(value, (int, F)):
        raise TypeError('Expected an exact integer or Fraction')
    return F(value)


def add(poly, key, value):
    poly[key] = poly.get(key, 0) + value


def multiply(left, right):
    result = {}
    for (a, b), x in left.items():
        for (c, d), y in right.items():
            add(result, (a + c, b + d), x * y)
    return result


@lru_cache(None)
def gaps():
    """All 16 square-free coefficients of pairing minus determinant."""
    orbit = sp.matrix_sum([Z] * 4)
    columns = [dict() for _ in PERMS]
    for size in range(5):
        for indices in combinations(range(4), size):
            subset = frozenset(indices)
            outside = sorted(set(range(4)) - subset)
            target = sp.leibniz([[orbit[i][j] for j in outside] for i in outside])
            for column, p in zip(columns, PERMS):
                column[subset] = sp.product(b for b, i in zip(sp.B_SPECTRUM, p) if i not in subset) - target
    return columns


def load():
    return json.loads((HERE / 'certificate.json').read_text())


def numerators(data):
    fixed_second = data.get('fixedSecond', 0)
    if type(fixed_second) is not int or fixed_second not in (0, 1) or data.get('schemaVersion') != 1 or data.get('family') != f'a=(x,{fixed_second},y,s+i*t)':
        raise ValueError('Unknown certificate schema or family')
    if data.get('U') != [[v.encode() for v in row] for row in sp.U] or data.get('b') != [v.encode() for v in sp.B_SPECTRUM]:
        raise ValueError('Certificate fixture does not match archived U,b')
    out = [dict() for _ in PERMS]
    seen = set()
    for term in data['terms']:
        p, vector, text = term['pairing'], term['linear'], term['mass']
        if type(p) is not int or not 0 <= p < 24 or len(vector) != 3 or any(type(v) is not int or v not in (-1, 0, 1) for v in vector) or not any(vector):
            raise ValueError('Invalid pairing or linear square')
        if type(text) is not str:
            raise ValueError('Exact rational masses must be strings')
        mass = F(text)
        key = (p, tuple(vector))
        if mass <= 0 or key in seen:
            raise ValueError('Coefficient must be positive and square term unique')
        seen.add(key)
        linear = dict(zip(((0, 0), (1, 0), (0, 1)), vector))
        for exponent, coefficient in multiply(linear, linear).items():
            add(out[p], exponent, mass * coefficient)
    return out


def weights_at(s, t, data=None):
    s, t = exact(s), exact(t)
    polys = numerators(load() if data is None else data)
    denominator = 1 + s*s + t*t
    return tuple(sum(c * s**a * t**b for (a, b), c in poly.items()) / denominator for poly in polys)


def verify_plane(data):
    polys = numerators(data)
    if sp.e.multiply(sp.U, sp.e.adjoint(sp.U)) != sp.e.identity(4) or sp.e.multiply(sp.e.adjoint(sp.U), sp.U) != sp.e.identity(4):
        raise ArithmeticError('Archived U is not unitary')
    total = {}
    for poly in polys:
        for key, value in poly.items():
            add(total, key, value)
    if {k: v for k, v in total.items() if v} != {(0, 0): F(1), (2, 0): F(1), (0, 2): F(1)}:
        raise ArithmeticError('Polynomial masses do not sum to 1+s²+t²')
    for free in FREE:
        residual = {}
        for poly, column in zip(polys, gaps()):
            fixed_second = data.get('fixedSecond', 0)
            constant = column[free] + fixed_second * column[free | {1}]
            slope = column[free | {3}] + fixed_second * column[free | {1, 3}]
            gap = {(0, 0): constant, (1, 0): slope, (0, 1): Q(0, 1) * slope}
            for key, value in multiply(poly, gap).items():
                residual[key] = residual.get(key, Z) + value
        if any(value != Z for value in residual.values()):
            raise ArithmeticError('A complex polynomial gap coefficient is nonzero')
    return {
        'positiveSquareCoefficients': len(data['terms']),
        'supportedPermutations': sum(bool(any(poly.values())) for poly in polys),
        'minimumSquareCoefficient': str(min(F(term['mass']) for term in data['terms'])),
        'realCoefficientEquations': 90,
        'denominator': '1+s^2+t^2',
        'domain': 'all real s,t; all complex x,y',
        'fixedSecond': data.get('fixedSecond', 0),
        'strictMassClaim': False,
    }


def discrepancy(weights):
    return {key: sum((w * column[key] for w, column in zip(weights, gaps())), Z) for key in gaps()[0]}


def singleton_marginals(weights):
    return [[sum(w for w, p in zip(weights, PERMS) if p[j] == i) for j in range(4)] for i in range(4)]


def verify_repeated(data):
    texts = data['repeatedMiddleWeights']
    if len(texts) != 24 or any(type(v) is not str for v in texts):
        raise ValueError('Expected 24 exact rational repeated-middle masses')
    weights = tuple(map(F, texts))
    if min(weights) < 0 or sum(weights) != 1:
        raise ArithmeticError('Repeated-middle law is not a probability law')
    expected = [[v.real*v.real + v.imag*v.imag for v in row] for row in sp.U]
    if singleton_marginals(weights) != expected:
        raise ArithmeticError('Repeated-middle singleton marginals differ')
    delta = discrepancy(weights)
    q = -delta[frozenset((0, 3))]
    expected_delta = {frozenset(pair): sign*q for pair, sign in [((0, 1), 1), ((0, 3), -1), ((1, 2), -1), ((2, 3), 1)]}
    if q == Z or any(value != expected_delta.get(key, Z) for key, value in delta.items()):
        raise ArithmeticError('Repeated-middle discrepancy is not q(a1-a3)(a2-a4) with q nonzero')
    return {'nonzeroQ': q.encode(), 'positiveMasses': sum(w > 0 for w in weights), 'domain': 'a=(x,d,y,d), arbitrary complex x,y,d'}


def direct_check(a, weights):
    if min(weights) < 0 or sum(weights) != 1:
        raise ArithmeticError('Invalid pointwise masses')
    actual = sp.e.determinant(sp.matrix_sum(a))
    if actual != sp.leibniz(sp.matrix_sum(a)):
        raise ArithmeticError('Two determinant methods disagree')
    if actual != sum((w * sp.vertex(a, p) for w, p in zip(weights, PERMS)), Z):
        raise ArithmeticError('Direct determinant mixture mismatch')
    return actual


def report(data):
    positive = verify_plane(data)
    second_data = json.loads((HERE / 'certificate-one.json').read_text())
    second_plane = verify_plane(second_data)
    for s, t in ((0, 0), (0, 1), (100, -100), (F(7, 3), F(-11, 5))):
        direct_check([Q(-3, 2), O, Q(5, -7), Q(s, t)], weights_at(s, t, second_data))
    repeated = verify_repeated(data)
    previous = archived.check_certificate(json.loads((HERE.parent / 'four-cycle-spectral/fiber_certificate.json').read_text()))
    samples = []
    for s, t in ((0, 0), (1, 0), (0, 1), (F(7, 3), F(-11, 5)), (100, -100), (-2, -2)):
        weights = weights_at(s, t, data)
        a = [Q(F(-7, 3), 2), Z, Q(3, F(-5, 7)), Q(s, t)]
        samples.append({'s': str(s), 't': str(t), 'determinant': direct_check(a, weights).encode(), 'minimumMass': str(min(weights))})
    repeated_weights = tuple(map(F, data['repeatedMiddleWeights']))
    for d in (Z, Q(2, -3), Q(100, 100)):
        direct_check([Q(-3, 2), d, Q(7, -11), d], repeated_weights)
    if weights_at(0, 0, data) == weights_at(1, 1, data):
        raise ArithmeticError('Expected genuinely parameter-dependent weights')
    expected = [[v.real*v.real + v.imag*v.imag for v in row] for row in sp.U]
    if singleton_marginals(weights_at(1, 0, data)) == expected:
        raise ArithmeticError('Unexpected violation of singleton obstruction')
    return {'schemaVersion': 1, 'plane': positive, 'nonzeroSecondPlane': second_plane, 'repeatedMiddle': repeated, 'archivedAllFourConstantLawSeparator': previous, 'directChecks': samples}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-report', action='store_true', help='write deterministic report.json after all exact checks pass')
    args = parser.parse_args()
    output = report(load())
    encoded = json.dumps(output, indent=2) + '\n'
    if args.write_report:
        (HERE / 'report.json').write_text(encoded)
    else:
        if json.loads((HERE / 'report.json').read_text()) != output:
            raise ArithmeticError('Archived report differs from exact replay')
    print(encoded, end='')


if __name__ == '__main__':
    main()
