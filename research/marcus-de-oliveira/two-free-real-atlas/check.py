"""Exact global real-line certificate. Python standard library, no optimizer.

The proof replay constructs the original determinant and all 24 vertices as
polynomials in three independent complex variables, then checks coefficients.
"""
from fractions import Fraction as F
from itertools import permutations
from pathlib import Path
import argparse
import json

HERE = Path(__file__).resolve().parent
PERMS = tuple(permutations(range(4)))
ZERO = (F(0), F(0))
ONE = (F(1), F(0))
U_NUMERATORS = ((-31, 492, -240, 300), (-300, 225, 500, 0),
                (240, -180, 225, 500), (-492, -256, -180, 225))
SPECTRUM = ((F(3), F(-2)), (F(2), F(0)), (F(1), F(2)), (F(0), F(-1)))


def add(a, b):
    return a[0]+b[0], a[1]+b[1]


def mul(a, b):
    return a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0]


def scale(q, a):
    return q*a[0], q*a[1]


def constant(a):
    return {} if a == ZERO else {(0, 0, 0): a}


def poly_add(a, b):
    out = dict(a)
    for monomial, value in b.items():
        out[monomial] = add(out.get(monomial, ZERO), value)
        if out[monomial] == ZERO:
            del out[monomial]
    return out


def poly_mul(a, b):
    out = {}
    for m, x in a.items():
        for n, y in b.items():
            key = tuple(i+j for i, j in zip(m, n))
            out[key] = add(out.get(key, ZERO), mul(x, y))
    return {key: value for key, value in out.items() if value != ZERO}


def poly_scale(q, a):
    return {key: scale(q, value) for key, value in a.items() if q != 0}


def poly_product(polys):
    out = constant(ONE)
    for poly in polys:
        out = poly_mul(out, poly)
    return out


def determinant(matrix):
    out = {}
    for perm in PERMS:
        inversions = sum(perm[i] > perm[j] for i in range(4) for j in range(i+1, 4))
        out = poly_add(out, poly_scale((-1)**inversions,
            poly_product(matrix[i][perm[i]] for i in range(4))))
    return out


def fixed_data():
    u = [[F(x, 625) for x in row] for row in U_NUMERATORS]
    for i in range(4):
        for j in range(4):
            if sum(u[i][k]*u[j][k] for k in range(4)) != (1 if i == j else 0):
                raise ValueError('Archived U is not orthogonal')
    b = [[ZERO for _ in range(4)] for _ in range(4)]
    for i in range(4):
        for j in range(4):
            for k in range(4):
                b[i][j] = add(b[i][j], scale(u[i][k]*u[j][k], SPECTRUM[k]))
    variables = [{(1, 0, 0): ONE}, {(0, 1, 0): ONE}, {(0, 0, 1): ONE}, {}]
    matrix = [[poly_add(constant(b[i][j]), variables[i] if i == j else {})
               for j in range(4)] for i in range(4)]
    vertices = [poly_product(poly_add(variables[p[j]], constant(SPECTRUM[j]))
                             for j in range(4)) for p in PERMS]
    return determinant(matrix), vertices


def load_certificate():
    return json.loads((HERE / 'certificate.json').read_text())


def decode(certificate):
    if certificate.get('fixedFourth') != ['0', '0'] or certificate.get('axis') != 1:
        raise ValueError('Wrong fixed family')
    terms = []
    seen = set()
    for item in certificate['terms']:
        perm, center, mass = tuple(item['permutation']), F(item['center']), F(item['mass'])
        if perm not in PERMS or mass <= 0 or (perm, center) in seen:
            raise ValueError('Invalid permutation, nonpositive mass, or duplicate term')
        seen.add((perm, center))
        terms.append((perm, center, mass))
    if not terms:
        raise ValueError('Empty certificate')
    return terms


def verify(certificate):
    terms = decode(certificate)
    target, vertices = fixed_data()
    numerators = {perm: {} for perm in PERMS}
    denominator = {}
    for perm, center, mass in terms:
        linear = poly_add({(0, 1, 0): ONE}, constant((-center, F(0))))
        square = poly_scale(mass, poly_mul(linear, linear))
        numerators[perm] = poly_add(numerators[perm], square)
        denominator = poly_add(denominator, square)
    a, b, c = [denominator.get((0, degree, 0), ZERO)[0] for degree in range(3)]
    if [str(a), str(b), str(c)] != certificate['denominator']:
        raise ValueError('Denominator does not equal the numerator sum')
    discriminant_gap = 4*a*c-b*b
    if c <= 0 or discriminant_gap <= 0:
        raise ValueError('Denominator is not certified positive on the projective line')
    residual = poly_scale(-1, poly_mul(denominator, target))
    for perm, vertex in zip(PERMS, vertices):
        residual = poly_add(residual, poly_mul(numerators[perm], vertex))
    if residual:
        raise ValueError('The full three-variable polynomial identity fails')
    return {'positiveSquareTerms': len(terms),
            'usedPermutations': sum(bool(p) for p in numerators.values()),
            'allPermutationsChecked': len(PERMS),
            'complexCoefficientSlots': 16,
            'denominatorAscending': [str(a), str(b), str(c)],
            'fourACMinusBSquared': str(discriminant_gap),
            'realLineCovered': True, 'projectiveInfinityCovered': True}


def exact_real(value):
    if isinstance(value, bool) or not isinstance(value, (int, F)):
        raise TypeError('Use an integer or Fraction, not a floating-point value')
    return F(value)


def weights_homogeneous(u, v, certificate=None):
    """Weights at t=u/v; (1,0) is the projective limit, not a finite eigenvalue."""
    u, v = exact_real(u), exact_real(v)
    if u == v == 0:
        raise ValueError('A projective coordinate pair cannot be zero')
    certificate = load_certificate() if certificate is None else certificate
    verify(certificate)
    raw = {perm: F(0) for perm in PERMS}
    for perm, center, mass in decode(certificate):
        raw[perm] += mass*(u-center*v)**2
    denominator = sum(raw.values())
    return [(perm, value/denominator) for perm, value in raw.items() if value]


def weights_at(t, certificate=None):
    return weights_homogeneous(t, 1, certificate)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--weights', type=F)
    parser.add_argument('--infinity', action='store_true')
    args = parser.parse_args()
    if args.weights is not None and args.infinity:
        parser.error('Choose a finite parameter or infinity')
    if args.weights is not None or args.infinity:
        pairs = weights_homogeneous(1, 0) if args.infinity else weights_at(args.weights)
        print(json.dumps([{'columnToRow': perm, 'mass': str(mass)} for perm, mass in pairs], indent=2))
    else:
        print(json.dumps(verify(load_certificate()), indent=2))
