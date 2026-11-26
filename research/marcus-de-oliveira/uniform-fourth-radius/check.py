"""Exact uniform fourth-eigenvalue radius; standard library, no optimizer or float."""
from collections import defaultdict
from fractions import Fraction as F
from itertools import permutations
from math import factorial
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import resource

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


uniform = load('radius_uniform', ROOT / 'projective-fiber-stability/uniform/verify.py')
center = load('radius_center', ROOT / 'projective-fiber-stability/check.py')
sp = uniform.sp
PERMS = tuple(permutations(range(4)))
RADIUS = F(1, 10**100)
T = F(25)
S = F(1, 4096)
MU = uniform.EPSILON / 121


def poly_mul(a, b):
    result = defaultdict(F)
    for (i, j), x in a.items():
        for (k, l), y in b.items():
            result[i+k, j+l] += x*y
    return {key: value for key, value in result.items() if value}


def infinity_polynomial():
    result = {(0, 0): F(25000)}
    for factor in ({(1, 0): 18, (0, 1): -1, (0, 0): 4},
                   {(2, 0): 1, (1, 0): -4, (0, 2): 1, (0, 0): -1},
                   {(2, 0): 2, (1, 0): -1, (0, 2): 2, (0, 1): 1, (0, 0): -1}):
        result = poly_mul(result, factor)
    return result


def verify_infinity(rank):
    poly = infinity_polynomial()
    # The determinant has six affine rows and hence total degree at most six.
    # The 28-point triangular grid proves the complete polynomial identity.
    for x in range(7):
        for y in range(7-x):
            columns = uniform.vertex_rows(sp.Q(x, y), reciprocal=True)
            actual = uniform.determinant([columns[i] for i in rank['infinityBasis']])
            expected = sum(c*x**i*y**j for (i, j), c in poly.items())
            if actual != expected:
                raise ArithmeticError('Reciprocal determinant polynomial identity failed')
    sums = [sum(abs(c) for (i, j), c in poly.items() if i+j == degree) for degree in range(7)]
    if sums != [100000, 1075000, 3600000, 4500000, 9050000, 3800000, 0]:
        raise ArithmeticError('Unexpected reciprocal coefficient bounds')
    return sum(sums[d] / T**d for d in range(1, 7))


def augmented_patches(patches, corrections):
    result = []
    for patch, correction in zip(patches, corrections):
        terms = defaultdict(F)
        top = max(len(centers) for _, _, centers, _ in patch['terms'])
        max_center = max((abs(F(c)) for _, _, cs, _ in patch['terms'] for c in cs), default=F(0))
        if top not in (1, 2) or (1+max_center)**(2*top) > 121:
            raise ArithmeticError('Radial denominator bound exceeds121')
        for (j, word, centers, weight), amount in zip(patch['terms'], correction):
            mass = F(weight) - uniform.EPSILON*F(amount)
            if mass <= 0 or len(centers) > top:
                raise ArithmeticError('Base law exceeds its positivity/radial bounds')
            terms[j, word, tuple(centers)] += mass
        for j in range(patch['degree']+1):
            for p in PERMS:
                word = ''.join(str(i+1) for i in p)
                for centers in ((), ('0',)*top):
                    terms[j, word, centers] += uniform.EPSILON
        for j in range(patch['degree']+1):
            if sum(mass for (control, _, _), mass in terms.items() if control == j) != 1:
                raise ArithmeticError('Bernstein control coefficient mass is not one')
        result.append(dict(patch, terms=[[j, word, list(c), str(m)] for (j, word, c), m in terms.items()]))
    return result


def bounds(radius=RADIUS):
    if not isinstance(radius, F) or radius <= 0 or 2*radius > 1:
        raise ValueError('Use a positive exact radius with complex displacement at most one')
    eta = S/1000
    gamma = eta**5
    # Near one of Y=-2,0,1,2, |P'(Y)| <=663<700 for |Y|<=3.
    if 6*S-700*eta < gamma or S+radius > F(1, 1024):
        raise ArithmeticError('Pole/separation bound failed')
    determinant_lower = gamma/T**2
    inverse_finite = 9*factorial(8)*(100*T)**8/determinant_lower
    inverse_infinity = 9*factorial(8)*F(1000)**8/50000
    variation_finite, variation_infinity = 200*T*radius, 2000*radius
    result = {'finiteDeterminantLower': determinant_lower, 'finiteInverseRowBound': inverse_finite,
              'infinityInverseRowBound': inverse_infinity, 'baseMassFloor': MU, 'massFloorAfterMotion': MU/2}
    for name, inverse, variation in (('finite', inverse_finite, variation_finite),
                                     ('infinity', inverse_infinity, variation_infinity)):
        neumann = 9*inverse*variation
        mass_change = 50*inverse*variation
        if neumann >= F(1, 2) or mass_change >= MU/2:
            raise ArithmeticError('Radius fails inverse/positivity bounds: '+name)
        result[name+'NeumannBound'] = neumann
        result[name+'MassChangeRelativeToFloor'] = mass_change/MU
    return result


def signed_masses():
    determinant = sp.leibniz(sp.U)
    if determinant.imag != 0 or abs(determinant.real) != 1:
        raise ArithmeticError('Fixture is not real orthogonal with determinant +/-1')
    result = []
    for p in PERMS:  # column-to-row convention
        sign = (-1)**sum(p[i] > p[j] for i in range(4) for j in range(i+1, 4))
        mass = determinant*sign*sp.product(sp.U[p[j]][j] for j in range(4))
        if mass.imag:
            raise ArithmeticError('The fixed signed law must be real')
        result.append(mass.real)
    if sum(result) != 1 or sum(map(abs, result)) > 24:
        raise ArithmeticError('Signed mass bound failed')
    return result


def raw_rows(w, delta, reciprocal=False):
    b = [z+delta for z in sp.B_SPECTRUM]
    columns = []
    for i, j, l, k in PERMS:  # row-to-column convention
        middle = sp.O+w*b[j] if reciprocal else w+b[j]
        values = [b[i]*b[l]*b[k], middle*b[l]*b[k], middle*b[i]*b[k], middle*b[k]]
        columns.append([F(1)]+[getattr(v, part) for part in ('real', 'imag') for v in values])
    return list(map(list, zip(*columns)))


def real_det(matrix):
    a = [row.copy() for row in matrix]
    result = F(1)
    for j in range(len(a)):
        k = next((k for k in range(j, len(a)) if a[k][j]), None)
        if k is None:
            return F(0)
        if j != k:
            a[j], a[k] = a[k], a[j]
            result = -result
        pivot = a[j][j]
        result *= pivot
        for k in range(j+1, len(a)):
            ratio = a[k][j]/pivot
            for l in range(j+1, len(a)):
                a[k][l] -= ratio*a[j][l]
    return result


def dot(a, b):
    return sum(x*y for x, y in zip(a, b))


def regular_weights(w, delta, patches, rank, signed, reciprocal=False):
    law = uniform.cover.weights_at(w.real, w.imag, patches)
    inverses = [tuple(q.index(j) for j in range(4)) for q in PERMS]
    base = [law[''.join(str(i+1) for i in p)] for p in inverses]
    signs = [signed[PERMS.index(p)] for p in inverses]
    coordinate = sp.O/w if reciprocal else w
    original, moved = raw_rows(coordinate, sp.Z, reciprocal), raw_rows(coordinate, delta, reciprocal)
    if any(dot(row, base) != dot(row, signs) for row in original):
        raise ArithmeticError('Base law does not solve normalized raw equations')
    choices = [rank['infinityBasis']] if reciprocal else rank['bases']
    basis = max(choices, key=lambda b: abs(real_det([[row[j] for j in b] for row in original])))
    residual = [dot([a-b for a, b in zip(new, old)], [m-s for m, s in zip(base, signs)])
                for old, new in zip(original, moved)]
    change = sp.exact_solve([[row[j] for j in basis] for row in moved], [-r for r in residual])
    for j, amount in zip(basis, change):
        base[j] += amount
    if min(base) <= MU/2 or any(dot(row, base) != dot(row, signs) for row in moved):
        raise ArithmeticError('Moved positive law failed exact coefficient equations')
    return [base[inverses.index(p)] for p in PERMS]


def regression(patches, rank, signed, centers):
    tests = [('zero', sp.Z, sp.Q(RADIUS, -RADIUS), False),
             ('ordinary', sp.Q(7, 5), sp.Q(-RADIUS, RADIUS), False),
             ('reciprocal', sp.Q(T, T), sp.Q(-RADIUS, -RADIUS), True)]
    records = []
    for name, w, delta, reciprocal in tests:
        masses = regular_weights(w, delta, patches, rank, signed, reciprocal)
        a = [sp.Q(2, 3), w+delta, sp.Q(-1, 2), delta]
        value = sum((m*sp.vertex(a, p) for m, p in zip(masses, PERMS)), sp.Z)
        if value != sp.leibniz(sp.matrix_sum(a)):
            raise ArithmeticError('Dense determinant regression failed')
        records.append({'case': name, 'support': len(masses), 'massAboveHalfFloor': min(masses) > MU/2,
                        'weightsSha256': hashlib.sha256('|'.join(map(str, masses)).encode()).hexdigest()})
    delta = sp.Q(RADIUS, RADIUS)
    for j, pole in enumerate(center.POLES):
        eta = sp.Q(S, -S)+delta
        masses = center.weights_at(j, delta, eta, centers)
        a = [sp.Q(2, 3), pole+eta, sp.Q(-1, 2), delta]
        if min(masses) <= MU/2 or sum((m*sp.vertex(a, p) for m, p in zip(masses, PERMS)), sp.Z) != sp.leibniz(sp.matrix_sum(a)):
            raise ArithmeticError('Collision boundary regression failed')
        records.append({'case': 'pole-'+str(j+1), 'support': len(masses), 'massAboveHalfFloor': True,
                        'weightsSha256': hashlib.sha256('|'.join(map(str, masses)).encode()).hexdigest()})
    return records


def verify():
    path = ROOT / 'projective-fiber-stability'
    patches = json.loads((ROOT / 'complex-direction-cover/patches.json').read_text())
    corrections = json.loads((path / 'uniform/corrections.json').read_text())
    rank = json.loads((path / 'uniform/rank.json').read_text())
    centers = json.loads((path / 'bases.json').read_text())
    atlas = uniform.cover.atlas
    if [[sp.Q(F(v, 625)) for v in row] for row in atlas.U_NUMERATORS] != sp.U or [sp.Q(*z) for z in atlas.SPECTRUM] != sp.B_SPECTRUM:
        raise ArithmeticError('Archived proof fixtures differ')
    if any(z.real*z.real+z.imag*z.imag > 16 for z in sp.B_SPECTRUM):
        raise ArithmeticError('Opposing spectrum exceeds modulus bound4')
    uniform.strict_patches(patches, corrections)
    uniform.rank_identities(rank)
    center.verify(centers)
    augmented = augmented_patches(patches, corrections)
    multiplier_bounds = [sum(abs(F(v)) for _, _, _, v in terms) for terms in rank['multipliers']]
    if max(multiplier_bounds) > 1:
        raise ArithmeticError('Bezout coefficient sums exceed one')
    loss = verify_infinity(rank)
    if loss >= 50000:
        raise ArithmeticError('Infinity determinant lower bound failed')
    signed = signed_masses()
    inequalities = bounds()
    records = regression(augmented, rank, signed, centers)
    return {'radiusPowerOfTen': -100, 'radius': str(RADIUS), 'domain': 'max(abs(Re delta),abs(Im delta)) <= radius; fixed archived real U,b; all x,z,y complex',
            'T': str(T), 'poleHalfwidth': str(S), 'bezoutCoefficientL1': list(map(str, multiplier_bounds)),
            'infinityDeterminantLoss': str(loss), 'signedTotalVariation': str(sum(map(abs, signed))),
            'bounds': {key: str(value) for key, value in inequalities.items()}, 'exactRegressions': records,
            'scope': 'Exact rational replay plus written continuous bounds; not Lean, not universal MDO, independent-owner review pending'}


if __name__ == '__main__':
    resource.setrlimit(resource.RLIMIT_CPU, (60, 60))
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--check', type=Path)
    args = parser.parse_args()
    report = verify()
    if args.check and json.loads(args.check.read_text()) != report:
        raise SystemExit('Archived certificate differs from exact replay')
    if args.output:
        args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({'radiusPowerOfTen': report['radiusPowerOfTen'], 'exactRegressions': len(report['exactRegressions']), 'verified': True}))
