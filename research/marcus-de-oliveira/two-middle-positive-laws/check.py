#!/usr/bin/env python3
"""Uniform complex second-eigenvalue square with an unrestricted complex fourth.

All proof replay and evaluation use the Python standard library.
"""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from pathlib import Path
import argparse
import importlib.util
import json
import sys

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('archived_full_four_cycle', HERE.parent / 'full-four-cycle/check.py')
prior = importlib.util.module_from_spec(spec)
spec.loader.exec_module(prior)
sp, Q, Z, O, PERMS = prior.sp, prior.Q, prior.Z, prior.O, prior.PERMS
MONOMIALS = tuple((a, b) for a in range(4) for b in range(4-a))
KEYS = tuple((row, exponent) for row in range(9) for exponent in MONOMIALS)
RHS = tuple(F(row == 0 and exponent in ((0, 0), (2, 0), (0, 2))) for row, exponent in KEYS)


def load(name):
    return json.loads((HERE / name).read_text())


def exact(value):
    if isinstance(value, bool) or not isinstance(value, (int, F)):
        raise TypeError('Exact integer or Fraction required')
    return F(value)


def rational_matrix(value, rows, cols):
    if len(value) != rows or any(len(row) != cols for row in value) or any(type(x) is not str for row in value for x in row):
        raise ValueError('Incorrect exact matrix dimensions or scalar type')
    return [list(map(F, row)) for row in value]


def transpose(a):
    return list(map(list, zip(*a)))


def matmul(a, b):
    columns = transpose(b)
    return [[sum((x*y for x, y in zip(row, col) if x and y), F(0)) for col in columns] for row in a]


def matvec(a, v):
    return [sum((x*y for x, y in zip(row, v) if x and y), F(0)) for row in a]


def solve(a, rhs):
    """Gauss-Jordan with multiple right-hand sides, exact rational arithmetic."""
    n = len(a)
    if len(rhs) != n or any(len(row) != n for row in a):
        raise ValueError('Expected square coefficient matrix')
    rows = [list(map(F, row)) + list(map(F, target)) for row, target in zip(a, rhs)]
    for j in range(n):
        pivot = next((i for i in range(j, n) if rows[i][j]), None)
        if pivot is None:
            raise ArithmeticError('Singular coefficient matrix')
        rows[j], rows[pivot] = rows[pivot], rows[j]
        factor = rows[j][j]
        rows[j] = [x/factor for x in rows[j]]
        for i in range(n):
            if i != j and rows[i][j]:
                factor = rows[i][j]
                rows[i] = [x-factor*y for x, y in zip(rows[i], rows[j])]
    return [row[n:] for row in rows]


def inverse(a):
    return solve(a, [[F(i == j) for j in range(len(a))] for i in range(len(a))])


def square(vector):
    return prior.multiply(dict(zip(((0, 0), (1, 0), (0, 1)), vector)), dict(zip(((0, 0), (1, 0), (0, 1)), vector)))


@lru_cache(None)
def columns():
    generators = []
    for vector in product((-1, 0, 1), repeat=3):
        if any(vector) and next(x for x in vector if x) > 0:
            generators.append((vector, square(vector)))
    return tuple((j, vector, polynomial) for j in range(24) for vector, polynomial in generators
                 if PERMS[j][3] != 1 or vector[0]+vector[2] == 0)


def coefficient_matrix(eta):
    """Nine real rows, expanded in ten monomials of Re/Im a4: 90 equations."""
    affine = []
    for column in prior.gaps():
        entries = []
        for free in prior.FREE:
            constant = column[free] + eta * column[free | {1}]
            slope = column[free | {3}] + eta * column[free | {1, 3}]
            entries.append((constant, slope, Q(0, 1)*slope))
        affine.append(entries)
    output = []
    for row, (a, b) in KEYS:
        line = []
        for pairing, vector, polynomial in columns():
            if row == 0:
                value = F(polynomial.get((a, b), 0))
            else:
                part, index = ('real', row-1) if row <= 4 else ('imag', row-5)
                constant, real_slope, imag_slope = affine[pairing][index]
                value = (getattr(constant, part)*polynomial.get((a, b), 0)
                         + getattr(real_slope, part)*polynomial.get((a-1, b), 0)
                         + getattr(imag_slope, part)*polynomial.get((a, b-1), 0))
            line.append(F(value))
        output.append(line)
    return output


def difference(a, b):
    return [[x-y for x, y in zip(row, other)] for row, other in zip(a, b)]


def indices(value, length, maximum):
    if len(value) != length or len(set(value)) != length or any(type(i) is not int or not 0 <= i < maximum for i in value):
        raise ValueError('Invalid distinct matrix indices')
    return value


def verify(center, witness):
    if sp.e.multiply(sp.U, sp.e.adjoint(sp.U)) != sp.e.identity(4) or sp.e.multiply(sp.e.adjoint(sp.U), sp.U) != sp.e.identity(4):
        raise ArithmeticError('Archived fixture is not unitary')
    expected_columns = [{'pairing': j, 'linear': list(v)} for j, v, p in columns()]
    if center['columns'] != expected_columns:
        raise ValueError('Certificate columns differ from the declared square basis')
    selected = indices(center['selectedRows'], 60, 90)
    support = indices(center['support'], 60, len(columns()))
    if len(center['coefficients']) != len(columns()) or any(type(x) is not str for x in center['coefficients']):
        raise ValueError('Expected exact rational center coefficients')
    coefficients = list(map(F, center['coefficients']))
    if min(coefficients) < 0 or min(coefficients[j] for j in support) <= 0:
        raise ArithmeticError('Center coefficients lack required positivity')
    for j, (pairing, vector, poly) in enumerate(columns()):
        if PERMS[pairing][3] == 1 and coefficients[j] != 0:
            raise ArithmeticError('Forced pairings must have identically zero mass')
        if PERMS[pairing][3] != 1 and sum(abs(x) for x in vector) == 1 and coefficients[j] < F(1, 10000):
            raise ArithmeticError('The eighteen center Gram floors are missing')
    if any(PERMS[columns()[j][0]][3] == 1 for j in support):
        raise ArithmeticError('Perturbation basis changes a forced zero pairing')
    a0 = coefficient_matrix(Z)
    ar = difference(coefficient_matrix(O), a0)
    ai = difference(coefficient_matrix(Q(0, 1)), a0)
    if matvec(a0, coefficients) != list(RHS):
        raise ArithmeticError('Center polynomial coefficient identities fail')
    active = lambda a: [[row[j] for j in support] for row in a]
    a0b, arb, aib = map(active, (a0, ar, ai))
    b0 = [a0b[i] for i in selected]
    br, bi = [arb[i] for i in selected], [aib[i] for i in selected]
    inv = inverse(b0)
    xr, xi = matmul(inv, br), matmul(inv, bi)
    residual_r, residual_i = matvec(ar, coefficients), matvec(ai, coefficients)
    vr = matvec(inv, [-residual_r[i] for i in selected])
    vi = matvec(inv, [-residual_i[i] for i in selected])
    if matvec(a0b, vr) != [-x for x in residual_r] or matvec(a0b, vi) != [-x for x in residual_i]:
        raise ArithmeticError('The first-order full-system compatibility identities fail')
    w = rational_matrix(witness['basis'], 60, 26)
    wr = rational_matrix(witness['actionReal'], 26, 26)
    wi = rational_matrix(witness['actionImag'], 26, 26)
    sr = rational_matrix(witness['seedReal'], 26, 1)
    si = rational_matrix(witness['seedImag'], 26, 1)
    basis_rows = indices(witness['basisRows'], 26, 60)
    inverse([w[i] for i in basis_rows])  # The archived subspace has dimension exactly 26.
    if matmul(w, sr) != [[x] for x in vr] or matmul(w, si) != [[x] for x in vi]:
        raise ArithmeticError('Correction seeds are outside the archived subspace')
    xrw, xiw = matmul(xr, w), matmul(xi, w)
    if xrw != matmul(w, wr) or xiw != matmul(w, wi):
        raise ArithmeticError('Archived subspace is not invariant')
    if matmul(arb, w) != matmul(a0b, xrw) or matmul(aib, w) != matmul(a0b, xiw):
        raise ArithmeticError('Residual operators do not annihilate the subspace')
    y0 = [coefficients[j] for j in support]
    L = max(sum((abs(xr[i][j])+abs(xi[i][j]))*y0[j]/y0[i] for j in range(60)) for i in range(60))
    M = max((abs(vr[i])+abs(vi[i]))/y0[i] for i in range(60))
    radius = F(witness['radius'])
    rho, nu = radius*L, radius*M
    if (str(L), str(M)) != (witness['L'], witness['M']) or radius <= 0 or rho >= F(3, 4) or nu/(1-rho) >= F(1, 6):
        raise ArithmeticError('The uniform positivity bounds fail')
    return dict(center=center, witness=witness, a0=a0, ar=ar, ai=ai, a0b=a0b, arb=arb, aib=aib,
                b0=b0, br=br, bi=bi, selected=selected, support=support, coefficients=coefficients,
                y0=y0, vr=vr, vi=vi, radius=radius, L=L, M=M, rho=rho, nu=nu, relative=nu/(1-rho))


@lru_cache(None)
def prepared():
    return verify(load('center.json'), load('invariance.json'))

def coefficients_at(real_eta, imag_eta):
    real_eta, imag_eta = exact(real_eta), exact(imag_eta)
    data = prepared()
    if max(abs(real_eta), abs(imag_eta)) > data['radius']:
        raise ValueError('Second eigenvalue lies outside the certified closed square')
    matrix = [[x+real_eta*r+imag_eta*i for x, r, i in zip(row, rr, ii)]
              for row, rr, ii in zip(data['b0'], data['br'], data['bi'])]
    # y = y0 + (B0 + u Br + v Bi)^-1 B0 (u vr + v vi).
    seed = [real_eta*r+imag_eta*i for r, i in zip(data['vr'], data['vi'])]
    correction = solve(matrix, [[x] for x in matvec(data['b0'], seed)])
    result = data['coefficients'].copy()
    for j, y, delta in zip(data['support'], data['y0'], correction):
        result[j] = y+delta[0]
    return result


def weights_at(real_eta, imag_eta, real_fourth, imag_fourth):
    s, t = exact(real_fourth), exact(imag_fourth)
    coeff = coefficients_at(real_eta, imag_eta)
    masses = [F(0)]*24
    for value, (pairing, vector, poly) in zip(coeff, columns()):
        masses[pairing] += value*(vector[0]+vector[1]*s+vector[2]*t)**2/(1+s*s+t*t)
    return masses


def collision_check():
    delta, eta = Q(0, 1), Z
    a = [delta, eta, delta, delta]
    actual = sp.e.determinant(sp.matrix_sum(a))
    constant = (eta-delta)*sp.product(delta+b for b in sp.B_SPECTRUM[:3])
    if actual != Z or constant == Z:
        raise ArithmeticError('Collision certificate is degenerate')
    for p in PERMS:
        if sp.vertex(a, p) != (constant if p[3] == 1 else Z):
            raise ArithmeticError('Collision pairing classification fails')
    return {'forcedPairings': [j for j, p in enumerate(PERMS) if p[3] == 1], 'fourthEigenvalue': ['0', '1'], 'nonzeroVertexAtSecondZero': constant.encode()}


def report():
    data = prepared()
    collision = collision_check()
    samples = []
    r = data['radius']
    for u, v, s, t in ((0, 0, 0, 0), (r, r, 100, -100), (-r, r, 0, 1), (r, -r, F(7, 3), F(-11, 5)), (-r, -r, -3, 2)):
        masses = weights_at(u, v, s, t)
        if sum(masses) != 1 or any(m < F(1, 12000) for m, p in zip(masses, PERMS) if p[3] != 1) or any(m != 0 for m, p in zip(masses, PERMS) if p[3] == 1):
            raise ArithmeticError('Pointwise uniform mass bound fails')
        a = [Q(-7, 3), Q(u, v), Q(11, -5), Q(s, t)]
        determinant = prior.direct_check(a, masses)
        samples.append({'second': a[1].encode(), 'fourth': a[3].encode(), 'determinant': determinant.encode()})
    return {'schemaVersion': 1, 'family': 'a=(x,u+i*v,y,s+i*t); x,y complex and s,t real unrestricted',
            'secondSquareRadius': str(r), 'positivePermutations': 18, 'uniformMassLowerBound': '1/12000',
            'zeroPermutations': 6, 'invariantSubspaceDimension': 26, 'selectedVariables': 60,
            'realCoefficientEquations': 90, 'L': str(data['L']), 'M': str(data['M']),
            'rho': str(data['rho']), 'nu': str(data['nu']), 'relativeDisplacementBound': str(data['relative']),
            'collision': collision, 'directChecks': samples}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-report', action='store_true')
    args = parser.parse_args()
    result = report()
    encoded = json.dumps(result, indent=2)+'\n'
    if args.write_report:
        (HERE/'report.json').write_text(encoded)
    elif json.loads((HERE/'report.json').read_text()) != result:
        raise ArithmeticError('Archived report does not match exact replay')
    print(encoded, end='')


if __name__ == '__main__':
    main()
