#!/usr/bin/env python3
"""Bounded numerical MDO baseline. Standard library only; no proof claim.

The planar hull uses exact integer orientation of the stored binary64 values.
The input unitary, matrix products and determinants remain approximate.
"""
from __future__ import annotations
import argparse
import cmath
from collections import Counter
from itertools import permutations
import json
import math
from pathlib import Path
import platform
import random
import sys
import time

EPS = sys.float_info.epsilon
MAX_TRIALS = 200
MAX_DIMENSION = 6
CPU_LIMIT_SECONDS = 60.0
FAMILIES = ('generic', 'near_monomial', 'clustered', 'repeated', 'structured_control')


def dot(left, right):
    terms = [a.conjugate() * b for a, b in zip(left, right)]
    return complex(math.fsum(z.real for z in terms), math.fsum(z.imag for z in terms))


def adjoint(matrix):
    return [[matrix[j][i].conjugate() for j in range(len(matrix))] for i in range(len(matrix))]


def multiply(left, right):
    n = len(left)
    return [[sum(left[i][k] * right[k][j] for k in range(n)) for j in range(n)] for i in range(n)]


def identity(n):
    return [[complex(i == j) for j in range(n)] for i in range(n)]


def frobenius(matrix):
    return math.sqrt(math.fsum(abs(z) ** 2 for row in matrix for z in row))


def difference(left, right):
    return [[a - b for a, b in zip(x, y)] for x, y in zip(left, right)]


def complex_qr(rng, n):
    """Complex Gaussian columns, two-pass modified Gram--Schmidt QR.

    This is a numerical unitary construction, not an exact certificate or a
    claim that the finite generated sample exhausts the unitary group.
    """
    columns = []
    for _ in range(n):
        column = [complex(rng.gauss(0, 1), rng.gauss(0, 1)) for _ in range(n)]
        for _ in range(2):
            for old in columns:
                coefficient = dot(old, column)
                column = [value - coefficient * basis for value, basis in zip(column, old)]
        norm = math.sqrt(max(0.0, dot(column, column).real))
        if norm <= 1e-12:
            raise ValueError('Numerically rank-deficient QR sample; do not silently repair it.')
        columns.append([value / norm for value in column])
    return [[columns[j][i] for j in range(n)] for i in range(n)]


def determinant_lu(matrix):
    """Partial pivoting; no tolerance-based declaration of singularity."""
    work = [row[:] for row in matrix]
    determinant = 1 + 0j
    for column in range(len(work)):
        pivot = max(range(column, len(work)), key=lambda row: abs(work[row][column]))
        if work[pivot][column] == 0:
            return 0j
        if pivot != column:
            work[pivot], work[column] = work[column], work[pivot]
            determinant = -determinant
        diagonal = work[column][column]
        determinant *= diagonal
        for row in range(column + 1, len(work)):
            factor = work[row][column] / diagonal
            for index in range(column + 1, len(work)):
                work[row][index] -= factor * work[column][index]
            work[row][column] = 0j
    return determinant


def determinant_permutations(matrix):
    """Independent Leibniz determinant with fsum real/imag accumulation."""
    values = []
    n = len(matrix)
    for permutation in permutations(range(n)):
        sign = -1 if sum(permutation[i] > permutation[j] for i in range(n) for j in range(i + 1, n)) % 2 else 1
        value = complex(sign)
        for i, j in enumerate(permutation):
            value *= matrix[i][j]
        values.append(value)
    return complex(math.fsum(value.real for value in values), math.fsum(value.imag for value in values))


def permutation_points(a, b):
    result = []
    for permutation in permutations(range(len(a))):
        product = 1 + 0j
        for i, j in enumerate(permutation):
            product *= a[i] + b[j]
        result.append((permutation, product))
    return result


def cross(origin, left, right):
    return (left[0] - origin[0]) * (right[1] - origin[1]) - (left[1] - origin[1]) * (right[0] - origin[0])


def integer_coordinates(points):
    """One exact power-of-two scale for all stored finite float coordinates."""
    ratios = [[float(z.real).as_integer_ratio(), float(z.imag).as_integer_ratio()] for z in points]
    scale = max(denominator for pair in ratios for _, denominator in pair)
    return [tuple(numerator * (scale // denominator) for numerator, denominator in pair) for pair in ratios], scale


def hull_indices(points):
    """Monotone chain with exact integer cross products, no geometric epsilon."""
    unique = sorted({point: index for index, point in enumerate(points)}.items())
    if len(unique) <= 1:
        return [unique[0][1]] if unique else []
    lower, upper = [], []
    for point in unique:
        while len(lower) >= 2 and cross(lower[-2][0], lower[-1][0], point[0]) <= 0:
            lower.pop()
        lower.append(point)
    for point in reversed(unique):
        while len(upper) >= 2 and cross(upper[-2][0], upper[-1][0], point[0]) <= 0:
            upper.pop()
        upper.append(point)
    return [index for _, index in lower[:-1] + upper[:-1]]


def hull_gap(vertices, query):
    """Signed support gap for polygon; distance for point/segment.

    All edge orientations are exact for stored floats. The returned distances
    and separating normal are floating point; neither certifies the underlying
    mathematical input. Degenerate hulls are handled in their true dimension.
    """
    if not vertices:
        raise ValueError('A hull needs at least one point.')
    if not all(math.isfinite(z.real) and math.isfinite(z.imag) for z in [*vertices, query]):
        raise ValueError('Hull values must be finite.')
    integers, scale = integer_coordinates([*vertices, query])
    hull = hull_indices(integers[:-1])
    q = integers[-1]
    if len(hull) == 1:
        p = integers[hull[0]]
        gap = math.hypot((q[0] - p[0]) / scale, (q[1] - p[1]) / scale)
        return {'dimension': 0, 'vertex_indices': hull, 'signed_support_gap': gap, 'outside_stored_float_hull': q != p, 'normal': None}
    if len(hull) == 2:
        p, r = (integers[index] for index in hull)
        dx, dy = r[0] - p[0], r[1] - p[1]
        squared = dx * dx + dy * dy
        projection = (q[0] - p[0]) * dx + (q[1] - p[1]) * dy
        orientation = cross(p, r, q)
        if projection < 0:
            gap = math.hypot((q[0] - p[0]) / scale, (q[1] - p[1]) / scale)
        elif projection > squared:
            gap = math.hypot((q[0] - r[0]) / scale, (q[1] - r[1]) / scale)
        else:
            gap = abs(orientation) / math.hypot(dx, dy) / scale
        return {'dimension': 1, 'vertex_indices': hull, 'signed_support_gap': gap, 'outside_stored_float_hull': orientation != 0 or projection < 0 or projection > squared, 'normal': None}
    strongest = None
    outside = False
    for offset, index in enumerate(hull):
        successor = hull[(offset + 1) % len(hull)]
        p, r = integers[index], integers[successor]
        dx, dy = r[0] - p[0], r[1] - p[1]
        orientation = cross(p, r, q)
        length = math.hypot(dx, dy)
        gap = -orientation / length / scale
        outside = outside or orientation < 0
        if strongest is None or gap > strongest['signed_support_gap']:
            strongest = {'signed_support_gap': gap, 'normal': [dy / length, -dx / length], 'edge_vertex_indices': [index, successor]}
    return {'dimension': 2, 'vertex_indices': hull, 'outside_stored_float_hull': outside, **strongest}


def finite_complex(value):
    return isinstance(value, (int, float, complex)) and math.isfinite(complex(value).real) and math.isfinite(complex(value).imag)


def validate_sample(a, b, u):
    n = len(a)
    if not 1 <= n <= MAX_DIMENSION or len(b) != n or len(u) != n or any(len(row) != n for row in u):
        raise ValueError('Mismatched dimensions or dimension outside 1..6.')
    if not all(finite_complex(z) for z in [*a, *b, *(z for row in u for z in row)]):
        raise ValueError('Sample values must be finite complex scalars.')
    residual = frobenius(difference(multiply(adjoint(u), u), identity(n)))
    if residual > 1e-10:
        raise ValueError(f'Not approximately unitary: residual {residual}.')
    return residual


def analyze_case(a, b, u):
    n = len(a)
    residual = validate_sample(a, b, u)
    b_matrix = [[sum(u[i][k] * b[k] * u[j][k].conjugate() for k in range(n)) for j in range(n)] for i in range(n)]
    total = [[b_matrix[i][j] + (a[i] if i == j else 0) for j in range(n)] for i in range(n)]
    det_lu = determinant_lu(total)
    det_check = determinant_permutations(total)
    points = permutation_points(a, b)
    geometric = hull_gap([z for _, z in points], det_lu)
    point_scale = max(1.0, abs(det_lu), *(abs(z) for _, z in points))
    matrix_scale = max(1.0, max(map(abs, a)) + max(map(abs, b)))
    disagreement = abs(det_lu - det_check)
    # Deliberately heuristic screening allowance, NOT a rigorous interval bound.
    allowance = (256 * EPS * n ** 3 + 8 * residual) * matrix_scale ** n + 8 * disagreement
    normal_residual = frobenius(difference(multiply(b_matrix, adjoint(b_matrix)), multiply(adjoint(b_matrix), b_matrix))) / max(1.0, frobenius(b_matrix) ** 2)
    prescribed = [[u[i][j] * b[j] for j in range(n)] for i in range(n)]
    spectral_residual = frobenius(difference(multiply(b_matrix, u), prescribed)) / max(1.0, frobenius(b_matrix))
    gap = geometric['signed_support_gap']
    if gap > allowance:
        classification = 'unverified_outside_candidate'
    elif gap >= -allowance:
        classification = 'numerically_near_boundary'
    else:
        classification = 'numerically_inside'
    return {
        'determinant_lu': encode(det_lu), 'determinant_leibniz': encode(det_check),
        'determinant_disagreement': disagreement, 'unitarity_frobenius_residual': residual,
        'normality_relative_frobenius_residual': normal_residual,
        'spectral_equation_relative_frobenius_residual': spectral_residual,
        'max_U_imaginary_component': max(abs(z.imag) for row in u for z in row),
        'permutation_count': len(points), 'point_scale': point_scale,
        'heuristic_roundoff_allowance': allowance, 'heuristic_relative_allowance': allowance / point_scale,
        'relative_signed_support_gap': gap / point_scale,
        'classification': classification, 'hull': geometric,
        'hull_vertices': [{'permutation': list(points[index][0]), 'point': encode(points[index][1])} for index in geometric['vertex_indices']],
    }


def encode(value):
    return [float(value.real), float(value.imag)]


def complex_uniform(rng, radius=2):
    return complex(rng.uniform(-radius, radius), rng.uniform(-radius, radius))


def generate_case(n, family, index, seed):
    rng = random.Random(seed)
    a = [complex_uniform(rng) for _ in range(n)]
    b = [complex_uniform(rng) for _ in range(n)]
    u = complex_qr(rng, n)
    detail = 'complex Gaussian QR mixing and independent uniform complex spectra'
    if family == 'clustered':
        ca, cb = complex_uniform(rng), complex_uniform(rng)
        spread = (1e-2, 1e-4, 1e-6, 1e-8)[index % 4]
        a = [ca + spread * complex_uniform(rng, 1) for _ in range(n)]
        b = [cb + spread * complex_uniform(rng, 1) for _ in range(n)]
        detail = f'complex clusters with spread {spread!r}; generic QR mixing'
    elif family == 'repeated':
        pools = [[complex(rng.randint(-4, 4) / 2, rng.randint(-4, 4) / 2) for _ in range(2)] for _ in range(2)]
        a = [pools[0][i % 2] for i in range(n)]
        b = [pools[1][i % 2] for i in range(n)]
        detail = 'two-value repeated complex spectra on half-integer grid; generic QR mixing'
    elif family == 'near_monomial':
        vector = [complex(rng.gauss(0, 1), rng.gauss(0, 1)) for _ in range(n)]
        norm_squared = dot(vector, vector).real
        angle = (1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9)[index]
        coefficient = cmath.exp(1j * angle) - 1
        perturbation = [[complex(i == j) + coefficient * vector[i] * vector[j].conjugate() / norm_squared for j in range(n)] for i in range(n)]
        permutation = list(range(n)); rng.shuffle(permutation)
        phases = [cmath.exp(1j * rng.uniform(-math.pi, math.pi)) for _ in range(n)]
        u = [[phases[i] * perturbation[permutation[i]][j] for j in range(n)] for i in range(n)]
        detail = f'complex monomial times unitary rank-one phase perturbation, angle {angle!r}'
    elif family == 'structured_control':
        if index < 3:
            a = [complex(1, 0.5)] * n
            b = [complex(rng.randint(-3, 3), rng.randint(-3, 3)) for _ in range(n)]
            detail = 'point-hull control: scalar A; binary-exact bounded spectra'
        elif index < 6:
            a = [complex(rng.randint(-5, 5), 0) for _ in range(n)]
            b = [complex(rng.randint(-5, 5), 0) for _ in range(n)]
            detail = 'segment-hull control: real spectra, complex QR; Hermitian pair'
        else:
            a = [0j] * n
            b = [0j, *(complex(rng.randint(-3, 3), rng.randint(-3, 3)) for _ in range(n - 1))]
            detail = 'singular point-hull control: A=0 and one zero b eigenvalue'
    return a, b, u, detail


def summarize(cases, errors):
    valid = [case for case in cases if 'analysis' in case]
    analyses = [case['analysis'] for case in valid]
    def group(predicate):
        selected = [case for case in valid if predicate(case)]
        return {'count': len(selected), 'classification_counts': dict(Counter(case['analysis']['classification'] for case in selected)), 'stored_float_outside_count': sum(case['analysis']['hull']['outside_stored_float_hull'] for case in selected), 'max_relative_signed_support_gap': max((case['analysis']['relative_signed_support_gap'] for case in selected), default=None)}
    return {
        'completed_trials': len(valid), 'errors': errors,
        'classification_counts': dict(Counter(item['classification'] for item in analyses)),
        'hull_dimension_counts': dict(Counter(str(item['hull']['dimension']) for item in analyses)),
        'stored_float_outside_count': sum(item['hull']['outside_stored_float_hull'] for item in analyses),
        'max_unitarity_residual': max((item['unitarity_frobenius_residual'] for item in analyses), default=None),
        'max_normality_residual': max((item['normality_relative_frobenius_residual'] for item in analyses), default=None),
        'max_spectral_equation_residual': max((item['spectral_equation_relative_frobenius_residual'] for item in analyses), default=None),
        'max_determinant_disagreement': max((item['determinant_disagreement'] for item in analyses), default=None),
        'min_max_U_imaginary_component': min((item['max_U_imaginary_component'] for item in analyses), default=None),
        'by_dimension': {str(n): group(lambda case, n=n: case['dimension'] == n) for n in range(2, 7)},
        'by_family': {family: group(lambda case, family=family: case['family'] == family) for family in FAMILIES},
        'unverified_candidate_ids': [case['id'] for case in valid if case['analysis']['classification'] == 'unverified_outside_candidate'],
        'strongest_nondegenerate_ids': [case['id'] for case in sorted((case for case in valid if case['analysis']['hull']['dimension'] == 2), key=lambda case: case['analysis']['relative_signed_support_gap'], reverse=True)[:10]],
        'interpretation': 'Finite numerical observations only. A floating point hull violation is not a counterexample. Heuristic screening allowances are not interval certificates. No universal statement is proved.',
    }


def run(seed, output):
    start_cpu, start_wall = time.process_time(), time.monotonic()
    cases, errors = [], []
    for n in range(2, 7):
        for family_index, family in enumerate(FAMILIES):
            for index in range(8):
                if time.process_time() - start_cpu >= CPU_LIMIT_SECONDS:
                    raise RuntimeError('CPU budget reached; rerun is not authorized without a new budget.')
                sample_seed = seed + n * 10000 + family_index * 100 + index
                case = {'id': f'n{n}-{family}-{index:02d}', 'dimension': n, 'family': family, 'family_index': index, 'seed': sample_seed}
                try:
                    a, b, u, detail = generate_case(n, family, index, sample_seed)
                    case.update({'construction': detail, 'a': list(map(encode, a)), 'b': list(map(encode, b)), 'U': [[encode(z) for z in row] for row in u], 'analysis': analyze_case(a, b, u)})
                except (ValueError, ArithmeticError) as error:
                    case['error'] = str(error); errors.append({'id': case['id'], 'error': str(error)})
                cases.append(case)
    report = {
        'schema': 'mdo.approximate-complex-pair-search.v1', 'approximate': True,
        'numeric_format': 'IEEE-754 binary64 complex components, JSON numbers; NOT exact Gaussian-rational witnesses',
        'base_seed': seed, 'bounds': {'dimensions': [2, 3, 4, 5, 6], 'trials_per_dimension': 40, 'max_total_trials': MAX_TRIALS, 'max_cpu_seconds': CPU_LIMIT_SECONDS},
        'environment': {'python': sys.version, 'platform': platform.platform(), 'dependencies': 'Python standard library only', 'float_mantissa_bits': sys.float_info.mant_dig},
        'timing': {'cpu_seconds': time.process_time() - start_cpu, 'wall_seconds': time.monotonic() - start_wall},
        'summary': summarize(cases, errors), 'cases': cases,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
    candidate_ids = set(report['summary']['strongest_nondegenerate_ids'] + report['summary']['unverified_candidate_ids'])
    best = {**{key: value for key, value in report.items() if key not in ('cases', 'summary')}, 'purpose': 'Complete strongest nondegenerate numerical cases for independent replay; these are not certified violations.', 'cases': [case for case in cases if case['id'] in candidate_ids]}
    output.with_name('strongest-cases.json').write_text(json.dumps(best, indent=2, allow_nan=False) + '\n')
    print(json.dumps({'output': str(output), 'timing': report['timing'], 'summary': report['summary']}, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seed', type=int, default=20261124)
    parser.add_argument('--output', type=Path, default=Path(__file__).with_name('report.json'))
    args = parser.parse_args()
    run(args.seed, args.output)


if __name__ == '__main__':
    main()
