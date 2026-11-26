"""Small reproducible floating screen, not a mathematical certificate."""
import argparse
import itertools
import json
import platform
from pathlib import Path

import numpy as np
import scipy
from scipy.spatial import ConvexHull, QhullError


def screen(seed=20261126, cases=100):
    if cases < 1 or cases > 1000:
        raise ValueError('cases must be between 1 and 1000 per dimension')
    rng = np.random.default_rng(seed)
    report = {'seed': seed, 'cases_per_dimension': cases, 'screening_allowance': 1e-9,
              'python': platform.python_version(), 'numpy': np.__version__,
              'scipy': scipy.__version__, 'dimensions': []}
    for n in (5, 6):
        permutations = np.array(list(itertools.permutations(range(n))))
        records = []
        for index in range(cases):
            a = rng.standard_normal(n)
            b = rng.standard_normal(n).astype(complex)
            imaginary = rng.standard_normal(4)
            if np.any(imaginary == 0):
                raise ValueError('unexpected zero imaginary input; no case silently skipped')
            b[:4] += 1j * imaginary
            raw = rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))
            unitary, triangular = np.linalg.qr(raw)
            diagonal = np.diag(triangular)
            unitary *= diagonal / np.abs(diagonal)
            matrix_b = (unitary * b) @ unitary.conj().T
            target = np.linalg.det(np.diag(a) + matrix_b)
            vertices = np.prod(a[None, :] + b[permutations], axis=1)
            scale = max(1.0, float(np.max(np.abs(vertices))), abs(target))
            points = np.column_stack((vertices.real, vertices.imag)) / scale
            point = np.array([target.real, target.imag]) / scale
            try:
                hull = ConvexHull(points)
                margin = float(-np.max(hull.equations[:, :2] @ point + hull.equations[:, 2]))
            except QhullError as error:
                raise ValueError('degenerate case requires explicit separate analysis') from error
            records.append({'index': index, 'margin': margin,
                            'unitarity_residual': float(np.linalg.norm(unitary.conj().T @ unitary - np.eye(n))),
                            'normality_residual': float(np.linalg.norm(matrix_b @ matrix_b.conj().T - matrix_b.conj().T @ matrix_b))})
        report['dimensions'].append({'n': n, 'completed': len(records),
                                     'outside_candidates': sum(x['margin'] < -1e-9 for x in records),
                                     'minimum_margin': min(x['margin'] for x in records),
                                     'records': records})
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seed', type=int, default=20261126)
    parser.add_argument('--cases', type=int, default=100)
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    result = screen(args.seed, args.cases)
    content = json.dumps(result, indent=2) + '\n'
    if args.write:
        Path(__file__).with_name('numerical-report.json').write_text(content)
    print(content, end='')
