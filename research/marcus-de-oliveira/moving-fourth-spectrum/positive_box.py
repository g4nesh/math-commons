"""Uniform convex weights on a two-dimensional neighborhood of a₃ = 1 − i."""

import argparse
import json
import sys
from fractions import Fraction as F
from itertools import product
from pathlib import Path

PARENT = Path(__file__).resolve().parents[1] / 'four-cycle-spectral'
sys.path.insert(0, str(PARENT))
from fiber_certificate import target
from spectral import Q, source_check, vertex

CENTER = Q(1, -1)
RADIUS = F(1, 1024)


def parent_weights():
    fixture = json.loads((PARENT / 'fiber_certificate.json').read_text())
    terms = fixture['weightsAtFixedFourthEigenvalue']
    return [tuple(t['columnToRow']) for t in terms], [F(t['mass']) for t in terms]


def corner_system(z, support):
    rows = [[F(1)] * (len(support) + 1)]
    for bits in product((0, 1), repeat=3):
        diagonal = [Q(bit) for bit in bits] + [z]
        values = [vertex(diagonal, p) for p in support] + [target(diagonal)]
        rows.extend([[value.real for value in values], [value.imag for value in values]])
    return rows


def subtract(a, b):
    return [[x - y for x, y in zip(r, s)] for r, s in zip(a, b)]


def pivots(matrix):
    """Indices of independent original rows, by exact elimination."""
    work = [list(row) for row in matrix]
    order = list(range(len(work)))
    rank = 0
    for j in range(len(work[0])):
        i = next((i for i in range(rank, len(work)) if work[i][j]), None)
        if i is None:
            continue
        work[rank], work[i] = work[i], work[rank]
        order[rank], order[i] = order[i], order[rank]
        for i in range(rank + 1, len(work)):
            if work[i][j]:
                ratio = work[i][j] / work[rank][j]
                for k in range(j + 1, len(work[i])):
                    work[i][k] -= ratio * work[rank][k]
        rank += 1
    return order[:rank]


def left_solve(a, c):
    """Solve A X = C with a square nonsingular rational A."""
    n = len(a)
    work = [list(row) + list(rhs) for row, rhs in zip(a, c)]
    for j in range(n):
        i = next(i for i in range(j, n) if work[i][j])
        work[i], work[j] = work[j], work[i]
        work[j] = [x / work[j][j] for x in work[j]]
        for i in range(n):
            if i != j and work[i][j]:
                multiplier = work[i][j]
                work[i] = [x - multiplier * y for x, y in zip(work[i], work[j])]
    return [row[n:] for row in work]


def verify():
    source_check()
    support, w = parent_weights()
    n = len(w)
    if n != 12 or len(set(support)) != n or min(w) <= 0 or sum(w) != 1:
        raise ValueError("Expected the twelve positive archived masses")

    base = corner_system(CENTER, support)
    dx = subtract(corner_system(CENTER + Q(1), support), base)
    dy = subtract(corner_system(CENTER + Q(0, 1), support), base)
    if any(sum(mass * c for mass, c in zip(w, row)) != row[-1] for row in base):
        raise ValueError("The archived masses fail a base equation")
    selected = pivots([r[:-1] for r in base])
    if len(selected) != n:
        raise ValueError("Base coefficient matrix is rank deficient")
    a = [base[i][:-1] for i in selected]

    # Every entry is real affine in x,y, where z = 1-i+x+iy. Every
    # (n+1)-minor of the augmented system has degree at most n+1 in each.
    # Vanishing on an (n+2)-by-(n+2) grid forces each minor to vanish identically.
    for x, y in product(range(n + 2), repeat=2):
        grid = [[c + x * u + y * v for c, u, v in zip(r, s, t)] for r, s, t in zip(base, dx, dy)]
        if len(pivots(grid)) > n:
            raise ValueError(f"Augmented rank exceeds {n} at {(x, y)}")

    # Write weights as w + diag(w)v. Precondition the selected equations
    # by A(0)^{-1} so v solves (I + x Mx + y My)v = x rx + y ry.
    systems = []
    for derivative in (dx, dy):
        columns = []
        for i in selected:
            row = derivative[i]
            scaled = [mass * c for mass, c in zip(w, row)]
            columns.append(scaled + [row[-1] - sum(scaled)])
        systems.append([[value / w[i] for value in row]
                        for i, row in enumerate(left_solve(a, columns))])
    h = max(sum(abs(u) + abs(v) for u, v in zip(rx[:-1], ry[:-1])) for rx, ry in zip(*systems))
    k = max(abs(rx[-1]) + abs(ry[-1]) for rx, ry in zip(*systems))
    radius = F(1)
    while radius * (h + k) >= 1:
        radius /= 2
    deviation = radius * k / (1 - radius * h)
    if radius != RADIUS or deviation >= 1:
        raise ValueError("Perturbation bound cannot certify strictly positive weights")
    return {'radius': str(radius), 'maxRelativeWeightChange': str(deviation),
            'matrixNormBound': str(h), 'residualNormBound': str(k),
            'selectedRows': selected, 'rankTestGridPoints': (n + 2) ** 2}


def weights_at(x, y):
    """Return exact weights for a rational displacement (x,y) from 1−i."""
    if max(abs(x), abs(y)) > RADIUS:
        raise ValueError("Displacement is outside the certified closed square")
    support, _ = parent_weights()
    selected = pivots([r[:-1] for r in corner_system(CENTER, support)])
    system = corner_system(CENTER + Q(x, y), support)
    square = [system[i] for i in selected]
    solution = left_solve([r[:-1] for r in square], [[r[-1]] for r in square])
    weights = [r[0] for r in solution]
    if min(weights) <= 0 or any(sum(v * c for v, c in zip(weights, row)) != row[-1] for row in system):
        raise ValueError("Weights failed positivity or exact reconstruction")
    return [{'columnToRow': p, 'weight': str(w)} for p, w in zip(support, weights)]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--offset', nargs=2, type=F, metavar=('X', 'Y'),
                        help='return the weights at fourth eigenvalue 1-i+X+iY')
    args = parser.parse_args()
    print(json.dumps(verify() if args.offset is None else weights_at(*args.offset), indent=2))
