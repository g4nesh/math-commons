"""Replay Claude's proposed PMSH inequality on this specific real four-cycle.

Equation (★), its one-based indices and the relabeling group come from
pmsh/equivalence.md at 38eb4ff58c9b3d1caefb3ab372f9dd37cca15ab2 in PR 18.
"""

import json
from fractions import Fraction
from itertools import permutations
from pathlib import Path


def slack(u):
    def magnitude(i, j, k, l):
        return abs(u[i][k] * u[j][l] - u[i][l] * u[j][k])

    right = (magnitude(0, 2, 0, 2) ** 2 + magnitude(0, 3, 0, 3) ** 2
             + magnitude(0, 2, 1, 2) ** 2
             + (magnitude(0, 1, 0, 1) + magnitude(0, 2, 0, 3)) ** 2)
    return right - magnitude(0, 1, 2, 3) ** 2


def check():
    fixture = Path(__file__).resolve().parents[1] / 'single-cycle' / 'four-cycle-obstruction.json'
    u = [[Fraction(q) for q in row] for row in json.loads(fixture.read_text())['U']]
    orders = list(permutations(range(4)))
    values = [slack([[v[i][j] for j in cols] for i in rows])
              for v in (u, list(zip(*u))) for rows in orders for cols in orders]
    if (len(values), slack(u), min(values)) != (1152, Fraction(5779, 15625), Fraction(15136, 390625)):
        raise ValueError("Unexpected value of the proposed minor inequality on this fixture")
    return {'relabellings': len(values), 'fixtureSlack': str(slack(u)), 'minimumSlack': str(min(values))}


if __name__ == '__main__':
    print(json.dumps(check()))
