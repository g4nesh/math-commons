"""Replay a fixed distribution on C³ and a separator against a fixed distribution on C⁴.

Both certificates use exactly the same unitary and fixed second spectrum.
See the accompanying proof for the multi-affine interpolation argument.
"""

import json
from fractions import Fraction
from itertools import permutations, product
from pathlib import Path

from spectral import BASE_A, Q, Z, e, leibniz, matrix_sum, source_check, vertex


def target(diagonal):
    matrix = matrix_sum(diagonal)
    value = leibniz(matrix)
    if value != e.determinant(matrix):
        raise ValueError("Elimination and the Leibniz formula disagree")
    return value


def check_certificate(certificate):
    source_check()
    pairings = list(permutations(range(4)))
    weights = {}
    for entry in certificate['weightsAtFixedFourthEigenvalue']:
        perm = tuple(entry['columnToRow'])
        mass = Fraction(entry['mass'])
        if perm not in pairings or perm in weights or mass <= 0:
            raise ValueError("Invalid permutation or mass")
        weights[perm] = mass
    if sum(weights.values()) != 1:
        raise ValueError("Mass does not sum to one")

    for corner in product((0, 1), repeat=3):
        diagonal = [Q(bit) for bit in corner] + [BASE_A[3]]
        mixture = sum((mass * vertex(diagonal, perm) for perm, mass in weights.items()), Z)
        if mixture != target(diagonal):
            raise ValueError(f"Fixed weights fail at interpolation corner {corner}")

    # Apply one real linear functional to each of the 24 vertex functions
    # and then to the actual determinant function, at the same test diagonals.
    scores = dict.fromkeys(pairings, Fraction(0))
    target_score = Fraction(0)
    for term in certificate['separatingFunctional']:
        bits = term['diagonal']
        if len(bits) != 4 or any(type(bit) is not int or bit not in (0, 1) for bit in bits):
            raise ValueError("Expected a Boolean test diagonal of length four")
        diagonal = [Q(bit) for bit in bits]
        x, y = Fraction(term['realCoefficient']), Fraction(term['imagCoefficient'])
        actual = target(diagonal)
        target_score += x * actual.real + y * actual.imag
        for perm in pairings:
            value = vertex(diagonal, perm)
            scores[perm] += x * value.real + y * value.imag

    ceiling = max(scores.values())
    if target_score <= ceiling:
        raise ValueError("Functional does not strictly separate the target from every pairing")
    return {
        'positiveWeights': len(weights),
        'complexVariables': 3,
        'interpolationCorners': 8,
        'separatorDiagonals': len(certificate['separatingFunctional']),
        'pairingsChecked': len(scores),
        'maximumPairingScore': str(ceiling),
        'determinantScore': str(target_score),
        'strictGap': str(target_score - ceiling),
    }


if __name__ == '__main__':
    certificate = json.loads(Path(__file__).with_suffix('.json').read_text())
    print(json.dumps(check_certificate(certificate), indent=2))
