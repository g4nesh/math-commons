"""Certify the whole complex two-eigenvalue family by a quadratic identity."""
from fractions import Fraction as F
from itertools import combinations_with_replacement, permutations
from pathlib import Path
import json
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "four-cycle-spectral"))
import spectral as s

PERMS = tuple(permutations(range(4)))
MONOMIALS = tuple(combinations_with_replacement(range(8), 2))


def gap_rows(third=s.Q(2, F(-1, 4)), fourth=s.Z):
    """Each row is Re(sum L_k * (vertex-target)_k) in the eight real coordinates."""
    s.source_check()
    samples = []
    for x, y in ((0, 0), (1, 0), (0, 1), (1, 1)):
        diagonal = [s.Q(x), s.Q(y), third, fourth]
        target = s.e.determinant(s.matrix_sum(diagonal))
        samples.append([s.vertex(diagonal, p)-target for p in PERMS])
    rows = []
    for constant, x_value, y_value, both in zip(*samples):
        coeff = [constant, x_value-constant, y_value-constant, both-x_value-y_value+constant]
        rows.append([c for z in coeff for c in (z.real, -z.imag)])
    return rows


def quadratic_product(a, b):
    return [a[i]*b[j] + (a[j]*b[i] if i != j else 0) for i, j in MONOMIALS]


def determinant_parts():
    """Coefficients of Re and Im of L_0 L_3 - L_1 L_2, with complex entries."""
    real = {(0, 6):1, (1, 7):-1, (2, 4):-1, (3, 5):1}
    imag = {(0, 7):1, (1, 6):1, (2, 5):-1, (3, 4):-1}
    return [real.get(m, 0) for m in MONOMIALS], [imag.get(m, 0) for m in MONOMIALS]


def modular_determinant(columns, prime=1_000_000_007):
    """A nonzero determinant modulo the prime certifies rational nonsingularity."""
    rows = [[F(x).numerator*pow(F(x).denominator, -1, prime) % prime for x in row]
            for row in zip(*columns)]
    assert len(rows) == len(columns)
    det = 1
    for j in range(len(rows)):
        k = next((i for i in range(j, len(rows)) if rows[i][j]), None)
        if k is None:
            return 0
        if k != j:
            rows[k], rows[j] = rows[j], rows[k]
            det = -det
        pivot = rows[j][j]
        det = det*pivot % prime
        inv = pow(pivot, -1, prime)
        for row in rows[j+1:]:
            ratio = row[j]*inv % prime
            row[j:] = [(a-ratio*b) % prime for a, b in zip(row[j:], rows[j][j:])]
    return det


def verify(certificate):
    assert certificate["fixedThird"] == ["2", "-1/4"] and certificate["fixedFourth"] == ["0", "0"]
    rows = gap_rows()
    witness = list(map(F, certificate["strictSeparator"]))
    assert len(witness) == 8
    margins = [sum(a*b for a, b in zip(witness, row)) for row in rows]
    assert min(margins) == F(certificate["minimumGap"]) > 0
    identifiers = {"".join(str(i+1) for i in p): j for j, p in enumerate(PERMS)}
    products = [F(0)]*len(MONOMIALS)
    terms = certificate["products"]
    assert 0 < len(terms) <= 36
    pairs = set()
    total = F(0)
    columns = []
    for left, right, value in terms:
        i, j = identifiers[left], identifiers[right]
        assert i <= j and (i, j) not in pairs
        pairs.add((i, j))
        weight = F(value)
        assert weight > 0
        total += weight
        q = quadratic_product(rows[i], rows[j])
        columns.append(q+[1])
        products = [x+weight*y for x, y in zip(products, q)]
    assert total == 1
    real_scale, imag_scale = map(F, certificate["phase"])
    assert real_scale != 0 or imag_scale != 0
    real, imag = determinant_parts()
    assert products == [real_scale*a+imag_scale*b for a, b in zip(real, imag)]
    columns.extend([[-a for a in real]+[0], [-b for b in imag]+[0]])
    residue = modular_determinant(columns)
    assert residue != 0
    return {"strictlyPositiveProducts": len(terms), "unrestrictedRealCoefficients": len(products),
            "uniformLawSeparationMargin": str(min(margins)), "permutations": len(margins),
            "determinantMod1000000007": residue}


if __name__ == "__main__":
    print(json.dumps(verify(json.loads((HERE / "certificate.json").read_text())), indent=2))
