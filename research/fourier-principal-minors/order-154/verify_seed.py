"""Contributor comparison checker for F_154. Python 3.8+, standard library."""
from functools import lru_cache
from math import comb

P = 11
ORDER = 14
# Low degree first, including each monic leading coefficient.
FACTORS = ((1, 4, 6, 1), (1, 6, 4, 1))
ZERO = (0, 0, 0)
ONE = (1, 0, 0)
X = (0, 1, 0)


def verify_factorization():
    product = [0] * 7
    for i, a in enumerate(FACTORS[0]):
        for j, b in enumerate(FACTORS[1]):
            product[i + j] += a * b
    if [c % P for c in product] != [(-1) ** i % P for i in range(7)]:
        raise ValueError("Wrong factorization of Phi_14")
    for factor in FACTORS:
        # A cubic over a field is irreducible exactly when it has no root there.
        if any(sum(c * t ** i for i, c in enumerate(factor)) % P == 0 for t in range(P)):
            raise ValueError(f"Reducible cubic: {factor}")


def multiply(left, right, factor):
    a, b, c = left
    d, e, f = right
    coeffs = [a*d, a*e + b*d, a*f + b*e + c*d, b*f + c*e, c*f]
    for k in (4, 3):
        for i in range(3):
            coeffs[k - 3 + i] -= coeffs[k] * factor[i]
    return tuple(x % P for x in coeffs[:3])


def check_chart(factor):
    def mul(a, b):
        return multiply(a, b, factor)

    def power(value, exponent):
        result = ONE
        while exponent:
            if exponent % 2:
                result = mul(result, value)
            value = mul(value, value)
            exponent //= 2
        return result

    @lru_cache(maxsize=None)
    def inverse(value):
        if value == ZERO:
            raise ZeroDivisionError("zero pivot")
        answer = power(value, P**3 - 2)
        if mul(value, answer) != ONE:
            raise ArithmeticError("field inversion failed")
        return answer

    if power(X, ORDER) != ONE or any(power(X, k) == ONE for k in (2, 7)):
        raise ValueError("root does not have exact order 14")
    powers = [power(X, k) for k in range(ORDER)]
    initial = [[powers[(i*j) % ORDER] for j in range(ORDER)] for i in range(ORDER)]
    by_size = [0] * (ORDER + 1)

    def visit(matrix, index, selected):
        # Decisions on earlier coordinates have already been eliminated. The
        # determinant of a selected submatrix is the product of its pivots.
        if not matrix:
            by_size[len(selected)] += 1
            return
        tail = [row[1:] for row in matrix[1:]]
        visit(tail, index + 1, selected)  # Skip this index.
        pivot = matrix[0][0]
        extended = selected + (index,)
        if pivot == ZERO:
            raise ArithmeticError(f"singular minor at {extended} in chart {factor}")
        inv_pivot = inverse(pivot)
        reduced = []
        for i, row in enumerate(tail):
            scale = mul(matrix[i+1][0], inv_pivot)
            new_row = []
            for j, value in enumerate(row):
                change = mul(scale, matrix[0][j+1])
                new_row.append(tuple((a - b) % P for a, b in zip(value, change)))
            reduced.append(new_row)
        visit(reduced, index + 1, extended)  # Include it and take its Schur complement.

    visit(initial, 0, ())
    if by_size != [comb(ORDER, k) for k in range(ORDER + 1)]:
        raise RuntimeError(f"incomplete subset counts: {by_size}")
    return by_size


if __name__ == "__main__":
    verify_factorization()
    for chart in FACTORS:
        counts = check_chart(chart)
        print(f"cubic {chart}: {sum(counts)} nonsingular principal minors; counts {counts}")
