"""Check the two-coordinate obstruction using original normal matrices and Schur complements."""
from fractions import Fraction as F
from check import v


def add(a, b):
    return [[x + y for x, y in zip(row, other)] for row, other in zip(a, b)]


def effective(block, j, outside):
    return [[block[i][k] - block[i][outside] * block[outside][k] / block[outside][outside]
             for k in j] for i in j]


def zero_extend(block, indices):
    result = [[v.ZERO for _ in range(4)] for _ in range(4)]
    for i, p in enumerate(indices):
        for j, q in enumerate(indices):
            result[p][q] = block[i][j]
    return result


def check():
    shift = [[v.QComplex(int(j == (i+1) % 3)) for j in range(3)] for i in range(3)]
    identity = v.identity(3)
    square = v.multiply(shift, shift)
    f = [[(identity[i][j] - shift[i][j] + square[i][j]) / 2 for j in range(3)] for i in range(3)]
    assert v.multiply(f, add(identity, shift)) == identity
    a = zero_extend(f, [0, 1, 2])
    b = zero_extend(v.adjoint(f), [1, 2, 3])
    for matrix in (a, b):
        assert v.multiply(matrix, v.adjoint(matrix)) == v.multiply(v.adjoint(matrix), matrix)
    left, right = effective(a, [1, 2], 0), effective(b, [1, 2], 3)
    for actual, target in zip((left, right), ([[1, -1], [0, 1]], [[1, 0], [-1, 1]])):
        assert actual == [[v.QComplex(x) for x in row] for row in target]
    scale = a[0][0] * b[3][3]
    actual = v.determinant(add(a, b))
    assert actual == scale * v.determinant(add(left, right)) == v.QComplex(F(3, 4))
    assert scale * v.determinant(add(left, v.identity(2))) == v.ONE
    assert F(1, 3) * F(1, 4) + F(2, 3) == F(3, 4)
    aligned = zero_extend(f, [1, 2, 3])
    assert effective(aligned, [1, 2], 3) == left
    assert v.determinant(add(a, aligned)) == scale * v.determinant(add(left, v.identity(2))) == v.ONE
    return 'normal blocks, opposite triangular Schur complements, direct determinant 3/4; spectral replacement 1'


if __name__ == '__main__':
    print(check())
