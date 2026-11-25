"""Finite regression suite; author-side tests, not independent mathematical review."""
from copy import deepcopy
from fractions import Fraction as F
from itertools import permutations
import json
from pathlib import Path
import random
import subprocess
import sys
import tempfile
import unittest

from verify import (InvalidWitness, QComplex as C, ZERO, ONE, adjoint,
                    certify_hull, check_hull_certificate, determinant, diagonal,
                    identity, load_witness, multiply, parse_complex, parse_witness,
                    verify_witness)

HERE = Path(__file__).resolve().parent


def leibniz(matrix):
    """Independent determinant algorithm for small regression matrices."""
    result = ZERO
    for permutation in permutations(range(len(matrix))):
        inversions = sum(permutation[i] > permutation[j] for i in range(len(matrix)) for j in range(i + 1, len(matrix)))
        term = ONE if inversions % 2 == 0 else -ONE
        for row, col in enumerate(permutation):
            term *= matrix[row][col]
        result += term
    return result


def fixture(name):
    return load_witness(HERE / 'fixtures' / f'{name}.json')


def encoded_witness(a, b, unitary):
    return {'schemaVersion': 1, 'a': [z.encode() for z in a], 'b': [z.encode() for z in b], 'U': [[z.encode() for z in row] for row in unitary]}


class ExactArithmeticTests(unittest.TestCase):
    def test_complex_arithmetic_and_conjugation(self):
        z, w = C(F(2, 3), F(-4, 7)), C(F(5, 11), F(3, 5))
        self.assertEqual((z * w) / w, z)
        self.assertEqual(z * z.conjugate(), C(F(4, 9) + F(16, 49)))
        self.assertEqual(C(1) / C(3), C(F(1, 3)))
        with self.assertRaises(TypeError):
            C(0.1)
        with self.assertRaises(ZeroDivisionError):
            z / ZERO

    def test_elimination_against_leibniz_with_exact_complex_pivots(self):
        rng = random.Random(2611)
        for n in range(1, 5):
            for _ in range(8):
                matrix = [[C(rng.randint(-2, 2), rng.randint(-2, 2)) for _ in range(n)] for _ in range(n)]
                self.assertEqual(determinant(matrix), leibniz(matrix))
        self.assertEqual(determinant([[ZERO, C(0, 1)], [C(2), C(3)]]), C(0, -2))
        self.assertEqual(determinant([[C(1), C(0, 1)], [C(2), C(0, 2)]]), ZERO)


class HullTests(unittest.TestCase):
    def assert_certificate(self, points, target, membership):
        certificate = certify_hull(points, target)
        self.assertEqual(certificate['membership'], membership)
        self.assertTrue(check_hull_certificate(points, target, certificate))
        if membership == 'inside':
            self.assertTrue(all(isinstance(weight, F) for _, weight in certificate['weights']))
            self.assertLessEqual(len(certificate['weights']), 3)
        else:
            self.assertGreater(certificate['gap'], 0)
        return certificate

    def test_point_hull_inside_and_outside(self):
        self.assert_certificate([(2, 3), (2, 3)], (2, 3), 'inside')
        self.assert_certificate([(2, 3)], (2, 4), 'outside')

    def test_segments_horizontal_vertical_oblique_and_duplicates(self):
        cases = [([(0, 0), (2, 0), (1, 0)], (1, 0), (1, 1), (3, 0)),
                 ([(2, -3), (2, 5)], (2, F(1, 3)), (3, 0), (2, -4)),
                 ([(-1, -1), (3, 3), (0, 0)], (F(1, 3), F(1, 3)), (0, 1), (4, 4))]
        for points, inside, off_line, beyond in cases:
            self.assert_certificate(points, inside, 'inside')
            self.assert_certificate(points, off_line, 'outside')
            self.assert_certificate(points, beyond, 'outside')
            self.assert_certificate(points, points[0], 'inside')
            self.assert_certificate(points, points[1], 'inside')

    def test_polygon_interior_boundary_vertices_and_strict_outside(self):
        square = [(0, 0), (2, 0), (2, 2), (0, 2), (1, 0), (1, 1), (0, 0)]
        for target in [(1, 1), (0, 0), (2, F(1, 7)), (F(3, 5), 2)]:
            self.assert_certificate(square, target, 'inside')
        for target in [(-1, 1), (3, 1), (1, -1), (1, 3), (2 + F(1, 10**30), 1)]:
            self.assert_certificate(square, target, 'outside')

    def test_direct_certificate_checker_rejects_corrupted_weights_and_gap(self):
        points, target = [(0, 0), (2, 0), (0, 2)], (F(1, 2), F(1, 2))
        certificate = certify_hull(points, target)
        certificate['weights'][0] = (certificate['weights'][0][0], F(-1))
        with self.assertRaises(ArithmeticError):
            check_hull_certificate(points, target, certificate)
        outside = certify_hull(points, (3, 3))
        outside['gap'] = F(0)
        with self.assertRaises(ArithmeticError):
            check_hull_certificate(points, (3, 3), outside)
        with self.assertRaises(TypeError):
            certify_hull(points, (0.1, 0.2))

    def test_direct_checker_rejects_float_rounding_and_all_nonexact_scalar_fields(self):
        # A float weight rounds 2**53+3 up to 2**53+4 in ordinary Python arithmetic.
        # These are different exact singleton points; accepting that equality is unsound.
        point, target = (2**53 + 3, 0), (2**53 + 4, 0)
        with self.assertRaises(TypeError):
            check_hull_certificate([point], target, {'membership': 'inside', 'weights': [(point, 1.0)]})
        with self.assertRaises(ArithmeticError):
            check_hull_certificate([point], target, {'membership': 'inside', 'weights': [(point, F(1))]})
        valid_inside = {'membership': 'inside', 'weights': [((0, 0), F(1))]}
        for invalid in (0.0, False):
            with self.assertRaises(TypeError):
                check_hull_certificate([(invalid, 0)], (0, 0), valid_inside)
            with self.assertRaises(TypeError):
                check_hull_certificate([(0, 0)], (invalid, 0), valid_inside)
            with self.assertRaises(TypeError):
                check_hull_certificate([(0, 0)], (0, 0), {'membership': 'inside', 'weights': [((invalid, 0), F(1))]})
        for invalid_weight in (1.0, True):
            with self.assertRaises(TypeError):
                check_hull_certificate([(0, 0)], (0, 0), {'membership': 'inside', 'weights': [((0, 0), invalid_weight)]})
        outside = certify_hull([(0, 0)], (1, 0))
        for field in ('bound', 'gap', 'targetValue'):
            for value in (float(outside[field]), True):
                corrupted = deepcopy(outside); corrupted[field] = value
                with self.assertRaises(TypeError):
                    check_hull_certificate([(0, 0)], (1, 0), corrupted)
        for value in (1.0, True):
            corrupted = deepcopy(outside); corrupted['normal'] = (value, 0)
            with self.assertRaises(TypeError):
                check_hull_certificate([(0, 0)], (1, 0), corrupted)
        corrupted = deepcopy(outside); corrupted['targetValue'] = F(99)
        with self.assertRaises(ArithmeticError):
            check_hull_certificate([(0, 0)], (1, 0), corrupted)

    def test_deterministic_exact_grid_against_square_characterization(self):
        points = [(-2, -2), (2, -2), (2, 2), (-2, 2)]
        for x in range(-6, 7):
            for y in range(-6, 7):
                target = (F(x, 2), F(y, 2))
                expected = 'inside' if -2 <= target[0] <= 2 and -2 <= target[1] <= 2 else 'outside'
                self.assert_certificate(points, target, expected)


class WitnessTests(unittest.TestCase):
    def test_published_determinants_and_full_permutation_multisets(self):
        examples = [('goyal-2026', ['2/3', '2/9'], {(6, -2), (2, 0), (2, 1), (0, -1), (0, 2), (2, -2)}),
                    ('zhang-ma-2026', ['2864/49', '1840/49'], {(-30, -50), (280, 280), (-130, 10), (-28, -56), (232, 376), (24, -20)})]
        for name, expected, points in examples:
            result = verify_witness(fixture(name))
            self.assertEqual(result['determinant'], expected)
            self.assertEqual(result['membership'], 'inside')
            self.assertEqual(result['permutationCount'], 6)
            self.assertEqual({tuple(F(x) for x in entry['point']) for entry in result['permutationPoints']}, points)
            matrix = [[parse_complex(z) for z in row] for row in result['AplusB']]
            self.assertEqual(leibniz(matrix).encode(), expected)
            self.assertEqual(result, json.loads((HERE / 'results' / f'{name}.json').read_text()))

    def test_genuinely_complex_n4_householder_and_independent_determinant(self):
        witness = fixture('complex-householder-n4')
        a, b, unitary, _ = parse_witness(witness)
        self.assertTrue(any(z.imag for row in unitary for z in row))
        self.assertEqual(adjoint(unitary), unitary)
        self.assertNotEqual([list(row) for row in zip(*unitary)], adjoint(unitary))
        self.assertEqual(multiply(unitary, unitary), identity(4))
        result = verify_witness(witness)
        self.assertEqual(result['determinant'], ['-880/27', '800/81'])
        self.assertEqual(result['membership'], 'inside')
        self.assertEqual(result['permutationCount'], 24)
        self.assertEqual(leibniz([[parse_complex(z) for z in row] for row in result['AplusB']]).encode(), result['determinant'])
        self.assertEqual(result, json.loads((HERE / 'results/complex-householder-n4.json').read_text()))

    def test_n1_singular_sum_and_point_hull(self):
        witness = encoded_witness([C(1, 1)], [C(-1, -1)], [[C(0, 1)]])
        result = verify_witness(witness)
        self.assertEqual(result['determinant'], ['0', '0'])
        self.assertTrue(result['singularSum'])
        self.assertEqual(result['hullDimension'], 0)
        self.assertEqual(result['convexCombination'][0]['weight'], '1')

    def test_n2_real_rotation_segment_and_repeated_spectra(self):
        rotation = [[C(F(3, 5)), C(F(4, 5))], [C(F(-4, 5)), C(F(3, 5))]]
        result = verify_witness(encoded_witness([C(0), C(2)], [C(1), C(3)], rotation))
        self.assertEqual(result['hullDimension'], 1)
        self.assertEqual(result['determinant'], ['189/25', '0'])
        self.assertEqual({entry['weight'] for entry in result['convexCombination']}, {'9/25', '16/25'})
        repeated = verify_witness(encoded_witness([C(0), C(0)], [C(1, 1), C(2, -1)], rotation))
        self.assertEqual(repeated['uniquePointCount'], 1)
        self.assertEqual(repeated['membership'], 'inside')

    def test_dimension_six_enumerates_all_720_permutations(self):
        result = verify_witness(encoded_witness([C(i) for i in range(6)], [C(i + 1) for i in range(6)], identity(6)))
        self.assertEqual(result['permutationCount'], 720)
        self.assertEqual(result['membership'], 'inside')

    def test_reject_nonunitary_without_rounding(self):
        witness = fixture('goyal-2026')
        witness['U'][0][0] = ['-333333333333/1000000000000', '0']
        with self.assertRaisesRegex(InvalidWitness, 'not exactly unitary'):
            verify_witness(witness)

    def test_reject_malformed_dimensions_values_unknown_keys_and_wrong_target(self):
        original = fixture('goyal-2026')
        bad_inputs = []
        for entry in [1, 0.5, True, '0.5', '1e-3', '1/0', '1/-2', '9' * 65, ['1']]:
            bad = deepcopy(original)
            bad['a'][0] = entry if isinstance(entry, list) else [entry, '0']
            bad_inputs.append(bad)
        for field, value in [('a', []), ('b', []), ('U', [[['1', '0']]]), ('schemaVersion', True), ('unexpected', 'ignored?')]:
            bad = deepcopy(original); bad[field] = value; bad_inputs.append(bad)
        bad_inputs.append(encoded_witness([C(1)] * 7, [C(2)] * 7, identity(7)))
        bad = deepcopy(original); bad['expectedDeterminant'] = ['0', '0']; bad_inputs.append(bad)
        for bad in bad_inputs:
            with self.assertRaises(InvalidWitness):
                verify_witness(bad)

    def test_duplicate_keys_invalid_json_and_cli_failure(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / 'bad.json'
            path.write_text('{"schemaVersion":1,"schemaVersion":1}')
            with self.assertRaisesRegex(InvalidWitness, 'Duplicate'):
                load_witness(path)
            path.write_text('{"schemaVersion":NaN}')
            with self.assertRaises(InvalidWitness):
                load_witness(path)
            path.write_text('{}')
            result = subprocess.run([sys.executable, str(HERE / 'verify.py'), str(path)], capture_output=True, text=True)
            self.assertEqual(result.returncode, 2)
            self.assertFalse(json.loads(result.stderr)['certified'])
            self.assertEqual(result.stdout, '')


if __name__ == '__main__':
    unittest.main()
