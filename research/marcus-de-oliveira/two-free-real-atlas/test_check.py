import copy
from fractions import Fraction as F
from pathlib import Path
import sys
import unittest

import check


class RealLineCertificateTests(unittest.TestCase):
    def setUp(self):
        self.certificate = check.load_certificate()

    def test_full_polynomial_identity(self):
        report = check.verify(self.certificate)
        self.assertEqual(report['positiveSquareTerms'], 29)
        self.assertEqual(report['usedPermutations'], 14)
        self.assertTrue(report['projectiveInfinityCovered'])

    def test_original_matrices_at_nonboolean_complex_spectra(self):
        # A separate archived implementation uses Gaussian elimination and
        # independently constructs the original matrix, not polynomial evaluation.
        sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'four-cycle-spectral'))
        import spectral as s
        for t in (-1000000, -10, -2, 0, F(1, 2), 1, 3, 1000000):
            a = [s.Q(F(7, 3), F(-2, 5)), s.Q(t), s.Q(F(-11, 4), F(13, 7)), s.Z]
            target = s.e.determinant(s.matrix_sum(a))
            value = sum((weight*s.vertex(a, perm) for perm, weight in
                         check.weights_at(t, self.certificate)), s.Z)
            self.assertEqual(target, value)

    def test_rank_drop_points_and_infinity(self):
        for u, v in ((-2, 1), (0, 1), (1, 0)):
            weights = check.weights_homogeneous(u, v, self.certificate)
            self.assertEqual(sum(w for _, w in weights), 1)
            self.assertTrue(all(w > 0 for _, w in weights))

    def test_projective_scaling_and_reciprocal_chart(self):
        expected = check.weights_homogeneous(3, 2, self.certificate)
        self.assertEqual(expected, check.weights_homogeneous(-9, -6, self.certificate))
        self.assertEqual(expected, check.weights_homogeneous(1, F(2, 3), self.certificate))

    def test_tampered_mass_and_center_are_rejected(self):
        for field in ('mass', 'center'):
            bad = copy.deepcopy(self.certificate)
            bad['terms'][0][field] = str(F(bad['terms'][0][field])+F(1, 1000))
            with self.assertRaises(ValueError):
                check.verify(bad)

    def test_wrong_permutation_and_negative_mass_are_rejected(self):
        for field, value in (('permutation', [0, 0, 1, 2]), ('mass', '-1')):
            bad = copy.deepcopy(self.certificate)
            bad['terms'][0][field] = value
            with self.assertRaises(ValueError):
                check.verify(bad)

    def test_denominator_tamper_and_unsafe_inputs_are_rejected(self):
        bad = copy.deepcopy(self.certificate)
        bad['denominator'][0] = '0'
        with self.assertRaises(ValueError):
            check.verify(bad)
        with self.assertRaises(ValueError):
            check.weights_homogeneous(0, 0)
        for t in (0.5, True, complex(1, 1)):
            with self.assertRaises(TypeError):
                check.weights_at(t)


if __name__ == '__main__':
    unittest.main()
