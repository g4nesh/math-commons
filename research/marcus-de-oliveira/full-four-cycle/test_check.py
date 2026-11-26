"""Focused exact regressions; run with unittest discovery from repository root."""
import copy
from fractions import Fraction as F
import importlib.util
from pathlib import Path
import unittest

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('full_four_cycle_check', HERE / 'check.py')
c = importlib.util.module_from_spec(spec)
spec.loader.exec_module(c)


class CertificateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = c.load()

    def test_all_coefficients_and_report(self):
        import json
        self.assertEqual(c.report(self.data), json.loads((HERE / 'report.json').read_text()))

    def test_corrupted_coefficient_rejected(self):
        changed = copy.deepcopy(self.data)
        changed['terms'][0]['mass'] = str(F(changed['terms'][0]['mass']) + F(1, 10**9))
        with self.assertRaises(ArithmeticError):
            c.verify_plane(changed)

    def test_negative_and_duplicate_terms_rejected(self):
        changed = copy.deepcopy(self.data)
        changed['terms'][0]['mass'] = '-1'
        with self.assertRaises(ValueError):
            c.verify_plane(changed)
        changed = copy.deepcopy(self.data)
        changed['terms'].append(changed['terms'][0])
        with self.assertRaises(ValueError):
            c.verify_plane(changed)

    def test_exact_api_rejects_rounded_parameters(self):
        for bad in (1.0, True, complex(1, 2), '1/2'):
            with self.assertRaises(TypeError):
                c.weights_at(bad, 0, self.data)
        changed = copy.deepcopy(self.data)
        changed['terms'][0]['mass'] = 0.5
        with self.assertRaises(ValueError):
            c.verify_plane(changed)

    def test_far_complex_motion_and_free_variables(self):
        for s, t in ((10**6, -(10**5)), (F(3, 17), F(11, 23)), (0, 1)):
            w = c.weights_at(s, t, self.data)
            for x, y in ((c.Z, c.O), (c.Q(-20, 31), c.Q(97, -103))):
                c.direct_check([x, c.Z, y, c.Q(s, t)], w)
        self.assertNotEqual(c.weights_at(0, 0, self.data), c.weights_at(1, 1, self.data))

    def test_nonzero_second_entire_plane(self):
        import json
        data = json.loads((HERE / 'certificate-one.json').read_text())
        self.assertEqual(c.verify_plane(data)['positiveSquareCoefficients'], 56)
        for s, t in ((0, 0), (1000, -2000), (F(2, 7), F(-19, 11))):
            c.direct_check([c.Q(11, 5), c.O, c.Q(-7, 3), c.Q(s, t)], c.weights_at(s, t, data))

    def test_repeated_middle_hypothesis_is_sharp(self):
        c.verify_repeated(self.data)
        w = tuple(map(F, self.data['repeatedMiddleWeights']))
        c.direct_check([c.Q(2, 3), c.Q(-11, 5), c.Q(7, -1), c.Q(-11, 5)], w)
        # These same singleton-preserving weights cannot cover an unequal pair.
        with self.assertRaises(ArithmeticError):
            c.direct_check([c.O, c.Z, c.Z, c.O], w)

    def test_translation_invariant_quadratic_basis(self):
        # Six pair coefficients have zero vertex sums iff they have the stated q,r form.
        for q, r in ((c.Q(2, 3), c.Q(-7, 11)), (c.O, c.Z), (c.Z, c.O)):
            coeff = {(0, 1): q+r, (0, 2): -r, (0, 3): -q, (1, 2): -q, (1, 3): -r, (2, 3): q+r}
            self.assertTrue(all(sum((v for ij, v in coeff.items() if i in ij), c.Z) == c.Z for i in range(4)))
            a = [c.Q(2, 1), c.Q(-3, 2), c.Q(5, -7), c.Q(11, 13)]
            def value(a):
                return sum((v*a[i]*a[j] for (i, j), v in coeff.items()), c.Z)
            self.assertEqual(value(a), q*(a[0]-a[2])*(a[1]-a[3])+r*(a[0]-a[3])*(a[1]-a[2]))
            self.assertEqual(value(a), value([v+c.Q(9, -17) for v in a]))


if __name__ == '__main__':
    unittest.main()
