import copy
from fractions import Fraction as F
import importlib.util
from pathlib import Path
import unittest

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('split_pair_check', HERE / 'check.py')
check = importlib.util.module_from_spec(spec)
spec.loader.exec_module(check)


class SplitPairTests(unittest.TestCase):
    def test_full_exact_fixture(self):
        result = check.verify()
        self.assertEqual(result['squared_minors_checked'], 924)
        self.assertEqual(result['zonotope_normal_directions'], 31)

    def test_perturbed_complex_unitary_stays_positive(self):
        data = check.load()
        u = check.fixture(data)
        t = F(data['radius']) / 10
        c, s = (1 - t*t)/(1 + t*t), 2*t/(1 + t*t)
        rotation = check.exact.identity(6)
        rotation[0][0] = rotation[5][5] = check.Q(c)
        rotation[0][5] = rotation[5][0] = check.Q(0, s)
        moved = check.exact.multiply(rotation, u)
        weights = check.weights_at(moved, data)
        self.assertGreaterEqual(min(weights), F(1, 2000))
        canonical = check.canonical(moved)
        self.assertTrue(all(check.dot(row, [x-y for x,y in zip(weights, canonical)]) == 0
                            for row in check.structure()['matrix']))

    def test_outside_ball_is_rejected(self):
        with self.assertRaisesRegex(ValueError, 'outside certified ball'):
            check.weights_at(check.exact.identity(6))

    def test_missing_normal_is_rejected(self):
        data = check.load()
        data['normals'].pop()
        with self.assertRaisesRegex(ValueError, 'all facet normals'):
            check.verify(data, minors=False)

    def test_mass_bound_violation_is_rejected(self):
        data = check.load()
        imbalances = list(map(F, data['imbalances']))
        imbalances[0] = F(10)
        with self.assertRaisesRegex(ValueError, 'split-pair bounds'):
            check.lift(list(map(F, data['canonical'])), list(map(F, data['even_repairs'])), imbalances)

    def test_corrupt_correction_is_rejected(self):
        data = copy.deepcopy(check.load())
        data['weights'][0] = str(F(data['weights'][0]) + F(1, 100))
        with self.assertRaises(ValueError):
            check.verify(data, minors=False)


if __name__ == '__main__':
    unittest.main()
