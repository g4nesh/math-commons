"""Bounded regressions for exact scalar bounds and certificate failure modes."""
from fractions import Fraction as F
from pathlib import Path
import copy
import importlib.util
import json
import unittest

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('radius_check', HERE/'check.py')
c = importlib.util.module_from_spec(spec)
spec.loader.exec_module(c)


class RadiusTests(unittest.TestCase):
    def test_claimed_radius_satisfies_both_regions(self):
        result = c.bounds()
        self.assertEqual(result['baseMassFloor'], F(1, 1210000000))
        for region in ('finite', 'infinity'):
            self.assertLess(result[region+'NeumannBound'], F(1, 2))
            self.assertLess(result[region+'MassChangeRelativeToFloor'], F(1, 2))

    def test_unsupported_radius_is_rejected(self):
        for radius in (0, True, 1e-100, F(0), F(-1), F(1)):
            with self.assertRaises(ValueError):
                c.bounds(radius)
        with self.assertRaises(ArithmeticError):
            c.bounds(F(1, 10**80))

    def test_reciprocal_identity_and_changed_basis(self):
        rank = json.loads((c.ROOT/'projective-fiber-stability/uniform/rank.json').read_text())
        self.assertEqual(c.verify_infinity(rank), F(153348616, 3125))
        changed = copy.deepcopy(rank)
        changed['infinityBasis'][1] = changed['infinityBasis'][0]
        with self.assertRaises(ArithmeticError):
            c.verify_infinity(changed)

    def test_radial_center_bound_is_checked(self):
        patches = json.loads((c.ROOT/'complex-direction-cover/patches.json').read_text())
        corrections = json.loads((c.ROOT/'projective-fiber-stability/uniform/corrections.json').read_text())
        self.assertEqual(len(c.augmented_patches(patches, corrections)), 8)
        changed = copy.deepcopy(patches)
        term = next(term for term in changed[0]['terms'] if term[2])
        term[2][0] = '20'
        with self.assertRaises(ArithmeticError):
            c.augmented_patches(changed, corrections)

    def test_signed_weights_are_real_normalized_but_not_assumed_positive(self):
        masses = c.signed_masses()
        self.assertEqual(sum(masses), 1)
        self.assertEqual(sum(map(abs, masses)), F(10718023, 9765625))
        self.assertTrue(any(m < 0 for m in masses))

    def test_rational_determinant_handles_pivot_and_singularity(self):
        self.assertEqual(c.real_det([[F(0), F(2)], [F(3), F(4)]]), -6)
        self.assertEqual(c.real_det([[F(1, 2), F(1)], [F(1), F(2)]]), 0)


if __name__ == '__main__':
    unittest.main()
