#!/usr/bin/env python3
"""Fixed mathematical fixtures, not an extra random search."""
import math
from pathlib import Path
import random
import sys
import unittest
sys.path.insert(0, str(Path(__file__).parent))
from search import (analyze_case, complex_qr, determinant_lu, determinant_permutations,
                    hull_gap, identity, permutation_points, validate_sample)


class HullFixtures(unittest.TestCase):
    def test_square_inside_boundary_outside(self):
        vertices = [0j, 1+0j, 1+1j, 1j, 0j]
        inside = hull_gap(vertices, .5+.5j)
        self.assertEqual(inside['dimension'], 2)
        self.assertEqual(inside['signed_support_gap'], -.5)
        self.assertFalse(inside['outside_stored_float_hull'])
        self.assertEqual(hull_gap(vertices, .5+0j)['signed_support_gap'], 0)
        outside = hull_gap(vertices, 1.25+.5j)
        self.assertAlmostEqual(outside['signed_support_gap'], .25)
        self.assertTrue(outside['outside_stored_float_hull'])

    def test_point_hull(self):
        self.assertEqual(hull_gap([2+3j]*4, 2+3j)['signed_support_gap'], 0)
        gap = hull_gap([2+3j]*4, 5+7j)
        self.assertEqual(gap['dimension'], 0)
        self.assertEqual(gap['signed_support_gap'], 5)

    def test_segment_endpoints_and_perpendicular(self):
        vertices = [0j, 1+0j, 2+0j, 1+0j]
        self.assertEqual(hull_gap(vertices, 1+0j)['dimension'], 1)
        self.assertEqual(hull_gap(vertices, 1+0j)['signed_support_gap'], 0)
        self.assertEqual(hull_gap(vertices, 1+.25j)['signed_support_gap'], .25)
        self.assertEqual(hull_gap(vertices, 3+0j)['signed_support_gap'], 1)
        self.assertEqual(hull_gap(vertices, -1+0j)['signed_support_gap'], 1)

    def test_thin_polygon_not_collapsed_by_epsilon(self):
        result = hull_gap([0j, 1+0j, 1+1e-30j, 1e-30j], .5+2e-30j)
        self.assertEqual(result['dimension'], 2)
        self.assertTrue(result['outside_stored_float_hull'])
        self.assertAlmostEqual(result['signed_support_gap'] / 1e-30, 1)

    def test_invalid_hulls_rejected(self):
        for points, query in [([], 0j), ([complex(float('nan'), 0)], 0j), ([0j], complex(float('inf'), 0))]:
            with self.assertRaises(ValueError): hull_gap(points, query)


class MatrixFixtures(unittest.TestCase):
    def test_two_by_two_complex_mixing_segment_formula(self):
        # U = [[sqrt(t), i*sqrt(1-t)], [i*sqrt(1-t), sqrt(t)]]
        # det = t*z_identity + (1-t)*z_swap, even for complex spectra.
        t = .25
        u = [[complex(math.sqrt(t)), 1j*math.sqrt(1-t)], [1j*math.sqrt(1-t), complex(math.sqrt(t))]]
        a, b = [1+2j, -2+.5j], [3-1j, .25+4j]
        result = analyze_case(a, b, u)
        points = permutation_points(a, b)
        expected = t*points[0][1] + (1-t)*points[1][1]
        actual = complex(*result['determinant_lu'])
        self.assertLess(abs(actual - expected), 2e-14)
        self.assertEqual(result['hull']['dimension'], 1)
        self.assertNotEqual(result['classification'], 'unverified_outside_candidate')

    def test_monomial_permutation_point(self):
        a, b = [1+1j, 2-1j, -3+2j], [4-2j, -1+3j, .5+1j]
        result = analyze_case(a, b, identity(3))
        expected = permutation_points(a, b)[0][1]
        self.assertEqual(complex(*result['determinant_lu']), expected)
        self.assertLessEqual(result['hull']['signed_support_gap'], 0)
        self.assertFalse(result['hull']['outside_stored_float_hull'])

    def test_exact_singular_sum(self):
        result = analyze_case([1+1j, 2-1j], [-1-1j, -2+1j], identity(2))
        self.assertEqual(complex(*result['determinant_lu']), 0j)
        self.assertNotEqual(result['classification'], 'unverified_outside_candidate')

    def test_lu_matches_direct_determinant_with_pivot(self):
        matrix = [[0j, 1+2j, 3+0j], [2-1j, 1j, 4+1j], [1+0j, -2+0j, .5-.5j]]
        self.assertLess(abs(determinant_lu(matrix) - determinant_permutations(matrix)), 1e-12)

    def test_qr_is_complex_and_approximately_unitary(self):
        u = complex_qr(random.Random(101), 4)
        self.assertGreater(max(abs(z.imag) for row in u for z in row), .1)
        self.assertLess(validate_sample([1j]*4, [2j]*4, u), 1e-14)

    def test_malformed_dimensions_nonfinite_and_nonunitary(self):
        invalid = [([1j], [1j, 2j], [[1+0j]]), ([1j, 2j], [1j, 2j], [[1j]]),
                   ([complex(float('nan'),0)], [0j], [[1+0j]]), ([1j], [2j], [[2+0j]])]
        for a, b, u in invalid:
            with self.assertRaises(ValueError): validate_sample(a,b,u)


if __name__ == '__main__':
    unittest.main(verbosity=2)
