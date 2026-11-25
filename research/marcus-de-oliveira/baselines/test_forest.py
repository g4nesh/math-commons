"""Bounded exact regressions and rejection tests for the forest construction."""
from fractions import Fraction as F
import unittest

import forest


class ForestTests(unittest.TestCase):
    def test_all_fixed_compounds_and_determinants(self):
        checked = [forest.check_fixture(fixture) for fixture in forest.fixtures()]
        self.assertEqual(len(checked), 11)
        self.assertEqual(sum(item["compoundEntriesChecked"] for item in checked), 3332)

    def test_complex_phases_are_present(self):
        fixture = next(f for f in forest.fixtures() if f["name"] == "complex-two-coordinate")
        self.assertTrue(any(z.imag for row in forest.build_unitary(2, fixture["gates"]) for z in row))

    def test_single_block_explicit_weights(self):
        gates = [(0, 1, forest.phased_block())]
        self.assertEqual(forest.permutation_weights(2, gates),
                         {(0, 1): F(9, 25), (1, 0): F(16, 25)})

    def test_cycle_obstruction_is_not_moc_violation(self):
        result = forest.cycle_obstruction()
        self.assertEqual(result["difference"], "1152/3125")
        self.assertTrue(0 < F(result["actualSquaredModulusAndDeterminant"]) < 1)

    def test_repeated_edge_rejected(self):
        gate = (0, 1, forest.phased_block())
        with self.assertRaisesRegex(ValueError, "Repeated edge or cycle"):
            forest.build_unitary(2, [gate, gate])

    def test_cycle_rejected(self):
        block = forest.phased_block()
        with self.assertRaisesRegex(ValueError, "Repeated edge or cycle"):
            forest.permutation_weights(3, [(0, 1, block), (1, 2, block), (0, 2, block)])

    def test_bad_endpoints_rejected(self):
        for edge in ((0, 0), (1, 0), (-1, 1), (0, 2), (False, 1)):
            with self.subTest(edge=edge), self.assertRaises(ValueError):
                forest.build_unitary(2, [(*edge, forest.phased_block())])

    def test_bad_dimensions_rejected(self):
        for n in (0, 7, True, 2.0):
            with self.subTest(n=n), self.assertRaises(ValueError):
                forest.build_unitary(n, [])

    def test_nonunitary_and_malformed_block_rejected(self):
        for block in ([[forest.O, forest.O], [forest.Z, forest.O]], [[forest.O]]):
            with self.subTest(block=block), self.assertRaises(ValueError):
                forest.build_unitary(2, [(0, 1, block)])

    def test_inexact_scalar_rejected(self):
        with self.assertRaises(TypeError):
            forest.Q(0.6)
        with self.assertRaises(TypeError):
            forest.build_unitary(2, [(0, 1, [[0.6, 0.8], [-0.8, 0.6]])])

    def test_zero_and_singular_determinants(self):
        self.assertEqual(forest.det([]), forest.O)
        self.assertEqual(forest.determinant_leibniz([]), forest.O)
        matrix = [[forest.O, forest.Q(0, 1)], [forest.O, forest.Q(0, 1)]]
        self.assertEqual(forest.det(matrix), forest.Z)
        self.assertEqual(forest.determinant_leibniz(matrix), forest.Z)


if __name__ == "__main__":
    unittest.main()
