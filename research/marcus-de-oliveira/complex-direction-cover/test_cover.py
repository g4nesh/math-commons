"""Independent evaluations of the explicit probability formula at boundaries and interiors."""
import copy
from fractions import Fraction as F
import json
from pathlib import Path
import sys
import unittest
import check as c

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "four-cycle-spectral"))
import spectral as s


class CoverTests(unittest.TestCase):
    def test_actual_matrix_determinants_at_patch_boundaries_and_negative_radii(self):
        for a, b in [(0, 0), (2, 0), (0, -2), (2, 1), (-2, 2), (1, 2),
                     (F(-3, 2), F(1, 2)), (F(-1, 3), F(-2, 1))]:
            weights = c.weights_at(a, b)
            self.assertEqual(sum(weights.values()), 1)
            self.assertTrue(all(w > 0 for w in weights.values()))
            for x, y in [(s.Q(-1, 2), s.Q(3, -1)), (s.Z, s.Z), (s.Q(3), s.Q(3))]:
                diag = [x, s.Q(a, b), y, s.Z]
                direct = s.e.determinant(s.matrix_sum(diag))
                mix = sum((w*s.vertex(diag, tuple(int(i)-1 for i in word))
                           for word, w in weights.items()), s.Z)
                self.assertEqual(direct, mix)

    def test_inexact_input_and_uncovered_direction(self):
        with self.assertRaises(TypeError):
            c.weights_at(0.5, 2)
        with self.assertRaises(ValueError):
            c.weights_at(1, 1)

    def test_changed_mixed_weight_breaks_the_certificate(self):
        patch = json.loads((c.HERE / "patches.json").read_text())[0]
        changed = copy.deepcopy(patch)
        changed["terms"][0][3] = str(F(changed["terms"][0][3])+1)
        with self.assertRaises(AssertionError):
            c.verify_patch(changed)


if __name__ == "__main__":
    unittest.main()
