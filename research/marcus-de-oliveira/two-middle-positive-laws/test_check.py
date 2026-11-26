import copy
from fractions import Fraction as F
import importlib.util
from pathlib import Path
import unittest

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('two_middle_check', HERE/'check.py')
c = importlib.util.module_from_spec(spec)
spec.loader.exec_module(c)


class UniformTubeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = c.prepared()

    def test_full_proof_and_archived_report(self):
        self.assertEqual(c.report(), c.load('report.json'))
        self.assertLess(self.data['relative'], F(1, 6))

    def test_nonzero_second_and_unbounded_fourth(self):
        r = self.data['radius']
        for u, v, s, t in ((r, -r, 100000, -1000000), (r/3, r/5, F(13, 7), F(-23, 11))):
            weights = c.weights_at(u, v, s, t)
            self.assertEqual(sum(weights), 1)
            for p, w in zip(c.PERMS, weights):
                if p[3] == 1:
                    self.assertEqual(w, 0)
                else:
                    self.assertGreaterEqual(w, F(1, 12000))
            for x, y in ((c.Z, c.O), (c.Q(19, -3), c.Q(-17, 29))):
                c.prior.direct_check([x, c.Q(u, v), y, c.Q(s, t)], weights)

    def test_boundary_included_and_outside_rejected(self):
        r = self.data['radius']
        self.assertTrue(all(v >= 0 for v in c.coefficients_at(r, r)))
        for u, v in ((r+F(1,10**9), 0), (0, -r-F(1,10**9))):
            with self.assertRaises(ValueError):
                c.coefficients_at(u, v)

    def test_exact_input_contract(self):
        for value in (True, 0.0, complex(0), '0'):
            with self.assertRaises(TypeError):
                c.weights_at(value, 0, 1, 0)
            with self.assertRaises(TypeError):
                c.weights_at(0, 0, value, 0)

    def test_center_corruption_is_detected(self):
        center = copy.deepcopy(c.load('center.json'))
        j = center['support'][0]
        center['coefficients'][j] = str(F(center['coefficients'][j])+F(1,10**9))
        with self.assertRaises(ArithmeticError):
            c.verify(center, c.load('invariance.json'))

    def test_invariant_action_corruption_is_detected(self):
        witness = copy.deepcopy(c.load('invariance.json'))
        witness['actionReal'][0][0] = str(F(witness['actionReal'][0][0])+1)
        with self.assertRaises(ArithmeticError):
            c.verify(c.load('center.json'), witness)

    def test_collision_face_for_general_second_entry(self):
        delta = c.Q(0,1)
        null = [row[3] for row in c.sp.U]
        for eta in (c.Z, c.O, c.Q(17,-23)):
            a = [delta, eta, delta, delta]
            matrix = c.sp.matrix_sum(a)
            self.assertTrue(all(sum((x*y for x,y in zip(row,null)), c.Z) == c.Z for row in matrix))
            self.assertEqual(c.sp.e.determinant(matrix), c.Z)
            nonzero = (eta-delta)*c.sp.product(delta+b for b in c.sp.B_SPECTRUM[:3])
            self.assertNotEqual(nonzero, c.Z)
            for p in c.PERMS:
                self.assertEqual(c.sp.vertex(a,p), nonzero if p[3] == 1 else c.Z)


if __name__ == '__main__':
    unittest.main()
