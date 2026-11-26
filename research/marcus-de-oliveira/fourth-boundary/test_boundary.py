import copy
from fractions import Fraction as F
from pathlib import Path
import json
import unittest
import boundary as b
import discover


class BoundaryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.positive=json.loads((b.HERE/'positive.json').read_text())
        cls.dual=json.loads((b.HERE/'dual.json').read_text())

    def test_full_exact_replay_matches_report(self):
        self.assertEqual(b.verify(),json.loads((b.HERE/'verification.json').read_text()))

    def test_rational_weights_cover_all_three_formulas(self):
        for t in (F(3,20),F(1,6),F(1,2),F(97,100),F(1)):
            weights=dict(b.weights_at(t));A,rhs=b.system(t)
            self.assertTrue(all(sum(row[j]*weights.get(p,0)for j,p in enumerate(b.P))==v for row,v in zip(A,rhs)))
            a=[b.m.s.Q(F(2,3),F(-7,5)),b.m.s.Q(-2,4),b.m.s.Q(F(1,9),F(1,2)),b.m.s.Q(t,-t)]
            target=b.m.s.leibniz(b.m.s.matrix_sum(a))
            self.assertEqual(target,sum((w*b.m.s.vertex(a,p)for p,w in weights.items()),b.m.s.Z))

    def test_forbidden_parameters_and_inexact_inputs(self):
        for t in (F(0),F(1,10),F(1,7),F(-1),F(2)):
            with self.assertRaises(ValueError):b.weights_at(t)
        for t in (0.5,True,'1/2'):
            with self.assertRaises(TypeError):b.weights_at(t)

    def test_changed_positive_mass_is_rejected(self):
        data=copy.deepcopy(self.positive);data['bases'][1]['numerators'][3][0]+=1
        with self.assertRaises(ValueError):b.verify_positive(data)

    def test_float_polynomial_is_rejected(self):
        data=copy.deepcopy(self.positive);data['bases'][0]['numerators'][0][0]=26244.0
        with self.assertRaises(ValueError):b.verify_positive(data)

    def test_dual_target_tampering_is_rejected(self):
        data=copy.deepcopy(self.dual);data['targetPolynomial'][0]+=1
        with self.assertRaises(ValueError):b.verify_dual(data)

    def test_missing_permutation_and_changed_normalization_rejected(self):
        data=copy.deepcopy(self.dual);data['columnPolynomials'].pop()
        with self.assertRaises(ValueError):b.verify_dual(data)
        data=copy.deepcopy(self.dual);data['commonFactor'][3]='-50176/5'
        with self.assertRaises(ValueError):b.verify_dual(data)

    def test_degree_bound_is_checked_before_interpolation(self):
        data=copy.deepcopy(self.dual);data['columnPolynomials'][0]=[0]*9+[1]
        with self.assertRaisesRegex(ValueError,'degree bound'):b.verify_dual(data)

    def test_zero_endpoint_needs_its_separate_dual(self):
        data=copy.deepcopy(self.dual);data['zeroDual']=['0']*17
        with self.assertRaises(ValueError):b.verify_dual(data)

    def test_exact_phase_one_discovery_matches_replayed_cases(self):
        zero=discover.phase(F(0));positive=discover.phase(F(3,20))
        self.assertEqual(F(zero['objective']),F(3472,33125))
        self.assertEqual(F(positive['objective']),0)
        self.assertEqual(positive['support'],self.positive['bases'][0]['support'])


if __name__=='__main__':unittest.main()
