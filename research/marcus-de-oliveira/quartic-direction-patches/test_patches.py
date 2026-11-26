"""Proof replay and direct physical regressions for quartic angular patches."""
import copy
from fractions import Fraction as F
import json
from pathlib import Path
import unittest
import patches as p

HERE=Path(__file__).resolve().parent


class QuarticPatchTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.certificate=json.loads((HERE/'certificate.json').read_text())
        cls.patches=cls.certificate['patches']

    def test_full_exact_certificate_and_saved_report(self):
        self.assertEqual(p.verify(self.certificate),json.loads((HERE/'verification.json').read_text()))

    def test_direct_determinants_at_every_patch_endpoint(self):
        s=p.s
        for patch in self.patches:
            for k,t in ((F(patch['interval'][0]),F(-2)),(F(patch['center']),F(0)),(F(patch['interval'][1]),F(3,7))):
                weights=p.weights_at(patch,k,t)
                omega=s.Q(1,k)if patch['chart']=='real'else s.Q(k,1)
                diagonal=[s.Q(F(1,3),2),omega*s.Q(t),s.Q(-3,F(1,2)),s.Z]
                target=s.leibniz(s.matrix_sum(diagonal))
                reconstructed=sum((s.vertex(diagonal,perm)*s.Q(w)for perm,w in zip(p.PERMS,weights)),s.Z)
                self.assertEqual(target,reconstructed)
                self.assertEqual(sum(weights),1)
                self.assertGreaterEqual(min(weights),0)

    def test_positive_normalizer_at_all_primitive_roots(self):
        for patch in self.patches:
            k=F(patch['center'])
            for t in (-2,-1,0,1,2):
                self.assertEqual(sum(p.weights_at(patch,k,t)),1)

    def test_bernstein_conversion_exact(self):
        # 1+2x+3x² on [0,1] has Bernstein coordinates 1,2,6.
        self.assertEqual(p.bernstein([F(1),F(2),F(3)],F(0),F(1)),[F(1),F(2),F(6)])
        self.assertEqual(p.bernstein([F(0),F(0),F(1)],F(-1),F(1)),[F(1),F(-1),F(1)])

    def test_altered_identity_and_duplicate_support_rejected(self):
        altered=copy.deepcopy(self.patches[0])
        altered['terms'][0]['numerator'][0]=str(F(altered['terms'][0]['numerator'][0])+1)
        with self.assertRaisesRegex(ValueError,'Polynomial identity'):
            p.check_patch(altered)
        repeated=copy.deepcopy(self.patches[0]);repeated['terms'].append(repeated['terms'][0])
        with self.assertRaisesRegex(ValueError,'repeated support'):
            p.check_patch(repeated)

    def test_negative_sign_and_bad_chart_rejected(self):
        altered=copy.deepcopy(self.patches[0])
        altered['terms'][0]['numerator']=[str(-F(v))for v in altered['terms'][0]['numerator']]
        with self.assertRaisesRegex(ValueError,'Nonpositive Bernstein'):
            p.check_patch(altered)
        altered=copy.deepcopy(self.patches[0]);altered['chart']='unclassified'
        with self.assertRaisesRegex(ValueError,'Unknown angular chart'):
            p.check_patch(altered)

    def test_missing_stated_patch_rejected(self):
        altered=copy.deepcopy(self.certificate);altered['patches'].pop()
        with self.assertRaisesRegex(ValueError,'both be present'):
            p.verify(altered)

    def test_exact_inputs_and_closed_interval(self):
        patch=self.patches[0]
        for k in map(F,patch['interval']):self.assertEqual(sum(p.weights_at(patch,k,F(1))),1)
        for value in (True,0.5,'1/2',complex(1)):
            with self.assertRaises(TypeError):p.weights_at(patch,value,F(0))
        with self.assertRaises(ValueError):p.weights_at(patch,F(patch['interval'][0])-1,F(0))


if __name__=='__main__':unittest.main()
