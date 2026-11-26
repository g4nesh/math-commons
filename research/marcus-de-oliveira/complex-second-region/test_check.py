import copy
from fractions import Fraction as F
import importlib.util
from pathlib import Path
import unittest

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('complex_region_check',HERE/'check.py')
c=importlib.util.module_from_spec(spec);spec.loader.exec_module(c)


class RectangleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.data=c.prepared()

    def test_exact_proof_and_report(self):
        self.assertEqual(c.report(),c.load('report.json'))
        self.assertEqual(self.data['numeratorDisplacement'],F(785,989))

    def test_arbitrary_fourth_and_free_entries(self):
        r=self.data['height']
        for u,v,s,t in ((c.LOWER,-r,100000,-200000),(c.UPPER,r,F(13,7),F(-19,11))):
            weights=c.weights_at(u,v,s,t)
            self.assertEqual(sum(weights),1)
            self.assertTrue(all(weights[p]>=self.data['massFloor'] for p in c.ALLOWED))
            self.assertTrue(all(weights[p]==0 for p in range(24) if p not in c.ALLOWED))
            for x,y in ((c.Z,c.O),(c.Q(17,-3),c.Q(-23,31))):
                c.g.prior.direct_check([x,c.Q(u,v),y,c.Q(s,t)],weights)

    def test_closed_boundaries_and_outside(self):
        r=self.data['height']
        for u,v in ((c.LOWER,r),(c.UPPER,-r)):c.weights_at(u,v,0,0)
        for u,v in ((c.LOWER-F(1,10**9),0),(c.UPPER+F(1,10**9),0),(F(1,2),r+F(1,10**9))):
            with self.assertRaises(ValueError):c.weights_at(u,v,0,0)

    def test_exact_input_contract(self):
        for bad in (True,0.0,'0',complex(0)):
            with self.assertRaises(TypeError):c.weights_at(F(1,2),bad,0,0)
            with self.assertRaises(TypeError):c.weights_at(F(1,2),0,bad,0)

    def test_center_floor_corruption(self):
        law=copy.deepcopy(c.load('real-law.json'));law['terms'][0][3]='-1'
        with self.assertRaises(ValueError):c.center_controls(law)

    def test_polynomial_relation_corruption(self):
        rank=copy.deepcopy(c.load('rank.json'))
        for term in rank['relations'][-1]:
            row=term[1]
            if any(self.data['a0'][row]) or any(self.data['ar'][row]) or any(self.data['ai'][row]):
                term[2]=str(F(term[2])+1);break
        with self.assertRaises(ArithmeticError):
            c.verify_relations(rank,self.data['a0'],self.data['ar'],self.data['ai'])


if __name__=='__main__':unittest.main()
