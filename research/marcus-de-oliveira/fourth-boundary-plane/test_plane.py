import copy,json,unittest
from fractions import Fraction as F
import plane as p

class PlaneTests(unittest.TestCase):
    def test_complete_exact_replay_matches_report(self):
        self.assertEqual(p.verify(),json.loads((p.HERE/'verification.json').read_text()))

    def test_nonconstant_weights_reconstruct_complex_determinants(self):
        all_weights=[]
        for s,t in ((F(31,200),F(-1,50)),(F(3,20),F(0)),(F(31,200),F(1,50))):
            weights=p.weights_at(s,t);all_weights.append(weights)
            a=[p.m.s.Q(F(-2,3),F(7,5)),p.m.s.Q(3,-2),p.m.s.Q(F(1,9),F(-1,7)),p.m.s.Q(s+t,t-s)]
            target=p.m.s.leibniz(p.m.s.matrix_sum(a))
            self.assertEqual(target,p.m.s.e.determinant(p.m.s.matrix_sum(a)))
            self.assertEqual(target,sum((w*p.m.s.vertex(a,perm)for perm,w in weights),p.m.s.Z))
        self.assertNotEqual(all_weights[0],all_weights[-1])

    def test_excluded_side_has_a_real_separating_functional(self):
        data=p.load_data();H=list(map(p.decode,data['dualColumns']));critical=p.decode(data['numerators'][1])
        for t in p.TI+(F(0),):
            self.assertLess(p.evaluate(critical,p.SI[0],t),0)
            self.assertTrue(all(p.evaluate(h,p.SI[0],t)>=0 for h in H))
            with self.assertRaises(ValueError):p.weights_at(p.SI[0],t)

    def test_boundary_brackets_cover_transverse_endpoints(self):
        for t in p.TI+(F(0),):
            lo,hi=p.boundary_bracket(t,40)
            self.assertLess(p.evaluate(p.Q,lo,t),0)
            self.assertGreaterEqual(p.evaluate(p.Q,hi,t),0)
            self.assertEqual(hi-lo,(p.SI[1]-p.SI[0])/2**40)

    def test_modified_critical_polynomial_is_rejected(self):
        data=p.load_data();data['numerators'][1][0][2]=str(F(data['numerators'][1][0][2])+1)
        with self.assertRaises(ValueError):p.verify(data)

    def test_modified_noncritical_mass_and_missing_dual_rejected(self):
        data=p.load_data();data['numerators'][3][0][2]=str(F(data['numerators'][3][0][2])+1)
        with self.assertRaises(ValueError):p.verify(data)
        data=p.load_data();data['dualColumns'].pop()
        with self.assertRaises(ValueError):p.verify(data)

    def test_inexact_inputs_domains_and_degree_rejected(self):
        for s,t in ((0.15,F(0)),(True,F(0))):
            with self.assertRaises(TypeError):p.weights_at(s,t)
        for s,t in ((F(1),F(0)),(F(3,20),F(1))):
            with self.assertRaises(ValueError):p.weights_at(s,t)
        for poly in ([[13,0,'1']],[[0,0,1.0]],[[0,0,'1'],[0,0,'2']]):
            with self.assertRaises(ValueError):p.decode(poly)

if __name__=='__main__':unittest.main()
