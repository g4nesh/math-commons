"""Exact regressions for the quantitative two-free-eigenvalue certificate."""
from fractions import Fraction as F
from pathlib import Path
import copy
import json
import unittest
import radius as r

HERE=Path(__file__).resolve().parent


class QuantitativeRadiusTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.bound=json.loads((HERE/'bound.json').read_text())
        cls.radius=F(cls.bound['fourthRealInterval'][1])
        cls.point=(cls.radius,-cls.radius)
        cls.values=r.coefficients_at(*cls.point)

    def test_complete_exact_replay_and_recorded_output(self):
        report=r.verify(self.bound)
        self.assertEqual(report,json.loads((HERE/'verification.json').read_text()))
        self.assertEqual(self.radius,F(1,32768))
        self.assertLess(F(self.bound['perturbationNormBound']),F(1,5))
        self.assertLess(F(self.bound['actionNormBound']),F(3,100))
        self.assertLess(F(self.bound['relativeSolutionChangeBound']),F(3,80))
        self.assertGreater(F(self.bound['minimumProductThroughoutSquare']),F(1,4000))
        self.assertEqual(len(report['cornerCertificates']),5)

    def test_offcenter_quadratic_identity_all_coefficients(self):
        u,v=self.point
        rows=r.c.gap_rows(fourth=r.c.s.Q(u,v))
        lhs=[F(0)]*36
        for weight,(i,j)in zip(self.values,r.PAIRS):
            lhs=[a+weight*b for a,b in zip(lhs,r.c.quadratic_product(rows[i],rows[j]))]
        real,imag=r.c.determinant_parts()
        self.assertEqual(lhs,[self.values[35]*a+self.values[36]*b for a,b in zip(real,imag)])
        self.assertEqual(sum(self.values[:35]),1)
        self.assertGreater(min(self.values[:35]),F(1,4000))
        self.assertNotEqual(tuple(self.values),r.BASE)

    def test_rankone_specialization_matches_direct_determinants(self):
        s=r.c.s;u,v=self.point;z=s.Q(u,v)
        rows=r.c.gap_rows(fourth=z)
        for x,y,phase in ((s.Q(2,3),s.Q(-4,1),s.Q(3,-2)),(s.Z,s.Q(F(7,3),F(-2,5)),s.Q(-1,4))):
            complex_coordinates=[phase,phase*x,phase*y,phase*x*y]
            coordinates=[component for a in complex_coordinates for component in(a.real,a.imag)]
            diagonal=[x,y,s.Q(2,F(-1,4)),z]
            target=s.e.determinant(s.matrix_sum(diagonal))
            gaps=[(phase*(s.vertex(diagonal,p)-target)).real for p in r.c.PERMS]
            self.assertEqual(gaps,[sum(a*b for a,b in zip(row,coordinates))for row in rows])
            self.assertEqual(complex_coordinates[0]*complex_coordinates[3]-complex_coordinates[1]*complex_coordinates[2],s.Z)
            self.assertEqual(sum(w*gaps[i]*gaps[j]for w,(i,j)in zip(self.values,r.PAIRS)),0)

    def test_separator_margin_for_every_permutation(self):
        base,dx,dy=r.gap_data();separator=list(map(F,r.SOURCE['strictSeparator']))
        lower=[]
        for row,a,b in zip(base,dx,dy):
            value=sum(x*y for x,y in zip(separator,row))
            slope=abs(sum(x*y for x,y in zip(separator,a)))+abs(sum(x*y for x,y in zip(separator,b)))
            lower.append(value-self.radius*slope)
        self.assertEqual(len(lower),24)
        self.assertEqual(min(lower),F(27177805772189,25622138880000000))
        self.assertGreaterEqual(min(lower),F(self.bound['uniformLawMarginThroughoutSquare']))
        self.assertGreater(min(lower),0)

    def test_closed_endpoints_and_exact_input_validation(self):
        self.assertEqual(r.exact_parameters(self.radius,-self.radius,self.radius),self.point)
        for bad in (0.0,True,False,'0',complex(0)):
            with self.assertRaises(TypeError):r.coefficients_at(bad,F(0))
        for u,v in ((2*self.radius,0),(0,-2*self.radius)):
            with self.assertRaises(ValueError):r.coefficients_at(u,v)

    def test_modified_bound_is_rejected(self):
        altered=copy.deepcopy(self.bound)
        altered['fourthRealInterval'][1]='1/16384'
        with self.assertRaisesRegex(ValueError,'differs from exact derivation'):
            r.verify(altered)

    def test_exact_solver_residual_and_singular_rejection(self):
        a=[[F(2),F(1)],[F(1),F(3)]]
        b=[[F(1),F(0)],[F(0),F(1)]]
        inverse=r.solve_many(a,b)
        self.assertEqual(r.matmul(inverse,a),b)
        with self.assertRaisesRegex(ValueError,'Singular'):
            r.solve_many([[F(1),F(2)],[F(2),F(4)]],[[F(1)],[F(2)]])


if __name__=='__main__':unittest.main()
