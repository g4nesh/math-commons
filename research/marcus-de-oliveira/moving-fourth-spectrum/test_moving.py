import copy
from fractions import Fraction as F
from itertools import permutations
import json
from math import comb
from pathlib import Path
import unittest
import moving as m
import zero_fiber
import positive_box
import claude_star

HERE=Path(__file__).parent


class MovingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data=json.loads((HERE/'derived.json').read_text())

    def test_saved_exact_replay_and_base(self):
        report=m.verify(self.data)
        self.assertEqual(report,json.loads((HERE/'verification.json').read_text()))
        self.assertEqual(m.weights_at(self.data,0),m.BASE_WEIGHTS)
        self.assertEqual(report['realCoefficientIdentities'],17)

    def test_independent_cramer_derivation_reproduces_coefficients(self):
        self.assertEqual(m.derive(),self.data)

    def test_bernstein_conversion_and_strict_signs(self):
        left,right=m.INTERVAL
        for p in self.data['numerators']:
            controls=m.bernstein(p,left,right);n=len(p)-1
            self.assertGreater(min(controls),0)
            for j in range(n+1):
                x=F(j,n)
                reconstructed=sum(c*comb(n,i)*x**i*(1-x)**(n-i)for i,c in enumerate(controls))
                self.assertEqual(reconstructed,m.evaluate(p,left+(right-left)*x))

    def test_arbitrary_complex_triples_and_closed_endpoints(self):
        triples=((m.s.Q(1,2),m.s.Q(-3,4),m.s.Q(F(1,7),F(-2,9))),
                 (m.s.Q(-100),m.s.Q(0,11),m.s.Q(17,-21)))
        for t in m.INTERVAL+(F(0),F(1,3)):
            w=m.weights_at(self.data,t)
            self.assertGreater(min(w),0);self.assertEqual(sum(w),1)
            for triple in triples:
                a=list(triple)+[m.s.Q(1,-1)+t]
                target=m.s.leibniz(m.s.matrix_sum(a))
                mixed=sum((mass*m.s.vertex(a,p)for p,mass in zip(m.SUPPORT,w)),m.s.Z)
                self.assertEqual(target,mixed)
        self.assertNotEqual(m.weights_at(self.data,0),m.weights_at(self.data,1))

    def test_invalid_parameter_and_corrupted_coefficients(self):
        for t in (0.0,True,m.s.Q(0)):
            with self.assertRaises(TypeError):m.weights_at(self.data,t)
        for t in (F(-1),F(3)):
            with self.assertRaises(ValueError):m.weights_at(self.data,t)
        bad=copy.deepcopy(self.data);bad['numerators'][2][1]+=1
        with self.assertRaises(ArithmeticError):m.verify(bad)
        bad=copy.deepcopy(self.data);bad['numerators'][0][0]=float(bad['numerators'][0][0])
        with self.assertRaises(ValueError):m.verify(bad)

    def test_singular_support_is_not_automatically_inconsistent(self):
        statuses=m.verify(self.data)['specifiedSupportAtFormulaPoles']
        self.assertEqual([(x['coefficientRank'],x['augmentedRank'])for x in statuses],[(11,12),(11,11),(11,12)])
        self.assertEqual([x['specifiedSupportConsistent']for x in statuses],[False,True,False])

    def test_zero_fiber_published_separator(self):
        data=json.loads((HERE/'zero_fiber.json').read_text())
        self.assertEqual(zero_fiber.replay(data),{'triples':6,'permutations':24,
            'minimumGapAtZero':'1210344/15625','openSquareRadius':'26085/58260907'})

    def test_zero_fiber_by_direct_affine_coefficients(self):
        data=json.loads((HERE/'zero_fiber.json').read_text());pairings=list(permutations(range(4)))
        gaps={p:[F(0),F(0),F(0)]for p in pairings}
        for probe in data['probes']:
            a=[m.s.Q(bit)for bit in probe['triple']]+[m.s.Z]
            matrix=m.s.matrix_sum(a)
            target=(m.s.leibniz(matrix),m.s.leibniz([row[:3]for row in matrix[:3]]))
            points={p:(m.s.vertex(a,p),m.s.product(a[p[j]]+m.s.B_SPECTRUM[j]
                                               for j in range(4) if p[j]!=3)) for p in pairings}
            weights=[(tuple(item['p']),F(item['w']))for item in probe['pointwiseWeights']]
            self.assertEqual(tuple(sum((w*points[p][k]for p,w in weights),m.s.Z)for k in range(2)),target)
            x,y=F(probe['real']),F(probe['imag'])
            functional=lambda z:x*z.real+y*z.imag
            for p in pairings:
                constant=target[0]-points[p][0];slope=target[1]-points[p][1]
                for i,value in enumerate((constant,slope,m.s.I*slope)):gaps[p][i]+=functional(value)
        self.assertGreater(min(v[0]for v in gaps.values()),0)
        radius=min(c/(abs(dx)+abs(dy))for c,dx,dy in gaps.values()if dx or dy)
        self.assertEqual(radius,F(26085,58260907))
        # A strictly smaller closed square has a positive verified margin at
        # every corner, checking the open-square quantifier and signs.
        for c,dx,dy in gaps.values():
            for x in (-radius/2,radius/2):
                for y in (-radius/2,radius/2):self.assertGreater(c+dx*x+dy*y,0)

    def test_zero_fiber_rejects_corrupted_pointwise_mass(self):
        data=json.loads((HERE/'zero_fiber.json').read_text());data['probes'][0]['pointwiseWeights'][0]['w']='0'
        with self.assertRaises(ValueError):zero_fiber.replay(data)

    def test_closed_positive_box_via_different_coefficient_basis(self):
        report=positive_box.verify()
        self.assertEqual(report['radius'],'1/1024')
        self.assertEqual(report['maxRelativeWeightChange'],'2753800000/592463385761')
        radius=F(1,1024)
        for x,y in ((-radius,-radius),(-radius,radius),(radius,-radius),(radius,radius),(F(0),F(0))):
            result=positive_box.weights_at(x,y)
            self.assertEqual([tuple(row['columnToRow'])for row in result],list(m.SUPPORT))
            weights=[F(row['weight'])for row in result]
            self.assertGreater(min(weights),0)
            rows,rhs=m.coefficient_system(m.s.Q(1+x,-1+y))
            for row,target in zip(rows,rhs):
                self.assertEqual(sum(w*c for w,c in zip(weights,row)),target)

    def test_positive_box_rejects_outside_displacement(self):
        with self.assertRaises(ValueError):positive_box.weights_at(F(1,512),F(0))

    def test_claude_inequality_is_only_a_fixture_check(self):
        self.assertEqual(claude_star.check(),{'relabellings':1152,'fixtureSlack':'5779/15625','minimumSlack':'15136/390625'})


if __name__=='__main__':unittest.main()
