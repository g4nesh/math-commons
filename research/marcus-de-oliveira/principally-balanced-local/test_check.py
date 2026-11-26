import unittest
from fractions import Fraction as F
from itertools import permutations
import check as c


class ExactLawTests(unittest.TestCase):
    def test_quadratic_field_order_and_inverse(self):
        for a,b,sign in [(0,0,0),(1,1,1),(-1,-1,-1),(1,-1,-1),(-1,1,1),(2,-1,1),(-2,1,-1)]:
            self.assertEqual(c.Q2(a,b).sign(),sign)
        self.assertEqual(c.Q2(3,2)*c.Q2(3,-2),c.ONE)
        self.assertEqual(c.Q2(3,2)/c.Q2(3,2),c.ONE)
        with self.assertRaises(ZeroDivisionError):
            c.ONE/c.ZERO
        for value in (0.5,True,'1'):
            with self.assertRaises(TypeError):
                c.Q2(value)

    def test_invalid_data_is_rejected(self):
        data={():c.ONE,(0,):c.Q2(F(1,2)),(1,):c.Q2(F(1,2)),(0,1):c.ZERO}
        with self.assertRaises(ValueError):
            c.construct([0,0],data)
        with self.assertRaises(ValueError):
            c.construct([0,1],{():c.ONE})
        bad=dict(data);bad[(0,)]=c.ONE
        with self.assertRaises(ValueError):
            c.construct([0,1],bad)
        outside=dict(data);outside[(0,)]=c.Q2(2);outside[(1,)]=c.Q2(-1)
        with self.assertRaisesRegex(ValueError,'positivity'):
            c.construct([0,1],outside)

    def test_exact_off_center_report(self):
        result=c.report()
        self.assertEqual(result['correctionL1'],['0','2999997000/1000002000001'])
        self.assertEqual(len(result['weights']),24)

    def test_direct_complex_determinant(self):
        result=c.report()
        spectrum=[c.Q2(F(a),F(b)) for a,b in result['spectrum']]
        weights=[c.Q2(F(row['weight'][0]),F(row['weight'][1])) for row in result['weights']]
        pairings=[row['permutation'] for row in result['weights']]
        skew=[[F(x) for x in row] for row in result['skewPartK']]
        zero=(c.ZERO,c.ZERO)
        one=(c.ONE,c.ZERO)
        def add(x,y):return (x[0]+y[0],x[1]+y[1])
        def mul(x,y):return (x[0]*y[0]-x[1]*y[1],x[0]*y[1]+x[1]*y[0])
        def scale(w,x):return (w*x[0],w*x[1])
        for a in [[(c.Q2(2),c.Q2(-3)),(c.Q2(0),c.Q2(7)),(c.Q2(F(-2,3)),c.Q2(F(5,7))),(c.Q2(100),c.Q2(-100))], [zero]*4]:
            matrix=[[(a[i] if i==j else (c.ZERO,c.Q2(skew[i][j]))) for j in range(4)] for i in range(4)]
            actual=zero
            for p in permutations(range(4)):
                term=one
                for i in range(4):term=mul(term,matrix[i][p[i]])
                sign=(-1)**sum(p[i]>p[j] for i in range(4) for j in range(i+1,4))
                actual=add(actual,scale(sign,term))
            mixture=zero
            for p,w in zip(pairings,weights):
                term=one
                for i in range(4):term=mul(term,add(a[i],(spectrum[p[i]],c.ZERO)))
                mixture=add(mixture,scale(w,term))
            self.assertEqual(actual,mixture)


if __name__=='__main__':
    unittest.main()
