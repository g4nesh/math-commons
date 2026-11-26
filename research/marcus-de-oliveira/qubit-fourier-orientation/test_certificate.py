import copy
import json
from fractions import Fraction as F
from pathlib import Path
import unittest
from model import FACETS,ORBIT_REPS,ROWS,orbit_column,require,solve,target
from verify import construct,verify
from fixtures import replay


class CertificateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data=json.loads(Path(__file__).with_name('certificate.json').read_text())
        cls.laws=verify(cls.data)

    def test_complete_vertex_count_and_exact_complex_bridge(self):
        self.assertEqual(len(self.laws),66)
        self.assertNotEqual(replay()[1]['triple_product'],'0')

    def test_missing_and_duplicate_vertex_rejected(self):
        bad=copy.deepcopy(self.data);bad['vertices'].pop()
        with self.assertRaises(ValueError):verify(bad)
        bad=copy.deepcopy(self.data);bad['vertices'].append(bad['vertices'][0])
        with self.assertRaises(ValueError):verify(bad)

    def test_float_and_negative_mass_rejected(self):
        bad=copy.deepcopy(self.data);bad['vertices'][0]['weights'][0]=0.0
        with self.assertRaises(ValueError):verify(bad)
        bad=copy.deepcopy(self.data);bad['vertices'][0]['weights'][0]='-1'
        with self.assertRaises(ValueError):verify(bad)

    def test_wrong_nonnegative_law_rejected(self):
        bad=copy.deepcopy(self.data);bad['vertices'][0]['weights']=['1']+['0']*29
        with self.assertRaises(ValueError):verify(bad)

    def test_construct_boundary_and_interior(self):
        columns=[orbit_column(c) for c in ORBIT_REPS]
        for point in [(1,1,1,0),(0,0,0,F(7,4)),(F(1,5),F(2,7),F(-1,8),F(1,2))]:
            weights=construct(point,self.laws)
            self.assertEqual(sum(weights),1);self.assertTrue(all(w>=0 for w in weights))
            self.assertEqual(tuple(sum(c[j]*w for c,w in zip(columns,weights)) for j in range(len(ROWS))),target(point))
        for point in [(2,0,0,0),(0.0,0,0,0),(True,0,0,0)]:
            with self.assertRaises(ValueError):construct(point,self.laws)

    def test_degenerate_linear_system_and_area_identity(self):
        self.assertIsNone(solve([(1,0),(2,0)],(1,2)))
        # Independent polynomial coefficient replay of the three-variable SOS.
        def mul(a,b):
            out={}
            for x,c in a.items():
                for y,d in b.items():
                    key=tuple(i+j for i,j in zip(x,y));out[key]=out.get(key,F(0))+c*d
            return out
        def square(a):return mul(a,a)
        left=square({(1,0,0):F(1),(0,1,0):F(1),(0,0,1):F(-1)})
        left[(1,1,0)]-=3;left[(0,0,2)]+=3
        right=square({(1,0,0):F(1,2),(0,1,0):F(1,2),(0,0,1):F(-2)})
        for k,v in square({(1,0,0):F(1),(0,1,0):F(-1)}).items():right[k]=right.get(k,F(0))+F(3,4)*v
        self.assertEqual(left,right)
        with self.assertRaises(ValueError):require(False,'must fail under -O')


if __name__=='__main__':unittest.main()
