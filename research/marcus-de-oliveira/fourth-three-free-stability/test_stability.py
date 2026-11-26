from fractions import Fraction as F
from pathlib import Path
import copy,json,sys,unittest
import check as proof
import neighborhood as n
c=n.c
sys.path.insert(0,str(c.HERE.parent/'four-cycle-spectral'))
import spectral as s


class StabilityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.certificate=json.loads((c.HERE/'certificate.json').read_text())
        cls.patches=json.loads(proof.SOURCE.read_text())

    def test_complete_exact_replay(self):
        actual=proof.verify(self.certificate)
        self.assertEqual(actual,json.loads((c.HERE/'verification.json').read_text()))
        self.assertEqual(actual['commonRadius'],'1/4194304')
        self.assertEqual(actual['invariantKernelDimensions'],[46]*6)

    def test_original_fixed_B_determinants_at_spectral_and_angular_endpoints(self):
        for item in self.certificate['patches']:
            index=item['index'];r=F(item['radius']);patch=self.patches[index]
            for real,imag,slope,t in ((r,-r,F(patch['interval'][0]),F(0)),(-r,r,F(patch['interval'][1]),F(3,7))):
                weights=proof.weights_at(index,real,imag,slope,t)
                delta=s.Q(real,imag);omega=s.Q(1,slope)if patch['chart']=='h'else s.Q(slope,1)
                diagonal=[delta+s.Q(F(1,3),2),delta+omega*s.Q(t),delta+s.Q(-3,F(1,2)),delta]
                target=s.e.determinant(s.matrix_sum(diagonal))
                mixture=sum((s.vertex(diagonal,p)*s.Q(w)for p,w in zip(c.a.PERMS,weights)),s.Z)
                self.assertEqual(target,mixture)
                self.assertEqual(sum(weights),1)
                self.assertGreaterEqual(min(weights),0)

    def test_quadratic_matrix_extraction_at_independent_complex_parameter(self):
        patch=self.patches[0];keys,base,matrices=n.data(patch);u,v=F(1,3),F(-2,7)
        monomials=(u,v,u*u,u*v,v*v)
        columns=[c.column(patch,t,(u,v))for t in patch['terms']]
        expected=[[col.get(key,F(0))for col in columns]for key in keys]
        actual=[[base[i][j]+sum(factor*matrix[i][j]for factor,matrix in zip(monomials,matrices))for j in range(64)]for i in range(len(keys))]
        self.assertEqual(actual,expected)
        self.assertTrue(all(gap.get((0,0,0),c.a.ZERO)==c.a.ZERO for gap in c.gaps(u,v)))

    def test_source_tamper_and_missing_sector_rejected(self):
        changed=copy.deepcopy(self.certificate);changed['sourceSha256']='0'*64
        with self.assertRaisesRegex(ValueError,'source changed'):proof.verify(changed)
        changed=copy.deepcopy(self.certificate);changed['patches'].pop()
        with self.assertRaisesRegex(ValueError,'six quadratic'):proof.verify(changed)

    def test_exact_input_and_domain_guards(self):
        for value in (0.0,True,'0',complex(0)):
            with self.assertRaises(TypeError):proof.coefficients_at(0,value,F(0))
        with self.assertRaises(ValueError):proof.coefficients_at(3,F(0),F(0))
        with self.assertRaises(ValueError):proof.coefficients_at(0,F(1),F(0))
        with self.assertRaises(ValueError):proof.weights_at(0,F(0),F(0),F(1),F(0))

    def test_fraction_free_pivoted_multiple_rhs(self):
        matrix=[[F(0),F(1,3),F(2)],[F(2),F(-1),F(1,7)],[F(1,2),F(2),F(-3)]]
        rhs=[[F(1),F(2)],[F(-2,3),F(1,7)],[F(4),F(0)]]
        self.assertEqual(c.matmul(matrix,c.solve_many(matrix,rhs)),rhs)
        with self.assertRaisesRegex(ValueError,'Singular'):
            c.solve_many([[F(1),F(2)],[F(2),F(4)]],[[F(1)],[F(2)]])


if __name__=='__main__':unittest.main()
