import unittest
import sympy as s
from check import equal, hermitian_matrix, require


class ExactContract(unittest.TestCase):
    def test_rational_identity_and_changed_sign(self):
        x=s.symbols('x')
        equal((x*x-1)/(x-1),x+1,'factor')
        with self.assertRaises(ValueError): equal((x*x-1)/(x-1),x-1,'altered sign')

    def test_hermitian_complex_evaluation(self):
        v,V=s.symbols('v V')
        M=hermitian_matrix(2+3*(v+V)+7*v*V,v,V)
        self.assertEqual(M,s.Matrix([[2,3],[3,7]]))
        z=1+2*s.I
        equal((s.Matrix([1,s.conjugate(z)]).T*M*s.Matrix([1,z]))[0],43,'complex value')

    def test_reject_nonhermitian_and_higher_degree(self):
        v,V=s.symbols('v V')
        for p in [v+2*V,v*v+V*V,1/(1+v*V)]:
            with self.assertRaises(ValueError): hermitian_matrix(p,v,V)

    def test_semidefinite_boundary_used_in_proof(self):
        M=s.Matrix([[-1,-2],[-2,-4]])
        self.assertEqual(M.det(),0)
        self.assertTrue(all(value<=0 for value in M.eigenvals()))
        with self.assertRaises(ValueError): require(False,'guards survive optimization')


if __name__=='__main__': unittest.main()
