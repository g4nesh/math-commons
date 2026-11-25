import json
from fractions import Fraction as F
from pathlib import Path
import unittest
import spectral as s
import fiber_certificate as fiber
from itertools import combinations, permutations


class SpectralTests(unittest.TestCase):
    def test_exact_coefficient_identities(self):
        target,points,area,nums,den = s.interval_algebra()
        self.assertEqual(area,(F(3000),F(-2000),F(500)))
        self.assertEqual(nums,s.NUMERATORS)
        self.assertEqual(den,s.DENOMINATOR)
        self.assertEqual(target[0],s.Q(F(269802,390625),F(-5707764,390625)))

    def test_closed_endpoints_are_positive_and_reconstruct(self):
        for t in s.INTERVAL:
            weights = s.interval_weights(t)
            self.assertGreater(min(weights),0)
            self.assertEqual(sum(weights),1)
            value = s.check_at(s.BASE_A[0]+s.I*t,s.TRIANGLE,weights)
            target = s.affine_data(s.BASE_A[0],s.I,s.TRIANGLE)[0]
            self.assertEqual(value,s.evaluate(target,t))

    def test_weights_match_separate_real_linear_system(self):
        for t in (F(-40),F(-13,7),F(0),F(1)):
            target,points = s.affine_data(s.BASE_A[0],s.I,s.TRIANGLE)
            vertices = [s.evaluate(p,t) for p in points]
            value = s.evaluate(target,t)
            solved = s.exact_solve([[1]*3,[z.real for z in vertices],[z.imag for z in vertices]],
                                   [1,value.real,value.imag])
            self.assertEqual(solved,s.interval_weights(t))

    def test_continuous_sign_proof_and_isolating_brackets(self):
        proof = s.positivity()
        self.assertEqual(proof['n2Discriminant'],'-244316944388864')
        self.assertEqual(proof['n1ConcaveEndpointValues'],['3403196','8044970'])
        self.assertLess(s.evaluate(s.NUMERATORS[0],F(-41)),0)
        self.assertLess(s.evaluate(s.NUMERATORS[2],F(3,2)),0)
        for t in (F(-100),F(0),F(2),F(100)):
            self.assertEqual(s.evaluate(s.DENOMINATOR,t),7812500*((t-2)**2+2))

    def test_interval_parameter_contract(self):
        for value in (F(-41),F(3,2)):
            with self.assertRaises(ValueError):s.interval_weights(value)
        for value in (0.0,True,s.Q(0)):
            with self.assertRaises(TypeError):s.interval_weights(value)
        self.assertNotEqual(s.interval_weights(0),s.interval_weights(1))

    def test_all_complex_plane_coefficients_and_off_grid_values(self):
        target,points,weights = s.plane_algebra()
        self.assertGreater(min(weights),0)
        for z in (s.Q(F(1,7),F(2,9)),s.Q(-100,13),s.Q(2),s.Q(0,-3)):
            self.assertEqual(s.check_at(z,s.PLANE,weights),s.evaluate(target,z))

    def test_normality_and_noncommuting_family(self):
        b = s.e.multiply(s.e.multiply(s.U,s.e.diagonal(s.B_SPECTRUM)),s.e.adjoint(s.U))
        self.assertEqual(s.e.multiply(b,s.e.adjoint(b)),s.e.multiply(s.e.adjoint(b),b))
        # The (0,1) commutator entry is (-1+i(t-1))*B01; both factors are
        # nonzero for every real t, so the whole family is noncommuting.
        self.assertNotEqual(b[0][1],s.Z)
        for t in s.INTERVAL+(F(0),):
            a=s.BASE_A.copy();a[0]+=s.I*t
            da=s.e.diagonal(a)
            self.assertEqual(s.e.multiply(da,s.e.adjoint(da)),s.e.multiply(s.e.adjoint(da),da))
            self.assertNotEqual(s.e.multiply(da,b),s.e.multiply(b,da))

    def test_rejects_a_weight_perturbation(self):
        w=list(s.interval_weights(0));w[0]+=F(1,100);w[1]-=F(1,100)
        with self.assertRaises(AssertionError):s.check_at(s.BASE_A[0],s.TRIANGLE,w)

    def test_saved_artifact_and_source_identity(self):
        s.source_check()
        saved=json.loads((Path(__file__).parent/'certificate.json').read_text())
        self.assertEqual(saved,s.report())
        self.assertNotIn('pythonVersion',saved)

    def test_three_variable_fiber_by_direct_principal_coefficients(self):
        certificate=json.loads((Path(__file__).parent/'fiber_certificate.json').read_text())
        terms=[(tuple(r['columnToRow']),F(r['mass'])) for r in certificate['weightsAtFixedFourthEigenvalue']]
        base=s.matrix_sum([s.Z,s.Z,s.Z,s.BASE_A[3]])
        for degree in range(4):
            for chosen in combinations(range(3),degree):
                # Independent of Boolean-corner interpolation: compare each
                # monomial coefficient via the complementary principal minor.
                remaining=[i for i in range(4) if i not in chosen]
                minor=[[base[i][j] for j in remaining] for i in remaining]
                actual=s.leibniz(minor)
                expected=sum((w*s.product(s.B_SPECTRUM[j]+(s.BASE_A[3] if p[j]==3 else s.Z)
                                          for j in range(4) if p[j] not in chosen)
                              for p,w in terms),s.Z)
                self.assertEqual(actual,expected)

    def test_fiber_separator_and_no_single_diagonal_counterexample(self):
        certificate=json.loads((Path(__file__).parent/'fiber_certificate.json').read_text())
        result=fiber.check_certificate(certificate)
        self.assertEqual(result['strictGap'],'11201152/15625')
        self.assertEqual(result['maximumPairingScore'],'0')
        self.assertEqual(result['positiveWeights'],12)
        for term in certificate['separatingFunctional']:
            a=[s.Q(bit) for bit in term['diagonal']]
            value=s.leibniz(s.matrix_sum(a))
            points=[s.vertex(a,p).pair() for p in permutations(range(4))]
            hull=s.e.certify_hull(points,value.pair())
            self.assertEqual(hull['membership'],'inside')
            self.assertTrue(s.e.check_hull_certificate(points,value.pair(),hull))

    def test_fiber_checker_rejects_corrupt_mass_and_empty_separator(self):
        path=Path(__file__).parent/'fiber_certificate.json'
        certificate=json.loads(path.read_text())
        certificate['weightsAtFixedFourthEigenvalue'][0]['mass']='0'
        with self.assertRaises(ValueError):fiber.check_certificate(certificate)
        certificate=json.loads(path.read_text())
        certificate['separatingFunctional']=[]
        with self.assertRaises(ValueError):fiber.check_certificate(certificate)


if __name__=='__main__':
    unittest.main()
