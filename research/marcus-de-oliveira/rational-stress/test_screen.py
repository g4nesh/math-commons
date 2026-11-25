import copy
from fractions import Fraction as F
from pathlib import Path
import sys
import unittest

HERE = Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
import screen as S
import replay as R
import holdout as H


class ScreenTests(unittest.TestCase):
    def test_declared_schedule(self):
        specs = list(S.case_specs())
        self.assertEqual(len(specs),900)
        self.assertEqual(len({spec['id'] for spec in specs}),900)
        self.assertEqual([sum(spec['n'] == n for spec in specs) for n in (4,5,6)],[360,360,180])

    def test_complex_givens_and_all_spectrum_families(self):
        # Tightest angle and largest dimension exercise the file parser's bounds.
        for family in S.FAMILIES:
            for mixing in S.MIXING:
                spec = dict(id='test-'+family+'-'+mixing,n=6,family=family,mixing=mixing,repetition=2,seed=12)
                witness = S.make_witness(spec)
                row = S.analyze(witness)
                self.assertTrue(row['noncommuting'])
                self.assertTrue(row['complexUnitary'])
                self.assertTrue(R.audit(witness,row))
                self.assertEqual(row['membership'],'inside')
                if family == 'singular': self.assertTrue(row['singularSum'])
                if family == 'hermitian': self.assertLessEqual(row['hullDimension'],1)

    def test_zheng_functional_exact(self):
        rows = ((1,4,4,4),(-4,1,4,-4),(-4,-4,1,4),(-4,4,-4,1))
        U = [[S.Q(F(value,7)) for value in row] for row in rows]
        self.assertEqual(S.V.multiply(U,S.V.adjoint(U)),S.V.identity(4))
        self.assertEqual(S.exotic_psi(U),F(-158,2401))
        for p in __import__('itertools').permutations(range(4)):
            permutation = [[S.Q(int(p[i] == j)) for j in range(4)] for i in range(4)]
            self.assertGreaterEqual(S.exotic_psi(permutation),0)

    def test_independent_replay_detects_tampered_certificate(self):
        witness = S.make_witness(next(S.case_specs()))
        row = S.analyze(witness)
        broken = copy.deepcopy(row)
        broken['convexCombination'][0]['weight'] = str(F(broken['convexCombination'][0]['weight'])+1)
        with self.assertRaisesRegex(ArithmeticError,'Weight sum'): R.audit(witness,broken)
        broken = copy.deepcopy(row); broken['determinant'][0] = str(F(broken['determinant'][0])+1)
        with self.assertRaisesRegex(ArithmeticError,'Determinant mismatch'): R.audit(witness,broken)

    def test_independent_determinant_row_swap_and_singular(self):
        self.assertEqual(R.determinant([[R.ZERO,R.ONE],[R.ONE,R.ZERO]]),(F(-1),F(0)))
        self.assertEqual(R.determinant([[R.ONE,R.ONE],[R.ONE,R.ONE]]),R.ZERO)

    def test_n7_holdout_is_separate_and_independently_checked(self):
        self.assertEqual(len(list(H.specs())),24)
        witness = S.make_witness(next(H.specs()))
        with self.assertRaises(S.V.InvalidWitness): S.V.verify_witness(witness)
        record = H.analyze(witness)
        self.assertEqual(record['permutationCount'],5040)
        self.assertTrue(record['noncommuting'])
        self.assertTrue(R.audit(witness,record))
        with self.assertRaises(ValueError): H.analyze(S.make_witness(next(S.case_specs())))


if __name__ == '__main__': unittest.main()
