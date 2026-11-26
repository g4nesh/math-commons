import copy,json,unittest
from fractions import Fraction as F
import exchanges as x

class ExchangeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.data=json.loads((x.HERE/'basis.json').read_text())

    def test_exact_report_and_all_repair_coordinates(self):
        self.assertEqual(x.verify(self.data),json.loads((x.HERE/'verification.json').read_text()))

    def test_exact_regeneration(self):
        self.assertEqual(x.generate(),self.data)

    def test_each_move_has_three_positive_and_three_negative_atoms(self):
        moves,_=x.moves((2,2,2))
        for move in moves:
            self.assertEqual(move.count(1),3);self.assertEqual(move.count(-1),3);self.assertEqual(sum(move),0)

    def test_tampered_coefficient_is_rejected(self):
        data=copy.deepcopy(self.data);data['interiorRepairCoefficients'][0]=str(F(data['interiorRepairCoefficients'][0])+1)
        with self.assertRaises(ValueError):x.verify(data)

    def test_missing_or_wrong_basis_is_rejected(self):
        data=copy.deepcopy(self.data);data['basisMoveIndices'].pop()
        with self.assertRaises(ValueError):x.verify(data)

    def test_non_kernel_vector_cannot_be_reconstructed(self):
        moves,_=x.moves((2,2,2));_,indices=x.rank_and_rows(moves)
        with self.assertRaises(ValueError):x.solve_coefficients([moves[i]for i in indices],[F(1)]+[F(0)]*89)

    def test_four_color_cone_dual_separates_known_obstruction(self):
        result=x.four_color_dual()
        self.assertEqual(result['canonicalDualValue'],'-672/15625')
        self.assertEqual(result['exchangeDirectionsAnnihilated'],16)

    def test_degenerate_counts_and_exact_types(self):
        report,_=x.check_pattern((3,0,1));self.assertEqual(report['kernelDimension'],0)
        for pattern in ((2,-1,2),(True,1,1),(),(8,0,0)):
            with self.assertRaises(ValueError):x.colorings(pattern)
        for vector in ([0.0],[True]):
            with self.assertRaises(TypeError):x.solve_coefficients([[1]],vector)

if __name__=='__main__':unittest.main()
