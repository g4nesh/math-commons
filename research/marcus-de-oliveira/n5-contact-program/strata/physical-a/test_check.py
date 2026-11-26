import importlib.util
from pathlib import Path
import unittest

spec=importlib.util.spec_from_file_location('physical_a_check',Path(__file__).with_name('check.py'))
c=importlib.util.module_from_spec(spec)
spec.loader.exec_module(c)


class PhysicalFacetA(unittest.TestCase):
    def test_facet_and_orbit(self):
        result=c.facet_data()
        self.assertEqual(result['physical_symmetry_orbit'],900)
        self.assertEqual(result['face_dimension'],40)

    def test_separate_equality_branches_and_strict_point(self):
        result=c.verify()['fixtures']
        self.assertEqual(result['block_zero_contact']['rhs'],'9/16')
        self.assertEqual(result['column_supported_contact']['rhs'],'10/121')
        self.assertEqual(result['zero_rhs_contact']['tau'],'75/121')
        self.assertEqual(result['strict_dense']['gap'],'17008/75625')

    def test_complex_row_and_column_phase_invariance(self):
        u=c.fixtures()['strict_dense']
        phase=c.Q(c.F(3,5),c.F(4,5))
        changed=[[x*phase if i==3 or j==2 else x for j,x in enumerate(row)]for i,row in enumerate(u)]
        # Apply the two phase changes separately at their intersection.
        changed[3][2]*=phase
        before=c.gap_data(u);after=c.gap_data(changed)
        self.assertEqual(before,after)

    def test_nonunitary_rejected(self):
        u=c.fixtures()['strict_dense'];u[0][0]+=c.Q(c.F(1,1000))
        with self.assertRaisesRegex(ValueError,'unitarity'):c.gap_data(u)


if __name__=='__main__':unittest.main()
