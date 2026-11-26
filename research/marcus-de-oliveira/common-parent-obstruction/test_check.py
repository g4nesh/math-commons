import copy
import importlib.util
from pathlib import Path
import unittest

spec=importlib.util.spec_from_file_location('obstruction_check',Path(__file__).with_name('check.py'))
c=importlib.util.module_from_spec(spec)
spec.loader.exec_module(c)
F=c.F


class CommonParentObstruction(unittest.TestCase):
    def test_full_physical_replay(self):
        result=c.verify()
        self.assertFalse(result['common_parent_for_these_two_directions'])
        self.assertTrue(result['separate_positive_laws_for_both'])
        self.assertEqual(result['squared_minors_checked'],1848)

    def test_dual_tampering(self):
        data=copy.deepcopy(c.load())
        data['dual'][30]='0'
        with self.assertRaises(ValueError):
            c.verify(data,minors=False)

    def test_invalid_norm_lower_bound(self):
        data=copy.deepcopy(c.load())
        data['norm_lower_bounds'][9]='100'
        with self.assertRaisesRegex(ValueError,'norm lower bounds'):
            c.verify(data,minors=False)

    def test_wrong_direction_rejected(self):
        data=copy.deepcopy(c.load())
        data['directions'][0]['direction'][0]='1'
        with self.assertRaisesRegex(ValueError,'unit Bloch'):
            c.verify(data,minors=False)

    def test_weight_redistribution_cannot_hide_coefficient_failure(self):
        data=copy.deepcopy(c.load())
        weights=data['directions'][0]['weights']
        # Preserve normalization and positivity; the polynomial identity must detect the change.
        high=sorted(range(33),key=lambda j:F(weights[j]),reverse=True)[:2]
        epsilon=min(F(weights[j]) for j in high)/10
        weights[high[0]]=str(F(weights[high[0]])+epsilon)
        weights[high[1]]=str(F(weights[high[1]])-epsilon)
        with self.assertRaises(ValueError):
            c.verify(data,minors=False)


if __name__=='__main__':
    unittest.main()
