from pathlib import Path
import copy
import importlib.util
import json
import unittest

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('coarsening',HERE/'check.py')
c = importlib.util.module_from_spec(spec)
spec.loader.exec_module(c)


class CoarseningTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.spaces,cls.records = c.construct()
        cls.certificate = json.loads((HERE/'certificate.json').read_text())

    def test_complete_certificate(self):
        self.assertTrue(c.verify(self.certificate,self.spaces,self.records))

    def test_orbit_three_requires_actual_row_space(self):
        r = self.records[3]
        self.assertEqual([k for k,v in r['outcomes'].items() if v['kind']=='member'],['row12'])
        self.assertFalse(r['coefficient'])
        c.validate_outcome(r['ray'],self.spaces['row12'],r['outcomes']['row12'])

    def test_altered_member_combination_rejected(self):
        r=self.records[3];bad=copy.deepcopy(r['outcomes']['row12'])
        bad['combination'][0][1]=str(c.F(bad['combination'][0][1])+1)
        with self.assertRaisesRegex(ValueError,'row-space certificate'):
            c.validate_outcome(r['ray'],self.spaces['row12'],bad)

    def test_zero_separator_rejected(self):
        r=self.records[10];bad=copy.deepcopy(r['outcomes']['row12']);bad['dot']=0
        with self.assertRaisesRegex(ValueError,'separating pairing'):
            c.validate_outcome(r['ray'],self.spaces['row12'],bad)

    def test_missing_coarsening_rejected(self):
        bad=copy.deepcopy(self.certificate);del bad['records'][0]['outcomes']['row23']
        with self.assertRaisesRegex(ValueError,'nine coarsenings'):
            c.verify(bad,self.spaces,self.records)

    def test_coefficient_exception_not_called_new_gap(self):
        r=self.records[38]
        self.assertTrue(r['coefficient']);self.assertFalse(r['covered'])
        self.assertEqual(r['orbit_size'],6)


if __name__=='__main__':
    unittest.main()
