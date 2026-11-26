import copy
import importlib.util
import json
from pathlib import Path
import unittest

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('paired_cone_check', HERE / 'check.py')
check = importlib.util.module_from_spec(spec)
spec.loader.exec_module(check)


class CertificateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = json.loads((HERE / 'certificate.json').read_text())

    def test_exact_obstruction_and_all_repairs(self):
        report = check.verify(self.data)
        self.assertEqual(report['separation'], -20832)
        self.assertEqual(report['violations_per_pr81_orbit'], [0, 0, 0, 2, 1, 2, 2])
        self.assertEqual(len(report['single_ray_repairs']), 7)

    def test_dual_corruption_is_rejected(self):
        data = copy.deepcopy(self.data)
        data['dual_numerators'][0] += 1
        with self.assertRaises(ValueError):
            check.verify(data)

    def test_spectral_corruption_is_rejected(self):
        data = copy.deepcopy(self.data)
        data['a'][0][0] += 1
        with self.assertRaises(ValueError):
            check.verify(data)

    def test_repair_coefficient_corruption_is_rejected(self):
        data = copy.deepcopy(self.data)
        data['repairs'][0]['positive_terms'][0][1] = '1'
        with self.assertRaises(ValueError):
            check.verify(data)

    def test_missing_repair_is_rejected(self):
        data = copy.deepcopy(self.data)
        data['repairs'].pop()
        with self.assertRaises(ValueError):
            check.verify(data)


if __name__ == '__main__':
    unittest.main()
