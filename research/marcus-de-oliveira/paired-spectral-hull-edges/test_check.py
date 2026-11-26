import importlib.util
import json
from pathlib import Path
import unittest

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('hull_certificate', HERE / 'check.py')
check = importlib.util.module_from_spec(spec)
spec.loader.exec_module(check)


class EdgeCertificates(unittest.TestCase):
    def setUp(self):
        self.data = json.loads((HERE / 'certificate.json').read_text())

    def test_exact_certificate_and_report(self):
        self.assertEqual(check.verify(self.data), json.loads((HERE / 'report.json').read_text()))

    def test_changed_positive_weight_rejected(self):
        row = self.data['edges'][0]['weights']
        j = next(j for j, x in enumerate(row) if check.F(x) > 0)
        row[j] = str(check.F(row[j]) + 1)
        with self.assertRaisesRegex(ValueError, 'exact edge certificate identity'):
            check.verify(self.data)

    def test_missing_edge_rejected(self):
        self.data['edges'].pop()
        with self.assertRaisesRegex(ValueError, 'every polygon edge'):
            check.verify(self.data)

    def test_reverse_orientation_rejected(self):
        self.data['edges'][0]['ends'].reverse()
        with self.assertRaisesRegex(ValueError, 'edge order and orientation'):
            check.verify(self.data)

    def test_positive_coefficient_outside_face_rejected(self):
        row = self.data['edges'][0]['weights']
        j = next(j for j, x in enumerate(row) if check.F(x) == 0)
        row[j] = '1'
        with self.assertRaisesRegex(ValueError, 'strictly positive face representation'):
            check.verify(self.data)

    def test_hull_degenerate_and_duplicate_cases(self):
        self.assertEqual(check.hull([(0, 0), (1, 0), (2, 0), (0, 0)]), [0, 2])
        self.assertEqual(check.hull([(1, 2), (1, 2)]), [0])
        self.assertEqual(check.hull([(0, 0), (1, 0), (1, 1), (0, 1), (1, 0)]), [0, 1, 2, 3])


if __name__ == '__main__':
    unittest.main()
