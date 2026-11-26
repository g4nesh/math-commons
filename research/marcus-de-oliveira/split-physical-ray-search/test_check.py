import importlib.util
import json
from pathlib import Path
import unittest

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('ray_search_check',HERE/'check.py')
c=importlib.util.module_from_spec(spec)
spec.loader.exec_module(c)


class ExactEndpoints(unittest.TestCase):
    def test_identity_chart(self):
        data=json.loads((HERE/'search.json').read_text())
        vectors=data['start_vectors'][0]
        self.assertEqual(c.rational_unitary(vectors,[0]*36),c.base.fixture({'vectors':vectors}))

    def test_complex_chart_is_exactly_unitary(self):
        data=json.loads((HERE/'search.json').read_text())
        theta=[(-1)**i*(i+1)*1000 for i in range(36)]
        u=c.rational_unitary(data['start_vectors'][0],theta)
        self.assertEqual(c.exact.multiply(c.exact.adjoint(u),u),c.exact.identity(6))

    def test_float_or_boolean_chart_rejected(self):
        data=json.loads((HERE/'search.json').read_text())
        for value in (0.0,True):
            theta=[0]*36;theta[4]=value
            with self.assertRaises(ValueError):c.rational_unitary(data['start_vectors'][0],theta)

    def test_all_exact_endpoints(self):
        result=c.verify()
        self.assertEqual(result['completed_runs'],140)
        self.assertEqual(result['exact_positive'],35)
        self.assertFalse(result['universal_sign_claim'])


if __name__=='__main__':unittest.main()
