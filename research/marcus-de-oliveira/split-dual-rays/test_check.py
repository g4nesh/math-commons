"""Adversarial local checks of the ray and transcript validation boundaries."""
from pathlib import Path
import importlib.util
import tempfile
import unittest

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("split_checker", HERE / "check.py")
c = importlib.util.module_from_spec(spec)
spec.loader.exec_module(c)


class CertificateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tables, _, _, cls.kernel = c.a.structure()
        cls.rays = sorted(c.a.read_lrs(HERE / "lrs-output.txt"))

    def test_full_certificate(self):
        self.assertEqual(c.check()["physical_orbits"], 62)

    def test_scaled_and_negative_rejected(self):
        for bad in ([2*x for x in self.rays[0]], [-x for x in self.rays[0]]):
            with self.assertRaises(ValueError):
                c.a.verify_ray(bad, self.kernel)

    def test_off_space_rejected(self):
        bad = [0]*33
        bad[0] = 1
        with self.assertRaisesRegex(ValueError, "row-space"):
            c.a.verify_ray(bad, self.kernel)

    def test_non_extreme_sum_rejected(self):
        combined = c.a.primitive([x+y for x,y in zip(self.rays[0], self.rays[1])])
        with self.assertRaisesRegex(ValueError, "support kernel"):
            c.a.verify_ray(combined, self.kernel)

    def test_bool_coordinate_rejected(self):
        bad = list(self.rays[0])
        bad[0] = bool(bad[0])
        with self.assertRaises(ValueError):
            c.a.verify_ray(bad, self.kernel)

    def test_missing_vertex_and_completion_rejected(self):
        text = (HERE / "lrs-output.txt").read_text()
        lines = text.splitlines()
        first = next(s for s in lines if len(s.split()) == 34 and s.split()[0] == "1")
        for bad in (text.replace(first + "\n", "", 1), text.replace("*Totals:", "*Partial:")):
            with tempfile.TemporaryDirectory() as temp:
                path = Path(temp) / "bad.txt"
                path.write_text(bad)
                with self.assertRaises(ValueError):
                    c.validate_output(path)


if __name__ == "__main__":
    unittest.main()
