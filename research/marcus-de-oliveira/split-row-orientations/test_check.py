import copy
import importlib.util
from pathlib import Path
import unittest

spec = importlib.util.spec_from_file_location('orientation_check', Path(__file__).with_name('check.py'))
c = importlib.util.module_from_spec(spec)
spec.loader.exec_module(c)
F, Q = c.F, c.Q


class OrientationCertificate(unittest.TestCase):
    def test_full_replay(self):
        report = c.verify()
        self.assertEqual(report['uniform_child_mass_floor'], '1/20000')
        self.assertEqual(len(report['mixed_exact_SU2_directions']), 3)

    def test_axis_and_south_pole_lifts(self):
        u = c.base.fixture(c.base.load())
        for n in ((1,0,0),(0,1,0),(0,0,1),(0,0,-1)):
            q = c.canonical_frame(c.rotate(c.direction_frame(n),u))
            self.assertEqual(q,c.canonical_at(n))
            weights = c.weights_at(n)
            self.assertGreaterEqual(min(weights),F(1,20000))

    def test_response_tampering_rejected(self):
        data = copy.deepcopy(c.load())
        data['odd_response'][0][1] = str(F(data['odd_response'][0][1])+F(1,1000000))
        with self.assertRaisesRegex(ValueError,'odd response'):
            c.verify(data,fixtures=False)

    def test_parent_tampering_rejected(self):
        data = copy.deepcopy(c.load())
        data['repairs'][0] = str(F(data['repairs'][0])+1)
        with self.assertRaisesRegex(ValueError,'parent masses'):
            c.verify(data,fixtures=False)

    def test_nonexact_and_nonunit_direction_rejected(self):
        for n in ((1.0,0,0),(True,0,0),(1,1,0),(0,0),('1',0,0)):
            with self.assertRaises(ValueError):
                c.weights_at(n)

    def test_nonzero_coarse_perturbation_all_direction_formula(self):
        u = c.base.fixture(c.base.load())
        t = F(1,10**10)
        phase = Q((1-t*t)/(1+t*t),2*t/(1+t*t))
        changed = [row[:] for row in u]
        changed[0] = [phase*x for x in changed[0]]
        rotation = list(c.mixed_rotations())[-1]
        weights = c.nearby_weights(changed,rotation)
        self.assertGreaterEqual(min(weights),F(1,40000))
        self.assertEqual(sum(weights),1)
        # Also test a non-gauge rotation coupling row0 to row4, changing the coarse plane.
        a,b = (1-t*t)/(1+t*t),2*t/(1+t*t)
        changed = [row[:] for row in u]
        changed[0] = [a*x+b*y for x,y in zip(u[0],u[4])]
        changed[4] = [-b*x+a*y for x,y in zip(u[0],u[4])]
        self.assertNotEqual(c.base.canonical(changed),c.base.canonical(u))
        self.assertGreaterEqual(min(c.nearby_weights(changed,rotation)),F(1,40000))


if __name__ == '__main__':
    unittest.main()
