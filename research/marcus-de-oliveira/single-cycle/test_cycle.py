import unittest
from fractions import Fraction as F
from itertools import permutations
import cycle as c
import cover_certificate as cover


class CycleTests(unittest.TestCase):
    def test_every_triangle_order_complex_and_glued_fixtures(self):
        fixtures = list(c.fixtures())
        self.assertEqual(len(fixtures), 14)
        total = 0
        for name, u, w in fixtures:
            with self.subTest(name=name):
                count, mismatch = c.compound_check(u, w)
                self.assertIsNone(mismatch)
                self.assertTrue(c.determinant_check(u, w)['equal'])
                total += count
        self.assertEqual(total, 1416)

    def test_formula_endpoints_and_midpoint(self):
        u, _ = c.block3(True)
        d = [[c.f.norm2(z) for z in row] for row in u]
        _, (lo, hi) = c.triangle_weights(d)
        for x in (lo, hi, (lo+hi)/2):
            w, _ = c.triangle_weights(d, x)
            self.assertIsNone(c.compound_check(u, w)[1])

    def test_all_permutation_boundary_matrices(self):
        for p in permutations(range(3)):
            d = [[F(int(p[j] == i)) for j in range(3)] for i in range(3)]
            w, interval = c.triangle_weights(d)
            self.assertEqual(w, {p:F(1)})
            self.assertEqual(interval[0], interval[1])

    def test_exact_triangle_repair(self):
        u, gates = c.block3()
        w, interval = c.triangle_weights([[c.f.norm2(z) for z in row] for row in u])
        self.assertEqual(interval, (F(1881,15625), F(81,625)))
        old = c.f.permutation_weights(3, gates, require_forest=False)
        self.assertEqual(sum(v for p,v in old.items() if p[0] == 0), F(6121,15625))
        self.assertEqual(sum(v for p,v in w.items() if p[0] == 0), F(11881,15625))

    def test_four_cycle_singletons_do_not_determine_middle_layer(self):
        block = [[c.Q(F(3,5)), c.Q(F(4,5))], [c.Q(F(-4,5)), c.Q(F(3,5))]]
        u = c.f.build_unitary(4, [(0,1,block),(1,2,block),(2,3,block),(0,3,block)], require_forest=False)
        w = c.greedy_birkhoff([[c.f.norm2(z) for z in row] for row in u])
        for i in range(4):
            for j in range(4):
                self.assertEqual(sum(v for p,v in w.items() if p[j] == i), c.f.norm2(u[i][j]))
        self.assertEqual(c.compound_check(u,w)[1], {'rows':[0,1], 'columns':[0,1], 'actual':'81/625', 'predicted':'3321/15625'})
        result = c.determinant_check(u,w)
        self.assertFalse(result['equal'])
        self.assertEqual(result['originalHullMembership'], 'inside')

    def test_gluing_rejects_shared_pair_and_false_support(self):
        block, _ = c.block3(True)
        left = c.embed(4,[0,1,2],block,c.certificate(block))
        right = c.embed(4,[1,2,3],block,c.certificate(block))
        with self.assertRaises(ValueError):
            c.glue(left,right)
        fake = (left[0], left[1], frozenset({0}))
        with self.assertRaises(ValueError):
            c.glue(fake,c.embed(4,[2,3],c.eye(2),c.certificate(c.eye(2))))

    def test_rejects_inexact_invalid_and_out_of_interval(self):
        for d in ([[1.,0,0],[0,1,0],[0,0,1]], [[True,0,0],[0,1,0],[0,0,1]], [[1,0,0],[0,1,0],[0,1,0]]):
            with self.assertRaises(ValueError):
                c.triangle_weights(d)
        with self.assertRaises(ValueError):
            c.triangle_weights([[F(1,3)]*3 for _ in range(3)],F(-1))
        with self.assertRaises(ValueError):
            c.certificate([[c.Q(2),c.Z],[c.Z,c.O]])
        with self.assertRaises(ValueError):
            c.certificate(c.eye(4))

    def test_singular_and_repeated_spectra(self):
        u, _ = c.block3(True)
        w = c.certificate(u)
        result = c.determinant_check(u,w,[c.Q(2,1)]*3,[c.Q(-2,-1)]*3)
        self.assertEqual(result['actual'], ['0','0'])
        self.assertTrue(result['equal'])
        result = c.determinant_check(u,w,[c.Z,c.Z,c.O],[c.O,c.Z,c.Z])
        self.assertTrue(result['equal'])

    def test_standalone_cover_proves_no_common_certificate(self):
        result = cover.certify()
        self.assertEqual(result['unitaryCoverValue'], '14953/15625')
        self.assertEqual(result['strictGap'], '672/15625')
        self.assertEqual(len(result['permutationChecks']), 24)
        self.assertGreaterEqual(min(r['hits'] for r in result['permutationChecks']), 1)
        self.assertEqual(sum(F(r['squaredMinor']) for r in result['cover']), F(14953,15625))
        # Independent Fraction-only matrix construction agrees with Q(i) code.
        block = [[c.Q(F(3,5)), c.Q(F(4,5))], [c.Q(F(-4,5)), c.Q(F(3,5))]]
        u = c.f.build_unitary(4, [(0,1,block),(1,2,block),(2,3,block),(0,3,block)], require_forest=False)
        self.assertEqual(u, [[c.Q(z) for z in row] for row in cover.matrix()])

    def test_cover_checker_rejects_missing_witness_and_duplicates(self):
        with self.assertRaises(ValueError):
            cover.check_cover(cover.COVER[1:])
        with self.assertRaises(ValueError):
            cover.check_cover(cover.COVER+cover.COVER[:1])


if __name__ == '__main__':
    unittest.main()
