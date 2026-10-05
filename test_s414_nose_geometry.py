"""Fore-nose reconstruction alternatives of stage V4 (geometry only; no flow solve)."""
from pathlib import Path
import sys
import unittest
import numpy as np
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT/'validation'))
import s414_nose_geometry as nose
from multifoil import repanel
from s414_geometry import elements as declared
from viscous import new_model


def solver_polygon(points):
    return new_model('fore', points, 1e6, .1, 9., (.25, 0.), 1e-5).foil.x.T[::-1]


def nose_turns(polygon):
    turn = nose.turns(polygon)
    return turn[polygon[1:-1, 0] <= nose.NOSE]


class NoseAlternatives(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.xy, cls.le = nose.source(); cls.s = nose.chord_parameter(cls.xy)
        cls.inputs = {v: nose.fore_input(v) for v in nose.VARIANTS}
        cls.parameters = nose.cosine_parameters(cls.s, cls.le, 340)

    def test_nominal_is_the_production_input_bit_for_bit(self):
        np.testing.assert_array_equal(self.inputs['nominal'], repanel(self.xy, 340))
        np.testing.assert_array_equal(self.inputs['nominal'], declared('nominal', [340, 240])[0][1])
        for variant in nose.VARIANTS:
            np.testing.assert_array_equal(nose.elements(variant)[1][1], declared('nominal', [340, 240])[1][1])

    def test_reported_polygon_is_convex_in_the_support_with_one_nearly_flat_vertex(self):
        i0, i1 = nose.support(self.xy, self.le); turn = nose.turns(self.xy)[i0-1:i1]
        self.assertEqual((i0, i1, i1-i0+1), (33, 47, 15)); self.assertTrue(np.all(turn > 0))
        self.assertEqual(int(np.argmin(turn))+i0, 45); self.assertAlmostEqual(turn.min(), 1.0077, places=3)

    def test_only_the_declared_samples_change_and_every_reported_vertex_is_kept(self):
        i0, i1 = nose.support(self.xy, self.le)
        for variant, lo, hi in (('single_curve', self.le-1, self.le+1), ('convex_nose', i0, i1)):
            changed = np.any(self.inputs[variant] != self.inputs['nominal'], axis=1)
            inside = (self.parameters > self.s[lo]) & (self.parameters < self.s[hi])
            self.assertTrue(changed.any()); self.assertFalse(np.any(changed & ~inside))
            curve = nose.Reconstruction(variant)
            for i in range(len(self.xy)):
                np.testing.assert_allclose(curve.point(float(self.s[i])), self.xy[i], rtol=0, atol=1e-14)
            np.testing.assert_array_equal(self.inputs[variant][[0, -1]], self.inputs['nominal'][[0, -1]])

    def test_arc_chain_is_convex_and_tangent_continuous(self):
        curve = nose.Reconstruction('convex_nose'); i0, i1 = curve.support
        self.assertTrue(all(k > 0 for pair in curve.arc_curvatures.values() for k in pair))
        self.assertEqual(curve.fairing['at_bound'], 0)
        self.assertLess(curve.fairing['energy_minimum'], curve.fairing['energy_three_point_start'])
        step = 1e-7
        for i in range(i0, i1+1):
            expected = np.array([np.cos(curve.tangent_angles[i]), np.sin(curve.tangent_angles[i])])
            for side in (-1, 1):
                if (i == i0 and side < 0) or (i == i1 and side > 0):
                    # Outside the support the production curve continues with the same tangent.
                    d = side*(curve.whole(self.s[i]+side*step)-self.xy[i])
                else:
                    d = side*(curve.point(float(self.s[i]+side*step*(self.s[i+side]-self.s[i])*side))-self.xy[i])
                np.testing.assert_allclose(d/np.linalg.norm(d), expected, rtol=0, atol=2e-5)

    def test_reversal_is_absent_only_in_the_convex_alternative(self):
        smallest = {v: float(nose_turns(solver_polygon(p)).min()) for v, p in self.inputs.items()}
        self.assertLess(smallest['nominal'], -.6); self.assertLess(smallest['single_curve'], -.6)
        self.assertGreater(smallest['convex_nose'], 0)
        self.assertTrue(np.all(nose_turns(self.inputs['convex_nose']) > 0))
        self.assertLess(np.linalg.norm(self.inputs['convex_nose']-self.inputs['nominal'], axis=1).max()*457.2, .13)

    def test_construction_refuses_data_that_are_not_convex(self):
        p0, p1 = np.array([0., 0.]), np.array([1., 0.]); up = np.array([np.cos(.2), np.sin(.2)]); down = np.array([np.cos(.2), -np.sin(.2)])
        sample, curvatures, _ = nose.biarc(p0, down, p1, up)
        self.assertTrue(curvatures[0] > 0 and curvatures[1] > 0)
        np.testing.assert_allclose(sample(0.), p0, atol=1e-14); np.testing.assert_allclose(sample(1.), p1, atol=1e-14)
        with self.assertRaisesRegex(ValueError, 'not convex-consistent'): nose.biarc(p0, up, p1, up)
        bent = self.xy.copy(); bent[45, 1] -= .002  # push the nearly flat vertex outward: the polygon becomes concave there
        with self.assertRaisesRegex(ValueError, 'not convex in the nose support'):
            nose.fair_convex_tangents(bent, 33, 47, np.array([1., 0.]), np.array([1., 0.]))


if __name__ == '__main__':
    unittest.main()
