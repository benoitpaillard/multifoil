"""Geometry/metric/wake invariants for the bounded local-resolution study."""
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent/'validation'))
from s414_nose_resolution import Solver, build_elements, subdivide
from viscous import CoupledSolver


class Subdivision(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.elements = build_elements()
        # No potential-flow matrix or nonlinear solve is required to check geometry.
        with patch.object(CoupledSolver, '_build_inviscid'):
            cls.base = Solver(cls.elements, subdivision_factor=1)

    def test_identity_preserves_all_panel_data(self):
        m = self.base.models[0]
        x, s, t, meta = subdivide(m.foil.x, m.foil.s, m.foil.t, 1)
        for actual, expected in ((x, m.foil.x), (s, m.foil.s), (t, m.foil.t)):
            np.testing.assert_array_equal(actual, expected)
        self.assertEqual(meta['original_nodes'], meta['refined_nodes'])
        self.assertEqual(self.base.discretization['wake_nodes'], {})

    def test_subdivision_is_nested_and_does_not_move_original_nodes(self):
        m = self.base.models[0]
        for factor in (2, 4):
            x, s, t, meta = subdivide(m.foil.x, m.foil.s, m.foil.t, factor)
            mapping = meta['original_to_refined']
            np.testing.assert_array_equal(x[:, mapping], m.foil.x)
            np.testing.assert_array_equal(s[mapping], m.foil.s)
            np.testing.assert_array_equal(t[:, mapping], m.foil.t)
            self.assertEqual(len(s), m.foil.N+(factor-1)*len(meta['selected_parent_panels']))
            for j in meta['selected_parent_panels']:
                a = mapping[j]
                for k in range(factor+1):
                    np.testing.assert_allclose(x[:, a+k],
                        m.foil.x[:, j]+(k/factor)*(m.foil.x[:, j+1]-m.foil.x[:, j]), atol=2e-16, rtol=0)
                    self.assertAlmostEqual(s[a+k],m.foil.s[j]+(k/factor)*(m.foil.s[j+1]-m.foil.s[j]),places=15)
            if factor == 2:
                half_x = x
            else:
                for point in half_x.T:
                    self.assertLess(np.min(np.linalg.norm(x.T-point, axis=1)), 2e-16)

    def test_refined_solver_freezes_unaffected_geometry_and_wake_controls(self):
        for factor in (2, 4):
            with patch.object(CoupledSolver, '_build_inviscid'):
                solver = Solver(self.elements, subdivision_factor=factor)
            for a, b in zip(self.base.models, solver.models):
                meta = solver.resolution_controls[a.geom.name]
                self.assertEqual(meta['wake_nodes'],self.base.resolution_controls[a.geom.name]['wake_nodes'])
                self.assertEqual(meta['wake_first_spacing'],self.base.resolution_controls[a.geom.name]['wake_first_spacing'])
                for attr in ('x', 's', 't'):
                    actual, expected = getattr(b.foil, attr), getattr(a.foil, attr)
                    mapping = meta['original_to_refined']
                    np.testing.assert_array_equal(actual[..., mapping], expected)
                if a.geom.name == 'aft':
                    np.testing.assert_array_equal(a.foil.x, b.foil.x)
                np.testing.assert_array_equal(a.foil.x[:, [0, -1]], b.foil.x[:, [0, -1]])

    def test_scope_rejects_invalid_factor_and_external_sampling(self):
        for factor in (0, 3, True, 2.0):
            with self.assertRaises(ValueError):
                Solver(self.elements, subdivision_factor=factor)
        with self.assertRaises(ValueError):
            Solver(self.elements, wake_nodes={'fore':99})
        with self.assertRaises(ValueError):
            Solver(list(reversed(self.elements)))


if __name__ == '__main__':
    unittest.main()
