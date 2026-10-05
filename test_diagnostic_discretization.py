"""Numerical compatibility and independence of the opt-in controls."""
import unittest
from unittest.mock import patch
import numpy as np
from multifoil import naca4
from viscous import CoupledSolver
from diagnostic_discretization import DiscretizedSolver


class DiscretizationTests(unittest.TestCase):
    def setUp(self):
        self.elements = [('main', naca4('0012', 81))]

    def test_no_options_preserve_actual_surfaces_and_inviscid_system(self):
        old = CoupledSolver(self.elements)
        new = DiscretizedSolver(self.elements)
        for key in ('aic', 'rhs', 'gref'):
            np.testing.assert_array_equal(getattr(old, key), getattr(new, key))
        for key in ('x', 's', 't'):
            np.testing.assert_array_equal(getattr(old.models[0].foil, key), getattr(new.models[0].foil, key))

    def test_explicit_auto_wake_reproduces_original_tracing_and_coupling(self):
        old = CoupledSolver(self.elements)
        old._build_wakes(2.)
        m = old.models[0]
        step = .5*((m.foil.s[1]-m.foil.s[0])+(m.foil.s[-1]-m.foil.s[-2]))
        new = DiscretizedSolver(self.elements, surface_nodes={'main': m.foil.N},
            wake_nodes={'main': m.wake.N}, wake_first_spacing={'main': step})
        for angle in (2., 3.):
            old._build_wakes(angle); new._build_wakes(angle)
            for key in ('x', 's', 't'):
                np.testing.assert_array_equal(getattr(old.models[0].wake, key), getattr(new.models[0].wake, key))
            np.testing.assert_array_equal(old.models[0].isol.uewiref, new.models[0].isol.uewiref)
            self.assertEqual(old.ss, new.ss); self.assertEqual(old.sources, new.sources)

    def test_surface_count_preserves_source_and_wake_controls(self):
        settings = dict(wake_nodes={'main': 30}, wake_first_spacing={'main': .001})
        low = DiscretizedSolver(self.elements, surface_nodes={'main': 61}, **settings)
        high = DiscretizedSolver(self.elements, surface_nodes={'main': 101}, **settings)
        np.testing.assert_array_equal(low.models[0].geom.xpoint, high.models[0].geom.xpoint)
        np.testing.assert_array_equal(low.models[0].foil.x[:, [0, -1]], high.models[0].foil.x[:, [0, -1]])
        self.assertEqual(low.nfoil, 61); self.assertEqual(high.nfoil, 101)
        # Straight velocity isolates station spacing from integration curvature.
        for solver in (low, high):
            with patch.object(solver, '_velocity', side_effect=lambda p,a: np.tile([[1.], [0.]], (1,p.shape[1]))):
                solver._build_wakes(0.)
            self.assertEqual(solver.models[0].wake.N, 30)
            self.assertAlmostEqual(np.diff(solver.models[0].wake.x[0])[0], .001, places=13)
        # TE-panel tangents can shift the tiny starting offset when surface
        # sampling changes; the controlled station increments must agree.
        np.testing.assert_array_equal(np.diff(low.models[0].wake.x), np.diff(high.models[0].wake.x))

    def test_invalid_controls_are_rejected_before_solve(self):
        for kwargs in ({'surface_nodes': {'main': 40}}, {'surface_nodes': {'main': 1201}},
                       {'surface_nodes': {'other': 81}}, {'wake_nodes': {'main': 15}},
                       {'wake_nodes': {'main': 16.5}}, {'wake_nodes': {'main': True}},
                       {'wake_first_spacing': {'main': float('nan')}},
                       {'wake_first_spacing': {'main': 0}}):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                DiscretizedSolver(self.elements, **kwargs)
        solver = DiscretizedSolver(self.elements, wake_nodes={'main': 30}, wake_first_spacing={'main': 1.})
        with self.assertRaisesRegex(ValueError, 'increasing geometric spacing'):
            solver._build_wakes(0.)


class DiscretizationFlowTests(unittest.TestCase):
    def test_explicit_auto_settings_repeat_converged_force_and_newton_state(self):
        elements = [('main', naca4('0012', 81))]
        old = CoupledSolver(elements)
        old._build_wakes(2.)
        m = old.models[0]
        step = .5*((m.foil.s[1]-m.foil.s[0])+(m.foil.s[-1]-m.foil.s[-2]))
        new = DiscretizedSolver(elements, surface_nodes={'main': m.foil.N},
            wake_nodes={'main': m.wake.N}, wake_first_spacing={'main': step})
        a, b = old.solve(2.), new.solve(2.)
        self.assertTrue(a.converged, a.diagnostics['failure_reason'])
        self.assertTrue(b.converged, b.diagnostics['failure_reason'])
        self.assertEqual((a.cl, a.cd, a.cm, a.residual), (b.cl, b.cd, b.cm, b.residual))
        np.testing.assert_array_equal(old.models[0].glob.U, new.models[0].glob.U)


if __name__ == '__main__':
    unittest.main()
