"""Contract tests for the fixed-lift objective; these are not aerodynamic validation."""
import unittest
from types import SimpleNamespace
import numpy as np
from optimization import solve_fixed_lift, LiftTargetError, camber_bump, minimum_clearance


class AnalyticSolver:
    def __init__(self, fail=False, supersonic=False):
        self.calls = []
        self.fail, self.supersonic = fail, supersonic

    def solve(self, alpha):
        self.calls.append(alpha)
        return SimpleNamespace(alpha=alpha, cl=.2*alpha+.003*alpha**2,
            cd=.01+.0002*alpha**2, cm=-.1, converged=not self.fail,
            diagnostics={'local_supersonic_detected': {'main': self.supersonic}})


class FixedLiftTests(unittest.TestCase):
    def test_nonlinear_target_and_callback(self):
        solver = AnalyticSolver()
        trace = []
        result, samples = solve_fixed_lift(solver, 1., (0., 8.), 1e-7, on_evaluation=trace.append)
        self.assertLess(abs(result.cl-1.), 1e-7)
        self.assertEqual(trace, samples)
        self.assertLessEqual(len(samples), 16)
        self.assertTrue(all(0 <= s.alpha <= 8 for s in samples))

    def test_failed_flow_cannot_be_an_objective(self):
        with self.assertRaises(LiftTargetError):
            solve_fixed_lift(AnalyticSolver(fail=True), 1., (0., 8.))

    def test_supersonic_cannot_be_an_objective(self):
        with self.assertRaises(LiftTargetError):
            solve_fixed_lift(AnalyticSolver(supersonic=True), 1., (0., 8.))

    def test_no_extrapolation(self):
        solver = AnalyticSolver()
        with self.assertRaises(LiftTargetError):
            solve_fixed_lift(solver, 3., (0., 8.))
        self.assertEqual(solver.calls, [0., 8.])

    def test_budget_and_input(self):
        solver = AnalyticSolver()
        with self.assertRaises(LiftTargetError):
            solve_fixed_lift(solver, 1., (0., 8.), 1e-14, max_evaluations=3)
        self.assertEqual(len(solver.calls), 3)
        with self.assertRaises(ValueError):
            solve_fixed_lift(solver, 1., (8., 0.))

    def test_camber_keeps_thickness_and_endpoints(self):
        xy = np.array([[1., .001], [.5, .05], [0., 0.], [.5, -.05], [1., -.001]])
        changed = camber_bump(xy, .01)
        np.testing.assert_allclose(changed[[0, 2, 4]], xy[[0, 2, 4]])
        self.assertAlmostEqual(changed[1, 1]-changed[3, 1], .1)
        self.assertAlmostEqual(changed[1, 1]-xy[1, 1], .01)

    def test_clearance_uses_segments(self):
        a = np.array([[0., 0.], [1., 0.], [1., 1.], [0., 1.]])
        b = np.array([[.4, 1.2], [.6, 1.2], [.6, 1.4], [.4, 1.4]])
        self.assertAlmostEqual(minimum_clearance(a, b), .2)


if __name__ == '__main__':
    unittest.main()
