"""Helpers of the S414 inviscid cross-formulation check (small linear problems only)."""
from pathlib import Path
import sys
from types import SimpleNamespace
import unittest
import numpy as np
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT/'validation'))
import report_s414_formulation_check as check
from multifoil import signed_area
from viscous import CoupledSolver, mf


class FormulationCheckHelpers(unittest.TestCase):
    def test_declared_schedule_and_identical_pressure_correction(self):
        self.assertEqual(check.SCHEDULE, ((340, 240), (480, 340), (680, 480)))
        model = SimpleNamespace(geom=mf.Geom(), oper=mf.Oper(), param=mf.Param())
        model.geom.chord = 1.; model.oper.Re = 1e6; model.oper.Ma = check.MACH; mf.init_thermo(model)
        speed = np.array([0., .4, 1., 1.6])
        np.testing.assert_allclose(check.cp_kt(speed), mf.get_cp(speed, model.param)[0], rtol=0, atol=1e-13)

    def test_operator_on_prescribed_nodes_is_the_production_operator(self):
        solver = CoupledSolver([('a', check.naca0012_open(80))])
        x = solver.models[0].foil.x
        direct = check.production_operator([('a', x)])
        np.testing.assert_array_equal(direct.gref, solver.gref)
        np.testing.assert_array_equal(direct.aic, solver.aic)

    def test_exact_control_is_closed_counterclockwise_and_satisfies_kutta(self):
        xy, theta, speed, lift = check.karman_trefftz(120)
        np.testing.assert_array_equal(xy[0], xy[-1]); self.assertGreater(signed_area(xy), 0)
        self.assertAlmostEqual(xy[:, 0].max()-xy[:, 0].min(), 1., places=12)
        self.assertGreater(xy[1, 1], xy[-2, 1])  # upper surface is traversed first
        # A finite-angle trailing edge is a (weak, power-law) stagnation point: speed falls toward it.
        near = [speed(theta[:1]+step, 4.)[0] for step in (1e-3, 1e-6, 1e-9)]
        self.assertTrue(np.isfinite(near).all() and near[0] > near[1] > near[2] > 0)
        self.assertGreater(lift(4.), lift(0.)); self.assertGreater(lift(0.), 0.)

    def test_both_methods_use_one_clockwise_midpoint_convention(self):
        pair = check.Pair([('a', check.clockwise(check.naca0012_open(120)))], far=False)
        e = pair.solve(4.)['elements']['a']; n = e['x'].shape[1]
        self.assertEqual((len(e['production']), len(e['hess'])), (n-1, n-1))
        inner = (e['xm'][0] > .05) & (e['xm'][0] < .9); lower = np.arange(n-1) < (n-1)//2
        for label in ('production', 'hess'):
            self.assertTrue(np.all(e[label][inner & lower] < 0) and np.all(e[label][inner & ~lower] > 0))
        self.assertLess(np.max(abs(e['hess'][inner]-e['production'][inner])), .02)
        s = check.crossing(e['s_mid'], e['production'])
        self.assertLess(np.interp(s, e['s_nodes'], e['x'][0]), .01)

    def test_taps_outside_the_sampled_branch_are_not_extrapolated(self):
        pair = check.Pair([('fore', check.clockwise(check.naca0012_open(120)))], far=False)
        solution = pair.solve(2.)
        values = check.tap_values(solution, 'hess', [('fore_upper', .5), ('fore_lower', .5), ('fore_upper', 1.5)])
        self.assertTrue(np.isfinite(values[:2]).all()); self.assertTrue(np.isnan(values[2]))


if __name__ == '__main__':
    unittest.main()
