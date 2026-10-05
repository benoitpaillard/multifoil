"""V5 tranche 2: Orr-Sommerfeld solver, inverse boundary layer and e^N bookkeeping (no coupled solve)."""
from pathlib import Path
import json
import sys
import unittest
import numpy as np
from scipy.interpolate import PchipInterpolator
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT/'validation'))
import s414_laminar_fd as fd
import s414_bubble_fd as bf
import s414_stability as stability
import s414_amplification as amp

STUDY = ROOT/'validation/results/s414_bubble_stability_20261005'
NU = 1e-6


class OrrSommerfeld(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.eta = bf.eta_grid()
        cls.profile = fd.similarity(0., cls.eta)
        cls.d1, cls.d2 = bf.measures(cls.eta, *cls.profile)[:2]
        cls.grid = stability.Grid(140, 60., 3.)
        cls.u, cls.upp = stability.on_grid(cls.grid, cls.eta/cls.d1, cls.profile[1],
                                           np.gradient(cls.profile[2], cls.eta)*cls.d1**2)

    def test_blasius_eigenvalues_match_published_values(self):
        frequency, shape = stability.temporal(.308, 998., self.u, self.upp, self.grid)
        self.assertAlmostEqual((frequency/.308).real, .3641, places=3); self.assertAlmostEqual((frequency/.308).imag, .0080, places=3)
        alpha, _, ok = stability.spatial(.1122, 998., self.u, self.upp, self.grid, .30-.004j, shape)
        self.assertTrue(ok); self.assertAlmostEqual(alpha.real, .3086, places=3); self.assertAlmostEqual(alpha.imag, -.0057, places=3)

    def test_seed_from_a_scan_leads_to_the_same_wave(self):
        found = stability.scan(998., self.u, self.upp, self.grid, np.linspace(.2, .4, 11))
        start, shape = stability.seed(.1122, found)
        alpha, _, ok = stability.spatial(.1122, 998., self.u, self.upp, self.grid, start, shape)
        self.assertTrue(ok); self.assertAlmostEqual(alpha.real, .30859, places=4)
        self.assertEqual(stability.seed(5., found), (None, None))

    def test_operator_derivative_is_consistent(self):
        alpha, omega, step = .31-.005j, .1122, 1e-3  # the operator is a quartic in alpha; a small step only adds round-off
        op, op_alpha = stability.operator(alpha, omega, 998., self.u, self.upp, self.grid)
        plus = stability.operator(alpha+step, omega, 998., self.u, self.upp, self.grid)[0]
        minus = stability.operator(alpha-step, omega, 998., self.u, self.upp, self.grid)[0]
        difference = (plus-minus)/(2*step)
        self.assertLess(np.max(abs(difference-op_alpha))/np.max(abs(op_alpha)), 1e-6)


class InverseBoundaryLayer(unittest.TestCase):
    def test_reversed_similar_profiles_match_tabulated_wall_shear(self):
        eta = bf.eta_grid()
        members = bf.family(np.round(np.arange(3.5, 7.01, .25), 2), eta)
        table = np.array([(2*m/(m+1), bf.measures(eta, *p)[2]/np.sqrt((m+1)/2)) for _, m, p in members])
        self.assertAlmostEqual(table[0, 0], -.1988, places=3)  # separation profile
        for beta, reference in ((-.18, -.0976), (-.15, -.1334)):
            self.assertAlmostEqual(np.interp(beta, table[:, 0], table[:, 1]), reference, delta=1e-3)

    def test_direct_mode_equals_tranche_one_and_round_trip_returns_unit_ratio(self):
        nodes = np.linspace(.02, 1.5, 38); speed = lambda x: 2*np.sin(x)
        new, old = bf.march(nodes, speed, NU), fd.march(nodes, speed, NU, m_start=1.)
        np.testing.assert_allclose(new['theta'], old['theta'], rtol=1e-8)
        back = bf.march(nodes, speed, NU, switch=.6, displacement=PchipInterpolator(new['xi'], new['delta_star']))
        self.assertEqual(back['status'], 'completed'); self.assertLess(np.max(abs(back['ratio']-1)), 1e-8)
        self.assertTrue(back['inverse'][-1] and not back['inverse'][0])

    def test_inverse_mode_passes_separation_and_levels_the_edge_velocity(self):
        nodes = np.linspace(.02, 2.1, 53); speed = lambda x: 2*np.sin(min(x, 1.9))
        direct = bf.march(np.linspace(.02, 1.6, 41), speed, NU)
        own = PchipInterpolator(direct['xi'], direct['delta_star'])
        d0, slope = float(own(1.6)), float(own.derivative()(1.6))/float(own(1.6))
        target = lambda x: float(own(x)) if x <= 1.6 else d0*np.exp(slope*(x-1.6)*(1+2*(x-1.6)))
        result = bf.march(nodes, speed, NU, switch=1.3, displacement=target)
        self.assertEqual(result['status'], 'completed')
        self.assertLess(result['cf'][-1], 0); self.assertLess(result['reverse'][-1], 0); self.assertGreater(result['H'][-1], 5)
        self.assertGreater(result['ue'][-1], 2*np.sin(2.1))  # plateau: above the attached potential value

    def test_join_hook_starts_the_prescription_from_the_solver_s_own_value(self):
        import run_s414_bubble_stability as study
        target = study.Displacement(np.array([1., 1.5, 2.]), np.array([2e-4, 3e-4, 5e-4]))
        target.join(1., 1.9e-4)
        self.assertAlmostEqual(target(1.), 1.9e-4); self.assertAlmostEqual(target(2.), 5e-4)
        scaled = study.Displacement(np.array([1., 1.5, 2.]), np.array([2e-4, 3e-4, 5e-4]), 'linear', True)
        scaled.join(1., 1.9e-4)
        self.assertAlmostEqual(scaled(2.), 5e-4*.95); self.assertAlmostEqual(scaled(1.25), .5*(1.9e-4+3e-4*.95))


class Envelope(unittest.TestCase):
    def test_integration_starts_at_the_neutral_point_and_drops_lost_waves(self):
        xi = np.array([0., 1., 2., 3., 4.])
        rate = np.array([[-1., 1., 1., 1., 1.],          # neutral point halfway through the first interval
                         [np.nan, np.nan, 2., 2., np.nan],  # appears at station 2, lost at station 4
                         [-1., -1., -1., -1., -1.]])        # never amplified
        each, total = amp.envelope(xi, rate)
        np.testing.assert_allclose(each[0], [0., .25, 1.25, 2.25, 3.25])
        np.testing.assert_allclose(each[1, :4], [0., 0., 0., 2.]); self.assertTrue(np.isnan(each[1, 4]))
        np.testing.assert_allclose(each[2], 0.)
        np.testing.assert_allclose(total, [0., .25, 1.25, 2.25, 3.25])


class SavedStudy(unittest.TestCase):
    def test_saved_summary_respects_the_declared_limits(self):
        summary = json.loads((STUDY/'summary.json').read_text())
        self.assertEqual(summary['state'], 'finished'); self.assertTrue(summary['part_a_gate_passed'])
        self.assertLessEqual(summary['counts']['marches'], 15); self.assertLessEqual(summary['counts']['envelopes'], 8)
        self.assertEqual(summary['counts']['coupled_flow_solves'], 0)
        self.assertEqual(set(summary['closure_screen']), {'cl0056', 'cl0175'})


if __name__ == '__main__':
    unittest.main()
