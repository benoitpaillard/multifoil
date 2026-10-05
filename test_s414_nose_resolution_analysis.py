"""Pure diagnostic-screen checks; no flow solver or archived result is executed."""
import sys
from pathlib import Path
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent/'validation'))
from report_s414_nose_resolution import (assess, common_profiles, compare_topology,
    contraction, location_on_polygon, supported_interpolation)


def state(factor, pressure_error=0., wake_shift=0.):
    x = np.linspace(0., .04, 2*factor+1); n = len(x); xy = np.array([x, np.zeros(n)])
    u = np.array([1e-4*(1+x), 2e-4*(1+x), 100*x, 1+x])
    surface = np.zeros(n, dtype=[('Cp', float), ('cf', float)])
    surface['Cp'] = (-1+x+pressure_error*np.sin(np.pi*x/.04))[::-1]; surface['cf'] = .001
    arrays = dict(xy=xy, U=u, lower=np.arange(n), upper=np.array([], dtype=int), turbulent=np.zeros(n, dtype=bool),
        wake_xy=np.array([[1., 1.1, 1.2], [wake_shift, wake_shift, wake_shift]]))
    layers = {e: {s: dict(transition=None, separation_regions=[]) for s in ('lower', 'upper')} for e in ('fore', 'aft')}
    result = dict(CL=.1, CD=.01, boundary_layers=layers)
    return dict(result=result, arrays={e: {k: v.copy() for k, v in arrays.items()} for e in ('fore', 'aft')},
        surfaces={e: surface.copy() for e in ('fore', 'aft')}, surface_s={e: x.tolist() for e in ('fore', 'aft')})


class ResolutionAnalysis(unittest.TestCase):
    def test_no_extrapolation_and_no_turbulent_shear_as_amplification(self):
        values, mask = supported_interpolation(np.array([0., 1., 2.]), np.array([0., 1., 90.]),
            np.arange(3), np.array([-.1, 0., .5, 1., 1.5, 2., 2.1]), np.array([False, False, True]))
        self.assertEqual(mask.tolist(), [False, True, True, True, False, False, False])
        self.assertEqual(values[2], .5)

    def test_all_fine_nodes_are_compared(self):
        states = {f: state(f) for f in (1, 2, 4)}; common = common_profiles(states)
        self.assertEqual(common['candidate_stations'], 9)
        self.assertEqual(common['common_supported_stations'], 9)
        self.assertEqual(common['common_laminar_stations'], 9)
        states[4]['surfaces']['fore']['Cp'][-2] += .05  # new fine node, absent on coarse/medium
        result = assess(states)
        self.assertGreater(result['criteria']['Cp']['value'], .049)
        self.assertEqual(result['numerical_screen'], 'fail')

    def test_laminar_mask_is_joint_across_all_grids(self):
        states = {f: state(f) for f in (1, 2, 4)}
        states[1]['arrays']['fore']['turbulent'][-1] = True
        common = common_profiles(states)
        self.assertEqual(common['common_laminar_stations'], 5)

    def test_identical_profiles_pass(self):
        result = assess({f: state(f) for f in (1, 2, 4)})
        self.assertEqual(result['numerical_screen'], 'pass')
        self.assertTrue(result['isolated_nose_refinement'])

    def test_missing_case_is_unknown_not_zero(self):
        result = assess({1: state(1), 2: state(2)})
        self.assertEqual(result['numerical_screen'], 'unknown')
        self.assertIsNone(result['isolated_nose_refinement'])
        self.assertEqual(result['pairs'], [])

    def test_wake_shift_confounds_attribution_without_invalidating_screen(self):
        result = assess({1: state(1), 2: state(2), 4: state(4, wake_shift=2e-5)})
        self.assertEqual(result['numerical_screen'], 'pass')
        self.assertFalse(result['isolated_nose_refinement'])

    def test_transition_width_is_physical_not_legacy_parameter(self):
        a, b = state(1), state(4)
        a['surface_s']['fore'] = (10*np.asarray(a['surface_s']['fore'])).tolist()
        a['result']['boundary_layers']['fore']['lower']['transition'] = dict(xy_over_cref=[.02, 0.])
        b['result']['boundary_layers']['fore']['lower']['transition'] = dict(xy_over_cref=[.028, 0.])
        result = compare_topology(a, b); row = result['sides'][0]
        self.assertAlmostEqual(row['physical_transition_shift'], .008)
        self.assertAlmostEqual(row['fine_cell_width'], .005)
        self.assertFalse(result['transition_within_fine_cell'])

    def test_coarse_topology_change_fails_even_when_finest_pair_matches(self):
        for element, side, changed in (
                ('fore', 'upper', dict(transition=dict(xy_over_cref=[.02, 0.]), separation_regions=[])),
                ('aft', 'lower', dict(transition=None, separation_regions=[dict(reaches_trailing_edge=True)]))):
            with self.subTest(element=element, side=side):
                states = {f: state(f) for f in (1, 2, 4)}
                states[1]['result']['boundary_layers'][element][side] = changed
                result = assess(states)
                self.assertFalse(result['pairs'][0]['topology']['same_topology'])
                self.assertTrue(result['pairs'][1]['topology']['same_topology'])
                self.assertFalse(result['topology']['same_topology'])
                self.assertEqual(result['numerical_screen'], 'fail')

    def test_coarse_transition_shift_is_reported_finest_shift_is_screened(self):
        states = {f: state(f) for f in (1, 2, 4)}
        for factor, position in ((1, .001), (2, .031), (4, .031)):
            states[factor]['result']['boundary_layers']['fore']['lower']['transition'] = dict(xy_over_cref=[position, 0.])
        result = assess(states)
        self.assertFalse(result['pairs'][0]['topology']['transition_within_fine_cell'])
        self.assertAlmostEqual(result['pairs'][0]['topology']['sides'][0]['physical_transition_shift'], .03)
        self.assertTrue(result['pairs'][1]['topology']['transition_within_fine_cell'])
        self.assertTrue(result['topology']['same_topology'])
        self.assertEqual(result['topology']['topology_screen_pairs'], [[1, 2], [2, 4]])
        self.assertEqual(result['topology']['transition_movement_screen_pair'], [2, 4])
        self.assertEqual(result['numerical_screen'], 'pass')

    def test_contraction_noise_floor_and_missing_support(self):
        self.assertTrue(contraction(1e-8, 5e-6, 1e-5))
        self.assertFalse(contraction(1e-8, 2e-5, 1e-5))
        self.assertTrue(contraction(.03, .02, 1e-5))
        self.assertFalse(contraction(.02, .03, 1e-5))
        self.assertIsNone(contraction(None, .01, 1e-5))

    def test_off_polygon_location_rejected(self):
        with self.assertRaisesRegex(ValueError, 'off the frozen polygon'):
            location_on_polygon(np.array([[0., 1.], [0., 0.]]), [.5, .01])


if __name__ == '__main__': unittest.main()
