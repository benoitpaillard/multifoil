"""Numerical contracts for incompressible shape search, not flow validation."""
import copy
import json
import contextlib
import io
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
from multifoil import naca4, transform, read_dat
from shape_optimization import ShapeFamily, camber_mode, paired_score, pattern_search, check_config, run


def config():
    return json.loads((Path(__file__).parent/'examples/s702_incompressible_design.json').read_text())


class ShapeTests(unittest.TestCase):
    def test_frozen_frame_and_thickness_under_rotation(self):
        xy = naca4('0012', 201)
        modes = [dict(element='foil', support=[.05, .95], bounds=[-.01, .01])]
        changed = ShapeFamily([('foil', xy)], modes).deform([.003])[0][1]
        # Same x on the upper/lower surface receives exactly the same shift.
        np.testing.assert_allclose(changed[:, 0], xy[:, 0], atol=2e-16)
        np.testing.assert_allclose((changed-xy)[:, 1], (changed-xy)[::-1, 1], atol=2e-16)
        np.testing.assert_array_equal(changed[[0, 100, -1]], xy[[0, 100, -1]])
        moved = transform(xy, scale=.7, angle_deg=23, translation=[.6, -.1])
        deformed = ShapeFamily([('foil', moved)], modes).deform([.003])[0][1]
        np.testing.assert_allclose(deformed, transform(changed, scale=.7, angle_deg=23, translation=[.6, -.1]), atol=3e-16)
        np.testing.assert_array_equal(ShapeFamily([('foil', xy)], modes).deform([0.])[0][1], xy)

    def test_compact_shape_is_c2_and_bounded(self):
        self.assertEqual(camber_mode(np.array([.5]), .1, .9)[0], 1.)
        np.testing.assert_array_equal(camber_mode(np.array([-.1, .1, .9, 1.1]), .1, .9), 0.)
        h = 1e-6
        # C2 contact: deformation scales as h^3 near the fixed end.
        ratio = camber_mode(np.array([.1+h, .1+2*h]), .1, .9)
        self.assertAlmostEqual(ratio[1]/ratio[0], 8., delta=4e-5)

    def test_mach_zero_is_required(self):
        cfg = config()
        check_config(cfg)
        for value in (.1, None, '0'):
            cfg['mach'] = value
            with self.assertRaises(ValueError):
                check_config(cfg)

    def test_shipped_case_paths_and_bounded_geometry(self):
        cfg = config()
        root = Path(__file__).parent/'examples'
        elements = [(e['name'], read_dat(root/e['file'])) for e in cfg['elements']]
        family = ShapeFamily(elements, cfg['modes'])
        for p in (-.002, 0., .002):
            changed = family.deform([p])
            np.testing.assert_array_equal(changed[1][1], elements[1][1])
            np.testing.assert_array_equal(changed[0][1][[0, -1]], elements[0][1][[0, -1]])


class ObjectiveTests(unittest.TestCase):
    def setUp(self):
        self.cfg = config()
        self.base = [dict(scenario=s['name'], point=p['name'], feasible=True, CL=p['cl'], CD=.01, CM=-.1)
                     for s in self.cfg['scenarios'] for p in self.cfg['points']]

    def test_transition_tradeoff_is_not_hidden(self):
        trial = copy.deepcopy(self.base)
        for r in trial:
            r['CD'] += -.001 if r['scenario'] == 'ncrit9' else .00008
        score = paired_score(trial, self.base, self.cfg)
        self.assertTrue(score['feasible'])
        self.assertAlmostEqual(score['score'], .00008)

    def test_off_design_moment_and_missing_flow_rejected(self):
        trial = copy.deepcopy(self.base)
        trial[0]['CD'] += .0002
        self.assertFalse(paired_score(trial, self.base, self.cfg)['feasible'])
        trial = copy.deepcopy(self.base)
        trial[0]['CM'] -= .02
        self.assertFalse(paired_score(trial, self.base, self.cfg)['feasible'])
        self.assertFalse(paired_score(trial[:-1], self.base, self.cfg)['feasible'])
        trial = copy.deepcopy(self.base)
        trial[0]['feasible'] = False
        self.assertFalse(paired_score(trial, self.base, self.cfg)['feasible'])
        trial = copy.deepcopy(self.base)
        trial[0]['CL'] += .001
        self.assertFalse(paired_score(trial, self.base, self.cfg)['feasible'])

    def test_search_finds_known_bounded_minimum_with_unique_calls(self):
        calls = []
        def objective(p, number):
            calls.append(tuple(p))
            return dict(feasible=True, score=(p[0]-.5)**2-.25, parameters=list(p))
        result = pattern_search(objective, [dict(bounds=[-1., 1.], step=.5)], 10, .5, 1e-5)
        self.assertEqual(result['parameters'], [.5])
        self.assertEqual(result['status'], 'step_limit')
        self.assertEqual(len(calls), len(set(calls)))
        self.assertTrue(all(-1 <= p[0] <= 1 for p in calls))

    def test_budget_is_not_claimed_as_convergence(self):
        def objective(p, number):
            return dict(feasible=True, score=-p[0])
        result = pattern_search(objective, [dict(bounds=[-2., 2.], step=1.)], 3, .5, .01)
        self.assertEqual(result['status'], 'budget_exhausted')
        self.assertEqual(len(result['history']), 3)


class DeliveryTests(unittest.TestCase):
    """The exported shape must depend on the mandatory fine-grid decision."""
    def exercise(self, reject_fine):
        cfg = config()
        examples = Path(__file__).parent/'examples'
        for e in cfg['elements']:
            e['file'] = str((examples/e['file']).resolve())
        family = ShapeFamily([(e['name'], read_dat(e['file'])) for e in cfg['elements']], cfg['modes'])
        original, basis = family.elements[0][1], family.basis[0]
        fine_calls = []

        class AnalyticPool:
            def __init__(self, workers):
                pass
            def __enter__(self):
                return self
            def __exit__(self, *args):
                pass
            def map(self, function, tasks):
                rows = []
                for folder, raw, conf, panels, scenario, point in tasks:
                    amplitude = float(np.sum((raw[0][1]-original)*basis)/np.sum(basis*basis))
                    fine = panels == cfg['verification_panels']
                    if fine:
                        fine_calls.append(amplitude)
                    # Both grids independently require their own baseline:
                    # a constant fine-grid offset must cancel in paired changes.
                    cd = .01 + (.0003 if fine else 0.) - .1*amplitude
                    if reject_fine and fine:
                        cd += .2*amplitude  # improvement reverses at fine resolution
                    rows.append(dict(scenario=scenario['name'], point=point['name'],
                                     feasible=True, CL=point['cl'], CD=cd, CM=-.1))
                return rows

        with tempfile.TemporaryDirectory() as folder:
            folder = Path(folder)
            path = folder/'case.json'
            path.write_text(json.dumps(cfg))
            out = folder/'out'
            with patch('shape_optimization.ProcessPoolExecutor', AnalyticPool), contextlib.redirect_stdout(io.StringIO()):
                run(path, out, 1)
            summary = json.loads((out/'summary.json').read_text())
            self.assertGreater(summary['candidate_parameters'][0], 0)
            self.assertEqual(summary['numerical_screen_passed'], not reject_fine)
            self.assertEqual(len(fine_calls), 8)
            selected = read_dat(out/'selected_geometry/fore.dat')
            if reject_fine:
                np.testing.assert_allclose(selected, original, atol=1e-15)
                self.assertEqual(summary['delivered_parameters'], [0.])
            else:
                self.assertGreater(summary['delivered_parameters'][0], 0.)
                self.assertGreater(np.max(abs(selected-original)), 1e-4)

    def test_fine_grid_own_baseline_allows_verified_export(self):
        self.exercise(False)

    def test_reversal_on_fine_grid_restores_baseline_export(self):
        self.exercise(True)


if __name__ == '__main__':
    unittest.main()
