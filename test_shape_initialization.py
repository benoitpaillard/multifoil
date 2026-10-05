"""Design integration contracts for selected states, failed paths and configuration."""
import copy
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace as NS
import unittest
from unittest.mock import patch

import numpy as np
from multifoil import naca4
from shape_optimization import check_config, flow_job
from test_shape_optimization import config


def policy_point():
    return dict(name='cl08', cl=.795, weight=1, start_alpha=0., alpha_bounds=[-2., 4.],
        initialization=dict(center=1., paths=[dict(name='cold', start_alpha=1.),
                                              dict(name='warm', start_alpha=0.)]))


class DesignInitializationTests(unittest.TestCase):
    def test_outside_domain_and_ambiguous_policy_rejected(self):
        cfg = config(); cfg['points'] = [policy_point()]
        check_config(cfg)
        for change in ('start', 'seed', 'center', 'duplicate', 'unknown', 'budget'):
            broken = copy.deepcopy(cfg); p = broken['points'][0]['initialization']
            if change == 'start': p['paths'][1]['start_alpha'] = -3
            if change == 'seed': p['paths'][1]['seed_alphas'] = [5]
            if change == 'center': p['center'] = 5
            if change == 'duplicate': p['paths'][1]['name'] = 'cold'
            if change == 'unknown': p['warm_start'] = 0
            if change == 'budget': p['max_brackets'] = 3
            with self.subTest(change=change), self.assertRaises(ValueError): check_config(broken)

    def exercise(self, fail_all=False, legacy=False):
        raw = [('foil', naca4('0012', 41))]
        point = policy_point()
        if legacy:
            point.pop('initialization'); point['alpha_bounds'] = [.95, 1.05]
        cfg = config()
        created, exports = [], []

        class Solver:
            def __init__(self, elements, **settings):
                self.start = settings['start_alpha']; created.append(self.start)
                self.elements = elements
                self.models = [NS(geom=NS(name='foil'), glob=NS(U=np.zeros((4, 1))),
                    foil=NS(x=np.zeros((2, 1))), isol=NS(xi=np.zeros(1)),
                    vsol=NS(turb=np.zeros(1), Is=[[], []]), wake=NS(x=np.zeros((2, 1))))]

            def solve(self, alpha):
                ok = not fail_all and self.start == 0.
                self.models[0].glob.U[:] = alpha
                return NS(alpha=alpha, cl=.7+.1*alpha if ok else float('nan'),
                    cd=.01 if ok else float('nan'), cm=0. if ok else float('nan'),
                    converged=ok, residual=1e-9 if ok else 1.,
                    surfaces=[dict(name=n, xy=x) for n, x in self.elements],
                    diagnostics=dict(boundary_layers={'foil': {'upper': {'separated_arclength_fraction': 0.}}},
                                     element_forces={'foil': {'CD_wake': .01}}))

        with tempfile.TemporaryDirectory() as tmp, patch('shape_optimization.CoupledSolver', Solver), \
                patch('shape_optimization.repanel', side_effect=lambda x, p: x), \
                patch('shape_optimization.save_results', side_effect=lambda f, s, e, c, m: exports.append(copy.deepcopy(c))):
            folder = Path(tmp)
            result = flow_job((folder, raw, cfg, [40], cfg['scenarios'][0], point))
            json.dumps(result, allow_nan=False)
            self.assertEqual(json.loads((folder/'result.json').read_text()), result)
            if fail_all:
                self.assertFalse(result['feasible'])
                self.assertTrue(all(result[k] is None for k in ('CL', 'CD', 'CM', 'alpha_deg')))
                self.assertNotIn('selected_flow_index', result)
                self.assertEqual(len(result['initialization_attempts']), 2)
                self.assertEqual(result['flow_evaluations'], 2)
            elif legacy:
                self.assertTrue(result['feasible'])
                self.assertEqual(result['selected_path'], 'legacy')
                self.assertEqual(created, [0.])
                self.assertEqual(result['selected_flow_index'], 0)
            else:
                self.assertTrue(result['feasible'])
                self.assertEqual(created, [1., 0.])
                self.assertEqual(result['selected_path'], 'warm')
                self.assertEqual(result['selected_flow_index'], 1)
                self.assertEqual(result['flow_evaluations'], 3)
                self.assertFalse(result['records'][0]['admissible'])
                self.assertIsNone(result['records'][0]['CL'])
                with np.load(folder/'foil_state.npz') as z:
                    np.testing.assert_allclose(z['U'], result['alpha_deg'])
                self.assertEqual(exports[-1]['initialization_history'],
                                 [r['initialization'] for r in result['records']])

    def test_warm_fallback_saves_selected_endpoint_not_last_solver(self):
        self.exercise()

    def test_failure_keeps_trace_without_finite_final_forces(self):
        self.exercise(fail_all=True)

    def test_existing_configuration_keeps_legacy_solve(self):
        self.exercise(legacy=True)


if __name__ == '__main__':
    unittest.main()
