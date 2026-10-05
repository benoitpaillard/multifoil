"""Run-deck guards for the S414 fixed-incidence transition-boundary paths (no flow solves)."""
import json
import math
from pathlib import Path
import shutil
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch
import numpy as np
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT/'validation'))
import run_s414_transition_boundary as run
from multifoil import Solution


class FakeSolver:
    """Minimal stand-in exposing only what the worker reads; fails at chosen incidences."""
    fail = ()

    def __init__(self, elements, **kwargs):
        self.kwargs = kwargs; self.amplification_jacobian = kwargs['amplification_jacobian']
        self.initialized = False; self.last_alpha = None; self.last_solution = None
        model = lambda name: SimpleNamespace(geom=SimpleNamespace(name=name),
            glob=SimpleNamespace(U=np.ones((4, 6))), foil=SimpleNamespace(x=np.zeros((2, 4))),
            isol=SimpleNamespace(xi=np.zeros(6)), vsol=SimpleNamespace(turb=np.zeros(6), Is=[[1, 0], [2, 3]]),
            wake=SimpleNamespace(x=np.array([[0., 1.], [0., 0.]]), s=np.zeros(2), t=np.zeros((2, 2)), N=2))
        self.models = [model('fore'), model('aft')]

    def solve(self, alpha):
        ok = not any(abs(alpha-a) < 1e-12 for a in self.fail)
        diag = dict(iterations=3, continuation=[dict(alpha_deg=alpha, converged=ok, CL=.1 if ok else math.nan)],
                    wall_seconds=.1, failure_reason='' if ok else 'Newton iteration limit reached',
                    boundary_layers={}, element_forces=dict(fore=dict(CD_wake=.004), aft=dict(CD_wake=.003)),
                    local_supersonic_detected={})
        sol = Solution(alpha, .1, -.05, .007, None, True, 1e-9, [], diag) if ok else \
            Solution(alpha, math.nan, math.nan, math.nan, None, False, .3, [], diag)
        if ok: self.initialized, self.last_alpha, self.last_solution = True, alpha, sol
        return sol


class TransitionBoundaryDeck(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.spec = run.specification()
        cls.archive = {s: json.loads((ROOT/run.ARCHIVE/j/'result.json').read_text()) for s, j in run.ARCHIVED_JOBS.items()}

    def test_bracket_is_the_archived_full_precision_pair(self):
        b = self.spec['bracket']; angles = self.spec['angles_deg']
        for result in self.archive.values():
            self.assertEqual([r['alpha_deg'] for r in result['records'][:2]], [b['lower_alpha_deg'], b['upper_alpha_deg']])
            self.assertEqual(result['initialization_paths'][0]['start_alpha'], b['archived_start_alpha_deg'])
        self.assertEqual((angles[0], angles[-1], len(angles)), (b['lower_alpha_deg'], b['upper_alpha_deg'], 5))
        np.testing.assert_allclose(np.diff(angles), (angles[-1]-angles[0])/4, rtol=0, atol=1e-14)
        # Full precision, not the rounded values quoted in the plan.
        self.assertNotIn(-2.56888, angles); self.assertNotIn(-2.46888, angles)

    def test_six_jobs_twenty_path_points_and_four_legacy_controls(self):
        jobs = self.spec['jobs']; b = self.spec['bracket']; angles = self.spec['angles_deg']
        self.assertEqual([j['name'] for j in jobs], ['nominal_legacy_endpoints', 'fore_source_cubic_legacy_endpoints',
            'nominal_corrected_increasing', 'nominal_corrected_decreasing',
            'fore_source_cubic_corrected_increasing', 'fore_source_cubic_corrected_decreasing'])
        path = [r for j in jobs for r in j['requests'] if r['role'] == 'path']
        self.assertEqual((len(path), self.spec['requested_path_points']), (20, 20))
        for job in jobs:
            self.assertEqual(job['start_alpha'], b['archived_start_alpha_deg']); self.assertEqual(job['wake_length'], 1.)
            requested = [r['alpha_deg'] for r in job['requests'] if r['role'] != 'endpoint_seed']
            if job['control']:
                self.assertEqual(job['amplification_jacobian'], 'legacy')
                self.assertEqual(requested, [b['lower_alpha_deg'], b['upper_alpha_deg']])
                self.assertEqual({r['role'] for r in job['requests']}, {'control_endpoint'})
            else:
                self.assertEqual(job['amplification_jacobian'], 'corrected_af')
                self.assertEqual(requested, angles if job['direction'] == 'increasing' else angles[::-1])
                seeds = [r for r in job['requests'] if r['role'] == 'endpoint_seed']
                self.assertEqual([s['alpha_deg'] for s in seeds],
                                 [] if job['direction'] == 'increasing' else [b['lower_alpha_deg']])
        self.assertEqual(sum(len(j['requests']) for j in jobs if j['control']), 4)

    def test_physical_and_numerical_settings_are_the_archived_ones(self):
        old = json.loads((ROOT/run.ARCHIVE/'study.json').read_text())
        for key in ('reynolds', 'mach', 'ncrit', 'panels', 'iterations', 'residual_tolerance', 'reference_chord'):
            self.assertEqual(self.spec[key], old[key])
        self.assertEqual((self.spec['ncrit'], self.spec['panels'], self.spec['iterations']), (9., [340, 240], 100))
        self.assertEqual(self.spec['expected_versions'], old['expected_versions'])
        self.assertEqual(self.spec['workers_limit'], 2)
        self.assertEqual(self.spec['equivalence_screens'], dict(CL=2e-5, CD_counts=.01, Cp=1e-4, transition_x_over_c=1e-4))
        for name in ('viscous.py', 'amplification_jacobian.py', 'vendor/mfoil.py', run.PROTOCOL, run.INTEGRATION,
                     'validation/run_s414_transition_boundary.py'):
            self.assertEqual(self.spec['sha256'][name], run.sha(ROOT/name))

    def test_only_the_fore_input_differs_between_geometries(self):
        nominal = run.geometry_elements('nominal', [340, 240]); cubic = run.geometry_elements('fore_source_cubic', [340, 240])
        np.testing.assert_array_equal(nominal[1][1], cubic[1][1])
        self.assertGreater(np.max(abs(nominal[0][1]-cubic[0][1])), 1e-6)
        with self.assertRaises(ValueError): run.geometry_elements('source_cubic', [340, 240])

    def test_changed_archive_or_unpassed_integration_gate_is_rejected(self):
        original = run.sha
        with patch.object(run, 'sha', side_effect=lambda p: 'changed' if str(p).endswith('cl0175_fore_source_cubic/flow/surface_001_fore.csv') else original(p)):
            with self.assertRaisesRegex(ValueError, 'Archived endpoint changed'): run.specification()
        real = run.read
        def failed(path):
            value = real(path)
            return dict(value, gate_passed=False) if str(path).endswith(run.INTEGRATION) else value
        with patch.object(run, 'read', failed), self.assertRaisesRegex(ValueError, 'V1 integration gate'):
            run.specification()

    def test_start_of_each_point_is_stated_without_calling_a_warm_start_cold(self):
        trail = lambda *steps: [dict(alpha_deg=a, converged=c) for a, c in steps]
        self.assertEqual(run.reached_by(False, None, -2.5, trail((-2.5, True))), 'cold solve at requested incidence')
        self.assertIn('continuation from cold solve at start incidence', run.reached_by(False, None, -2.6, trail((-2.5, True), (-2.6, True))))
        self.assertIn('after failed cold solve at start incidence', run.reached_by(False, None, -2.6, trail((-2.5, False), (-2.6, True))))
        self.assertIn('continuation from previous converged incidence -2.6', run.reached_by(True, -2.6, -2.5, trail((-2.5, True))))

    def run_fake(self, job, fail):
        FakeSolver.fail = fail
        with tempfile.TemporaryDirectory(dir=ROOT/'validation/results') as tmp, \
             patch.object(run, 'CoupledSolver', FakeSolver), patch.object(run, 'save_results', lambda *a, **k: None), \
             patch.object(run, 'transition_audit', lambda solver: []), \
             patch.object(run, 'geometry_elements', lambda *a: []):
            result = run.worker(Path(tmp), job, self.spec)
            states = sorted(p.name for p in (Path(tmp)/job['name']/'states').iterdir())
        return result, states

    def test_failed_seed_leaves_the_path_unevaluated_and_failures_have_no_force(self):
        b = self.spec['bracket']; angles = self.spec['angles_deg']
        job = next(j for j in self.spec['jobs'] if j['name'] == 'fore_source_cubic_corrected_decreasing')
        result, states = self.run_fake(job, fail=(b['upper_alpha_deg'],))
        self.assertEqual([r['role'] for r in result['records']], ['endpoint_seed', 'path'])
        self.assertEqual(result['error'], 'Endpoint seed failed; path not started')
        self.assertEqual([result['records'][1][k] for k in ('CL', 'CD', 'CM')], [None]*3)
        self.assertEqual(states, ['point_00_aft.npz', 'point_00_fore.npz'])
        result, _ = self.run_fake(job, fail=(b['lower_alpha_deg'],))
        self.assertEqual(len(result['records']), 1)
        # A mid-path failure keeps the direction and names the state it continued from.
        job = next(j for j in self.spec['jobs'] if j['name'] == 'nominal_corrected_increasing')
        result, states = self.run_fake(job, fail=(angles[2],))
        self.assertEqual([r['converged'] for r in result['records']], [True, True, False, True, True])
        self.assertEqual(result['records'][2]['CD'], None); self.assertEqual(result['error'], '')
        self.assertEqual(result['records'][3]['previous_converged_alpha_deg'], angles[1])
        self.assertEqual([r['approach'] for r in result['records']], ['cold']+['increasing']*4)
        self.assertEqual(len(states), 8)
        self.assertTrue(result['complete'])

    def test_control_gate_requires_identical_archived_exports(self):
        job = self.spec['jobs'][1]; archive = ROOT/job['archive']
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)/job['name']; shutil.copytree(archive/'flow', folder/'flow')
            (folder/'result.json').write_text(json.dumps(dict(records=[dict(converged=True), dict(converged=True)])))
            self.assertTrue(run.compare_control(folder, job)['passed'])
            path = folder/'flow/wake_001_aft.csv'; path.write_text(path.read_text().replace('0', '1', 1))
            outcome = run.compare_control(folder, job)
            self.assertFalse(outcome['passed'])
            self.assertEqual([c['item'] for c in outcome['comparisons'] if not c['identical']], ['flow/wake_001_aft.csv'])
            (folder/'result.json').write_text(json.dumps(dict(records=[dict(converged=True)])))
            self.assertFalse(run.compare_control(folder, job)['passed'])


if __name__ == '__main__':
    unittest.main()
