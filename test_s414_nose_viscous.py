"""Run-deck guards for the S414 convex-nose viscous check (no flow solves)."""
import json
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT/'validation'))
import run_s414_nose_viscous as run
import s414_nose_geometry as nose
from viscous import new_model


class ConvexNoseDeck(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.spec = run.specification(); cls.jobs = {j['name']: j for j in cls.spec['jobs']}
        cls.old = json.loads((ROOT/run.ARCHIVE/'study.json').read_text())

    def test_six_jobs_one_exact_control_three_targets_two_paths(self):
        self.assertEqual(list(self.jobs), ['nominal_cl0056_control', 'convex_nose_cl0056', 'convex_nose_cl0175',
                                           'convex_nose_cl0284', 'convex_nose_increasing', 'convex_nose_decreasing'])
        self.assertEqual([j['control'] for j in self.spec['jobs']], [True]+[False]*5)
        self.assertEqual({j['amplification_jacobian'] for j in self.spec['jobs']}, {'corrected_af'})
        self.assertEqual([(j['scenario'], j.get('target_cl')) for j in self.spec['jobs'][:4]],
                         [('nominal', .056), ('convex_nose', .056), ('convex_nose', .175), ('convex_nose', .284)])
        angles = self.spec['path_angles_deg']
        self.assertEqual((len(angles), angles), (8, sorted(angles)))
        self.assertEqual([r['alpha_deg'] for r in self.jobs['convex_nose_increasing']['requests']], angles)
        self.assertEqual([r['alpha_deg'] for r in self.jobs['convex_nose_decreasing']['requests']], angles[::-1])
        # Increasing: cold at its first incidence. Decreasing: cold at the archived warm anchor, above every path incidence.
        self.assertEqual((self.jobs['convex_nose_increasing']['start_alpha'], self.jobs['convex_nose_decreasing']['start_alpha']),
                         (angles[0], self.old['warm_start_alpha']))
        self.assertGreater(self.old['warm_start_alpha'], angles[-1])
        t2 = json.loads((ROOT/run.T2/'study.json').read_text())
        self.assertEqual(set(angles), set(t2['nominal_new_incidences_deg']) | set(t2['tranche1_angles_deg']))

    def test_policy_centers_and_settings_are_the_archived_ones(self):
        archived = {j['name']: j for j in self.old['jobs']}
        self.assertEqual(self.jobs['nominal_cl0056_control']['center'], archived['cl0056_corrected_af']['center'])
        self.assertEqual(self.jobs['convex_nose_cl0056']['center'], archived['cl0056_corrected_af']['center'])
        self.assertEqual(self.jobs['convex_nose_cl0284']['center'], archived['cl028_corrected_af']['center'])
        for key in ('reynolds', 'mach', 'ncrit', 'panels', 'iterations', 'cl_tolerance', 'residual_tolerance', 'initial_half_width',
                    'alpha_bounds', 'max_evaluations', 'max_brackets', 'warm_start_alpha', 'entry_offset', 'expected_versions'):
            self.assertEqual(self.spec[key], self.old[key])
        result = json.loads((ROOT/run.ARCHIVE/'cl0056_corrected_af/result.json').read_text())
        from dataclasses import asdict
        declared = [json.loads(json.dumps(asdict(p))) for p in run.paths_for(self.jobs['nominal_cl0056_control'], self.spec)]
        self.assertEqual(declared, result['initialization_paths'])

    def test_geometry_is_the_frozen_admitted_input_with_the_nominal_aft(self):
        nominal = run.geometry_elements('nominal', [340, 240]); convex = run.geometry_elements('convex_nose', [340, 240])
        # The frozen file is the definition; regenerating it may differ in the last bits on another SciPy build.
        np.testing.assert_allclose(convex[0][1], nose.fore_input('convex_nose'), rtol=0, atol=1e-12)
        np.testing.assert_array_equal(convex[1][1], nominal[1][1])
        self.assertGreater(np.max(abs(convex[0][1]-nominal[0][1])), 1e-5)
        polygon = new_model('fore', convex[0][1], 1e6, .1, 9., (.25, 0.), 1e-5).foil.x.T[::-1]
        self.assertGreater(nose.turns(polygon)[polygon[1:-1, 0] <= nose.NOSE].min(), 0)
        with self.assertRaises(ValueError): run.geometry_elements('single_curve', [340, 240])
        with self.assertRaises(ValueError): run.geometry_elements('convex_nose', [480, 340])

    def test_changed_geometry_or_missing_admission_is_rejected(self):
        original = run.sha
        with patch.object(run, 'sha', side_effect=lambda p: 'changed' if str(p).endswith(run.FORE) else original(p)):
            with self.assertRaisesRegex(ValueError, 'Frozen convex-nose input changed'): run.specification()
        real = run.read
        def unadmitted(path):
            value = real(path)
            return dict(value, convex_nose_material=False) if str(path).endswith(run.GEOMETRY+'/verification.json') else value
        with patch.object(run, 'read', unadmitted), self.assertRaisesRegex(ValueError, 'not admitted'):
            run.specification()
        with patch.object(run, 'sha', side_effect=lambda p: 'changed' if str(p).endswith('cl0056_corrected_af/fore_state.npz') else original(p)):
            with self.assertRaisesRegex(ValueError, 'Archived control changed'): run.specification()

    def test_control_gate_needs_identical_arrays_incidence_and_path(self):
        job = self.jobs['nominal_cl0056_control']; reference = (ROOT/job['reference_result']).parent
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)/job['name']; folder.mkdir()
            for name in ('result.json', 'fore_state.npz', 'aft_state.npz'): shutil.copyfile(reference/name, folder/name)
            self.assertTrue(run.compare_control(folder, job, .01)['passed'])
            with np.load(reference/'aft_state.npz') as saved: arrays = {k: saved[k].copy() for k in saved.files}
            arrays['U'][0, 5] *= 1+1e-15; np.savez_compressed(folder/'aft_state.npz', **arrays)
            outcome = run.compare_control(folder, job, .01); self.assertFalse(outcome['passed'])
            self.assertEqual([(a['element'], a['field']) for a in outcome['array_comparisons'] if not a['exactly_equal']], [('aft', 'U')])
            shutil.copyfile(reference/'aft_state.npz', folder/'aft_state.npz')
            result = json.loads((reference/'result.json').read_text()); result['selected_path'] = 'warm_anchor'
            (folder/'result.json').write_text(json.dumps(result)); self.assertFalse(run.compare_control(folder, job, .01)['passed'])


if __name__ == '__main__':
    unittest.main()
