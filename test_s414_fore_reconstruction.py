"""Guards for an isolated fore reconstruction and unchanged matched-lift protocol."""
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT/'validation'))
import run_s414_fore_reconstruction as run
from viscous import new_model, te_info


class ForeReconstructionTests(unittest.TestCase):
    def test_only_fore_input_changes_and_matches_existing_cubic(self):
        nominal = run.geometry_elements('nominal', [340, 240])
        changed = run.geometry_elements('fore_source_cubic', [340, 240])
        cubic = run.original_geometry('source_cubic', [340, 240])
        np.testing.assert_array_equal(changed[1][1], nominal[1][1])
        np.testing.assert_array_equal(changed[0][1], cubic[0][1])
        self.assertGreater(np.max(abs(changed[0][1]-nominal[0][1])), 1e-6)
        self.assertEqual(changed[0][1].shape, nominal[0][1].shape)
        for (name, a), (_, b) in zip(nominal, changed):
            # Independent interpolators differ by one ulp at a source endpoint.
            np.testing.assert_allclose(a[[0, -1]], b[[0, -1]], rtol=0, atol=2e-15)
            self.assertGreater(np.linalg.norm(a[0]-a[-1]), 1e-5)

    def test_final_aft_geometry_and_finite_te_are_preserved(self):
        saved = {}
        for scenario in run.SCENARIOS:
            saved[scenario] = [new_model(n, xy, 1e6, .1, 9., (.25, 0.), 1e-5)
                               for n, xy in run.geometry_elements(scenario, [340, 240])]
            for m in saved[scenario]: self.assertGreater(te_info(m.foil.x)[1], 1e-5)
        np.testing.assert_array_equal(saved['nominal'][1].foil.x, saved['fore_source_cubic'][1].foil.x)
        self.assertEqual(saved['nominal'][0].foil.N, saved['fore_source_cubic'][0].foil.N)

    def test_same_two_lifts_settings_and_initialization(self):
        spec = run.specification(); old = run.read(ROOT/run.OLD/'study.json')
        self.assertEqual(len(spec['jobs']), 4)
        self.assertEqual({(j['target_cl'], j['scenario']) for j in spec['jobs']},
                         {(.056, s) for s in run.SCENARIOS} | {(.175, s) for s in run.SCENARIOS})
        for key in ('reynolds', 'mach', 'ncrit', 'panels', 'iterations', 'cl_tolerance',
                    'residual_tolerance', 'initial_half_width', 'alpha_bounds', 'max_evaluations',
                    'max_brackets', 'warm_start_alpha', 'entry_offset', 'expected_versions'):
            self.assertEqual(spec[key], old[key])
        for job in spec['jobs']:
            original = next(j for j in old['jobs'] if j['target_cl'] == job['target_cl'])
            self.assertEqual(job['amplification_jacobian'], 'legacy')
            self.assertEqual(job['control'], job['scenario'] == 'nominal')
            for key in ('center', 'target_cl', 'wake_length'): self.assertEqual(job[key], original[key])

    def test_changed_archived_control_is_rejected(self):
        original = run.sha
        with patch.object(run, 'sha', side_effect=lambda p: 'changed' if str(p).endswith('cl0175_nominal/fore_state.npz') else original(p)):
            with self.assertRaisesRegex(ValueError, 'Archived control changed'): run.specification()

    def test_both_controls_must_pass(self):
        spec = dict(jobs=[dict(name='a', control=True, reference_result='a/result.json'),
                          dict(name='b', control=True, reference_result='b/result.json')], acceptance_control_CD_counts=.01)
        with patch.object(run, 'compare_control', side_effect=[dict(passed=True), dict(passed=False)]):
            self.assertFalse(run.controls(Path('/tmp'), spec)['passed'])

    def test_failed_target_cannot_have_finite_force(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); (root/'bad').mkdir()
            run.write(root/'bad/result.json', dict(feasible=False, CL=None, CD=.01, CM=None, alpha_deg=None))
            with self.assertRaisesRegex(ValueError, 'finite selected force'):
                run.summarize(root, dict(jobs=[dict(name='bad')]))

    def test_unknown_geometry_is_rejected(self):
        with self.assertRaises(ValueError): run.geometry_elements('source_cubic', [340, 240])


if __name__ == '__main__': unittest.main()
