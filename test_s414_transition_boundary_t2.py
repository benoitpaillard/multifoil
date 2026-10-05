"""Run-deck guards for tranche 2 of the S414 transition-boundary study (no flow solves)."""
import json
from pathlib import Path
import sys
import unittest
from unittest.mock import patch
import numpy as np
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT/'validation'))
import run_s414_transition_boundary_t2 as run


class TrancheTwoDeck(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.spec = run.specification(); cls.jobs = {j['name']: j for j in cls.spec['jobs']}
        cls.angles = cls.spec['tranche1_angles_deg']
        cls.selection = json.loads((ROOT/run.ANALYSIS/'tranche2_selection.json').read_text())

    def test_eighteen_evaluations_use_exactly_the_forty_two_point_allowance(self):
        kinds = [j['kind'] for j in self.spec['jobs']]
        self.assertEqual((kinds.count('restore_control'), kinds.count('continuation'), kinds.count('cold_check')), (2, 14, 2))
        self.assertEqual(self.spec['requested_points']+self.spec['tranche1_requested_points'], self.spec['allowance'])
        self.assertEqual(self.spec['allowance'], 42)
        self.assertEqual({j['amplification_jacobian'] for j in self.spec['jobs']}, {'corrected_af'})
        self.assertEqual([j['kind'] for j in self.spec['jobs'][:2]], ['restore_control']*2)

    def test_cubic_incidences_are_the_audited_selection_at_full_precision(self):
        chosen = sorted(l['alpha_deg'] for l in self.selection['fore_source_cubic']['new_locations'])
        self.assertEqual(self.spec['cubic_new_incidences_deg'], chosen)
        for k, alpha in enumerate(chosen):
            for direction in ('increasing', 'decreasing'):
                self.assertEqual(self.jobs[f'fore_source_cubic_{direction}_{k}']['alpha_deg'], alpha)
        self.assertTrue(self.angles[0] < chosen[0] < self.angles[1] < chosen[1] < chosen[2] < self.angles[2] < chosen[3] < self.angles[3])

    def test_bases_are_the_nearest_converged_states_on_each_directions_own_path(self):
        a = self.angles; base = lambda name: (self.jobs[name]['base']['job'], self.jobs[name]['base']['alpha_deg'])
        up, down = 'fore_source_cubic_corrected_increasing', 'fore_source_cubic_corrected_decreasing'
        self.assertEqual([base(f'fore_source_cubic_increasing_{k}') for k in range(4)], [(up, a[0]), (up, a[1]), (up, a[1]), (up, a[2])])
        # The decreasing path failed at a[2] in tranche 1, so its next converged state above is used.
        self.assertEqual([base(f'fore_source_cubic_decreasing_{k}') for k in range(4)], [(down, a[1]), (down, a[3]), (down, a[3]), (down, a[3])])
        for job in self.spec['jobs']:
            if job['kind'] == 'cold_check': self.assertNotIn('base', job); continue
            for name in job['base']['state_files'].values():
                self.assertEqual(self.spec['sha256'][name], run.sha(ROOT/name))

    def test_nominal_fallback_bracket_midpoint_and_quarter_points(self):
        b = self.spec['nominal_fallback_bracket']; lo, hi = b['lower_alpha_deg'], b['upper_alpha_deg']
        self.assertEqual((lo, hi), (-3.07, self.angles[0])); self.assertGreater(b['fore_lower_transition_change'], .05)
        np.testing.assert_allclose(self.spec['nominal_new_incidences_deg'], [lo+f*(hi-lo) for f in (.25, .5, .75)], rtol=0, atol=1e-14)
        for k in range(3):
            self.assertEqual(self.jobs[f'nominal_increasing_{k}']['base']['job'], run.ANCHOR_JOB)
            self.assertEqual(self.jobs[f'nominal_increasing_{k}']['base']['study'], run.ANCHOR)
            self.assertEqual(self.jobs[f'nominal_decreasing_{k}']['base']['job'], 'nominal_corrected_decreasing')
            self.assertEqual(self.jobs[f'nominal_decreasing_{k}']['base']['alpha_deg'], hi)

    def test_cold_checks_and_the_recorded_deviation(self):
        cubic, nominal = self.jobs['fore_source_cubic_cold_check'], self.jobs['nominal_cold_check']
        self.assertEqual(cubic['alpha_deg'], self.angles[1])
        self.assertEqual(cubic['deviation']['declared_alpha_deg'], self.selection['fore_source_cubic']['cold_check']['alpha_deg'])
        self.assertEqual(cubic['deviation']['declared_alpha_deg'], self.angles[2])
        self.assertIn('user', cubic['deviation']['approved_by'])
        self.assertEqual(nominal['alpha_deg'], self.spec['nominal_new_incidences_deg'][1]); self.assertNotIn('deviation', nominal)

    def test_changed_base_state_or_incomplete_audit_is_rejected(self):
        original = run.sha
        target = self.jobs['fore_source_cubic_increasing_1']['base']['state_files']['fore']
        with patch.object(run, 'sha', side_effect=lambda p: 'changed' if str(p).endswith(target) else original(p)):
            with self.assertRaisesRegex(ValueError, 'Audited tranche-1 file changed'): run.specification()
        real = run.read
        def incomplete(path):
            value = real(path)
            return dict(value, complete=False) if str(path).endswith(run.ANALYSIS+'/residual_recheck.json') else value
        with patch.object(run, 'read', incomplete), self.assertRaisesRegex(ValueError, 'not completely audited'):
            run.specification()

    def test_restore_gate_applies_the_declared_screens(self):
        job = self.jobs['nominal_restore_control']; screens = self.spec['equivalence_screens']
        ref = ROOT/job['reference']['study']/job['reference']['job']
        old = next(r for r in json.loads((ref/'result.json').read_text())['records'] if r['index'] == 1)
        import shutil, tempfile
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)/job['name']; (folder/'flow').mkdir(parents=True)
            for e in ('fore', 'aft'): shutil.copyfile(ref/f'flow/surface_001_{e}.csv', folder/f'flow/surface_000_{e}.csv')
            def outcome(shift):
                (folder/'result.json').write_text(json.dumps(dict(error='', records=[dict(old, CL=old['CL']+shift)])))
                return run.compare_reference(folder, job, screens)
            self.assertTrue(outcome(0.)['passed']); self.assertTrue(outcome(1e-5)['passed']); self.assertFalse(outcome(3e-5)['passed'])
            (folder/'result.json').write_text(json.dumps(dict(error='Base state not admitted: x', records=[])))
            self.assertFalse(run.compare_reference(folder, job, screens)['passed'])


if __name__ == '__main__':
    unittest.main()
