"""Guards for the four-case derivative-only comparison; no new flow solves."""
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT/'validation'))
import run_s414_jacobian as run


class JacobianStudyTests(unittest.TestCase):
    def test_declared_pairs_change_only_derivative_mode(self):
        spec = run.specification(); old = run.read(ROOT/run.OLD/'study.json')
        self.assertEqual(len(spec['jobs']), 4)
        self.assertEqual({(j['target_cl'], j['amplification_jacobian']) for j in spec['jobs']},
            {(.056, 'legacy'), (.284, 'legacy'), (.056, 'corrected_af'), (.284, 'corrected_af')})
        for key in ('reynolds', 'mach', 'ncrit', 'panels', 'reference_chord', 'iterations',
                    'cl_tolerance', 'residual_tolerance', 'initial_half_width', 'alpha_bounds',
                    'max_evaluations', 'max_brackets', 'warm_start_alpha', 'entry_offset', 'expected_versions'):
            self.assertEqual(spec[key], old[key])
        for job in spec['jobs']:
            original = next(j for j in old['jobs'] if j['target_cl'] == job['target_cl'])
            for key in ('center', 'scenario', 'wake_length', 'target_cl'):
                self.assertEqual(job[key], original[key])
            self.assertEqual(job['control'], job['amplification_jacobian'] == 'legacy')

    def test_modified_verified_adapter_is_rejected(self):
        original = run.sha
        with patch.object(run, 'sha', side_effect=lambda p: 'changed' if str(p).endswith('/amplification_jacobian.py') else original(p)):
            with self.assertRaisesRegex(ValueError, 'Jacobian check input changed'):
                run.specification()

    def test_both_legacy_controls_must_pass(self):
        spec = dict(jobs=[dict(name='a', control=True, reference_result='one/result.json'),
                          dict(name='b', control=True, reference_result='two/result.json')],
                    acceptance_control_CD_counts=.01)
        with patch.object(run, 'compare_control', side_effect=[dict(passed=True), dict(passed=False)]):
            result = run.controls(Path('/tmp'), spec)
        self.assertFalse(result['passed']); self.assertEqual(len(result['jobs']), 2)

    def test_failed_target_cannot_have_finite_force(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); (root/'bad').mkdir()
            run.write(root/'bad/result.json', dict(feasible=False, CL=None, CD=.01, CM=None, alpha_deg=None))
            with self.assertRaisesRegex(ValueError, 'finite selected force'):
                run.summarize(root, dict(jobs=[dict(name='bad')]))


if __name__ == '__main__':
    unittest.main()
