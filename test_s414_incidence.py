"""Keep measured-incidence evidence separate from matched-lift validation."""
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch
sys.path.insert(0, str(Path(__file__).resolve().parent/'validation'))
import run_s414_incidence as run


class IncidenceTests(unittest.TestCase):
    def solution(self, **kwargs):
        fields = dict(converged=True, cl=.15, cd=.007, cm=0., alpha=-3.07, residual=1e-9, diagnostics={})
        fields.update(kwargs); return SimpleNamespace(**fields)

    def test_fixed_angle_accepts_and_reports_lift_disagreement(self):
        job = dict(target_alpha=-3.07, reference_experimental_cl=.056)
        self.assertTrue(run.accepted_fixed_angle(self.solution(), job, dict(residual_tolerance=1e-7)))

    def test_wrong_incidence_failure_and_residual_are_rejected(self):
        for changed in (dict(alpha=-3.06), dict(residual=1e-6), dict(converged=False), dict(cd=float('nan'))):
            self.assertFalse(run.accepted_fixed_angle(self.solution(**changed), dict(target_alpha=-3.07), dict(residual_tolerance=1e-7)))

    def test_four_jobs_use_reported_angles_without_a_hidden_lift_target(self):
        spec = run.specification()
        self.assertEqual(len(spec['jobs']), 4)
        self.assertTrue(all(j['amplification_jacobian'] == 'corrected_af' for j in spec['jobs']))
        new = [j for j in spec['jobs'] if not j['control']]
        self.assertEqual([j['target_alpha'] for j in new], [-3.07, -1.04])
        self.assertEqual([j['reference_experimental_cl'] for j in new], [.056, .284])
        self.assertTrue(all('target_cl' not in j and j['comparison_basis']=='fixed_alpha' for j in new))
        self.assertTrue(all(j['comparison_basis']=='fixed_lift' and 'target_alpha' not in j for j in spec['jobs'] if j['control']))

    def test_prior_analysis_tampering_is_rejected(self):
        original = run.sha
        with patch.object(run, 'sha', side_effect=lambda p:'bad' if str(p).endswith(run.ANALYSIS+'/summary.json') else original(p)):
            with self.assertRaisesRegex(ValueError, 'Audited analysis changed'): run.specification()

    def test_failed_case_cannot_export_finite_force(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); (root/'bad').mkdir()
            run.write(root/'bad/result.json', dict(feasible=False, CL=None, CD=.01, CM=None, alpha_deg=None))
            with self.assertRaisesRegex(ValueError, 'finite selected force'):
                run.summarize(root, dict(jobs=[dict(name='bad')]))


if __name__ == '__main__': unittest.main()
