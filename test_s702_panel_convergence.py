"""Third-grid reference integrity, coverage and fixed-protocol checks; no CFD."""
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).parent/'validation'))
import run_s702_panel_convergence as run


class PanelConvergenceTests(unittest.TestCase):
    def test_same_conditions_and_incidence_paths_with_correct_grid_references(self):
        spec=run.specification();prior=run.previous.specification()
        primary=[j for j in spec['jobs'] if not j['control']]
        self.assertEqual({(j['flap'],j['target_cl']) for j in primary},{(f,c) for f in (5,10) for c in (.8,1.2)})
        self.assertEqual(len(primary),4)
        for job in spec['jobs']:
            old=next(j for j in prior['jobs'] if not j['control'] and j['flap']==job['flap'] and j['target_cl']==job['target_cl'])
            self.assertEqual(job['reference_result'],str(Path(run.FINE)/old['name']/'result.json'))
            self.assertEqual(run.flow.paths_for(job,spec),run.flow.paths_for(old,prior))
            self.assertEqual(run.previous.effective_spec(job,spec)['panels'],[340,240] if job['control'] else [440,320])
        for k in ('reynolds','mach','ncrit','wake_length','alpha_bounds','iterations','cl_tolerance','residual_tolerance'):
            self.assertEqual(spec[k],prior[k])
        self.assertTrue(spec['last_automatic_grid_extension'])

    def test_modified_fine_reference_is_rejected(self):
        spec=run.specification();path=run.ROOT/spec['jobs'][0]['reference_result'];real_sha=run.sha
        with patch.object(run,'sha',side_effect=lambda p:'changed' if Path(p)==path else real_sha(p)):
            with self.assertRaisesRegex(ValueError,'Audited fine-grid input changed'):run.specification()

    def test_failed_target_has_no_delta_or_passing_coverage(self):
        spec=run.specification()
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp)
            for job in spec['jobs']:
                folder=out/job['name'];folder.mkdir()
                r=run.read(run.ROOT/job['reference_result'])
                if job==spec['jobs'][0]:r.update(feasible=False,CL=None,CD=None,CM=None,alpha_deg=None)
                run.write(folder/'result.json',r)
            s=run.summarize(out,spec)
            self.assertEqual(s['feasible_jobs'],4)
            self.assertFalse(s['execution_checks_pass'])
            self.assertIsNone(s['jobs'][0]['difference_from_fine_CD_counts'])


if __name__=='__main__':unittest.main()
