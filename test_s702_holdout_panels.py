"""Protocol composition and failure coverage for the reserved panel check."""
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).parent/'validation'))
import run_s702_holdout_panels as run


class HoldoutPanelsTests(unittest.TestCase):
    def test_only_declared_panels_change_and_manifest_is_not_mutated(self):
        spec=run.specification()
        primary=[j for j in spec['jobs'] if not j['control']]
        self.assertEqual({(j['flap'],j['target_cl']) for j in primary},{(f,c) for f in (5,10) for c in (.8,1.2)})
        self.assertEqual(len(primary),4)
        for job in spec['jobs']:
            settings=run.effective_spec(job,spec)
            self.assertEqual(settings['panels'],[260,180] if job['control'] else [340,240])
            self.assertEqual(job['wake_length'],1.)
            self.assertEqual({k for k in settings if settings[k]!=spec[k]}, {'panels'} if job['control'] else set())
        self.assertEqual(spec['panels'],[340,240])
        self.assertEqual((spec['reynolds'],spec['mach'],spec['ncrit']),(700000.,.07,6.25))

    def test_worker_gets_the_declared_per_job_configuration(self):
        spec=run.specification()
        for job in (spec['jobs'][0],spec['jobs'][-1]):
            with patch.object(run.baseline,'worker',return_value='sentinel') as worker:
                self.assertEqual(run.worker(Path('/tmp'),job,spec),'sentinel')
                self.assertEqual(worker.call_args.args[2]['panels'],job['panels'])
        with self.assertRaisesRegex(ValueError,'Job settings'):
            run.effective_spec(dict(spec['jobs'][0],wake_length=2.),spec)

    def test_failed_fine_target_cannot_pass_by_omission(self):
        spec=run.specification()
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp)
            for job in spec['jobs']:
                folder=out/job['name'];folder.mkdir()
                r=run.read(run.ROOT/job['reference_result'])
                if job==spec['jobs'][0]:r.update(feasible=False,CL=None,CD=None,CM=None,alpha_deg=None)
                run.write(folder/'result.json',r)
            summary=run.summarize(out,spec)
            self.assertEqual(summary['feasible_jobs'],4)
            self.assertFalse(summary['execution_checks_pass'])
            self.assertIsNone(summary['jobs'][0]['difference_from_base_CD_counts'])


if __name__=='__main__':unittest.main()
