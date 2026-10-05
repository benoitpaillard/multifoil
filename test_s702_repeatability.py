"""Plan, restart and coverage contracts; no flow solves."""
import copy
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'validation'))
import run_s702_repeatability as run


class RepeatabilityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.spec=run.specification()

    def test_matrix_has_all_six_conditions_on_three_paths_and_one_warm_control(self):
        spec=self.spec
        self.assertEqual(len(spec['jobs']),11)
        self.assertEqual(spec['nominal_fixed_lift_targets'],19)
        expected={(f,cl) for f in (5,10) for cl in (.5,.8,1.2)}
        for path in ('cold','increasing','decreasing'):
            actual={(j['flap'],cl) for j in spec['jobs'] if j['path']==path for cl in j['targets']}
            self.assertEqual(actual,expected)
        for job in spec['jobs']:
            if job['path']=='decreasing':
                self.assertEqual(job['targets'],[1.2,.8,.5])

    def test_original_failed_targets_seed_from_computed_lift(self):
        for ref in self.spec['references']:
            if ref['original_failure']:
                result=run.read(ROOT/ref['original_failure'])
                self.assertFalse(result['feasible'])
                result['CD']=10000.
                for record in result['records']:
                    record['CD']=-999.
                self.assertEqual(ref['local_seed_alpha'],run.local_seed(result))

    def test_worker_passes_target_order_without_duplicate_output_keyword(self):
        job=next(j for j in self.spec['jobs'] if j['path']=='decreasing')
        with patch.object(run,'run_group',return_value={'complete_success':True}) as mock:
            run.worker(Path('/unused'),job,self.spec)
        args=mock.call_args.args
        self.assertNotIn('targets',args[1])
        self.assertEqual(args[2]['targets'],[1.2,.8,.5])
        self.assertTrue(args[2]['start_at_seed'])

    def test_warm_control_preserves_original_prefix_and_default_start(self):
        job=next(j for j in self.spec['jobs'] if j['path']=='warm_replay')
        with patch.object(run,'run_group') as mock:
            run.worker(Path('/unused'),job,self.spec)
        spec=mock.call_args.args[2]
        self.assertFalse(spec['start_at_seed'])
        self.assertEqual(spec['warm_seed_alphas'],[-0.8402477019000163,-0.24024770190001643])

    def test_prepare_rejects_changed_manifest_and_nonempty_unmanifested_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp)/'study'
            run.prepare(out,self.spec)
            run.prepare(out,self.spec)
            changed=copy.deepcopy(self.spec); changed['residual_tolerance']=1e-5
            with self.assertRaisesRegex(ValueError,'Study or sources changed'):
                run.prepare(out,changed)
            other=Path(tmp)/'other'; other.mkdir(); (other/'partial.txt').write_text('keep')
            with self.assertRaisesRegex(ValueError,'Nonempty'):
                run.prepare(other,self.spec)


if __name__=='__main__':
    unittest.main()
