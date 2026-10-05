"""Combined-policy composition and coverage checks; no flow simulations."""
import copy
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).parent/'validation'))
import run_design_entry_seed as run
from shape_optimization import initialization_settings
from fixed_lift_policy import solve_with_initialization


class CombinedPolicyTests(unittest.TestCase):
    def test_same_ordered_rule_covers_six_shapes_and_an_independent_repeat(self):
        spec = run.specification(True)
        self.assertEqual(len(spec['jobs']), 7)
        self.assertEqual({(j['flap'], j['amplitude']) for j in spec['jobs']},
                         {(f, a) for f in (5, 10) for a in (-.001, 0., .001)})
        for job in spec['jobs']:
            point = run.configuration(job, spec)['points'][0]
            policy = initialization_settings(point)
            self.assertEqual([p.name for p in policy['paths']], ['cold_local', 'warm_anchor', 'warm_entry'])
            self.assertEqual(policy['paths'][0].start_alpha, job['center'])
            self.assertEqual(policy['paths'][1].start_alpha, -1.)
            self.assertEqual(policy['paths'][2].seed_alphas, (job['center']-.4,))
        baseline = next(j for j in spec['jobs'] if j['name']=='d10_baseline_combined')
        repeat = next(j for j in spec['jobs'] if 'repeat_of' in j)
        self.assertEqual(run.configuration(baseline,spec),run.configuration(repeat,spec))
        self.assertNotEqual(baseline['name'], repeat['name'])

    def test_third_path_keeps_failed_prefix_and_gets_a_fresh_solver(self):
        spec = run.specification(True)
        job = dict(spec['jobs'][0], center=1., seed_alpha=.6)
        point = run.configuration(job, spec)['points'][0]
        created = []
        class Solver:
            def __init__(self, start):
                self.start = start; self.entered = False; self.history = []
            def solve(self, alpha):
                self.history.append(alpha)
                if self.start == -1. and abs(alpha-.6) < 1e-12: self.entered = True
                return SimpleNamespace(alpha=alpha, cl=.1*alpha if self.entered else float('nan'),
                    cd=.01 if self.entered else float('nan'), cm=0., converged=self.entered,
                    residual=1e-9 if self.entered else 1., diagnostics={})
        def factory(start):
            solver=Solver(start);created.append(solver);return solver
        result=solve_with_initialization(factory,.1,alpha_bounds=point['alpha_bounds'], **initialization_settings(point))
        self.assertTrue(result.feasible)
        self.assertEqual(result.selected_path,'warm_entry')
        self.assertEqual(len(created),3)
        self.assertTrue(all(not a['feasible'] for a in result.attempts[:2]))
        self.assertFalse(result.evaluations[0]['solution'].converged)
        self.assertFalse(result.evaluations[1]['solution'].converged)
        self.assertEqual(result.evaluations[2]['phase'],'warm_seed')
        self.assertEqual(result.selected_flow_index,5)

    def test_missing_repeat_coverage_and_repeat_disagreement_cannot_pass(self):
        spec = run.specification(True)
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp)
            for job in spec['jobs']:
                folder=out/job['name'];folder.mkdir()
                value=run.read(run.ROOT/job['reference_result'])
                run.write_json(folder/'result.json',value)
            self.assertTrue(run.summarize(out,spec)['regression_pass'])
            repeat=next(j for j in spec['jobs'] if 'repeat_of' in j)
            path=out/repeat['name']/'result.json'; original=run.read(path)
            failed=dict(original,feasible=False,CL=None,CD=None,CM=None)
            run.write_json(path,failed)
            summary=run.summarize(out,spec)
            self.assertFalse(summary['regression_pass'])
            self.assertEqual(summary['passed_jobs'],6)
            self.assertIsNone(summary['repeats'][0]['difference_CD_counts'])
            run.write_json(path,original)
            # Each comparison to the legacy reference can pass while repeats
            # disagree with one another by more than the declared tolerance.
            other=out/repeat['repeat_of']/'result.json'
            run.write_json(path,dict(original,CD=original['CD']+.8e-6))
            run.write_json(other,dict(run.read(other),CD=original['CD']-.8e-6))
            summary=run.summarize(out,spec)
            self.assertEqual(summary['passed_jobs'],7)
            self.assertFalse(summary['regression_pass'])
            self.assertFalse(summary['repeats'][0]['passed'])

    def test_monitor_routes_to_the_requested_campaign(self):
        with patch.object(run,'remote_status',return_value=2) as monitor, \
                patch('sys.argv',['run_design_entry_seed.py','--combined','--remote-status','--watch']):
            self.assertEqual(run.main(),2)
            monitor.assert_called_once_with(True,False,True)


if __name__ == '__main__':
    unittest.main()
