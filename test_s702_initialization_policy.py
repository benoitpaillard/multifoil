"""Frozen baseline and follow-up matrix contracts; no flow solves."""
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT/'validation'))
import run_s702_initialization_policy as run


class PolicyMatrixTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.spec = run.specification()

    def test_six_primary_conditions_and_separate_replay_control(self):
        primary = [j for j in self.spec['jobs'] if not j['control']]
        controls = [j for j in self.spec['jobs'] if j['control']]
        self.assertEqual({(j['flap'], j['target_cl']) for j in primary},
                         {(f, c) for f in (5, 10) for c in (.5, .8, 1.2)})
        self.assertEqual(len(primary), 6)
        self.assertEqual(len(controls), 1)
        self.assertEqual(controls[0]['seed_alphas'], [-.8402477019000163, -.24024770190001643])

    def test_same_ordered_policy_on_every_primary_condition(self):
        for job in self.spec['jobs']:
            paths = run.paths_for(job, self.spec)
            if job['control']:
                self.assertEqual([p.name for p in paths], ['archived_replay'])
            else:
                self.assertEqual([p.name for p in paths], ['cold_local', 'warm_anchor'])
                self.assertEqual(paths[0].start_alpha, job['center'])
                self.assertEqual(paths[1].start_alpha, -1.)
                self.assertTrue(all(not p.seed_alphas for p in paths))

    def test_original_physics_and_declared_budgets_are_preserved(self):
        s = self.spec
        self.assertEqual((s['reynolds'], s['mach'], s['ncrit'], s['panels']),
                         (700000., .07, 6.25, [260, 180]))
        self.assertEqual(s['max_evaluations'], 20)
        self.assertEqual(s['max_brackets'], 2)
        self.assertEqual(s['residual_tolerance'], 1e-7)
        for path, expected in s['sha256'].items():
            self.assertEqual(run.sha(ROOT/path), expected)


if __name__ == '__main__':
    unittest.main()
