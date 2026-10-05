"""Initialization, failure provenance and selected-state contracts; no CFD."""
from types import SimpleNamespace
import unittest

from fixed_lift_policy import InitializationPath, lift_seed, solve_with_initialization


class StatefulSolver:
    def __init__(self, start_alpha, fail_cold=False, residual=1e-9, slope=.1):
        self.start_alpha = start_alpha
        self.fail_cold = fail_cold
        self.residual = residual
        self.slope = slope
        self.history = []
        self.last_alpha = None
        self.coupling = {'state': None}

    def solve(self, alpha):
        self.history.append(alpha)
        self.last_alpha = alpha
        self.coupling['state'] = [alpha, len(self.history)]
        ok = not self.fail_cold or self.start_alpha == 0.
        self.last_solution = SimpleNamespace(alpha=alpha, cl=self.slope*alpha if ok else float('nan'),
            cd=.01 if ok else float('nan'), cm=0. if ok else float('nan'),
            converged=ok, residual=self.residual, diagnostics={})
        return self.last_solution


class InitializationPolicyTests(unittest.TestCase):
    def run_policy(self, factory=StatefulSolver, target=.1, center=1., **kwargs):
        return solve_with_initialization(factory, target, center, (-2., 4.),
            kwargs.pop('paths', [InitializationPath('cold', center), InitializationPath('warm', 0.)]),
            **kwargs)

    def test_first_success_stops_without_trying_other_path(self):
        created = []
        def factory(a):
            created.append(a)
            return StatefulSolver(a)
        out = self.run_policy(factory)
        self.assertTrue(out.feasible)
        self.assertEqual(created, [1.])
        self.assertEqual(out.selected_path, 'cold')
        self.assertLess(abs(out.solution.cl-.1), 2e-5)

    def test_failed_cold_flow_is_kept_before_warm_success(self):
        out = self.run_policy(lambda a: StatefulSolver(a, fail_cold=True))
        self.assertTrue(out.feasible)
        self.assertEqual(out.selected_path, 'warm')
        self.assertFalse(out.evaluations[0]['solution'].converged)
        self.assertFalse(out.attempts[0]['feasible'])
        self.assertTrue(out.attempts[1]['feasible'])
        self.assertEqual(out.selected_flow_index, len(out.evaluations)-1)

    def test_selected_left_endpoint_restores_entire_solver(self):
        out = self.run_policy(target=.095)
        self.assertEqual(len(out.evaluations), 2)
        self.assertEqual(out.selected_flow_index, 0)
        self.assertAlmostEqual(out.solution.alpha, .95)
        self.assertEqual(out.solver.last_alpha, out.solution.alpha)
        self.assertEqual(out.solver.coupling['state'], [.95, 1])
        self.assertEqual(out.solver.history, [.95])
        out.solver.solve(1.2)
        self.assertEqual(out.solution.alpha, .95)

    def test_failed_policy_has_no_selected_solver_or_force(self):
        out = self.run_policy(lambda a: StatefulSolver(a, residual=.01))
        self.assertFalse(out.feasible)
        self.assertIsNone(out.solution)
        self.assertIsNone(out.solver)
        self.assertIsNone(out.selected_flow_index)
        self.assertEqual(len(out.attempts), 2)
        self.assertEqual(len(out.evaluations), 2)
        self.assertTrue(all(not x['feasible'] for x in out.attempts))

    def test_warm_seed_is_logged_but_cannot_be_selected_as_target(self):
        seen = []
        out = self.run_policy(paths=[InitializationPath('replay', 0., (1.,))],
                              on_evaluation=lambda sol, context: seen.append(context))
        self.assertTrue(out.feasible)
        self.assertEqual(seen[0]['phase'], 'warm_seed')
        self.assertGreater(out.selected_flow_index, 0)
        self.assertEqual(len(seen), len(out.evaluations))

    def test_unbracketed_pair_recenters_using_computed_lift(self):
        out = self.run_policy(target=.15)
        self.assertTrue(out.feasible)
        self.assertEqual(len(out.attempts), 2)
        self.assertAlmostEqual(out.attempts[1]['center'], 1.5)
        self.assertEqual(out.attempts[1]['path'], 'cold')

    def test_invalid_branch_cannot_recenter_or_succeed(self):
        out = self.run_policy(lambda a: StatefulSolver(a, slope=-.1))
        self.assertFalse(out.feasible)
        self.assertEqual(len(out.attempts), 2)
        self.assertTrue(all(a['bracket'] == 0 for a in out.attempts))

    def test_callback_failure_propagates_without_retry(self):
        def broken(sol, context):
            raise RuntimeError('disk failure')
        with self.assertRaisesRegex(RuntimeError, 'disk failure'):
            self.run_policy(on_evaluation=broken)

    def test_seed_rejects_extrapolation_outside_domain(self):
        solver = StatefulSolver(0.)
        with self.assertRaisesRegex(RuntimeError, 'incidence domain'):
            lift_seed([solver.solve(0.), solver.solve(1.)], .4, (-2., 4.))

    def test_input_validation_precedes_any_solve(self):
        for kwargs in ({'half_width': 0.}, {'residual_tolerance': float('nan')},
                       {'max_brackets': 3}, {'paths': []},
                       {'paths': [InitializationPath('outside', -3.)]}):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                self.run_policy(**kwargs)


if __name__ == '__main__':
    unittest.main()
