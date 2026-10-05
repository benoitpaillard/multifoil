"""Tests for diagnostic validity; the original growth-Jacobian defect is retained."""
import ast
import csv
import inspect
import json
from pathlib import Path
import sys
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT/'validation'))
import report_s414_fore_nose as report
from report_s702_wake_extent import thermo


class ForeNoseTests(unittest.TestCase):
    def test_integer_flags_are_not_used_as_bitwise_boolean_masks(self):
        a = {'lower': np.array([3, 2, 1, 0]), 'turbulent': np.array([1, 1, 0, 0])}
        np.testing.assert_array_equal(report.laminar_indices(a), [3, 2])
        a['turbulent'] = np.array([1, 0, 1, 0])
        with self.assertRaises(ValueError): report.laminar_indices(a)

    def test_corrected_derivative_matches_finite_difference_without_changing_value(self):
        spec = json.loads((report.STUDY/'study.json').read_text())
        p = thermo(spec)
        with np.load(report.STUDY/'cl0056_nominal/fore_state.npz') as a:
            j = report.laminar_indices(a)[-1]; u = a['U'][:, j].copy()
        before = u.copy(); original = report.mf.get_damp(u, p)[0]
        rows = report.rate_derivatives(u, p)
        self.assertGreater(max(r['scaled_error'] for r in rows), .3)
        self.assertLess(max(r['corrected_scaled_error'] for r in rows), 2e-5)
        self.assertEqual(report.mf.get_damp(u, p)[0], original)
        np.testing.assert_array_equal(u, before)

    def test_transition_bracket_rejects_unphysical_root(self):
        spec = json.loads((report.STUDY/'study.json').read_text())
        p = thermo(spec)
        u = np.array([[.0001, .0001], [.0002, .0002], [0., .03], [.8, .8]])
        with self.assertRaises(ValueError): report.transition_root(u, np.array([.01, .011]), p)

    def test_turbulent_shear_is_never_exported_as_amplification(self):
        with (report.OUT/'fore_lower_stations.csv').open() as stream: rows = list(csv.DictReader(stream))
        self.assertTrue(any(r['laminar'] == 'False' for r in rows))
        for r in rows:
            if r['laminar'] == 'False':
                self.assertEqual(r['N'], ''); self.assertEqual(r['dN_dxi'], '')
                self.assertNotEqual(r['shear'], '')
            else:
                self.assertEqual(r['shear'], ''); self.assertNotEqual(r['N'], '')

    def test_failed_derivative_check_and_successful_balance_checks_remain_distinct(self):
        s = json.loads((report.OUT/'summary.json').read_text())
        self.assertFalse(s['original_growth_derivative_check_pass'])
        self.assertTrue(s['diagnostic_corrected_derivative_check_pass'])
        self.assertTrue(s['amplification_balances_and_transition_roots_pass'])
        self.assertEqual(s['new_nonlinear_flow_solves'], 0)
        self.assertEqual(s['archived_incidence_inviscid_evaluations'], 3)
        for c in s['cases']:
            self.assertTrue(c['state_and_flags_restored_exactly'])
            self.assertTrue(c['frozen_march_same_transition_cell'])
            self.assertLess(c['maximum_amplification_balance'], 1e-7)
        self.assertIsNone(s['cases'][0]['N_at_x0034'])

    def test_no_nonlinear_solver_or_vendor_mutation_is_called(self):
        tree = ast.parse(inspect.getsource(report))
        methods = {n.func.attr for n in ast.walk(tree) if isinstance(n, ast.Call) and isinstance(n.func, ast.Attribute)}
        self.assertFalse({'solve', '_newton', '_initialize', 'solve_viscous'} & methods)
        assignments = [n for n in ast.walk(tree) if isinstance(n, ast.Assign)]
        self.assertFalse(any(isinstance(t, ast.Attribute) and isinstance(t.value, ast.Name)
            and t.value.id == 'mf' for n in assignments for t in n.targets))


if __name__ == '__main__': unittest.main()
