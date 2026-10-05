"""Coverage, source support and audit contracts for the S414 diagnostic."""
import csv
import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT/'validation'))
import report_s414_lower_bucket as report


def rows(name):
    with (report.OUT/name).open() as stream:
        return list(csv.DictReader(stream))


class LowerBucketAnalysisTests(unittest.TestCase):
    def test_missing_reference_is_not_clamped_or_extrapolated(self):
        curve = [{'CL': 0, 'CD': .01}, {'CL': 1, 'CD': .02}]
        self.assertIsNone(report.reference_at(curve, -.001))
        self.assertIsNone(report.reference_at(curve, 1.001))
        self.assertAlmostEqual(report.reference_at(curve, .25), .0125)
        with self.assertRaises(ValueError):
            report.reference_at(curve[::-1], .25)

    def test_physical_taps_keep_identity_across_cases_coordinates_and_methods(self):
        pressure = rows('pressure_at_taps.csv')
        masks = []
        for case in ('clm0072', 'cl0056', 'cl0175', 'cl028'):
            for method in (('multifoil', 'mses') if case == 'cl028' else ('multifoil',)):
                for coordinate in ('table', 'plotted'):
                    mask = {(r['surface'], r['table_index']) for r in pressure if r['case'] == case
                        and r['method'] == method and r['coordinate'] == coordinate and r['included'] == 'True'}
                    self.assertEqual(len(mask), 87)
                    self.assertEqual(sum(s == 'fore_lower' for s, _ in mask), 35)
                    masks.append(mask)
        self.assertTrue(all(mask == masks[0] for mask in masks))
        excluded = rows('excluded_taps.csv')
        self.assertEqual(sum(r['reason'] == 'duplicate_endpoint' for r in excluded), 16)
        self.assertEqual(sum(r['reason'] == 'unsupported_source_center' for r in excluded), 12)

    def test_published_mses_gaps_remain_explicit(self):
        targets = rows('targets.csv')
        low = next(r for r in targets if float(r['target_cl']) < 0)
        self.assertEqual(low['mses_CD'], '')
        self.assertEqual(low['mses_error_CD_counts'], '')
        self.assertEqual(low['experiment_rounding_fully_supported'], 'False')
        self.assertEqual(low['experiment_rounding_CD_low'], '')
        self.assertEqual({r['case'] for r in rows('pressure_at_taps.csv') if r['method'] == 'mses'}, {'cl028'})

    def test_every_original_source_recovery_interval_is_preserved(self):
        with (report.SOURCE/'pressure_recovery_intervals.csv').open() as stream:
            reference = list(csv.DictReader(stream))
        jobs = json.loads((report.STUDY/'study.json').read_text())['jobs']
        mapping = {j['pressure_case']: j['pressure_figure'] for j in jobs}
        actual = [r for r in rows('recovery_intervals.csv') if r['method'] == 'multifoil']
        self.assertEqual(len(actual), len(reference))
        for row in actual:
            matches = [r for r in reference if r['figure'] == mapping[row['case']]
                and r['coordinate'] == row['coordinate']+'_x' and r['left_x'] == row['left_x'] and r['right_x'] == row['right_x']]
            self.assertEqual(len(matches), 1)
            self.assertAlmostEqual(float(row['experiment_delta_Cp']), float(matches[0]['delta_Cp']), places=14)
            self.assertAlmostEqual(float(row['experiment_slope']), float(matches[0]['slope']), places=12)

    def test_report_requires_full_audit_and_preserves_all_requests(self):
        summary = json.loads((report.OUT/'summary.json').read_text())
        receipt = json.loads((report.OUT/'residual_recheck.json').read_text())
        self.assertTrue(summary['full_state_audit_pass'] and receipt['complete'])
        self.assertEqual(len(receipt['states']), 4)
        self.assertEqual(summary['requested_evaluations'], 22)
        self.assertEqual(summary['failed_requests'], 0)
        self.assertEqual(summary['failed_internal_points'], 0)
        self.assertLessEqual(max(r['max_residual'] for r in receipt['states']), 1e-7)
        self.assertEqual(receipt['new_flow_solves'], 0)
        self.assertTrue(all(r['exactly_equal'] for r in receipt['controls']['array_comparisons']))

    def test_drag_contrasts_sum_to_total_without_adding_friction(self):
        diagnosis = json.loads((report.OUT/'diagnosis.json').read_text())
        for row in diagnosis['adjacent_lift_contrasts']:
            self.assertAlmostEqual(row['fore_delta_CD_counts']+row['aft_delta_CD_counts'], row['model_delta_CD_counts'], places=12)
        # Independent benchmark reference check at the primary discrepant condition.
        t = next(r for r in rows('targets.csv') if r['case'] == 'cl0056')
        expected = .006177+(.056-.0556)/(.0851-.0556)*(.006370-.006177)
        self.assertAlmostEqual(float(t['experiment_CD']), expected, places=14)
        self.assertGreater(float(t['error_CD_counts']), 50)


if __name__ == '__main__':
    unittest.main()
