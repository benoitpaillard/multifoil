"""Checks for calibration leakage, missing cases and independent extraction."""
import hashlib
import json
from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT/'validation'))
from s702_sweep_metrics import score, select, leave_one_flap_out
from s702_ncrit_sweep import local_seed


class SweepMetricsTests(unittest.TestCase):
    def test_local_retry_uses_model_lift_and_rejects_reversed_branch(self):
        result = dict(target_cl=.8, CD_experiment=99., records=[
            dict(admissible=True,alpha_deg=-1.,CL=.7,CD=88.),
            dict(admissible=True,alpha_deg=0.,CL=.9,CD=-77.)])
        self.assertAlmostEqual(local_seed(result),-.5)
        result['records'][1]['CL']=.6
        with self.assertRaises(ValueError):
            local_seed(result)

    def row(self, n, flap, cl, error=0., feasible=True):
        return dict(ncrit=n, flap=flap, target_cl=cl, CD_experiment=.01,
                    CD=.01*(1+error) if feasible else None, feasible=feasible)

    def test_missing_difficult_case_cannot_win(self):
        rows = [self.row(4, 0, .5), self.row(4, 5, .5, feasible=False),
                self.row(9, 0, .5, .1), self.row(9, 5, .5, .1)]
        self.assertEqual(select(rows, [(0, .5), (5, .5)])['ncrit'], 9)

    def test_heldout_error_does_not_choose_its_own_parameter(self):
        rows = [self.row(4, 0, .5, 0), self.row(4, 5, .5, .8),
                self.row(9, 0, .5, .2), self.row(9, 5, .5, 0)]
        folds, overall = leave_one_flap_out(rows, [0, 5], [.5], [4, 9])
        self.assertEqual(folds[1]['heldout_flap'], 5)
        self.assertEqual(folds[1]['selected_ncrit'], 4)
        self.assertAlmostEqual(folds[1]['test_mape_percent'], 80)
        self.assertAlmostEqual(overall['mape_percent'], 50)

    def test_duplicate_case_is_not_extra_weight(self):
        row = self.row(6, 0, .5)
        with self.assertRaises(ValueError):
            score([row, row.copy()], [(0, .5)])

    def test_absolute_relative_errors_cannot_cancel(self):
        result = score([self.row(6, 0, .5, .2), self.row(6, 5, .5, -.2)], [(0, .5), (5, .5)])
        self.assertAlmostEqual(result['mape_percent'], 20)
        self.assertAlmostEqual(result['rmse_counts'], 20)


class HoldoutExtractionTests(unittest.TestCase):
    def test_pinned_source_and_export_hashes(self):
        folder = ROOT/'validation/data/s702/re700k_holdout'
        meta = json.loads((folder/'source.json').read_text())
        pdf = folder.parent/'NASA-CR-20220014964.pdf'
        self.assertEqual(hashlib.sha256(pdf.read_bytes()).hexdigest(), meta['pdf_sha256'])
        for name, expected in meta['sha256'].items():
            self.assertEqual(hashlib.sha256((folder/name).read_bytes()).hexdigest(), expected)
        self.assertEqual((meta['reynolds'], meta['mach'], meta['flaps']), (700000., .07, [0, 5, 10]))

    def test_drag_markers_match_independent_lift_curve(self):
        # Independent drawing paths provide CL twice. This catches swapped
        # marker series and axes, not merely a round-trip through one formula.
        for flap, expected_count in [(0, 24), (5, 23), (10, 24)]:
            folder = ROOT/f'validation/data/s702/re700k_holdout/d{flap:+d}'
            drag = np.genfromtxt(folder/'experiment_drag.csv', delimiter=',', names=True)
            lift = np.genfromtxt(folder/'experiment_lift.csv', delimiter=',', names=True)
            self.assertEqual(len(drag), expected_count)
            mismatch = np.min(abs(drag['CL'][:, None]-lift['CL'][None, :]), axis=1)
            self.assertLess(float(mismatch.max()), 1e-4)
            self.assertLess(drag['CL'].min(), .5)
            self.assertGreater(drag['CL'].max(), 1.2)


if __name__ == '__main__':
    unittest.main()
