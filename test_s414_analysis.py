"""Scientific comparison contracts for the S414 saved-state report."""
import ast
import csv
import inspect
import json
from pathlib import Path
import sys
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'validation'))
import report_s414_reconstruction as report
import audit_s414_state as audit


def rows(path):
    with path.open() as stream:return list(csv.DictReader(stream))


class S414AnalysisTests(unittest.TestCase):
    def test_lower_endpoint_repeats_do_not_shift_table_associations(self):
        ex=rows(report.QUAL/'pressure/experiment.csv')
        coords=rows(report.QUAL/'pressure/tap_coordinate_checks.csv')
        associated=report.associate_taps(ex,coords)
        self.assertEqual(len(associated),188)
        self.assertEqual(sum(r['table_index'] is None for r in associated),8)
        for case in report.CASES:
            lower=[r for r in associated if r['case']==case and r['surface']=='fore_lower']
            self.assertIsNone(lower[0]['table_index'])
            self.assertEqual((int(lower[1]['index']),lower[1]['table_index']),(1,0))
            self.assertIsNone(lower[-1]['table_index'])

    def test_common_support_uses_every_curve_and_both_coordinates(self):
        tap=dict(endpoint_duplicate=0,surface='aft_upper',table_index=2,table_x=.75,x_over_c=.751)
        curves={'nominal':(np.array([.7,1.]),np.zeros(2)),'shift':(np.array([.7505,1.]),np.zeros(2))}
        self.assertEqual(report.common_support(tap,curves,set()),'outside_common_curve_support')
        curves['shift']=curves['nominal']
        self.assertEqual(report.common_support(tap,curves,set()),'')
        self.assertEqual(report.common_support(tap,curves,{('aft','upper',2)}),'unsupported_source_center')
        tap['endpoint_duplicate']=1
        self.assertEqual(report.common_support(tap,curves,set()),'duplicate_endpoint')

    def test_extrapolation_and_nonmonotone_curves_rejected(self):
        for x,y,t in [([0,1],[0,1],[-.001]),([1,0],[0,1],[.5]),([0,0],[0,1],[0]),([0,1],[0,np.nan],[.5])]:
            with self.assertRaises(ValueError):report.interpolate(x,y,t)
        np.testing.assert_allclose(report.interpolate([0,1],[2,4],[0,.25,1]),[2,2.5,4])

    def test_missing_scores_cannot_look_like_perfect_agreement(self):
        empty=report.pressure_metric([dict(included=False,error_Cp=None)])
        self.assertEqual(empty['scored_taps'],0)
        self.assertIsNone(empty['rmse'])
        metric=report.pressure_metric([dict(included=True,error_Cp=.1),dict(included=True,error_Cp=-.1)])
        self.assertAlmostEqual(metric['rmse'],.1);self.assertAlmostEqual(metric['bias'],0.)

    def test_recovery_gradient_respects_unequal_tap_spacing(self):
        intervals=report.recovery_intervals([.8,.85,.86],[-.8,-.7,-.65])
        self.assertGreater(intervals[0]['delta_Cp'],intervals[1]['delta_Cp'])
        self.assertGreater(intervals[1]['slope'],intervals[0]['slope'])
        with self.assertRaises(ValueError):report.recovery_intervals([.8,.8],[1,2])

    def test_published_recovery_interval_comes_from_source_markers(self):
        taps=report.associate_taps(rows(report.QUAL/'pressure/experiment.csv'),rows(report.QUAL/'pressure/tap_coordinate_checks.csv'))
        for case in report.CASES:
            selected=[r for r in taps if r['case']==case and r['surface']=='aft_upper' and r['table_x'] is not None and .8<=r['table_x']<=.98]
            intervals=report.recovery_intervals([r['table_x'] for r in selected],[float(r['Cp']) for r in selected])
            peak=max(intervals,key=lambda r:r['slope'])
            self.assertEqual((peak['left_x'],peak['right_x']),(.88215,.89313))

    def test_residual_auditor_never_launches_nonlinear_solver(self):
        tree=ast.parse(inspect.getsource(audit.recompute))
        calls=[n.func.attr for n in ast.walk(tree) if isinstance(n,ast.Call) and isinstance(n.func,ast.Attribute)]
        self.assertTrue({'_system','_build_wakes','_build_source_coupling'}<=set(calls))
        self.assertFalse({'solve','_newton','_initialize'}&set(calls))

    def test_report_scores_identical_supported_taps_for_every_method(self):
        pressure=rows(report.OUT/'pressure_at_taps.csv')
        masks=[]
        for method in ('mses',)+report.SCENARIOS:
            for coordinate in ('table','plotted'):
                mask={(r['case'],r['surface'],r['table_index']) for r in pressure
                    if r['method']==method and r['coordinate']==coordinate and r['included']=='True'}
                self.assertEqual(len(mask),174);masks.append(mask)
        self.assertTrue(all(mask==masks[0] for mask in masks))
        summary=json.loads((report.OUT/'summary.json').read_text())
        self.assertTrue(summary['full_residual_recomputation_pass'])
        self.assertEqual((summary['converged_targets'],summary['failed_requests']),(9,0))


if __name__=='__main__':unittest.main()
