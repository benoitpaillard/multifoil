"""Guard audit/report distinctions between fixed angle and fixed lift."""
import copy
import csv
from pathlib import Path
import sys
import unittest
sys.path.insert(0, str(Path(__file__).resolve().parent/'validation'))
from audit_s414_incidence import check_fixed_angle_identity, OUT
from report_s414_incidence import condition_metrics


class IncidenceAnalysisTests(unittest.TestCase):
    def setUp(self):
        self.job = dict(comparison_basis='fixed_alpha', target_alpha=-3.07, reference_experimental_cl=.056)
        self.ref = dict(target_cl='.056', experiment_CD='.00618', experiment_caption_alpha='-3.07',
            appendix_CL='.0556', appendix_alpha='-3.074', appendix_CD='.006177')
        self.result = dict(feasible=True, comparison_basis='fixed_alpha', target_alpha=-3.07,
            reference_experimental_cl=.056, CL=.1072, CD=.01142, alpha_deg=-3.07,
            records=[dict(initialization=dict(index=0, path='cold_measured_incidence', phase='fixed_alpha'))],
            selected_flow_index=0, selected_path='cold_measured_incidence',
            attempts=[dict(path='cold_measured_incidence', start_alpha=-3.07, requested_alpha=-3.07, feasible=True)])

    def test_fixed_angle_lift_disagreement_is_reported_without_matched_lift_score(self):
        check_fixed_angle_identity(self.result, self.job, {})
        metric = condition_metrics(self.result, self.job, self.ref)
        self.assertAlmostEqual(metric['CL_minus_caption'], .0512)
        self.assertAlmostEqual(metric['CL_minus_appendix'], .0516)
        self.assertIsNone(metric['matched_lift_CD_error_counts'])

    def test_fixed_lift_score_keeps_count_units(self):
        job = dict(comparison_basis='fixed_lift', target_cl=.056)
        metric = condition_metrics(dict(self.result, CD=.00628), job, self.ref)
        self.assertAlmostEqual(metric['matched_lift_CD_error_counts'], 1.)

    def test_hidden_lift_target_is_rejected(self):
        with self.assertRaisesRegex(ValueError, 'lift target'):
            check_fixed_angle_identity(self.result, dict(self.job, target_cl=.056), {})
        with self.assertRaisesRegex(ValueError, 'lift target'):
            condition_metrics(self.result, dict(self.job, target_cl=.056), self.ref)

    def test_wrong_angle_and_extra_evaluation_are_rejected(self):
        for field, value in [('alpha_deg', -3.06), ('records', self.result['records']*2)]:
            with self.assertRaises(ValueError): check_fixed_angle_identity(dict(self.result, **{field:value}), self.job, {})

    def test_seed_or_alternative_path_cannot_be_presented_as_measured_case(self):
        modified = copy.deepcopy(self.result); modified['records'][0]['initialization']['phase'] = 'warm_seed'
        with self.assertRaisesRegex(ValueError, 'phase'): check_fixed_angle_identity(modified, self.job, {})
        with self.assertRaisesRegex(ValueError, 'warm path'):
            check_fixed_angle_identity(dict(self.result, selected_path='warm_anchor'), self.job, {})

    def test_failed_state_and_reassociated_measurement_are_rejected(self):
        with self.assertRaisesRegex(ValueError, 'Failed state'):
            condition_metrics(dict(self.result, feasible=False), self.job, self.ref)
        with self.assertRaisesRegex(ValueError, 'angle changed'):
            condition_metrics(self.result, self.job, dict(self.ref, experiment_caption_alpha='-3.074'))

    def test_exported_turbulent_shear_is_not_amplification(self):
        with (OUT/'fore_amplification.csv').open() as stream: rows = list(csv.DictReader(stream))
        self.assertTrue(any(r['laminar']=='False' for r in rows))
        for r in rows:
            self.assertEqual(bool(r['N']), r['laminar']=='True')
            self.assertEqual(bool(r['shear']), r['laminar']=='False')


if __name__ == '__main__': unittest.main()
