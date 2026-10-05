"""Regressions for the declared three-grid decision, coverage and reporting."""
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).parent/'validation'))
from report_s702_panel_convergence import assess_contraction, step_metrics


class PanelAuditTests(unittest.TestCase):
    def assess(self, individual, paired, **overrides):
        options = dict(coverage=True, control=True, audited=True)
        options.update(overrides)
        return assess_contraction(individual, paired, **options)

    def test_smaller_paired_maximum_cannot_hide_growing_individual_response(self):
        result = self.assess([[.02,.06], [.37,.27], [-.49,-.70], [-.68,.07]],
                             [[-.50,-.75], [-1.05,-.21]])
        self.assertFalse(result['numerical_screen_pass'])
        self.assertTrue(result['paired']['maximum_contracts'])
        self.assertFalse(result['individual']['maximum_contracts'])

    def test_coverage_control_and_audit_are_each_required(self):
        for gate in ('coverage', 'control', 'audited'):
            result = self.assess([], [], **{gate:False})
            self.assertFalse(result['decision_eligible'])
            self.assertFalse(result['numerical_screen_pass'])
            self.assertNotIn('paired', result)
        for bad in ([[1.,.1]]*3, [[1.,float('nan')]]*4, [[1.,None]]*4):
            with self.assertRaisesRegex(ValueError, 'Complete finite'):
                self.assess(bad, [[1.,.1]]*2)

    def test_contracting_changes_must_also_fall_strictly_below_threshold(self):
        self.assertTrue(self.assess([[1.,.1]]*4, [[1.,.2]]*2)['numerical_screen_pass'])
        self.assertFalse(self.assess([[1.,.5]]*4, [[1.,.2]]*2)['numerical_screen_pass'])

    def test_nonmonotonicity_is_reported_even_when_maximum_screen_passes(self):
        row = step_metrics(-.68, .07)
        self.assertTrue(row['reverses_direction'])
        self.assertTrue(row['contracts'])
        self.assertTrue(self.assess([[-.68,.07]]*4, [[1.,.2]]*2)['numerical_screen_pass'])
        self.assertFalse(step_metrics(.02,.06)['contracts'])
        self.assertIsNone(step_metrics(0.,.06)['absolute_step_ratio'])


if __name__ == '__main__':
    unittest.main()
