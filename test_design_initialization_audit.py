"""Audit rejection checks against the preserved completed integration study."""
import copy
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).parent/'validation'))
import report_design_initialization as report


class IntegrationAuditTests(unittest.TestCase):
    def setUp(self):
        self.study = report.ROOT/'validation/results/design_initialization_20260921'
        self.inputs = report.Inputs()
        self.spec = self.inputs.json(self.study/'study.json')
        self.runtime = self.inputs.json(self.study/'runtime.json')

    def check_job(self, name, mutation=None):
        job = next(j for j in self.spec['jobs'] if j['name'] == name)
        folder = self.study/name
        inputs = report.Inputs()
        original = inputs.json
        def read(path):
            value = original(path)
            if mutation is not None and Path(path) == folder/'result.json':
                value = copy.deepcopy(value); mutation(value)
            return value
        inputs.json = read
        return report.audit(folder, job, self.spec, self.runtime, inputs)

    def test_successful_selected_state_and_failed_target_are_distinct(self):
        self.assertTrue(self.check_job('d5_plus_policy')['result']['feasible'])
        self.assertFalse(self.check_job('d10_baseline_policy')['result']['feasible'])

    def test_wrong_selected_index_is_rejected(self):
        with self.assertRaisesRegex(ValueError, 'Selected flow index'):
            self.check_job('d5_plus_policy', lambda r:r.update(selected_flow_index=0))

    def test_failed_target_cannot_acquire_finite_final_drag(self):
        with self.assertRaisesRegex(ValueError, 'Failed flow has final forces'):
            self.check_job('d10_baseline_policy', lambda r:r.update(CD=.01))

    def test_deleted_failure_cannot_hide_missing_work(self):
        with self.assertRaisesRegex(ValueError, 'Work count mismatch'):
            self.check_job('d10_baseline_policy', lambda r:r['records'].pop())


if __name__ == '__main__':
    unittest.main()
