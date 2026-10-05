"""Operational checks for the exact baseline classification and live monitor."""
import contextlib
import io
import unittest
from unittest.mock import patch

from handover import resume_s702_policy as resume
from handover import policy_status as status


class PolicyOperationsTests(unittest.TestCase):
    def good_log(self, count=46):
        return f'FAIL: {resume.FAILURE}\n{resume.ASSERTION}\nRan {count} tests in 1.0s\n\nFAILED (failures=1)\n'

    def test_only_recorded_failure_is_accepted(self):
        resume.check_failure(self.good_log(), 46)
        resume.check_failure(self.good_log(1), 1)
        for text in (self.good_log()+'ERROR: another_test\n',
                     self.good_log()+'FAIL: another_test\n',
                     self.good_log().replace(resume.ASSERTION, 'AssertionError: another failure'),
                     self.good_log().replace('Ran 46', 'Ran 45')):
            with self.subTest(text=text), self.assertRaises(ValueError):
                resume.check_failure(text, 46)

    def test_finished_failure_exits_watch_with_failure_code(self):
        data = dict(launch={'phase': 'finished_with_regression_failures', 'exit_code': 2},
                    launcher_alive=False)
        with patch.object(status, 'fetch', return_value=data), patch.object(status, 'render'), \
                patch('sys.argv', ['policy_status.py', '--watch']), patch.object(status.time, 'sleep') as sleep:
            self.assertEqual(status.main(), 2)
            sleep.assert_not_called()

    def test_dead_launcher_cannot_look_like_running_job(self):
        data = dict(launch={'phase': 'running_study'}, launcher_alive=False)
        with patch.object(status, 'fetch', return_value=data), patch.object(status, 'render'), \
                patch('sys.argv', ['policy_status.py', '--watch']):
            self.assertEqual(status.main(), 1)

    def test_completed_is_separate_from_converged(self):
        data = dict(launch={'phase': 'finished_with_regression_failures', 'exit_code': 2},
                    observed_utc='test', known_baseline_failure=True, launcher_alive=False,
                    study={'completed_jobs': 7, 'total_jobs': 7}, jobs=[],
                    summary={'converged_primary': 5, 'converged_control': 1, 'force_regression_pass': False})
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            status.render(data)
        self.assertIn('Completed: 7/7', output.getvalue())
        self.assertIn('Primary converged: 5/6', output.getvalue())
        self.assertIn('Force regression: FAIL', output.getvalue())


if __name__ == '__main__':
    unittest.main()
