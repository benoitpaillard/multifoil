"""Coordinator safety/regression checks. Fake local commands only; no CFD/SSH/Codex."""
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

spec = importlib.util.spec_from_file_location('pilot', Path(__file__).with_name('pilot.py'))
pilot = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pilot)


FAKE = r'''
import hashlib,json,pathlib,sys,time
operation,result,job_id,run=sys.argv[1:]
result=pathlib.Path(result)
run=pathlib.Path(run)
artifact=run/(operation+'.txt')
artifact.write_text('independent evidence\n')
evidence=[{'path':str(artifact),'sha256':hashlib.sha256(artifact.read_bytes()).hexdigest()}]
statuses={'prepare':'prepared','submit':'running','poll':'complete','collect':'collected'}
receipt={'status':statuses.get(operation,'audited'),'job_id':job_id,'evidence':evidence}
if operation=='submit':
    count=run/'submission_count.txt'
    count.write_text(str(int(count.read_text())+1) if count.exists() else '1')
    if (run/'simulate_lost_reply').exists():
        sys.exit(7)
if operation=='audit':
    receipt.update(evidence_verified=True,criteria=[{'id':'resolution','outcome':'fail','evidence':[str(artifact)]}])
if operation=='slow':
    time.sleep(5)
result.write_text(json.dumps(receipt))
'''


class PilotTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        (self.root / 'fake.py').write_text(FAKE)
        self.plan = {'schema': 1, 'id': 'fake-pilot', 'project_root': str(self.root),
                     'resources': {'remote_workers': 2, 'max_final_cases': 6},
                     'files': ['fake.py'], 'criteria': [{'id': 'resolution'}], 'commands': {}}
        for step in ('prepare', 'submit', 'poll', 'collect', 'audit'):
            self.plan['commands'][step] = {
                'argv': [sys.executable, '{snapshot}/fake.py', step, '{result}', '{job_id}', '{run}'],
                'timeout_seconds': 10}
        self.path = self.root / 'plan.yaml'
        self.path.write_text(json.dumps(self.plan))
        self.run = self.root / 'run'
        pilot.freeze(self.path, self.run)

    def tearDown(self):
        self.temporary.cleanup()

    def accepted(self, run=None):
        run = run or self.run
        pilot.operate(run, 'prepare')
        # Stub only the independent model process; exercise gate interpretation.
        result = {'status': 'PASS', 'summary': 'protocol verified', 'findings': [], 'unresolved': [],
                  'criteria': [{'id': 'resolution', 'outcome': 'not_evaluated', 'evidence': []}]}
        with patch.object(pilot, 'command', return_value=result):
            pilot.operate(run, 'prelaunch_review')

    def test_complete_scientific_failure_is_audited_and_reviewable(self):
        self.accepted()
        self.assertEqual(pilot.operate(self.run, 'submit')['status'], 'REMOTE_ACTIVE')
        self.assertEqual(pilot.operate(self.run, 'poll')['status'], 'REMOTE_FINISHED')
        pilot.operate(self.run, 'collect')
        state = pilot.operate(self.run, 'audit')
        self.assertEqual(state['status'], 'AUDITED')
        self.assertEqual(state['evidence_verification'], 'passed')
        self.assertEqual(state['scientific_criteria'][0]['outcome'], 'fail')
        self.assertEqual(json.loads((self.root/'orchestration/.remote_resources.json').read_text())['claims'], {})

    def test_crash_after_remote_submission_cannot_duplicate_job(self):
        self.accepted()
        (self.run / 'simulate_lost_reply').touch()
        with self.assertRaises(pilot.PilotError):
            pilot.operate(self.run, 'submit')
        self.assertEqual(pilot.operate(self.run, 'status')['status'], 'SUBMISSION_PENDING')
        with self.assertRaises(pilot.PilotError):
            pilot.operate(self.run, 'submit')
        self.assertEqual((self.run/'submission_count.txt').read_text(), '1')
        self.assertEqual(pilot.operate(self.run, 'advance')['status'], 'REMOTE_FINISHED')

    def test_mutation_of_frozen_source_blocks_next_step(self):
        (self.run/'snapshot/fake.py').chmod(0o644)
        (self.run/'snapshot/fake.py').write_text('changed')
        with self.assertRaisesRegex(pilot.PilotError, 'Frozen input changed'):
            pilot.operate(self.run, 'prepare')

    def test_failed_review_blocks_launch(self):
        pilot.operate(self.run, 'prepare')
        result = {'status': 'FAIL', 'summary': 'bad protocol', 'findings': ['missing control'], 'unresolved': [],
                  'criteria': [{'id': 'resolution', 'outcome': 'not_evaluated', 'evidence': []}]}
        with patch.object(pilot, 'command', return_value=result):
            state = pilot.operate(self.run, 'prelaunch_review')
        self.assertEqual(state['status'], 'PRELAUNCH_REJECTED')
        with self.assertRaises(pilot.PilotError):
            pilot.operate(self.run, 'submit')

    def test_worker_budget_shared_across_run_directories(self):
        self.accepted()
        pilot.operate(self.run, 'submit')
        other = self.root/'separate-results/run2'
        pilot.freeze(self.path, other)
        self.accepted(other)
        with self.assertRaisesRegex(pilot.PilotError, 'two-worker allowance'):
            pilot.operate(other, 'submit')
        self.assertFalse((other/'submission_count.txt').exists())

    def test_altered_completed_receipt_is_rejected(self):
        pilot.operate(self.run, 'prepare')
        (self.run/'attempts/prepare/001/result.json').write_text('{}')
        with self.assertRaisesRegex(pilot.PilotError, 'Recorded result changed'):
            pilot.operate(self.run, 'status')

    def test_external_evidence_hash_checked(self):
        (self.run/'evidence').write_text('content')
        with self.assertRaisesRegex(pilot.PilotError, 'hash mismatch'):
            pilot.check_evidence(self.run, {'evidence': [{'path': 'evidence', 'sha256': 'wrong'}]})

    def test_changed_prepared_payload_blocks_submission(self):
        self.accepted()
        (self.run/'prepare.txt').write_text('different payload')
        with self.assertRaisesRegex(pilot.PilotError, 'hash mismatch'):
            pilot.operate(self.run, 'submit')
        self.assertFalse((self.run/'submission_count.txt').exists())

    def test_unexpected_new_file_in_snapshot_is_rejected(self):
        (self.run/'snapshot/unreviewed.py').write_text('pass')
        with self.assertRaisesRegex(pilot.PilotError, 'Unexpected path'):
            pilot.operate(self.run, 'prepare')

    def test_declared_generated_directory_allowed_but_input_still_protected(self):
        self.plan['writable_snapshot_paths'] = ['generated/results']
        self.path.write_text(json.dumps(self.plan))
        run = self.root/'generated-run'
        pilot.freeze(self.path, run)
        output = run/'snapshot/generated/results'
        output.mkdir(parents=True)
        (output/'artifact').write_text('permitted')
        self.assertEqual(pilot.operate(run, 'status')['status'], 'FROZEN')

    def test_watch_cannot_launch_unsubmitted_job(self):
        self.assertEqual(pilot.watch(self.run, 10, 1)['status'], 'FROZEN')
        self.assertFalse((self.run/'submission_count.txt').exists())

    def test_status_remains_readable_during_mutating_operation(self):
        with pilot.locked(self.run/'.lock'):
            self.assertEqual(pilot.operate(self.run, 'status')['status'], 'FROZEN')

    def test_bad_snapshot_destination_rejected(self):
        self.plan['files'] = [{'source': 'fake.py', 'destination': '../escape.py'}]
        self.path.write_text(json.dumps(self.plan))
        with self.assertRaises(pilot.PilotError):
            pilot.freeze(self.path, self.root/'bad-run')

    def test_timeout_is_durable_and_does_not_mark_completion(self):
        self.plan['commands']['prepare']['argv'][2] = 'slow'
        self.plan['commands']['prepare']['timeout_seconds'] = .05
        self.path.write_text(json.dumps(self.plan))
        run = self.root/'timeout-run'
        pilot.freeze(self.path, run)
        with self.assertRaises(pilot.subprocess.TimeoutExpired):
            pilot.operate(run, 'prepare')
        state = pilot.operate(run, 'status')
        self.assertEqual(state['status'], 'FROZEN')
        self.assertEqual(state['steps']['prepare'][0]['status'], 'failed')

    def test_criterion_failure_not_mislabeled_as_execution_failure(self):
        pilot.check_criteria(self.plan, {'criteria': [
            {'id': 'resolution', 'outcome': 'fail', 'evidence': ['measurement.csv']}]})
        with self.assertRaises(pilot.PilotError):
            pilot.check_criteria(self.plan, {'criteria': []})


if __name__ == '__main__':
    unittest.main()
