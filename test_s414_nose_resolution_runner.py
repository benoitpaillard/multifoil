"""Short contract tests with synthetic state arrays and fake solver; no flow solves."""
from pathlib import Path
import contextlib
import io
import sys
import tempfile
import time
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT/'validation'))
import run_s414_nose_resolution as run


def archived():
    return dict(feasible=True, selected_flow_index=4, selected_path='cold_local', alpha_deg=4.,
        initialization_paths=[dict(name='cold_local', start_alpha=-1., seed_alphas=[])],
        records=[dict(converged=True, alpha_deg=float(i), initialization=dict(path='cold_local', index=i)) for i in range(5)])


class Histories(unittest.TestCase):
    def test_archived_history_ends_at_fixed_selected_incidence(self):
        start, requests = run.archived_requests(archived())
        self.assertEqual(start, -1.)
        self.assertEqual([r['alpha_deg'] for r in requests], list(map(float, range(5))))
        self.assertEqual([r['role'] for r in requests], ['initialization']*4+['selected'])

    def test_failed_or_different_history_rejected(self):
        for key, value in [('feasible', False), ('selected_path', 'warm_anchor'), ('selected_flow_index', 3), ('alpha_deg', 2.)]:
            data=archived();data[key]=value
            with self.assertRaises(ValueError):run.archived_requests(data)
        data=archived();data['records'][0]['converged']=False
        with self.assertRaises(ValueError):run.archived_requests(data)


class BaselineGate(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.root=Path(self.tmp.name)
        self.reference=self.root/'reference';self.actual=self.root/'actual'
        self.reference.mkdir();self.actual.mkdir()
        self.result=dict(feasible=True,alpha_deg=-3.,CL=.056,CD=.01,CM=-.1)
        for folder in (self.reference,self.actual):
            run.write(folder/'result.json',self.result)
            for element in ('fore','aft'):
                np.savez(folder/(element+'_state.npz'),**{f:np.arange(3.) for f in run.STATE_FIELDS})
        self.job=dict(reference_result='reference/result.json')

    def tearDown(self):self.tmp.cleanup()

    def test_gate_demands_all18_arrays_and_exact_forces(self):
        gate=run.compare_control(self.actual,self.job,self.root)
        self.assertTrue(gate['passed']);self.assertEqual(len(gate['array_comparisons']),18)
        different=dict(self.result,CD=self.result['CD']+1e-15)
        run.write(self.actual/'result.json',different)
        self.assertFalse(run.compare_control(self.actual,self.job,self.root)['passed'])

    def test_changed_array_rejects_control(self):
        arrays={f:np.arange(3.) for f in run.STATE_FIELDS};arrays['U'][1]+=1e-15
        np.savez(self.actual/'fore_state.npz',**arrays)
        self.assertFalse(run.compare_control(self.actual,self.job,self.root)['passed'])


class Restart(unittest.TestCase):
    def test_completed_artifact_change_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder=Path(tmp);(folder/'data').write_text('original')
            run.write(folder/'result.json',dict(complete=True,artifact_sha256={'data':run.sha(folder/'data')}))
            run.verify_result(folder,{})
            (folder/'data').write_text('changed')
            with self.assertRaisesRegex(ValueError,'artifact changed'):run.verify_result(folder,{})

    def test_interrupted_attempt_is_preserved_without_relaunch(self):
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp);(out/'case').mkdir();(out/'case'/'run.log').write_text('partial attempt')
            with patch.object(run.subprocess,'Popen') as launch:
                with self.assertRaisesRegex(ValueError,'Interrupted attempt'):
                    run.run_phase(out,{},[dict(name='case')],2,float('inf'),lambda *a,**k:None)
                launch.assert_not_called()
            self.assertEqual((out/'case'/'run.log').read_text(),'partial attempt')


class FailedReplay(unittest.TestCase):
    def test_failure_stops_history_and_carries_no_selected_forces(self):
        import s414_nose_resolution as adapter
        import multifoil
        solution=SimpleNamespace(alpha=-3.,cl=float('nan'),cd=float('nan'),cm=float('nan'),converged=False,residual=1.,
            diagnostics=dict(iterations=100,continuation=[dict(alpha_deg=-3.,converged=False,CL=float('nan'))],
                             wall_seconds=.1,failure_reason='synthetic failure'))
        fake=SimpleNamespace(resolution_controls={},models=[],initialized=False,last_alpha=None,elements=[])
        calls=[]
        def solve(alpha):calls.append(alpha);return solution
        fake.solve=solve
        spec=dict(qualification='fake',reynolds=1e6,mach=.1,ncrit=9.,iterations=100,residual_tolerance=1e-7,
            wake_length=1.,minimum_te_gap=1e-5,continuation_step=1.,min_step=.25)
        job=dict(name='fake',subdivision_factor=1,start_alpha=-3.,
            requests=[dict(alpha_deg=-3.,role='initialization'),dict(alpha_deg=-2.,role='selected')])
        with tempfile.TemporaryDirectory() as tmp, patch.object(adapter,'Solver',return_value=fake), \
                patch.object(adapter,'build_elements',return_value=[]),patch.object(multifoil,'save_results'):
            result=run.worker(Path(tmp),job,spec)
        self.assertEqual(calls,[-3.]);self.assertFalse(result['feasible']);self.assertTrue(result['complete'])
        self.assertEqual(result['execution_state'],'complete')
        self.assertTrue(all(result[k] is None for k in ('CL','CD','CM','alpha_deg')))
        self.assertEqual(len(result['records']),1)


class RemoteIdentity(unittest.TestCase):
    def test_remote_programs_compile_without_executing(self):
        import s414_nose_resolution_remote as remote
        for name in ('SUBMIT','POLL','COLLECT'):
            compile(getattr(remote,name),name,'exec')

    def test_poll_recovers_lost_submit_receipt_only_for_prepared_identity(self):
        import s414_nose_resolution_remote as remote
        with tempfile.TemporaryDirectory() as tmp:
            launch=Path(tmp);preparation=dict(archive_sha256='archive',study_sha256='study')
            run.write(launch/'preparation.json',preparation)
            receipt=dict(job_id='pilot',**preparation)
            response=dict(job_id='pilot',status='running',launch=receipt)
            with patch.object(remote,'LAUNCH',launch),patch.object(remote,'remote',return_value=response), \
                    patch.object(sys,'argv',['remote','poll','--job-id','pilot']),contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(remote.main(),0)
            self.assertEqual(run.read(launch/'launch_receipt.json'),receipt)

    def test_poll_cannot_adopt_a_different_job(self):
        import s414_nose_resolution_remote as remote
        with tempfile.TemporaryDirectory() as tmp:
            launch=Path(tmp);preparation=dict(archive_sha256='archive',study_sha256='study')
            run.write(launch/'preparation.json',preparation)
            response=dict(launch=dict(job_id='different',**preparation))
            with patch.object(remote,'LAUNCH',launch),patch.object(remote,'remote',return_value=response), \
                    patch.object(sys,'argv',['remote','poll','--job-id','pilot']):
                with self.assertRaisesRegex(ValueError,'prepared identity'):remote.main()
            self.assertFalse((launch/'launch_receipt.json').exists())


class CompleteLaunchBudget(unittest.TestCase):
    def test_timeout_kills_supervisor_and_detached_recorded_worker(self):
        import run_s414_nose_resolution_launch as guard
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);(root/'validation').mkdir()
            worker=root/'validation/run_s414_nose_resolution.py'
            worker.write_text('import signal,time\nsignal.signal(signal.SIGTERM,signal.SIG_IGN)\ntime.sleep(60)\n')
            attempt=root/run.STUDY/'fake';attempt.mkdir(parents=True)
            driver=root/'fake_supervisor.py'
            driver.write_text('import json,pathlib,signal,subprocess,sys,time\n'
                'signal.signal(signal.SIGTERM,signal.SIG_IGN)\n'
                "p=subprocess.Popen([sys.executable,'validation/run_s414_nose_resolution.py','--worker','fake'],start_new_session=True)\n"
                f"pathlib.Path({str(attempt/'attempt.json')!r}).write_text(json.dumps(dict(job='fake',pid=p.pid)))\n"
                'time.sleep(60)\n')
            with (root/'log').open('wb') as log:
                with self.assertRaises(guard.subprocess.TimeoutExpired):
                    guard.execute([sys.executable,str(driver)],root,log,.4,grace_seconds=.1)
            until=time.monotonic()+1
            while guard.worker_groups(root) and time.monotonic()<until:time.sleep(.01)
            self.assertEqual(guard.worker_groups(root),[])

    def test_timeout_also_kills_unregistered_detached_preflight_child(self):
        import run_s414_nose_resolution_launch as guard
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);driver=root/'fake_preflight.py'
            driver.write_text('import pathlib,signal,subprocess,sys,time\n'
                'signal.signal(signal.SIGTERM,signal.SIG_IGN)\n'
                "p=subprocess.Popen([sys.executable,'-c','import time;time.sleep(60)'],start_new_session=True)\n"
                "pathlib.Path('child.pid').write_text(str(p.pid))\ntime.sleep(60)\n")
            with (root/'log').open('wb') as log:
                with self.assertRaises(guard.subprocess.TimeoutExpired):
                    guard.execute([sys.executable,str(driver)],root,log,.4,grace_seconds=.1)
            pid=int((root/'child.pid').read_text());proc=Path('/proc')/str(pid)/'stat';until=time.monotonic()+1
            while proc.exists() and proc.read_text().split()[2]!='Z' and time.monotonic()<until:time.sleep(.01)
            self.assertTrue(not proc.exists() or proc.read_text().split()[2]=='Z')

    def test_no_phase_launch_after_budget_expires(self):
        import run_s414_nose_resolution_launch as guard
        with patch.object(guard.subprocess,'Popen') as process:
            with self.assertRaises(TimeoutError):guard.execute(['unused'],Path('/tmp'),None,0.)
            process.assert_not_called()


if __name__=='__main__':unittest.main()
