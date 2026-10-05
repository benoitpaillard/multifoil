"""Frozen lower-bucket protocol and archived-control gates; no flow solves."""
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
import numpy as np

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'validation'))
import run_s414_lower_bucket as run
import run_s414_reconstruction as previous


class LowerBucketTests(unittest.TestCase):
    def test_four_targets_preserve_baseline_settings_and_worker(self):
        spec=run.specification();old=run.read(ROOT/run.OLD/'study.json')
        self.assertEqual(len(spec['jobs']),4)
        self.assertEqual(sorted(j['target_cl'] for j in spec['jobs']),[-.072,.056,.175,.284])
        for key in ('reynolds','mach','ncrit','panels','reference_chord','iterations','cl_tolerance',
                    'residual_tolerance','initial_half_width','alpha_bounds','max_evaluations','max_brackets',
                    'warm_start_alpha','entry_offset','expected_versions'):
            self.assertEqual(spec[key],old[key])
        self.assertTrue(all(j['scenario']=='nominal' and j['wake_length']==1. for j in spec['jobs']))
        self.assertEqual(sum(j['control'] for j in spec['jobs']),1)
        self.assertIs(run.worker,previous.worker)

    def test_control_initialization_is_exact_and_missing_mses_is_not_extrapolated(self):
        spec=run.specification();old=run.read(ROOT/run.OLD/'study.json')
        control=next(j for j in spec['jobs'] if j['control'])
        archived=next(j for j in old['jobs'] if j['name']=='cl028_nominal')
        for key in ('center','target_cl','experiment_alpha','scenario','wake_length'):
            self.assertEqual(control[key],archived[key])
        negative=next(j for j in spec['jobs'] if j['target_cl']<0)
        self.assertEqual(negative['center'],-4.09)
        self.assertFalse(negative['mses_drag_available']);self.assertFalse(negative['mses_pressure_available'])
        self.assertEqual(sum(j['mses_pressure_available'] for j in spec['jobs']),1)
        self.assertEqual(sum(j['mses_drag_available'] for j in spec['jobs']),3)

    def test_qualified_reference_tampering_is_rejected(self):
        original=run.sha
        with patch.object(run,'sha',side_effect=lambda p:'bad' if str(p).endswith(run.SOURCE+'/free_boundary_pressure.csv') else original(p)):
            with self.assertRaisesRegex(ValueError,'Qualified measurement changed'):run.specification()

    def test_all_state_arrays_and_selected_path_are_required(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);a=root/'actual';b=root/'reference';a.mkdir();b.mkdir()
            r=dict(feasible=True,CD=.007,selected_path='cold_local')
            for p in (a,b):
                run.write(p/'result.json',r)
                for element in ('fore','aft'):np.savez(p/(element+'_state.npz'),U=[1.,2.],turbulent=[False,True])
            self.assertTrue(run.compare_control(a,b,.01)['passed'])
            np.savez(a/'aft_state.npz',U=[1.,2.+1e-12],turbulent=[False,True])
            self.assertFalse(run.compare_control(a,b,.01)['passed'])
            np.savez(a/'aft_state.npz',U=[1.,2.],turbulent=[False,True])
            run.write(a/'result.json',dict(r,selected_path='warm_anchor'))
            self.assertFalse(run.compare_control(a,b,.01)['passed'])

    def test_failed_control_cannot_open_new_case_gate(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);(root/run.CONTROL).mkdir();ref=root/'reference';ref.mkdir()
            run.write(root/run.CONTROL/'result.json',dict(feasible=False,CL=None,CD=None))
            run.write(ref/'result.json',dict(feasible=True,CL=.284,CD=.007))
            with patch.object(run,'ROOT',root):
                with self.assertRaisesRegex(ValueError,'not launched'):
                    run.check_control(root,dict(reference_control='reference',acceptance_control_CD_counts=.01))
            self.assertFalse(run.read(root/'controls.json')['passed'])

    def test_failed_target_cannot_export_finite_force(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);(root/'bad').mkdir()
            run.write(root/'bad/result.json',dict(feasible=False,CL=None,CD=.01,CM=None,alpha_deg=None))
            with self.assertRaisesRegex(ValueError,'finite forces'):
                run.summarize(root,dict(jobs=[dict(name='bad')]))

    def test_manifest_refuses_setting_changes_on_resume(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);spec=dict(ncrit=9.,sha256={'viscous.py':run.sha(ROOT/'viscous.py')})
            run.prepare(root,spec);run.prepare(root,spec)
            with self.assertRaisesRegex(ValueError,'changed'):run.prepare(root,dict(spec,ncrit=8.))


if __name__=='__main__':unittest.main()
