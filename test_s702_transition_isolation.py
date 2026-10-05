"""Threshold assignment, geometry-source checks and launch/resume safeguards."""
from pathlib import Path
from types import SimpleNamespace
import sys
import tempfile
import unittest
from unittest.mock import patch
import numpy as np

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'validation'))
import run_s702_transition_isolation as run
from qualify_s702_geometry import register_outline


class TransitionIsolationTests(unittest.TestCase):
    def test_registration_preserves_x_inconsistency_and_relative_placement(self):
        baseline=np.array([[0.,0.],[1.,.1]])
        fore=np.array([[10.02,30.],[30.02,28.]])
        aft=fore+np.array([14.,1.])
        (a,b),_=register_outline(fore,aft,10.,20.,baseline)
        np.testing.assert_allclose(a[:,0]-baseline[:,0],.001,atol=1e-15)
        np.testing.assert_allclose(a[:,1],baseline[:,1],atol=1e-15)
        np.testing.assert_allclose(b-a,[[.7,-.05],[.7,-.05]],atol=1e-15)

    def test_published_outlines_do_not_resolve_tap_support(self):
        s=run.read(ROOT/run.GEO/'summary.json')
        self.assertTrue(s['original_figure1_reextraction_matches'])
        self.assertEqual(s['pressure_outline_comparisons'],8)
        self.assertLess(s['maximum_chordwise_difference'],2e-5)
        self.assertLess(s['maximum_ordinate_difference_after_fore_registration'],5e-6)
        self.assertGreater(min(s['aft_nose_minus_tap_range']),.0008)
        self.assertFalse(s['correction_justified'])

    def test_complete_factorial_and_original_control_centers(self):
        spec=run.specification();old=run.read(ROOT/run.PARENT/'study.json')
        mixed=[j for j in spec['jobs'] if not j['control']]
        self.assertEqual(len(spec['jobs']),10);self.assertEqual(len(mixed),8)
        for case in {j['pressure_case'] for j in mixed}:
            jobs=[j for j in mixed if j['pressure_case']==case]
            self.assertEqual({tuple(j['ncrit_by_element'][e] for e in ('fore','aft')) for j in jobs},{(6.25,9.),(9.,6.25)})
            endpoints=[j for j in old['jobs'] if j['pressure_case']==case]
            self.assertEqual({j['center'] for j in jobs},{sum(j['center'] for j in endpoints)/2})
        for j in spec['jobs'][:2]:
            self.assertTrue(j['control'])
            original=next(v for v in old['jobs'] if v['name']==Path(j['reference_result']).parent.name)
            self.assertEqual(j['center'],original['center'])
        self.assertEqual(spec['panels'],[340,240]);self.assertEqual(spec['reynolds'],1e6)

    def test_every_fresh_solver_receives_named_thresholds(self):
        created=[]
        def constructor(*args,**kwargs):
            value=SimpleNamespace(models=[SimpleNamespace(geom=SimpleNamespace(name=e),param=SimpleNamespace(ncrit=kwargs['ncrit'])) for e in ('aft','fore')])
            created.append(value);return value
        spec=dict(reynolds=1e6,mach=.1,iterations=100,residual_tolerance=1e-7)
        j=dict(wake_length=1.,ncrit_by_element=dict(fore=9.,aft=6.25))
        with patch.object(run,'CoupledSolver',side_effect=constructor):
            first=run.make_solver([],j,spec,0.);first.models[0].param.ncrit=3.
            second=run.make_solver([],j,spec,-1.)
        self.assertEqual({m.geom.name:m.param.ncrit for m in second.models},j['ncrit_by_element'])
        self.assertEqual(len(created),2)

    def test_invalid_thresholds_rejected_before_constructing(self):
        for values in ({'fore':9.},{'fore':9.,'aft':float('nan')},{'fore':9.,'aft':0.}):
            with patch.object(run,'CoupledSolver') as constructor,self.assertRaises(ValueError):
                run.make_solver([],dict(ncrit_by_element=values),{},0.)
            constructor.assert_not_called()

    def test_common_threshold_keeps_original_constructor(self):
        for n in (9.,6.25):
            model=SimpleNamespace(models=[SimpleNamespace(geom=SimpleNamespace(name=e),param=SimpleNamespace(ncrit=n)) for e in ('fore','aft')])
            with patch.object(run,'CoupledSolver',return_value=model) as constructor:
                run.make_solver([],dict(wake_length=1.,ncrit_by_element=dict(fore=n,aft=n)),dict(reynolds=1e6,mach=.1,iterations=100,residual_tolerance=1e-7),1.2)
            self.assertEqual(constructor.call_args.kwargs['ncrit'],n)
            self.assertEqual([m.param.ncrit for m in model.models],[n,n])

    def test_control_rejects_changed_state_even_with_identical_drag(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);out=root/'out';out.mkdir();jobs=[]
            for i in range(2):
                reference=root/f'ref{i}';reference.mkdir();folder=out/f'control{i}';folder.mkdir()
                r=dict(feasible=True,CD=.01,selected_path='cold_local')
                run.write(reference/'result.json',r);run.write(folder/'result.json',r)
                for e in ('fore','aft'):
                    np.savez(reference/(e+'_state.npz'),U=np.array([1.]))
                    np.savez(folder/(e+'_state.npz'),U=np.array([1.]))
                jobs.append(dict(name=folder.name,control=True,reference_result=str(reference/'result.json')))
            spec=dict(jobs=jobs,acceptance_control_CD_counts=.01)
            self.assertTrue(run.check_controls(out,spec)['passed'])
            np.savez(out/'control0/aft_state.npz',U=np.array([1.00001]))
            self.assertFalse(run.check_controls(out,spec)['passed'])

    def test_changed_settings_cannot_resume(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder=Path(tmp);spec=dict(sha256={},jobs=[dict(ncrit_by_element=dict(fore=9.,aft=6.25))])
            run.prepare(folder,spec)
            with self.assertRaises(ValueError):run.prepare(folder,dict(spec,jobs=[dict(ncrit_by_element=dict(fore=6.25,aft=9.))]))


if __name__=='__main__':unittest.main()
