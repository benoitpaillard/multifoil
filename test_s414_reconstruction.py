"""S414 extraction, reconstruction and execution guards; no flow solves."""
import csv
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT/'validation'))
import extract_s414_drag as drag
import run_s414_reconstruction as run
from s414_geometry import elements, SCENARIOS, SHIFT, cubic_repanel
from multifoil import read_dat, repanel


def rows(name):
    with (ROOT/run.DRAG/name).open() as stream:
        return [{k: float(v) if k in ('CL', 'CD', 'alpha_deg') else v for k,v in r.items()}
                for r in csv.DictReader(stream)]


class S414ReconstructionTests(unittest.TestCase):
    def test_drag_sign_and_caption_associations(self):
        source = run.read(ROOT/run.DRAG/'source.json')
        self.assertEqual((source['experiment_drag_markers'], source['mses_drag_vertices']), (22,24))
        expected = [.00655,.01449]
        for target, caption in zip(source['targets'], expected):
            value = drag.interpolate(rows('experiment_drag.csv'), target['CL'])
            self.assertAlmostEqual(value, caption, delta=6e-6)
            self.assertEqual(value, target['experiment_CD_interpolated'])
        np.testing.assert_allclose([t['mses_minus_experiment_counts'] for t in source['targets']],
                                   [5.331863388738315,-13.894908115366038], atol=1e-10)

    def test_no_extrapolation_or_nonmonotone_interpolation(self):
        curve = [dict(CL=0.,CD=.01),dict(CL=1.,CD=.02)]
        with self.assertRaises(ValueError): drag.interpolate(curve,1.1)
        with self.assertRaises(ValueError): drag.interpolate(curve[::-1],.5)
        with self.assertRaises(ValueError): drag.interpolate([dict(CL=0.,CD=float('nan')),curve[1]],.5)
        self.assertEqual(drag.interpolate(curve,.5),.015)

    def test_lift_seed_stops_at_stall_and_retains_original(self):
        curve = [dict(CL=0.,alpha_deg=0.),dict(CL=1.,alpha_deg=4.),dict(CL=.8,alpha_deg=8.)]
        self.assertEqual(drag.ascending_branch(curve),curve[:2])
        self.assertEqual(len(curve),3)
        with self.assertRaises(ValueError):
            drag.ascending_branch([dict(CL=.3,alpha_deg=0.),dict(CL=.2,alpha_deg=1.),dict(CL=.4,alpha_deg=2.)])

    def test_original_drag_lift_plot_offset_not_fitted_away(self):
        source = run.read(ROOT/run.DRAG/'source.json')
        self.assertGreater(source['maximum_drag_lift_CL_offset'],8e-5)
        self.assertLess(source['maximum_drag_lift_CL_offset'],1e-4)
        for target in source['targets']:
            self.assertNotEqual(target['experimental_marker_CL'],target['CL'])

    def test_nominal_uses_production_paneling_and_finite_edges(self):
        nominal = elements('nominal')
        for (name, xy),count in zip(nominal,(340,240)):
            source=read_dat(ROOT/run.QUAL/f'report_test_{name}.dat')
            np.testing.assert_array_equal(xy,repanel(source,count))
            self.assertGreater(np.linalg.norm(xy[0]-xy[-1]),.00099)
            self.assertEqual(xy.shape,(count+1,2))

    def test_displacement_is_rigid_and_only_on_aft(self):
        nominal = dict(elements('nominal'))
        for scenario,sign in (('aft_z_minus',-1),('aft_z_plus',1)):
            changed = dict(elements(scenario))
            np.testing.assert_array_equal(changed['fore'],nominal['fore'])
            np.testing.assert_allclose(changed['aft']-nominal['aft'],
                np.tile([0.,sign*SHIFT],(len(nominal['aft']),1)),atol=1e-16,rtol=0)
        self.assertAlmostEqual(SHIFT*457.2,.13)

    def test_cubic_preserves_edges_and_has_bounded_contour_change(self):
        for (name,xy),(_,base) in zip(elements('source_cubic'),elements('nominal')):
            source=read_dat(ROOT/run.QUAL/f'report_test_{name}.dat')
            np.testing.assert_array_equal(xy[[0,-1]],source[[0,-1]])
            shift=np.linalg.norm(xy-base,axis=1).max()
            self.assertGreater(shift,1e-5)
            self.assertLess(shift,.00014)
        with self.assertRaises(ValueError): elements('unregistered')
        with self.assertRaises(ValueError): cubic_repanel(source,10)

    def test_study_matrix_conditions_and_frozen_sources(self):
        spec = run.specification()
        self.assertEqual(len(spec['jobs']),9)
        self.assertEqual(sum(j['control'] for j in spec['jobs']),3)
        self.assertEqual((spec['reynolds'],spec['mach'],spec['ncrit'],spec['reference_chord']), (1e6,.1,9.,1.))
        for scenario in SCENARIOS:
            jobs=[j for j in spec['jobs'] if j['scenario']==scenario and 'repeat_of' not in j]
            self.assertEqual({j['target_cl'] for j in jobs},{.284,1.041})
            np.testing.assert_allclose([j['center'] for j in jobs],[-1.5717960133708968,5.063277366434972])

    def test_resume_refuses_changed_protocol(self):
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp); spec=dict(sha256={},ncrit=9.)
            run.prepare(out,spec)
            with self.assertRaises(ValueError): run.prepare(out,dict(spec,ncrit=6.25))

    def test_repeat_gate_requires_arrays_and_path_not_just_drag(self):
        spec=dict(jobs=[dict(name='repeat',repeat_of='original')],acceptance_control_CD_counts=.01)
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp)
            for name in ('original','repeat'):
                folder=out/name; folder.mkdir()
                run.write(folder/'result.json',dict(feasible=True,CD=.01,selected_path='cold_local'))
                for element in ('fore','aft'): np.savez(folder/(element+'_state.npz'),U=[1.,2.])
            self.assertTrue(run.check_controls(out,spec)['passed'])
            np.savez(out/'repeat/aft_state.npz',U=[1.,2.000000001])
            self.assertFalse(run.check_controls(out,spec)['passed'])
            np.savez(out/'repeat/aft_state.npz',U=[1.,2.])
            run.write(out/'repeat/result.json',dict(feasible=True,CD=.01,selected_path='warm_anchor'))
            self.assertFalse(run.check_controls(out,spec)['passed'])

    def test_failed_target_cannot_contain_selected_force(self):
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp); (out/'failed').mkdir()
            run.write(out/'failed/result.json',dict(feasible=False,CL=None,CD=.01,CM=None,alpha_deg=None))
            with self.assertRaisesRegex(ValueError,'Failed target'):
                run.summarize(out,dict(jobs=[dict(name='failed')]))


if __name__=='__main__': unittest.main()
