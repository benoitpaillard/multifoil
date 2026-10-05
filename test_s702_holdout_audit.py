"""Panel-report geometry, resolution and job-configuration audit regressions."""
from pathlib import Path
import sys
import unittest
import numpy as np

sys.path.insert(0,str(Path(__file__).parent/'validation'))
import report_s702_holdout_panels as report


class HoldoutAuditTests(unittest.TestCase):
    def test_extra_vertices_do_not_create_a_geometry_change(self):
        a=np.array([[0.,0.],[1.,0.],[1.,1.],[0.,1.]])
        b=np.array([[0.,0.],[.5,0.],[1.,0.],[1.,.5],[1.,1.],[.5,1.],[0.,1.]])
        self.assertEqual(report.contour_sample_distance(a,b),0.)
        self.assertAlmostEqual(report.contour_sample_distance(a,b+[0.,.1]),.1)

    def test_resolution_excludes_laminar_step_and_preserves_station_location(self):
        v=np.array([(0,10.,.01,.5),(1,6.,1.2,.7),(1,5.,1.16,.8),(1,4.5,1.18,.9)],
                   dtype=[('turbulent',int),('Hk',float),('ue_over_vinf',float),('x_over_local_chord',float)])
        r=report.resolution(v)
        self.assertEqual(r['max_Hk_drop'],1.)
        self.assertAlmostEqual(r['max_velocity_step'],.04)
        self.assertEqual((r['velocity_step_x0'],r['velocity_step_x1']),(.7,.8))
        self.assertEqual(report.resolution(v[:1])['max_Hk_drop'],0.)

    def test_control_requires_its_declared_panel_override(self):
        inputs=report.Inputs();folder=report.ROOT/'validation/results/s702_holdout_panels_20260922'
        spec=inputs.json(folder/'study.json');rt=inputs.json(folder/'runtime.json')
        job=next(j for j in spec['jobs'] if j['control'])
        self.assertTrue(report.audit(folder/job['name'],job,report.effective_spec(job,spec),rt,inputs)['result']['feasible'])
        with self.assertRaisesRegex(ValueError,'Export settings'):
            report.audit(folder/job['name'],job,spec,rt,inputs)


if __name__=='__main__':unittest.main()
