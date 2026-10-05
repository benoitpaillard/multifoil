"""Regression for the TN808 slot-entry failure and unchanged smooth inputs."""
import sys
import unittest
from pathlib import Path
import numpy as np
from multifoil import naca4,read_dat,repanel,validate_geometry
from viscous import new_model
from paneling import curvature_paneling
from vendor import mfoil as mf
sys.path.insert(0,str(Path(__file__).parent/'validation'))
from tn808_geometry import geometry,geometry_audit


class PanelingTests(unittest.TestCase):
    def test_smooth_inputs_are_bitwise_unchanged(self):
        inputs=[naca4('0012',141).T]
        data=Path(__file__).parent/'validation/data/nlr7301'
        inputs += [read_dat(data/name).T for name in ('NLR7301MainECARP.dat','NLR7301FlapECARP.dat')]
        for xy in inputs:
            original=mf.spline_curvature(xy,181,2.,.1,None)
            fixed=curvature_paneling(xy,181)
            self.assertFalse(fixed[3])
            for before,after in zip(original,fixed[:3]):np.testing.assert_array_equal(before,after)

    def test_slot_entry_has_ordered_nonintersecting_solver_panels(self):
        raw,_=geometry(20)
        xy=repanel(raw[0][1],260)
        te=.5*(xy[0]+xy[-1]);le=xy[np.argmax(np.linalg.norm(xy-te,axis=1))]
        chord=np.linalg.norm(te-le);axis=(te-le)/chord
        frame=np.column_stack([axis,[-axis[1],axis[0]]])
        local=frame.T@(xy[::-1].T-le[:,None])/chord
        original=mf.spline_curvature(local,len(xy),2.,.1,None)
        self.assertTrue(np.any(np.diff(original[1])<=0))
        model=new_model('main',xy,2190000,.105,9,np.array([.25,0]),1e-5)
        self.assertTrue(model.paneling_recovery)
        self.assertTrue(np.all(np.diff(model.foil.s)>0))
        validate_geometry([('main',model.foil.x.T)])
        np.testing.assert_allclose(model.foil.x[:,[0,-1]],xy[[-1,0]].T,atol=1e-12)

    def test_reconstruction_table_and_chord(self):
        self.assertLess(geometry_audit()['corrected_main_table_max_abs_error'],1e-7)
        for angle in (10,20,30):
            raw,audit=geometry(angle)
            self.assertAlmostEqual(audit['physical_chord_angle_deg'],-angle,places=10)
            self.assertAlmostEqual(audit['nose'][0],.8975,places=12)
            upper=raw[1][1][:301][::-1]
            self.assertAlmostEqual(.0159-np.interp(.9,upper[:,0],upper[:,1]),.02,places=12)


if __name__=='__main__':unittest.main()
