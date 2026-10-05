"""Numerical/physical checks for the coupled viscous backend (no paid tools)."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
import copy
import math
import unittest
import numpy as np

from multifoil import Solution, naca4
from viscous import (CoupledSolver, analyze_stall, dense_boundary_system,
                     source_stream_matrix, source_velocity_matrix, mf)


class ViscousPhysics(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.xy=naca4('0012',141)
        cls.solver=CoupledSolver([('main',cls.xy)])
        cls.zero=cls.solver.solve(0)
        cls.result=cls.solver.solve(2)
        assert cls.result.converged,cls.result.diagnostics

    def test_symmetry_positive_drag_and_residuals(self):
        self.assertLess(abs(self.zero.cl),1e-8)
        self.assertGreater(self.result.cl,.15)
        self.assertGreater(self.result.cd,.003)
        self.assertLess(self.result.cd,.01)
        self.assertLess(self.result.residual,1e-7)
        self.assertTrue(np.all(self.result.surfaces[0]['theta']>0))

    def test_boundary_layer_residual_and_jacobian_match_upstream(self):
        m=copy.deepcopy(self.solver.models[0])
        r,j=dense_boundary_system(m)
        m.glob.realloc=True
        mf.build_glob_sys(m)
        np.testing.assert_allclose(r,m.glob.R,atol=1e-12)
        np.testing.assert_allclose(j,m.glob.R_U.toarray(),atol=1e-10,rtol=1e-11)

    def test_single_element_against_native_mfoil(self):
        # Same actual panel geometry and closures; the reference retains the
        # original wake-source approximation and original Newton algorithm.
        m=copy.deepcopy(self.solver.models[0])
        m.oper.initbl=True;m.param.rtol=1e-7;m.param.niglob=80
        original=mf.build_glob_sys
        def compatible(model):
            model.glob.realloc=True # upstream sparse shape comparison fix
            original(model)
        try:
            mf.build_glob_sys=compatible
            mf.solve_viscous(m)
        finally:mf.build_glob_sys=original
        self.assertTrue(m.glob.conv)
        self.assertLess(abs(self.result.cl/m.post.cl-1),.005)
        self.assertLess(abs(self.result.cd/m.post.cd-1),.002)

    def test_cross_body_sources_against_independent_velocity_integral(self):
        # Integrate analytic normal velocity along a second closed contour.
        # This catches source branch jumps AND missing wake half-panel terms.
        m=self.solver.models[0]
        targets=m.foil.x+np.array([[0.],[.4]])
        psi=source_stream_matrix(m,targets,cross_body=True)
        dx=np.diff(targets,axis=1);length=np.linalg.norm(dx,axis=0)
        normals=np.vstack((-dx[1],dx[0]))/length
        q,w=np.polynomial.legendre.leggauss(4)
        integral=np.zeros((len(length),psi.shape[1]))
        for a,b in zip(q,w):
            pts=targets[:,:-1]+.5*(a+1)*dx
            integral-=.5*b*length[:,None]*source_velocity_matrix(m,pts,normals)
        np.testing.assert_allclose(np.diff(psi,axis=0),integral,atol=2e-10,rtol=2e-7)

    def test_reynolds_and_transition_sensitivity(self):
        earlier=CoupledSolver([('main',self.xy)],ncrit=5).solve(2)
        higher=CoupledSolver([('main',self.xy)],reynolds=2e6).solve(2)
        self.assertTrue(earlier.converged and higher.converged)
        def transition(r):
            return r.diagnostics['boundary_layers']['main']['upper']['transition']['x_over_local_chord']
        self.assertLess(transition(earlier),transition(self.result)-.05)
        self.assertGreater(earlier.cd,self.result.cd)
        self.assertLess(higher.cd,self.result.cd)

    def test_multiple_elements_order_and_far_limit(self):
        xy=naca4('0012',81)
        solo=CoupledSolver([('a',xy)]).solve(2)
        e=[('a',xy),('b',xy+[0,100])]
        forward=CoupledSolver(e,tolerance=1e-9).solve(2)
        reverse=CoupledSolver(e[::-1],tolerance=1e-9).solve(2)
        self.assertTrue(solo.converged and forward.converged and reverse.converged)
        self.assertLess(abs(forward.cl/(2*solo.cl)-1),.002)
        self.assertLess(abs(forward.cd/(2*solo.cd)-1),.002)
        # At 100 chords the streamfunction system loses digits to its large
        # constant offsets. Require agreement well below one force-count.
        self.assertLess(abs(forward.cl-reverse.cl),5e-7)
        self.assertAlmostEqual(forward.cd,reverse.cd,places=8)

    def test_local_chord_and_reynolds_scaling(self):
        half=CoupledSolver([('main',.5*self.xy)],moment_reference=(.125,0),reynolds=2e6,wake_length=.5).solve(2)
        self.assertTrue(half.converged)
        self.assertLess(abs(half.cl/(.5*self.result.cl)-1),.001)
        self.assertLess(abs(half.cd/(.5*self.result.cd)-1),.001)

    def test_failure_is_not_a_force_or_a_stall_prediction(self):
        failed=CoupledSolver([('main',naca4('0012',81))],iterations=1).solve(0)
        self.assertFalse(failed.converged)
        self.assertTrue(math.isnan(failed.cl) and math.isnan(failed.cd))
        self.assertFalse(failed.surfaces)
        def point(a,cl,ok=True):
            return Solution(a,cl,0,.01,None,ok,0,[],{})
        good=[point(0,.5),point(1,1.),point(2,.8)]
        self.assertEqual(analyze_stall(good)['status'],'peak_bracketed')
        self.assertEqual(analyze_stall(good[:2]+[point(2,math.nan,False)])['status'],'not_bracketed')
        self.assertEqual(analyze_stall(good[:2])['status'],'not_bracketed')


if __name__=='__main__':unittest.main()
