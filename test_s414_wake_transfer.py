"""Wake-transfer diagnostic contracts; no flow solves."""
import ast
import inspect
from pathlib import Path
import sys
import unittest
from types import SimpleNamespace
import numpy as np

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'validation'))
import report_s414_wake_transfer as report
from vendor import mfoil as mf


class WakeTransferTests(unittest.TestCase):
    def test_nearest_open_segments_and_interpolation(self):
        wake=np.array([[0,0],[2,0],[2,2]],float)
        distance,index,fraction,projected=report.nearest_polyline([[1,1],[3,2],[-1,0]],wake)
        np.testing.assert_allclose(distance,[1,1,1])
        self.assertEqual(index[0],0)
        np.testing.assert_allclose(projected,[[1,0],[2,2],[0,0]])
        np.testing.assert_allclose(report.segment_values([10,20,40],index,fraction),[15,40,10])
        # Closing the path would falsely put this point on a wake segment.
        d,*_=report.nearest_polyline([[1,1]],[[0,0],[0,2],[2,2]])
        self.assertEqual(d[0],1.)

    def test_zero_length_wake_segment_rejected(self):
        with self.assertRaises(ValueError):report.nearest_polyline([[0,1]],[[0,0],[0,0]])

    def test_source_column_partition_preserves_sign_and_total(self):
        models=[]
        for name in ('fore','aft'):
            # 2 surface + 2 wake nodes give 2 source columns, one in each group.
            models.append(SimpleNamespace(geom=SimpleNamespace(name=name),foil=SimpleNamespace(N=2),
                wake=SimpleNamespace(N=2),isol=SimpleNamespace(sgnue=np.array([-1.,1.])),
                glob=SimpleNamespace(U=np.array([[1,1,1,1],[1,2,3,4],[1,1,1,1],[2,2,2,2]],float)),
                vsol=SimpleNamespace(sigma_m=np.array([[-1,1,0,0],[0,0,-1,1]],float))))
        operator=np.arange(32,dtype=float).reshape(8,4)/10
        mass=np.zeros((4,8));mass[:2,:4]=models[0].vsol.sigma_m;mass[2:,4:]=models[1].vsol.sigma_m
        signs=np.array([-1,1,1,1,-1,1,1,1])
        solver=SimpleNamespace(models=models,nsys=8,ss=[slice(0,4),slice(4,8)],
            sources=[slice(0,2),slice(2,4)],ue_sigma=operator,D=signs[:,None]*(operator@mass))
        parts=report.source_parts(solver)
        self.assertEqual(set(parts),{'fore_surface','fore_wake','aft_surface','aft_wake'})
        np.testing.assert_allclose(parts['fore_wake'],signs*operator[:,1]*2)
        np.testing.assert_allclose(sum(parts.values()),solver.D@np.tile([2,4,6,8],2))

    def test_wake_closure_flag_is_an_assumption_not_TE_observation(self):
        model=SimpleNamespace(param=mf.Param())
        self.assertFalse(mf.build_param(model,0).turb)
        self.assertFalse(mf.build_param(model,1).turb)
        self.assertTrue(mf.build_param(model,2).turb)
        self.assertTrue(mf.build_param(model,2).wake)

    def test_wake_initialization_uses_only_its_own_TE_states(self):
        code=inspect.getsource(mf.wake_sys)
        self.assertIn('get_cttr(Ul',code)
        self.assertIn('get_cttr(Uu',code)
        # Both source layers and the target wake belong to the same model M.
        self.assertIn('M.vsol.Is[0][-1]',code)
        self.assertIn('M.vsol.Is[1][-1]',code)
        self.assertIn('M.vsol.Is[2][ 0]',code)

    def test_diagnostics_never_solve_or_compute_perturbed_forces(self):
        tree=ast.parse(inspect.getsource(report))
        calls={n.func.attr for n in ast.walk(tree) if isinstance(n,ast.Call) and isinstance(n.func,ast.Attribute)}
        self.assertFalse({'solve','_newton','_initialize','calc_force'}&calls)
        self.assertIn('_system',calls)


if __name__=='__main__':unittest.main()
