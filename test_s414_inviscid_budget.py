"""Scope and pressure-convention regression checks for the inviscid diagnostic."""
from pathlib import Path
import sys
import unittest
from types import SimpleNamespace
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent/'validation'))
import report_s414_inviscid_budget as p


class InviscidBudgetTests(unittest.TestCase):
    def test_component_identity_with_opposing_errors(self):
        a,b,c=p.split_error([2,-1],[1,-2],[.5,-1.5])
        np.testing.assert_allclose(a+b,c)
        np.testing.assert_allclose(a,[1.5,.5]);np.testing.assert_allclose(b,[-1,-1])

    def test_split_rejects_broadcast_and_nonfinite(self):
        with self.assertRaises(ValueError):p.split_error([1,2],1,[0,0])
        with self.assertRaises(ValueError):p.split_error([np.nan],[0],[0])

    def test_pressure_uses_raw_speed_and_is_sign_invariant(self):
        u=np.array([0,.4,1,1.4]);m=.1
        np.testing.assert_allclose(p.cp_from_raw_speed(u,0),1-u*u)
        np.testing.assert_allclose(p.cp_from_raw_speed(u,m),p.cp_from_raw_speed(-u,m))
        beta=np.sqrt(1-m*m)
        self.assertAlmostEqual(p.cp_from_raw_speed([0],m)[0],1/(beta+m*m/(2*(1+beta))))

    def test_independent_formula_matches_vendor_convention(self):
        for m in [0,.1,.4]:
            param=SimpleNamespace(Vinf=1.,Minf=m,KTb=np.sqrt(1-m*m),KTl=m*m/(1+np.sqrt(1-m*m))**2)
            u=np.array([-.3,0.,.7,1.,1.5])
            np.testing.assert_allclose(p.cp_from_raw_speed(u,m),p.mf.get_cp(u,param)[0],atol=1e-14)

    def test_matched_lift_cannot_substitute_for_measured_angle(self):
        t=dict(comparison_basis='fixed_alpha',alpha_deg=-3.07,reference_caption_alpha=-3.07,matched_lift_CD_error_counts=None)
        p.admit_angle(t)
        for changed in [dict(comparison_basis='fixed_lift'),dict(alpha_deg=-3.46),dict(matched_lift_CD_error_counts=1)]:
            with self.assertRaises(ValueError):p.admit_angle(dict(t,**changed))

    def test_invalid_pressure_inputs_rejected(self):
        for u,m in [([np.nan],.1),([1],.7),([1],-.1),([100],.6)]:
            with self.assertRaises(ValueError):p.cp_from_raw_speed(u,m)


if __name__=='__main__':unittest.main()
