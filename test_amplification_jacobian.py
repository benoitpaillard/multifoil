"""Derivative-only correctness, isolation and branch coverage without flow solves."""
import copy
from pathlib import Path
import sys
from types import SimpleNamespace
import unittest
from unittest.mock import patch
import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT/'validation'))
import amplification_jacobian as option


class AmplificationJacobianTests(unittest.TestCase):
    def setUp(self):
        model = SimpleNamespace(geom=option.mf.Geom(), oper=option.mf.Oper(), param=option.mf.Param())
        model.geom.chord = 1.; model.oper.Re = 1e6; model.oper.Ma = .1; model.param.ncrit = 9
        option.mf.init_thermo(model)
        self.param = model.param

    def state(self):
        with np.load(ROOT/'validation/results/s414_lower_bucket_20261003/cl0056_nominal/fore_state.npz') as saved:
            ids = saved['lower']; lam = ids[saved['turbulent'][ids] == 0]
            return saved['U'][:, lam[-1]].copy()

    def test_legacy_and_untagged_are_exact_even_after_opt_in(self):
        U = self.state(); a = SimpleNamespace(param=copy.deepcopy(self.param)); b = SimpleNamespace(param=copy.deepcopy(self.param))
        option.configure([a], option.CORRECTED); option.configure([b], option.LEGACY)
        expected = option._original_get_damp(U, self.param)
        for p in (self.param, b.param):
            actual = option.mf.get_damp(U, p)
            self.assertEqual(actual[0], expected[0]); np.testing.assert_array_equal(actual[1], expected[1])
        self.assertNotEqual(option.mf.get_damp(U, a.param)[1][1], expected[1][1])

    def test_missing_factor_corrected_without_value_or_state_change(self):
        U = self.state(); before = U.copy(); self.param.amplification_jacobian = option.CORRECTED
        value, derivative = option.get_damp(U, self.param)
        self.assertEqual(value, option._original_get_damp(U, self.param)[0])
        for relative in (1e-5, 1e-6):
            for j in range(4):
                h = relative*max(abs(U[j]), 1e-6); a, b = U.copy(), U.copy(); a[j] += h; b[j] -= h
                fd = (option.get_damp(a, self.param)[0]-option.get_damp(b, self.param)[0])/(2*h)
                self.assertLess(abs(fd-derivative[j])/max(1, abs(derivative[j])), 2e-5)
        np.testing.assert_array_equal(U, before)

    def test_copy_and_side_parameters_preserve_mode(self):
        model = SimpleNamespace(param=self.param)
        option.configure([model], option.CORRECTED)
        for si in (0, 1, 2):
            self.assertEqual(option.mf.build_param(copy.deepcopy(model), si).amplification_jacobian, option.CORRECTED)

    def test_inactive_and_hk_floor_preserve_derivative(self):
        self.param.amplification_jacobian = option.CORRECTED
        for U in (np.array([1e-7, 2e-7, 0., .8]), np.array([.001, .00101, 8.9, .8])):
            expected = option._original_get_damp(U, self.param); actual = option.get_damp(U, self.param)
            self.assertEqual(actual[0], expected[0]); np.testing.assert_array_equal(actual[1], expected[1])

    def test_activation_ramp_and_full_growth(self):
        self.param.amplification_jacobian = option.CORRECTED
        for rn in (.25, .75, 1.5):
            hmi = 1/5.; critical = 2.492*hmi**.43+.7*(np.tanh(14*hmi-9.24)+1)
            U = np.array([.001, .006, 2., .8]); ret = 10**(critical-.1+.2*rn)
            with patch.object(option.mf, 'get_Hk', return_value=(6., np.array([0., 1., 0., 0.]))), \
                 patch.object(option.mf, 'get_Ret', return_value=(ret, np.zeros(4))):
                old = option._original_get_damp(U, self.param); new = option.get_damp(U, self.param)
            rf = 1 if rn >= 1 else 3*rn**2-2*rn**3
            expected = rf*(.028*5-.0345*np.exp(-(3.87*hmi-2.52)**2))/U[0]*np.exp(-20*hmi)*hmi**2
            self.assertEqual(new[0], old[0]); self.assertAlmostEqual(new[1][1]-old[1][1], expected, places=12)

    def test_invalid_mode_rejected_before_constructing_solver(self):
        with self.assertRaises(ValueError): option.AmplificationJacobianSolver([], amplification_jacobian='typo')
        with self.assertRaises(ValueError): option.configure([], 'typo')


if __name__ == '__main__':
    unittest.main()
