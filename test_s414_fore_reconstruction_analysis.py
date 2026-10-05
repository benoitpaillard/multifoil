"""Coverage, compressibility and transition-accounting regression checks."""
from pathlib import Path
import json
import sys
import unittest
import numpy as np
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT/'validation'))
import report_s414_fore_reconstruction as report


class ForeAnalysisTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.folder = report.STUDY/'cl0056_fore_source_cubic'
        cls.result = json.loads((cls.folder/'result.json').read_text())
        cls.surface = np.genfromtxt(cls.folder/'flow/surface_004_fore.csv', delimiter=',', names=True, dtype=None, encoding='utf-8')
        cls.model = report.new_model('fore', report.geometry_elements('fore_source_cubic', [340,240])[0][1], 1e6,.1,9,(.25,0),1e-5)
        cls.layer = cls.result['boundary_layers']['fore']['lower']

    def test_missing_target_never_supplies_paired_drag(self):
        with self.assertRaisesRegex(ValueError, 'Incomplete target'):
            report.paired_drag(dict(feasible=True,CD=.01), dict(feasible=False,CD=None))
        with self.assertRaisesRegex(ValueError, 'Incomplete target'):
            report.paired_drag(dict(feasible=False,CD=.01), dict(feasible=True,CD=.007))

    def test_pair_sign_is_new_minus_old(self):
        self.assertAlmostEqual(report.paired_drag(dict(feasible=True,CD=.01),dict(feasible=True,CD=.009)), -10)

    def test_velocity_inverse_matches_solver_for_zero_and_nonzero_mach(self):
        for mach in (0.,.1):
            model = report.new_model('fore', report.geometry_elements('nominal',[340,240])[0][1],1e6,mach,9,(.25,0),1e-5)
            param = report.mf.build_param(model,0)
            u = np.array([1e-4,.2,1.,1.8])
            np.testing.assert_allclose(report.raw_speed(report.mf.get_uk(u,param)[0],param),u,rtol=3e-15,atol=1e-15)

    def test_export_growth_closes_and_recovers_recorded_transition(self):
        original = self.surface.copy()
        rows, stats = report.lower_history(self.model, self.surface, self.layer)
        self.assertAlmostEqual(stats['base_N']+stats['regularization_N'],9,places=7)
        self.assertLess(stats['maximum_growth_balance'],1e-7)
        self.assertEqual(stats['fore_transition_x'],self.layer['transition']['xy_over_cref'][0])
        self.assertTrue(all(r['N'] is None for r in rows if not r['laminar']))
        np.testing.assert_array_equal(original,self.surface)

    def test_corrupt_amplification_is_detected(self):
        changed = self.surface.copy()
        j = np.flatnonzero((changed['side']=='lower') & (changed['turbulent']==0))[-1]
        changed['amplification_or_shear'][j] += .1
        with self.assertRaises(ValueError): report.lower_history(self.model,changed,self.layer)

    def test_corrupt_velocity_cannot_masquerade_as_matching_pressure(self):
        changed = self.surface.copy()
        j = np.flatnonzero(changed['side']=='lower')[5]
        changed['ue_over_vinf'][j] *= 1.01
        with self.assertRaisesRegex(ValueError,'Pressure inversion mismatch'):
            report.lower_history(self.model,changed,self.layer)


if __name__ == '__main__': unittest.main()
