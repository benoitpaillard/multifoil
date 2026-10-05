"""Independent analytic and admission checks for the pressure-load diagnostic."""
import sys
from pathlib import Path
import unittest
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent/'validation'))
import report_s414_pressure_budget as p


class PressureBudgetTests(unittest.TestCase):
    def test_constant_pressure_on_closed_contour_has_zero_force(self):
        f = p.loads([1, 0, 0, 1], [1, 1, 0, 0], [3, 3, 3, 3], closed=True)
        np.testing.assert_allclose(f, 0., atol=1e-14)

    def test_flat_plate_lift_and_side_sign(self):
        upper = p.integrate([0, 1], [0, 0], [0, 1], [-2, -2], 0, 1, 'upper')
        lower = p.integrate([0, 1], [0, 0], [0, 1], [.5, .5], 0, 1, 'lower')
        np.testing.assert_allclose(upper+lower, [2.5, 0])
        self.assertAlmostEqual(p.lift(upper+lower, 60), 1.25)

    def test_sloped_geometry_exact_linear_cp(self):
        # Cp=2x+1, y=3x, lower branch: integral Cp dx=2, -integral Cp dy=-6.
        np.testing.assert_allclose(p.integrate([0,.2,.8,1], [0,.6,2.4,3],
            [0,.7,1], [1,2.4,3], 0, 1, 'lower'), [2,-6], atol=1e-13)

    def test_geometry_kink_is_preserved_between_taps(self):
        # y=x on [0,.5], y=1-x on [.5,1]; Cp=x => CC=.25.
        np.testing.assert_allclose(p.integrate([0,.5,1], [0,.5,0], [0,1], [0,1],
            0,1,'lower'), [.5,.25], atol=1e-14)

    def test_partial_and_missing_partition_closes(self):
        x=[0,.2,.7,1]; y=[0,.1,.03,0]; cp=[1,-2,-.5,0]
        parts=sum((p.integrate(x,y,x,cp,a,b,'upper') for a,b in [(0,.13),(.13,.91),(.91,1)]), np.zeros(2))
        np.testing.assert_allclose(parts,p.integrate(x,y,x,cp,0,1,'upper'),atol=1e-14)

    def test_report_lift_conversion_roundtrip(self):
        alpha=-3.074; cl=.0556; cd=.006177
        cn=p.experimental_normal(cl,cd,alpha); a=np.deg2rad(alpha)
        self.assertAlmostEqual(cn/np.cos(a)-cd*np.tan(a),cl)

    def test_upper_negative_weight_conditional_range(self):
        np.testing.assert_allclose(p.signed_bound(-.01,-2,1),[-.01,.02])

    def test_unsupported_extrapolation_rejected(self):
        with self.assertRaises(ValueError): p.interp([.1,1],[1,0],[0])
        with self.assertRaises(ValueError): p.integrate([0,1],[0,0],[.1,1],[0,0],0,1,'upper')

    def test_reversed_and_nonfinite_curves_rejected(self):
        for x,y in [([1,0],[0,0]),([0,0],[0,0]),([0,1],[0,np.nan])]:
            with self.assertRaises(ValueError): p.curve(x,y)

    def test_duplicate_endpoint_excluded_without_creating_fake_coverage(self):
        rows=[dict(included='True',table_index=str(i),x_over_c=str(x),experiment_Cp='0') for i,x in enumerate([.1,.9])]
        rows.append(dict(included='False',table_index='',x_over_c='1',experiment_Cp='0'))
        kept,x,_=p.physical_taps(rows)
        self.assertEqual(len(kept),2); self.assertEqual(x[-1],.9)

    def test_interior_missing_tap_cannot_be_silently_bridged(self):
        rows=[dict(included='True',table_index=str(i),x_over_c=str(x),experiment_Cp='0') for i,x in [(0,.1),(2,.9)]]
        with self.assertRaises(ValueError): p.physical_taps(rows)

    def test_normal_force_independent_of_ordinate_mapping(self):
        a=p.integrate([0,1],[0,0],[0,.4,1],[-1,.2,0],0,1,'upper')
        b=p.integrate([0,.7,1],[0,.1,0],[0,.4,1],[-1,.2,0],0,1,'upper')
        self.assertAlmostEqual(a[0],b[0]); self.assertNotEqual(a[1],b[1])


if __name__ == '__main__':
    unittest.main()
