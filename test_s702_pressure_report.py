"""Guard geometric pressure mapping, unsupported taps and surface metrics."""
import sys
from pathlib import Path
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent/'validation'))
from report_s702_matched_pressure import pressure_branch, at_taps, metric


class PressureReportTests(unittest.TestCase):
    def contour(self):
        # Stagnation-based labels deliberately disagree with geometric surfaces
        # near the nose. A tap belongs to the same wall as incidence changes.
        return np.array([(1.,-.4,'upper'),(.5,-.8,'upper'),(0.,1.,'lower'),
                         (.5,.3,'lower'),(1.,-.4,'lower')],
            dtype=[('x_over_cref',float),('Cp',float),('side','U5')])

    def test_geometric_upper_includes_nose_even_if_bl_label_is_lower(self):
        c=self.contour();upper=pressure_branch(c,'upper')
        self.assertEqual(upper['Cp'].tolist(),[1.,-.8,-.4])
        c['side']='lower'
        np.testing.assert_array_equal(pressure_branch(c,'upper')['Cp'],upper['Cp'])

    def test_correct_wall_in_overlap_region(self):
        c=self.contour();tap=[dict(x_over_cref=.5,Cp=-.8,endpoint_duplicate=0,drawing_index=2)]
        self.assertAlmostEqual(at_taps(pressure_branch(c,'upper'),tap)[0]['error_Cp'],0.)
        self.assertAlmostEqual(at_taps(pressure_branch(c,'lower'),tap)[0]['error_Cp'],1.1)

    def test_interpolation_is_linear_and_does_not_extrapolate(self):
        c=self.contour();tap=[dict(x_over_cref=x,Cp=-.6,endpoint_duplicate=0,drawing_index=i)
            for i,x in enumerate((.75,1.00001))]
        rows=at_taps(pressure_branch(c,'upper'),tap)
        self.assertAlmostEqual(rows[0]['model_Cp'],-.6)
        self.assertIsNone(rows[1]['model_Cp']);self.assertIsNone(rows[1]['error_Cp'])
        self.assertEqual(rows[1]['exclusion'],'outside_model_support')

    def test_duplicate_cannot_count_towards_accuracy_or_coverage(self):
        c=self.contour();tap=[dict(x_over_cref=.5,Cp=-.8,endpoint_duplicate=d,drawing_index=d) for d in (0,1)]
        rows=at_taps(pressure_branch(c,'upper'),tap);m=metric(rows)
        self.assertEqual((m['scored_taps'],m['excluded_taps'],m['total_plotted_taps']),(1,1,2))
        self.assertEqual(rows[1]['exclusion'],'duplicate_endpoint')

    def test_nonmonotone_branch_is_not_silently_sorted(self):
        c=self.contour();c=np.insert(c,1,np.array((1.1,-.4,'upper'),dtype=c.dtype))
        with self.assertRaises(ValueError):pressure_branch(c,'upper')

    def test_empty_score_is_not_zero_error(self):
        with self.assertRaises(ValueError):metric([dict(included=False,error_Cp=None)])

    def test_bias_and_rmse_do_not_cancel_opposite_errors(self):
        m=metric([dict(included=True,error_Cp=e) for e in (-.2,.2)])
        self.assertAlmostEqual(m['bias'],0.);self.assertAlmostEqual(m['rmse'],.2)


if __name__=='__main__':unittest.main()
