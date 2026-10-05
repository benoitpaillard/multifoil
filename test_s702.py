"""Independent geometry and reference contracts for the S702 reconstruction."""
import json
from pathlib import Path
import sys
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT/'validation'))
from s702_geometry import geometry, DATA
from optimization import minimum_clearance


class S702DataTests(unittest.TestCase):
    def test_hinge_and_deflection_preserve_upstream_elements(self):
        audit = json.loads((DATA/'geometry_extraction.json').read_text())
        np.testing.assert_allclose(audit['hinge_from_cross'], [.9001, -.0114], atol=6e-5)
        baseline = geometry(0)
        hinge = np.array([.9001, -.0114])
        for angle in [-2, 0, 5, 10]:
            g = geometry(angle)
            np.testing.assert_array_equal(g[0][1], baseline[0][1])
            upstream = g[1][1][g[1][1][:, 0] < .85]
            reference = baseline[1][1][baseline[1][1][:, 0] < .85]
            # Repeated drawing translations accumulate up to a few 0.001-pt
            # coordinate roundoffs; keep these plotted differences unmodified.
            np.testing.assert_allclose(upstream, reference, atol=5e-6)
            a = np.deg2rad(-angle)
            rotation = np.array([[np.cos(a), -np.sin(a)], [np.sin(a), np.cos(a)]])
            expected = (baseline[1][1][[0, -1]]-hinge) @ rotation.T+hinge
            np.testing.assert_allclose(g[1][1][[0, -1]], expected, atol=1.5e-5)
            self.assertGreater(minimum_clearance(g[0][1], g[1][1]), .019)

    def test_reference_calibration_matches_printed_maximum_lift(self):
        # Values independently stated in the report, not fitted to the extraction.
        for angle, printed_clmax in [(0, 1.80), (10, 2.01)]:
            rows = np.genfromtxt(DATA/f'd{angle:+d}'/'experiment_lift.csv', delimiter=',', names=True)
            self.assertLess(abs(rows['CL'].max()-printed_clmax), .001)
        experiment = np.genfromtxt(DATA/'d+0'/'experiment_lift.csv', delimiter=',', names=True)
        mses = np.genfromtxt(DATA/'d+0'/'mses_lift.csv', delimiter=',', names=True)
        overprediction = mses['CL'].max()/experiment['CL'].max()-1
        self.assertTrue(.10 < overprediction < .13)

    def test_drag_markers_correspond_to_lift_measurement_sequence(self):
        # A wrong axis or symbol block can still look like a plausible polar.
        # Match the independently drawn lift values for the same measurement order.
        for angle in [-2, 0, 5, 10]:
            folder = DATA/f'd{angle:+d}'
            drag = np.genfromtxt(folder/'experiment_drag.csv', delimiter=',', names=True)
            lift = np.genfromtxt(folder/'experiment_lift.csv', delimiter=',', names=True)
            # Separate copies of the symbols have subpixel translation offsets
            # (up to 1.03e-4 CL in this PDF), not distinct force measurements.
            np.testing.assert_allclose(drag['CL'], lift['CL'][:len(drag)], atol=2e-4)
            self.assertTrue(np.all(np.diff(drag['CL']) > 0))
            self.assertTrue(np.all(drag['CD'] > 0))


if __name__ == '__main__':
    unittest.main()
