"""Pressure-extraction checks against independent report information."""
import hashlib
import json
from pathlib import Path
import unittest

import numpy as np

DATA = Path(__file__).resolve().parent/'validation/data/s702'


class S702PressureDataTests(unittest.TestCase):
    def test_fore_upper_taps_match_table_ii(self):
        # These are the printed Table II(a) locations, independent of the Cp
        # figure and its axis transformation. All 36 fore upper taps are shown.
        printed = np.array([.00026, .00198, .00682, .01462, .02564, .03961,
            .05661, .07680, .09926, .12465, .15099, .18048, .21118, .24364,
            .27736, .31247, .34791, .38394, .41894, .45569, .49061, .52582,
            .55998, .59183, .62247, .65136, .67858, .70403, .72663, .74583,
            .76395, .77766, .78828, .79672, .80154, .80270])
        rows = np.genfromtxt(DATA/'pressure_cl1.04/experiment_fore_upper.csv',
                            delimiter=',', names=True)
        np.testing.assert_allclose(rows['x_over_cref'], printed, atol=5e-5, rtol=0)

    def test_aft_upper_markers_match_subset_of_printed_taps(self):
        # Some pressure orifices are omitted from this plotted comparison.
        printed = np.array([.69941, .70194, .70930, .72177, .73844, .75893,
            .78231, .80843, .83542, .84940, .86303, .87212, .88138, .88928,
            .90395, .91666, .93969, .95997, .97649, .98892, .99660, 1.])
        rows = np.genfromtxt(DATA/'pressure_cl1.04/experiment_aft_upper.csv',
                            delimiter=',', names=True)
        distance = abs(rows['x_over_cref'][:, None]-printed[None, :]).min(axis=1)
        self.assertEqual(len(rows), 18)
        self.assertLess(float(distance.max()), 8e-5)

    def test_caption_and_pressure_provenance(self):
        source = json.loads((DATA/'pressure_cl1.04/source.json').read_text())
        self.assertEqual((source['target_cl'], source['reynolds'], source['mach']),
                         (1.04, 1e6, .1))
        self.assertEqual(hashlib.sha256((DATA/'NASA-CR-20220014964.pdf').read_bytes()).hexdigest(),
                         source['pdf_sha256'])
        for name, sha in source['files_sha256'].items():
            self.assertEqual(hashlib.sha256((DATA/'pressure_cl1.04'/name).read_bytes()).hexdigest(), sha)
        # The separately plotted lift at the reported pressure-test incidence
        # rounds to the figure's CL=1.04; no extra precision is invented.
        lift = np.genfromtxt(DATA/'d+0/experiment_lift.csv', delimiter=',', names=True)
        row = lift[np.argmin(abs(lift['alpha_deg']-6.10))]
        self.assertEqual(round(float(row['CL']), 2), source['target_cl'])


if __name__ == '__main__':
    unittest.main()
