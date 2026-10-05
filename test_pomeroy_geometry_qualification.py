"""Extraction, transformation and placement-consistency regressions; no CFD."""
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
import xml.etree.ElementTree as ET

import numpy as np

sys.path.insert(0,str(Path(__file__).parent/'validation'))
import qualify_pomeroy_geometry as q


class PomeroyGeometryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.bbox=q.ROOT/'validation/results/pomeroy_geometry_qualification_20261001/coordinates_bbox.html'
        cls.local,cls.rows=q.extract(cls.bbox)

    def test_split_tables_signs_and_endpoints_match_printed_rows(self):
        self.assertEqual(len(self.rows),243)
        for n,te in zip(q.NAMES,(.0009,.0029,.0033)):
            self.assertEqual(self.local[n].shape,(81,2))
            np.testing.assert_array_equal(self.local[n][[0,-1]],[[1.,te],[1.,-te]])
        np.testing.assert_array_equal(self.local['main'][73:75],[[.8951,.0108],[.9219,.0088]])
        np.testing.assert_array_equal(self.local['flap1'][59:61],[[.3327,-.0565],[.3700,-.0506]])
        np.testing.assert_array_equal(self.local['flap2'][45:47],[[.0110,-.0310],[.0206,-.0388]])

    def test_table1_uses_clockwise_absolute_deflections_without_renormalization(self):
        assembled=q.assemble(self.local)
        np.testing.assert_allclose(.5*(assembled['main'][0]+assembled['main'][-1]),[.7028,0.],atol=1e-15)
        tip=.5*(assembled['flap2'][0]+assembled['flap2'][-1])
        self.assertAlmostEqual(tip[0],1.000005249185088,places=12)
        self.assertLess(tip[1],-.27)
        self.assertNotEqual(float(tip[0]),1.)

    def test_gap_uses_distance_to_segment_not_nearest_node(self):
        distance,foot=q.point_to_contour(np.array([2.,.5]),np.array([[0.,0.],[1.,0.],[1.,1.],[0.,1.]]))
        self.assertEqual(distance,1.)
        np.testing.assert_array_equal(foot,[1.,.5])

    def test_published_gap_mismatch_survives_TE_choice_and_resampling(self):
        rows=q.rigging_checks(self.local)
        self.assertEqual(len(rows),18)
        self.assertTrue(all(r['gap']>.015 for r in rows))
        a=next(r for r in rows if r['first']=='main' and r['target_panels']==80 and r['TE_definition']=='midpoint')
        self.assertAlmostEqual(a['gap'],.023140176738932247,places=12)
        self.assertAlmostEqual(a['nominal_overhang'],-.0168,places=12)
        self.assertTrue(all(r['gap_minus_table2']>.005 for r in rows))

    def test_missing_coordinate_word_cannot_shift_pairing_silently(self):
        tree=ET.parse(self.bbox)
        word=next(w for w in tree.findall('.//h:word',q.NS) if w.text=='1.0000')
        word.text='unreadable'
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'bad.html';tree.write(path)
            with self.assertRaisesRegex(ValueError,'Coordinate block length differs'):q.extract(path)

    def test_changed_primary_source_is_rejected(self):
        with patch.object(q,'sha',return_value='changed'):
            with self.assertRaisesRegex(ValueError,'Primary PDF changed'):q.qualify(Path('/unused'))


if __name__=='__main__':unittest.main()
