"""Scientific extraction contracts; no aerodynamic solves."""
import json
from pathlib import Path
import sys
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'validation'))
from qualify_s414_geometry import (OUT, parse_pts, comparison_geometry, geometry_metrics,
                                   point_distances, symmetric_vertex_distance, sha)
from extract_s414_pressure import pressure_coordinates


class S414QualificationTests(unittest.TestCase):
    def test_boundary_blocks_are_not_all_solid_airfoils(self):
        a=parse_pts(OUT/'sources/s414.pts')
        b=parse_pts(OUT/'sources/s414-1.pts')
        self.assertEqual([len(x['points']) for x in a],[39,50,32,23,27])
        self.assertEqual([len(x['points']) for x in b],[53,201,86,51,48,99])
        self.assertGreater(np.ptp(b[0]['points'][:,0]),19.)

    def test_global_chord_normalization_preserves_relative_placement(self):
        blocks=parse_pts(OUT/'sources/s414-1.pts')
        elems,scale=comparison_geometry(blocks,'s414-1.pts')
        self.assertAlmostEqual(scale,1.249952)
        np.testing.assert_allclose(elems['fore']*scale,blocks[1]['points'],atol=1e-15)
        np.testing.assert_allclose(elems['aft']*scale,blocks[2]['points'],atol=1e-15)
        self.assertAlmostEqual(elems['aft'][:,0].max(),1.)
        self.assertLess(elems['fore'][:,0].max(),.81)

    def test_published_finite_edges_are_not_public_file_sharp_edges(self):
        elements={name:np.loadtxt(OUT/f'report_test_{name}.dat',skiprows=1) for name in ('fore','aft')}
        metrics=geometry_metrics(elements)
        for name in ('fore','aft'):
            self.assertGreater(metrics[name+'_TE_thickness'],.00099)
        for file in ('s414.pts','s414-1.pts'):
            public,_=comparison_geometry(parse_pts(OUT/'sources'/file),file)
            self.assertEqual(geometry_metrics(public)['fore_TE_thickness'],0.)

    def test_distance_uses_segments_not_nearest_vertices(self):
        square=np.array([[0,0],[1,0],[1,1],[0,1],[0,0]],float)
        np.testing.assert_allclose(point_distances([[.5,-.2],[1.1,.5]],square),[.2,.1])
        self.assertAlmostEqual(symmetric_vertex_distance(square,square),0.)

    def test_pressure_affine_calibration_has_correct_sign_and_units(self):
        cal=dict(x0=20.,scale=400.,cp0=300.,cp_scale=100.)
        self.assertEqual(pressure_coordinates(220.,200.,cal),(.5,-1.))
        self.assertEqual(pressure_coordinates(420.,400.,cal),(1.,1.))

    def test_unresolved_taps_and_source_mismatches_stay_visible(self):
        geometry=json.loads((OUT/'summary.json').read_text())
        pressure=json.loads((OUT/'pressure/source.json').read_text())
        self.assertFalse(geometry['tested_geometry_benchmark_admitted'])
        self.assertEqual(geometry['tap_count'],90)
        self.assertEqual(len(geometry['unsupported_tap_centers']),3)
        self.assertEqual(pressure['experiment_markers'],188)
        self.assertEqual(pressure['repeated_endpoints'],8)
        for case in pressure['cases']:
            self.assertFalse(case['exact_table_coordinate_check_pass'])
            self.assertTrue(case['ordered_nearest_table_association_pass'])
            self.assertGreater(case['maximum_tap_x_error'],.0002)
        self.assertFalse(pressure['mses_executable_run'])
        self.assertFalse(pressure['measured_transition_endpoints_available'])

    def test_source_and_derived_hashes_match_receipts(self):
        for name in ('input_manifest.json','pressure/source.json'):
            receipt=json.loads((OUT/name).read_text())
            for path,digest in receipt['input_sha256'].items():
                self.assertEqual(sha(ROOT/path),digest,path)
            base=OUT if name=='input_manifest.json' else OUT/'pressure'
            for path,digest in receipt.get('outputs_sha256',receipt.get('output_sha256',{})).items():
                self.assertEqual(sha(base/path),digest,path)


if __name__=='__main__':unittest.main()
