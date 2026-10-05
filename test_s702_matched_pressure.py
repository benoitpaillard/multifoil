"""Independent data checks and physical-condition/seed safeguards."""
import csv
import json
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'validation'))
import extract_s702_flap_pressure as extract
import run_s702_matched_pressure as run


class MatchedPressureTests(unittest.TestCase):
    def test_independent_fore_and_prehinge_taps(self):
        for c in extract.CASES:
            rows=run.csv_rows(extract.OUT/c['name']/'fore_upper.csv')
            np.testing.assert_allclose([float(r['x_over_cref']) for r in rows],extract.FORE_UPPER,atol=8e-5,rtol=0)
            aft=run.csv_rows(extract.OUT/c['name']/'aft_upper.csv')
            x=np.array([float(r['x_over_cref']) for r in aft if float(r['x_over_cref'])<.89])
            self.assertLess(float(abs(x[:,None]-extract.AFT_UPPER_PREHINGE).min(axis=1).max()),8e-5)

    def test_duplicated_endpoints_not_lower_taps(self):
        for c in extract.CASES:
            for e in ('fore','aft'):
                up=run.csv_rows(extract.OUT/c['name']/(e+'_upper.csv'))
                lo=run.csv_rows(extract.OUT/c['name']/(e+'_lower.csv'))
                self.assertEqual(sum(int(r['endpoint_duplicate']) for r in lo),2)
                for i in (0,-1):
                    self.assertEqual(lo[i]['endpoint_duplicate'],'1')
                    self.assertLess(abs(float(lo[i]['Cp'])-float(up[i]['Cp'])),.001)

    def test_nearest_lift_guard_and_source_association(self):
        source=run.read(extract.OUT/'source.json')
        for c in source['cases']:
            rows=run.csv_rows(extract.DATA/f'd+{c["flap"]}'/'experiment_lift.csv')
            self.assertEqual(extract.match_lift(rows,c['alpha'])['target_cl'],c['target_cl'])
        with self.assertRaises(ValueError):
            extract.match_lift([dict(alpha_deg=1.,CL=.5,drawing_index=1)],1.1)

    def test_job_conditions_and_frozen_sources(self):
        s=run.specification()
        self.assertEqual(len(s['jobs']),8)
        self.assertEqual((s['reynolds'],s['mach'],s['panels']),(1e6,.1,[340,240]))
        for c in extract.CASES:
            jobs=[j for j in s['jobs'] if j['pressure_case']==c['name']]
            self.assertEqual({j['ncrit'] for j in jobs},{9.,6.25})
            self.assertEqual(len({j['target_cl'] for j in jobs}),1)
            self.assertTrue(all(-6<j['center']<8 for j in jobs))

    def test_seed_extrapolation_guard(self):
        rows=run.csv_rows(ROOT/run.SEEDS)
        with self.assertRaises(ValueError):run.seed(rows,5,9.,3.)

    def test_resume_rejects_changed_protocol(self):
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp);spec=dict(sha256={},panels=[340,240])
            run.prepare(out,spec)
            with self.assertRaises(ValueError):run.prepare(out,dict(spec,panels=[440,320]))


if __name__=='__main__':unittest.main()
