"""Coverage and tamper regressions for the saved-state discretization audit."""
import copy
import json
from pathlib import Path
import shutil
import sys
import tempfile
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).parent/'validation'))
import report_s702_discretization as report


class DiscretizationAuditTests(unittest.TestCase):
    def setUp(self):
        self.study = report.ROOT/'validation/results/s702_discretization_20261001'
        self.spec = json.loads((self.study/'study.json').read_text())
        self.results = {(j['flap'], j['stage']):json.loads((self.study/j['name']/'result.json').read_text())
                        for j in self.spec['jobs'] if not j['control']}

    def test_missing_or_failed_pair_cannot_produce_attribution(self):
        absent = dict(self.results); del absent[10, 'surface']
        with self.assertRaisesRegex(ValueError, 'coverage'):
            report.paired_increments(absent)
        for mutation in ({'feasible':False}, {'CD':None}, {'CD':float('nan')}):
            failed = copy.deepcopy(self.results); failed[10, 'surface'].update(mutation)
            with self.assertRaisesRegex(ValueError, 'coverage'):
                report.paired_increments(failed)

    def test_own_baseline_and_element_changes_telescope(self):
        pairs, increments, total = report.paired_increments(self.results)
        self.assertAlmostEqual(total, -.754786244392871, places=10)
        self.assertAlmostEqual(increments[1]['paired_change_counts'], -.5612800411452548, places=10)
        self.assertAlmostEqual(sum(r['aft_paired_change_counts']+r['fore_paired_change_counts'] for r in increments), total, places=10)
        shifted = copy.deepcopy(self.results)
        # A common offset to both flap drags cannot change any paired result.
        for r in shifted.values():
            r['CD'] += .003; r['element_forces']['fore']['CD_wake'] += .003
        pairs2, _, total2 = report.paired_increments(shifted)
        np.testing.assert_allclose([v['flap_delta_CD_counts'] for v in pairs2], [v['flap_delta_CD_counts'] for v in pairs], atol=1e-12, rtol=0)
        self.assertAlmostEqual(total2, total, places=10)

    def test_hidden_second_control_change_is_rejected(self):
        setups = {stage:json.loads((self.study/('d10_'+stage)/'numerical_setup.json').read_text()) for stage in report.STAGES}
        report.check_factor_controls(setups)
        for field in ('source_sha256', 'wake_first_spacing'):
            bad = copy.deepcopy(setups)
            if field == 'source_sha256': bad['surface'][field]['aft'] = 'altered'
            else: bad['surface']['requested'][field]['aft'] *= 1.01
            with self.assertRaisesRegex(ValueError, 'Unexpected factor change'):
                report.check_factor_controls(bad)

    def test_missing_control_array_cannot_be_hidden_by_intersection(self):
        arrays = {e:dict(U=np.ones((4,3)), xy=np.ones((2,3))) for e in report.ELEMENTS}
        missing = copy.deepcopy(arrays); del missing['aft']['U']
        with self.assertRaisesRegex(ValueError, 'Missing reference array'):
            report.compare_arrays(missing, arrays, 'missing state')

    def test_saved_source_hash_and_geometry_tampering_are_detected(self):
        job = next(j for j in self.spec['jobs'] if j['name'] == 'd10_contour')
        runtime = json.loads((self.study/'runtime.json').read_text())
        expected = dict(self.spec, panels=job['panels'], source_panels=job['source_panels'])
        state = report.audit(self.study/job['name'], job, expected, runtime, report.Inputs())
        with tempfile.TemporaryDirectory(dir=self.study.parent, prefix='discretization_audit_test_') as tmp:
            folder = Path(tmp)/job['name']; shutil.copytree(self.study/job['name'], folder)
            report.audit_setup(folder, job, self.spec, state, report.Inputs())
            path = folder/'source_aft.npy'; original = path.read_bytes()
            xy = np.load(path); xy[20,1] += 1e-4; np.save(path, xy)
            with self.assertRaisesRegex(ValueError, 'Source array hash differs'):
                report.audit_setup(folder, job, self.spec, state, report.Inputs())
            path.write_bytes(original)
            bad = copy.deepcopy(state); bad['arrays']['aft']['source_xy'][1,20] += 1e-4
            with self.assertRaisesRegex(ValueError, 'Regularized source differs'):
                report.audit_setup(folder, job, self.spec, bad, report.Inputs())


if __name__ == '__main__':
    unittest.main()
