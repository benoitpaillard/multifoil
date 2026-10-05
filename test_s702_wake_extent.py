"""Wake-study coverage, source integrity and interpretation safeguards; no CFD."""
import copy
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).parent/'validation'))
import run_s702_wake_extent as run


class WakeExtentTests(unittest.TestCase):
    def setUp(self):
        self.spec = run.specification()

    def populate(self, folder):
        for job in self.spec['jobs']:
            target = folder/job['name']; target.mkdir()
            run.write(target/'result.json', run.read(run.ROOT/job['reference_result']))

    def test_matched_matrix_retains_one_baseline_control_and_common_policy(self):
        primary = [j for j in self.spec['jobs'] if not j['control']]
        controls = [j for j in self.spec['jobs'] if j['control']]
        self.assertEqual({(j['flap'], j['target_cl']) for j in primary}, {(f, c) for f in (5, 10) for c in (.8, 1.2)})
        self.assertEqual(len(primary), 4)
        self.assertEqual([j['wake_length'] for j in primary], [2.]*4)
        self.assertEqual([(j['flap'], j['target_cl'], j['wake_length']) for j in controls], [(10, .8, 1.)])
        old = run.read(run.ROOT/run.OLD/'study.json')
        for key in ('mach', 'reynolds', 'ncrit', 'panels', 'reference_chord', 'cl_tolerance', 'residual_tolerance'):
            self.assertEqual(self.spec[key], old[key])
        for job in self.spec['jobs']:
            paths = run.paths_for(job, self.spec)
            self.assertEqual([p.name for p in paths], ['cold_local', 'warm_anchor', 'warm_entry'])
            self.assertEqual(paths[-1].seed_alphas, (job['center']-.4,))

    def test_audited_reference_change_is_rejected(self):
        path = self.spec['jobs'][0]['reference_result']
        real_sha = run.sha
        with patch.object(run, 'sha', side_effect=lambda p: 'changed' if Path(p) == run.ROOT/path else real_sha(p)):
            with self.assertRaisesRegex(ValueError, 'Audited reference changed'):
                run.specification()

    def test_failed_case_prevents_complete_comparison_and_keeps_null_delta(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp); self.populate(out)
            self.assertTrue(run.summarize(out, self.spec)['execution_checks_pass'])
            path = out/self.spec['jobs'][0]['name']/'result.json'
            result = run.read(path)
            result.update(feasible=False, CL=None, CD=None, CM=None, alpha_deg=None)
            run.write(path, result)
            summary = run.summarize(out, self.spec)
            self.assertFalse(summary['execution_checks_pass'])
            self.assertEqual(summary['feasible_jobs'], 4)
            self.assertIsNone(summary['jobs'][0]['difference_from_1c_CD_counts'])
            result['CD'] = .01
            run.write(path, result)
            with self.assertRaisesRegex(ValueError, 'Failed target'):
                run.summarize(out, self.spec)

    def test_large_extent_response_is_observable_but_control_disagreement_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp); self.populate(out)
            def change(job):
                path = out/job['name']/'result.json'
                r = copy.deepcopy(run.read(path))
                r['CD'] += .002
                r['element_forces']['aft']['CD_wake'] += .002
                r['records'][r['selected_flow_index']]['CD'] += .002
                run.write(path, r)
            change(self.spec['jobs'][0])
            summary = run.summarize(out, self.spec)
            self.assertTrue(summary['execution_checks_pass'])
            self.assertAlmostEqual(summary['jobs'][0]['difference_from_1c_CD_counts'], 20.)
            change(next(j for j in self.spec['jobs'] if j['control']))
            summary = run.summarize(out, self.spec)
            self.assertTrue(summary['coverage_complete'])
            self.assertFalse(summary['control_pass'])
            self.assertFalse(summary['execution_checks_pass'])


if __name__ == '__main__':
    unittest.main()
