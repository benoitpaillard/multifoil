"""Declared factor path, frozen references, and missing-coverage gates."""
import copy
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
import numpy as np

sys.path.insert(0,str(Path(__file__).parent/'validation'))
import run_s702_discretization as run


class DiscretizationProtocolTests(unittest.TestCase):
    def test_matrix_changes_one_declared_factor_per_step(self):
        spec = run.specification()
        self.assertEqual(len(spec['jobs']),11)
        self.assertEqual(sum(j['control'] for j in spec['jobs']),1)
        factors = ('source_panels','panels','wake_count_recipe','wake_spacing_recipe')
        for flap in (5,10):
            jobs = [next(j for j in spec['jobs'] if j['flap']==flap and j['stage']==stage) for stage in run.STAGES]
            for i,(a,b) in enumerate(zip(jobs[:-1],jobs[1:])):
                self.assertEqual([k for k in factors if a[k]!=b[k]],[factors[i]])
                for key in ('center','target_cl','wake_length','reference_result'):
                    self.assertEqual(a[key],b[key])
            for job in jobs:
                self.assertEqual(run.flow.paths_for(job,spec),run.flow.paths_for(jobs[0],spec))
        self.assertEqual((spec['reynolds'],spec['mach'],spec['ncrit']),(700000.,.07,6.25))

    def test_numerical_recipes_and_surface_source_are_independent(self):
        spec = run.specification()
        setups = {}
        for stage in run.STAGES:
            job = next(j for j in spec['jobs'] if j['name']=='d10_'+stage)
            setups[stage] = run.numerical_setup(job,spec)
        for key in ('wake_nodes','wake_first_spacing'):
            self.assertEqual(setups['fine'][1][key],setups['contour'][1][key])
            self.assertEqual(setups['fine'][1][key],setups['surface'][1][key])
        for a,b in zip(setups['contour'][0],setups['surface'][0]):
            np.testing.assert_array_equal(a[1],b[1])
        self.assertEqual(setups['surface'][1]['wake_first_spacing'],setups['wake_count'][1]['wake_first_spacing'])
        self.assertEqual(setups['wake_count'][1]['wake_nodes'],setups['wake_spacing'][1]['wake_nodes'])

    def test_modified_reference_is_rejected(self):
        spec = run.specification()
        target = run.ROOT/next(j['finer_reference_result'] for j in spec['jobs'])
        real = run.sha
        with patch.object(run,'sha',side_effect=lambda p:'changed' if Path(p)==target else real(p)):
            with self.assertRaisesRegex(ValueError,'Audited reference changed'):
                run.specification()

    def test_failed_case_cannot_pass_coverage_or_pair(self):
        spec = run.specification()
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            for job in spec['jobs']:
                folder = out/job['name']; folder.mkdir()
                r = copy.deepcopy(run.read(run.ROOT/job['reference_result']))
                # Fail the fine flap-10 case so no saved-array control can be asserted.
                if job['name']=='d10_fine':
                    r.update(feasible=False,CL=None,CD=None,CM=None,alpha_deg=None)
                run.write(folder/'result.json',r)
            summary = run.summarize(out,spec)
            self.assertFalse(summary['coverage_complete'])
            self.assertFalse(summary['execution_checks_pass'])
            self.assertFalse(summary['paired_stages'][0]['complete'])
            self.assertIsNone(summary['paired_stages'][0]['flap_delta_CD_counts'])


if __name__ == '__main__':unittest.main()
