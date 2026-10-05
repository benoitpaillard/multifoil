"""Coverage, interaction accounting and actual mixed-threshold audit checks."""
from pathlib import Path
import sys
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'validation'))
import audit_s702_transition_state as state
import report_s702_wake_extent as common
from report_s702_transition_isolation import factorial,KEYS


class TransitionReportTests(unittest.TestCase):
    def test_additive_effects_and_units(self):
        f=factorial(dict(zip(KEYS,[100.,110.,96.,106.])))
        self.assertEqual((f['fore_effect'],f['aft_effect'],f['interaction'],f['combined_change']),(10.,-4.,0.,6.))

    def test_interaction_changes_both_conditional_effects(self):
        f=factorial(dict(zip(KEYS,[100.,110.,96.,109.])))
        self.assertEqual(f['interaction'],3.)
        self.assertEqual(f['conditional_fore_effect'],13.)
        self.assertEqual(f['conditional_aft_effect'],-1.)
        self.assertEqual(f['fore_effect']+f['aft_effect']+f['interaction'],f['combined_change'])

    def test_missing_endpoint_cannot_receive_complete_score(self):
        with self.assertRaisesRegex(ValueError,'complete four-state'):factorial(dict(zip(KEYS[:3],[1.,2.,3.])))

    def test_nonfinite_endpoint_cannot_receive_score(self):
        with self.assertRaisesRegex(ValueError,'Nonfinite'):factorial(dict(zip(KEYS,[1.,2.,3.,float('nan')])))

    def test_common_and_mixed_threshold_maps(self):
        self.assertEqual(state.thresholds(dict(ncrit=6.25)),dict(fore=6.25,aft=6.25))
        expected=dict(fore=9.,aft=6.25)
        self.assertEqual(state.thresholds(dict(ncrit_by_element=expected)),expected)
        state.check_threshold_provenance(dict(ncrit=9.,ncrit_by_element=expected),expected)
        with self.assertRaises(ValueError):state.check_threshold_provenance(dict(ncrit=6.25,ncrit_by_element=expected),expected)

    def test_complete_residual_receipt_retains_actual_threshold_maps(self):
        inputs=state.Inputs();receipt=inputs.json(state.OUT/'residual_recheck.json')
        self.assertTrue(receipt['complete']);self.assertEqual(len(receipt['states']),18)
        jobs=inputs.json(state.STUDY/'study.json')['jobs']+inputs.json(state.PARENT/'study.json')['jobs']
        for job in jobs:
            row=next(r for r in receipt['states'] if r['job']==job['name'])
            self.assertEqual(row['ncrit_by_element'],state.thresholds(job))

    def test_wrong_aft_provenance_rejected_on_actual_mixed_state(self):
        spec=state.Inputs().json(state.STUDY/'study.json');runtime=state.Inputs().json(state.STUDY/'runtime.json')
        job=next(j for j in spec['jobs'] if j['name']=='d5_a2p03_aft_only')
        class CorruptInputs(state.Inputs):
            def json(self,path):
                value=super().json(path)
                if Path(path).name=='metadata.json':value['diagnostics'][0]['ncrit_by_element']['aft']=9.
                return value
        with self.assertRaisesRegex(ValueError,'Per-element Ncrit differs'):
            state.audit(state.STUDY/job['name'],job,spec,runtime,CorruptInputs())

    def test_original_common_threshold_audits_are_preserved(self):
        inputs=state.Inputs();spec=inputs.json(state.PARENT/'study.json');runtime=inputs.json(state.PARENT/'runtime.json')
        for j in spec['jobs']:
            settings=dict(spec,ncrit=j['ncrit'])
            original=common.audit(state.PARENT/j['name'],j,settings,runtime,state.Inputs())
            adapted=state.audit(state.PARENT/j['name'],j,settings,runtime,state.Inputs())
            self.assertEqual(original['result'],adapted['result']);self.assertEqual(original['residual'],adapted['residual'])
            for e in ('fore','aft'):
                for field in original['arrays'][e]:
                    np.testing.assert_array_equal(original['arrays'][e][field],adapted['arrays'][e][field])


if __name__=='__main__':unittest.main()
