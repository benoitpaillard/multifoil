"""Saved-state regression gates for factorial attribution and reference labeling."""
import copy
import json
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

import numpy as np

sys.path.insert(0,str(Path(__file__).parent/'validation'))
import report_s702_surface_factorial as report


class SurfaceFactorialAuditTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.states={};cls.setups={}
        for flap in (5,10):
            for cell,study,name in (
                ('B','s702_discretization_20261001','contour'),('C','s702_discretization_20261001','surface'),
                ('F','s702_surface_factorial_20261001','fore_only'),('A','s702_surface_factorial_20261001','aft_only')):
                folder=report.ROOT/'validation/results'/study/f'd{flap}_{name}'
                state=dict(result=json.loads((folder/'result.json').read_text()),arrays={})
                for e in report.ELEMENTS:
                    with np.load(folder/f'{e}_state.npz') as saved:state['arrays'][e]={k:saved[k].copy() for k in saved.files}
                cls.states[flap,cell]=state
                cls.setups[flap,cell]=json.loads((folder/'numerical_setup.json').read_text())

    def test_no_attribution_from_missing_or_failed_cell(self):
        missing=dict(self.states);del missing[10,'A']
        with self.assertRaisesRegex(ValueError,'coverage'):report.contrasts(missing)
        for update in ({'feasible':False},{'CD':None},{'CD':float('nan')}):
            failed=copy.deepcopy(self.states);failed[10,'A']['result'].update(update)
            with self.assertRaisesRegex(ValueError,'coverage'):report.contrasts(failed)

    def test_same_path_comparison_does_not_hide_interaction_path_change(self):
        rows=report.contrasts(self.states)
        lookup={r['contrast']:r for r in rows if r['scope']=='paired' and r['component']=='total'}
        self.assertTrue(lookup['aft_at_high_fore']['same_selected_path_within_each_flap'])
        self.assertTrue(lookup['fore_at_low_aft']['same_selected_path_within_each_flap'])
        self.assertFalse(lookup['interaction']['same_selected_path_within_each_flap'])
        self.assertFalse(lookup['aft_at_low_fore']['same_selected_path_within_each_flap'])
        self.assertAlmostEqual(lookup['fore_at_low_aft']['change_counts']+lookup['aft_at_high_fore']['change_counts'],
                               lookup['total']['change_counts'],places=10)
        self.assertAlmostEqual(lookup['aft_at_high_fore']['change_counts'],-.5915064730406541,places=10)

    def test_hidden_wake_or_geometry_change_is_rejected(self):
        states={c:self.states[10,c] for c in report.CELLS};setups={c:self.setups[10,c] for c in report.CELLS}
        report.check_fixed_controls(states,setups)
        modified=copy.deepcopy(setups);modified['F']['requested']['wake_first_spacing']['aft']*=.9
        with self.assertRaisesRegex(ValueError,'Fixed wake control'):report.check_fixed_controls(states,modified)
        modified=copy.deepcopy(states);modified['F']['arrays']['aft']['xy'][1,20]+=1e-5
        with self.assertRaisesRegex(ValueError,'Fixed element surface'):report.check_fixed_controls(modified,setups)
        modified=copy.deepcopy(states);modified['C']['arrays']['fore']['source_xy'][1,20]+=1e-5
        with self.assertRaisesRegex(ValueError,'Fixed regularized source'):report.check_fixed_controls(modified,setups)

    def test_changed_comparison_table_is_rejected_before_scoring(self):
        real=report.sha
        with patch.object(report,'sha',side_effect=lambda p:'changed' if Path(p)==report.ROOT/report.POINTS else real(p)):
            with self.assertRaisesRegex(ValueError,'Archived comparison table changed'):
                report.reference_comparison(report.Inputs(),{}, {},self.states)


if __name__=='__main__':unittest.main()
