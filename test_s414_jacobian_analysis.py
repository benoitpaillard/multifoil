"""Protect paired comparisons from missing coverage, changed geometry and unit errors."""
import copy
from pathlib import Path
import sys
import unittest
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent/'validation'))
from report_s414_jacobian import paired_metrics


class PairedAuditTests(unittest.TestCase):
    def state(self):
        layers = {e: {s: dict(transition=None) for s in ('lower', 'upper')} for e in ('fore', 'aft')}
        return dict(result=dict(feasible=True, CD=.01, CL=.056, alpha_deg=-3., selected_path='cold_local',
            element_forces={e: dict(CD_wake=.005) for e in layers}, boundary_layers=layers),
            arrays={e: dict(xy=np.array([[0., 1.], [0., 0.]]), turbulent=np.array([0, 1])) for e in layers},
            surfaces={e: dict(Cp=np.array([-.2, .1])) for e in layers})

    def test_count_units_and_no_transition_are_preserved(self):
        a = self.state(); b = copy.deepcopy(a); b['result']['CD'] += .0001
        row = paired_metrics(a, b)
        self.assertAlmostEqual(row['delta_CD_counts'], 1.)
        self.assertIsNone(row['fore_lower_delta_transition_x'])
        self.assertTrue(row['fore_lower_transition_presence_matches'])

    def test_failed_candidate_cannot_enter_paired_score(self):
        a = self.state(); b = copy.deepcopy(a); b['result']['feasible'] = False
        with self.assertRaisesRegex(ValueError, 'Both modes'): paired_metrics(a, b)

    def test_unannounced_geometry_change_is_rejected(self):
        a = self.state(); b = copy.deepcopy(a); b['arrays']['fore']['xy'][1, 0] = .00001
        with self.assertRaisesRegex(ValueError, 'geometry changed'): paired_metrics(a, b)


if __name__ == '__main__': unittest.main()
