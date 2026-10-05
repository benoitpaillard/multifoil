"""Scientific state/provenance checks for the archived drag-mechanism report."""
import copy
import csv
from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT/'validation'))
from report_s702_drag_mechanism import (
    Inputs, load_state, selected_index, pair_deltas, check_published,
    RESULTS, DATA, load_study, from_analysis,
)


class StateSelectionTests(unittest.TestCase):
    def setUp(self):
        self.result = dict(CL=.8, CD=.01, CM=-.2, alpha_deg=1., selected_flow_index=0)
        self.polar = [dict(CL=.8, CD=.01, CM=-.2, alpha_deg=1., converged=1),
                      dict(CL=1.2, CD=.02, CM=-.3, alpha_deg=5., converged=1)]

    def test_recorded_state_is_selected_even_when_export_has_later_solutions(self):
        self.assertEqual(selected_index(self.result, self.polar), 0)
        self.result['selected_flow_index'] = 1
        with self.assertRaisesRegex(ValueError, 'index disagrees'):
            selected_index(self.result, self.polar)

    def test_matching_forces_do_not_override_nonconvergence(self):
        self.polar[0]['converged'] = 0
        with self.assertRaisesRegex(ValueError, 'did not converge'):
            selected_index(self.result, self.polar)

    def test_legacy_export_with_duplicate_force_states_is_ambiguous(self):
        del self.result['selected_flow_index']
        self.polar.append(self.polar[0].copy())
        with self.assertRaisesRegex(ValueError, 'ambiguous'):
            selected_index(self.result, self.polar)

    def test_nonfinite_selected_drag_is_rejected(self):
        self.result['CD'] = float('nan')
        with self.assertRaises(ValueError):
            selected_index(self.result, self.polar)

    def test_changed_archived_bytes_are_rejected(self):
        inputs = Inputs()
        path = 'validation/results/s702_flap_effect/s702_flap_5_to_10.csv'
        inputs.original[path]['sha256'] = '0'*64
        with self.assertRaisesRegex(ValueError, 'Archived input changed'):
            inputs.add(ROOT/path)


class ArchivedMechanismTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.inputs = Inputs()
        cls.refs = from_analysis(RESULTS/'s702_ncrit_analysis')
        cls.states, specs = [], {}
        for point in cls.refs:
            if not point['series'].startswith('multifoil'):
                continue
            folder = (ROOT/point['result_path']).parents[2]
            if folder not in specs:
                specs[folder], _, _ = load_study(folder,
                    DATA/'re700k_holdout' if point['dataset'] == 'reserved' else DATA)
            cls.states.append(load_state(point, specs[folder], cls.inputs))
        cls.reproduced, cls.deltas, _, _ = pair_deltas(cls.states, cls.refs)

    def test_complete_coverage_selected_states_and_actual_exports(self):
        states = [s['summary'] for s in self.states]
        self.assertEqual(len(states), 42)
        self.assertEqual(len({s['state_id'] for s in states}), 42)
        self.assertEqual(sum(s['primary'] for s in states), 16)
        self.assertEqual(sum(s['reused'] for s in states), 12)
        self.assertLessEqual(max(abs(s['lift_error']) for s in states), 2e-5)
        self.assertLessEqual(max(s['residual'] for s in states), 1e-7)

    def test_legacy_alias_and_two_retry_selections_preserve_original_failures(self):
        states = [s['summary'] for s in self.states]
        reused = [s for s in states if s['reused']]
        self.assertTrue(all(s['flow_directory'].endswith(f"cl{s['target_cl']:.4f}") for s in reused))
        retried = [s for s in states if s['numerical_retry']]
        self.assertEqual({(s['dataset'], s['ncrit'], s['flap'], s['target_cl']) for s in retried},
                         {('reserved', 6.25, 10, .5), ('reserved', 6.25, 10, .8)})
        for s in retried:
            original = self.inputs.json(ROOT/s['retry_of'])
            self.assertFalse(original['feasible'])
            self.assertIsNone(original['CD'])
        failed_cold = self.inputs.json(RESULTS/'s702_ncrit_holdout_local_cl08/n6p25_d+10/cl0p8/result.json')
        self.assertFalse(failed_cold['feasible'])
        self.assertIsNone(failed_cold['CD'])

    def test_all_six_reference_deltas_reproduced(self):
        original = self.inputs.csv(RESULTS/'s702_flap_effect/s702_flap_5_to_10.csv')
        check_published(self.reproduced, original)
        with self.assertRaisesRegex(ValueError, 'Missing published'):
            check_published(self.reproduced[:-1], original)
        changed = copy.deepcopy(self.reproduced)
        changed[0]['multifoil_n9_delta_10_minus_5_counts'] += .01
        with self.assertRaisesRegex(ValueError, 'published'):
            check_published(changed, original)

    def test_friction_is_never_added_to_wake_drag(self):
        for state in self.states:
            s = state['summary']
            self.assertAlmostEqual(s['CD'], s['fore_CD_wake']+s['aft_CD_wake'], places=12)
            self.assertGreater(s['fore_CD_friction']+s['aft_CD_friction'], .001)
        for d in self.deltas:
            self.assertAlmostEqual(d['delta_CD_counts'], d['fore_delta_wake_counts']+d['aft_delta_wake_counts'], places=8)

    def test_exported_geometry_is_identical_across_reynolds_lift_and_ncrit(self):
        by_flap = {}
        fore = None
        for state in self.states:
            s = state['summary']
            for element in ('fore', 'aft'):
                g = np.loadtxt(ROOT/s['flow_directory']/f'geometry_{element}.dat', skiprows=1)
                key = (s['flap'], element)
                if key in by_flap:
                    np.testing.assert_array_equal(g, by_flap[key])
                else:
                    by_flap[key] = g
                if element == 'fore':
                    if fore is not None:
                        np.testing.assert_array_equal(g, fore)
                    fore = g

    def test_generated_tables_match_authoritative_selected_states(self):
        with (RESULTS/'s702_drag_mechanism/selected_states.csv').open(newline='') as stream:
            exported = {r['state_id']: r for r in csv.DictReader(stream)}
        self.assertEqual(set(exported), {s['summary']['state_id'] for s in self.states})
        for state in self.states:
            s = state['summary']
            for field in ('CL', 'CD', 'CM', 'alpha_deg', 'residual', 'fore_CD_wake', 'aft_CD_wake'):
                self.assertEqual(float(exported[s['state_id']][field]), s[field])


if __name__ == '__main__':
    unittest.main()
