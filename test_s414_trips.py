"""Source, branch-selection and matched-lift contracts for S414 trip data."""
import csv
import json
from pathlib import Path
import sys
import unittest
import numpy as np
import fitz

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'validation'))
import qualify_s414_trips as q


def read_rows(name):
    with (q.OUT/name).open() as stream:return list(csv.DictReader(stream))


class S414TripTests(unittest.TestCase):
    def test_missing_drag_and_poststall_rows_are_preserved(self):
        text='R = 0.70 × 106, M = 0.07, transition free\n-1.0 -.1 .01 -.1\n0.0 .5 .02 -.1\n1.0 .3 — -.1'
        rows=q.parse_table(text,q.CONDITIONS[0],'free',88)
        self.assertEqual(len(rows),3)
        self.assertIsNone(rows[-1]['CD']);self.assertFalse(rows[-1]['drag_available'])
        self.assertEqual([r['prestall_branch'] for r in rows],[True,True,False])
        with self.assertRaises(ValueError):q.parse_table(text,q.CONDITIONS[0],'fore_fixed',88)

    def test_extrapolation_or_reordering_folded_polars_forbidden(self):
        rows=[dict(CL=0.,CD=.01,alpha_deg=0.),dict(CL=1.,CD=.02,alpha_deg=2.),dict(CL=.8,CD=.1,alpha_deg=3.)]
        with self.assertRaises(ValueError):q.interpolate(rows,.9)
        with self.assertRaises(ValueError):q.interpolate(q.ascending_branch(rows),1.1)
        self.assertAlmostEqual(q.interpolate(q.ascending_branch(rows),.5),.015)
        self.assertEqual(len(rows),3)

    def test_published_table_values_from_independent_pages(self):
        with fitz.open(q.PDF) as doc:
            free=q.parse_table(doc[90].get_text(sort=True),q.CONDITIONS[1],'free',91)
            fore=q.parse_table(doc[91].get_text(sort=True),q.CONDITIONS[1],'fore_fixed',92)
        self.assertEqual((free[7]['CL'],free[7]['CD']),(.2837,.006549))
        self.assertEqual((fore[3]['CL'],fore[3]['CD']),(.2726,.010300))
        self.assertEqual((free[16]['CL'],free[16]['CD']),(1.0410,.014489))
        self.assertEqual(sum(r['prestall_branch'] for r in free),29)
        self.assertFalse(fore[-1]['prestall_branch'])

    def test_target_delta_uses_each_own_lift_grid(self):
        free=.006549+(.284-.2837)/(.3903-.2837)*(.006745-.006549)
        fore=.010300+(.284-.2726)/(.3746-.2726)*(.010393-.010300)
        row=next(r for r in read_rows('matched_lift_targets.csv') if r['reynolds']=='1000000' and r['CL']=='0.284')
        self.assertAlmostEqual(float(row['experiment_fore_minus_free_counts']),1e4*(fore-free),places=10)
        # Comparing same acquisition index would omit the measured lift change.
        self.assertNotAlmostEqual(float(row['experiment_fore_minus_free_counts']),1e4*(.010300-.006549),places=3)

    def test_all_source_rows_and_missing_values_retained(self):
        rows=read_rows('experimental_tables.csv');self.assertEqual(len(rows),169)
        missing=[r for r in rows if r['CD']=='']
        self.assertEqual(len(missing),1)
        self.assertEqual(missing[0]['prestall_branch'],'False')
        self.assertEqual(missing[0]['drag_available'],'False')
        for c in q.CONDITIONS:
            for state,count in zip(q.STATES,c['table_counts']):
                self.assertEqual(len([r for r in rows if r['reynolds']==str(c['reynolds']) and r['configuration']==state]),count)

    def test_nominal_trip_labels_do_not_pass_causal_admission(self):
        e=json.loads((q.OUT/'evidence_qualification.json').read_text())
        self.assertTrue(e['admission']['matched_global_force_context'])
        self.assertFalse(e['admission']['clean_wake_disturbance_intervention'])
        self.assertFalse(e['admission']['matched_tripped_aft_pressure_benchmark'])
        self.assertEqual(e['trips'][0]['x_over_total_chord'],.02)
        self.assertEqual(e['trips'][2]['x_over_aft_chord'],.05)
        self.assertNotIn('x_over_total_chord',e['trips'][2])
        self.assertIn('1.50',f"{e['excluded_conditions'][1]['reynolds']/1e6:.2f}")

    def test_strict_copy_failure_and_original_offsets_are_retained(self):
        r=json.loads((q.OUT/'summary.json').read_text())
        self.assertFalse(r['strict_repeated_pressure_copy_check_pass'])
        self.assertGreater(r['maximum_repeated_pressure_Cp_difference'],3e-4)
        self.assertLess(r['maximum_repeated_pressure_Cp_difference'],4e-4)
        pressure=read_rows('free_boundary_pressure.csv')
        self.assertEqual(len(pressure),376)
        self.assertEqual(sum(r['endpoint_duplicate']=='True' for r in pressure),16)
        for case in q.PRESSURES:
            rr=[r for r in pressure if r['figure']==case['figure'] and r['surface']=='aft_upper']
            self.assertEqual(len(rr),20)
            self.assertTrue(any(float(r['plotted_x'])!=float(r['table_x']) for r in rr))

    def test_pressure_sequence_does_not_imply_monotone_bubble_shift(self):
        intervals=read_rows('pressure_recovery_intervals.csv')
        peaks=[]
        for case in q.PRESSURES:
            rr=[r for r in intervals if r['figure']==case['figure'] and r['coordinate']=='table_x']
            peak=max(rr,key=lambda r:float(r['delta_Cp']))
            peaks.append(float(peak['left_x']))
        np.testing.assert_array_equal(peaks,[.88215,.88215,.89313,.88215])


if __name__=='__main__':unittest.main()
