"""The declared tranche-2 selection rule of the S414 transition-boundary study (synthetic inputs)."""
from pathlib import Path
import sys
import unittest
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT/'validation'))
import report_s414_transition_boundary as report

ANGLES = [-2.6, -2.575, -2.55, -2.525, -2.5]


def point(x, drag=.007):
    layer = dict(transition=None if x is None else dict(xy_over_cref=[x, 0.]))
    return dict(record=dict(CD=drag, boundary_layers=dict(fore=dict(lower=layer))))


def paths(increasing, decreasing):
    make = lambda values: {a: point(x) for a, x in zip(ANGLES, values) if x != 'failed'}
    return dict(increasing=make(increasing), decreasing=make(decreasing))


def pairs(p, discordant=()):
    out = {}
    for a in ANGLES:
        both = a in p['increasing'] and a in p['decreasing']
        out[a] = dict(status='both_converged' if both else 'only_one', within_screens=both and a not in discordant)
    return out


class SelectionRule(unittest.TestCase):
    def test_late_transition_everywhere_selects_nothing_and_names_the_fallback(self):
        p = paths([.75, .76, .77, .78, .79], [.75, .76, .77, .78, .79])
        s = report.select_refinement(ANGLES, p, pairs(p))
        self.assertEqual((s['new_locations'], s['cold_check'], s['ranked_eligible']), ([], None, []))
        self.assertIn('archived CL 0.056 to 0.175 history', s['fallback'])

    def test_path_difference_outranks_transition_change_and_fills_four_locations(self):
        p = paths([.08, .39, .66, .74, .77], [.08, .39, 'failed', .74, .77])
        s = report.select_refinement(ANGLES, p, pairs(p))
        self.assertEqual(s['ranked_eligible'], [[ANGLES[1], ANGLES[2]], [ANGLES[2], ANGLES[3]], [ANGLES[0], ANGLES[1]]])
        self.assertEqual([round(l['alpha_deg'], 6) for l in s['new_locations']], [-2.5625, -2.5375, -2.5875, -2.56875])
        self.assertEqual([l['rule'] for l in s['new_locations']], ['midpoint']*3+['quarter point of top-ranked interval'])
        self.assertEqual(s['cold_check']['alpha_deg'], ANGLES[2])
        self.assertTrue(all(i['transition_change_unbounded'] for i in s['intervals'][1:3]))
        self.assertFalse(s['intervals'][3]['eligible'])

    def test_transition_change_alone_ranks_by_size_then_lower_angle(self):
        p = paths([.10, .40, .50, .75, .77], [.10, .40, .50, .75, .77])
        s = report.select_refinement(ANGLES, p, pairs(p))
        self.assertEqual(s['ranked_eligible'], [[ANGLES[0], ANGLES[1]], [ANGLES[2], ANGLES[3]], [ANGLES[1], ANGLES[2]]])
        self.assertEqual(len(s['new_locations']), 4)
        self.assertEqual(s['cold_check']['rule'], 'midpoint of top-ranked interval')

    def test_drag_never_enters_the_selection(self):
        a = paths([.10, .40, .50, .75, .77], [.10, .40, .50, .75, .77])
        b = paths([.10, .40, .50, .75, .77], [.10, .40, .50, .75, .77])
        for direction in b.values():
            for k, value in enumerate(direction.values()): value['record']['CD'] = .02-.003*k
        self.assertEqual(report.select_refinement(ANGLES, a, pairs(a)), report.select_refinement(ANGLES, b, pairs(b)))

    def test_screen_violation_at_a_converged_angle_is_a_path_difference(self):
        p = paths([.75, .76, .77, .78, .79], [.75, .76, .77, .78, .79])
        s = report.select_refinement(ANGLES, p, pairs(p, discordant=(ANGLES[3],)))
        self.assertEqual(s['ranked_eligible'], [[ANGLES[2], ANGLES[3]], [ANGLES[3], ANGLES[4]]])
        self.assertEqual(s['cold_check'], dict(alpha_deg=ANGLES[3], rule='angle with a path difference'))


if __name__ == '__main__':
    unittest.main()
