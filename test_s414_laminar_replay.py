"""V5 tranche 1: independent laminar solver, reconstructions and the two-root station (no coupled solve)."""
from pathlib import Path
import csv
import json
import sys
import unittest
import warnings
import numpy as np
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT/'validation'))
import s414_laminar_fd as fd
import run_s414_laminar_replay as replay

STUDY = ROOT/'validation/results/s414_laminar_replay_20261004'


class IndependentSolver(unittest.TestCase):
    def test_similarity_profiles_match_tabulated_values(self):
        eta = fd.eta_grid()
        d1, d2, shear = fd.measures(eta, *fd.similarity(0., eta))
        self.assertAlmostEqual(shear, .332057, places=4); self.assertAlmostEqual(d1/d2, 2.5911, places=3)
        self.assertAlmostEqual(fd.measures(eta, *fd.similarity(1., eta))[2], 1.232588, places=3)

    def test_marching_preserves_a_similarity_flow(self):
        for m in (1., 1/3, -.05):
            r = fd.march(np.geomspace(1e-3, 1., 21), lambda x, m=m: x**m, 1e-6, m_start=m)
            self.assertEqual(r['status'], 'completed')
            self.assertLess(np.max(abs(r['wall_shear']-r['wall_shear'][0])), 1e-8)

    def test_howarth_separation_with_both_schemes(self):
        nodes = np.r_[np.geomspace(1e-4, 2e-3, 8)[:-1], np.linspace(2e-3, .13, 101)]
        for scheme in ('bdf2', 'box'):
            r = fd.march(nodes, lambda x: 1-x, 1e-6, m_start=0., scheme=scheme)
            self.assertEqual(r['status'], 'separated'); self.assertAlmostEqual(r['xi_separation'], .1198, delta=5e-4)

    def test_backward_scheme_has_no_sawtooth_at_a_kink_and_the_box_scheme_does(self):
        kinked = lambda x: 1. if x <= .5 else 1-.1*(x-.5)
        nodes = np.r_[np.linspace(1e-3, .5, 51)[:-1], np.linspace(.5, 1., 51)]
        index = {}
        for scheme in ('bdf2', 'box'):
            r = fd.march(nodes, kinked, 1e-6, m_start=0., scheme=scheme)
            self.assertEqual(r['status'], 'completed')
            index[scheme] = fd.oscillation(r['wall_shear'][r['xi'] > .5])
        self.assertLess(index['bdf2'], 1e-3); self.assertGreater(index['box'], 1e-2)

    def test_thwaites_on_a_flat_plate_and_on_the_howarth_flow(self):
        x = np.linspace(1e-6, 1., 2001)
        flat = fd.thwaites(x, np.ones_like(x), np.zeros_like(x), 1e-6)
        self.assertAlmostEqual(flat['theta'][-1]/np.sqrt(1e-6), .671, places=3)
        self.assertIsNone(flat['xi_separation'])
        x = np.linspace(1e-6, .13, 2001)
        self.assertAlmostEqual(fd.thwaites(x, 1-x, -np.ones_like(x), 1e-6)['xi_separation'], .1232, delta=1e-3)


class Reconstructions(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model, cls.param = replay.parameters(replay.MACH)
        cls.pressure = replay.Pressure(cls.param)
        cls.surfaces = {case: replay.Surface(np.load(replay.STATES/directory/'fore_state.npz'), cls.pressure)
                        for case, (directory, _) in replay.CASES.items()}
        with replay.TAPS.open() as stream:
            rows = [r for r in csv.DictReader(stream) if r["coordinate"] == "table"]
        cls.taps = {case: np.array([(float(r['x']), float(r['experiment_Cp']), float(r['viscous_Cp']))
                                    for r in rows if r['case'] == case][:replay.TAP_COUNT]) for case in replay.CASES}

    def test_pressure_relations_invert_the_inherited_ones(self):
        u = np.array([.1, .9, 1.5])
        q = self.pressure.physical(u)
        for value, speed in zip(u, q):
            self.assertAlmostEqual(replay.mf.get_uk(value, self.param)[0], speed, places=14)
            self.assertAlmostEqual(replay.mf.get_cp(value, self.param)[0], float(self.pressure.cp(speed)), places=13)
        np.testing.assert_allclose(self.pressure.speed(self.pressure.cp(q)), q, rtol=1e-11)  # cancellation near stagnation

    def test_model_tap_values_agree_with_the_archived_tap_table(self):
        for case, s in self.surfaces.items():
            own = s.cp(s.xi_of_x(self.taps[case][:, 0]))
            self.assertLess(np.max(abs(own-self.taps[case][:, 2])), .03)

    def test_every_reconstruction_passes_through_the_measured_taps(self):
        s = self.surfaces['cl0056']; data = self.taps['cl0056']; tap_xi = s.xi_of_x(data[:, 0])
        for family in replay.FAMILIES:
            for upstream in replay.UPSTREAM:
                speed = replay.reconstruction(s, tap_xi, data[:, 1], family, upstream)
                np.testing.assert_allclose(self.pressure.cp(speed(tap_xi)), data[:, 1], atol=1e-12)
                self.assertTrue(np.all(speed(s.nodes(tap_xi)) > 0))

    def test_shape_families_reduce_to_the_native_profile_for_the_model_s_own_taps(self):
        s = self.surfaces['cl0175']; tap_xi = s.xi_of_x(self.taps['cl0175'][:, 0]); nodes = s.nodes(tap_xi)
        for family in ('additive', 'ratio'):
            speed = replay.reconstruction(s, tap_xi, s.cp(tap_xi), family, 'late')
            self.assertLess(np.max(abs(speed(nodes)-s.q(nodes))), 1e-12)

    def test_one_station_has_two_roots_for_the_same_upstream_state_and_edge_velocity(self):
        s = self.surfaces['cl0056']; j = 26
        left, right = s.U[:, j-1], s.U[:, j]
        roots = []
        for guess in (right, left):
            with warnings.catch_warnings():
                warnings.simplefilter('error', RuntimeWarning)
                state, error, _ = replay.station(left, guess.copy(), right[3], s.xi[j-1:j+1], self.param)
            self.assertLess(error, 1e-9); roots.append(float(replay.mf.get_Hk(state, self.param)[0]))
        self.assertAlmostEqual(roots[0], 4.2124, places=3); self.assertAlmostEqual(roots[1], 3.6707, places=3)

    def test_closure_march_reproduces_the_archive_where_the_root_is_unique(self):
        s = self.surfaces['cl028']; end = int(np.sum(s.xi < s.xi_of_x(replay.X_END)))
        r = replay.closure_march(s.xi[:end], s.U[3, :end], self.param)
        self.assertEqual(r['status'], 'completed')
        np.testing.assert_allclose(r['states'][0, 12:], s.U[0, 12:end], rtol=1e-5)
        np.testing.assert_allclose(r['Hk'][12:], [replay.mf.get_Hk(u, self.param)[0] for u in s.U[:, 12:end].T], atol=1e-4)


class SavedStudy(unittest.TestCase):
    def test_saved_outcomes_and_limits(self):
        summary = json.loads((STUDY/'summary.json').read_text())
        self.assertEqual(summary['state'], 'finished'); self.assertTrue(summary['part_b_gate_passed'])
        self.assertEqual({k: (v['outcome'], v['admitted']) for k, v in summary['outcomes'].items()},
                         dict(cl0056=('not bounded', 0), cl0175=('not bounded', 0), cl028=('separation retained', 3)))
        self.assertEqual((summary['controls']['fd_nontrivial_passed'], summary['controls']['fd_nontrivial']), (10, 54))
        self.assertEqual((summary['counts']['fd_marches'], summary['counts']['closure_marches'],
                          summary['counts']['coupled_flow_solves']), (114, 121, 0))
        self.assertLess(summary['fd_max_sawtooth_index'], 1e-2)

    def test_first_execution_is_preserved(self):
        first = json.loads((STUDY/'initial_run_box_scheme/summary.json').read_text())
        self.assertEqual(first['counts']['fd_marches'], 114)
        manifest = json.loads((STUDY/'initial_run_box_scheme/input_manifest.json').read_text())['files']
        self.assertEqual(replay.digest(STUDY/'initial_run_box_scheme/source_snapshot/s414_laminar_fd.py'),
                         manifest['validation/s414_laminar_fd.py'])


if __name__ == '__main__':
    unittest.main()
