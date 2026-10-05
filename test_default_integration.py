"""Integration of the corrected amplification Jacobian as the CoupledSolver default.

Covers dispatch, propagation, provenance and isolated branch derivatives. The
archived-state replay is validation/verify_s414_default_integration.py. The
small NACA 0012 solves here are integration checks, not aerodynamic evidence.
"""
import copy
import csv
import inspect
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np

ROOT = Path(__file__).resolve().parent
import amplification_jacobian as option
import multifoil
import shape_optimization
import validate
import viscous
from viscous import CoupledSolver, mf

STATE = ROOT/'validation/results/s414_lower_bucket_20261003/cl0056_nominal/fore_state.npz'


def thermo_param():
    model = SimpleNamespace(geom=mf.Geom(), oper=mf.Oper(), param=mf.Param())
    model.geom.chord = 1.; model.oper.Re = 1e6; model.oper.Ma = .1; model.param.ncrit = 9
    mf.init_thermo(model)
    return model.param


def regime(U, param):
    """Branch of the vendored growth law at one station; mirrors its arithmetic."""
    hk = mf.get_Hk(U, param)[0]
    clipped = bool(hk < 1.05)
    hk = max(hk, 1.05); hmi = 1/(hk-1)
    critical = 2.492*hmi**.43+.7*(np.tanh(14*hmi-9.24)+1)
    rn = (np.log(mf.get_Ret(U, param)[0])/np.log(10)-(critical-.1))/.2
    return ('inactive' if rn < 0 else 'ramp' if rn < 1 else 'full'), clipped, float(rn)


def state_for(param, hk, rn, amplification=2., ue=.8):
    """Station with a requested kinematic shape factor and activation coordinate."""
    hmi = 1/(max(hk, 1.05)-1)
    critical = 2.492*hmi**.43+.7*(np.tanh(14*hmi-9.24)+1)
    lo, hi = 1.0001, 40.
    for _ in range(200):  # get_Hk is monotone in delta*/theta at fixed edge speed
        mid = .5*(lo+hi)
        if mf.get_Hk(np.array([1e-3, mid*1e-3, amplification, ue]), param)[0] < hk: lo = mid
        else: hi = mid
    shape = .5*(lo+hi)
    unit = mf.get_Ret(np.array([1., shape, amplification, ue]), param)[0]
    theta = 10**(critical-.1+.2*rn)/unit
    return np.array([theta, shape*theta, amplification, ue])


class DefaultAndDispatch(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        xy = multifoil.naca4('0012', 81)
        cls.standard = CoupledSolver([('main', xy)])
        cls.legacy = CoupledSolver([('main', xy)], amplification_jacobian='legacy')
        cls.adapter = option.AmplificationJacobianSolver([('main', xy)])
        cls.adapter_legacy = option.AmplificationJacobianSolver([('main', xy)], amplification_jacobian=option.LEGACY)
        with np.load(STATE) as saved:
            ids = saved['lower']; lam = ids[saved['turbulent'][ids] == 0]
            cls.U = saved['U'][:, lam[-1]].copy()

    def test_declared_defaults_and_single_implementation(self):
        default = inspect.signature(CoupledSolver.__init__).parameters['amplification_jacobian'].default
        self.assertEqual(default, 'corrected_af')
        self.assertEqual(inspect.signature(option.AmplificationJacobianSolver.__init__)
                         .parameters['amplification_jacobian'].default, option.CORRECTED)
        self.assertTrue(issubclass(option.AmplificationJacobianSolver, CoupledSolver))
        # The adapter must not carry a second copy of the solve or provenance logic.
        self.assertEqual(set(vars(option.AmplificationJacobianSolver)) - {'__module__', '__doc__', '__init__',
                         '__firstlineno__', '__static_attributes__'}, set())
        self.assertEqual(option.MODES, ('legacy', 'corrected_af'))

    def test_models_copies_and_side_parameters_carry_mode(self):
        for solver, mode in ((self.standard, 'corrected_af'), (self.legacy, 'legacy'),
                             (self.adapter, 'corrected_af'), (self.adapter_legacy, 'legacy')):
            self.assertEqual(solver.amplification_jacobian, mode)
            for model in list(solver.models)+copy.deepcopy(solver.models):
                self.assertEqual(model.param.amplification_jacobian, mode)
                for side in (0, 1, 2):
                    self.assertEqual(mf.build_param(model, side).amplification_jacobian, mode)

    def test_corrected_and_legacy_coexist_in_one_process(self):
        expected = option._original_get_damp(self.U, thermo_param())
        self.assertIs(mf.get_damp, option.get_damp)
        for solver in (self.legacy, self.adapter_legacy):
            actual = mf.get_damp(self.U, mf.build_param(solver.models[0], 0))
            self.assertEqual(actual[0], option._original_get_damp(self.U, mf.build_param(solver.models[0], 0))[0])
            np.testing.assert_array_equal(actual[1], option._original_get_damp(self.U, mf.build_param(solver.models[0], 0))[1])
        untagged = mf.get_damp(self.U, thermo_param())
        self.assertEqual(untagged[0], expected[0]); np.testing.assert_array_equal(untagged[1], expected[1])
        a = mf.get_damp(self.U, mf.build_param(self.standard.models[0], 0))
        b = mf.get_damp(self.U, mf.build_param(self.adapter.models[0], 0))
        c = mf.get_damp(self.U, mf.build_param(self.legacy.models[0], 0))
        self.assertEqual(a[0], c[0]); np.testing.assert_array_equal(a[1], b[1])
        self.assertNotEqual(a[1][1], c[1][1])
        # Constructing another legacy solver later must not downgrade existing corrected models.
        CoupledSolver([('main', multifoil.naca4('0012', 81))], amplification_jacobian='legacy')
        np.testing.assert_array_equal(mf.get_damp(self.U, mf.build_param(self.standard.models[0], 0))[1], a[1])

    def test_invalid_mode_is_rejected_at_every_entry(self):
        for build in (lambda: CoupledSolver([], amplification_jacobian='typo'),
                      lambda: option.AmplificationJacobianSolver([], amplification_jacobian='typo'),
                      lambda: option.configure([], 'typo'),
                      lambda: shape_optimization.check_config({'amplification_jacobian': 'typo'})):
            with self.assertRaisesRegex(ValueError, 'Unknown amplification Jacobian'): build()
        bad = SimpleNamespace(param=SimpleNamespace(amplification_jacobian='typo'))
        with self.assertRaisesRegex(ValueError, 'Unknown amplification Jacobian'):
            option.get_damp(self.U, bad.param)

    def test_both_import_orders_construct_without_circular_failure(self):
        code = ("import {first}\nimport {second}\nfrom multifoil import naca4\n"
                "import viscous, amplification_jacobian as o\n"
                "s={entry}([('m',naca4('0012',61))])\n"
                "assert s.amplification_jacobian=='corrected_af' and viscous.mf.get_damp is o.get_damp\n"
                "assert s.models[0].param.amplification_jacobian=='corrected_af'\nprint('ok')")
        for first, second in (('viscous', 'amplification_jacobian'), ('amplification_jacobian', 'viscous')):
            for entry in ('viscous.CoupledSolver', 'o.AmplificationJacobianSolver'):
                done = subprocess.run([sys.executable, '-c', code.format(first=first, second=second, entry=entry)],
                                      cwd=ROOT, text=True, capture_output=True)
                self.assertEqual((done.returncode, done.stdout.strip()), (0, 'ok'), done.stderr)


class GrowthDerivativeBranches(unittest.TestCase):
    """Value is untouched in every branch; the derivative is consistent on smooth branches."""
    def setUp(self):
        self.param = thermo_param(); self.param.amplification_jacobian = option.CORRECTED
        self.legacy = thermo_param(); self.legacy.amplification_jacobian = option.LEGACY

    def finite_difference(self, U, relative):
        base = regime(U, self.param); out = np.zeros(4)
        for j in range(4):
            h = relative*max(abs(U[j]), 1e-6); a, b = U.copy(), U.copy(); a[j] += h; b[j] -= h
            # A central difference is only meaningful when both samples stay on one branch.
            self.assertEqual((regime(a, self.param)[:2], regime(b, self.param)[:2]), (base[:2], base[:2]))
            out[j] = (option.get_damp(a, self.param)[0]-option.get_damp(b, self.param)[0])/(2*h)
        return out

    def test_inactive_ramp_and_full_branches(self):
        worst_legacy = 0.
        for hk in (2.3, 2.6, 3.5, 6., 11.):
            for rn, name in ((-.5, 'inactive'), (.25, 'ramp'), (.5, 'ramp'), (.75, 'ramp'), (1.5, 'full'), (4., 'full')):
                for amplification in (2., 8.9):
                    U = state_for(self.param, hk, rn, amplification)
                    self.assertEqual(regime(U, self.param)[:2], (name, False))
                    value, derivative = option.get_damp(U, self.param)
                    old = option._original_get_damp(U, self.legacy)
                    self.assertEqual(value, old[0])
                    if name == 'inactive': np.testing.assert_array_equal(derivative, old[1])
                    scale = max(1., float(abs(derivative).max()))
                    for relative in (1e-5, 1e-6):
                        fd = self.finite_difference(U, relative)
                        self.assertLess(abs(fd-derivative).max()/scale, 2e-5, (hk, rn, amplification, relative))
                        worst_legacy = max(worst_legacy, abs(fd-old[1]).max()/scale)
        # The historical derivative must visibly fail somewhere, or this test checks nothing.
        self.assertGreater(worst_legacy, 1e-3)

    def test_shape_factor_floor_inactive_and_active(self):
        for rn, name in ((-.5, 'inactive'), (.5, 'ramp'), (2., 'full')):
            U = state_for(self.param, 1.02, rn)
            self.assertEqual(regime(U, self.param)[:2], (name, True))
            value, derivative = option.get_damp(U, self.param)
            old = option._original_get_damp(U, self.legacy)
            # With Hk held at its floor the missing term multiplies a zero Hk derivative.
            self.assertEqual(value, old[0]); np.testing.assert_array_equal(derivative, old[1])
            fd = self.finite_difference(U, 1e-6)
            self.assertLess(abs(fd-derivative).max()/max(1., float(abs(derivative).max())), 2e-5)

    def test_branch_boundaries_are_continuous_in_value_and_one_sided_in_slope(self):
        # Activation start and saturation: the cubic ramp is C1, so value and derivative join.
        for hk in (2.6, 6.):
            for boundary in (0., 1.):
                below, above = (state_for(self.param, hk, boundary+s*1e-9) for s in (-1, 1))
                self.assertNotEqual(regime(below, self.param)[0], regime(above, self.param)[0])
                a, b = option.get_damp(below, self.param), option.get_damp(above, self.param)
                inside = option.get_damp(state_for(self.param, hk, .5), self.param)
                self.assertAlmostEqual(a[0], b[0], delta=1e-6*max(1., abs(inside[0])))
                # Compare the jump with the slope inside the ramp, not with a vanishing endpoint slope.
                self.assertLess(abs(a[1]-b[1]).max(), 1e-6*abs(inside[1]).max())
        # Shape-factor floor: value is continuous but the slope jumps (vendored TODO, unchanged here).
        below, above = state_for(self.param, 1.05-1e-9, 2.), state_for(self.param, 1.05+1e-9, 2.)
        self.assertEqual((regime(below, self.param)[1], regime(above, self.param)[1]), (True, False))
        a, b = option.get_damp(below, self.param), option.get_damp(above, self.param)
        self.assertAlmostEqual(a[0], b[0], delta=1e-5*max(1., abs(b[0])))
        for state, mode_value in ((below, a), (above, b)):
            old = option._original_get_damp(state, self.legacy)
            self.assertEqual(mode_value[0], old[0])


class PropagationAndProvenance(unittest.TestCase):
    def case(self, folder, **extra):
        path = Path(folder)/'case.json'
        path.write_text(json.dumps(dict(alpha_deg=[0.], reynolds=1e6, mach=0., ncrit=9,
            elements=[dict(name='main', naca='0012', panels=120)], **extra)))
        return path

    def metadata(self, out):
        meta = json.loads((Path(out)/'metadata.json').read_text())
        with (Path(out)/'polar.csv').open() as f: polar = list(csv.DictReader(f))
        return meta, polar

    def test_cli_config_and_failed_output_record_mode(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp); plain = self.case(tmp)
            (tmp/'legacy').mkdir(); configured = self.case(tmp/'legacy', amplification_jacobian='legacy')
            runs = [('default', [str(plain)], 'corrected_af', 0),
                    ('flag', [str(plain), '--amplification-jacobian', 'legacy'], 'legacy', 0),
                    ('config', [str(configured)], 'legacy', 0),
                    ('override', [str(configured), '--amplification-jacobian', 'corrected_af'], 'corrected_af', 0),
                    ('failed', [str(plain), '--iterations', '1'], 'corrected_af', 2),
                    ('failed_legacy', [str(plain), '--iterations', '1', '--amplification-jacobian', 'legacy'], 'legacy', 2)]
            forces = {}
            for name, arguments, mode, code in runs:
                with patch('builtins.print'):
                    self.assertEqual(multifoil.main(arguments+['--out', str(tmp/name)]), code, name)
                meta, polar = self.metadata(tmp/name)
                self.assertEqual(meta['case']['amplification_jacobian'], mode)
                self.assertEqual([d['amplification_jacobian'] for d in meta['diagnostics']], [mode])
                if code:
                    self.assertEqual(meta['diagnostics'][0]['status'], 'nonconverged')
                    self.assertEqual(polar[0]['converged'], '0')
                    # A failed solve keeps its mode but never a finite force.
                    self.assertTrue(all(not np.isfinite(float(polar[0][k])) for k in ('CL', 'CD', 'CM')))
                else:
                    self.assertEqual(meta['diagnostics'][0]['status'], 'converged')
                    forces[name] = [float(polar[0][k]) for k in ('CL', 'CD', 'CM')]
            # Derivative-only: both modes converge to the same symmetric solution here.
            np.testing.assert_allclose(forces['default'], forces['flag'], atol=1e-7)
            np.testing.assert_array_equal(forces['flag'], forces['config'])
            np.testing.assert_array_equal(forces['default'], forces['override'])

    def test_cli_rejects_bad_or_misplaced_mode_without_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp); plain = self.case(tmp)
            (tmp/'bad').mkdir(); bad = self.case(tmp/'bad', amplification_jacobian='typo')
            with patch('sys.stderr'), self.assertRaises(SystemExit) as stop:
                multifoil.main([str(plain), '--amplification-jacobian', 'typo', '--out', str(tmp/'a')])
            self.assertEqual(stop.exception.code, 2)
            with patch('builtins.print'):
                self.assertEqual(multifoil.main([str(bad), '--out', str(tmp/'b')]), 1)
                self.assertEqual(multifoil.main([str(plain), '--backend', 'panel', '--amplification-jacobian',
                                                 'legacy', '--out', str(tmp/'c')]), 1)
            self.assertFalse(any((tmp/name).exists() for name in 'abc'))

    def test_validation_runner_passes_declared_default(self):
        seen = {}
        class Stop(RuntimeError): pass
        def recorder(*args, **kwargs):
            seen.update(kwargs); raise Stop('captured')
        with tempfile.TemporaryDirectory() as tmp, patch.object(viscous, 'CoupledSolver', recorder), patch('builtins.print'):
            try: validate.main(['--case', 'naca0012_ladson', '--out', tmp])
            except Stop: pass
        self.assertEqual(seen.get('amplification_jacobian'), 'corrected_af')

    def test_optimizer_flow_job_passes_configured_mode(self):
        cfg = json.loads((ROOT/'examples/s702_incompressible_design.json').read_text())
        raw = [(e['name'], multifoil.read_dat(ROOT/'examples'/e['file'])) for e in cfg['elements']]
        for configured, expected in ((None, 'corrected_af'), ('legacy', 'legacy'), ('corrected_af', 'corrected_af')):
            seen = []
            def recorder(*args, **kwargs):
                seen.append(kwargs); raise ValueError('captured before any flow solve')
            trial = dict(cfg)
            if configured: trial['amplification_jacobian'] = configured
            with tempfile.TemporaryDirectory() as tmp, patch.object(shape_optimization, 'CoupledSolver', recorder):
                result = shape_optimization.flow_job((tmp, raw, trial, cfg['panels'], cfg['scenarios'][0], cfg['points'][0]))
            self.assertFalse(result['feasible']); self.assertIsNone(result['CD'])
            self.assertEqual(seen[0]['amplification_jacobian'], expected)

    def test_optimizer_manifest_hashes_correction_and_guards_resume(self):
        self.assertIn('amplification_jacobian.py', shape_optimization.CODE)
        config = ROOT/'examples/s702_incompressible_design.json'
        class Captured(Exception): pass
        class Reached(Exception): pass
        captured = {}
        def capture(path, value):
            captured['spec'] = copy.deepcopy(value); raise Captured()
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            with patch.object(shape_optimization, 'write_json', capture), self.assertRaises(Captured):
                shape_optimization.run(config, tmp/'fresh', 1)
            spec = captured['spec']
            self.assertEqual(spec['config']['amplification_jacobian'], 'corrected_af')
            self.assertEqual(set(spec['code_sha256']), set(shape_optimization.CODE))
            def study(name, value):
                (tmp/name).mkdir(); (tmp/name/'study.json').write_text(json.dumps(value)); return tmp/name
            def search(*args, **kwargs): raise Reached()
            # An identical manifest resumes: execution proceeds past the guard to the search.
            with patch.object(shape_optimization, 'pattern_search', search), self.assertRaises(Reached):
                shape_optimization.run(config, study('same', spec), 1)
            changed = copy.deepcopy(spec); changed['code_sha256']['amplification_jacobian.py'] = '0'*64
            other = copy.deepcopy(spec); other['config']['amplification_jacobian'] = 'legacy'
            for name, value in (('code', changed), ('mode', other)):
                with patch.object(shape_optimization, 'pattern_search', search), \
                     self.assertRaisesRegex(ValueError, 'Code, input or configuration changed'):
                    shape_optimization.run(config, study(name, value), 1)


if __name__ == '__main__':
    unittest.main()
