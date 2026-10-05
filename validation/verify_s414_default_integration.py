#!/usr/bin/env python3
"""V1 replay: the integrated default Jacobian at archived S414 CL 0.056/0.284 states.

No nonlinear flow solve. Saved states are restored through the standard
CoupledSolver constructor and compared with a frozen pre-integration snapshot
evaluated in a separate process. Historical manifests are not relaxed: the
current-source differences are enumerated and must equal the declared edit set.
"""
import os
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
import argparse
import importlib.metadata
import json
from pathlib import Path
import platform
import subprocess
import sys
import tempfile
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from report_design_initialization import ROOT, Inputs, sha, write, csv_write
from report_s702_wake_extent import audit, require, array_close
from s414_geometry import elements
from viscous import CoupledSolver, mf, te_info
import amplification_jacobian as option

OUT = ROOT/'validation/results/s414_default_integration_20261004'
STUDY = ROOT/'validation/results/s414_lower_bucket_20261003'
LAUNCH = ROOT/'validation/results/s414_lower_bucket_launch_20261003'
STUDY_HASH = '2c1d57eba8101dcbc8283f23032447bd18e656dc3f18812a8d0b4b17c82f8cec'
FROZEN = ROOT/'validation/results/s414_fore_reconstruction_20261004'
FROZEN_HASH = None  # filled from the launch preparation receipt below
PRIOR = ROOT/'validation/results/s414_jacobian_check_20261003'
CHILD = ROOT/'validation/s414_default_integration_frozen_child.py'
CASES = ('cl0056_nominal', 'cl028_control')
# Files intentionally edited by the 4 October default integration.
INTEGRATED = {'multifoil.py', 'viscous.py', 'shape_optimization.py', 'validate.py', 'amplification_jacobian.py'}
UNCHANGED_CORE = ('paneling.py', 'vendor/mfoil.py', 'optimization.py', 'fixed_lift_policy.py')
CRITERION, RESIDUAL_ALLOWANCE = 2e-5, 1e-10
VARIABLES = ('theta', 'delta_star', 'N', 'ue')
STEPS = (1e-4, 1e-5)  # the two relative steps of the archived optional-mode check
LADDER = (1e-4, 5e-5, 2e-5, 1e-5, 5e-6, 2e-6, 1e-6)


def restore(solver, state):
    """Same restoration as report_s414_wake_transfer.restore, for a supplied solver."""
    alpha = state['result']['alpha_deg']
    direction = np.array([np.cos(np.deg2rad(alpha)), np.sin(np.deg2rad(alpha))])
    for model in solver.models:
        model.oper.alpha = alpha; model.isol.gam = model.isol.gamref@direction
    solver._build_wakes(alpha)
    for model in solver.models:
        a = state['arrays'][model.geom.name]
        array_close(model.foil.x, a['xy'], 'Rebuilt foil geometry', 2e-12)
        array_close(model.wake.x, a['wake_xy'], 'Rebuilt inviscid wake geometry', 2e-11)
        array_close(model.wake.t, a['wake_t'], 'Rebuilt wake tangents', 2e-11)
        array_close(model.wake.s, a['wake_s'], 'Rebuilt wake arclength', 2e-11)
        model.glob.U = a['U'].copy(); model.glob.Nsys = a['U'].shape[1]
        model.vsol.turb = a['turbulent'].copy()
        model.isol.Istag = [int(a['lower'][0]), int(a['upper'][0])]
        model.isol.sgnue = np.ones(model.foil.N); model.isol.sgnue[a['lower']] = -1
        model.isol.sstag = float(model.foil.s[a['lower'][0]]+a['xi'][a['lower'][0]])
        mf.identify_surfaces(model); mf.stagpoint_move(model)
        array_close(model.isol.xi, a['xi'], 'Restored stagnation distances', 2e-12)
        for k, label in enumerate(('lower', 'upper')):
            require(np.array_equal(model.vsol.Is[k], a[label]), 'Restored branch indices differ')
        if abs(te_info(model.foil.x)[1]) < 1e-12: model.vsol.wgap = np.zeros(model.wake.N)
        else: mf.set_wake_gap(model)
    solver._build_source_coupling()
    return solver


def regime(U, param):
    """Growth-law branch at one station, using the vendored arithmetic."""
    hk = mf.get_Hk(U, param)[0]; clipped = bool(hk < 1.05)
    hk = max(hk, 1.05); hmi = 1/(hk-1)
    critical = 2.492*hmi**.43+.7*(np.tanh(14*hmi-9.24)+1)
    rn = (np.log(mf.get_Ret(U, param)[0])/np.log(10)-(critical-.1))/.2
    return ('inactive' if rn < 0 else 'ramp' if rn < 1 else 'full')+('_hk_floor' if clipped else ''), float(rn), float(hk)


def source_classification(inputs):
    """Enumerate, never waive, the differences between frozen and current sources."""
    require(sha(inputs.add(STUDY/'study.json')) == STUDY_HASH, 'Unexpected archived study')
    spec = inputs.json(STUDY/'study.json'); differs = []
    for name, digest in spec['sha256'].items():
        require(sha(inputs.add(STUDY/'source_snapshot'/name)) == digest, 'Archived snapshot changed: '+name)
        if sha(ROOT/name) != digest: differs.append(name)
    require(set(differs) <= INTEGRATED, 'Undeclared source difference: '+str(sorted(set(differs)-INTEGRATED)))
    receipt = inputs.json(LAUNCH/'retrieved_files_sha256.json')
    for name, digest in receipt.items():
        require(sha(inputs.add(STUDY/name)) == digest, 'Retrieved output changed: '+name)
    frozen_spec = inputs.json(FROZEN/'study.json'); snapshot = FROZEN/'source_snapshot'
    preparation = inputs.json(ROOT/'validation/results/s414_fore_reconstruction_launch_20261004/preparation.json')
    require(sha(FROZEN/'study.json') == preparation['study_sha256'], 'Frozen optional-implementation study changed')
    frozen_differs = []
    for name in sorted(INTEGRATED | set(UNCHANGED_CORE)):
        require(sha(inputs.add(snapshot/name)) == frozen_spec['sha256'][name], 'Frozen snapshot changed: '+name)
        if sha(inputs.add(ROOT/name)) != frozen_spec['sha256'][name]: frozen_differs.append(name)
    require(set(frozen_differs) == INTEGRATED, 'Integration edit set differs from declaration: '+str(frozen_differs))
    other = sorted(name for name, digest in frozen_spec['sha256'].items()
                   if name not in INTEGRATED and (ROOT/name).exists() and sha(ROOT/name) != digest)
    require(not other, 'Other frozen inputs changed since the optional-implementation study: '+str(other))
    return spec, dict(archived_state_study=str(STUDY.relative_to(ROOT)), archived_study_sha256=STUDY_HASH,
        archived_inputs_checked=len(spec['sha256']), archived_snapshot_intact=True,
        archived_inputs_differing_from_current_source=sorted(differs),
        retrieved_outputs_verified=len(receipt),
        frozen_optional_implementation=str(snapshot.relative_to(ROOT)),
        frozen_study_sha256=preparation['study_sha256'], frozen_inputs_checked=len(frozen_spec['sha256']),
        files_changed_by_integration=sorted(frozen_differs), other_frozen_inputs_changed=other,
        unchanged_numerical_core=list(UNCHANGED_CORE),
        note='Historical load_states guards correctly reject the current source; this receipt does not relax them.')


def frozen_reference(spec, states, raw, folder):
    cases = [dict(name=name, alpha_deg=states[name]['result']['alpha_deg'],
                  wake_length=next(j for j in spec['jobs'] if j['name'] == name)['wake_length']) for name in CASES]
    (folder/'cases.json').write_text(json.dumps(dict(cases=cases, settings={k: spec[k] for k in
        ('reynolds', 'mach', 'ncrit', 'iterations', 'residual_tolerance')})))
    for name in CASES:
        payload = {'raw_'+element: xy for element, xy in raw}
        for element, arrays in states[name]['arrays'].items():
            payload.update({element+'_'+key: value for key, value in arrays.items()})
        np.savez(folder/(name+'_input.npz'), **payload)
    done = subprocess.run([sys.executable, str(CHILD), str(FROZEN/'source_snapshot'), str(folder)],
                          cwd=folder, text=True, capture_output=True)
    require(done.returncode == 0, 'Frozen reference evaluation failed:\n'+done.stderr[-2000:])
    return json.loads((folder/'frozen_record.json').read_text())


def compare(label, a, b, rows, case, exact=True):
    same = bool(np.array_equal(a, b)); difference = float(np.max(abs(a-b)))
    rows.append(dict(case=case, comparison=label, exactly_equal=same, maximum_absolute_difference=difference,
                     required='exact' if exact else f'<= {RESIDUAL_ALLOWANCE:g}', shape='x'.join(map(str, a.shape))))
    return same if exact else difference <= RESIDUAL_ALLOWANCE


def probe_nodes(solver):
    """Archived selection rule plus one representative of every growth branch present."""
    fore = solver.models[0]; p = mf.build_param(fore, 0)
    ids = np.asarray(fore.vsol.Is[0]); lam = ids[fore.vsol.turb[ids] == 0]
    active = [int(j) for j in lam if mf.get_damp(fore.glob.U[:, j], p)[0] > 1e-4]
    original = sorted(set([max(active, key=lambda j: mf.get_Hk(fore.glob.U[:, j], p)[0]), int(lam[-1])]))
    chosen = {node: 'archived_rule' for node in original}
    interior = [int(j) for j in lam[1:-1] if int(j) not in fore.isol.Istag]
    groups = {}
    for node in interior:
        groups.setdefault(regime(fore.glob.U[:, node], p)[0], []).append(node)
    for name, nodes in sorted(groups.items()):
        # The middle station of each branch stays away from its own boundaries.
        rn = np.array([regime(fore.glob.U[:, n], p)[1] for n in nodes])
        target = -1. if name.startswith('inactive') else .5 if name.startswith('ramp') else float(np.median(rn))
        node = nodes[int(np.argmin(abs(rn-target)))]
        chosen.setdefault(node, 'branch_'+name)
    return chosen, original, {name: len(nodes) for name, nodes in groups.items()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=OUT)
    args = parser.parse_args(); out = args.out.resolve(); out.mkdir(parents=True, exist_ok=True)
    require(not (out/'replay_summary.json').exists(), 'Preserve completed replay; choose a new output directory')
    inputs = Inputs()
    spec, sources = source_classification(inputs)
    prior = inputs.json(PRIOR/'summary.json')
    require(prior['complete'] and prior['passed'] and prior['comparisons'] == 32, 'Earlier derivative evidence changed')
    prior_rows = inputs.csv(PRIOR/'assembled_derivatives.csv')
    runtime = inputs.json(STUDY/'runtime.json')
    states = {name: audit(STUDY/name, next(j for j in spec['jobs'] if j['name'] == name), spec, runtime, inputs)
              for name in CASES}
    raw = elements('nominal', spec['panels'])
    settings = dict(reynolds=spec['reynolds'], mach=spec['mach'], ncrit=spec['ncrit'],
                    iterations=spec['iterations'], tolerance=spec['residual_tolerance'])
    replay, derivatives, census, growth, summaries, refinement = [], [], [], [], [], []
    with tempfile.TemporaryDirectory(prefix='s414_default_integration_') as tmp:
        tmp = Path(tmp); frozen = frozen_reference(spec, states, raw, tmp)
        require(frozen['optional_default'] == option.LEGACY, 'Frozen snapshot is not the opt-in implementation')
        for name in CASES:
            job = next(j for j in spec['jobs'] if j['name'] == name); state = states[name]
            kwargs = dict(settings, wake_length=job['wake_length'])
            # Standard entry point, no mode argument: this is what every new run now uses.
            standard = restore(CoupledSolver(raw, **kwargs), state)
            require(standard.amplification_jacobian == option.CORRECTED
                    and all(m.param.amplification_jacobian == option.CORRECTED for m in standard.models),
                    'Standard constructor did not select the corrected derivative')
            legacy = restore(CoupledSolver(raw, amplification_jacobian=option.LEGACY, **kwargs), state)
            adapter = restore(option.AmplificationJacobianSolver(raw, **kwargs), state)
            adapter_legacy = restore(option.AmplificationJacobianSolver(raw, amplification_jacobian=option.LEGACY, **kwargs), state)
            saved = [m.glob.U.copy() for m in standard.models]; flags = [m.vsol.turb.copy() for m in standard.models]
            r1, j1 = standard._system(); r0, j0 = legacy._system()
            ra, ja = adapter._system(); rb, jb = adapter_legacy._system()
            load = lambda kind, part: np.load(tmp/f'{name}_{kind}_{part}.npy')
            ok = [compare('standard default residual vs frozen optional corrected', r1, load('corrected_af', 'R'), replay, name),
                  compare('standard default Jacobian vs frozen optional corrected', j1, load('corrected_af', 'J'), replay, name),
                  compare('compatibility adapter default residual vs standard default', ra, r1, replay, name),
                  compare('compatibility adapter default Jacobian vs standard default', ja, j1, replay, name),
                  compare('explicit legacy residual vs frozen optional legacy', r0, load('legacy', 'R'), replay, name),
                  compare('explicit legacy Jacobian vs frozen optional legacy', j0, load('legacy', 'J'), replay, name),
                  compare('explicit legacy residual vs pre-option historical solver', r0, load('plain', 'R'), replay, name),
                  compare('explicit legacy Jacobian vs pre-option historical solver', j0, load('plain', 'J'), replay, name),
                  compare('compatibility adapter legacy residual vs explicit legacy', rb, r0, replay, name),
                  compare('compatibility adapter legacy Jacobian vs explicit legacy', jb, j0, replay, name),
                  compare('corrected residual vs legacy residual (derivative-only allowance)', r1, r0, replay, name, exact=False)]
            require(all(ok), 'State replay differs from frozen implementation: '+name)
            blnorm = float(np.linalg.norm(r1[:3*standard.nsys])); coupling = float(np.max(abs(r1[3*standard.nsys:])))
            require(max(blnorm, coupling) <= spec['residual_tolerance'], 'Restored residual fails tolerance')
            jacobian_changed = int(np.count_nonzero(j1 != j0))
            require(jacobian_changed > 0, 'Corrected and legacy Jacobians are identical; nothing was tested')
            # Growth values (the equation) must be identical; only its derivative may differ.
            changed_derivatives = 0
            for model, old in zip(standard.models, legacy.models):
                for side, label in enumerate(('lower', 'upper')):
                    p1, p0 = mf.build_param(model, side), mf.build_param(old, side)
                    counts = {}
                    for node in model.vsol.Is[side]:
                        if model.vsol.turb[node]: continue
                        U = model.glob.U[:, node]; a, b = mf.get_damp(U, p1), mf.get_damp(U, p0)
                        c = option._original_get_damp(U, p0)
                        require(a[0] == b[0] == c[0] and np.array_equal(b[1], c[1]), 'Growth value or legacy derivative changed')
                        changed_derivatives += int(not np.array_equal(a[1], b[1]))
                        kind = regime(U, p1)[0]; counts[kind] = counts.get(kind, 0)+1
                    for kind, count in sorted(counts.items()):
                        census.append(dict(case=name, element=model.geom.name, side=label, branch=kind, laminar_stations=count))
            growth.append(dict(case=name, laminar_growth_values_identical=True, legacy_derivative_identical_to_vendored=True,
                               stations_with_changed_derivative=changed_derivatives))
            chosen, original, groups = probe_nodes(standard); fore = standard.models[0]; p = mf.build_param(fore, 0)
            def central(node, variable, step, size):
                try:
                    fore.glob.U[variable, node] = saved[0][variable, node]+step*size
                    upper = regime(fore.glob.U[:, node], p)[0]; plus = standard._system()[0]
                    fore.glob.U[variable, node] = saved[0][variable, node]-step*size
                    lower = regime(fore.glob.U[:, node], p)[0]; minus = standard._system()[0]
                finally:
                    fore.glob.U[:] = saved[0]
                return (plus-minus)/(2*step), upper, lower
            for node, reason in sorted(chosen.items()):
                require(node not in fore.isol.Istag, 'Probe must avoid moving stagnation coordinate')
                base = regime(saved[0][:, node], p)
                for variable in range(4):
                    direction = np.zeros(4*standard.nsys)
                    direction[4*node+variable] = size = max(abs(saved[0][variable, node]), 1e-6)
                    expected = j1@direction; old = j0@direction
                    scale = max(float(abs(expected).max()), 1e-10); exceeded = False
                    for step in STEPS:
                        numerical, upper, lower = central(node, variable, step, size)
                        derivatives.append(dict(case=name, node=node, selection=reason, x=float(fore.foil.x[0, node]),
                            branch=base[0], activation_coordinate=base[1], Hk=base[2],
                            # A central difference across a branch change is reported, never scored as smooth.
                            same_branch=bool(upper == lower == base[0]),
                            variable=VARIABLES[variable], relative_step=step,
                            corrected_relative_error=float(abs(numerical-expected).max()/scale),
                            legacy_relative_error=float(abs(numerical-old).max()/scale),
                            direction_scale=float(size), derivative_scale=scale))
                        exceeded = exceeded or derivatives[-1]['corrected_relative_error'] >= CRITERION
                    if exceeded:
                        # Exceeding the criterion triggers investigation, never a relaxed threshold:
                        # separate finite-difference truncation from a Jacobian error by step refinement.
                        previous = None
                        for step in LADDER:
                            numerical, upper, lower = central(node, variable, step, size)
                            error = float(abs(numerical-expected).max()/scale); worst = int(np.argmax(abs(numerical-expected)))
                            refinement.append(dict(case=name, node=node, selection=reason, branch=base[0],
                                variable=VARIABLES[variable], relative_step=step, same_branch=bool(upper == lower == base[0]),
                                corrected_relative_error=error, legacy_relative_error=float(abs(numerical-old).max()/scale),
                                observed_order=None if previous is None else float(np.log(previous[1]/error)/np.log(previous[0]/step)),
                                worst_residual_row=worst))
                            previous = (step, error)
            back, jback = standard._system()
            require(np.array_equal(back, r1) and np.array_equal(jback, j1), 'Standard system not restored exactly')
            require(all(np.array_equal(m.glob.U, u) and np.array_equal(m.vsol.turb, f)
                        for m, u, f in zip(standard.models, saved, flags)), 'Saved states changed')
            for model in standard.models:
                array_close(model.glob.U, state['arrays'][model.geom.name]['U'], 'Replay must not change the archived solution', 0.)
            summaries.append(dict(case=name, alpha_deg=state['result']['alpha_deg'], nsys=int(standard.nsys),
                restored_BL_residual_l2=blnorm, restored_coupling_residual_max=coupling,
                recorded_residual=state['residual'], jacobian_entries_changed_by_correction=jacobian_changed,
                archived_rule_nodes=original, probe_nodes={str(k): v for k, v in sorted(chosen.items())},
                fore_lower_interior_branch_counts=groups, saved_arrays_restored_exactly=True))
            print(name, 'replay exact; max corrected derivative error',
                  max(r['corrected_relative_error'] for r in derivatives if r['case'] == name and r['same_branch']), flush=True)
    smooth = [r for r in derivatives if r['same_branch']]; crossing = [r for r in derivatives if not r['same_branch']]
    archived = [r for r in smooth if r['selection'] == 'archived_rule']
    require(len(archived) == len(prior_rows) == 32, 'Archived 32-probe set was not reproduced')
    # Same probes as the earlier optional-mode check: errors must match that receipt closely.
    drift = max(abs(a['corrected_relative_error']-float(b['corrected_relative_error'])) for a, b in zip(
        sorted(archived, key=lambda r: (r['case'], r['node'], r['variable'], r['relative_step'])),
        sorted(prior_rows, key=lambda r: (r['case'], int(r['node']), r['variable'], float(r['relative_step'])))))
    csv_write(out/'state_replay.csv', replay); csv_write(out/'assembled_derivatives.csv', derivatives)
    csv_write(out/'growth_branch_census.csv', census)
    if refinement: csv_write(out/'step_refinement.csv', refinement)
    for name in ('amplification_jacobian.py', 'viscous.py', 'multifoil.py', 'validation/s414_geometry.py',
                 'validation/verify_s414_default_integration.py', 'validation/s414_default_integration_frozen_child.py',
                 'validation/report_s702_wake_extent.py', 'validation/report_design_initialization.py',
                 'validation/report_s702_drag_mechanism.py'):
        inputs.add(ROOT/name)
    write(out/'replay_input_manifest.json', dict(files=inputs.files))
    branches = sorted({r['branch'] for r in smooth})
    exceed = [r for r in smooth if r['corrected_relative_error'] >= CRITERION]
    # The archived criterion is applied unchanged to the archived probe set at both steps.
    archived_pass = all(r['corrected_relative_error'] < CRITERION for r in archived)
    investigations = []
    for key in sorted({(r['case'], r['node'], r['variable']) for r in exceed}):
        ladder = [r for r in refinement if (r['case'], r['node'], r['variable']) == key]
        orders = [r['observed_order'] for r in ladder if r['observed_order'] is not None]
        finest = ladder[-1]; fine = next(r for r in ladder if r['relative_step'] == 1e-5)
        truncation = bool(all(r['same_branch'] for r in ladder) and all(1.8 <= o <= 2.2 for o in orders)
                          and fine['corrected_relative_error'] < CRITERION and finest['corrected_relative_error'] < CRITERION/100)
        investigations.append(dict(case=key[0], node=key[1], variable=key[2],
            coarse_step_error=ladder[0]['corrected_relative_error'], legacy_coarse_step_error=ladder[0]['legacy_relative_error'],
            error_at_1e_5=fine['corrected_relative_error'], finest_step_error=finest['corrected_relative_error'],
            legacy_finest_step_error=finest['legacy_relative_error'],
            observed_orders=orders, classified_as_finite_difference_truncation=truncation))
    passed = bool(smooth and archived_pass and all(i['classified_as_finite_difference_truncation'] for i in investigations)
                  and not crossing and all(r['exactly_equal'] or r['required'] != 'exact' for r in replay))
    summary = dict(complete=True, passed=passed, criterion=CRITERION, residual_allowance=RESIDUAL_ALLOWANCE,
        entry_point='viscous.CoupledSolver default constructor (amplification_jacobian omitted)',
        sources=sources, frozen_child_modules=frozen['loaded'], cases=summaries, growth=growth,
        state_replay_comparisons=len(replay), state_replay_exact=sum(r['exactly_equal'] for r in replay),
        maximum_corrected_vs_legacy_residual_change=max(r['maximum_absolute_difference'] for r in replay if r['required'] != 'exact'),
        derivative_comparisons=len(derivatives), smooth_branch_comparisons=len(smooth),
        branch_crossing_comparisons=len(crossing),
        archived_rule_comparisons=len(archived), archived_rule_all_below_criterion=archived_pass,
        archived_rule_maximum_corrected_relative_error=max(r['corrected_relative_error'] for r in archived),
        archived_rule_maximum_drift_from_prior_receipt=drift,
        added_branch_comparisons=len(smooth)-len(archived),
        declared_step_comparisons_at_or_above_criterion=len(exceed),
        exceedance_investigations=investigations,
        exceedance_rule='Criterion unchanged. A declared-step exceedance passes only if step refinement shows '
            'second-order convergence (1.8-2.2) on one branch to below the criterion at 1e-5 and below criterion/100 at 1e-6.',
        assembled_branches_covered=branches,
        branches_not_present_in_saved_states='Hk floor (Hk < 1.05) does not occur at any saved laminar station; '
            'covered by isolated component tests in test_default_integration.py',
        maximum_corrected_relative_error_all_declared_steps=max(r['corrected_relative_error'] for r in smooth),
        maximum_corrected_relative_error_at_1e_5=max(r['corrected_relative_error'] for r in smooth if r['relative_step'] == 1e-5),
        maximum_legacy_relative_error=max(r['legacy_relative_error'] for r in smooth),
        prior_optional_mode_maximum=prior['maximum_corrected_relative_error'],
        new_nonlinear_flow_solves=0, runtime=dict(python=platform.python_version(),
            versions={k: importlib.metadata.version(k) for k in ('numpy', 'scipy')}),
        qualification='Implementation verification only: no aerodynamic benefit or solution uniqueness is claimed.')
    write(out/'replay_summary.json', summary)
    require(passed, 'Default-integration replay failed; preserve receipt and stop nonlinear interpretation')
    print('PASS', len(replay), 'replay comparisons;', len(smooth), 'derivative comparisons;', len(exceed),
          'declared-step exceedances resolved as truncation; max error at 1e-5',
          summary['maximum_corrected_relative_error_at_1e_5'])


if __name__ == '__main__':
    main()
