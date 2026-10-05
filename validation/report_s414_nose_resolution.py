#!/usr/bin/env python3
"""Frozen numerical screens and graphs for the audited nose-resolution pilot.

All differences compare identical incidences and common physical stations on
the same reconstructed polygon. Lift may change: these are not matched-lift
drag errors. A failed convergence screen is a scientific result, not an excuse
to retry a case or change acceptance criteria.
"""
import argparse
import csv
import json
from pathlib import Path
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import run_s414_nose_resolution as run

OUT = ROOT/'validation/results/s414_nose_resolution_analysis_20261005'
SCREENS = dict(Cp=0.01, theta_relative=0.02, delta_star_relative=0.02, N=0.1, CD_counts=1.)
NOISE_FLOORS = dict(Cp=1e-5, theta_relative=1e-4, delta_star_relative=1e-4, N=1e-3, CD_counts=.01)
WAKE_MOVEMENT_LIMIT = 1e-5
COLORS = {1: '#335c81', 2: '#d66b31', 4: '#267b58'}


def require(condition, message):
    if not condition: raise ValueError(message)


def maximum(values):
    return float(np.max(np.abs(values))) if len(values) else None


def supported_interpolation(s, values, branch, query, turbulent=None):
    """Linear interpolation with explicit domain and two-node laminar support.

At an exact node, that node alone supports its value. Outside the node range
the returned mask is false; callers must never treat np.interp's clamp as data.
    """
    s, values, branch, query = map(np.asarray, (s, values, branch, query))
    order = branch[np.argsort(s[branch])]; domain = s[order]
    require(len(domain) >= 2 and np.all(np.diff(domain) > 0), 'Invalid branch coordinates')
    supported = (query >= domain[0]) & (query <= domain[-1])
    right = np.searchsorted(domain, query, side='left').clip(0, len(domain)-1)
    left = np.maximum(right-1, 0)
    exact = domain[right] == query; left[exact] = right[exact]
    if turbulent is not None:
        flags = np.asarray(turbulent).astype(bool)
        supported &= ~flags[order[left]] & ~flags[order[right]]
    return np.interp(query, domain, values[order]), supported


def contraction(first, second, floor):
    if first is None or second is None: return None
    return bool(second <= (floor if first <= floor else first))


def location_on_polygon(xy, point):
    """Physical (chord-sum) arc and enclosing cell width; independent of foil.s."""
    xy, point = np.asarray(xy), np.asarray(point)
    delta = np.diff(xy, axis=1).T; squared = np.sum(delta*delta, axis=1)
    require(np.all(squared > 0), 'Degenerate polygon cell')
    fraction = np.clip(np.sum((point-xy[:, :-1].T)*delta, axis=1)/squared, 0., 1.)
    projected = xy[:, :-1].T+fraction[:, None]*delta
    cell = int(np.argmin(np.linalg.norm(projected-point, axis=1)))
    require(np.linalg.norm(projected[cell]-point) < 1e-9, 'Reported location is off the frozen polygon')
    lengths = np.sqrt(squared)
    return dict(physical_s=float(np.sum(lengths[:cell])+fraction[cell]*lengths[cell]),
                cell_width=float(lengths[cell]), cell=cell)


def compare_topology(a, b):
    sides = []; same = True; within = True
    for element in ('fore', 'aft'):
        for side in ('lower', 'upper'):
            aa = a['result']['boundary_layers'][element][side]; bb = b['result']['boundary_layers'][element][side]
            ta, tb = aa['transition'], bb['transition']
            topology = ((ta is None) == (tb is None) and
                [r['reaches_trailing_edge'] for r in aa['separation_regions']] ==
                [r['reaches_trailing_edge'] for r in bb['separation_regions']])
            row = dict(element=element, side=side, same_topology=topology,
                physical_transition_shift=None, fine_cell_width=None, transition_within_fine_cell=None)
            if ta is not None and tb is not None:
                la = location_on_polygon(a['arrays'][element]['xy'], ta['xy_over_cref'])
                lb = location_on_polygon(b['arrays'][element]['xy'], tb['xy_over_cref'])
                shift = abs(la['physical_s']-lb['physical_s'])
                row.update(physical_transition_shift=shift, fine_cell_width=lb['cell_width'],
                    transition_within_fine_cell=bool(shift < lb['cell_width']))
                within &= row['transition_within_fine_cell']
            elif ta is not None or tb is not None:
                within = False
            sides.append(row); same &= topology
    return dict(same_topology=bool(same), transition_within_fine_cell=bool(within), sides=sides)


def common_profiles(states):
    """Use every fine lower nose station; support and N mask common to all grids."""
    fine = states[4]; fa = fine['arrays']['fore']; lower = fa['lower']
    ids = lower[fa['xy'][0, lower] <= .04]
    query = np.asarray(fine['surface_s']['fore'])[ids]
    order = np.argsort(query); query, ids = query[order], ids[order]
    profiles = {}; supported = np.ones(len(query), dtype=bool); laminar = supported.copy()
    for factor, state in states.items():
        arr = state['arrays']['fore']; s = np.asarray(state['surface_s']['fore']); branch = arr['lower']
        values = {'Cp': state['surfaces']['fore']['Cp'][::-1], 'theta': arr['U'][0],
                  'delta_star': arr['U'][1], 'N': arr['U'][2], 'cf': state['surfaces']['fore']['cf'][::-1]}
        profiles[factor] = {}
        for key, value in values.items():
            interp, mask = supported_interpolation(s, value, branch, query, arr['turbulent'] if key == 'N' else None)
            profiles[factor][key] = interp
            if key == 'N': laminar &= mask
            else: supported &= mask
    laminar &= supported
    return dict(s=query, xy=fa['xy'][:, ids], profiles=profiles, supported=supported, laminar=laminar,
        candidate_stations=len(query), common_supported_stations=int(np.sum(supported)),
        common_laminar_stations=int(np.sum(laminar)))


def pair_metrics(a, b, common, first, second):
    aa, bb = common['profiles'][first], common['profiles'][second]
    mask = common['supported']; laminar = common['laminar']
    values = dict(Cp=maximum((aa['Cp']-bb['Cp'])[mask]),
        theta_relative=maximum(((aa['theta']-bb['theta'])/np.maximum(abs(bb['theta']), 1e-12))[mask]),
        delta_star_relative=maximum(((aa['delta_star']-bb['delta_star'])/np.maximum(abs(bb['delta_star']), 1e-12))[mask]),
        N=maximum((aa['N']-bb['N'])[laminar]), CD_counts=abs(a['result']['CD']-b['result']['CD'])*1e4)
    wake = {e: maximum(np.linalg.norm(a['arrays'][e]['wake_xy']-b['arrays'][e]['wake_xy'], axis=0)) for e in ('fore', 'aft')}
    global_cp = {}
    for e in ('fore', 'aft'):
        sa, sb = np.asarray(a['surface_s'][e]), np.asarray(b['surface_s'][e])
        qa = a['surfaces'][e]['Cp'][::-1]; qb = b['surfaces'][e]['Cp'][::-1]
        require(sb[0] >= sa[0] and sb[-1] <= sa[-1], 'Surface parameter extent changed')
        global_cp[e] = maximum(np.interp(sb, sa, qa)-qb)
    return dict(from_factor=first, to_factor=second, metrics=values, topology=compare_topology(a, b),
        wake_maximum_movement_by_element=wake, maximum_wake_movement=max(wake.values()),
        all_surface_maximum_Cp_difference=global_cp,
        delta_CL=b['result']['CL']-a['result']['CL'], delta_CD_counts=(b['result']['CD']-a['result']['CD'])*1e4)


def assess(states):
    if set(states) != {1, 2, 4}:
        return dict(numerical_screen='unknown', isolated_nose_refinement=None,
            reason='At least one selected state is missing; no refinement order is inferred', pairs=[])
    common = common_profiles(states)
    pairs = [pair_metrics(states[a], states[b], common, a, b) for a, b in ((1, 2), (2, 4))]
    criteria = {}
    for key, threshold in SCREENS.items():
        first, last = (p['metrics'][key] for p in pairs)
        criteria[key] = dict(value=last, threshold=threshold, noise_floor=NOISE_FLOORS[key],
            within_screen=None if last is None else bool(last <= threshold),
            contraction=contraction(first, last, NOISE_FLOORS[key]))
    # Freeze the distinction explicitly: topology must agree on all four
    # surfaces across BOTH adjacent grid pairs. Transition movement is a
    # finest-grid convergence screen (2 -> 4); 1 -> 2 movement is reported.
    topology = dict(pairs[-1]['topology'])
    topology['same_topology'] = all(p['topology']['same_topology'] for p in pairs)
    topology['topology_screen_pairs'] = [[1, 2], [2, 4]]
    topology['transition_movement_screen_pair'] = [2, 4]
    topology['coarse_transition_movement_role'] = 'reported_only'
    statuses = [v[k] for v in criteria.values() for k in ('within_screen', 'contraction')]
    statuses += [topology['same_topology'], topology['transition_within_fine_cell']]
    outcome = 'fail' if False in statuses else 'unknown' if None in statuses else 'pass'
    # Baseline-to-fine movement can exceed either adjacent increment; include it.
    wake = {e: maximum(np.linalg.norm(states[1]['arrays'][e]['wake_xy']-states[4]['arrays'][e]['wake_xy'], axis=0))
            for e in ('fore', 'aft')}
    maximum_wake = max(max(wake.values()), *(p['maximum_wake_movement'] for p in pairs))
    return dict(numerical_screen=outcome, isolated_nose_refinement=maximum_wake <= WAKE_MOVEMENT_LIMIT,
        maximum_wake_movement=maximum_wake, wake_movement_limit=WAKE_MOVEMENT_LIMIT,
        baseline_to_fine_wake_movement=wake, criteria=criteria, pairs=pairs,
        sampling={k: common[k] for k in ('candidate_stations', 'common_supported_stations', 'common_laminar_stations')},
        topology=topology, interpretation='Fixed-incidence grid sensitivity; different CL, not matched-lift accuracy')


def csv_write(path, rows):
    if not rows: return
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)


def figures(out, cases):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    names = ('cl0056', 'cl0284'); artifacts = []
    fig, axes = plt.subplots(5, 2, figsize=(10, 13), layout='constrained', sharex='col')
    labels = [('Cp', 'Cp'), ('theta', 'Momentum thickness / c'), ('delta_star', 'Displacement thickness / c'),
              ('N', 'Laminar amplification N'), ('cf', 'Skin friction Cf')]
    for col, name in enumerate(names):
        states = cases[name]
        for factor, state in sorted(states.items()):
            a = state['arrays']['fore']; ids = a['lower']; ids = ids[a['xy'][0, ids] <= .04]
            values = dict(Cp=state['surfaces']['fore']['Cp'][::-1], theta=a['U'][0], delta_star=a['U'][1],
                N=a['U'][2], cf=state['surfaces']['fore']['cf'][::-1])
            for row, (key, label) in enumerate(labels):
                use = ids[~a['turbulent'][ids].astype(bool)] if key == 'N' else ids
                axes[row, col].plot(a['xy'][0, use], values[key][use], '.-', color=COLORS[factor], ms=3, lw=1,
                    label=f'nose spacing / {factor}')
                axes[row, col].set_ylabel(label); axes[row, col].grid(alpha=.2)
            axes[0, col].set_title(f"{name}: fixed α={state['result']['alpha_deg']:.7f}°")
        axes[0, col].invert_yaxis(); axes[-1, col].set_xlabel('x / reference chord')
        axes[0, col].legend(fontsize=8)
    fig.suptitle('Same polygon and fixed incidence; lift may differ between grids')
    for suffix in ('png', 'pdf'):
        p = out/f'nose_profiles.{suffix}'; fig.savefig(p, dpi=180); artifacts.append(p)
    plt.close(fig)
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), layout='constrained')
    for col, name in enumerate(names):
        for row, element in enumerate(('fore', 'aft')):
            for factor, state in sorted(cases[name].items()):
                surface = state['surfaces'][element]
                axes[row, col].plot(surface['x_over_cref'], surface['Cp'], '-', color=COLORS[factor],
                    lw=1, label=f'factor {factor}')
            axes[row, col].set_title(f'{name}: {element}'); axes[row, col].set_ylabel('Cp')
            axes[row, col].set_xlabel('x / reference chord'); axes[row, col].invert_yaxis()
            axes[row, col].grid(alpha=.2); axes[row, col].legend(fontsize=8)
    for suffix in ('png', 'pdf'):
        p = out/f'full_surface_pressure.{suffix}'; fig.savefig(p, dpi=180); artifacts.append(p)
    plt.close(fig)
    fig, axes = plt.subplots(3, 2, figsize=(9, 8), layout='constrained')
    for col, name in enumerate(names):
        states = cases[name]
        for row, (key, label, multiplier) in enumerate((('CL', 'Lift coefficient', 1.), ('CD', 'Drag [counts]', 1e4))):
            factors = sorted(states)
            axes[row, col].plot(factors, [states[f]['result'][key]*multiplier for f in factors], 'o-', color='#335c81')
            axes[row, col].set_ylabel(label)
        for factor, state in sorted(states.items()):
            tr = state['result']['boundary_layers']['fore']['lower']['transition']
            if tr: axes[2, col].plot(factor, tr['xy_over_cref'][0], 'o', color=COLORS[factor])
        axes[0, col].set_title(name); axes[2, col].set_ylabel('Fore-lower transition x / c')
        for ax in axes[:, col]: ax.set_xticks([1, 2, 4]); ax.grid(alpha=.2)
        axes[2, col].set_xlabel('Subdivision factor; failed cases omitted, never assigned forces')
    for suffix in ('png', 'pdf'):
        p = out/f'fixed_incidence_response.{suffix}'; fig.savefig(p, dpi=180); artifacts.append(p)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(9, 4), layout='constrained')
    available = next((cases[n] for n in names if cases[n]), {})
    for factor, state in sorted(available.items()):
        xy = state['arrays']['fore']['xy']; mask = xy[0] <= .04
        ax.plot(xy[0, mask], xy[1, mask], '.-', ms=3, lw=.7, color=COLORS[factor], label=f'factor {factor}')
    ax.set_aspect('equal'); ax.set_xlabel('x / c'); ax.set_ylabel('y / c'); ax.grid(alpha=.2)
    ax.set_title('Identical nose polygon; collinear inserted stations'); ax.legend()
    p = out/'nose_geometry.png'; fig.savefig(p, dpi=180); artifacts.append(p); plt.close(fig)
    return artifacts


def report(study, out):
    audit_path = out/'residual_recheck.json'; audit = run.read(audit_path)
    require(audit['complete'] and audit['evidence_verified'] and len(audit['states']) == audit['saved_states'], 'Saved-state audit incomplete')
    for name, expected in audit['input_sha256'].items(): require(run.sha(ROOT/name) == expected, 'Audited input changed: '+name)
    spec = run.read(study/'study.json'); cases = {'cl0056': {}, 'cl0284': {}}; rows = []; extrema = []
    for job in spec['jobs']:
        folder = study/job['name']; path = folder/'result.json'
        result = run.read(path) if path.exists() else None
        row = dict(job=job['name'], case=job['pressure_case'], subdivision_factor=job['subdivision_factor'],
            requested_alpha_deg=job['alpha_deg'], selected_converged=bool(result and result['feasible']),
            CL=None, CD_counts=None, fore_CD_counts=None, aft_CD_counts=None, fore_lower_transition_x=None)
        if result and result['feasible']:
            index = result['selected_flow_index']; arrays = {}; surfaces = {}
            for e in ('fore', 'aft'):
                with np.load(folder/(e+'_state.npz'), allow_pickle=False) as saved: arrays[e] = {k: saved[k].copy() for k in saved.files}
                surfaces[e] = np.atleast_1d(np.genfromtxt(folder/f'flow/surface_{index:03d}_{e}.csv', names=True, delimiter=',', dtype=None, encoding='utf-8'))
            recheck = next(r for r in audit['states'] if r['job'] == job['name'] and r['index'] == index)
            cases[job['pressure_case']][job['subdivision_factor']] = dict(result=result, arrays=arrays, surfaces=surfaces, surface_s=recheck['surface_s'])
            tr = result['boundary_layers']['fore']['lower']['transition']
            row.update(CL=result['CL'], CD_counts=result['CD']*1e4,
                fore_CD_counts=result['element_forces']['fore']['CD_wake']*1e4,
                aft_CD_counts=result['element_forces']['aft']['CD_wake']*1e4,
                fore_lower_transition_x=tr['xy_over_cref'][0] if tr else None)
            for e in ('fore', 'aft'):
                for side in ('lower', 'upper'):
                    values = surfaces[e][surfaces[e]['side'] == side]
                    layer = result['boundary_layers'][e][side]
                    extrema.append(dict(job=job['name'], element=e, side=side, nodes=len(values),
                        minimum_Cp=float(np.min(values['Cp'])), maximum_Cp=float(np.max(values['Cp'])),
                        minimum_Cf=float(np.min(values['cf'])), maximum_Hk=float(np.max(values['Hk'])),
                        maximum_theta=float(np.max(values['theta'])), maximum_delta_star=float(np.max(values['delta_star'])),
                        transition_present=layer['transition'] is not None, separation_intervals=len(layer['separation_regions']),
                        element_CL=result['element_forces'][e]['CL'], element_CD_counts=result['element_forces'][e]['CD_wake']*1e4))
        rows.append(row)
    screens = {name: assess(states) for name, states in cases.items()}
    csv_write(out/'selected_cases.csv', rows); artifacts = [out/'selected_cases.csv', audit_path]
    if extrema:
        path = out/'full_surface_extrema.csv'; csv_write(path, extrema); artifacts.append(path)
    for name, states in cases.items():
        if set(states) != {1, 2, 4}: continue
        common = common_profiles(states); profile_rows = []
        for i, value in enumerate(common['s']):
            if not common['supported'][i]: continue
            row = dict(foil_parameter_s=float(value), x_over_c=float(common['xy'][0, i]), y_over_c=float(common['xy'][1, i]),
                joint_laminar_support=bool(common['laminar'][i]))
            for factor, profiles in common['profiles'].items():
                for key, data in profiles.items(): row[f'factor{factor}_{key}'] = float(data[i]) if key != 'N' or common['laminar'][i] else None
            profile_rows.append(row)
        path = out/(name+'_common_profiles.csv'); csv_write(path, profile_rows); artifacts.append(path)
    artifacts += figures(out, cases)
    numerical = [s['numerical_screen'] for s in screens.values()]
    overall = 'fail' if 'fail' in numerical else 'unknown' if 'unknown' in numerical else 'pass'
    gates = dict(execution_complete=audit['execution_state'] == 'finished', evidence_verified=True,
        exact_control_gate=audit['control_gate_passed'], selected_coverage=f'{sum(r["selected_converged"] for r in rows)}/6',
        numerical_screen=overall, isolated_nose_refinement=(None if any(s['isolated_nose_refinement'] is None for s in screens.values())
            else all(s['isolated_nose_refinement'] for s in screens.values())))
    lines = ['# S414 local coupled nose-resolution pilot', '',
        f"Execution: {audit['execution_state']}. Evidence verified: yes ({audit['saved_states']} saved states).",
        f"Selected coverage: {gates['selected_coverage']}. Exact baseline controls: {audit['control_gate_passed']}.", '',
        'The comparison fixes incidence and the nominal reconstructed polygon. Refined CL values may differ; these drag differences are not matched-lift experimental errors.', '',
        '| Archived case label | Numerical screen | Nose-only attribution |', '|---|---|---|']
    for name, screen in screens.items(): lines.append(f"| {name} | {screen['numerical_screen']} | {screen['isolated_nose_refinement']} |")
    lines += ['', 'All fine-grid fore-lower stations at x/c ≤ 0.04 enter the declared common-domain pressure/thickness screen. Amplification uses only stations with laminar bracketing support in all three grids; turbulent shear is never interpreted as N.', '',
        'Interpolation uses the frozen solver surface parameter. Transition shifts and enclosing-cell widths use physical polygon arc. Traceable wake coordinates may move despite fixed wake controls; movement above 1e-5c marks causal attribution confounded.', '',
        'No physical closure, Ncrit or geometry correction follows automatically. A failed numerical screen or missing state is retained; it does not authorize additional simulations. Independent scientific review and the next decision remain pending.', '',
        '[Nose profiles](nose_profiles.png) · [Full surface pressure](full_surface_pressure.png) · [Fixed-incidence response](fixed_incidence_response.png) · [Geometry](nose_geometry.png)', '']
    findings = out/'FINDINGS.md'; findings.write_text('\n'.join(lines)); artifacts.append(findings)
    result = dict(schema=1, evidence_verified=True, criteria=gates, scientific_cases=screens,
        qualification='Fixed-incidence diagnostic; no matched-lift accuracy claim; no physical correction selected',
        independent_review='pending', new_flow_solves=0, study_sha256=run.sha(study/'study.json'),
        audit_sha256=run.sha(audit_path), artifact_sha256={str(p.relative_to(out)): run.sha(p) for p in artifacts},
        artifacts=[str(p) for p in artifacts])
    run.write(out/'report.json', result); return result


def orchestration_result(result, out):
    """Scientific failures are completed evidence, never a request for a retry."""
    gates = result['criteria']; evidence = [dict(path=str(Path(p).resolve()), sha256=run.sha(p)) for p in result['artifacts']]
    evidence.append(dict(path=str((out/'report.json').resolve()), sha256=run.sha(out/'report.json')))
    references = [str((out/'residual_recheck.json').resolve()), str((out/'report.json').resolve())]
    flag = lambda value: 'inconclusive' if value is None else ('pass' if value else 'fail')
    rows = [('execution_complete', flag(gates['execution_complete'])), ('evidence_integrity', 'pass'),
        ('archived_controls', flag(gates['exact_control_gate'])),
        ('selected_coverage', 'pass' if gates['selected_coverage'] == '6/6' else 'inconclusive'),
        ('numerical_screen', 'inconclusive' if gates['numerical_screen'] == 'unknown' else gates['numerical_screen']),
        ('isolated_nose_refinement', flag(gates['isolated_nose_refinement']))]
    return dict(evidence_verified=True, criteria=[dict(id=key, outcome=value, evidence=references) for key, value in rows],
        evidence=evidence, scientific_interpretation='Independent review pending; no correction or additional runs authorized')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--study', type=Path, default=ROOT/run.STUDY)
    parser.add_argument('--out', type=Path, default=OUT)
    parser.add_argument('--launch', type=Path, default=ROOT/run.LAUNCH)
    parser.add_argument('--audit', action='store_true', help='Perform fresh saved-state audit before reporting; never solves a flow')
    parser.add_argument('--result', type=Path, help='Structured orchestrator result; scientific failures keep exit status zero')
    args = parser.parse_args()
    if args.audit:
        from audit_s414_nose_resolution import audit
        audit(args.study.resolve(), args.launch.resolve(), args.out.resolve())
    result = report(args.study.resolve(), args.out.resolve())
    if args.result: run.write(args.result, orchestration_result(result, args.out.resolve()))
    print(json.dumps(dict(evidence_verified=result['evidence_verified'], criteria=result['criteria'], artifacts=result['artifacts'])))
    return 0


if __name__ == '__main__': raise SystemExit(main())
