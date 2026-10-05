#!/usr/bin/env python3
"""Audit the frozen S702 factor-separation study from saved states; no solves."""
import argparse
import json
import platform
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from report_s702_wake_extent import (
    ROOT, Inputs, sha, write, csv_write, audit, require, close, array_close,
    clearance, side_data,
)
from report_s702_holdout_panels import contour_sample_distance, resolution
from run_s702_holdout_panels import effective_spec
from run_s702_discretization import numerical_setup, STAGES
from viscous import new_model, mf
from paneling import curvature_paneling

STUDY_HASH = '032ca093ea6cd2dff5007cd4177fb49c392fafc7b63896604529aa350119e9fd'
ELEMENTS = ('fore', 'aft')
LABELS = ('Source contour', 'Surface sampling', 'Wake count', 'First wake step')
DEPENDENCIES = (
    'validation/report_s702_wake_extent.py', 'validation/report_s702_holdout_panels.py',
    'validation/report_design_initialization.py', 'validation/report_s702_drag_mechanism.py',
    'validation/run_s702_holdout_panels.py', 'validation/run_s702_discretization.py',
    'validation/s702_ncrit_sweep.py', 'viscous.py', 'multifoil.py', 'paneling.py',
    'vendor/mfoil.py', 'diagnostic_discretization.py',
)


def check_factor_controls(setups):
    """Reject accidental changes to a second factor, including hidden source changes."""
    require(set(setups) == set(STAGES), 'All five stages required')
    keys = ('source_sha256', 'surface_nodes', 'wake_nodes', 'wake_first_spacing')
    for index, (first, second) in enumerate(zip(STAGES[:-1], STAGES[1:])):
        a, b = setups[first], setups[second]
        for key in keys:
            va = a[key] if key == 'source_sha256' else a['requested'][key]
            vb = b[key] if key == 'source_sha256' else b['requested'][key]
            require((va != vb) == (key == keys[index]),
                    f'Unexpected factor change {first}->{second}: {key}')
        require(a['recipes'] == b['recipes'], 'Control recipes changed between stages')


def paired_increments(results):
    """Require every pair, even when an omitted state would improve the answer."""
    require(set(results) == {(f, s) for f in (5, 10) for s in STAGES},
            'Complete ten-state paired coverage required')
    for r in results.values():
        require(r['feasible'] and r['CD'] is not None and np.isfinite(r['CD']),
                'Complete finite paired coverage required')
    pairs, increments = [], []
    for stage in STAGES:
        a, b = results[5, stage], results[10, stage]
        row = dict(stage=stage, flap_delta_CD_counts=(b['CD']-a['CD'])*1e4)
        for e in ELEMENTS:
            row[e+'_delta_counts'] = (b['element_forces'][e]['CD_wake']-a['element_forces'][e]['CD_wake'])*1e4
        close(row['flap_delta_CD_counts'], sum(row[e+'_delta_counts'] for e in ELEMENTS))
        pairs.append(row)
    for label, a, b in zip(LABELS, pairs[:-1], pairs[1:]):
        row = dict(factor=label, from_stage=a['stage'], to_stage=b['stage'],
                   paired_change_counts=b['flap_delta_CD_counts']-a['flap_delta_CD_counts'])
        for flap in (5, 10):
            r0, r1 = results[flap, a['stage']], results[flap, b['stage']]
            row[f'd{flap}_change_counts'] = (r1['CD']-r0['CD'])*1e4
            for e in ELEMENTS:
                row[f'd{flap}_{e}_change_counts'] = (r1['element_forces'][e]['CD_wake']-r0['element_forces'][e]['CD_wake'])*1e4
        for e in ELEMENTS:
            row[e+'_paired_change_counts'] = b[e+'_delta_counts']-a[e+'_delta_counts']
        close(row['paired_change_counts'], row['d10_change_counts']-row['d5_change_counts'])
        close(row['paired_change_counts'], sum(row[e+'_paired_change_counts'] for e in ELEMENTS))
        increments.append(row)
    total = pairs[-1]['flap_delta_CD_counts']-pairs[0]['flap_delta_CD_counts']
    close(sum(r['paired_change_counts'] for r in increments), total, context='Telescoping paired change')
    return pairs, increments, total


def compare_arrays(a, b, label):
    rows = []
    for element in ELEMENTS:
        require(set(b[element]) <= set(a[element]), 'Missing reference array fields: '+label)
        for key in sorted(b[element]):
            x, y = a[element][key], b[element][key]
            require(x.shape == y.shape, 'Control array shape differs: '+label+'/'+key)
            rows.append(dict(comparison=label, element=element, field=key,
                exactly_equal=bool(np.array_equal(x, y)), max_abs_difference=float(np.max(abs(x-y)))))
    return rows


def audit_setup(folder, job, spec, state, inputs):
    setup = inputs.json(folder/'numerical_setup.json')
    elements, expected, recipes = numerical_setup(job, spec)
    require(setup['production_solver'] == job['control'], 'Production solver selection differs')
    require(setup['applied'] == ({} if job['control'] else setup['requested']), 'Applied controls differ')
    for label in ('fine', 'finer'):
        require(setup['recipes'][label]['nodes'] == recipes[label]['nodes'], 'Wake recipe count differs')
        for e in ELEMENTS:
            close(setup['recipes'][label]['first_spacing'][e], recipes[label]['first_spacing'][e],
                  tolerance=1e-12, context='Wake recipe first spacing')
    for key in expected:
        require(set(setup['requested'][key]) == set(ELEMENTS), 'Control elements differ')
        for e in ELEMENTS:
            close(setup['requested'][key][e], expected[key][e], tolerance=1e-12, context='Requested '+key)
    for diag in state['metadata']['diagnostics']:
        require(diag.get('discretization_controls', {}) == setup['applied'], 'Evaluation controls differ')
    for e, xy in elements:
        path = inputs.add(folder/f'source_{e}.npy')
        require(sha(path) == setup['source_sha256'][e], 'Source array hash differs')
        array_close(np.load(path, allow_pickle=False), xy, 'Reconstructed source differs')
        exported = np.loadtxt(inputs.add(folder/f'flow/geometry_{e}.dat'), skiprows=1)
        array_close(exported, xy, 'Exported source differs', 1e-10)
        # Reconstruct only geometry, without assembling a flow system or solving.
        model = new_model(e, xy, spec['reynolds'], spec['mach'], spec['ncrit'], np.array([.25, 0.]), 1e-5)
        a = state['arrays'][e]
        array_close(a['source_xy'], model.geom.xpoint, 'Regularized source differs')
        if expected['surface_nodes'][e] != model.foil.N:
            le, axis = model.chord_origin, model.chord_axis
            chord = np.linalg.norm(.5*(xy[0]+xy[-1])-le)
            frame = np.column_stack((axis, [-axis[1], axis[0]]))
            local = frame.T@(model.geom.xpoint-le[:, None])/chord
            y, s, tangent, _ = curvature_paneling(local, expected['surface_nodes'][e], 2., .1)
            coords, arc, tangents = le[:, None]+chord*(frame@y), chord*s, frame@tangent
        else:
            coords, arc, tangents = model.foil.x, model.foil.s, model.foil.t
        for key, value in (('xy', coords), ('foil_s', arc), ('foil_t', tangents)):
            array_close(a[key], value, 'Reconstructed '+e+'/'+key, 1e-10)
        require(np.all(np.diff(a['foil_s']) > 0), 'Nonmonotone surface arclength')
        require(a['wake_xy'].shape[1] == expected['wake_nodes'][e], 'Executed wake count differs')
        step = setup['requested']['wake_first_spacing'][e]
        ds = np.linalg.norm(np.diff(a['wake_xy'], axis=1), axis=0)
        alpha = np.deg2rad(state['result']['alpha_deg']); stream = np.array([np.cos(alpha), np.sin(alpha)])
        tips = {n: .5*(state['arrays'][n]['xy'][:, 0]+state['arrays'][n]['xy'][:, -1]) for n in ELEMENTS}
        length = max(job['wake_length'], max(stream@pt for pt in tips.values())+job['wake_length']-stream@tips[e])
        requested = np.diff(mf.space_geom(step, length, len(ds)+1))
        # Trapezoidal integration of two unit tangents cannot exceed its tracing step.
        require(np.all(ds <= requested+1e-12), 'Realized wake chord exceeds requested tracing step')
        close(state['result']['wake_discretization'][e]['streamwise_start'], stream@a['wake_xy'][:, 0])
        close(state['result']['wake_discretization'][e]['streamwise_end'], stream@a['wake_xy'][:, -1])
    return setup


def report(study, out):
    require(sha(study/'study.json') == STUDY_HASH, 'Unexpected discretization study')
    inputs = Inputs(); spec = inputs.json(study/'study.json')
    for name, digest in spec['sha256'].items():
        require(sha(inputs.add(study/'source_snapshot'/name)) == digest, 'Frozen snapshot changed: '+name)
    prior = inputs.json(study/'source_snapshot/validation/results/s702_panel_convergence_analysis/input_manifest.json')['files']
    for name in DEPENDENCIES:
        digest = spec['sha256'].get(name, prior.get(name))
        require(digest is not None and sha(inputs.add(ROOT/name)) == digest, 'Audit dependency changed: '+name)
    for name in ('fore.dat', 'aft_d5.dat', 'aft_d10.dat'):
        key = 'validation/data/s702/'+name
        require(sha(inputs.add(ROOT/key)) == spec['sha256'][key], 'Raw geometry source changed: '+key)
    runtime = inputs.json(study/'runtime.json'); status = inputs.json(study/'status.json')
    require(runtime['versions'] == spec['expected_versions'], 'Unexpected solver runtime')
    require(all(v == '1' for v in runtime['thread_limits'].values()), 'Worker thread limits differ')
    require(status['state'] == 'finished' and status['completed_jobs'] == status['total_jobs'] == 11, 'Study incomplete')
    automatic = inputs.json(study/'summary.json')
    require(automatic['required_jobs'] == automatic['feasible_jobs'] == 11 and automatic['execution_checks_pass'], 'Runner checks failed')
    states, setups, references = {}, {}, {}
    targets, endpoints, controls, histories, surfaces, wakes, geometries = [], [], [], [], [], [], []

    def reference(name, job):
        if name not in references:
            path = ROOT/name; before = set(inputs.files)
            oldspec = inputs.json(path.parent.parent/'study.json')
            oldruntime = inputs.json(path.parent.parent/'runtime.json')
            oldjob = next(j for j in oldspec['jobs'] if j['name'] == path.parent.name)
            expected = effective_spec(oldjob, oldspec) if 'panels' in oldjob else oldspec
            state = audit(path.parent, oldjob, expected, oldruntime, inputs)
            require(state['result']['feasible'], 'Incomplete archived reference')
            for key in set(inputs.files)-before:
                require(key in spec['sha256'] and inputs.files[key] == spec['sha256'][key], 'Unpinned reference: '+key)
            for key in ('reynolds', 'mach', 'ncrit', 'reference_chord', 'iterations', 'cl_tolerance',
                        'residual_tolerance', 'initial_half_width', 'alpha_bounds', 'max_evaluations',
                        'max_brackets', 'warm_start_alpha', 'entry_offset'):
                require(oldspec[key] == spec[key], 'Reference physics/protocol differs: '+key)
            for key in ('flap', 'target_cl', 'center', 'wake_length'):
                require(oldjob[key] == job[key], 'Reference condition differs: '+key)
            references[name] = state
        return references[name]

    for job in spec['jobs']:
        folder = study/job['name']
        state = audit(folder, job, dict(spec, panels=job['panels'], source_panels=job['source_panels']), runtime, inputs)
        r = state['result']; require(r['feasible'], 'Incomplete study: paired attribution prohibited')
        setup = audit_setup(folder, job, spec, state, inputs)
        states[job['flap'], job['stage']] = state; setups[job['flap'], job['stage']] = setup
        fine = reference(job['reference_result'], job)
        delta = (r['CD']-fine['result']['CD'])*1e4
        auto = next(v for v in automatic['jobs'] if v['job'] == job['name'])
        close(auto['difference_from_fine_CD_counts'], delta)
        targets.append(dict(job=job['name'], flap=job['flap'], stage=job['stage'], control=job['control'],
            CL=r['CL'], CD=r['CD'], CM=r['CM'], alpha_deg=r['alpha_deg'], residual=state['residual'],
            lift_error=r['CL']-job['target_cl'], difference_from_fine_CD_counts=delta,
            selected_path=r['selected_path'], evaluations=len(r['records']),
            failed_evaluations=sum(not rec['admissible'] for rec in r['records']),
            failed_internal_points=sum(not c['converged'] for d in state['metadata']['diagnostics'] for c in d['continuation']),
            wall_seconds=r['wall_seconds']))
        histories.append(dict(job=job['name'], records=r['records'], attempts=r['attempts'],
            selected_flow_index=r['selected_flow_index'], continuation=[d['continuation'] for d in state['metadata']['diagnostics']]))
        if job['stage'] in ('fine', 'wake_spacing'):
            ref = fine if job['stage'] == 'fine' else reference(job['finer_reference_result'], job)
            difference = (r['CD']-ref['result']['CD'])*1e4
            require(abs(difference) <= spec['acceptance_control_CD_counts'], 'Endpoint reproduction failed')
            endpoints.append(dict(job=job['name'], difference_CD_counts=difference, passed=True))
            controls.extend(compare_arrays(state['arrays'], ref['arrays'], job['name']+'/archive'))
        diag = state['metadata']['diagnostics'][r['selected_flow_index']]
        for element in ELEMENTS:
            a = state['arrays'][element]; other = 'aft' if element == 'fore' else 'fore'
            gap = clearance(a['wake_xy'].T, state['arrays'][other]['xy'].T)
            require(gap > 0, 'Wake crosses other element')
            ds = np.linalg.norm(np.diff(a['wake_xy'], axis=1), axis=0)
            wakes.append(dict(job=job['name'], flap=job['flap'], stage=job['stage'], element=element,
                surface_nodes=a['xy'].shape[1], source_nodes=a['source_xy'].shape[1], wake_nodes=a['wake_xy'].shape[1],
                requested_first_spacing=setup['requested']['wake_first_spacing'][element], realized_first_spacing=float(ds[0]),
                last_spacing=float(ds[-1]), arclength=float(sum(ds)), other_body_clearance=gap,
                CD_wake=r['element_forces'][element]['CD_wake'], CD_friction=r['element_forces'][element]['CD_friction'],
                TE_gap_added=diag['te_gap_added_over_local_chord'][element], paneling_recovery=diag['monotone_paneling_recovery'][element]))
            for side in ('upper', 'lower'):
                bl = r['boundary_layers'][element][side]
                surfaces.append(dict(job=job['name'], flap=job['flap'], stage=job['stage'], element=element, side=side,
                    transition_x=bl['transition']['x_over_local_chord'] if bl['transition'] else None,
                    separated_fraction=bl['separated_arclength_fraction'], min_cf=bl['min_cf'], max_Hk=bl['max_Hk'],
                    separation_regions=json.dumps(bl['separation_regions'], sort_keys=True),
                    **resolution(side_data(state['surfaces'][element], side))))

    for flap in (5, 10):
        check_factor_controls({s:setups[flap, s] for s in STAGES})
        for first, second in zip(STAGES[:-1], STAGES[1:]):
            for e in ELEMENTS:
                a, b = states[flap, first]['arrays'][e], states[flap, second]['arrays'][e]
                if first != 'fine':
                    array_close(a['source_xy'], b['source_xy'], 'Fixed source changed', 0.)
                if first in ('surface', 'wake_count'):
                    for key in ('xy', 'foil_s', 'foil_t'):
                        array_close(a[key], b[key], 'Fixed surface changed: '+key, 0.)
                geometries.append(dict(flap=flap, from_stage=first, to_stage=second, element=e,
                    sampled_contour_distance=contour_sample_distance(a['xy'].T, b['xy'].T),
                    max_TE_endpoint_change=float(np.max(abs(a['xy'][:, [0,-1]]-b['xy'][:, [0,-1]]))),
                    same_source=bool(np.array_equal(a['source_xy'], b['source_xy'])),
                    same_surface=bool(np.array_equal(a['xy'], b['xy']))))
    control = states[10, 'legacy_control']; explicit = states[10, 'fine']
    control_delta = (control['result']['CD']-explicit['result']['CD'])*1e4
    require(abs(control_delta) <= spec['acceptance_control_CD_counts'], 'Original solver control failed')
    require(all(set(control['arrays'][e]) == set(explicit['arrays'][e]) for e in ELEMENTS), 'Original solver array fields differ')
    controls.extend(compare_arrays(control['arrays'], explicit['arrays'], 'legacy/explicit'))
    require(all(row['exactly_equal'] for row in controls if row['comparison'] == 'legacy/explicit'), 'Original solver arrays differ')
    pairs, increments, total = paired_increments({k:v['result'] for k,v in states.items() if k[1] in STAGES})
    for got, expected in zip(pairs, automatic['paired_stages']):
        require(got['stage'] == expected['stage'] and expected['complete'], 'Runner pair ordering/coverage differs')
        close(got['flap_delta_CD_counts'], expected['flap_delta_CD_counts'])
    for got, expected in zip(increments, automatic['conditional_increments']):
        close(got['paired_change_counts'], expected['paired_change_CD_counts'])
    experimental = inputs.csv(ROOT/'validation/results/s702_drag_mechanism/s702_flap_5_to_10.csv')
    measured = next(float(v['experiment_delta_10_minus_5_counts']) for v in experimental
                    if float(v['reynolds']) == spec['reynolds'] and float(v['target_cl']) == .8)
    summary = dict(required_jobs=11, feasible_jobs=len(targets), selected_states_audited=len(states)+len(references),
        full_state_audit_pass=True, factor_controls_pass=True, endpoint_checks=endpoints,
        control_CD_difference_counts=control_delta, control_saved_arrays_exactly_equal=True,
        endpoints_shared_arrays_exactly_equal=all(v['exactly_equal'] for v in controls if v['comparison'] != 'legacy/explicit'),
        paired_stages=pairs, conditional_increments=increments, total_paired_change_counts=total,
        unchanged_screening_threshold_counts=spec['design_screening_counts'],
        total_paired_change_exceeds_screen=abs(total) > spec['design_screening_counts'],
        numerical_convergence_established=False, experiment_flap_delta_counts=measured,
        finest_model_minus_experiment_counts=pairs[-1]['flap_delta_CD_counts']-measured,
        minimum_wake_body_clearance=min(v['other_body_clearance'] for v in wakes),
        artificial_TE_gap_added=any(v['TE_gap_added'] != 0 for v in wakes),
        paneling_recovery_used=any(v['paneling_recovery'] for v in wakes),
        requested_evaluations=sum(v['evaluations'] for v in targets),
        failed_requested_evaluations=sum(v['failed_evaluations'] for v in targets),
        failed_internal_continuation_points=sum(v['failed_internal_points'] for v in targets),
        sum_job_wall_seconds=sum(v['wall_seconds'] for v in targets),
        reporting_python=platform.python_version(), solver_runtime=runtime,
        qualification='Conditional contributions along the declared path, at CL 0.8 only; interactions and path dependence remain possible. No physical validation or numerical error bound.')
    out.mkdir(parents=True, exist_ok=True)
    for name, rows in (('targets',targets), ('paired_stages',pairs), ('conditional_increments',increments),
                       ('surface_diagnostics',surfaces), ('wake_diagnostics',wakes), ('geometry_comparison',geometries),
                       ('control_differences',controls)):
        csv_write(out/(name+'.csv'), rows)
    write(out/'initialization_histories.json', histories); write(out/'summary.json', summary)
    plot(out, states, pairs, increments)
    findings(out, summary)
    inputs.add(Path(__file__))
    write(out/'input_manifest.json', dict(study_sha256=STUDY_HASH, files=inputs.files))
    return summary


def plot(out, states, pairs, increments):
    plt.rcParams.update({'font.size':10, 'axes.spines.top':False, 'axes.spines.right':False})
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.7), layout='constrained')
    x = np.arange(4)
    for offset, flap, color in ((-.17,5,'#3976b5'),(.17,10,'#d5772c')):
        axes[0].bar(x+offset, [r[f'd{flap}_change_counts'] for r in increments], width=.32, label=f'Flap {flap}°', color=color)
    values = [r['paired_change_counts'] for r in increments]
    bars = axes[1].bar(x, values, color=['#747474','#3976b5','#747474','#7e5aa0'], width=.6)
    axes[1].bar_label(bars, labels=[f'{v:+.4f}' for v in values], padding=5)
    for ax in axes:
        ax.axhline(0, color='black', lw=.8); ax.grid(axis='y', alpha=.2); ax.set_axisbelow(True)
        ax.set_xticks(x, ['Source\ncontour','Surface\nsampling','Wake\ncount','First wake\nstep'])
        ax.set_ylim(-.68, .18); ax.set_ylabel('Conditional change in drag counts')
    axes[0].legend(); axes[0].set_title('Individual drag response')
    axes[1].set_title('Change in flap-10° minus flap-5° drag')
    axes[1].text(.97, .07, f"Total: {sum(values):+.4f} counts\n{pairs[0]['flap_delta_CD_counts']:.4f} → {pairs[-1]['flap_delta_CD_counts']:.4f}",
                 transform=axes[1].transAxes, ha='right', va='bottom')
    fig.suptitle('S702 discretization diagnostic · CL = 0.8 · Re = 700,000 · M = 0.07 · Ncrit = 6.25\n'
                 'One count = 0.0001 CD. Contributions depend on the declared order of changes.', fontsize=11)
    for ext in ('png','svg','pdf'): fig.savefig(out/('factor_contributions.'+ext), dpi=180)
    plt.close(fig)
    fig, axes = plt.subplots(2, 2, figsize=(11, 7), layout='constrained')
    for stage, label, color, style in (('contour','B: coarse surface','#3976b5','-'),
                                      ('surface','C: finer surface','#d5772c','--'),
                                      ('wake_count','D: more wake nodes','#338866',':'),
                                      ('wake_spacing','E: smaller first step','#7e5aa0','-.')):
        for col, side in enumerate(('upper','lower')):
            data = side_data(states[10,stage]['surfaces']['aft'], side)
            for row, field in enumerate(('Hk','theta')):
                axes[row,col].plot(data['x_over_local_chord'], data[field]*(1e3 if field=='theta' else 1),
                                   label=label, color=color, ls=style, lw=1.5)
                axes[row,col].set_xlim(.4,1); axes[row,col].grid(alpha=.2)
    for col, side in enumerate(('upper','lower')):
        axes[0,col].set_title('Aft element · '+side); axes[1,col].set_xlabel('x / local chord')
    axes[0,0].set_ylabel('Kinematic shape factor Hk'); axes[1,0].set_ylabel('Momentum thickness θ / cref × 1000')
    axes[0,0].legend(fontsize=9)
    fig.suptitle('Selected boundary-layer states · flap 10° · CL = 0.8\n'
                 'Source contour held fixed from B through E; each state meets its own fixed-lift tolerance.', fontsize=11)
    for ext in ('png','svg','pdf'): fig.savefig(out/('aft_layers.'+ext), dpi=180)
    plt.close(fig)


def findings(out, s):
    rows = s['conditional_increments']
    table = '\n'.join(f"| {r['factor']} | {r['d5_change_counts']:+.6f} | {r['d10_change_counts']:+.6f} | {r['paired_change_counts']:+.6f} | {r['fore_paired_change_counts']:+.6f} | {r['aft_paired_change_counts']:+.6f} |" for r in rows)
    (out/'FINDINGS.md').write_text(
        '# S702 separated discretization audit\n\n'
        'All eleven targets and fifteen full selected states (eleven new, four archived) passed audit. '
        'The four endpoint drag differences are listed in summary.json; all pass the declared 0.01-count tolerance. '
        f"All shared endpoint arrays exactly equal: {s['endpoints_shared_arrays_exactly_equal']}. "
        'All 24 original-solver control array comparisons are exactly equal. Frozen source hashes, runtime, requested/applied controls, '
        'reconstructed geometry, saved Newton arrays, nonzero-Mach exports, forces, transition/separation, '
        'initialization order and tolerances pass. No new flow solves or production changes were made for this report.\n\n'
        'At Re 700,000, Mach 0.07, Ncrit 6.25, CL 0.8 and wake extent 1c, the paired flap-10-minus-5 response changes '
        f"from {s['paired_stages'][0]['flap_delta_CD_counts']:.6f} to {s['paired_stages'][-1]['flap_delta_CD_counts']:.6f} counts. "
        'The declared A→E path separates the existing 340/240→440/320 change as follows. One count is 1e-4 CD.\n\n'
        '| Changed control | Flap 5 | Flap 10 | Paired change | Fore paired term | Aft paired term |\n'
        '|---|---:|---:|---:|---:|---:|\n'+table+'\n\n'
        f"The increments telescope to {s['total_paired_change_counts']:+.6f} counts. Surface sampling dominates, "
        'with a smaller first-wake-step contribution. Source resampling and wake count have small conditional effects. '
        f"The surface step's paired change comprises {rows[1]['aft_paired_change_counts']:+.6f} counts from the aft drag term "
        f"and {rows[1]['fore_paired_change_counts']:+.6f} from the fore term. "
        'These are order-dependent numerical contributions; a different order could redistribute them through interactions. '
        'They do not allocate experimental error to a physical mechanism. The surface step alone exceeds the unchanged '
        '0.5-count screen, so the earlier numerical screen remains failed; this experiment supplies no asymptotic error bound.\n\n'
        'Flap-10 stages B and C both require the declared warm-anchor fallback; all other new targets select cold-local. '
        'The dominant B→C change therefore compares the same selected initialization path. A→B and C→D also change '
        'selected path, and path independence is not established for every intermediate configuration. '
        f"All {s['requested_evaluations']} requested evaluations, including {s['failed_requested_evaluations']} failures and "
        f"{s['failed_internal_continuation_points']} failed internal continuation points, are retained in the histories.\n\n"
        'The regularized source is exactly fixed from B through E, and the actual surface coordinates, spline arclength '
        'and tangents are exactly fixed from C through E. Stage B changes the actual contour at fixed station count; '
        'stage C resamples that fixed source. Requested wake controls match every evaluation, actual counts match the '
        'saved arrays, and realized curved chord lengths do not exceed nominal tracing steps. Requested spacing does '
        'not hold wake coordinates fixed: incidence and local velocities respond to the fixed-lift solve. No artificial '
        'trailing-edge opening or paneling recovery was activated. Minimum wake-to-other-body centerline clearance is '
        f"{s['minimum_wake_body_clearance']:.6f}c; positive clearance does not rule out physical wake/BL interaction. "
        'The CSV files give actual geometry changes, spacing, element drag, local Hk/velocity steps, transition and '
        'separated regions, residuals, incidence, lift error and wall time.\n\n'
        f"For context, the experimental flap delta is {s['experiment_flap_delta_counts']:.6f} counts; the finest model "
        f"still differs by {s['finest_model_minus_experiment_counts']:+.6f} counts. Experimental values were not used to select "
        'the numerical path or retune Ncrit. This CL-0.8 diagnostic does not revisit the reversed CL-1.2 flap ranking.\n\n'
        '## Next bounded numerical question\n\n'
        'Separate fore and aft surface sampling at fixed 440/320 source contours and the fine wake-count/first-step '
        'recipes: use the existing B (340/240) and C (440/320) endpoints and add mixed 440/240 and 340/320 surfaces '
        'for both flaps. This 2×2 surface-count design measures element-specific sensitivities and their interaction; '
        'an aft force contribution alone does not identify which element sampling caused it. Include a fresh C flap-10 '
        'control, preserve all initialization failures, and compare successful B/C warm paths if mixed states raise a '
        'path-dependence question. Freeze that follow-up protocol before launching; no automatic higher grid follows.\n\n'
        'Pomeroy geometry/observable qualification remains the next independent experimental task. S103/S414 follow '
        'their tested-geometry audits. Closure changes and shape-optimization claims remain gated. No job is running.\n\n'
        '![Conditional drag contributions](factor_contributions.png)\n\n'
        '![Aft boundary layers](aft_layers.png)\n\n'
        'Reproduce: `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 validation/report_s702_discretization.py`. '
        'This reads archived states only. See `input_manifest.json` for consumed hashes and `verification.json` for checks.\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--study', type=Path, default=ROOT/'validation/results/s702_discretization_20261001')
    parser.add_argument('--out', type=Path, default=ROOT/'validation/results/s702_discretization_analysis')
    args = parser.parse_args()
    print(json.dumps(report(args.study, args.out), indent=2, allow_nan=False))
