#!/usr/bin/env python3
"""Audit the frozen three-grid S702 study and apply its declared screen; no solves."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from report_s702_wake_extent import (
    ROOT, Inputs, sha, write, csv_write, audit, require, close, clearance, side_data,
)
from report_s702_holdout_panels import contour_sample_distance, resolution
from run_s702_holdout_panels import effective_spec

STUDY_HASH = 'b522cead604443e3470eb7cb38df9efdd1863ea15a6bf8cc6440cbf6eeb43833'
GRIDS = ('base', 'fine', 'finer')
LABELS = ('260/180', '340/240', '440/320')
ELEMENTS = ('fore', 'aft')
DEPENDENCIES = (
    'validation/report_s702_wake_extent.py', 'validation/report_s702_holdout_panels.py',
    'validation/report_design_initialization.py', 'validation/report_s702_drag_mechanism.py',
    'validation/run_s702_holdout_panels.py', 'viscous.py', 'multifoil.py', 'vendor/mfoil.py',
)


def step_metrics(first, second):
    return dict(first_step_counts=float(first), second_step_counts=float(second),
                contracts=bool(abs(second) < abs(first)),
                reverses_direction=bool(first * second < 0),
                absolute_step_ratio=float(abs(second / first)) if first else None)


def assess_contraction(individual, paired, *, coverage, control, audited, threshold=.5):
    """No omitted/failed condition may produce a passing numerical screen."""
    result = dict(decision_eligible=bool(coverage and control and audited),
                  screening_threshold_counts=threshold, numerical_screen_pass=False)
    if not result['decision_eligible']:
        return result
    for name, steps, expected in (('individual', individual, 4), ('paired', paired, 2)):
        values = np.asarray(steps, dtype=float)
        require(values.shape == (expected, 2) and np.isfinite(values).all(),
                'Complete finite '+name+' steps required')
        maxima = np.max(abs(values), axis=0)
        result[name] = dict(first_max_counts=float(maxima[0]), second_max_counts=float(maxima[1]),
                            maximum_contracts=bool(maxima[1] < maxima[0]),
                            below_threshold=bool(maxima[1] < threshold))
    result['numerical_screen_pass'] = all(
        result[name]['maximum_contracts'] and result[name]['below_threshold']
        for name in ('individual', 'paired'))
    return result


def report(study, out):
    require(sha(study/'study.json') == STUDY_HASH, 'Unexpected third-grid study')
    inputs = Inputs(); spec = inputs.json(study/'study.json')
    runtime = inputs.json(study/'runtime.json')
    for name, digest in spec['sha256'].items():
        require(sha(inputs.add(study/'source_snapshot'/name)) == digest, 'Frozen source changed: '+name)
    # Older audit helpers are pinned transitively by the previous audit receipt.
    prior_audit = inputs.json(study/'source_snapshot/validation/results/s702_holdout_panels_analysis/input_manifest.json')['files']
    for name in DEPENDENCIES:
        digest = spec['sha256'].get(name, prior_audit.get(name))
        require(digest is not None and sha(ROOT/name) == digest, 'Audit dependency changed: '+name)
    require(runtime['versions'] == spec['expected_versions'], 'Unexpected numerical runtime')
    require(inputs.json(study/'status.json')['state'] == 'finished', 'Study unfinished')
    automatic = inputs.json(study/'summary.json')
    states, references, targets, controls, histories = {}, {}, [], [], []

    def reference(name, job):
        if name not in references:
            path = ROOT/name; before = set(inputs.files)
            oldspec = inputs.json(path.parent.parent/'study.json')
            oldruntime = inputs.json(path.parent.parent/'runtime.json')
            oldjob = next(j for j in oldspec['jobs'] if j['name'] == path.parent.name)
            expected = effective_spec(oldjob, oldspec) if 'panels' in oldjob else oldspec
            state = audit(path.parent, oldjob, expected, oldruntime, inputs)
            require(state['result']['feasible'], 'Incomplete archived grid coverage')
            for key in set(inputs.files)-before:
                require(key in spec['sha256'] and inputs.files[key] == spec['sha256'][key],
                        'Reference differs from frozen inputs: '+key)
            for key in ('reynolds', 'mach', 'ncrit', 'reference_chord', 'iterations', 'cl_tolerance',
                        'residual_tolerance', 'initial_half_width', 'alpha_bounds', 'max_evaluations',
                        'max_brackets', 'warm_start_alpha'):
                require(oldspec[key] == spec[key], 'Grid physics/protocol differs: '+key)
            for key in ('flap', 'target_cl', 'center'):
                require(oldjob[key] == job[key], 'Reference condition differs: '+key)
            for filename in ('fore.dat', f"aft_d{job['flap']}.dat"):
                key = 'validation/data/s702/'+filename
                require(oldspec['sha256'][key] == spec['sha256'][key], 'Raw geometry differs')
            references[name] = state
        return references[name]

    for job in spec['jobs']:
        current = audit(study/job['name'], job, effective_spec(job, spec), runtime, inputs)
        r = current['result']
        require(r['feasible'], 'Incomplete coverage: paired conclusions prohibited; inspect runner summary')
        fine = reference(job['reference_result'], job)
        delta = (r['CD']-fine['result']['CD'])*1e4
        auto = next(v for v in automatic['jobs'] if v['job'] == job['name'])
        close(auto['difference_from_fine_CD_counts'], delta)
        targets.append(dict(job=job['name'], control=job['control'], flap=job['flap'], target_cl=job['target_cl'],
            panels='/'.join(map(str, job['panels'])), difference_from_fine_CD_counts=delta,
            selected_path=r['selected_path'], evaluations=len(r['records']),
            failed_evaluations=sum(not v['admissible'] for v in r['records']),
            failed_internal_points=sum(not c['converged'] for d in current['metadata']['diagnostics'] for c in d['continuation']),
            wall_seconds=r['wall_seconds']))
        histories.append(dict(job=job['name'], records=r['records'], attempts=r['attempts'],
            selected_flow_index=r['selected_flow_index'],
            continuation=[d['continuation'] for d in current['metadata']['diagnostics']]))
        if job['control']:
            require(abs(delta) <= spec['acceptance_control_CD_counts'], 'Control drag tolerance failed')
            for e in ELEMENTS:
                require(set(fine['arrays'][e]) == set(current['arrays'][e]), 'Control array fields differ')
                for key, a in fine['arrays'][e].items():
                    b = current['arrays'][e][key]
                    controls.append(dict(element=e, field=key, exactly_equal=bool(np.array_equal(a,b)),
                                         max_abs_difference=float(np.max(abs(a-b)))))
        else:
            states[job['flap'], job['target_cl'], 'base'] = reference(job['coarse_reference_result'], job)
            states[job['flap'], job['target_cl'], 'fine'] = fine
            states[job['flap'], job['target_cl'], 'finer'] = current
    require(len(targets) == automatic['required_jobs'] == automatic['feasible_jobs'] == 5
            and automatic['execution_checks_pass'], 'Coverage mismatch')
    require(len(states) == 12 and len(references) == 8, 'Incomplete three-grid matrix')

    forces, surfaces, wakes, geometry, steps = [], [], [], [], []
    for (flap, cl, grid), state in states.items():
        r = state['result']; diag = state['metadata']['diagnostics'][r['selected_flow_index']]
        row = dict(flap=flap, target_cl=cl, grid=grid, panels=LABELS[GRIDS.index(grid)],
                   **{k:r[k] for k in ('CL', 'CD', 'CM', 'alpha_deg', 'selected_path', 'wall_seconds')},
                   lift_error=r['CL']-cl, residual=state['residual'])
        for e in ELEMENTS:
            for field in ('CL', 'CD_wake', 'CM'):
                row[e+'_'+field] = r['element_forces'][e][field]
            a = state['arrays'][e]; w = state['wakes'][e]; other = 'aft' if e == 'fore' else 'fore'
            gap = clearance(a['wake_xy'].T, state['arrays'][other]['xy'].T)
            require(gap > 0, 'Wake/body intersection')
            wakes.append(dict(flap=flap, target_cl=cl, grid=grid, element=e, nodes=len(w),
                arclength=state['stations'][e]['s'][-1],
                first_spacing=float(np.linalg.norm(a['wake_xy'][:,1]-a['wake_xy'][:,0])),
                last_spacing=float(np.linalg.norm(a['wake_xy'][:,-1]-a['wake_xy'][:,-2])),
                other_body_clearance=gap, CD_wake=r['element_forces'][e]['CD_wake'],
                terminal_theta=float(w['theta'][-1]), terminal_delta_star=float(w['delta_star'][-1]),
                terminal_Hk=float(w['Hk'][-1]), terminal_ue=float(w['ue_over_vinf'][-1]),
                TE_gap_added=diag['te_gap_added_over_local_chord'][e],
                paneling_recovery=diag['monotone_paneling_recovery'][e]))
            for side in ('upper', 'lower'):
                values = side_data(state['surfaces'][e], side)
                layer = r['boundary_layers'][e][side]; tr = layer['transition']
                sr = dict(flap=flap, target_cl=cl, grid=grid, element=e, side=side,
                    transition_x=tr['x_over_local_chord'] if tr else None,
                    separated_fraction=layer['separated_arclength_fraction'],
                    separation_regions=len(layer['separation_regions']), max_Hk=layer['max_Hk'],
                    min_cf=layer['min_cf'], **resolution(values))
                for field in ('Cp', 'Hk', 'theta', 'delta_star', 'ue_over_vinf', 'cf'):
                    sr['TE_'+field] = float(values[field][-1])
                surfaces.append(sr)
        forces.append(row)
    for flap in (5, 10):
        for cl in (.8, 1.2):
            rr = [states[flap,cl,g]['result'] for g in GRIDS]
            changes = np.diff([r['CD'] for r in rr])*1e4
            row = dict(flap=flap, target_cl=cl, **step_metrics(*changes))
            for e in ELEMENTS:
                parts = np.diff([r['element_forces'][e]['CD_wake'] for r in rr])*1e4
                row[e+'_first_step_counts'], row[e+'_second_step_counts'] = map(float, parts)
            steps.append(row)
            for g0, g1 in zip(GRIDS[:-1], GRIDS[1:]):
                for e in ELEMENTS:
                    a, b = [states[flap,cl,g]['arrays'][e]['xy'].T for g in (g0,g1)]
                    geometry.append(dict(flap=flap, target_cl=cl, from_grid=g0, to_grid=g1, element=e,
                        from_nodes=len(a), to_nodes=len(b),
                        max_vertex_to_other_polyline_distance=contour_sample_distance(a,b),
                        from_TE_gap=float(np.linalg.norm(a[0]-a[-1])), to_TE_gap=float(np.linalg.norm(b[0]-b[-1])),
                        max_TE_endpoint_change=float(np.max(abs(a[[0,-1]]-b[[0,-1]])))))

    experiment_path = 'validation/results/s702_drag_mechanism/s702_flap_5_to_10.csv'
    require(sha(ROOT/experiment_path) == prior_audit[experiment_path], 'Experimental context changed')
    experimental = inputs.csv(ROOT/experiment_path)
    pairs = []
    for cl in (.8, 1.2):
        ref = next(v for v in experimental if float(v['reynolds']) == spec['reynolds']
                   and float(v['mach']) == spec['mach'] and float(v['target_cl']) == cl)
        row = dict(target_cl=cl, experiment_flap_delta_counts=float(ref['experiment_delta_10_minus_5_counts']))
        for grid in GRIDS:
            a,b = [states[f,cl,grid]['result'] for f in (5,10)]
            row[grid+'_flap_delta_CD_counts'] = (b['CD']-a['CD'])*1e4
            row[grid+'_flap_delta_error_counts'] = row[grid+'_flap_delta_CD_counts']-row['experiment_flap_delta_counts']
            for e in ELEMENTS:
                row[grid+'_'+e+'_flap_delta_CD_counts'] = (b['element_forces'][e]['CD_wake']-a['element_forces'][e]['CD_wake'])*1e4
        row.update(step_metrics(*np.diff([row[g+'_flap_delta_CD_counts'] for g in GRIDS])))
        pairs.append(row)
    decision = assess_contraction(
        [[v['first_step_counts'],v['second_step_counts']] for v in steps],
        [[v['first_step_counts'],v['second_step_counts']] for v in pairs],
        coverage=True, control=automatic['control_pass'], audited=True, threshold=spec['design_screening_counts'])
    summary = dict(required_jobs=5, feasible_jobs=5, selected_states_audited=13, full_state_audit_pass=True,
        control_CD_difference_counts=next(v['difference_from_fine_CD_counts'] for v in targets if v['control']),
        control_saved_arrays_exactly_equal=all(v['exactly_equal'] for v in controls), decision=decision,
        nonmonotone_conditions=[dict(flap=v['flap'], target_cl=v['target_cl']) for v in steps if v['reverses_direction']],
        noncontracting_conditions=[dict(flap=v['flap'], target_cl=v['target_cl']) for v in steps if not v['contracts']],
        minimum_wake_body_clearance=min(v['other_body_clearance'] for v in wakes),
        artificial_TE_gap_added=any(v['TE_gap_added'] != 0 for v in wakes),
        paneling_recovery_used=any(v['paneling_recovery'] for v in wakes),
        maximum_absolute_lift_error=max(abs(v['lift_error']) for v in forces),
        maximum_selected_residual=max(v['residual'] for v in forces),
        transition_presence_and_separation_counts_unchanged=all(
            len({(v['transition_x'] is not None, v['separation_regions']) for v in surfaces
                 if (v['flap'],v['target_cl'],v['element'],v['side']) == (flap,cl,e,side)}) == 1
            for flap in (5,10) for cl in (.8,1.2) for e in ELEMENTS for side in ('upper','lower')),
        requested_evaluations=sum(v['evaluations'] for v in targets),
        failed_requested_evaluations=sum(v['failed_evaluations'] for v in targets),
        failed_internal_continuation_points=sum(v['failed_internal_points'] for v in targets),
        qualification='Operational refinement changes contours and derived wake discretization; no asymptotic error bound or independently validated design ranking.')
    for grid in GRIDS:
        group = [v for v in surfaces if v['grid'] == grid]
        summary[grid+'_max_post_transition_Hk_drop'] = max(v['max_Hk_drop'] for v in group)
        summary[grid+'_max_post_transition_velocity_step'] = max(v['max_velocity_step'] for v in group)
    for grid in GRIDS[1:]:
        group = [v for v in geometry if v['to_grid'] == grid]
        summary[grid+'_max_sampled_geometry_change'] = max(v['max_vertex_to_other_polyline_distance'] for v in group)
        summary[grid+'_max_TE_endpoint_change'] = max(v['max_TE_endpoint_change'] for v in group)
    out.mkdir(parents=True, exist_ok=True)
    for name, rows in (('targets',targets), ('forces',forces), ('individual_steps',steps), ('flap_deltas',pairs),
                       ('surface_diagnostics',surfaces), ('wake_diagnostics',wakes), ('geometry_comparison',geometry),
                       ('control_differences',controls)):
        csv_write(out/(name+'.csv'), rows)
    write(out/'summary.json', summary); write(out/'initialization_histories.json', histories)
    plots(out, states, pairs)
    for name in DEPENDENCIES + ('validation/report_s702_panel_convergence.py', 'test_s702_panel_audit.py'):
        inputs.add(ROOT/name)
    write(out/'input_manifest.json', dict(study_sha256=STUDY_HASH, files=inputs.files))
    return summary


def plots(out, states, pairs):
    def save(fig, name):
        fig.tight_layout()
        for ext in ('png', 'svg'):
            fig.savefig(out/(name+'.'+ext), dpi=170)
        plt.close(fig)
    fig, axes = plt.subplots(1,2,figsize=(11.5,4.7))
    for flap in (5,10):
        for cl in (.8,1.2):
            values = np.array([states[flap,cl,g]['result']['CD'] for g in GRIDS])*1e4
            axes[0].plot(LABELS,values-values[0], 'o-', label=f'Flap {flap}°, CL {cl:g}')
    axes[0].set(ylabel='CD − own 260/180 value [counts]', xlabel='Fore / aft panels', title='Grid sensitivity has not settled')
    for p,color in zip(pairs, ('#277da8','#cf661f')):
        axes[1].plot(LABELS,[p[g+'_flap_delta_CD_counts'] for g in GRIDS], 'o-', color=color, label=f"CL {p['target_cl']:g}, model")
        axes[1].axhline(p['experiment_flap_delta_counts'],color=color,ls=':',label=f"CL {p['target_cl']:g}, experiment")
    axes[1].set(ylabel='Flap 10° − 5° CD [counts]', xlabel='Fore / aft panels', title='Large model–measurement differences remain')
    for ax in axes:
        ax.axhline(0,color='0.55',lw=.6); ax.grid(alpha=.2); ax.legend(fontsize=8)
    save(fig, 'three_grid_drag')
    fig, axes = plt.subplots(2,2,figsize=(11.5,7))
    for col,side in enumerate(('upper','lower')):
        for cl,color in ((.8,'#277da8'), (1.2,'#cf661f')):
            for grid,label,style in zip(GRIDS,LABELS,('-o','--.','-.*')):
                v = side_data(states[10,cl,grid]['surfaces']['aft'],side)
                v = v[v['turbulent']==1]
                for row,field in enumerate(('Hk','ue_over_vinf')):
                    axes[row,col].plot(v['x_over_local_chord'],v[field],style,color=color,ms=3,label=f'CL {cl:g}, {label}')
        axes[0,col].set_title('Flap 10° aft '+side+' after transition')
        axes[1,col].set_xlabel('x / aft local chord')
    axes[0,0].set_ylabel('Hk'); axes[1,0].set_ylabel('Ue / Uinf'); axes[0,1].legend(fontsize=7)
    for ax in axes.flat: ax.grid(alpha=.2)
    save(fig, 'aft_layers')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--study', type=Path, default=ROOT/'validation/results/s702_panel_convergence_20260922')
    parser.add_argument('--out', type=Path, default=ROOT/'validation/results/s702_panel_convergence_analysis')
    args = parser.parse_args()
    print(json.dumps(report(args.study.resolve(), args.out.resolve()), indent=2))
