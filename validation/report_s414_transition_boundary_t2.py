#!/usr/bin/env python3
"""Tables, declared screens and figures for V2 tranche 2, combined with audited tranche 1.

Reads the retrieved studies only after the tranche-2 audit is complete. No flow solve.
Same-angle comparisons are numerical path diagnostics, never matched-lift drag scores.
"""
import argparse
import json
from pathlib import Path
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from report_design_initialization import ROOT, Inputs, sha, write, csv_write
from report_s702_drag_mechanism import require
import report_s414_transition_boundary as t1report
from report_s414_transition_boundary import describe, compare, transition_x, style, SURFACE, INK, SECONDARY, INCREASING, DECREASING
import run_s414_transition_boundary as deck
import run_s414_transition_boundary_t2 as t2

T1, STUDY = ROOT/deck.STUDY, ROOT/t2.STUDY
OUT = ROOT/'validation/results/s414_transition_boundary_t2_analysis_20261004'
COLD = '#1baf7a'  # validated categorical slot 3; direct-labelled because its contrast is below 3:1
TARGET_CL = .175
TITLES = {'nominal': 'Nominal reconstruction: archived fallback bracket', 'fore_source_cubic': 'Fore-only cubic reconstruction: refined bracket'}
MEASURES = [('CL', 'Lift coefficient'), ('CD_counts', 'Drag [counts]'), ('fore_lower_transition_x', 'Fore-lower transition, x / chord')]


def figure(out, name, scenario, rows, anchor=None):
    data = [r for r in rows if r['scenario'] == scenario]
    fig, axes = plt.subplots(len(MEASURES), 1, figsize=(8.4, 2.5*len(MEASURES)+.9), sharex=True, layout='constrained')
    fig.patch.set_facecolor(SURFACE); handles = {}
    for ax, (key, label) in zip(axes, MEASURES):
        style(ax)
        for direction, color, marker, size in (('increasing', INCREASING, '>', 11), ('decreasing', DECREASING, '<', 7)):
            series = sorted((r for r in data if r['direction'] == direction), key=lambda r: r['alpha_deg'])
            # A failed point is a gap in the line and a mark on the axis, never a value.
            values = [r[key] if r['converged'] and r.get(key) is not None else np.nan for r in series]
            handles[direction] = ax.plot([r['alpha_deg'] for r in series], values, color=color, linewidth=1.6, marker=marker,
                markersize=size, markeredgecolor=SURFACE, markeredgewidth=1.3, solid_capstyle='round',
                label=f'Reached from {"below" if direction == "increasing" else "above"} ({direction} incidence)')[0]
            for r in series:
                if not r['converged']:
                    handles['failed'] = ax.plot([r['alpha_deg']], [.03], transform=ax.get_xaxis_transform(), linestyle='none', marker='x',
                        markersize=8, markeredgewidth=1.8, color=INK, clip_on=False, label='Not converged (no value)')[0]
        for r in data:
            if r['direction'] == 'cold' and r['converged'] and r.get(key) is not None:
                handles['cold'] = ax.plot([r['alpha_deg']], [r[key]], linestyle='none', marker='D', markersize=7.5, color=COLD,
                    markeredgecolor=SURFACE, markeredgewidth=1.3, label='Cold start at that incidence', zorder=5)[0]
                if key == 'CL':
                    ax.annotate('cold start', (r['alpha_deg'], r[key]), xytext=(0, 11), textcoords='offset points', ha='center',
                                fontsize=8.5, color=SECONDARY)
            if r['direction'] == 'cold' and not r['converged']:
                handles['failed'] = ax.plot([r['alpha_deg']], [.03], transform=ax.get_xaxis_transform(), linestyle='none', marker='x',
                    markersize=8, markeredgewidth=1.8, color=INK, clip_on=False, label='Not converged (no value)')[0]
        if anchor and anchor.get(key) is not None:
            handles['anchor'] = ax.plot([anchor['alpha_deg']], [anchor[key]], linestyle='none', marker='o', markersize=8, color=INK,
                markeredgecolor=SURFACE, markeredgewidth=1.3, label='Archived audited state', zorder=4)[0]
        if key == 'CL' and scenario == 'fore_source_cubic':
            ax.axhline(TARGET_CL, color=SECONDARY, linewidth=.8); ax.annotate('CL 0.175', (.995, TARGET_CL), xycoords=ax.get_yaxis_transform(),
                xytext=(0, 3), textcoords='offset points', ha='right', fontsize=8.5, color=SECONDARY)
        ax.set_ylabel(label)
    axes[0].set_title(TITLES[scenario], fontsize=10.5, color=INK, loc='left'); axes[-1].set_xlabel('Incidence [deg]')
    fig.legend(list(handles.values()), [h.get_label() for h in handles.values()], frameon=False, fontsize=8.5,
               loc='outside lower center', ncols=2, labelcolor=SECONDARY)
    for suffix in ('png', 'pdf', 'svg'): fig.savefig(out/f'{name}.{suffix}', dpi=170, facecolor=SURFACE)
    plt.close(fig)


def shift(first, second):
    a, b = (transition_x(r['boundary_layers']['fore']['lower']) for r in (first, second))
    return float('inf') if a is None or b is None else abs(b-a)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=OUT)
    out = parser.parse_args().out.resolve()
    audit = json.loads((out/'residual_recheck.json').read_text())
    require(audit['complete'] and len(audit['states']) == audit['saved_states'], 'Independent audit is not complete')
    inputs = Inputs(); spec = inputs.json(STUDY/'study.json'); spec1 = inputs.json(T1/'study.json'); screens = spec['equivalence_screens']
    for name, digest in audit['input_sha256'].items():
        if name.startswith(t2.STUDY+'/') and 'source_snapshot' not in name:
            require(sha(ROOT/name) == digest, 'Audited output changed: '+name)
    rows, points = [], {}
    def add(job, rec, folder, source):
        row = describe(job, rec, folder); row['source'] = source; rows.append(row)
        points[(job['scenario'], job['direction'], rec['requested_alpha_deg'], source)] = dict(record=rec, folder=folder)
    for job in spec1['jobs']:
        if job['control']: continue
        for rec in inputs.json(T1/job['name']/'result.json')['records']:
            if rec['role'] == 'path': add(job, rec, T1/job['name'], 'tranche 1 live path')
    restore_rows = []
    for job in spec['jobs']:
        result = inputs.json(STUDY/job['name']/'result.json'); rec = result['records'][0]
        if job['kind'] == 'restore_control':
            restore_rows.append(dict(job=job['name'], alpha_deg=job['alpha_deg'], **{k: v for k, v in audit['controls']['jobs'][job['name']].items() if k != 'reference'},
                restored_residual=max(result['restore']['BL_residual_l2'], result['restore']['coupling_residual_max']), iterations=rec['iterations_final_step']))
            continue
        add(dict(job, direction='cold' if job['kind'] == 'cold_check' else job['direction']), rec, STUDY/job['name'],
            'tranche 2 cold check' if job['kind'] == 'cold_check' else 'tranche 2 restored continuation')
    for r in rows:
        if r['direction'] not in ('increasing', 'decreasing', 'cold'): r['direction'] = 'cold'
    keys = list(dict.fromkeys(k for r in rows for k in r)); csv_write(out/'points.csv', [{k: r.get(k) for k in keys} for r in rows])
    csv_write(out/'restore_controls.csv', restore_rows)
    get = lambda scenario, direction, alpha: next((p for (s, d, a, _), p in points.items()
        if s == scenario and d == direction and a == alpha and p['record']['converged']), None)
    directions, colds, response, lift, multiplicity = [], [], [], [], []
    for scenario in deck.ARCHIVED_JOBS:
        new = spec['cubic_new_incidences_deg'] if scenario == 'fore_source_cubic' else spec['nominal_new_incidences_deg']
        for alpha in new:
            a, b = get(scenario, 'increasing', alpha), get(scenario, 'decreasing', alpha)
            row = dict(scenario=scenario, alpha_deg=alpha, status='both_converged' if a and b else 'neither_converged' if not a and not b
                       else ('only_increasing_converged' if a else 'only_decreasing_converged'))
            row.update(compare(a, b, screens) if a and b else dict(within_screens=None)); directions.append(row)
        cold_alpha = next(j['alpha_deg'] for j in spec['jobs'] if j['kind'] == 'cold_check' and j['scenario'] == scenario)
        cold = get(scenario, 'cold', cold_alpha)
        for direction in ('increasing', 'decreasing'):
            other = get(scenario, direction, cold_alpha)
            row = dict(scenario=scenario, alpha_deg=cold_alpha, comparison=f'cold start minus state reached from {"below" if direction == "increasing" else "above"}',
                       status='both_converged' if cold and other else 'cold_failed' if not cold else 'path_state_missing')
            row.update(compare(cold, other, screens) if cold and other else dict(within_screens=None)); colds.append(row)
        angles = sorted({a for (s, _, a, _) in points if s == scenario})
        merged = []
        for alpha in angles:
            found = {d: get(scenario, d, alpha) for d in ('increasing', 'decreasing', 'cold')}
            states = {d: p for d, p in found.items() if p}
            entry = dict(scenario=scenario, alpha_deg=alpha, converged_from=' '.join(sorted(states)) or 'none')
            names = list(states); distinct = False
            for i, a in enumerate(names):
                for b in names[i+1:]:
                    if not compare(states[a], states[b], screens)['within_screens']:
                        distinct = True; multiplicity.append(dict(scenario=scenario, alpha_deg=alpha, first=a, second=b, **compare(states[a], states[b], screens)))
            entry['distinct_converged_states'] = distinct
            for d in ('increasing', 'decreasing', 'cold'):
                rec = states[d]['record'] if d in states else None
                entry[d+'_CL'] = rec['CL'] if rec else None; entry[d+'_CD_counts'] = rec['CD']*1e4 if rec else None
                entry[d+'_fore_lower_transition_x'] = transition_x(rec['boundary_layers']['fore']['lower']) if rec else None
            response.append(entry)
            if states and not distinct: merged.append((alpha, next(iter(states.values()))['record']))
        for (a0, r0), (a1, r1) in zip(merged, merged[1:]):
            lift.append(dict(scenario=scenario, lower_alpha_deg=a0, upper_alpha_deg=a1, lower_CL=r0['CL'], upper_CL=r1['CL'],
                delta_CL=r1['CL']-r0['CL'], lift_increases=bool(r1['CL'] > r0['CL']),
                brackets_target_CL=bool((r0['CL']-TARGET_CL)*(r1['CL']-TARGET_CL) < 0),
                lower_CD_counts=r0['CD']*1e4, upper_CD_counts=r1['CD']*1e4,
                fore_lower_transition_change=shift(r0, r1)))
    uniform = lambda data: [{k: r.get(k) for k in dict.fromkeys(k for r in data for k in r)} for r in data]
    for name, data in (('direction_pairs.csv', directions), ('cold_checks.csv', colds), ('combined_response.csv', response),
                       ('lift_response.csv', lift)):
        csv_write(out/name, uniform(data))
    if multiplicity: csv_write(out/'distinct_states.csv', uniform(multiplicity))
    anchor = inputs.json(ROOT/t2.ANCHOR/t2.ANCHOR_JOB/'result.json')
    anchor_row = dict(alpha_deg=anchor['alpha_deg'], CL=anchor['CL'], CD_counts=anchor['CD']*1e4,
                      fore_lower_transition_x=transition_x(anchor['boundary_layers']['fore']['lower']))
    figure(out, 'cubic_refined_response', 'fore_source_cubic', rows)
    figure(out, 'nominal_fallback_response', 'nominal', [r for r in rows if r['source'] != 'tranche 1 live path'
           or r['alpha_deg'] == spec1['angles_deg'][0]], anchor_row)
    both = [r for r in directions if r['status'] == 'both_converged']; new_rows = [r for r in rows if r['source'].startswith('tranche 2')]
    # Resolution counts converged states only: an incidence where nothing converged resolves nothing.
    spacing = {s: float(np.min(np.diff(sorted({r['alpha_deg'] for r in rows if r['scenario'] == s and r['converged']}))))
               for s in deck.ARCHIVED_JOBS}
    widest = {s: max(lift, key=lambda r: r['fore_lower_transition_change'] if r['scenario'] == s else -1) for s in deck.ARCHIVED_JOBS}
    summary = dict(complete=True, audited_saved_states=audit['saved_states'], maximum_fresh_residual=max(s['max_residual'] for s in audit['states']),
        restore_controls=restore_rows, restore_controls_passed=audit['controls']['passed'],
        requested_new_evaluations=len(new_rows), converged_new_evaluations=sum(r['converged'] for r in new_rows),
        failed_new_evaluations=[dict(job=r['job'], direction=r['direction'], alpha_deg=r['alpha_deg'], reached_by=r['reached_by'],
                                     failure_reason=r['failure_reason']) for r in new_rows if not r['converged']],
        direction_pairs=len(directions), direction_pairs_both_converged=len(both),
        direction_pairs_within_screens=sum(bool(r['within_screens']) for r in both),
        largest_direction_difference=(dict(CL=max(abs(r['delta_CL']) for r in both), CD_counts=max(abs(r['delta_CD_counts']) for r in both),
            Cp=max(r['max_abs_delta_Cp'] for r in both), transition_x=max(r['max_transition_shift_x'] for r in both)) if both else None),
        cold_checks=colds, distinct_converged_states_at_one_incidence=multiplicity,
        nonmonotone_lift_intervals=[r for r in lift if not r['lift_increases']],
        intervals_bracketing_target_CL={s: [[r['lower_alpha_deg'], r['upper_alpha_deg'], r['lower_CD_counts'], r['upper_CD_counts']]
                                            for r in lift if r['scenario'] == s and r['brackets_target_CL']] for s in deck.ARCHIVED_JOBS},
        largest_transition_change_between_adjacent_samples={s: max((r for r in lift if r['scenario'] == s),
            key=lambda r: r['fore_lower_transition_change']) for s in deck.ARCHIVED_JOBS},
        finest_converged_angular_spacing_deg=spacing,
        unresolved_transition_interval_deg={s: [widest[s]['lower_alpha_deg'], widest[s]['upper_alpha_deg'],
                                                widest[s]['upper_alpha_deg']-widest[s]['lower_alpha_deg']] for s in deck.ARCHIVED_JOBS}, total_requested_points_V2=spec['requested_points']+spec['tranche1_requested_points'],
        allowance=spec['allowance'], equivalence_screens=screens, new_flow_solves=0,
        qualification='Sampled-path numerical diagnostic at fixed incidence; not matched-lift validation, not physical hysteresis, '
                      'not solution uniqueness. V2 ends with this tranche')
    write(out/'summary.json', summary); write(out/'report_input_manifest.json', dict(files=inputs.files))
    print(json.dumps({k: v for k, v in summary.items() if k not in ('cold_checks', 'restore_controls')}, indent=1))


if __name__ == '__main__':
    main()
