#!/usr/bin/env python3
"""Tables, declared screens and figures for the V4 convex-nose viscous check.

Reads the retrieved study only after its independent audit is complete, plus
audited archived results for the nominal and cubic reconstructions. No flow solve.
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
from report_s414_transition_boundary import describe, compare, transition_x, style, surface, SURFACE, INK, SECONDARY, MUTED
import run_s414_nose_viscous as run

STUDY = ROOT/run.STUDY
OUT = ROOT/'validation/results/s414_nose_viscous_analysis_20261004'
LOWER = ROOT/'validation/results/s414_lower_bucket_20261003'
JACOBIAN = ROOT/run.ARCHIVE
FORE = ROOT/'validation/results/s414_fore_reconstruction_20261004'
REFERENCES = ROOT/'validation/results/s414_lower_bucket_analysis_20261003/targets.csv'
TAPS = ROOT/'validation/results/s414_lower_bucket_analysis_20261003/pressure_at_taps.csv'
COMBINED = ROOT/'validation/results/s414_transition_boundary_t2_analysis_20261004/combined_response.csv'
COLORS = {'nominal': '#2a78d6', 'convex_nose': '#eb6834', 'fore_source_cubic': '#1baf7a'}  # validated slots 1, 2, 3
NAMES = {'nominal': 'Nominal (PCHIP)', 'convex_nose': 'Convex nose (reversal removed)', 'fore_source_cubic': 'Fore-only cubic'}
CASES = {.056: 'cl0056', .175: 'cl0175', .284: 'cl028'}


def fixed(folder):
    r = json.loads((folder/'result.json').read_text())
    return None if not r['feasible'] else dict(alpha_deg=r['alpha_deg'], CL=r['CL'], CD_counts=r['CD']*1e4,
        fore_CD_counts=r['element_forces']['fore']['CD_wake']*1e4, aft_CD_counts=r['element_forces']['aft']['CD_wake']*1e4,
        fore_lower_transition_x=transition_x(r['boundary_layers']['fore']['lower']), selected_path=r['selected_path'],
        fore_lower_max_Hk=r['boundary_layers']['fore']['lower']['max_Hk'], index=r['selected_flow_index'], folder=folder)


def tap_cp(state, x):
    # Geometric lower surface (after the smallest-x node), as in the archived tap comparison.
    values = surface(state['folder'], state['index'], 'fore'); lower = values[int(np.argmin(values['x_over_cref'])):]
    return float(np.interp(x, lower['x_over_cref'], lower['Cp']))


def figures(out, response, targets, references):
    fig, axes = plt.subplots(3, 1, figsize=(8.4, 8.4), sharex=True, layout='constrained'); fig.patch.set_facecolor(SURFACE)
    measures = [('CL', 'Lift coefficient'), ('CD_counts', 'Drag [counts]'), ('fore_lower_transition_x', 'Fore-lower transition, x / chord')]
    handles = {}
    for ax, (key, label) in zip(axes, measures):
        style(ax)
        for scenario in ('nominal', 'fore_source_cubic', 'convex_nose'):
            data = sorted((r for r in response if r['scenario'] == scenario), key=lambda r: r['alpha_deg'])
            # A point where no route converged is a gap, never a value.
            handles[scenario] = ax.plot([r['alpha_deg'] for r in data], [r[key] if r[key] is not None else np.nan for r in data],
                color=COLORS[scenario], linewidth=1.6, marker='o', markersize=6.5, markeredgecolor=SURFACE, markeredgewidth=1.2,
                solid_capstyle='round', label=NAMES[scenario])[0]
        ax.set_ylabel(label)
    last = max((r for r in response if r['scenario'] == 'convex_nose'), key=lambda r: r['alpha_deg'])
    first = min((r for r in response if r['scenario'] == 'convex_nose'), key=lambda r: r['alpha_deg'])
    axes[1].annotate('convex nose: transition stays at the nose', (first['alpha_deg'], first['CD_counts']), xytext=(4, -22),
                     textcoords='offset points', ha='left', fontsize=8.5, color=SECONDARY)
    axes[-1].set_xlabel('Incidence [deg]')
    axes[0].set_title('Three reconstructions of the fore nose at the same incidences', fontsize=10.5, color=INK, loc='left')
    fig.legend(list(handles.values()), [h.get_label() for h in handles.values()], frameon=False, fontsize=8.5,
               loc='outside lower center', ncols=3, labelcolor=SECONDARY)
    for suffix in ('png', 'pdf', 'svg'): fig.savefig(out/f'incidence_response.{suffix}', dpi=170, facecolor=SURFACE)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(8.4, 5.2), layout='constrained'); fig.patch.set_facecolor(SURFACE); style(ax); handles = {}
    for scenario in ('nominal', 'fore_source_cubic', 'convex_nose'):
        data = sorted((r for r in response if r['scenario'] == scenario and r['CL'] is not None), key=lambda r: r['alpha_deg'])
        ax.plot([r['CL'] for r in data], [r['CD_counts'] for r in data], color=COLORS[scenario], linewidth=1.3, solid_capstyle='round')
        marks = [t for t in targets if t['geometry'] == scenario and t['CD_counts'] is not None and t['admitted']]
        handles[scenario] = ax.plot([t['CL'] for t in marks], [t['CD_counts'] for t in marks], linestyle='none', marker='o', markersize=8.5,
            color=COLORS[scenario], markeredgecolor=SURFACE, markeredgewidth=1.4, label=NAMES[scenario]+', matched-lift state')[0]
    lifts = sorted(references)
    handles['experiment'] = ax.plot(lifts, [references[c]['experiment'] for c in lifts], linestyle='none', marker='s', markersize=7.5, color=INK,
        markeredgecolor=SURFACE, markeredgewidth=1.3, label='Experiment (2010)')[0]
    handles['mses'] = ax.plot(lifts, [references[c]['mses'] for c in lifts], linestyle='none', marker='D', markersize=6.5, color=MUTED,
        markeredgecolor=SURFACE, markeredgewidth=1.3, label='Published MSES (2010)')[0]
    ax.annotate('fixed-incidence samples', (.1549, 112.2), xytext=(0, 9), textcoords='offset points', ha='center', fontsize=8.5, color=SECONDARY)
    ax.set_xlabel('Lift coefficient'); ax.set_ylabel('Drag [counts]'); ax.set_ylim(55, 125)
    ax.set_title('Matched-lift drag: the low-lift error is unchanged; the convex nose keeps it up to CL 0.175', fontsize=10.5, color=INK, loc='left')
    fig.legend(list(handles.values()), [h.get_label() for h in handles.values()], frameon=False, fontsize=8.5,
               loc='outside lower center', ncols=2, labelcolor=SECONDARY)
    for suffix in ('png', 'pdf', 'svg'): fig.savefig(out/f'drag_against_lift.{suffix}', dpi=170, facecolor=SURFACE)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=OUT)
    out = parser.parse_args().out.resolve()
    audit = json.loads((out/'residual_recheck.json').read_text())
    require(audit['complete'] and len(audit['states']) == audit['saved_states'], 'Independent audit is not complete')
    inputs = Inputs(); spec = inputs.json(STUDY/'study.json'); screens = spec['equivalence_screens']; angles = spec['path_angles_deg']
    for name, digest in audit['input_sha256'].items():
        if name.startswith(run.STUDY+'/') and 'source_snapshot' not in name:
            require(sha(ROOT/name) == digest, 'Audited output changed: '+name)
    jobs = {j['name']: j for j in spec['jobs']}
    # ---- fixed-incidence paths of the convex nose
    rows, points = [], {}
    for direction in ('increasing', 'decreasing'):
        job = jobs['convex_nose_'+direction]; folder = STUDY/job['name']
        for rec in inputs.json(folder/'result.json')['records']:
            rows.append(describe(job, rec, folder))
            if rec['converged']: points[(direction, rec['requested_alpha_deg'])] = dict(record=rec, folder=folder)
    keys = list(dict.fromkeys(k for r in rows for k in r)); csv_write(out/'path_points.csv', [{k: r.get(k) for k in keys} for r in rows])
    pairs, convex = [], []
    for alpha in angles:
        a, b = points.get(('increasing', alpha)), points.get(('decreasing', alpha))
        row = dict(alpha_deg=alpha, status='both_converged' if a and b else 'neither_converged' if not a and not b else
                   ('only_increasing_converged' if a else 'only_decreasing_converged'))
        row.update(compare(a, b, screens) if a and b else dict(within_screens=None)); pairs.append(row)
        state = a or b; rec = state['record'] if state else None
        convex.append(dict(scenario='convex_nose', alpha_deg=alpha, CL=rec['CL'] if rec else None, CD_counts=rec['CD']*1e4 if rec else None,
            fore_lower_transition_x=transition_x(rec['boundary_layers']['fore']['lower']) if rec else None,
            converged_from=' '.join(d for d, s in (('increasing', a), ('decreasing', b)) if s) or 'none'))
    csv_write(out/'direction_pairs.csv', [{k: r.get(k) for k in dict.fromkeys(k for r in pairs for k in r)} for r in pairs])
    # ---- same incidences for the other two reconstructions (audited V2 output)
    response = list(convex)
    for r in inputs.csv(COMBINED):
        if r['distinct_converged_states'] == 'True': continue
        value = lambda key: next((float(r[f'{d}_{key}']) for d in ('increasing', 'decreasing', 'cold') if r[f'{d}_{key}'] != ''), None)
        response.append(dict(scenario=r['scenario'], alpha_deg=float(r['alpha_deg']), CL=value('CL'), CD_counts=value('CD_counts'),
                             fore_lower_transition_x=value('fore_lower_transition_x'), converged_from=r['converged_from']))
    csv_write(out/'incidence_response.csv', response)
    other = {(r['scenario'], r['alpha_deg']): r for r in response}
    sensitivity = []
    for r in convex:
        base = other.get(('nominal', r['alpha_deg']))
        if r['CL'] is not None and base and base['CL'] is not None:
            sensitivity.append(dict(alpha_deg=r['alpha_deg'], comparison='convex nose minus nominal, same incidence, different lift',
                delta_CL=r['CL']-base['CL'], delta_CD_counts=r['CD_counts']-base['CD_counts'],
                nominal_fore_lower_transition_x=base['fore_lower_transition_x'], convex_fore_lower_transition_x=r['fore_lower_transition_x']))
    csv_write(out/'same_incidence_sensitivity.csv', sensitivity)
    # ---- one monotone bracket of CL 0.175 on the paths decides admission of that target
    lift = [(r['alpha_deg'], r['CL']) for r in convex if r['CL'] is not None]
    brackets = [(a0, a1) for (a0, c0), (a1, c1) in zip(lift, lift[1:]) if (c0-.175)*(c1-.175) < 0]
    monotone = all(c1 > c0 for (_, c0), (_, c1) in zip(lift, lift[1:]))
    comparable = [p for p in pairs if p['status'] == 'both_converged']
    admit_0175 = bool(len(brackets) == 1 and monotone and all(p['within_screens'] for p in comparable))
    # ---- matched-lift targets
    reference_rows = {r['case']: r for r in inputs.csv(REFERENCES)}
    references = {cl: dict(experiment=float(reference_rows[c]['experiment_CD'])*1e4, mses=float(reference_rows[c]['mses_CD'])*1e4) for cl, c in CASES.items()}
    states = {('nominal', .056): fixed(JACOBIAN/'cl0056_corrected_af'), ('nominal', .175): fixed(LOWER/'cl0175_nominal'),
              ('nominal', .284): fixed(JACOBIAN/'cl028_corrected_af'),
              ('fore_source_cubic', .056): fixed(FORE/'cl0056_fore_source_cubic'), ('fore_source_cubic', .175): fixed(FORE/'cl0175_fore_source_cubic'),
              ('convex_nose', .056): fixed(STUDY/'convex_nose_cl0056'), ('convex_nose', .175): fixed(STUDY/'convex_nose_cl0175'),
              ('convex_nose', .284): fixed(STUDY/'convex_nose_cl0284')}
    for name in ('cl0056_corrected_af', 'cl028_corrected_af'): inputs.add(JACOBIAN/name/'result.json')
    inputs.add(LOWER/'cl0175_nominal/result.json'); inputs.add(FORE/'cl0056_fore_source_cubic/result.json'); inputs.add(FORE/'cl0175_fore_source_cubic/result.json')
    control = fixed(STUDY/'nominal_cl0056_control')
    require(control['CD_counts'] == states[('nominal', .056)]['CD_counts'] and audit['controls']['passed'], 'Nominal control differs from the archive')
    targets = []
    for (geometry, cl), state in states.items():
        admitted = state is not None and (geometry, cl) != ('convex_nose', .175) or (state is not None and admit_0175)
        row = dict(geometry=geometry, target_cl=cl, reached=state is not None, admitted=bool(admitted),
            note='' if state is not None else 'target not reached; no force', experiment_CD_counts=references[cl]['experiment'],
            mses_CD_counts=references[cl]['mses'], mses_error_CD_counts=references[cl]['mses']-references[cl]['experiment'])
        row.update({k: (state[k] if state else None) for k in ('alpha_deg', 'CL', 'CD_counts', 'fore_CD_counts', 'aft_CD_counts',
                                                              'fore_lower_transition_x', 'fore_lower_max_Hk', 'selected_path')})
        row['error_CD_counts'] = state['CD_counts']-references[cl]['experiment'] if state and admitted else None
        targets.append(row)
    csv_write(out/'matched_lift_targets.csv', targets)
    find = lambda geometry, cl: next(t for t in targets if t['geometry'] == geometry and t['target_cl'] == cl)
    paired = []
    for geometry in ('nominal', 'fore_source_cubic', 'convex_nose'):
        low, high = find(geometry, .056), find(geometry, .175)
        paired.append(dict(series=geometry, CD_change_counts=high['CD_counts']-low['CD_counts'] if low['admitted'] and high['admitted'] else None,
                           note='' if low['admitted'] and high['admitted'] else 'CL 0.175 target missing'))
    paired += [dict(series='experiment', CD_change_counts=references[.175]['experiment']-references[.056]['experiment'], note=''),
               dict(series='published MSES 2010', CD_change_counts=references[.175]['mses']-references[.056]['mses'], note='')]
    csv_write(out/'paired_drag_change.csv', paired)
    # ---- supported nose taps at the two matched lifts
    taps = []
    tap_rows = inputs.csv(TAPS)
    for cl in (.056, .175):
        for r in tap_rows:
            if (r['case'] == CASES[cl] and r['method'] == 'multifoil' and r['coordinate'] == 'table' and r['included'] == 'True'
                    and r['surface'] == 'fore_lower' and float(r['x_over_c']) <= .06):
                x = float(r['x_over_c']); nominal, new = states[('nominal', cl)], states[('convex_nose', cl)]
                row = dict(target_cl=cl, x_over_c=x, experiment_Cp=float(r['experiment_Cp']), nominal_Cp=tap_cp(nominal, x), convex_nose_Cp=tap_cp(new, x))
                require(abs(row['nominal_Cp']-float(r['model_Cp'])) < 1e-9, 'Nominal tap pressure differs from the archived comparison')
                row.update(nominal_error=row['nominal_Cp']-row['experiment_Cp'], convex_nose_error=row['convex_nose_Cp']-row['experiment_Cp'])
                taps.append(row)
    csv_write(out/'nose_taps.csv', taps)
    rmse = lambda cl, key: float(np.sqrt(np.mean([t[key]**2 for t in taps if t['target_cl'] == cl])))
    figures(out, response, targets, references)
    path = [r for r in rows if r['role'] == 'path']
    summary = dict(complete=True, audited_saved_states=audit['saved_states'], maximum_fresh_residual=max(s['max_residual'] for s in audit['states']),
        nominal_control_identical_to_archive=audit['controls']['passed'],
        requested_path_points=len(path), converged_path_points=sum(r['converged'] for r in path),
        failed_path_points=[dict(job=r['job'], alpha_deg=r['alpha_deg'], reached_by=r['reached_by'], failure_reason=r['failure_reason'])
                            for r in path if not r['converged']],
        direction_pairs_both_converged=len(comparable), direction_pairs_within_screens=sum(bool(p['within_screens']) for p in comparable),
        largest_direction_difference=dict(CL=max(abs(p['delta_CL']) for p in comparable), CD_counts=max(abs(p['delta_CD_counts']) for p in comparable),
            Cp=max(p['max_abs_delta_Cp'] for p in comparable), transition_x=max(p['max_transition_shift_x'] for p in comparable)),
        convex_lift_monotone_on_paths=bool(monotone), convex_CL0175_brackets=brackets, convex_CL0175_admitted=admit_0175,
        convex_transition_range_on_paths=[min(r['fore_lower_transition_x'] for r in convex if r['CL'] is not None),
                                          max(r['fore_lower_transition_x'] for r in convex if r['CL'] is not None)],
        targets=[{k: v for k, v in t.items()} for t in targets], paired_drag_change=paired,
        nose_tap_rmse={str(cl): dict(nominal=rmse(cl, 'nominal_error'), convex_nose=rmse(cl, 'convex_nose_error'),
                                     taps=sum(t['target_cl'] == cl for t in taps)) for cl in (.056, .175)},
        same_incidence_sensitivity=sensitivity, equivalence_screens=screens, new_flow_solves=0,
        qualification='Reconstruction sensitivity on a figure-derived polygon; neither contour is the tested one. '
                      'Sampled incidences only; not solution uniqueness')
    write(out/'summary.json', summary); write(out/'report_input_manifest.json', dict(files=inputs.files))
    print(json.dumps({k: v for k, v in summary.items() if k not in ('targets', 'same_incidence_sensitivity')}, indent=1))


if __name__ == '__main__':
    main()
