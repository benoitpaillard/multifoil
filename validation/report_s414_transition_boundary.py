#!/usr/bin/env python3
"""Tables, declared equivalence screens, tranche-2 selection and figures for V2 tranche 1.

Reads the retrieved study only after its independent audit is complete. No flow solve.
Same-angle comparisons are numerical path diagnostics, never matched-lift drag scores.
"""
import argparse
import json
import math
from pathlib import Path
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from report_design_initialization import ROOT, Inputs, sha, write, csv_write
from report_s702_drag_mechanism import require
import run_s414_transition_boundary as deck

STUDY = ROOT/deck.STUDY
OUT = ROOT/'validation/results/s414_transition_boundary_analysis_20261004'
SURFACES = [(e, s) for e in deck.ELEMENTS for s in ('upper', 'lower')]
PROBE_X, TARGET_CL, TRANSITION_CHANGE = .034, .175, .05
INCREASING, DECREASING = '#2a78d6', '#eb6834'  # validated categorical slots 1 and 2
SURFACE, INK, SECONDARY, MUTED, GRID, AXIS = '#fcfcfb', '#0b0b0b', '#52514e', '#898781', '#e1e0d9', '#c3c2b7'
TITLES = {'nominal': 'Nominal reconstruction', 'fore_source_cubic': 'Fore-only cubic reconstruction'}


def surface(folder, index, element):
    return np.genfromtxt(folder/f'flow/surface_{index:03d}_{element}.csv', delimiter=',', names=True, dtype=None, encoding='utf8')


def transition_x(layer):
    return None if layer['transition'] is None else layer['transition']['xy_over_cref'][0]


def describe(job, rec, folder):
    """One requested evaluation as a flat row; a failed one keeps only its numerical record."""
    row = dict(job=job['name'], scenario=job['scenario'], mode=job['amplification_jacobian'], direction=job['direction'],
        index=rec['index'], role=rec['role'], alpha_deg=rec['requested_alpha_deg'], converged=rec['converged'],
        approach=rec['approach'], reached_by=rec['reached_by'], iterations_final_step=rec['iterations_final_step'],
        residual=rec['residual'], CL=rec['CL'], CD=rec['CD'], CM=rec['CM'], failure_reason=rec['failure_reason'])
    if not rec['converged']: return row
    row.update(CD_counts=rec['CD']*1e4, fore_CD_counts=rec['element_forces']['fore']['CD_wake']*1e4,
               aft_CD_counts=rec['element_forces']['aft']['CD_wake']*1e4)
    for element, side in SURFACES:
        layer = rec['boundary_layers'][element][side]; key = f'{element}_{side}'
        row[key+'_transition_x'] = transition_x(layer); row[key+'_min_cf'] = layer['min_cf']; row[key+'_max_Hk'] = layer['max_Hk']
        row[key+'_separation_regions'] = len(layer['separation_regions'])
        row[key+'_separated_fraction'] = layer['separated_arclength_fraction']
    values = surface(folder, rec['index'], 'fore'); lower = values[values['side'] == 'lower']
    laminar = lower[lower['turbulent'] == 0]
    # Amplification at a fixed fore-lower position, only where the layer is still laminar there.
    row['fore_lower_N_at_x0034'] = (float(np.interp(PROBE_X, laminar['x_over_cref'], laminar['amplification_or_shear']))
                                    if len(laminar) > 1 and laminar['x_over_cref'][0] <= PROBE_X <= laminar['x_over_cref'][-1] else None)
    return row


def compare(a, b, screens):
    """Declared same-angle equivalence screens between two converged evaluations."""
    ra, fa, rb, fb = a['record'], a['folder'], b['record'], b['folder']
    cp, topology = 0., True
    for element in deck.ELEMENTS:
        sa, sb = surface(fa, ra['index'], element), surface(fb, rb['index'], element)
        cp = max(cp, float(np.max(abs(sa['Cp']-sb['Cp'])))); topology = topology and bool(np.array_equal(sa['turbulent'], sb['turbulent']))
    shift = 0.
    for element, side in SURFACES:
        la, lb = ra['boundary_layers'][element][side], rb['boundary_layers'][element][side]
        xa, xb = transition_x(la), transition_x(lb)
        topology = topology and (xa is None) == (xb is None) and len(la['separation_regions']) == len(lb['separation_regions'])
        if xa is not None and xb is not None: shift = max(shift, abs(xa-xb))
    out = dict(delta_CL=ra['CL']-rb['CL'], delta_CD_counts=(ra['CD']-rb['CD'])*1e4, max_abs_delta_Cp=cp,
               max_transition_shift_x=shift, same_topology=topology)
    out['within_screens'] = bool(abs(out['delta_CL']) <= screens['CL'] and abs(out['delta_CD_counts']) <= screens['CD_counts']
        and cp <= screens['Cp'] and shift <= screens['transition_x_over_c'] and topology)
    return out


def select_refinement(angles, paths, pairs):
    """The tranche-2 rule declared in the tranche-1 protocol; drag is never an input."""
    intervals = []
    for k in range(len(angles)-1):
        ends = (angles[k], angles[k+1]); t = 0.
        discordant = any(pairs[a]['status'] != 'both_converged' or not pairs[a]['within_screens'] for a in ends)
        for direction in ('increasing', 'decreasing'):
            lo, hi = (paths[direction].get(a) for a in ends)
            if lo is None and hi is None: continue
            if lo is None or hi is None: t = math.inf; continue
            xa, xb = (transition_x(p['record']['boundary_layers']['fore']['lower']) for p in (lo, hi))
            t = max(t, math.inf if xa is None or xb is None else abs(xa-xb))
        intervals.append(dict(lower_alpha_deg=ends[0], upper_alpha_deg=ends[1], path_difference=bool(discordant),
            fore_lower_transition_change=None if math.isinf(t) else t, transition_change_unbounded=math.isinf(t),
            eligible=bool(discordant or t > TRANSITION_CHANGE), _t=t))
    ranked = sorted((i for i in intervals if i['eligible']), key=lambda i: (not i['path_difference'], -i['_t'], i['lower_alpha_deg']))
    locations = [dict(alpha_deg=.5*(i['lower_alpha_deg']+i['upper_alpha_deg']), rule='midpoint',
                      interval=[i['lower_alpha_deg'], i['upper_alpha_deg']]) for i in ranked][:4]
    if ranked and len(locations) < 4:
        top = ranked[0]; width = top['upper_alpha_deg']-top['lower_alpha_deg']
        for fraction in (.25, .75):
            if len(locations) < 4:
                locations.append(dict(alpha_deg=top['lower_alpha_deg']+fraction*width, rule='quarter point of top-ranked interval',
                                      interval=[top['lower_alpha_deg'], top['upper_alpha_deg']]))
    discordant_angles = [a for a in angles if pairs[a]['status'] != 'both_converged' or not pairs[a]['within_screens']]
    cold = (dict(alpha_deg=discordant_angles[0], rule='angle with a path difference') if discordant_angles else
            dict(alpha_deg=locations[0]['alpha_deg'], rule='midpoint of top-ranked interval') if locations else None)
    for i in intervals: del i['_t']
    return dict(intervals=intervals, ranked_eligible=[[i['lower_alpha_deg'], i['upper_alpha_deg']] for i in ranked],
        new_locations=locations, cold_check=cold,
        fallback='No eligible interval: apply the same rule to the bracket from the archived CL 0.056 to 0.175 history'
                 if not ranked else None)


def style(ax):
    ax.set_facecolor(SURFACE); ax.grid(color=GRID, linewidth=.6); ax.set_axisbelow(True)
    for side in ('top', 'right'): ax.spines[side].set_visible(False)
    for side in ('left', 'bottom'): ax.spines[side].set_color(AXIS)
    ax.tick_params(colors=MUTED, labelsize=8.5); ax.xaxis.label.set_color(SECONDARY); ax.yaxis.label.set_color(SECONDARY)


def figure(out, name, title, rows, measures, angles):
    scenarios = list(TITLES)
    # One vertical scale per measure across both geometries, so a flat response reads as flat.
    fig, axes = plt.subplots(len(measures), 2, figsize=(10, 2.35*len(measures)+.9), sharex=True, sharey='row',
                             layout='constrained', squeeze=False)
    fig.patch.set_facecolor(SURFACE); handles = {}
    for column, scenario in enumerate(scenarios):
        for line, (key, label, log) in enumerate(measures):
            ax = axes[line, column]; style(ax)
            for direction, color, marker, size in (('increasing', INCREASING, '>', 11), ('decreasing', DECREASING, '<', 7)):
                data = sorted((r for r in rows if r['scenario'] == scenario and r['direction'] == direction and r['role'] == 'path'),
                              key=lambda r: r['alpha_deg'])
                # A failed or unavailable point is a gap: the line never bridges it.
                values = [r[key] if r['converged'] and r.get(key) is not None else np.nan for r in data]
                handles[direction] = ax.plot([r['alpha_deg'] for r in data], values, color=color, linewidth=1.6,
                    marker=marker, markersize=size, markeredgecolor=SURFACE, markeredgewidth=1.3, solid_capstyle='round',
                    label=f'{direction.capitalize()} incidence path')[0]
                for r in data:
                    if not r['converged']:
                        # A failed point has no value: mark its incidence on the axis, never a number.
                        handles['failed'] = ax.plot([r['alpha_deg']], [0.03], transform=ax.get_xaxis_transform(), linestyle='none',
                            marker='x', markersize=8, markeredgewidth=1.8, color=INK, clip_on=False, label='Not converged (no value)')[0]
                        if line == 0:
                            ax.annotate(f'{direction} path:\nnot converged', (r['alpha_deg'], .03), xycoords=ax.get_xaxis_transform(),
                                xytext=(8, 10), textcoords='offset points', fontsize=8, color=SECONDARY)
            if log: ax.set_yscale('log')
            if column == 0: ax.set_ylabel(label)
            if line == 0: ax.set_title(TITLES[scenario], fontsize=10, color=INK, loc='left')
            if line == len(measures)-1:
                ax.set_xlabel('Incidence [deg]'); ax.set_xticks(angles, [f'{a:.3f}' for a in angles])
    fig.legend(list(handles.values()), [h.get_label() for h in handles.values()], frameon=False, fontsize=8.5,
               loc='outside lower center', ncols=3, labelcolor=SECONDARY)
    fig.suptitle(title, fontsize=11.5, color=INK, x=.01, ha='left')
    for suffix in ('png', 'pdf', 'svg'): fig.savefig(out/f'{name}.{suffix}', dpi=170, facecolor=SURFACE)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=OUT)
    out = parser.parse_args().out.resolve()
    audit = json.loads((out/'residual_recheck.json').read_text())
    require(audit['complete'] and len(audit['states']) == audit['saved_states'], 'Independent audit is not complete')
    inputs = Inputs(); spec = inputs.json(STUDY/'study.json'); screens = spec['equivalence_screens']; angles = spec['angles_deg']
    for name, digest in audit['input_sha256'].items():
        if name.startswith(deck.STUDY) and 'source_snapshot' not in name:
            require(sha(ROOT/name) == digest, 'Audited output changed: '+name)
    rows, points = [], {}
    for job in spec['jobs']:
        folder = STUDY/job['name']; result = inputs.json(folder/'result.json')
        for rec in result['records']:
            rows.append(describe(job, rec, folder)); points[(job['name'], rec['index'])] = dict(record=rec, folder=folder, job=job)
        for request in job['requests'][len(result['records']):]:
            rows.append(dict(job=job['name'], scenario=job['scenario'], mode=job['amplification_jacobian'], direction=job['direction'],
                index=None, role=request['role'], alpha_deg=request['alpha_deg'], converged=False, reached_by='not evaluated',
                failure_reason=result['error']))
    keys = list(dict.fromkeys(k for r in rows for k in r)); csv_write(out/'points.csv', [{k: r.get(k) for k in keys} for r in rows])
    path = lambda scenario, direction: {p['record']['requested_alpha_deg']: p for (name, _), p in points.items()
        if name == f'{scenario}_corrected_{direction}' and p['record']['role'] == 'path' and p['record']['converged']}
    directions, modes, trips, geometry, monotone, selection = [], [], [], [], [], {}
    for scenario in deck.ARCHIVED_JOBS:
        paths = dict(increasing=path(scenario, 'increasing'), decreasing=path(scenario, 'decreasing')); pairs = {}
        for alpha in angles:
            a, b = paths['increasing'].get(alpha), paths['decreasing'].get(alpha)
            row = dict(scenario=scenario, alpha_deg=alpha, status='both_converged' if a and b else
                       'neither_converged' if not a and not b else ('only_increasing_converged' if a else 'only_decreasing_converged'))
            row.update(compare(a, b, screens) if a and b else dict(within_screens=None))
            pairs[alpha] = row; directions.append(row)
        control = {r['requested_alpha_deg']: dict(record=r, folder=STUDY/f'{scenario}_legacy_endpoints')
                   for r in inputs.json(STUDY/f'{scenario}_legacy_endpoints/result.json')['records']}
        lo, hi = angles[0], angles[-1]
        for alpha, corrected, basis in ((lo, paths['increasing'].get(lo), 'recorded initialization of the lower endpoint'),
                                        (hi, paths['decreasing'].get(hi), 'recorded initialization of the upper endpoint (through the lower one)')):
            row = dict(scenario=scenario, alpha_deg=alpha, comparison='corrected_af minus legacy', initialization=basis)
            if corrected and control[alpha]['record']['converged']:
                row.update(compare(corrected, control[alpha], screens),
                    corrected_iterations=corrected['record']['iterations_final_step'], legacy_iterations=control[alpha]['record']['iterations_final_step'],
                    corrected_reached_by=corrected['record']['reached_by'], legacy_reached_by=control[alpha]['record']['reached_by'])
            modes.append(row)
        seed = points.get((f'{scenario}_corrected_decreasing', 0)); last = paths['decreasing'].get(lo)
        if seed and seed['record']['converged'] and last:
            trips.append(dict(scenario=scenario, alpha_deg=lo, comparison='return to the lower endpoint minus its seed', **compare(last, seed, screens)))
        merged = {a: paths['increasing'].get(a) or paths['decreasing'].get(a) for a in angles}
        ordered = [(a, merged[a]['record']) for a in angles if merged[a]]
        for (a0, r0), (a1, r1) in zip(ordered, ordered[1:]):
            monotone.append(dict(scenario=scenario, lower_alpha_deg=a0, upper_alpha_deg=a1, lower_CL=r0['CL'], upper_CL=r1['CL'],
                delta_CL=r1['CL']-r0['CL'], lift_increases=bool(r1['CL'] > r0['CL']),
                brackets_target_CL=bool((r0['CL']-TARGET_CL)*(r1['CL']-TARGET_CL) < 0)))
        selection[scenario] = select_refinement(angles, paths, pairs)
    nominal, cubic = (dict(increasing=path(s, 'increasing'), decreasing=path(s, 'decreasing')) for s in deck.ARCHIVED_JOBS)
    for alpha in angles:
        a = cubic['increasing'].get(alpha) or cubic['decreasing'].get(alpha); b = nominal['increasing'].get(alpha) or nominal['decreasing'].get(alpha)
        if a and b:
            ra, rb = a['record'], b['record']
            geometry.append(dict(alpha_deg=alpha, comparison='cubic minus nominal, same incidence, different lift',
                delta_CL=ra['CL']-rb['CL'], delta_CD_counts=(ra['CD']-rb['CD'])*1e4,
                fore_delta_CD_counts=(ra['element_forces']['fore']['CD_wake']-rb['element_forces']['fore']['CD_wake'])*1e4,
                aft_delta_CD_counts=(ra['element_forces']['aft']['CD_wake']-rb['element_forces']['aft']['CD_wake'])*1e4,
                nominal_fore_lower_transition_x=transition_x(rb['boundary_layers']['fore']['lower']),
                cubic_fore_lower_transition_x=transition_x(ra['boundary_layers']['fore']['lower'])))
    uniform = lambda data: [{k: r.get(k) for k in dict.fromkeys(k for r in data for k in r)} for r in data]
    for name, data in (('direction_pairs.csv', directions), ('mode_pairs.csv', modes), ('round_trips.csv', trips),
                       ('geometry_pairs.csv', geometry), ('lift_response.csv', monotone)):
        csv_write(out/name, uniform(data))
    write(out/'tranche2_selection.json', selection)
    figure(out, 'path_response', 'Fixed-incidence paths, corrected Jacobian: lift, drag and fore-lower transition', rows,
        [('CL', 'Lift coefficient', False), ('CD_counts', 'Drag [counts]', False),
         ('fore_lower_transition_x', 'Fore-lower transition, x / chord', False)], angles)
    figure(out, 'path_diagnostics', 'Fixed-incidence paths: element drag, fore-lower boundary layer and residual', rows,
        [('fore_CD_counts', 'Fore wake drag [counts]', False), ('aft_CD_counts', 'Aft wake drag [counts]', False),
         ('fore_lower_min_cf', 'Fore-lower minimum Cf', False), ('fore_lower_max_Hk', 'Fore-lower maximum Hk', False),
         ('fore_lower_N_at_x0034', 'N at x/c 0.034 (laminar)', False), ('residual', 'Newton residual', True)], angles)
    both = [r for r in directions if r['status'] == 'both_converged']
    path_rows = [r for r in rows if r['role'] == 'path']
    summary = dict(complete=True, audited_saved_states=audit['saved_states'], maximum_fresh_residual=max(s['max_residual'] for s in audit['states']),
        requested_path_points=len(path_rows), converged_path_points=sum(r['converged'] for r in path_rows),
        failed_path_points=[dict(job=r['job'], alpha_deg=r['alpha_deg'], reached_by=r['reached_by'], failure_reason=r['failure_reason'])
                            for r in path_rows if not r['converged']],
        direction_pairs=len(directions), direction_pairs_both_converged=len(both),
        direction_pairs_within_screens=sum(bool(r['within_screens']) for r in both),
        largest_direction_difference=dict(CL=max(abs(r['delta_CL']) for r in both), CD_counts=max(abs(r['delta_CD_counts']) for r in both),
            Cp=max(r['max_abs_delta_Cp'] for r in both), transition_x=max(r['max_transition_shift_x'] for r in both)),
        mode_pairs_within_screens=sum(bool(r.get('within_screens')) for r in modes), mode_pairs=len(modes),
        round_trips_within_screens=sum(bool(r['within_screens']) for r in trips), round_trips=len(trips),
        nonmonotone_lift_intervals=[r for r in monotone if not r['lift_increases']],
        intervals_bracketing_target_CL={s: [[r['lower_alpha_deg'], r['upper_alpha_deg']] for r in monotone
                                           if r['scenario'] == s and r['brackets_target_CL']] for s in deck.ARCHIVED_JOBS},
        achieved_angular_resolution_deg=angles[1]-angles[0], equivalence_screens=screens, tranche2_selection=selection,
        new_flow_solves=0, qualification='Sampled-path numerical diagnostic at fixed incidence; not matched-lift validation, '
                                         'not physical hysteresis, not solution uniqueness')
    write(out/'summary.json', summary); write(out/'report_input_manifest.json', dict(files=inputs.files))
    print(json.dumps({k: v for k, v in summary.items() if k != 'tranche2_selection'}, indent=1))


if __name__ == '__main__':
    main()
