#!/usr/bin/env python3
"""Audit archived S702 states and explain modeled element drag changes; no solves.

Run from the project: python3 validation/report_s702_drag_mechanism.py --out PATH
See S702_DRAG_MECHANISM_PROTOCOL.md for the predeclared scope and tolerances.
"""
import argparse
import csv
import importlib.metadata
import json
import platform
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from plot_s702_flap_effect import from_analysis, CASES, TARGETS
from report_s702_ncrit_sweep import (
    ROOT, RESULTS, DATA, read, sha, csv_write, flow_directory, load_study,
)

ELEMENTS = ('fore', 'aft')
SIDES = ('upper', 'lower')
NCRITS = (9., 6.25)
FORCE_ATOL = 1e-12
DELTA_ATOL = 1e-8
COLORS = {5: '#247AA5', 10: '#D16726'}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def close(actual, expected, tolerance=FORCE_ATOL, context='value'):
    require(np.isfinite(actual) and np.isfinite(expected)
            and abs(actual-expected) <= tolerance, f'{context}: {actual} != {expected}')


def relative(path):
    return str(Path(path).resolve().relative_to(ROOT))


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')


class Inputs:
    """Pin consumed files to the delivered archive and record report provenance."""
    def __init__(self):
        self.original = {r['path']: r for r in read(ROOT/'handover/archive_manifest.json')['files']}
        self.files = {}

    def add(self, path):
        path = Path(path).resolve()
        key = relative(path)
        if key not in self.files:
            digest = sha(path)
            original = self.original.get(key)
            if original:
                require(digest == original['sha256'], f'Archived input changed: {key}')
            self.files[key] = dict(sha256=digest, bytes=path.stat().st_size,
                                  verified_against_delivered_archive=original is not None)
        return path

    def json(self, path):
        return read(self.add(path))

    def csv(self, path):
        with self.add(path).open(newline='') as stream:
            return list(csv.DictReader(stream))

    def array(self, path):
        data = np.atleast_1d(np.genfromtxt(self.add(path), delimiter=',', names=True,
                                         dtype=None, encoding='utf-8'))
        for field in data.dtype.names:
            if data[field].dtype.kind in 'fiu':
                require(np.isfinite(data[field]).all(), f'Nonfinite {field}: {path}')
        return data

    def sources(self, spec, folder):
        for name, expected in spec['sha256'].items():
            path = ROOT/name
            if sha(path) != expected:
                path = folder/'source_snapshot'/name
            require(path.is_file() and sha(path) == expected, f'Study source changed: {name}')
            self.add(path)


def selected_index(result, polar):
    """Use the recorded index; require unambiguous force matching for old exports."""
    matches = [i for i, p in enumerate(polar) if all(
        abs(float(p[key])-result[key]) <= FORCE_ATOL for key in ('CL', 'CD', 'CM', 'alpha_deg'))]
    if 'selected_flow_index' in result:
        index = result['selected_flow_index']
        require(index in matches, 'Selected flow index disagrees with returned forces')
    else:
        require(len(matches) == 1, 'Legacy selected flow is missing or ambiguous')
        index = matches[0]
    require(int(polar[index]['converged']) == 1, 'Selected polar state did not converge')
    return index


def side_data(surface, side):
    values = surface[surface['side'] == side]
    # Export order is upper TE -> stagnation -> lower TE. Do not sort by x:
    # the flow may cross a surface whose projected x is locally nonmonotone.
    return values[::-1] if side == 'upper' else values


def separation_endpoints(values):
    """Reconstruct Cf-negative interval endpoints from exported nodes only."""
    xy = np.column_stack((values['x_over_cref'], values['y_over_cref']))
    s = np.r_[0., np.cumsum(np.linalg.norm(np.diff(xy, axis=0), axis=1))]
    cf = values['cf']
    breaks = [s[0]]
    for i in range(len(s)-1):
        if cf[i]*cf[i+1] < 0:
            breaks.append(s[i]-cf[i]*(s[i+1]-s[i])/(cf[i+1]-cf[i]))
    breaks.append(s[-1])
    regions = []
    for a, b in zip(breaks[:-1], breaks[1:]):
        if np.interp((a+b)/2, s, cf) < -1e-9:
            regions.append([np.array([np.interp(t, s, xy[:,k]) for k in (0, 1)]) for t in (a, b)])
    return regions


def audit_surface(values, layer):
    close(float(min(values['cf'])), layer['min_cf'], context='surface min Cf')
    close(float(max(values['Hk'])), layer['max_Hk'], context='surface max Hk')
    regions = separation_endpoints(values)
    require(len(regions) == len(layer['separation_regions']), 'Separation interval count mismatch')
    for xy, region in zip(regions, layer['separation_regions']):
        for point, end in zip(xy, ('start', 'end')):
            require(np.max(abs(point-region[end]['xy_over_cref'])) < 1e-10,
                    'Separation endpoint differs from selected surface export')
    turbulent = values['turbulent'].astype(bool)
    tr = layer['transition']
    require(bool(turbulent.any()) == (tr is not None), 'Transition flag mismatch')
    if tr is not None:
        first = np.flatnonzero(turbulent)[0]
        require(turbulent[first:].all(), 'Unexpected relaminarization flag')
        if first:
            points = np.column_stack((values['x_over_cref'], values['y_over_cref']))
            a, b = points[first-1:first+1]
            t = float(np.dot(np.array(tr['xy_over_cref'])-a, b-a)/np.dot(b-a, b-a))
            require(-1e-9 <= t <= 1+1e-9, 'Transition outside laminar/turbulent cell')
            require(np.max(abs(a+t*(b-a)-tr['xy_over_cref'])) < 1e-10,
                    'Transition differs from selected surface export')


def load_state(point, spec, inputs):
    path = ROOT/point['result_path']
    result = inputs.json(path)
    require(result['feasible'], 'Failed result cannot supply a force state')
    for key, expected in [('CL', point['achieved_CL']), ('CD', point['CD']),
                          ('ncrit', point['ncrit']), ('flap', point['flap_deg']),
                          ('target_cl', point['target_cl'])]:
        close(result[key], float(expected), context=key)
    flow = flow_directory(result)
    polar = inputs.csv(flow/'polar.csv')
    index = selected_index(result, polar)
    p = polar[index]
    metadata = inputs.json(flow/'metadata.json')
    diag = metadata['diagnostics'][index]
    require(diag['status'] == 'converged' and not diag['failure_reason'], 'Nonconverged diagnostics')
    if result.get('reused'):
        source = inputs.json(ROOT/result['source_path'])
        require(sha(ROOT/result['source_path']) == result['source_sha256'], 'Reused source mismatch')
        source_target = next(t for t in source['targets'] if t['target_cl'] == result['target_cl'])
        for key in ('CL', 'CD', 'CM', 'alpha_deg'):
            close(result[key], source_target[key], context='reused '+key)
        records = [r for r in result['records'] if all(
            abs(r[k]-result[k]) <= FORCE_ATOL for k in ('CL', 'CD', 'alpha_deg'))]
        require(len(records) == 1, 'Ambiguous reused record')
        record = records[0]
    else:
        record = result['records'][index]
    require(record['converged'] and record['admissible'], 'Selected record inadmissible')
    for key in ('CL', 'CD', 'CM', 'alpha_deg', 'residual'):
        close(record[key], float(p[key]), context='record/polar '+key)
    residual = float(p['residual'])
    require(residual <= spec['residual_tolerance'], 'Residual tolerance not met')
    require(abs(result['CL']-result['target_cl']) <= spec['cl_tolerance'], 'Lift tolerance not met')
    history = diag['history'][-1]
    close(residual, max(history['bl_residual_l2'], history['coupling_residual_max']), context='history residual')
    cfg = metadata['case']['study']
    inputs.sources(cfg, flow.parents[1] if result.get('reused') else path.parents[2])
    for key in ('reynolds', 'mach'):
        close(diag[key], point[key], context=key)
        close(cfg[key], spec[key], context='export config '+key)
    close(diag['ncrit'], float(point['ncrit']), context='Ncrit')
    require(all(n == float(point['ncrit']) for n in diag.get('ncrit_by_element', {}).values()),
            'Per-element Ncrit mismatch')
    close(cfg.get('reference_chord', 1.), 1., context='assembly chord')
    close(cfg['wake_length'], 1., context='wake length')
    require(spec['panels'] == [260, 180], 'Panel protocol mismatch')
    panels = cfg.get('panels', cfg.get('settings', {}).get('base'))
    require(panels == spec['panels'], 'Export panel configuration mismatch')
    job = metadata['case']['job']
    close(job.get('flap', job.get('deflection')), result['flap'], context='Export flap')
    for key in ('boundary_layers', 'element_forces'):
        if key in result:
            require(result[key] == diag[key], f'{key} not from selected diagnostics')
    forces = diag['element_forces']
    close(sum(forces[e]['CD_wake'] for e in ELEMENTS), result['CD'], context='wake CD sum')
    for key in ('CL', 'CM'):
        close(sum(forces[e][key] for e in ELEMENTS), result[key], context='element '+key+' sum')
    sid = f"{point['dataset']}_n{float(point['ncrit']):g}_d{result['flap']:+d}_cl{result['target_cl']:g}"
    summary = dict(state_id=sid, dataset=point['dataset'], reynolds=point['reynolds'], mach=point['mach'],
        ncrit=float(point['ncrit']), flap=result['flap'], target_cl=result['target_cl'],
        primary=result['flap'] in (5, 10) and result['target_cl'] in (.8, 1.2),
        CL=result['CL'], CD=result['CD'], CM=result['CM'], alpha_deg=result['alpha_deg'],
        lift_error=result['CL']-result['target_cl'], residual=residual,
        cl_tolerance=spec['cl_tolerance'], residual_tolerance=spec['residual_tolerance'],
        cd_sum_error=result['CD']-sum(forces[e]['CD_wake'] for e in ELEMENTS),
        selected_solve_wall_seconds=diag['wall_seconds'], target_wall_seconds=result.get('wall_seconds'),
        selected_iterations=diag['iterations'], reused=result['reused'], numerical_retry=point['numerical_retry'],
        retry_of=point['retry_of'], result_path=point['result_path'], flow_directory=relative(flow),
        selected_flow_index=index, protocol=spec['protocol'])
    surfaces, wakes, surface_rows, wake_rows, intervals = {}, {}, [], [], []
    for element, panel_count in zip(ELEMENTS, spec['panels']):
        inputs.add(flow/f'geometry_{element}.dat')
        surf = inputs.array(flow/f'surface_{index:03d}_{element}.csv')
        wake = inputs.array(flow/f'wake_{index:03d}_{element}.csv')
        require(len(surf) == panel_count+1, 'Actual surface panel count mismatch')
        require(set(surf['side']) == set(SIDES), 'Surface side labels mismatch')
        require(np.all(wake['theta'] > 0) and np.all(wake['ue_over_vinf'] > 0), 'Invalid wake')
        surfaces[element], wakes[element] = surf, wake
        summary[element+'_CD_wake'] = forces[element]['CD_wake']
        summary[element+'_CD_friction'] = forces[element]['CD_friction']
        wr = dict(state_id=sid, element=element, stations=len(wake))
        for endpoint, k in [('first', 0), ('last', -1)]:
            wr.update({endpoint+'_'+field: float(wake[field][k]) for field in wake.dtype.names})
        wake_rows.append(wr)
        for side in SIDES:
            values = side_data(surf, side)
            layer = diag['boundary_layers'][element][side]
            audit_surface(values, layer)
            tr = layer['transition']
            sr = dict(state_id=sid, dataset=point['dataset'], ncrit=float(point['ncrit']),
                flap=result['flap'], target_cl=result['target_cl'], element=element, side=side,
                transition_x_local=tr['x_over_local_chord'] if tr else None,
                transition_x_assembly=tr['xy_over_cref'][0] if tr else None,
                laminar_to_te=tr is None, separated_arclength_fraction=layer['separated_arclength_fraction'],
                separation_regions=len(layer['separation_regions']), min_cf=layer['min_cf'], max_Hk=layer['max_Hk'])
            sr.update({'te_'+field: float(values[field][-1]) for field in
                       ('x_over_local_chord', 'Cp', 'cf', 'ue_over_vinf', 'theta', 'delta_star', 'Hk')})
            surface_rows.append(sr)
            for j, region in enumerate(layer['separation_regions']):
                ir = dict(state_id=sid, element=element, side=side, interval=j,
                          reaches_trailing_edge=region['reaches_trailing_edge'])
                for end in ('start', 'end'):
                    loc = region[end]
                    ir.update({end+'_x_local': loc['x_over_local_chord'], end+'_x_assembly': loc['xy_over_cref'][0],
                               end+'_y_assembly': loc['xy_over_cref'][1], end+'_s_assembly': loc['s_over_cref']})
                intervals.append(ir)
    provenance = dict(state_id=sid, result_path=relative(path), flow_directory=relative(flow),
        selected_flow_index=index, study_path=relative(path.parents[2]/'study.json'),
        solver_versions=metadata['versions'], numerical_settings={k:v for k,v in cfg.items()
            if k not in ('jobs', 'sha256', 'seed_results')}, executed_source_sha256=cfg['sha256'],
        requested_records=result['records'], attempts=result.get('attempts', []),
        selected_initialization=diag['initialization'], selected_continuation=diag.get('continuation', []),
        selected_newton_history=diag['history'], retry_of=point['retry_of'])
    return dict(summary=summary, surfaces=surfaces, wakes=wakes, layers=diag['boundary_layers'],
                surface_rows=surface_rows, wake_rows=wake_rows, intervals=intervals, provenance=provenance)


def pair_deltas(states, references):
    lookup = {(s['summary']['dataset'], s['summary']['ncrit'], s['summary']['target_cl'], s['summary']['flap']):s
              for s in states}
    refs = {(r['dataset'], r['target_cl'], r['flap_deg'], r['series']):r for r in references}
    reproduced, contributions, surface_deltas = [], [], []
    for dataset, re, mach, _, _ in CASES:
        for cl in TARGETS:
            row = dict(dataset=dataset, reynolds=re, mach=mach, target_cl=cl)
            for series in ('experiment', 'mses_n9', 'multifoil_n9', 'multifoil_n6p25'):
                if series.startswith('multifoil'):
                    n = 9. if series == 'multifoil_n9' else 6.25
                    a, b = (lookup[(dataset, n, cl, f)]['summary'] for f in (5, 10))
                    value = (b['CD']-a['CD'])*1e4
                else:
                    value = (refs[(dataset, cl, 10, series)]['CD']-refs[(dataset, cl, 5, series)]['CD'])*1e4
                row[series+'_delta_10_minus_5_counts'] = value
            reproduced.append(row)
            for n in NCRITS:
                a, b = (lookup[(dataset, n, cl, f)] for f in (5, 10))
                sa, sb = a['summary'], b['summary']
                delta = dict(dataset=dataset, reynolds=re, mach=mach, target_cl=cl, ncrit=n,
                    state_5=sa['state_id'], state_10=sb['state_id'],
                    delta_CD_counts=(sb['CD']-sa['CD'])*1e4,
                    experiment_delta_counts=row['experiment_delta_10_minus_5_counts'],
                    mses_n9_delta_counts=row['mses_n9_delta_10_minus_5_counts'])
                delta['delta_error_counts'] = delta['delta_CD_counts']-delta['experiment_delta_counts']
                for element in ELEMENTS:
                    for kind in ('wake', 'friction'):
                        field = element+'_CD_'+kind
                        delta[element+'_delta_'+kind+'_counts'] = (sb[field]-sa[field])*1e4
                close(delta['fore_delta_wake_counts']+delta['aft_delta_wake_counts'],
                      delta['delta_CD_counts'], DELTA_ATOL, 'element delta sum')
                delta['largest_absolute_wake_change'] = max(ELEMENTS, key=lambda e: abs(delta[e+'_delta_wake_counts']))
                contributions.append(delta)
                for sr_a, sr_b in zip(a['surface_rows'], b['surface_rows']):
                    sd = dict(dataset=dataset, target_cl=cl, ncrit=n, element=sr_a['element'], side=sr_a['side'])
                    for field in ('transition_x_local', 'separated_arclength_fraction', 'max_Hk',
                                  'te_theta', 'te_delta_star', 'te_ue_over_vinf', 'te_Cp', 'te_cf', 'te_Hk'):
                        sd[field+'_5'], sd[field+'_10'] = sr_a[field], sr_b[field]
                        sd['delta_'+field] = sr_b[field]-sr_a[field] if sr_a[field] is not None and sr_b[field] is not None else None
                    surface_deltas.append(sd)
    return reproduced, contributions, surface_deltas, lookup


def check_published(reproduced, original):
    require(len(original) == len(reproduced) == 6, 'Missing published flap comparison')
    lookup = {(r['dataset'], float(r['target_cl'])):r for r in original}
    require(len(lookup) == 6, 'Duplicate published comparison')
    for row in reproduced:
        other = lookup[(row['dataset'], row['target_cl'])]
        for key, value in row.items():
            if key != 'dataset':
                close(value, float(other[key]), DELTA_ATOL, 'published '+key)


def save_figure(fig, stem):
    for ext in ('png', 'svg'):
        fig.savefig(stem.with_suffix('.'+ext), dpi=150)
    plt.close(fig)


def contribution_figure(rows, out):
    fig, axes = plt.subplots(2, 3, figsize=(13.5, 8), layout='constrained')
    for ax, (dataset, _, _, _, _), cl in [(axes[i,j], case, cl) for i,case in enumerate(CASES) for j,cl in enumerate(TARGETS)]:
        subset = [r for r in rows if r['dataset'] == dataset and r['target_cl'] == cl]
        for x, row in enumerate(subset):
            pos = neg = 0.
            for element, color in [('fore', '#247AA5'), ('aft', '#D99A35')]:
                value = row[element+'_delta_wake_counts']
                ax.bar(x, value, .5, bottom=pos if value >= 0 else neg, color=color,
                       label=element.title()+' wake contribution' if x == 0 else None)
                if value >= 0:
                    pos += value
                else:
                    neg += value
            ax.plot(x, row['delta_CD_counts'], 'D', color='#1D2935', ms=6,
                    label='MultiFoil total' if x == 0 else None)
        ax.axhline(subset[0]['experiment_delta_counts'], color='#1D2935', ls='--', label='Experiment')
        ax.axhline(subset[0]['mses_n9_delta_counts'], color='#8654A4', ls=':', label='Published MSES N=9')
        ax.axhline(0, color='#888', lw=.6)
        ax.set(xticks=[0,1], xticklabels=['Ncrit 9', 'Ncrit 6.25'],
               title=f"Re = {subset[0]['reynolds']:,.0f} · CL = {cl:g}",
               ylabel='CD(10°) − CD(5°), counts')
        ax.grid(axis='y', alpha=.15)
    handles, labels = axes[0,0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='outside lower center', ncols=3, fontsize=9)
    fig.suptitle('S702 | Which modeled element changes the flap drag ranking?\n'
                 'Element wake terms sum to total CD; friction is a separate diagnostic', fontsize=14)
    save_figure(fig, out/'element_drag_contributions')


def pair_figures(a, b, folder):
    sa = a['summary']
    stem = f"{sa['dataset']}_n{sa['ncrit']:g}_cl{sa['target_cl']:g}".replace('.', 'p')
    quantities = [('Cp', 'Cp'), ('cf', 'Cf'), ('ue_over_vinf', 'Ue / V∞'), ('Hk', 'Hk'),
                  ('theta', 'θ / assembly chord'), ('delta_star', 'δ* / assembly chord')]
    fig, axes = plt.subplots(6, 4, figsize=(16, 16), layout='constrained', sharex='col')
    for col, (element, side) in enumerate((e,s) for e in ELEMENTS for s in SIDES):
        for row, (field, label) in enumerate(quantities):
            ax = axes[row,col]
            for state in (a,b):
                flap = state['summary']['flap']
                color = COLORS[flap]
                values = side_data(state['surfaces'][element], side)
                layer = state['layers'][element][side]
                ax.plot(values['x_over_local_chord'], values[field], color=color, lw=1.25, label=f'{flap}° flap')
                if layer['transition'] is not None:
                    ax.axvline(layer['transition']['x_over_local_chord'], color=color, ls=':', lw=.8)
                for region in layer['separation_regions']:
                    ax.axvspan(region['start']['x_over_local_chord'], region['end']['x_over_local_chord'], color=color, alpha=.10)
            if field == 'Cp':
                ax.invert_yaxis()
            if field == 'cf':
                ax.axhline(0, color='#666', lw=.6)
            ax.set_ylabel(label)
            ax.grid(alpha=.15)
            if row == 0:
                ax.set_title(element.title()+' '+side)
            if row == 5:
                ax.set_xlabel('x / local element chord')
    axes[0,0].legend(fontsize=9)
    info = ' | '.join(f"{s['summary']['flap']}°: CL={s['summary']['CL']:.7f}, α={s['summary']['alpha_deg']:.5f}°"
                      + (' (numerical retry)' if s['summary']['numerical_retry'] else '') for s in (a,b))
    fig.suptitle(f"S702 · Re={sa['reynolds']:,.0f}, M={sa['mach']:g}, Ncrit={sa['ncrit']:g}, target CL={sa['target_cl']:g}\n"
                 f'{info}\nDotted: natural transition; shaded: Cf < 0. All curves are modeled selected states.', fontsize=12)
    save_figure(fig, folder/(stem+'_surfaces'))
    fig, axes = plt.subplots(4,2,figsize=(11,11),layout='constrained',sharex='col')
    for col, element in enumerate(ELEMENTS):
        for row,(field,label) in enumerate(quantities[2:]):
            for state in (a,b):
                flap = state['summary']['flap']
                wake = state['wakes'][element]
                axes[row,col].plot(wake['x_over_cref'], wake[field], '.-', ms=2.5, lw=1.2,
                                   color=COLORS[flap], label=f'{flap}° flap')
            axes[row,col].set_ylabel(label)
            axes[row,col].grid(alpha=.15)
            if row == 0:
                axes[row,col].set_title(element.title()+' wake')
            if row == 3:
                axes[row,col].set_xlabel('x / assembly chord')
    axes[0,0].legend()
    fig.suptitle(f"S702 wake exports · Re={sa['reynolds']:,.0f}, M={sa['mach']:g}, Ncrit={sa['ncrit']:g}, target CL={sa['target_cl']:g}\n"
                 +info+'\nSeparate non-merging modeled wakes; not an independent drag validation.', fontsize=10)
    save_figure(fig, folder/(stem+'_wakes'))


def report_text(summary, contributions):
    lines = ['# S702 archived drag mechanism diagnostic', '',
        f"Audited **{summary['coverage']['selected_states']}/42 selected states** at all 21 operating conditions, including "
        'all 16 primary states. All six archived flap-delta rows are reproduced. No new flow solves or model changes.', '',
        '## Element contributions at equal lift', '',
        'CD(10°) − CD(5°), in counts (1 count = 0.0001 CD). The element columns are modeled wake contributions; '
        'experimental element drag is unavailable. The error column is total modeled delta minus experimental delta.', '',
        '| Re | CL | Ncrit | Fore | Aft | Total | Experiment | Delta error |',
        '| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |']
    for r in contributions:
        lines.append(f"| {r['reynolds']:,.0f} | {r['target_cl']:g} | {r['ncrit']:g} | "
            f"{r['fore_delta_wake_counts']:+.2f} | {r['aft_delta_wake_counts']:+.2f} | {r['delta_CD_counts']:+.2f} | "
            f"{r['experiment_delta_counts']:+.2f} | {r['delta_error_counts']:+.2f} |")
    lines += ['', '![Element contributions](element_drag_contributions.png)', '',
        '## Audit and conventions', '',
        f"Maximum selected lift error: {summary['audit']['maximum_lift_error']:.3g} (tolerance 2e-5). "
        f"Maximum residual: {summary['audit']['maximum_residual']:.3g} (tolerance 1e-7). "
        f"Maximum total-CD minus element-wake sum: {summary['audit']['maximum_cd_sum_error']:.3g}. "
        'Selected records, polar rows and Newton residual history agree. Surface extrema, transition cells and '
        'separation endpoints were checked directly against the corresponding exports.', '',
        'The 42 states include 12 legacy reused Ncrit-9 calibration states and the two separately recovered '
        'reserved Ncrit-6.25 targets. Reserved coverage before retries remains 7/9; the original failures and '
        'the unsuccessful cold CL=0.8 retry remain archived. The first successful cold-then-warm initialization '
        'is retained, with no drag-based path selection. This report does not measure initialization repeatability.', '',
        'Friction is a diagnostic already represented within the wake drag estimate; **do not add it to CD**. '
        'There is no validated surface-by-surface CD decomposition. Surface TE thickness changes and separation '
        'are associations, not additive shares of the total drag error. Hk is kinematic shape factor. '
        'Thickness is normalized by assembly chord. A null transition means no natural transition before the TE. '
        'No flap-0/CL-1.04 pressure taps are presented as measurements of these states.', '',
        '## Files and reproduction', '',
        '- `selected_states.csv`: all 42 states, element forces, residuals, timing and actual selected export paths.',
        '- `surface_summary.csv`, `separation_intervals.csv`, `surface_pair_deltas.csv`: four surfaces per state and paired changes.',
        '- `wake_summary.csv`: first/last exported wake stations, with complete evolution in paired plots.',
        '- `element_drag_deltas.csv`: both Ncrit scenarios, all six conditions; friction deltas kept separate.',
        '- `s702_flap_5_to_10.csv`: reproduced six-row published comparison, including published MSES.',
        '- `comparison_points.csv`: all 21 conditions and four source series (84 rows).',
        '- `state_provenance.json`, `input_manifest.json`: executed settings, initialization records, versions and consumed hashes.',
        '- `pairs/`: paired surface and wake figures for all 12 scenarios, including CL=0.5 controls.', '',
        'Run from the project root with an unused output directory:', '',
        '```bash', 'OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 validation/report_s702_drag_mechanism.py --out validation/results/s702_drag_mechanism_recheck',
        '```', '',
        'The original archive manifest and per-study source guards are checked for consumed inputs. '
        'The output records original solver environments separately from this report environment. '
        'Python/NumPy/SciPy versions may differ because no solver is executed. '
        'Existing output directories are refused to preserve completed reports.', '',
        'See `FINDINGS.md` for the interpretation and next discriminating task. '
        'This is an implementation/data audit and model diagnosis, not new independent aerodynamic validation.', '']
    return '\n'.join(lines)


def build_report(out):
    require(not out.exists(), 'Output directory exists; use a new directory to preserve prior reports')
    inputs = Inputs()
    analysis = RESULTS/'s702_ncrit_analysis'
    for name in ('points.csv', 'holdout_after_retries.csv', 'holdout_points.csv', 'local_retry_points.csv',
                 'frozen_selection.json', 'summary.json', 'transfer_summary.json'):
        inputs.add(analysis/name)
    frozen = inputs.json(analysis/'frozen_selection.json')
    require(sha(analysis/'summary.json') == frozen['calibration_summary_sha256'], 'Frozen summary changed')
    require(sha(DATA/'re700k_holdout/source.json') == frozen['holdout_reference_sha256'], 'Reserved source changed')
    refs = from_analysis(analysis)
    specs = {}
    states = []
    for point in refs:
        if not point['series'].startswith('multifoil'):
            continue
        path = ROOT/point['result_path']
        folder = path.parents[2]
        if folder not in specs:
            spec, _, _ = load_study(folder, DATA/'re700k_holdout' if point['dataset'] == 'reserved' else DATA)
            inputs.add(folder/'study.json')
            inputs.sources(spec, folder)
            specs[folder] = spec
        state = load_state(point, specs[folder], inputs)
        states.append(state)
    require(len(states) == 42, 'Incomplete selected state coverage')
    reproduced, contributions, surface_deltas, lookup = pair_deltas(states, refs)
    check_published(reproduced, inputs.csv(RESULTS/'s702_flap_effect/s702_flap_5_to_10.csv'))
    summaries = [s['summary'] for s in states]
    summary = dict(schema_version=1, new_flow_solves=0, numerical_core_changed=False,
        coverage=dict(operating_conditions=21, selected_states=len(states), primary_states=sum(s['primary'] for s in summaries),
                      numerical_retries=sum(s['numerical_retry'] for s in summaries), reused_states=sum(s['reused'] for s in summaries),
                      original_reserved_selected_coverage='7/9', after_retries_reserved_selected_coverage='9/9'),
        audit=dict(published_delta_rows_reproduced=len(reproduced), published_delta_atol_counts=DELTA_ATOL,
            force_atol=FORCE_ATOL, maximum_lift_error=max(abs(s['lift_error']) for s in summaries),
            maximum_residual=max(s['residual'] for s in summaries),
            maximum_cd_sum_error=max(abs(s['cd_sum_error']) for s in summaries)),
        report_environment=dict(python=platform.python_version(), **{name:importlib.metadata.version(name)
            for name in ('numpy', 'scipy', 'matplotlib')}),
        qualification='Archived model associations; experimental element/surface drag unavailable; no validated contour gradients')
    # All scientific data are validated before any new report directory is made.
    out.mkdir(parents=True)
    (out/'pairs').mkdir()
    csv_write(out/'selected_states.csv', summaries)
    csv_write(out/'comparison_points.csv', refs)
    csv_write(out/'s702_flap_5_to_10.csv', reproduced)
    csv_write(out/'element_drag_deltas.csv', contributions)
    csv_write(out/'surface_pair_deltas.csv', surface_deltas)
    for name, key in [('surface_summary.csv', 'surface_rows'), ('separation_intervals.csv', 'intervals'), ('wake_summary.csv', 'wake_rows')]:
        csv_write(out/name, [row for s in states for row in s[key]])
    write_json(out/'state_provenance.json', [s['provenance'] for s in states])
    contribution_figure(contributions, out)
    for dataset, _, _, _, _ in CASES:
        for cl in TARGETS:
            for n in NCRITS:
                pair_figures(lookup[(dataset,n,cl,5)], lookup[(dataset,n,cl,10)], out/'pairs')
    summary['status'] = 'complete'
    write_json(out/'summary.json', summary)
    (out/'README.md').write_text(report_text(summary, contributions))
    for name in ('validation/report_s702_drag_mechanism.py', 'validation/S702_DRAG_MECHANISM_PROTOCOL.md',
                 'validation/report_s702_ncrit_sweep.py', 'validation/plot_s702_flap_effect.py',
                 'validation/s702_sweep_metrics.py', 'handover/archive_manifest.json'):
        inputs.add(ROOT/name)
    write_json(out/'input_manifest.json', dict(schema_version=1, files=inputs.files))
    return summary


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out', type=Path, default=RESULTS/'s702_drag_mechanism')
    args = ap.parse_args()
    summary = build_report(args.out.resolve())
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
