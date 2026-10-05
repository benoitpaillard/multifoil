#!/usr/bin/env python3
"""Audit saved S702 initialization paths and plot their differences; no solves."""
import argparse
from itertools import combinations
from pathlib import Path
import platform

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from report_s702_drag_mechanism import (
    Inputs, ROOT, read, sha, require, close, load_state, write_json, relative,
)
from report_s702_ncrit_sweep import csv_write

STUDY_HASH = 'b8248d34f040518483943295d56fc90ef3fe3369fecb443b630026789fa6adee'
PATHS = ('cold', 'increasing', 'decreasing', 'warm_replay')
COLORS = dict(zip(PATHS, ('#277da8', '#e07a27', '#32936f', '#9355a2')))


def state(path, spec, inputs, label):
    result = inputs.json(path)
    point = dict(result_path=relative(path), achieved_CL=result['CL'], CD=result['CD'],
                 ncrit=result['ncrit'], flap_deg=result['flap'], target_cl=result['target_cl'],
                 reynolds=spec['reynolds'], mach=spec['mach'], dataset=label,
                 numerical_retry=len(result['attempts']) > 1, retry_of='')
    loaded = load_state(point, spec, inputs)
    arrays = {}
    for element in ('fore', 'aft'):
        with np.load(inputs.add(path.parent/f'{element}_state.npz'), allow_pickle=False) as z:
            arrays[element] = {k: z[k].copy() for k in z.files}
        a = arrays[element]
        require(all(np.isfinite(v).all() for v in a.values()), 'Nonfinite selected state')
        surface, wake = loaded['surfaces'][element], loaded['wakes'][element]
        count = a['xy'].shape[1]
        for row, field in ((0, 'theta'), (1, 'delta_star')):
            require(np.max(abs(a['U'][row, :count][::-1]-surface[field])) < 1e-12,
                    'Saved Newton state differs from selected surface')
            require(np.max(abs(a['U'][row, count:]-wake[field])) < 1e-12,
                    'Saved Newton state differs from selected wake')
        require(np.array_equal(a['turbulent'][:count][::-1], surface['turbulent']),
                'Saved turbulence flags differ from selected surface')
    loaded.update(arrays=arrays, result=result, label=label)
    return loaded


def compare_states(a, b):
    """Pointwise differences at corresponding panel/wake indices, not interpolation."""
    rows = []
    for element in ('fore', 'aft'):
        aa, bb = a['arrays'][element], b['arrays'][element]
        for field in aa:
            require(aa[field].shape == bb[field].shape, 'State topology dimensions changed')
            pieces = range(4) if field == 'U' else (None,)
            for k in pieces:
                x, y = (aa[field][k], bb[field][k]) if k is not None else (aa[field], bb[field])
                rows.append(dict(element=element, field=field if k is None else f'U[{k}]',
                    max_abs_difference=float(np.max(abs(x-y))),
                    exactly_equal=bool(np.array_equal(x, y))))
    return rows


def plot(out, rows, states, references):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), sharey=True)
    for ax, flap in zip(axes, (5, 10)):
        for offset, path in zip((-.024, -.008, .008, .024), PATHS):
            group = [r for r in rows if r['flap'] == flap and r['path'] == path]
            good = [r for r in group if r['feasible']]
            ax.scatter([r['target_cl']+offset for r in good],
                       [r['difference_CD_counts'] for r in good], label=path,
                       color=COLORS[path], marker='o', s=38)
            for r in group:
                if not r['feasible']:
                    ax.text(.8, .93, 'cold failed', transform=ax.get_xaxis_transform(),
                            color=COLORS[path], ha='center')
        ax.axhline(0, color='0.5', lw=.7)
        ax.set(title=f'Flap {flap}°', xlabel='Target CL', xticks=[.5, .8, 1.2])
        ax.grid(alpha=.2)
    axes[0].set_ylabel('CD − archived selected CD [counts]')
    axes[1].legend(fontsize=8)
    fig.suptitle('18/19 converged; the missing cold result has no plotted force')
    fig.tight_layout()
    for ext in ('png', 'svg'): fig.savefig(out/f'drag_repeatability.{ext}', dpi=170)
    plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(11, 7))
    ref = references[10, .8]
    selected = [s for s in states if s['result']['flap'] == 10 and s['result']['target_cl'] == .8]
    for col, element in enumerate(('fore', 'aft')):
        for s in [ref]+selected:
            values = s['surfaces'][element]
            upper = values[values['side'] == 'upper']
            label = s['label']
            for ax, field in ((axes[0, col], 'Cp'), (axes[1, col], 'Hk')):
                ax.plot(upper['x_over_local_chord'], upper[field], label=label,
                        color=COLORS.get(label, 'black'), lw=1.2,
                        linestyle='--' if label == 'archive' else '-')
                ax.grid(alpha=.2)
        axes[0, col].set_title(element+' upper surface')
        axes[0, col].invert_yaxis()
        axes[1, col].set_xlabel('x / local chord')
    axes[0, 0].set_ylabel('Cp'); axes[1, 0].set_ylabel('Hk')
    axes[0, 1].legend(fontsize=8)
    fig.suptitle('Flap 10°, CL=0.8: successful warm paths overlap; cold path failed')
    fig.tight_layout()
    for ext in ('png', 'svg'): fig.savefig(out/f'difficult_state.{ext}', dpi=170)
    plt.close(fig)


def report(study, out):
    require(sha(study/'study.json') == STUDY_HASH, 'Unexpected original study manifest')
    inputs = Inputs()
    spec = inputs.json(study/'study.json')
    inputs.sources(spec, study)
    runtime = inputs.json(study/'runtime.json')
    require(runtime['versions'] == spec['expected_versions'], 'Numerical runtime mismatch')
    status = inputs.json(study/'status.json')
    summary = inputs.json(study/'summary.json')
    require(status['state'] == 'finished', 'Matrix has not finished')
    references = {}
    for ref in spec['references']:
        references[ref['flap'], ref['target_cl']] = state(ROOT/ref['result_path'], spec, inputs, 'archive')
    rows, states, failures, comparisons = [], [], [], []
    for job in spec['jobs']:
        group = inputs.json(study/job['name']/'result.json')
        require([r['target_cl'] for r in group['targets']] == job['targets'], 'Target order changed')
        for result in group['targets']:
            path = study/job['name']/('cl'+str(result['target_cl']).replace('.', 'p'))/'result.json'
            require(inputs.json(path) == result, 'Group and individual results disagree')
            row = dict(job=job['name'], path=job['path'], flap=job['flap'],
                       target_cl=result['target_cl'], feasible=result['feasible'],
                       first_attempt_feasible=result['attempts'][0]['feasible'],
                       attempts=len(result['attempts']), requested_evaluations=len(result['records']),
                       wall_seconds=result['wall_seconds'], CD=result['CD'], difference_CD_counts=None)
            if not result['feasible']:
                require(all(result[k] is None for k in ('CL', 'CD', 'CM', 'alpha_deg')),
                        'Failed target has finite reported forces')
                require('selected_flow_index' not in result, 'Failed target has selected state')
                meta = inputs.json(path.parent/'flow/metadata.json')
                failures.append(dict(**row, attempts_detail=result['attempts'], records=result['records'],
                    continuation=[d.get('continuation', []) for d in meta['diagnostics']],
                    failure_reason=result['error']))
            else:
                s = state(path, spec, inputs, job['path'])
                s['job'] = job['name']; states.append(s)
                ref = references[job['flap'], result['target_cl']]
                row['difference_CD_counts'] = (result['CD']-ref['result']['CD'])*1e4
                for k in ('CL', 'CM', 'alpha_deg'):
                    row['difference_'+k] = result[k]-ref['result'][k]
            rows.append(row)
    require(len(rows) == 19 and len(states) == summary['converged_targets'] == status['converged_targets'],
            'Coverage summary mismatch')
    require(sum(r['first_attempt_feasible'] for r in rows) == summary['first_attempt_successes'],
            'First-attempt coverage mismatch')
    require(sum(r['requested_evaluations'] for r in rows) == summary['requested_evaluations'],
            'Work summary mismatch')
    for old, new in zip(inputs.csv(study/'comparison.csv'), rows):
        require(old['job'] == new['job'], 'Comparison order mismatch')
        require((old['feasible'] == 'True') == new['feasible'], 'Comparison coverage mismatch')
        if new['feasible']: close(float(old['difference_CD_counts']), new['difference_CD_counts'])

    spans = []
    for (flap, cl), ref in references.items():
        group = [s for s in states if (s['result']['flap'], s['result']['target_cl']) == (flap, cl)]
        paths = [s for s in group if s['label'] != 'warm_replay']
        complete = len(paths) == 3
        span = (max(s['result']['CD'] for s in paths)-min(s['result']['CD'] for s in paths))*1e4
        spans.append(dict(flap=flap, target_cl=cl, complete=complete,
                          converged_paths=len(paths), required_paths=3,
                          CD_span_counts=span if complete else None,
                          observed_successful_only_CD_span_counts=span))
        for a, b in combinations([ref]+group, 2):
            for d in compare_states(a, b):
                comparisons.append(dict(flap=flap, target_cl=cl, first=a['label'], second=b['label'], **d))
    topology = [d for d in comparisons if d['field'] in ('turbulent', 'lower', 'upper')]
    require(all(d['exactly_equal'] for d in topology), 'Successful states change discrete topology')
    path_only = [d for d in comparisons if d['first'] != 'archive']
    maxima = {field: max(d['max_abs_difference'] for d in path_only if d['field'] == field)
              for field in sorted({d['field'] for d in path_only})}
    out.mkdir(parents=True, exist_ok=True)
    csv_write(out/'targets.csv', rows)
    csv_write(out/'path_spans.csv', spans)
    csv_write(out/'state_differences.csv', comparisons)
    write_json(out/'failures.json', failures)
    plot(out, rows, states, references)
    result = dict(coverage=dict(required=19, converged=len(states),
                              first_attempts=sum(r['first_attempt_feasible'] for r in rows),
                              fully_covered_conditions=sum(r['complete'] for r in spans)),
                  max_complete_CD_span_counts=max(r['CD_span_counts'] for r in spans if r['complete']),
                  max_absolute_CD_difference_from_archive_counts=max(abs(r['difference_CD_counts']) for r in rows if r['feasible']),
                  successful_paths_have_identical_turbulence_and_surface_indices=True,
                  max_pointwise_state_differences_between_successful_paths=maxima,
                  failed_jobs=[r['job'] for r in failures],
                  requested_evaluations=summary['requested_evaluations'],
                  failed_requested_evaluations=summary['failed_requested_evaluations'],
                  selected_states_audited=len(states)+len(references),
                  numerical_runtime=runtime, report_python=platform.python_version(),
                  qualification='Numerical repeatability at these six conditions, not independent aerodynamic validation; the cold CL=0.8 flap-10 condition remains incomplete.')
    write_json(out/'summary.json', result)
    inputs.add(Path(__file__))
    inputs.add(ROOT/'validation/report_s702_drag_mechanism.py')
    write_json(out/'input_manifest.json', dict(study_sha256=STUDY_HASH, files=inputs.files))
    return result


if __name__ == '__main__':
    import json
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--study', type=Path, default=ROOT/'validation/results/s702_repeatability_20260915')
    ap.add_argument('--out', type=Path, default=ROOT/'validation/results/s702_repeatability_analysis')
    args = ap.parse_args()
    print(json.dumps(report(args.study.resolve(), args.out.resolve()), indent=2))
