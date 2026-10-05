#!/usr/bin/env python3
"""Archived S414 pressure budgets; no solves, extrapolation, or fitted closure."""
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
PRIOR = ROOT/'validation/results/s414_incidence_analysis_20261003'
STUDY = ROOT/'validation/results/s414_incidence_20261003'
QUAL = ROOT/'validation/results/s414_geometry_qualification_20261002'
OUT = ROOT/'validation/results/s414_pressure_budget_20261003'
SURFACES = ('fore_upper', 'fore_lower', 'aft_upper', 'aft_lower')
REGIONS = {'fore': (('nose', 0., .06), ('middle', .06, .60), ('tail_slot', .60, 1.1)),
           'aft': (('nose', 0., .80), ('front_recovery', .80, .86),
                   ('bubble_recovery', .86, .94), ('tail', .94, 1.1))}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


class Inputs:
    def __init__(self):
        self.files = {}

    def add(self, path):
        path = Path(path)
        self.files[str(path.relative_to(ROOT))] = sha(path)
        return path

    def json(self, path):
        return json.loads(self.add(path).read_text())

    def csv(self, path):
        with self.add(path).open() as stream:
            return list(csv.DictReader(stream))


def write(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')


def csv_write(path, rows):
    require(bool(rows), 'Empty export')
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def curve(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    require(x.ndim == 1 and x.shape == y.shape and len(x) >= 2
            and np.isfinite(x).all() and np.isfinite(y).all()
            and np.all(np.diff(x) > 0), 'Invalid or nonmonotone curve')
    return x, y


def interp(x, y, q):
    x, y = curve(x, y); q = np.asarray(q, float)
    require(np.isfinite(q).all() and np.all((x[0] <= q) & (q <= x[-1])),
            'Extrapolation forbidden')
    return np.interp(q, x, y)


def branch(surface, side):
    require(side in ('upper', 'lower'), 'Invalid side')
    n = int(np.argmin(surface['x_over_cref']))
    b = surface[:n+1][::-1] if side == 'upper' else surface[n:]
    curve(b['x_over_cref'], b['Cp'])
    return b


def loads(x, y, cp, closed=False):
    """CCW contour: CN=integral Cp dx; CC=-integral Cp dy. Coordinates / cref."""
    xy = np.column_stack((x, y)); cp = np.asarray(cp, float)
    require(xy.shape == (len(cp), 2) and len(cp) >= 2 and
            np.isfinite(xy).all() and np.isfinite(cp).all(), 'Invalid contour')
    if closed:
        xy = np.vstack((xy, xy[0])); cp = np.r_[cp, cp[0]]
    force = np.sum(.5*(cp[:-1]+cp[1:])[:, None]*np.diff(xy, axis=0), axis=0)
    return np.array([force[0], -force[1]])


def integrate(gx, gy, px, cp, lo, hi, side):
    """Integrate linear-in-x Cp along the same polygon, splitting both knot sets."""
    gx, gy = curve(gx, gy); px, cp = curve(px, cp)
    require(side in ('upper', 'lower'), 'Invalid side')
    require(max(gx[0], px[0]) <= lo <= hi <= min(gx[-1], px[-1]),
            'Unsupported integration interval')
    if lo == hi:
        return np.zeros(2)
    x = np.unique(np.r_[lo, hi, gx[(lo < gx) & (gx < hi)], px[(lo < px) & (px < hi)]])
    answer = loads(x, interp(gx, gy, x), interp(px, cp, x))
    return answer * (-1 if side == 'upper' else 1)


def lift(force, alpha):
    a = np.deg2rad(alpha)
    return float(force[0]*np.cos(a)-force[1]*np.sin(a))


def experimental_normal(cl, cd, alpha):
    """2010 report's definition CL=CN/cos(alpha)-CD_wake*tan(alpha)."""
    a = np.deg2rad(alpha)
    return float(cl*np.cos(a)+cd*np.sin(a))


def signed_bound(weight, cp_min, cp_max):
    require(np.isfinite([weight, cp_min, cp_max]).all() and cp_min <= cp_max,
            'Invalid conditional Cp range')
    return sorted((weight*cp_min, weight*cp_max))


def physical_taps(rows):
    kept = [r for r in rows if r['included'] == 'True']
    kept.sort(key=lambda r: int(r['table_index']))
    require(len(kept) >= 2 and len({r['table_index'] for r in kept}) == len(kept),
            'Repeated or missing physical taps')
    require(all(int(b['table_index']) == int(a['table_index'])+1
                for a, b in zip(kept, kept[1:])), 'Interior missing tap requires a separate gap budget')
    x, cp = curve([float(r['x_over_c']) for r in kept], [float(r['experiment_Cp']) for r in kept])
    return kept, x, cp


def verify_prior(inputs):
    receipt = inputs.json(PRIOR/'verification.json')
    require(receipt['complete'] and receipt['selected_states_audited'] == 4, 'Prior audit incomplete')
    for name, digest in receipt['output_sha256'].items():
        require(sha(inputs.add(PRIOR/name)) == digest, 'Prior output changed: '+name)
    for name, digest in receipt['code_sha256'].items():
        require(sha(inputs.add(ROOT/name)) == digest, 'Prior verification code changed: '+name)
    manifest = inputs.json(PRIOR/'input_manifest.json')['files']
    for name, digest in manifest.items():
        require(sha(inputs.add(ROOT/name)) == digest, 'Prior input changed: '+name)
    return receipt


def main():
    OUT.mkdir(exist_ok=True)
    inputs = Inputs(); prior = verify_prior(inputs)
    targets = inputs.json(PRIOR/'summary.json')['targets']
    taps = inputs.csv(PRIOR/'pressure_at_taps.csv')
    old_taps = inputs.csv(ROOT/'validation/results/s414_lower_bucket_analysis_20261003/pressure_at_taps.csv')
    for name in ('source_manifest.json', 'sources/s414_report_2010.pdf', 'sources/s414_report_2010.txt'):
        inputs.add(QUAL/name)
    source = inputs.json(QUAL/'source_manifest.json')
    original = next(r for r in source['sources'] if r['name'] == 's414_report_2010.pdf')
    require(sha(QUAL/'sources/s414_report_2010.pdf') == original['sha256'], 'Primary source changed')
    surfaces_out, intervals, missing, budgets, dense_out, mses_out = [], [], [], [], [], []
    max_force_error = 0.
    for t in targets:
        job, case = t['job'], t['case']; folder = STUDY/job
        result = inputs.json(folder/'result.json'); i = result['selected_flow_index']
        require(result['feasible'], 'Failed state cannot enter budget')
        surfaces = {}
        for e in ('fore', 'aft'):
            s = np.genfromtxt(inputs.add(folder/f'flow/surface_{i:03d}_{e}.csv'),
                              delimiter=',', names=True, dtype=None, encoding='utf8')
            surfaces[e] = s
            f = loads(s['x_over_cref'], s['y_over_cref'], s['Cp'], closed=True)
            err = lift(f, t['alpha_deg'])-result['element_forces'][e]['CL']
            max_force_error = max(max_force_error, abs(err))
            require(abs(err) < 1e-12, 'Dense pressure integral fails archived element lift')
            dense_out.append(dict(job=job, element=e, CN=float(f[0]), CC_pressure=float(f[1]),
                CL_pressure=lift(f, t['alpha_deg']), archived_CL=result['element_forces'][e]['CL'], error=err))
        for coord in ('table', 'plotted'):
            base = dict(job=job, case=case, comparison_basis=t['comparison_basis'], coordinate=coord)
            exp_total = np.zeros(2); model_total = np.zeros(2); dense_supported = np.zeros(2)
            omitted_total = np.zeros(2); cn_bounds = np.zeros(2); cp_coverage = 0
            for surface in SURFACES:
                e, side = surface.split('_'); b = branch(surfaces[e], side)
                gx, gy, dense_cp = b['x_over_cref'], b['y_over_cref'], b['Cp']
                rows = [r for r in taps if r['job'] == job and r['coordinate'] == coord and r['surface'] == surface]
                supported, tx, ex = physical_taps(rows); cp_coverage += len(tx)
                mc = interp(gx, dense_cp, tx)
                require(np.max(np.abs(mc-np.array([float(r['model_Cp']) for r in supported]))) < 1e-13,
                        'Tap interpolation differs from previous audit')
                fe = integrate(gx, gy, tx, ex, tx[0], tx[-1], side)
                fm = integrate(gx, gy, tx, mc, tx[0], tx[-1], side)
                fd = integrate(gx, gy, gx, dense_cp, tx[0], tx[-1], side)
                exp_total += fe; model_total += fm; dense_supported += fd
                surfaces_out.append(dict(**base, surface=surface, taps=len(tx), x_first=float(tx[0]), x_last=float(tx[-1]),
                    experiment_CN=float(fe[0]), model_CN=float(fm[0]), delta_CN=float(fm[0]-fe[0]),
                    experiment_CC_pressure=float(fe[1]), model_CC_pressure=float(fm[1]),
                    delta_CL_common_angle=lift(fm-fe, t['reference_caption_alpha']),
                    model_dense_supported_CN=float(fd[0]), model_sampling_correction_CN=float(fd[0]-fm[0])))
                # Use all distinct physical readings for an explicitly conditional unmeasured-region range.
                physical = [float(r['experiment_Cp']) for r in rows if r['table_index'] != '']
                cp_min, cp_max = min(physical), max(physical)
                sign = -1 if side == 'upper' else 1
                for region, lo, hi in [('uncovered_nose', gx[0], tx[0]), ('uncovered_tail', tx[-1], gx[-1])]:
                    f = integrate(gx, gy, gx, dense_cp, lo, hi, side); omitted_total += f
                    bounds = signed_bound(sign*(hi-lo), cp_min, cp_max); cn_bounds += bounds
                    missing.append(dict(**base, surface=surface, region=region, left_x=float(lo), right_x=float(hi),
                        model_CN=float(f[0]), model_CC_pressure=float(f[1]), assumed_Cp_min=cp_min, assumed_Cp_max=cp_max,
                        conditional_experiment_CN_low=bounds[0], conditional_experiment_CN_high=bounds[1]))
                # Split at physical taps and predeclared region edges for cumulative and regional diagnostics.
                cuts = np.unique(np.r_[tx, [v for _, lo, hi in REGIONS[e] for v in (lo, hi) if tx[0] < v < tx[-1]]])
                for lo, hi in zip(cuts[:-1], cuts[1:]):
                    region = next(name for name, low, high in REGIONS[e] if low <= (lo+hi)/2 < high)
                    ef = integrate(gx, gy, tx, ex, lo, hi, side)
                    mf = integrate(gx, gy, tx, mc, lo, hi, side)
                    intervals.append(dict(**base, surface=surface, region=region, left_x=float(lo), right_x=float(hi),
                        experiment_CN=float(ef[0]), model_CN=float(mf[0]), delta_CN=float(mf[0]-ef[0]),
                        experiment_CC_pressure=float(ef[1]), model_CC_pressure=float(mf[1])))
                if case == 'cl028' and t['comparison_basis'] == 'fixed_lift':
                    mt = [r for r in old_taps if r['case'] == case and r['surface'] == surface and
                          r['coordinate'] == coord and r['method'] == 'mses' and r['included'] == 'True']
                    mt.sort(key=lambda r: int(r['table_index']))
                    require([r['table_index'] for r in mt] == [r['table_index'] for r in supported], 'MSES tap coverage differs')
                    require(np.array_equal(tx, np.array([float(r['x_over_c']) for r in mt])), 'MSES tap coordinates differ')
                    ms = integrate(gx, gy, tx, [float(r['model_Cp']) for r in mt], tx[0], tx[-1], side)
                    mses_out.append(dict(case=case, coordinate=coord, surface=surface, taps=len(tx),
                        experiment_CN=float(fe[0]), multifoil_matched_lift_CN=float(fm[0]),
                        published_mses_CN=float(ms[0]), mses_minus_experiment_CN=float(ms[0]-fe[0]),
                        mses_basis='Published nominal CL 0.28; no measured-angle MSES state'))
            # Include the blunt trailing-edge closure explicitly, not as a measured tap.
            for e, s in surfaces.items():
                f = loads(s['x_over_cref'][[-1, 0]], s['y_over_cref'][[-1, 0]], s['Cp'][[-1, 0]])
                omitted_total += f
                phys = [float(r['experiment_Cp']) for r in taps if r['job'] == job and r['coordinate'] == coord
                        and r['surface'].startswith(e+'_') and r['table_index'] != '']
                weight = s['x_over_cref'][0]-s['x_over_cref'][-1]
                bounds = signed_bound(weight, min(phys), max(phys)); cn_bounds += bounds
                missing.append(dict(**base, surface=e+'_TE', region='unmeasured_TE_closure',
                    left_x=float(s['x_over_cref'][-1]), right_x=float(s['x_over_cref'][0]),
                    model_CN=float(f[0]), model_CC_pressure=float(f[1]), assumed_Cp_min=min(phys), assumed_Cp_max=max(phys),
                    conditional_experiment_CN_low=bounds[0], conditional_experiment_CN_high=bounds[1]))
            full = np.array([sum(r[key] for r in dense_out if r['job'] == job) for key in ('CN', 'CC_pressure')])
            require(np.max(np.abs(full-dense_supported-omitted_total)) < 1e-12, 'Coverage partition does not close')
            require(cp_coverage == 87, 'Physical tap coverage changed')
            # Use the internally consistent appendix row; caption rounding is a separate diagnostic.
            exp_cn = experimental_normal(t['reference_appendix_CL'], t['reference_appendix_CD'], t['reference_appendix_alpha'])
            caption_cn = experimental_normal(t['reference_caption_CL'], t['reference_appendix_CD'], t['reference_caption_alpha'])
            a = np.deg2rad(t['alpha_deg'])
            source_style_cl = full[0]/np.cos(a)-t['CD_context']*np.tan(a)
            expected_gap = full[0]-exp_cn
            decomposition = (model_total[0]-exp_total[0]) + (full[0]-model_total[0]) + (exp_total[0]-exp_cn)
            require(abs(expected_gap-decomposition) < 1e-12, 'Normal-force budget does not close')
            budgets.append(dict(**base, taps=cp_coverage, alpha=t['alpha_deg'], model_CL=t['CL'],
                reported_appendix_CL=t['reference_appendix_CL'], reported_appendix_CN=exp_cn,
                caption_approx_CN=caption_cn, experiment_partial_CN=float(exp_total[0]),
                model_partial_CN=float(model_total[0]), delta_partial_CN=float(model_total[0]-exp_total[0]),
                model_full_CN=float(full[0]), model_sampling_correction_CN=float(dense_supported[0]-model_total[0]),
                model_uncovered_CN=float(omitted_total[0]), experimental_unresolved_CN=float(exp_cn-exp_total[0]),
                full_model_minus_reported_CN=float(expected_gap),
                conditional_missing_CN_low=float(cn_bounds[0]), conditional_missing_CN_high=float(cn_bounds[1]),
                conditional_completed_exp_CN_low=float(exp_total[0]+cn_bounds[0]),
                conditional_completed_exp_CN_high=float(exp_total[0]+cn_bounds[1]),
                model_CL_using_report_formula=float(source_style_cl),
                lift_definition_change=float(source_style_cl-t['CL']),
                experiment_partial_pressure_CL_at_caption=lift(exp_total, t['reference_caption_alpha']),
                delta_partial_pressure_CL_at_caption=lift(model_total-exp_total, t['reference_caption_alpha'])))
    for name, rows in [('surface_budgets', surfaces_out), ('interval_budgets', intervals),
                       ('uncovered_regions', missing), ('budgets', budgets), ('dense_force_checks', dense_out),
                       ('published_mses_partial', mses_out)]:
        csv_write(OUT/(name+'.csv'), rows)
    summary = dict(complete=True, cases=2, states=4, coordinate_conventions=['table', 'plotted'],
        supported_taps_per_case=87, dense_element_force_checks=8, maximum_dense_CL_error=max_force_error,
        prior_audited_states=prior['selected_states_audited'], new_nonlinear_flow_solves=0,
        total_force_validation=False, missing_region_bounds='Conditional on per-surface observed Cp range; not measurement uncertainty bounds',
        experimental_reduction='CN=CL*cos(alpha)+CD_wake*sin(alpha), consistent appendix rows',
        mses_basis='Published CL 0.28 only; partial identical-tap normal force; no new MSES run', budgets=budgets)
    write(OUT/'summary.json', summary)
    plot(budgets, surfaces_out, intervals, mses_out)
    inputs.add(Path(__file__)); inputs.add(ROOT/'test_s414_pressure_budget.py')
    write(OUT/'input_manifest.json', dict(created_utc=datetime.now(timezone.utc).isoformat(),
        python=platform.python_version(), numpy=np.__version__, matplotlib=matplotlib.__version__, files=inputs.files))
    for row in budgets:
        if row['coordinate'] == 'table':
            print(row['job'], 'partial ΔCN', row['delta_partial_CN'], 'full ΔCN', row['full_model_minus_reported_CN'], flush=True)


def save(fig, name):
    for suffix in ('png', 'pdf', 'svg'):
        fig.savefig(OUT/(name+'.'+suffix), dpi=170, bbox_inches='tight')
    plt.close(fig)


def plot(budgets, surfaces, intervals, mses):
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
    colors = {'fixed_lift': '#2364aa', 'fixed_alpha': '#d86a00'}
    names = {'fixed_lift': 'Matched lift', 'fixed_alpha': 'Measured angle'}
    cases = ('cl0056', 'cl028')
    titles = ('Experiment CL 0.056, α −3.07°', 'Experiment CL 0.284, α −1.04°')
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), layout='constrained')
    for ax, case, title in zip(axes, cases, titles):
        for offset, basis in [(-.19, 'fixed_lift'), (.19, 'fixed_alpha')]:
            rows = [r for r in surfaces if r['case'] == case and r['coordinate'] == 'table' and r['comparison_basis'] == basis]
            values = [next(r['delta_CN'] for r in rows if r['surface'] == s) for s in SURFACES]
            bars = ax.bar(np.arange(4)+offset, values, width=.36, color=colors[basis], label=names[basis])
            ax.bar_label(bars, fmt='%+.4f', fontsize=8, padding=3)
        ax.set_xticks(np.arange(4), ['Fore upper', 'Fore lower', 'Aft upper', 'Aft lower'], rotation=18)
        ax.axhline(0, color='k', lw=.8); ax.set(title=title, ylabel='Partial ΔCN: model − experiment', ylim=(-.051,.064))
        ax.grid(axis='y', alpha=.2); ax.legend(fontsize=8)
    fig.suptitle('Identical supported taps reveal pressure-load cancellation at matched lift\n87 physical taps per case; table coordinates; uncovered nose/TE contributions are excluded')
    save(fig, 'surface_load_errors')
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), layout='constrained')
    for col, (case, title) in enumerate(zip(cases, titles)):
        for row, e in enumerate(('fore', 'aft')):
            ax = axes[row, col]
            for basis in ('fixed_lift', 'fixed_alpha'):
                for side, style in [('upper', '-'), ('lower', '--')]:
                    values = [r for r in intervals if r['case'] == case and r['comparison_basis'] == basis
                              and r['coordinate'] == 'table' and r['surface'] == e+'_'+side]
                    x = [values[0]['left_x']]+[r['right_x'] for r in values]
                    y = np.r_[0., np.cumsum([r['delta_CN'] for r in values])]
                    ax.plot(x, y, style, color=colors[basis], label=names[basis]+': '+side)
            ax.axhline(0, color='k', lw=.8); ax.grid(alpha=.2)
            ax.set(title=title+' · '+e, xlabel='x/c', ylabel='Cumulative partial ΔCN')
            ax.legend(fontsize=7)
    fig.suptitle('Where pressure differences accumulate\nLinear tap interpolation; each surface begins at its first supported tap, with no endpoint completion')
    save(fig, 'cumulative_pressure_load')
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), layout='constrained')
    for ax, case, title in zip(axes, cases, titles):
        r = next(r for r in budgets if r['case'] == case and r['comparison_basis'] == 'fixed_alpha' and r['coordinate'] == 'table')
        values = [r['delta_partial_CN'], r['model_sampling_correction_CN'], r['model_uncovered_CN'], -r['experimental_unresolved_CN']]
        start = 0.
        for i, val in enumerate(values):
            ax.bar(i, val, bottom=start, color=('#d86a00' if i == 0 else '#777777'), width=.65)
            ax.annotate(f'{val:+.5f}', (i, start+val), xytext=(0,5 if val >= 0 else -13), textcoords='offset points', ha='center', fontsize=8)
            start += val
        ax.bar(4, r['full_model_minus_reported_CN'], color='#2364aa', width=.65)
        ax.text(4, start+.002, f'{start:+.5f}', ha='center', fontsize=8)
        ax.set_xticks(range(5), ['Same-tap\npressure gap', 'Model sampling\ncorrection', 'Model uncovered\nnose / TE', 'Experiment\nclosure remainder', 'Full CN\ndifference'], rotation=15)
        ax.set(title=title, ylabel='Normal-force coefficient difference', ylim=(-.003,.09)); ax.grid(axis='y', alpha=.2)
    fig.suptitle('Measured-angle excess normal force: exact arithmetic budget\nThe experimental remainder includes missing regions, digitization, quadrature and reduction; it is not measured endpoint loading')
    save(fig, 'normal_force_budget')
    fig, ax = plt.subplots(figsize=(9, 4.5), layout='constrained')
    for j, (key, label, color) in enumerate([('multifoil_matched_lift_CN', 'MultiFoil at matched lift', '#2364aa'), ('published_mses_CN', 'Published MSES nominal CL 0.28', '#168a6e')]):
        selected = [r for r in mses if r['coordinate'] == 'table']
        values = [next(r[key]-r['experiment_CN'] for r in selected if r['surface'] == s) for s in SURFACES]
        bars = ax.bar(np.arange(4)+(j-.5)*.36, values, width=.34, label=label, color=color)
        ax.bar_label(bars, fmt='%+.4f', fontsize=8, padding=3)
    ax.axhline(0, color='k', lw=.8); ax.set_xticks(range(4), [s.replace('_', ' ') for s in SURFACES])
    ax.set(ylabel='Partial ΔCN: calculation − experiment', ylim=(-.05,.035), title='Published MSES context on the identical supported taps\nExperiment CL 0.284; nominal theoretical lift is rounded; no measured-angle MSES run')
    ax.legend(fontsize=9); ax.grid(axis='y', alpha=.2); save(fig, 'published_mses_loads')


if __name__ == '__main__':
    main()
