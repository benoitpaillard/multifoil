#!/usr/bin/env python3
"""Render the completed pilot without running the flow solver."""
from pathlib import Path
import argparse
import hashlib
import itertools
import json
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from validate import verify_data, read_csv, at_angle, write_csv
from validation.optimization_pilot import summarize

COLORS = ['#7c8793', '#0875b5', '#237b52', '#d47713', '#9b4b98', '#633345']
SETTINGS = ('coarse', 'base', 'fine', 'long_wake')
LABELS = ['180 / 120 panels', '260 / 180 panels', '340 / 240 panels',
          '260 / 180; 2c wake', '260 / 180; Ncrit 6', '340 / 240; Ncrit 6']


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=ROOT/'validation/results/optimization_pilot')
    parser.add_argument('--transition-followup', type=Path)
    parser.add_argument('--transition-retry', type=Path, action='append', default=[])
    args=parser.parse_args()
    out = args.out.resolve()
    verify_data()
    spec = json.loads((out/'study.json').read_text())
    for filename, expected in spec['sha256'].items():
        if hashlib.sha256((ROOT/filename).read_bytes()).hexdigest() != expected:
            raise ValueError(f'Code/input differs from executed study: {filename}')
    summary = summarize(out, spec)
    if summary['completed'] != summary['requested']:
        raise ValueError('Study incomplete; resume it before making the final report')
    rows = [json.loads((out/f'{g}__{s}'/'result.json').read_text()) for g, s in spec['jobs']]
    records = {(r['geometry'], r['setting']): r for r in rows}
    followup_attempts=[]
    if args.transition_retry and not args.transition_followup:
        raise ValueError('A retry must be paired with its original follow-up')
    if args.transition_followup:
        directory=args.transition_followup.resolve()
        follow=json.loads((directory/'study.json').read_text())
        for key in ('target_cl','reynolds','mach','reference_chord','moment_reference','cl_tolerance'):
            if follow[key]!=spec[key]:
                raise ValueError(f'Follow-up condition mismatch: {key}')
        for filename, expected in follow['sha256'].items():
            if hashlib.sha256((ROOT/filename).read_bytes()).hexdigest()!=expected:
                raise ValueError(f'Follow-up code/input changed: {filename}')
        for g,s in follow['jobs']:
            r=json.loads((directory/f'{g}__{s}'/'result.json').read_text())
            records[(g,s)]=r
            followup_attempts.append(dict(attempt='original',start_alpha=follow['settings'][s]['start_alpha'],
                                         alpha_bounds=follow['alpha_bounds'],**r))
    for retry_path in args.transition_retry:
        retry_dir=retry_path.resolve()
        retry=json.loads((retry_dir/'study.json').read_text())
        if retry['source_study_sha256']!=hashlib.sha256((directory/'study.json').read_bytes()).hexdigest():
            raise ValueError('Retry source does not match the selected follow-up')
        for f,h in retry['sha256'].items():
            if hashlib.sha256((ROOT/f).read_bytes()).hexdigest()!=h:
                raise ValueError(f'Retry code/input changed: {f}')
        for g,s in retry['jobs']:
            for k,v in retry['settings'][s].items():
                if k!='start_alpha' and follow['settings'][s][k]!=v:
                    raise ValueError(f'Retry changes a physical/numerical setting: {k}')
            r=json.loads((retry_dir/f'{g}__{s}'/'result.json').read_text())
            r['retry_start_alpha']=retry['settings'][s]['start_alpha']
            r['retry_alpha_bounds']=retry['alpha_bounds']
            r['retry_directory']=retry_dir.name
            if not follow['alpha_bounds'][0]<=retry['alpha_bounds'][0]<retry['alpha_bounds'][1]<=follow['alpha_bounds'][1]:
                raise ValueError('Retry bracket extends outside the original domain')
            records[(g,s)]=r
            followup_attempts.append(dict(attempt=retry_dir.name,**r))
    if followup_attempts:
        (out/'followup_attempts.json').write_text(json.dumps(followup_attempts,indent=2,allow_nan=False)+'\n')

    def value(g, s, q='CD'):
        r = records.get((g, s))
        return r[q] if r and r['feasible'] else np.nan

    def delta(g, s):
        return 1e4*(value(g, s)-value('baseline', s))

    pairs = []
    geometries = ('baseline', 'camber_minus', 'camber_plus', 'flap_down', 'flap_up')
    for first, second in itertools.combinations(geometries, 2):
        d = [1e4*(value(second, s)-value(first, s)) for s in SETTINGS]
        valid = bool(np.isfinite(d).all())
        pairs.append(dict(first=first, second=second,
            min_second_minus_first_counts=min(d) if valid else None,
            max_second_minus_first_counts=max(d) if valid else None,
            order_stable=valid and (all(v > 0 for v in d) or all(v < 0 for v in d))))
    write_csv(out/'pairwise_rankings.csv', pairs, list(pairs[0]))
    reference = json.loads((ROOT/'validation/data/nlr7301_mses_official/reference.json').read_text())
    old = at_angle(read_csv(ROOT/'validation/results/nlr7301_mses/polar.csv'), reference['alpha_deg'])
    conditions = json.loads((ROOT/'validation/results/nlr7301_mses/metadata.json').read_text())['case']
    if old is None or old['converged'] != '1' or any(conditions[k] != reference[k] for k in ('reynolds','mach','ncrit')):
        raise ValueError('Missing matching archived baseline')
    context = []
    for name, alpha, cl, cd in [('same_alpha', float(old['alpha_deg']), float(old['CL']), float(old['CD'])),
            ('same_CL', value('baseline','base','alpha_deg'), value('baseline','base','CL'), value('baseline','base'))]:
        if not np.isfinite([alpha, cl, cd]).all():
            continue
        context.append(dict(comparison=name, alpha_deg=alpha, model_CL=cl, model_CD=cd,
            reference_CL=reference['CL'], reference_CD=reference['CD'], delta_cd_counts=(cd-reference['CD'])*1e4,
            qualification='Exact MSES geometry unconfirmed; context only'))
    write_csv(out/'official_mses_context.csv', context, list(context[0]))
    budgets=[]
    for g in geometries:
        if not records[(g,'base')]['feasible']:
            continue
        diag=json.loads((out/f'{g}__base'/'flow/metadata.json').read_text())['diagnostics'][-1]
        for element, forces in diag['element_forces'].items():
            budgets.append(dict(geometry=g, element=element, CD_wake=forces['CD_wake'],
                                CD_friction=forces['CD_friction']))
    if budgets:
        write_csv(out/'base_drag_budget.csv', budgets, list(budgets[0]))
    plot(out, records, value, delta, reference)
    report(out, rows, records, summary, pairs, context, value, delta)
    print(f'Saved report, PNG, PDF and comparison tables in {out}')


def plot(out, records, value, delta, reference):
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,
        'axes.spines.top':False,'axes.spines.right':False,'savefig.facecolor':'white'})
    fig, axes = plt.subplots(2, 2, figsize=(12, 8.5))
    for ax, mode, points, key, title in [
        (axes[0,0], 'camber', ('camber_minus','baseline','camber_plus'), 'camber_amplitude', 'A  Main-element camber'),
        (axes[0,1], 'placement', ('flap_down','baseline','flap_up'), 'flap_dy', 'B  Flap vertical position')]:
        settings=(*SETTINGS,'ncrit6')+ (('ncrit6_fine',) if np.isfinite(value('baseline','ncrit6_fine')) else ())
        for i, setting in enumerate(settings):
            if mode == 'placement' and setting.startswith('ncrit6'):
                continue
            names = points
            if mode == 'camber' and setting == 'base':
                names = ('camber_minus','camber_minus_half','baseline','camber_plus_half','camber_plus')
            x = [100*records[(g, setting)][key] for g in names]
            style='s:' if setting=='ncrit6_fine' else ('o--' if setting=='ncrit6' else 'o-')
            ax.plot(x, [delta(g, setting) for g in names], style,
                    color=COLORS[i], ms=4, lw=1.5, label=LABELS[i])
        ax.axhline(0, color='#4e5863', lw=.7)
        ax.set_title(title, loc='left', fontweight='semibold')
        ax.set_ylabel('Drag change from own baseline [counts]')
        ax.set_xlabel('Camber amplitude [% local chord]' if mode=='camber' else 'Flap vertical shift [% reference chord]')
        ax.grid(alpha=.2)
    handles,labels_top=axes[0,0].get_legend_handles_labels()
    fig.legend(handles,labels_top,loc='upper center',bbox_to_anchor=(.53,.925),
               ncol=3,fontsize=8.6,frameon=False)
    labels = ['180/120','260/180','340/240','260/180\n2c wake','260/180\nstart 0°']
    for ax, q, title, ylabel in [
        (axes[1,0], 'CD', 'C  Baseline absolute drag', '$C_D$ [drag counts]')]:
        vals = [value('baseline', s, q)*(1e4 if q=='CD' else 1) for s in (*SETTINGS,'start0')]
        ax.plot(range(len(vals)), vals, 'o-', color=COLORS[1], ms=5)
        ax.set_xticks(range(len(vals)), labels, fontsize=8.5)
        ax.set_title(title, loc='left', fontweight='semibold')
        ax.set_ylabel(ylabel)
        ax.grid(alpha=.2)
        if q=='CD':
            ax.axhline(reference['CD']*1e4, color=COLORS[3], ls='--', lw=1.2,
                       label='MSES at same lift; geometry unconfirmed')
            ax.legend(fontsize=8.3, frameon=False)
    ax=axes[1,1]
    derivative_settings=['coarse','base','fine','ncrit6']
    if ('camber_plus','ncrit6_fine') in records:
        derivative_settings.append('ncrit6_fine')
    derivatives=[(value('camber_plus',s)-value('camber_minus',s))/.002 for s in derivative_settings]
    derivative_labels=['180/120\nNcrit 9','260/180\nNcrit 9','340/240\nNcrit 9','260/180\nNcrit 6','340/240\nNcrit 6']
    ax.bar(range(len(derivatives)),derivatives,color=[COLORS[1] if not s.startswith('ncrit6') else COLORS[4] for s in derivative_settings],width=.6)
    ax.axhline(0,color='#4e5863',lw=.8)
    ax.set_xticks(range(len(derivatives)),derivative_labels[:len(derivatives)],fontsize=8.5)
    ax.set_title('D  Camber derivative vs transition',loc='left',fontweight='semibold')
    ax.set_ylabel('Central difference dCD/dp; p = camber/chord')
    ax.grid(axis='y',alpha=.2)
    ax.set_axisbelow(True)
    fig.suptitle('MultiFoil: first optimization-readiness pilot', x=.08, y=.99, ha='left', fontsize=18, fontweight='semibold')
    fig.text(.08,.945,'NLR7301  |  fixed CL = 2.484  |  Re = 2.51 million  |  M = 0.185',color='#505a66')
    fig.subplots_adjust(left=.08,right=.98,top=.82,bottom=.135,hspace=.48,wspace=.25)
    fig.text(.08,.071,'One drag count = 0.0001 CD. Numerical sensitivity checks do not establish experimental geometry rankings.',fontsize=8.7,color='#505a66')
    fig.text(.08,.046,'Independent non-merging wakes. Ncrit 6 tests natural-transition sensitivity; it is not a forced-trip model.',fontsize=8.7,color='#505a66')
    if ('baseline','ncrit6_fine') in records and 'retry_start_alpha' in records[('baseline','ncrit6_fine')]:
        fig.text(.08,.021,'Fine Ncrit 6 baseline uses a narrower-bracket retry; the original procedure failed at 4°. See report for all attempts.',fontsize=8.7,color='#505a66')
    for suffix in ('png','pdf'):
        fig.savefig(out/f'optimization_pilot.{suffix}',dpi=180)
    plt.close(fig)


def report(out, rows, records, summary, pairs, context, value, delta):
    lines = ['# First optimization-readiness pilot','',
        '**Numerical verification of design screening; independent geometry rankings remain unvalidated.**','',
        f"Completed {summary['completed']}/{summary['requested']} combinations; {summary['feasible']} satisfy the declared checks.",'',
        '![Pilot results](optimization_pilot.png)','',
        'Methods, assumptions and reproduction commands: [OPTIMIZATION_PILOT.md](../../OPTIMIZATION_PILOT.md).','',
        '## Fixed-lift drag changes','',
        'Differences use the baseline with the same numerical settings. Ranges span coarse, base, fine and doubled wake. '
        'They are observed sensitivity ranges, not statistical intervals or error bounds.','',
        '| Geometry | Base change, counts | Range, counts | Sign stable? |','|---|---:|---:|---|']
    if ('baseline','ncrit6_fine') in records:
        lines[5:5]=['Three additional fine-grid Ncrit 6 checks were selected after the initial reversal appeared. '
                     'They are an adaptive follow-up, not independent holdouts.','']
    for r in summary['geometry_comparison']:
        if r['complete']:
            lines.append(f"| {r['geometry']} | {delta(r['geometry'],'base'):+.3f} | {r['min_delta_cd_counts']:+.3f} to {r['max_delta_cd_counts']:+.3f} | {r['sign_stable']} |")
        else:
            lines.append(f"| {r['geometry']} | unavailable | incomplete | unknown |")
    base_drag=[value('baseline',s) for s in ('coarse','base','fine')]
    grid_span=(max(base_drag)-min(base_drag))*1e4
    lines += ['',f'The baseline absolute drag spans {grid_span:.3f} counts across the three panel densities. '
        'This study does not establish sub-count absolute accuracy. Difference consistency can still '
        'support a local screening direction because numerical bias can cancel between nearby geometries. '
        'It does not prove that the physical model predicts the same improvement.', '',
        'The perturbation endpoints are sampled candidates, not optima. No combined camber/placement '
        'search, off-design sweep, aircraft trim constraint or independent finalist verification was performed.']
    lines += ['',f"{sum(p['order_stable'] for p in pairs)}/{len(pairs)} pairwise orders retain their sign across all four settings. "
        'This is internal consistency, not agreement with MSES or experiment on changed geometries.','',
        '## Baseline results','', '| Setting | Alpha, degrees | CL | CD | CM |','|---|---:|---:|---:|---:|']
    baseline_settings=(*SETTINGS,'ncrit6','start0')+(('ncrit6_fine',) if ('baseline','ncrit6_fine') in records else ())
    for s in baseline_settings:
        r=records[('baseline',s)]
        lines.append(f"| {s} | {value('baseline',s,'alpha_deg'):.6f} | {value('baseline',s,'CL'):.7f} | {value('baseline',s):.8f} | {value('baseline',s,'CM'):.6f} |" if r['feasible'] else
                     f'| {s} | failed | unavailable | unavailable | unavailable |')
    budget_path=out/'base_drag_budget.csv'
    if budget_path.exists():
        budgets={(r['geometry'],r['element']):float(r['CD_wake']) for r in read_csv(budget_path)}
        lines += ['', '## Element contributions at base settings','',
            '| Geometry | Main wake-drag change, counts | Flap wake-drag change, counts |',
            '|---|---:|---:|']
        for g in ('camber_minus','camber_plus','flap_down','flap_up'):
            if all((g,e) in budgets and ('baseline',e) in budgets for e in ('main','flap')):
                d=[1e4*(budgets[(g,e)]-budgets[('baseline',e)]) for e in ('main','flap')]
                lines.append(f'| {g} | {d[0]:+.3f} | {d[1]:+.3f} |')
        lines += ['', 'The two wake contributions sum to total profile drag. Skin friction is separately '
            'recorded in `base_drag_budget.csv` and must not be added to wake drag. Opposing component '
            'changes can make a small total gain sensitive to modeling and discretization.','']
    lines += ['', '## Directional differences and robustness','',
        'Central differences below are dCD/dp, where p is camber amplitude/local chord or vertical shift/reference chord. '
        'Both endpoint geometries are re-solved at the required lift.','',
        '| Setting | Mode | Step | dCD/dp |','|---|---|---:|---:|']
    for d in read_csv(out/'directional_differences.csv'):
        lines.append(f"| {d['setting']} | {d['mode']} | {d['step']} | {float(d['dCD_dp']):.7f} |" if d['dCD_dp'] else
                     f"| {d['setting']} | {d['mode']} | {d['step']} | unavailable |")
    full=(value('camber_plus','base')-value('camber_minus','base'))/.002
    half=(value('camber_plus_half','base')-value('camber_minus_half','base'))/.001
    step_change=100*abs(full-half)/max(abs(full),abs(half),1e-15)
    path_counts=(value('baseline','start0')-value('baseline','base'))*1e4
    lines += ['',f'The full/half camber-step derivatives differ by {step_change:.3f}% '
        '(normalized by the larger magnitude). This is one local step-size check, not derivative convergence.', '',
        f'The alternate baseline initialization changes CD by {path_counts:+.6f} counts. '
        'Only one alternate path was checked.', '',
        f"For positive camber, the change is {delta('camber_plus','base'):+.3f} counts at Ncrit 9 and "
        f"{delta('camber_plus','ncrit6'):+.3f} counts at Ncrit 6. "
        'Ncrit sensitivity does not establish roughness or forced-transition robustness.','']
    if delta('camber_plus','base')*delta('camber_plus','ncrit6')<0:
        lines += ['**Transition sensitivity reverses the preferred camber direction.** '
            'The nominal Ncrit 9 gain is not robust across the two tested transition settings. '
            'This is a concrete reason to use an objective covering multiple operating and transition '
            'conditions and to verify changed geometries independently before accepting an optimum.','']
    if ('baseline','ncrit6_fine') in records:
        fine_delta=delta('camber_plus','ncrit6_fine')
        fine_derivative=(value('camber_plus','ncrit6_fine')-value('camber_minus','ncrit6_fine'))/.002
        outcome=(f'{fine_delta:+.3f} counts for positive camber' if np.isfinite(fine_delta) else
                 'no baseline-relative drag change because its baseline did not converge')
        lines += [f'The fine-grid Ncrit 6 follow-up gives {outcome} '
            f'and a central derivative dCD/dp = {fine_derivative:.7f}. '
            'Its separate study, raw solutions and comparison table are in `../transition_followup/`.','']
        if np.isfinite(fine_delta) and fine_delta*delta('camber_plus','ncrit6')<0:
            lines += ['The baseline-relative positive-camber change at Ncrit 6 changes sign between '
                'the base and fine calculations. The fine baseline required a different bracket and '
                'initialization, so this is combined resolution/procedure sensitivity, not an isolated '
                'grid-error estimate. The central difference still uses the two modified geometries '
                'with the same original procedure and remains positive on both grids. A small '
                'positive-camber benefit at Ncrit 6 is not resolved by these checks.','']
        attempts=json.loads((out/'followup_attempts.json').read_text())
        for attempt in attempts:
            if not attempt['feasible']:
                lines += [f"Recorded failed {attempt['attempt']} attempt: {attempt['geometry']} / {attempt['setting']}: "
                          f"{attempt['failure_reason']}. This failure is retained, not omitted from the evidence.",'']
        r=records[('baseline','ncrit6_fine')]
        if 'retry_start_alpha' in r:
            lines += [f"The fine-grid baseline was retried from alpha {r['retry_start_alpha']:g} degrees "
                f"using bracket {r['retry_alpha_bounds']} inside the original [4,7] domain, "
                'with unchanged physical conditions. The reported baseline uses the final retry; the modified '
                'geometries retain their original initialization. See `followup_attempts.json` and '
                f"`../{r['retry_directory']}/` for the attempts. The earlier failures occurred at the "
                'lower bracketing incidence before target lift was reached; they do not establish '
                'physical infeasibility at the requested CL. Branch uniqueness is not proven.','']
    lines += ['', '## Official MSES context','',
        '[MIT example](https://web.mit.edu/drela/Public/web/mses/nlf26.pdf): MSES 3.11, alpha 6 degrees, '
        'CL 2.4840 and CD 0.02264, with the same printed Re, Mach and Ncrit 9. Exact coordinates and full input deck remain unconfirmed.','',
        '| Matching choice | MultiFoil alpha | MultiFoil CL | MultiFoil CD | Difference, counts |','|---|---:|---:|---:|---:|']
    for r in context:
        lines.append(f"| {r['comparison']} | {r['alpha_deg']:.6f} | {r['model_CL']:.7f} | {r['model_CD']:.8f} | {r['delta_cd_counts']:+.3f} |")
    lines += ['', 'These remain unscored arithmetic comparisons. A single baseline MSES point cannot validate geometry rankings. '
        'The matching choice changes the drag comparison because the codes predict different lift at the same incidence.','',
        '## Audit files','',
        '`study.json` records settings and code/input hashes. `results.csv` records all final outcomes; '
        '`directional_differences.csv` contains central finite differences of re-solved fixed-lift drag; '
        '`pairwise_rankings.csv` compares geometry order; `official_mses_context.csv` records the qualified comparisons. '
        'Per-job folders retain modified coordinates, requested flow evaluations, actual panel/wake coordinates, '
        'boundary-layer and transition diagnostics, versions and failure reasons.','',
        'Small residuals and stable signs do not validate missing wake confluence, cove recirculation, massive separation or shocks. '
        'No independent MSES/RANS runs, new wind-tunnel data or continuous shape optimization were performed.','']
    (out/'REPORT.md').write_text('\n'.join(lines))


if __name__=='__main__':
    main()
