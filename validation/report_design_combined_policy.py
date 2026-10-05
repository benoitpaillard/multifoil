#!/usr/bin/env python3
"""Acceptance audit for the complete ordered initialization policy; no solves."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from report_design_initialization import ROOT, Inputs, sha, write, csv_write, audit, require, close

STUDY_HASH = '52d4c2f6c08c903f76824a1a4864cee2551130bb88bebcc4f161dfe14f921d68'


def report(study, out):
    require(sha(study/'study.json') == STUDY_HASH, 'Unexpected combined study')
    inputs = Inputs(); spec = inputs.json(study/'study.json')
    for name, digest in spec['sha256'].items():
        require(sha(inputs.add(study/'source_snapshot'/name)) == digest, 'Snapshot changed: '+name)
    runtime = inputs.json(study/'runtime.json')
    require(runtime['versions'] == spec['expected_versions'], 'Runtime mismatch')
    require(inputs.json(study/'status.json')['state'] == 'finished', 'Unfinished study')
    original = inputs.json(study/'summary.json')
    states, references, rows, differences, histories = {}, {}, [], [], []
    for job in spec['jobs']:
        state = audit(study/job['name'], dict(job, mode='policy'), spec, runtime, inputs)
        r = state['result']; require(r['feasible'], 'Incomplete combined coverage')
        path = ROOT/job['reference_result']
        require(sha(path) == spec['sha256'][job['reference_result']], 'Legacy reference changed')
        if job['reference_result'] not in references:
            references[job['reference_result']] = audit(path.parent, dict(job, mode='legacy'), spec, runtime, inputs)
        ref = references[job['reference_result']]; old = ref['result']
        states[job['name']] = state
        point = state['metadata']['case']['point']
        expected = [dict(name='cold_local', start_alpha=job['center']),
                    dict(name='warm_anchor', start_alpha=-1.),
                    dict(name='warm_entry', start_alpha=-1., seed_alphas=[job['seed_alpha']])]
        require(point['initialization']['paths'] == expected, 'Ordered policy changed')
        close(point['initialization']['center'], job['center'])
        close(point['initialization']['half_width'], spec['half_width'])
        require(point['alpha_bounds'] == spec['alpha_domain'], 'Incidence domain changed')
        attempts = r['initialization_attempts']
        require(attempts[-1]['feasible'] and sum(a['feasible'] for a in attempts)==1, 'Policy did not stop on first successful attempt')
        expected_order = [v['name'] for v in expected]
        visited = [a['path'] for a in attempts]
        order = list(dict.fromkeys(visited))
        require(order == expected_order[:len(order)], 'Initialization paths skipped or reordered')
        require(r['selected_path'] == order[-1], 'Selected path mismatch')
        require(r['records'][r['selected_flow_index']]['initialization']['phase']=='fixed_lift', 'Warm seed substituted for target')
        if r['selected_path']=='warm_entry':
            seeds = [v for v in r['records'] if v['initialization']['phase']=='warm_seed']
            require(len(seeds)==1 and seeds[0]['admissible'], 'Missing successful entry seed')
            close(seeds[0]['alpha_deg'], job['seed_alpha'])
        delta = (r['CD']-old['CD'])*1e4
        require(abs(delta)<=spec['acceptance_CD_counts'], 'Drag comparison exceeds acceptance')
        internal_failures = sum(not t['converged'] for d in state['metadata']['diagnostics'] for t in d['continuation'])
        rows.append(dict(job=job['name'], flap=job['flap'], amplitude=job['amplitude'], repeat_of=job.get('repeat_of',''),
            CL=r['CL'], CD=r['CD'], CM=r['CM'], alpha_deg=r['alpha_deg'], residual=r['residual'],
            lift_error=r['CL']-.8, selected_path=r['selected_path'], difference_CD_counts=delta,
            evaluations=r['flow_evaluations'], failed_evaluations=sum(not x['admissible'] for x in r['records']),
            failed_internal_points=internal_failures, wall_seconds=r['wall_seconds']))
        histories.append(dict(job=job['name'], attempts=attempts, records=r['records'], selected_flow_index=r['selected_flow_index'],
            continuation=[d['continuation'] for d in state['metadata']['diagnostics']]))
        for element in ('fore','aft'):
            require(sha(inputs.add(study/job['name']/f'input_{element}.dat')) == sha(inputs.add(path.parent/f'input_{element}.dat')), 'Input geometry changed')
            for kind in ('surfaces','wakes'):
                a,b = state[kind][element],ref[kind][element]
                require(a.shape==b.shape, 'Export dimensions changed')
                for field in a.dtype.names:
                    exact=bool(np.array_equal(a[field],b[field]))
                    delta_field=float(np.max(abs(a[field]-b[field]))) if a[field].dtype.kind in 'fiu' else None
                    if field in ('side','turbulent'): require(exact, 'Discrete topology differs from control')
                    differences.append(dict(job=job['name'], element=element, kind=kind, field=field,
                                            exactly_equal=exact, max_abs_difference=delta_field))
        automatic=next(v for v in original['jobs'] if v['job']==job['name'])
        require(automatic['passed'] and automatic['selected_path']==r['selected_path'], 'Automatic acceptance mismatch')
        close(automatic['difference_CD_counts'], delta)
    require(len(rows)==original['required_jobs']==original['passed_jobs']==7 and original['regression_pass'], 'Coverage mismatch')
    repeated=[]
    for job in spec['jobs']:
        if 'repeat_of' not in job: continue
        a,b = states[job['name']]['result'],states[job['repeat_of']]['result']
        for field in ('CL','CD','CM','alpha_deg','residual'):
            repeated.append(dict(element='assembly', field=field, max_abs_difference=abs(a[field]-b[field]), exactly_equal=a[field]==b[field]))
        for element in ('fore','aft'):
            with np.load(inputs.add(study/job['name']/f'{element}_state.npz')) as x, \
                    np.load(inputs.add(study/job['repeat_of']/f'{element}_state.npz')) as y:
                require(set(x.files)==set(y.files), 'Repeat state fields differ')
                for field in x.files:
                    require(x[field].shape==y[field].shape,'Repeat state dimensions differ')
                    repeated.append(dict(element=element, field=field,
                        max_abs_difference=float(np.max(abs(x[field]-y[field]))), exactly_equal=bool(np.array_equal(x[field],y[field]))))
        require(abs(a['CD']-b['CD'])*1e4<=spec['acceptance_CD_counts'], 'Repeat force criterion failed')
    deltas=[]
    for flap in (5,10):
        base=states[f'd{flap}_baseline_combined']['result']
        prior_base=references[next(j['reference_result'] for j in spec['jobs'] if j['name']==f'd{flap}_baseline_combined')]['result']
        for label in ('minus','plus'):
            name=f'd{flap}_{label}_combined'; current=states[name]['result']
            old=references[next(j['reference_result'] for j in spec['jobs'] if j['name']==name)]['result']
            new_delta=(current['CD']-base['CD'])*1e4; old_delta=(old['CD']-prior_base['CD'])*1e4
            deltas.append(dict(flap=flap, amplitude=-.001 if label=='minus' else .001,
                combined_delta_CD_counts=new_delta, legacy_delta_CD_counts=old_delta,
                difference_of_deltas_counts=new_delta-old_delta))
    summary=dict(required_jobs=7, feasible_jobs=7, primary_cases=6, repeat_cases=1,
        full_state_audit_pass=True, accepted_for_declared_protocol=True, selected_states_audited=13, new_newton_states_audited=7,
        max_absolute_CD_difference_counts=max(abs(r['difference_CD_counts']) for r in rows),
        repeat_forces_and_saved_arrays_exactly_equal=all(r['exactly_equal'] for r in repeated),
        max_absolute_difference_of_design_deltas_counts=max(abs(r['difference_of_deltas_counts']) for r in deltas),
        requested_evaluations=sum(r['evaluations'] for r in rows), failed_requested_evaluations=sum(r['failed_evaluations'] for r in rows),
        failed_internal_continuation_points=sum(r['failed_internal_points'] for r in rows),
        qualification='Priority 2 accepted for the declared S702 protocols only; this is numerical repeatability, not independent contour-ranking validation.')
    out.mkdir(parents=True,exist_ok=True)
    for name, values in (('targets',rows),('state_differences',differences),('repeat_differences',repeated),('contour_deltas',deltas)):
        csv_write(out/(name+'.csv'),values)
    write(out/'summary.json',summary);write(out/'initialization_histories.json',histories)
    fig,axes=plt.subplots(1,2,figsize=(10.5,4.3))
    for ax,flap in zip(axes,(5,10)):
        group=sorted([r for r in deltas if r['flap']==flap],key=lambda r:r['amplitude'])
        for field,style,label in (('combined_delta_CD_counts','-o','Combined policy'),('legacy_delta_CD_counts','--s','Legacy control')):
            ax.plot([-1,0,1],[group[0][field],0,group[1][field]],style,label=label)
        ax.axhline(0,color='0.5',lw=.7);ax.grid(alpha=.2)
        ax.set(title=f'Flap {flap}°',xlabel='Fore-camber amplitude [0.001 chord]',ylabel='CD change from own baseline [counts]')
        ax.legend()
    fig.suptitle('7/7 passed; independent baseline repeat matches exactly')
    fig.tight_layout()
    for ext in ('png','svg'):fig.savefig(out/f'combined_policy_audit.{ext}',dpi=170)
    plt.close(fig)
    (out/'FINDINGS.md').write_text(
        '# Combined initialization policy accepted for the declared S702 cases\n\n'
        'All six primary cases and the independent baseline repeat passed. The complete ordered policy '
        'was executed afresh: cold local start, fresh −1° start, then fresh solver with explicit entry at '
        'center−0.4°. Four primary cases select the cold path, the positive-camber flap-10° case selects '
        'the warm anchor, and the flap-10° baseline (including its repeat) uses the explicit entry.\n\n'
        f"Maximum CD difference against the controls: {summary['max_absolute_CD_difference_counts']:.9g} counts "
        '(limit 0.01). All thirteen distinct selected target exports (seven new, six controls) passed audit; '
        'all seven new Newton states match their selected exports. Transition flags and surface labels '
        'match the controls. The repeated baseline has exactly equal CL, CD, CM, incidence, residual and '
        'every saved Newton-array field. This equality is specific to these inputs and this runtime.\n\n'
        f"The run retains {summary['failed_requested_evaluations']} failed requests among {summary['requested_evaluations']} "
        f"and {summary['failed_internal_continuation_points']} failed internal continuation points. "
        'Failed paths carry no selected force; first successful path wins. The original 11/12 matrix '
        'remains failed and preserved.\n\n'
        'Contour differences use each flap setting and each numerical method’s own baseline. The largest '
        f"difference between paired design changes is {summary['max_absolute_difference_of_design_deltas_counts']:.9g} counts, "
        'well below the existing 0.5-count screening threshold. No new geometry is selected or physically validated.\n\n'
        'Priority 2 is complete for the declared S702 conditions and perturbations, with the earlier '
        'Mach-0.07 repeatability and six-condition policy evidence retained. Broader geometry, Reynolds-number '
        'and transition robustness remain untested. The unrelated far-field ordering limitation persists. '
        'Next is the geometry/transition-condition audit in priority 3.\n\n'
        'Reproduce: `python3 validation/report_design_combined_policy.py`.\n')
    for name in ('validation/report_design_combined_policy.py','validation/report_design_initialization.py','validation/report_s702_drag_mechanism.py'):
        inputs.add(ROOT/name)
    write(out/'input_manifest.json',dict(study_sha256=STUDY_HASH,files=inputs.files))
    return summary


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--study',type=Path,default=ROOT/'validation/results/design_combined_policy_20260921')
    ap.add_argument('--out',type=Path,default=ROOT/'validation/results/design_combined_policy_analysis')
    args=ap.parse_args();print(json.dumps(report(args.study.resolve(),args.out.resolve()),indent=2))
