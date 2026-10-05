#!/usr/bin/env python3
"""Audit saved S702 fore/aft sampling states and compare with experiment/MSES."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from report_s702_discretization import (
    ROOT, Inputs, sha, write, csv_write, audit, audit_setup, compare_arrays,
    require, close, array_close, clearance, side_data, resolution,
    contour_sample_distance, DEPENDENCIES, ELEMENTS,
)
from run_s702_surface_factorial import CELLS, effects

STUDY_HASH = '125e68020d67746ea67b70075f68d67a215e07576ef31ea73c1f2f494fbfaf0e'
POINTS = 'validation/results/s702_drag_mechanism/comparison_points.csv'
POINTS_HASH = 'aaf6a220d63d0d85202d90182e54f79e10fa37b78d67725672b592c3f2b475ac'
PARENT_AUDIT = 'validation/results/s702_discretization_analysis'
TERMS = ('fore_at_low_aft','fore_at_high_aft','aft_at_low_fore','aft_at_high_fore','interaction','total')
PAIRS = dict(fore_at_low_aft=('B','F'),fore_at_high_aft=('A','C'),
             aft_at_low_fore=('B','A'),aft_at_high_fore=('F','C'),
             interaction=('B','F','A','C'),total=('B','C'))


def check_fixed_controls(states, setups):
    require(set(states) == set(setups) == set(CELLS), 'Complete B/F/A/C states and setups required')
    for cell in CELLS:
        base, setup = setups['B'], setups[cell]
        for key in ('source_sha256','recipes'):
            require(setup[key] == base[key], 'Fixed source/recipe changed: '+key)
        for key in ('wake_nodes','wake_first_spacing'):
            require(setup['requested'][key] == base['requested'][key], 'Fixed wake control changed: '+key)
        require(setup['requested']['surface_nodes'] == dict(zip(ELEMENTS,[n+1 for n in CELLS[cell]])),
                'Factorial surface count differs: '+cell)
        for element in ELEMENTS:
            array_close(states['B']['arrays'][element]['source_xy'], states[cell]['arrays'][element]['source_xy'],
                        'Fixed regularized source changed', 0.)
    # Holding an element's count fixed must hold its realized surface exactly fixed.
    for element, pairs in (('fore',(('B','A'),('F','C'))), ('aft',(('B','F'),('A','C')))):
        for first, second in pairs:
            for key in ('xy','foil_s','foil_t'):
                array_close(states[first]['arrays'][element][key], states[second]['arrays'][element][key],
                            'Fixed element surface changed: '+element+'/'+key, 0.)


def contrasts(states):
    require(set(states) == {(f,c) for f in (5,10) for c in CELLS}, 'Complete eight-cell coverage required')
    for state in states.values():
        r=state['result']
        require(r['feasible'] and r['CD'] is not None and np.isfinite(r['CD']), 'Complete finite coverage required')
    rows=[]
    for component in ('total','fore','aft'):
        values={}
        for flap in (5,10):
            values[flap]={c:states[flap,c]['result']['CD']*1e4 if component=='total'
                          else states[flap,c]['result']['element_forces'][component]['CD_wake']*1e4 for c in CELLS}
        values['paired']={c:values[10][c]-values[5][c] for c in CELLS}
        for scope in (5,10,'paired'):
            e=effects(values[scope])
            for term in TERMS:
                flaps=(5,10) if scope=='paired' else (scope,)
                matched=all(len({states[f,c]['result']['selected_path'] for c in PAIRS[term]})==1 for f in flaps)
                rows.append(dict(scope=str(scope),component=component,contrast=term,change_counts=e[term],
                    same_selected_path_within_each_flap=matched))
    lookup={(r['scope'],r['component'],r['contrast']):r for r in rows}
    for term in TERMS:
        for scope in ('5','10','paired'):
            close(lookup[scope,'total',term]['change_counts'],sum(lookup[scope,e,term]['change_counts'] for e in ELEMENTS))
        close(lookup['paired','total',term]['change_counts'],
              lookup['10','total',term]['change_counts']-lookup['5','total',term]['change_counts'])
    return rows


def reference_comparison(inputs, spec, prior, states):
    """Reuse hash-locked audited context; do not silently promote its numerical protocol."""
    require(sha(inputs.add(ROOT/POINTS)) == POINTS_HASH, 'Archived comparison table changed')
    points=inputs.csv(ROOT/POINTS)
    require(sha(inputs.add(ROOT/spec['reference_comparison']))==spec['sha256'][spec['reference_comparison']],
            'Frozen experimental/MSES context changed')
    paired=inputs.csv(ROOT/spec['reference_comparison'])
    original_manifest=inputs.json(ROOT/'validation/results/s702_drag_mechanism/input_manifest.json')['files']
    third_manifest=inputs.json(ROOT/'validation/results/s702_panel_convergence_analysis/input_manifest.json')['files']
    require(sha(ROOT/'validation/results/s702_panel_convergence_analysis/input_manifest.json') ==
            spec['sha256']['validation/results/s702_panel_convergence_analysis/input_manifest.json'], 'Third-grid audit changed')
    def old_result(name,pinned):
        r=inputs.json(ROOT/name);digest=pinned[name]
        require(sha(ROOT/name)==(digest['sha256'] if isinstance(digest,dict) else digest), 'Archived context state changed: '+name)
        require(r['feasible'] and abs(r['CL']-cl)<=spec['cl_tolerance'], 'Archived context state inadmissible')
        return r
    rows=[]
    for cl in (.8,1.2):
        group={}
        for series in ('experiment','mses_n9','multifoil_n9'):
            p={f:next(v for v in points if float(v['reynolds'])==spec['reynolds'] and float(v['mach'])==spec['mach']
                      and float(v['target_cl'])==cl and float(v['flap_deg'])==f and v['series']==series) for f in (5,10)}
            if series=='multifoil_n9':
                for f in (5,10):close(old_result(p[f]['result_path'],original_manifest)['CD'],float(p[f]['CD']))
            group[series]=dict(series=series,ncrit=None if series=='experiment' else 9.,
                protocol={'experiment':'Digitized experimental curve', 'mses_n9':'Published MSES curve; exact deck unavailable',
                          'multifoil_n9':'Archived 260/180 surfaces; automatic wake'}[series],
                d5_CD=float(p[5]['CD']),d10_CD=float(p[10]['CD']),
                source='NASA/CR-20220014964 Figures '+p[5]['reference_figure']+' and '+p[10]['reference_figure'],
                source_url=p[5]['reference_url'])
        for series in ('experiment','mses_n9','multifoil_n9'):
            expected=next(v for v in paired if float(v['reynolds'])==spec['reynolds'] and float(v['mach'])==spec['mach'] and float(v['target_cl'])==cl)
            close((group[series]['d10_CD']-group[series]['d5_CD'])*1e4,float(expected[series+'_delta_10_minus_5_counts']))
        if cl==.8:
            group['multifoil_C']=dict(series='multifoil_C',ncrit=6.25,protocol='Factorial C: 440/320 surfaces, fine wake recipes',
                d5_CD=states[5,'C']['result']['CD'],d10_CD=states[10,'C']['result']['CD'],source='This audited factorial',source_url='')
            archived={f:old_result(f'validation/results/s702_discretization_20261001/d{f}_wake_spacing/result.json',prior) for f in (5,10)}
        else:
            archived={f:old_result(f'validation/results/s702_panel_convergence_20260922/d{f}_cl1p2_finer/result.json',third_manifest) for f in (5,10)}
        group['multifoil_finer']=dict(series='multifoil_finer',ncrit=6.25,protocol='Archived 440/320 surfaces and finer wake recipes',
            d5_CD=archived[5]['CD'],d10_CD=archived[10]['CD'],source='Previous audited finest operational grid',source_url='')
        for row in group.values():
            row=dict(target_cl=cl,reynolds=spec['reynolds'],mach=spec['mach'],**row)
            row['flap_delta_counts']=(row['d10_CD']-row['d5_CD'])*1e4
            for f in (5,10):row[f'd{f}_error_vs_experiment_counts']=(row[f'd{f}_CD']-group['experiment'][f'd{f}_CD'])*1e4
            row['delta_error_vs_experiment_counts']=row['d10_error_vs_experiment_counts']-row['d5_error_vs_experiment_counts']
            rows.append(row)
    return rows


def report(study,out):
    require(sha(study/'study.json')==STUDY_HASH,'Unexpected factorial study')
    inputs=Inputs();spec=inputs.json(study/'study.json')
    for name,digest in spec['sha256'].items():
        require(sha(inputs.add(study/'source_snapshot'/name))==digest,'Frozen snapshot changed: '+name)
    prior=inputs.json(study/'source_snapshot'/PARENT_AUDIT/'input_manifest.json')['files']
    for name in (*DEPENDENCIES,'validation/report_s702_discretization.py','validation/run_s702_surface_factorial.py'):
        digest=spec['sha256'].get(name,prior.get(name))
        require(digest is not None and sha(inputs.add(ROOT/name))==digest,'Audit dependency changed: '+name)
    for name in ('fore.dat','aft_d5.dat','aft_d10.dat'):
        key='validation/data/s702/'+name
        require(sha(inputs.add(ROOT/key))==spec['sha256'][key],'Raw geometry changed')
    runtime=inputs.json(study/'runtime.json');status=inputs.json(study/'status.json');automatic=inputs.json(study/'summary.json')
    require(runtime['versions']==spec['expected_versions'] and all(v=='1' for v in runtime['thread_limits'].values()),'Runtime differs')
    require(status['state']=='finished' and status['completed_jobs']==status['total_jobs']==5,'Incomplete execution')
    require(automatic['required_jobs']==automatic['feasible_jobs']==5 and automatic['execution_checks_pass'],'Runner gate failed')
    states,setups={},{}
    targets,histories,surfaces,wakes,geometry=[],[],[],[],[]
    def selected(folder,job,study_spec,run,cell,new,repeat=False):
        state=audit(folder,job,dict(study_spec,panels=job['panels'],source_panels=job['source_panels']),run,inputs)
        r=state['result'];require(r['feasible'],'Failed target: complete attribution prohibited')
        setup=audit_setup(folder,job,study_spec,state,inputs)
        targets.append(dict(job=job['name'],flap=job['flap'],cell=cell,new=new,repeat=repeat,
            panels='/'.join(map(str,job['panels'])),CL=r['CL'],CD=r['CD'],CM=r['CM'],alpha_deg=r['alpha_deg'],
            lift_error=r['CL']-job['target_cl'],residual=state['residual'],selected_path=r['selected_path'],
            evaluations=len(r['records']),failed_requests=sum(not rec['admissible'] for rec in r['records']),
            failed_internal_points=sum(not c['converged'] for d in state['metadata']['diagnostics'] for c in d['continuation']),
            wall_seconds=r['wall_seconds']))
        histories.append(dict(job=job['name'],new=new,repeat=repeat,records=r['records'],attempts=r['attempts'],
            selected_flow_index=r['selected_flow_index'],continuation=[d['continuation'] for d in state['metadata']['diagnostics']]))
        diag=state['metadata']['diagnostics'][r['selected_flow_index']]
        for e in ELEMENTS:
            a=state['arrays'][e];other='aft' if e=='fore' else 'fore'
            gap=clearance(a['wake_xy'].T,state['arrays'][other]['xy'].T);require(gap>0,'Wake crosses another element')
            ds=np.linalg.norm(np.diff(a['wake_xy'],axis=1),axis=0)
            wakes.append(dict(job=job['name'],cell=cell,flap=job['flap'],element=e,surface_nodes=a['xy'].shape[1],
                source_nodes=a['source_xy'].shape[1],wake_nodes=a['wake_xy'].shape[1],
                requested_first_spacing=setup['requested']['wake_first_spacing'][e],realized_first_spacing=float(ds[0]),
                arclength=float(sum(ds)),last_spacing=float(ds[-1]),clearance=gap,
                CD_wake=r['element_forces'][e]['CD_wake'],CD_friction=r['element_forces'][e]['CD_friction'],
                TE_gap_added=diag['te_gap_added_over_local_chord'][e],paneling_recovery=diag['monotone_paneling_recovery'][e]))
            for side in ('upper','lower'):
                bl=r['boundary_layers'][e][side]
                surfaces.append(dict(job=job['name'],cell=cell,flap=job['flap'],element=e,side=side,
                    transition_x=bl['transition']['x_over_local_chord'] if bl['transition'] else None,
                    separated_fraction=bl['separated_arclength_fraction'],min_cf=bl['min_cf'],max_Hk=bl['max_Hk'],
                    separation_regions=json.dumps(bl['separation_regions'],sort_keys=True),
                    **resolution(side_data(state['surfaces'][e],side))))
        return state,setup
    for name in ('study.json','runtime.json'):
        key=spec['parent_study']+'/'+name
        require(sha(inputs.add(ROOT/key))==spec['sha256'][key],'Frozen parent manifest/runtime changed')
    parent=inputs.json(ROOT/spec['parent_study']/'study.json');parent_runtime=inputs.json(ROOT/spec['parent_study']/'runtime.json')
    for key in ('reynolds','mach','ncrit','reference_chord','iterations','cl_tolerance','residual_tolerance',
                'initial_half_width','alpha_bounds','max_evaluations','max_brackets','warm_start_alpha','entry_offset'):
        require(parent[key]==spec[key],'Parent physics/protocol differs: '+key)
    require(parent_runtime['versions']==runtime['versions'] and parent_runtime['python']==runtime['python'],'Parent runtime differs')
    for ref in spec['references']:
        before=set(inputs.files);key=(ref['flap'],ref['cell'])
        require(ref['job'] in parent['jobs'],'Reference job absent from parent')
        states[key],setups[key]=selected((ROOT/ref['result']).parent,ref['job'],parent,parent_runtime,ref['cell'],False)
        for name in set(inputs.files)-before:
            require(name in spec['sha256'] and inputs.files[name]==spec['sha256'][name],'Unpinned reference: '+name)
    control=None
    for job in spec['jobs']:
        state,setup=selected(study/job['name'],job,spec,runtime,job['cell'],True,job['repeat_control'])
        if job['repeat_control']:control=state
        else:states[job['flap'],job['cell']],setups[job['flap'],job['cell']]=state,setup
        auto=next(v for v in automatic['jobs'] if v['job']==job['name'])
        for key in ('CD','CL'):close(state['result'][key],auto[key])
    require(len(targets)==9 and control is not None,'Nine selected states required')
    original=states[10,'C'];control_delta=(control['result']['CD']-original['result']['CD'])*1e4
    require(abs(control_delta)<=spec['acceptance_control_CD_counts'],'Repeat drag differs')
    require(all(set(control['arrays'][e])==set(original['arrays'][e]) for e in ELEMENTS),'Control fields differ')
    controls=compare_arrays(control['arrays'],original['arrays'],'C repeat/archive')
    require(len(controls)==24 and all(v['exactly_equal'] for v in controls),'Control saved arrays differ')
    for flap in (5,10):
        check_fixed_controls({c:states[flap,c] for c in CELLS},{c:setups[flap,c] for c in CELLS})
        for first,second in (('B','F'),('B','A'),('F','C'),('A','C')):
            for e in ELEMENTS:
                a,b=states[flap,first]['arrays'][e],states[flap,second]['arrays'][e]
                geometry.append(dict(flap=flap,first=first,second=second,element=e,
                    sampled_contour_distance=contour_sample_distance(a['xy'].T,b['xy'].T),
                    max_TE_endpoint_change=float(np.max(abs(a['xy'][:,[0,-1]]-b['xy'][:,[0,-1]])))))
    rows=contrasts(states)
    for row in rows:
        group=automatic['paired_effects_counts'] if row['scope']=='paired' else automatic['individual_effects_counts'][row['scope']]
        component='total_drag' if row['component']=='total' else row['component']+'_drag'
        close(row['change_counts'],group[component][row['contrast']])
    pairs=[dict(cell=c,panels='/'.join(map(str,CELLS[c])),flap_delta_counts=(states[10,c]['result']['CD']-states[5,c]['result']['CD'])*1e4,
                d5_path=states[5,c]['result']['selected_path'],d10_path=states[10,c]['result']['selected_path']) for c in CELLS]
    for p in pairs:close(p['flap_delta_counts'],next(v['flap_delta_CD_counts'] for v in automatic['paired_cells'] if v['cell']==p['cell']))
    comparison=reference_comparison(inputs,spec,prior,states)
    new=[r for r in targets if r['new']]
    summary=dict(required_jobs=5,feasible_jobs=5,selected_states_audited=9,full_state_audit_pass=True,
        factorial_controls_pass=True,control_CD_difference_counts=control_delta,control_saved_arrays_exactly_equal=True,
        paired_cells=pairs,paired_contrasts=[r for r in rows if r['scope']=='paired' and r['component']=='total'],
        same_path_attribution_terms=[r['contrast'] for r in rows if r['scope']=='paired' and r['component']=='total' and r['same_selected_path_within_each_flap']],
        requested_evaluations=sum(v['evaluations'] for v in new),failed_requested_evaluations=sum(v['failed_requests'] for v in new),
        failed_internal_continuation_points=sum(v['failed_internal_points'] for v in new),
        all_nine_requested_evaluations=sum(v['evaluations'] for v in targets),
        all_nine_failed_requests=sum(v['failed_requests'] for v in targets),
        all_nine_failed_internal_points=sum(v['failed_internal_points'] for v in targets),
        new_sum_job_wall_seconds=sum(v['wall_seconds'] for v in new),
        minimum_wake_body_clearance=min(v['clearance'] for v in wakes),
        maximum_absolute_lift_error=max(abs(v['lift_error']) for v in targets),maximum_selected_residual=max(v['residual'] for v in targets),
        artificial_TE_gap_added=any(v['TE_gap_added']!=0 for v in wakes),paneling_recovery_used=any(v['paneling_recovery'] for v in wakes),
        numerical_convergence_established=False,independent_validation_pass=False,comparison= comparison,
        qualification='Aft surface sampling dominates the measured numerical response; one mixed cell changes selected path. No error bound or physical validation.')
    out.mkdir(parents=True,exist_ok=True)
    for name,table in (('targets',targets),('contrasts',rows),('paired_cells',pairs),('comparison',comparison),
                       ('surface_diagnostics',surfaces),('wake_diagnostics',wakes),('geometry_comparison',geometry),('control_differences',controls)):
        csv_write(out/(name+'.csv'),table)
    write(out/'initialization_histories.json',histories);write(out/'summary.json',summary)
    plot(out,states,pairs,rows,comparison);findings(out,summary)
    inputs.add(Path(__file__));write(out/'input_manifest.json',dict(study_sha256=STUDY_HASH,files=inputs.files))
    return summary


def save(fig,out,name):
    for ext in ('png','svg','pdf'):fig.savefig(out/(name+'.'+ext),dpi=180)
    plt.close(fig)


def plot(out,states,pairs,rows,comparison):
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(1,2,figsize=(12,4.6),layout='constrained')
    values={r['cell']:r['flap_delta_counts'] for r in pairs}
    matrix=np.array([[values['B'],values['F']],[values['A'],values['C']]])
    axes[0].imshow(matrix,cmap='Blues_r',vmin=-2.5,vmax=-1.7)
    for (iy,ix),cell in zip(((0,0),(0,1),(1,0),(1,1)),('B','F','A','C')):
        axes[0].text(ix,iy,f'{cell}: {values[cell]:+.4f}\n'+('warm entry*' if cell=='A' else 'warm anchor'),ha='center',va='center',
                     color='white' if cell in ('A','C') else '#152c42')
    axes[0].set_xticks([0,1],['340','440']);axes[0].set_yticks([0,1],['240','320'])
    axes[0].set_xlabel('Fore surface panels');axes[0].set_ylabel('Aft surface panels')
    axes[0].set_title('Flap-10° minus flap-5° drag [counts]')
    terms=('fore_at_low_aft','aft_at_high_fore','interaction')
    selected=[next(r for r in rows if r['scope']=='paired' and r['component']=='total' and r['contrast']==t) for t in terms]
    bars=axes[1].bar(np.arange(3),[r['change_counts'] for r in selected],color=['#357cac','#c47729','#888888'])
    axes[1].bar_label(bars,labels=[f"{r['change_counts']:+.4f}" for r in selected],padding=5)
    axes[1].axhline(0,color='black',lw=.8);axes[1].set_ylim(-.73,.18);axes[1].grid(axis='y',alpha=.2)
    axes[1].set_xticks([0,1,2],['Fore: B → F\nsame paths','Aft: F → C\nsame paths','Interaction*\nnot an extra increment'])
    axes[1].set_ylabel('Change in paired drag [counts]');axes[1].set_title('Conditional effects and interaction')
    fig.suptitle('S702 · CL = 0.8 · Re = 700,000 · M = 0.07 · Ncrit = 6.25\n'
                 'Source and wake recipes fixed. *Flap-10 cell A needs warm entry; all flap-5 cells use cold local.',fontsize=10)
    save(fig,out,'surface_factors')
    fig,axes=plt.subplots(1,2,figsize=(12,5.4),layout='constrained')
    styles=(('experiment','Experiment','#222222','o','-'),
            ('mses_n9','Published MSES · N=9','#8854a5','s','--'),
            ('multifoil_n9','MultiFoil · N=9 · archived 260/180','#cf9b28','^',':'),
            ('multifoil_C','MultiFoil C · N=6.25 · fixed fine wake','#247ca4','D','-'),
            ('multifoil_finer','MultiFoil · N=6.25 · finer wake (prior)','#4d9a62','v','--'))
    for i,(series,label,color,marker,line) in enumerate(styles):
        r=next(v for v in comparison if v['target_cl']==.8 and v['series']==series)
        axes[0].plot([5,10],[r['d5_CD']*1e4,r['d10_CD']*1e4],color=color,marker=marker,ls=line,label=label)
        axes[1].barh(i,r['flap_delta_counts'],color=color)
        value=r['flap_delta_counts'];axes[1].text(value+(.35 if value>=0 else -.35),i,f'{value:+.2f}',ha='left' if value>=0 else 'right',va='center')
    axes[0].set_xticks([5,10]);axes[0].set_xlabel('Flap angle [degrees]');axes[0].set_ylabel('CD [counts = 0.0001]')
    axes[0].set_title('Absolute drag at CL = 0.8');axes[0].grid(alpha=.2)
    axes[1].set_yticks(range(5),['Experiment','Published MSES N=9','MultiFoil N=9 (archive)','MultiFoil C N=6.25','MultiFoil finer N=6.25'])
    axes[1].invert_yaxis();axes[1].axvline(0,color='black',lw=.8);axes[1].set_xlim(-25,20)
    axes[1].set_xlabel('CD(10°) − CD(5°) [counts]');axes[1].set_title('Flap effect: negative means less drag')
    fig.legend(*axes[0].get_legend_handles_labels(),loc='outside lower center',ncol=2,fontsize=8.5)
    fig.suptitle('S702 · Re = 700,000 · M = 0.07 · common target CL = 0.8\n'
                 'Transition settings and numerical protocols differ as labeled; published MSES is not a new run.',fontsize=11)
    save(fig,out,'experiment_mses_comparison')
    fig,axes=plt.subplots(1,2,figsize=(11,4),layout='constrained')
    for flap,color in ((5,'#357cac'),(10,'#c47729')):
        for e in ELEMENTS:
            xy=states[flap,'C']['arrays'][e]['xy']
            axes[0].plot(xy[0],xy[1],color=color,label=f'Flap {flap}°' if e=='fore' else None)
    axes[0].axis('equal');axes[0].set_xlabel('x / cref');axes[0].set_ylabel('y / cref');axes[0].legend()
    axes[0].set_title('Actual C solver contours')
    for cell,color in (('B','#777777'),('F','#357cac'),('A','#ab76bb'),('C','#c47729')):
        data=side_data(states[10,cell]['surfaces']['aft'],'upper')
        axes[1].plot(data['x_over_local_chord'],data['Hk'],label=cell+' '+str(CELLS[cell]),color=color,ls='--' if cell in ('A','C') else '-')
    axes[1].set_xlim(.70,.95);axes[1].set_xlabel('x / aft local chord');axes[1].set_ylabel('Kinematic shape factor Hk')
    axes[1].set_title('Flap 10°: aft upper boundary layer');axes[1].legend(fontsize=8)
    for ax in axes:ax.grid(alpha=.2)
    save(fig,out,'geometry_and_aft_layer')


def findings(out,s):
    contrasts=s['paired_contrasts'];lookup={r['contrast']:r for r in contrasts}
    table='\n'.join(f"| {r['contrast']} | {r['change_counts']:+.6f} | {'yes' if r['same_selected_path_within_each_flap'] else 'no'} |" for r in contrasts)
    comparisons='\n'.join(f"| {r['target_cl']:g} | {r['series']} | {r['d5_CD']:.8f} | {r['d10_CD']:.8f} | {r['flap_delta_counts']:+.5f} | {r['delta_error_vs_experiment_counts']:+.5f} |" for r in s['comparison'])
    (out/'FINDINGS.md').write_text(
        '# S702 fore/aft surface sampling: completed audit\n\n'
        'All five new targets and all nine selected states (five new, four archived endpoints) passed full audit. '
        'The repeated C control reproduces drag exactly and all 24 saved array comparisons are exactly equal. '
        'Frozen hashes, runtime, every requested/applied control, reconstructed sources/surfaces, Newton arrays, '
        'nonzero-Mach export conversion, residual/lift tolerances, force sums, transition/separation, '
        'initialization order and wake clearance were checked. No flow solve or production-core change was made by this report.\n\n'
        '## Numerical result\n\n'
        'With 440/320 source contours and fine wake recipes fixed, B=340/240, F=440/240, A=340/320 and C=440/320 '
        'surface panels. All cases use Re 700,000, M 0.07, CL 0.8 and Ncrit 6.25. The paired drag change means '
        'CD(flap 10)−CD(flap 5). One count is 1e-4 CD.\n\n'
        '| Contrast in paired response | Change [counts] | Same selected path within each flap? |\n|---|---:|---|\n'+table+'\n\n'
        'Aft surface sampling dominates. The F→C aft change is -0.591506 count, while B→F fore sampling changes '
        'the paired response by +0.030226. Those two increments telescope to -0.561280 and both compare cold-local '
        'flap-5 states and warm-anchor flap-10 states. This supports the diagnosis without using the mixed A state. '
        'The alternate path gives an aft effect of -0.542617 and a fore effect of -0.018663; the factorial interaction '
        'is -0.048890 count. The interaction is a difference of conditional effects, not an extra term to add to B→F→C.\n\n'
        'Cell A at flap 10 requires warm-entry after cold-local and warm-anchor fail. Therefore its conditional contrasts '
        'and the interaction retain a path-change qualification; universal path independence has not been shown. '
        'No matched-path rerun is needed to support the limited same-path F→C observation, and no extra solve is launched. '
        'Same path names do not themselves prove uniqueness of a nonlinear solution. The original 0.5-count screen '
        'remains failed and no numerical error bound or valid small design gain follows from this study.\n\n'
        'Regularized sources and wake recipes are exactly fixed across all four cells at each flap. Each element\'s '
        'surface coordinates, arclength and tangents are identical whenever that element\'s count is unchanged. '
        f"Minimum wake-to-other-body centerline clearance is {s['minimum_wake_body_clearance']:.6f}c; this is not proof of absent physical wake/BL interaction. "
        f"Maximum lift error is {s['maximum_absolute_lift_error']:.3g}; maximum selected residual is {s['maximum_selected_residual']:.3g}. "
        'No artificial TE opening or paneling recovery was used. Detailed geometry, surface, wake and cost tables are saved.\n\n'
        f"The five new jobs contain {s['requested_evaluations']} requested evaluations, {s['failed_requested_evaluations']} failed requests and "
        f"{s['failed_internal_continuation_points']} failed internal continuation points. Including the four archived endpoints, the audit retains "
        f"{s['all_nine_requested_evaluations']} requests, {s['all_nine_failed_requests']} failed requests and {s['all_nine_failed_internal_points']} failed internal points. "
        'All failures remain recorded, with no substituted force values.\n\n'
        '## Experiment and published MSES\n\n'
        '| CL | Series | CD flap 5 | CD flap 10 | Flap delta [counts] | Delta minus experiment [counts] |\n'
        '|---:|---|---:|---:|---:|---:|\n'+comparisons+'\n\n'
        'Experimental and MSES values are digitized from NASA/CR-20220014964, Figures 25(a) and 26(b), '
        'at Re 700,000 and M 0.07. The absolute values and all six historical flap differences were preserved; '
        'the rows above select CL 0.8/1.2. Their hashes and original audited result paths are traced in the input manifest. '
        'MultiFoil Ncrit 9 uses the archived 260/180 protocol. MultiFoil C uses Ncrit 6.25, 440/320 surfaces and fine wake '
        'recipes; the earlier finest operational grid uses finer wake recipes. The Ncrit-9 parameter agrees with published '
        'MSES, but exact geometry/input decks and transition formulations are not identical or fully qualified. '
        'These are explicitly labeled comparisons, not exact code-to-code reproductions. MSES curves are published references, '
        'not a new MSES run. The default remains Ncrit 9; 6.25 remains the frozen S702 calibration.\n\n'
        'At CL 0.8 the current C flap delta remains about 18.33 counts above experiment; the earlier finest operational '
        'grid remains 18.13 counts above. At CL 1.2 the finest prediction still has the wrong sign: +4.53 counts versus '
        '-8.76 measured. This factorial diagnoses a sub-count numerical sensitivity; it does not resolve the large '
        'experimental discrepancy. No uncertainty-based accuracy pass is assigned.\n\n'
        '## Next validation work\n\n'
        'Close this bounded numerical diagnostic and qualify independent experimental inputs next. Extract and visually '
        'check Pomeroy\'s coordinate/placement tables and map measured wake/pressure quantities to supported model outputs '
        'before benchmarking. Preserve S103/S414 as explicit experiment/MultiFoil/published-MSES comparisons, subject to '
        'their tested-geometry and transition-method audits. Do not launch another global grid, closure change or optimizer '
        'to compensate for the unresolved physical discrepancy. No simulation is running.\n\n'
        '![Sampling effects](surface_factors.png)\n\n![Experiment and published MSES](experiment_mses_comparison.png)\n\n'
        '![Geometry and aft layer](geometry_and_aft_layer.png)\n\n'
        'Reproduce with `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 validation/report_s702_surface_factorial.py`. '
        'This reads saved states only. See `verification.json` for checks and `input_manifest.json` for provenance.\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--study',type=Path,default=ROOT/'validation/results/s702_surface_factorial_20261001')
    parser.add_argument('--out',type=Path,default=ROOT/'validation/results/s702_surface_factorial_analysis')
    args=parser.parse_args();print(json.dumps(report(args.study,args.out),indent=2,allow_nan=False))
