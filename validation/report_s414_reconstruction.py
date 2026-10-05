#!/usr/bin/env python3
"""S414 pressure/BL/drag comparison after a complete saved-state audit."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key]='1'
import platform
import importlib.metadata
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from audit_s414_state import ROOT,STUDY,OUT,STUDY_HASH,Inputs,load_states,require,sha,write,close
from report_s702_wake_extent import csv_write,clearance,side_data
from report_s702_matched_pressure import pressure_branch,save
from s414_geometry import SCENARIOS,QUAL as QUAL_REL

QUAL=ROOT/QUAL_REL
DRAG=ROOT/'validation/results/s414_drag_reference_20261002'
CASES=('cl028','cl104')
SURFACES=('fore_upper','fore_lower','aft_upper','aft_lower')
COLORS={'nominal':'#2166ac','source_cubic':'#d6604d','aft_z_minus':'#9970ab','aft_z_plus':'#e6ab02','mses':'#168a6e'}
LABELS={'nominal':'MultiFoil nominal','source_cubic':'Cubic source','aft_z_minus':'Aft −0.13 mm',
        'aft_z_plus':'Aft +0.13 mm','mses':'Published MSES (2010)'}
REGIONS=(('aft_nose',-np.inf,.80),('aft_front_recovery',.80,.86),('aft_bubble_recovery',.86,.94),('aft_tail',.94,np.inf))


def interpolate(x,y,targets):
    x,y,targets=np.asarray(x,float),np.asarray(y,float),np.asarray(targets,float)
    require(len(x)>=2 and x.shape==y.shape and np.isfinite(x).all() and np.isfinite(y).all(), 'Invalid reference curve')
    require(np.all(np.diff(x)>0),'Nonmonotone reference curve')
    require(np.isfinite(targets).all() and np.all((x[0]<=targets)&(targets<=x[-1])),'Reference extrapolation forbidden')
    return np.interp(targets,x,y)


def associate_taps(experiment,coordinates):
    """Lower-curve endpoints are repeats: table index is physical, not plot index."""
    result=[]
    for case in CASES:
        for surface in SURFACES:
            points=[r for r in experiment if r['case']==case and r['surface']==surface]
            table=[r for r in coordinates if r['case']==case and r['surface']==surface]
            table=sorted(table,key=lambda r:int(r['index']));k=0
            for row in points:
                row=dict(row);duplicate=int(row['endpoint_duplicate'])==1
                row.update(table_x=None,table_index=None)
                if not duplicate:
                    require(k<len(table),'Missing physical tap association')
                    close(float(row['x_over_c']),float(table[k]['plotted_x']),context='Plotted/table tap identity')
                    row.update(table_x=float(table[k]['table_x']),table_index=int(table[k]['index']));k+=1
                result.append(row)
            require(k==len(table),'Unused physical tap coordinates')
    return result


def common_support(tap,branches,unsupported):
    if int(tap['endpoint_duplicate']):return 'duplicate_endpoint'
    element,side=tap['surface'].split('_')
    if (element,side,tap['table_index']) in unsupported:return 'unsupported_source_center'
    for coordinate in (tap['table_x'],float(tap['x_over_c'])):
        if any(not x[0]<=coordinate<=x[-1] for x,_ in branches.values()):return 'outside_common_curve_support'
    return ''


def pressure_metric(rows):
    selected=[r for r in rows if r['included']]
    if not selected:return dict(scored_taps=0,rmse=None,bias=None,maximum_absolute_error=None)
    error=np.array([r['error_Cp'] for r in selected])
    return dict(scored_taps=len(selected),rmse=float(np.sqrt(np.mean(error**2))),
        bias=float(error.mean()),maximum_absolute_error=float(abs(error).max()))


def recovery_intervals(x,cp):
    x,cp=np.asarray(x,float),np.asarray(cp,float)
    require(len(x)>=2 and np.all(np.diff(x)>0) and np.isfinite(cp).all(),'Invalid recovery taps')
    return [dict(left_x=float(a),right_x=float(b),delta_Cp=float(d),slope=float(d/(b-a)))
            for a,b,d in zip(x[:-1],x[1:],np.diff(cp))]


def report():
    OUT.mkdir(exist_ok=True);inputs=Inputs();spec,states,controls=load_states(inputs)
    receipt=inputs.json(OUT/'residual_recheck.json')
    require(receipt['complete'] and len(receipt['states'])==9,'Complete residual audit required')
    for name,digest in receipt['input_sha256'].items():
        require(sha(inputs.add(ROOT/name))==digest,'Residual receipt invalidated: '+name)
    checks={r['job']:r for r in receipt['states']}
    reference=inputs.json(DRAG/'source.json');targets_ref={r['case']:r for r in reference['targets']}
    curves={method:inputs.csv(DRAG/(method+'_drag.csv')) for method in ('experiment','mses')}
    pressure_source=inputs.json(QUAL/'pressure/source.json')
    taps=associate_taps(inputs.csv(QUAL/'pressure/experiment.csv'),inputs.csv(QUAL/'pressure/tap_coordinate_checks.csv'))
    mses=inputs.csv(QUAL/'pressure/mses_2010_envelope.csv')
    unsupported={(r['element'],r['surface'],int(r['index'])) for r in inputs.csv(QUAL/'pressure_tap_geometry_support.csv') if int(r['unsupported_center'])}
    require(len(unsupported)==3,'Original unsupported tap set changed')
    groups={c:{} for c in CASES};targets=[];layers=[];separations=[];forces=[];histories=[]
    for job in spec['jobs']:
        name=job['name'];state=states[name];result=state['result'];c=job['pressure_case'];ref=targets_ref[c]
        require(job['target_cl']==ref['CL'],'Measured lift association changed')
        require(checks[name]['max_residual']<=spec['residual_tolerance'],'Residual exceeds tolerance')
        context={}
        for method,curve in curves.items():
            x=[float(r['CL']) for r in curve];y=[float(r['CD']) for r in curve]
            values=interpolate(x,y,[ref['CL']-.0005,ref['CL'],ref['CL']+.0005,result['CL']])
            close(values[1],ref[method+'_CD_interpolated'],context='Reference CD interpolation')
            context.update({method+'_CD':float(values[1]),method+'_CD_at_actual_model_CL':float(values[3]),
                method+'_CD_counts_min_CL_rounding':float(min(values[:3])*1e4),
                method+'_CD_counts_max_CL_rounding':float(max(values[:3])*1e4)})
        nominal=states[c+'_nominal']['result']
        row=dict(job=name,case=c,scenario=job['scenario'],repeat='repeat_of' in job,
            target_cl=job['target_cl'],CL=result['CL'],CD=result['CD'],alpha_deg=result['alpha_deg'],
            error_CD_counts=(result['CD']-context['experiment_CD'])*1e4,
            difference_from_nominal_CD_counts=(result['CD']-nominal['CD'])*1e4,
            mses_minus_experiment_counts=(context['mses_CD']-context['experiment_CD'])*1e4,
            experiment_alpha=ref['alpha_deg'],alpha_error=result['alpha_deg']-ref['alpha_deg'],
            lift_error=result['CL']-job['target_cl'],residual=state['residual'],recomputed_residual=checks[name]['max_residual'],
            selected_path=result['selected_path'],path_matches_nominal=result['selected_path']==nominal['selected_path'],
            requested_evaluations=len(result['records']),failed_requests=sum(not r['admissible'] for r in result['records']),
            failed_internal_points=sum(not p['converged'] for d in state['metadata']['diagnostics'] for p in d['continuation']),
            wall_seconds=result['wall_seconds'],**context)
        targets.append(row);histories.append(dict(job=name,records=result['records'],attempts=result['attempts']))
        if row['repeat']:continue
        groups[c][job['scenario']]=dict(state=state,target=row)
        for element in ('fore','aft'):
            gap=clearance(state['arrays'][element]['wake_xy'].T,state['arrays']['aft' if element=='fore' else 'fore']['xy'].T)
            require(gap>0,'Wake intersects other element')
            forces.append(dict(job=name,case=c,scenario=job['scenario'],element=element,
                other_body_clearance_c=gap,**result['element_forces'][element]))
            for side in ('upper','lower'):
                bl=result['boundary_layers'][element][side];tr=bl['transition'];flow=side_data(state['surfaces'][element],side)
                layers.append(dict(job=name,case=c,scenario=job['scenario'],element=element,side=side,
                    transition_x=None if tr is None else tr['xy_over_cref'][0],min_cf=bl['min_cf'],max_Hk=bl['max_Hk'],
                    separation_regions=len(bl['separation_regions']),separated_fraction=bl['separated_arclength_fraction'],
                    TE_theta=float(flow['theta'][-1]),TE_cf=float(flow['cf'][-1])))
                for i,region in enumerate(bl['separation_regions']):
                    separations.append(dict(job=name,case=c,scenario=job['scenario'],element=element,side=side,interval=i,
                        start_x=region['start']['xy_over_cref'][0],end_x=region['end']['xy_over_cref'][0],
                        reaches_TE=region['reaches_trailing_edge']))
    pressures=[];metrics=[];features=[];recovery=[];excluded=[]
    for c in CASES:
        require(set(groups[c])==set(SCENARIOS),'Incomplete scenario coverage')
        for surface in SURFACES:
            element,side=surface.split('_');branches={}
            for scenario,group in groups[c].items():
                branch=pressure_branch(group['state']['surfaces'][element],side)
                branches[scenario]=(branch['x_over_cref'],branch['Cp'])
            curve=[r for r in mses if r['case']==c and r['surface']==surface]
            branches['mses']=(np.array([float(r['x_over_c']) for r in curve]),np.array([float(r['Cp']) for r in curve]))
            for tap in [r for r in taps if r['case']==c and r['surface']==surface]:
                reason=common_support(tap,branches,unsupported)
                if reason:excluded.append(dict(case=c,surface=surface,plot_index=int(tap['index']),table_index=tap['table_index'],
                    table_x=tap['table_x'],plotted_x=float(tap['x_over_c']),reason=reason))
                for coordinate in ('table','plotted'):
                    x=tap['table_x'] if coordinate=='table' else float(tap['x_over_c'])
                    for method,(bx,by) in branches.items():
                        cp=None if reason else float(interpolate(bx,by,[x])[0]);measured=float(tap['Cp'])
                        pressures.append(dict(case=c,surface=surface,plot_index=int(tap['index']),table_index=tap['table_index'],
                            coordinate=coordinate,x_over_c=x,method=method,experiment_Cp=measured,model_Cp=cp,
                            error_Cp=None if cp is None else cp-measured,included=not reason,exclusion=reason))
            regions=[('all',-np.inf,np.inf)]+(list(REGIONS) if surface=='aft_upper' else [])
            for coordinate in ('table','plotted'):
                for method in branches:
                    rows=[r for r in pressures if r['case']==c and r['surface']==surface and r['coordinate']==coordinate and r['method']==method]
                    for region,lo,hi in regions:
                        selected=[r for r in rows if r['x_over_c'] is not None and lo<=r['x_over_c']<hi]
                        metrics.append(dict(case=c,surface=surface,coordinate=coordinate,method=method,region=region,**pressure_metric(selected)))
                    if surface!='aft_upper':continue
                    selected=[r for r in rows if r['included'] and .80<=r['x_over_c']<=.98]
                    x=[r['x_over_c'] for r in selected]
                    cp=[r['model_Cp'] for r in selected];ex=[r['experiment_Cp'] for r in selected]
                    actual=recovery_intervals(x,cp);observed=recovery_intervals(x,ex)
                    for model,exp in zip(actual,observed):
                        features.append(dict(case=c,coordinate=coordinate,method=method,**model,
                            experiment_delta_Cp=exp['delta_Cp'],experiment_slope=exp['slope'],
                            error_delta_Cp=model['delta_Cp']-exp['delta_Cp']))
                    # Descriptive statistic of all predeclared tap intervals, not an experimental transition endpoint.
                    model_peak=max(actual,key=lambda r:r['slope']);exp_peak=max(observed,key=lambda r:r['slope'])
                    recovery.append(dict(case=c,coordinate=coordinate,method=method,
                        peak_left_x=model_peak['left_x'],peak_right_x=model_peak['right_x'],peak_slope=model_peak['slope'],
                        experimental_peak_left_x=exp_peak['left_x'],experimental_peak_right_x=exp_peak['right_x'],
                        experimental_peak_slope=exp_peak['slope'],
                        peak_interval_midpoint_shift_c=(model_peak['left_x']+model_peak['right_x']-exp_peak['left_x']-exp_peak['right_x'])/2))
    for name,rows in [('targets',targets),('pressure_at_taps',pressures),('pressure_metrics',metrics),
                     ('recovery_intervals',features),('strongest_sampled_recovery',recovery),('excluded_taps',excluded),
                     ('boundary_layers',layers),('separation_intervals',separations),('element_forces_and_wakes',forces)]:
        csv_write(OUT/(name+'.csv'),rows)
    primary=[r for r in targets if not r['repeat']]
    shifts={};diagnosis=[]
    for c in CASES:
        rows=[r for r in primary if r['case']==c]
        shifts[c]=dict(maximum_absolute_CD_change_counts=max(abs(r['difference_from_nominal_CD_counts']) for r in rows),
            maximum_alpha_change_deg=max(abs(r['alpha_deg']-groups[c]['nominal']['target']['alpha_deg']) for r in rows))
        values=[r for r in pressures if r['case']==c and r['surface']=='aft_upper' and r['coordinate']=='table'
                and r['table_index']==11 and r['included']]
        model=[r['model_Cp'] for r in values if r['method'] in SCENARIOS]
        nominal=next(r for r in values if r['method']=='nominal');theory=next(r for r in values if r['method']=='mses')
        transition=[r['transition_x'] for r in layers if r['case']==c and r['element']=='aft' and r['side']=='upper']
        separation=[r for r in separations if r['case']==c and r['element']=='aft' and r['side']=='upper']
        selected_metrics=[r for r in metrics if r['case']==c and r['surface']=='aft_upper' and r['region']=='aft_bubble_recovery']
        metric_delta=[abs(r['rmse']-next(v['rmse'] for v in selected_metrics if v['method']==r['method'] and v['coordinate']=='table'))
            for r in selected_metrics if r['coordinate']=='plotted']
        early=[r for r in features if r['case']==c and r['coordinate']=='table' and r['left_x']==.88215 and r['right_x']==.89313]
        diagnosis.append(dict(case=c,diagnostic_tap_x=nominal['x_over_c'],experiment_Cp=nominal['experiment_Cp'],
            nominal_Cp=nominal['model_Cp'],mses_Cp=theory['model_Cp'],scenario_Cp_min=min(model),scenario_Cp_max=max(model),
            maximum_Cp_change_from_nominal=max(abs(v-nominal['model_Cp']) for v in model),
            experiment_early_recovery_delta_Cp=early[0]['experiment_delta_Cp'],
            nominal_early_recovery_delta_Cp=next(r['delta_Cp'] for r in early if r['method']=='nominal'),
            mses_early_recovery_delta_Cp=next(r['delta_Cp'] for r in early if r['method']=='mses'),
            scenario_early_recovery_delta_Cp_min=min(r['delta_Cp'] for r in early if r['method'] in SCENARIOS),
            scenario_early_recovery_delta_Cp_max=max(r['delta_Cp'] for r in early if r['method'] in SCENARIOS),
            aft_upper_transition_x_min=min(transition),aft_upper_transition_x_max=max(transition),
            aft_upper_separation_start_x_min=min(r['start_x'] for r in separation),
            aft_upper_separation_start_x_max=max(r['start_x'] for r in separation),
            aft_upper_reattachment_x_min=min(r['end_x'] for r in separation),
            aft_upper_reattachment_x_max=max(r['end_x'] for r in separation),
            maximum_recovery_RMSE_change_from_coordinate_convention=max(metric_delta)))
    summary=dict(required_targets=9,converged_targets=9,full_state_audit_pass=True,full_residual_recomputation_pass=True,
        exact_repeat_pass=controls['passed'],exact_saved_array_comparisons=len(controls['array_comparisons']),new_flow_solves=0,
        maximum_lift_error=max(abs(r['lift_error']) for r in targets),maximum_recorded_residual=max(r['residual'] for r in targets),
        maximum_recomputed_residual=max(r['recomputed_residual'] for r in targets),
        requested_evaluations=sum(r['requested_evaluations'] for r in targets),failed_requests=sum(r['failed_requests'] for r in targets),
        failed_internal_points=sum(r['failed_internal_points'] for r in targets),selected_paths=sorted({r['selected_path'] for r in targets}),
        minimum_wake_body_clearance_c=min(r['other_body_clearance_c'] for r in forces),
        physical_tap_observations=180,common_supported_tap_observations=188-len(excluded),
        duplicated_endpoint_observations=sum(r['reason']=='duplicate_endpoint' for r in excluded),
        unsupported_source_observations=sum(r['reason']=='unsupported_source_center' for r in excluded),
        other_support_exclusions=sum(r['reason']=='outside_common_curve_support' for r in excluded),
        scenario_changes=shifts,nominal_targets=[r for r in primary if r['scenario']=='nominal'],
        strongest_sampled_recovery=[r for r in recovery if r['coordinate']=='table'],
        matched_mses_pressure_reference=True,mses_transition_method=pressure_source['mses_transition'],
        exact_tested_geometry_established=False,numerical_convergence_established=False,
        qualification='Robustness only to four declared reconstructions; no geometry fit, Ncrit fit, confidence bounds or proven unique physical cause.')
    require(summary['duplicated_endpoint_observations']==8 and summary['unsupported_source_observations']==6,'Expected exclusions changed')
    write(OUT/'summary.json',summary);write(OUT/'controls.json',controls);write(OUT/'initialization_histories.json',histories)
    write(OUT/'diagnosis.json',dict(cases=diagnosis,
        note='The diagnostic tap and largest-gradient summaries describe the predeclared full interval comparison; no exact experimental transition location is inferred.'))
    plot(groups,taps,mses,pressures,features,targets_ref)
    for name in ('validation/report_s414_reconstruction.py','validation/report_s702_matched_pressure.py'):
        inputs.add(ROOT/name)
    write(OUT/'input_manifest.json',dict(study_sha256=STUDY_HASH,files=dict(inputs.files),analysis_runtime=dict(
        python=platform.python_version(),versions={k:importlib.metadata.version(k) for k in ('numpy','scipy','matplotlib')})))
    print({k:v for k,v in summary.items() if k not in ('nominal_targets','strongest_sampled_recovery')})


def plot(groups,taps,mses,pressures,features,refs):
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    def cp(ax,c,element,side='upper',all_scenarios=True):
        observed=[r for r in taps if r['case']==c and r['surface']==element+'_'+side and not int(r['endpoint_duplicate'])]
        ax.plot([r['table_x'] for r in observed],[float(r['Cp']) for r in observed],
            'o',color='black',mfc='white',ms=4,zorder=5,label='Experiment (Table-II x)')
        theory=[r for r in mses if r['case']==c and r['surface']==element+'_'+side]
        ax.plot([float(r['x_over_c']) for r in theory],[float(r['Cp']) for r in theory],color=COLORS['mses'],lw=1.8,label=LABELS['mses'])
        for scenario in (SCENARIOS if all_scenarios else ('nominal',)):
            b=pressure_branch(groups[c][scenario]['state']['surfaces'][element],side)
            ax.plot(b['x_over_cref'],b['Cp'],color=COLORS[scenario],lw=1.8 if scenario=='nominal' else 1.,
                ls='-' if scenario=='nominal' else '--',label=LABELS[scenario])
        ax.grid(alpha=.2);ax.set(xlabel='x / reference chord',ylabel='Cp');ax.invert_yaxis()
    fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    for col,c in enumerate(CASES):
        cp(axes[0,col],c,'aft');axes[0,col].set(xlim=(.81,.95),ylim=(.05,-.95),title=f'CL {refs[c]["CL"]:.3f} · aft upper pressure')
        axes[0,col].axvspan(.88215,.89313,color='black',alpha=.045)
        for scenario in SCENARIOS:
            state=groups[c][scenario]['state'];flow=side_data(state['surfaces']['aft'],'upper')
            axes[1,col].plot(flow['x_over_cref'],flow['cf'],color=COLORS[scenario],ls='-' if scenario=='nominal' else '--',label=LABELS[scenario])
            tr=state['result']['boundary_layers']['aft']['upper']['transition']
            if tr:axes[1,col].plot(tr['xy_over_cref'][0],0.,'^',color=COLORS[scenario],ms=6)
        axes[1,col].axhline(0,color='.5',lw=.7);axes[1,col].set(xlim=(.81,.95),ylim=(-.0007,.0045),xlabel='x / reference chord',ylabel='Calculated Cf',title='Separation and transition (triangles)');axes[1,col].grid(alpha=.2)
    axes[0,0].legend(fontsize=8);axes[1,0].legend(fontsize=8)
    fig.suptitle('S414: aft recovery remains downstream across tested reconstructions\nExperiment/MultiFoil CL 0.284 / 1.041; published MSES CL 0.28 / 1.04\nRe 1,000,000 · M 0.10 · Ncrit 9 · grey band: strongest measured sampled recovery')
    save(fig,OUT,'aft_pressure_and_separation')
    fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    for row,c in enumerate(CASES):
        for col,element in enumerate(('fore','aft')):
            ax=axes[row,col];cp(ax,c,element,all_scenarios=False)
            ex=[r for r in taps if r['case']==c and r['surface']==element+'_lower' and not int(r['endpoint_duplicate'])]
            ax.plot([r['table_x'] for r in ex],[float(r['Cp']) for r in ex],'s',ms=3,mfc='white',color='black')
            theory=[r for r in mses if r['case']==c and r['surface']==element+'_lower']
            ax.plot([float(r['x_over_c']) for r in theory],[float(r['Cp']) for r in theory],'--',color=COLORS['mses'])
            b=pressure_branch(groups[c]['nominal']['state']['surfaces'][element],'lower');ax.plot(b['x_over_cref'],b['Cp'],'--',color=COLORS['nominal'])
            ax.set(title=f'CL {refs[c]["CL"]:.3f} · {element} (lower: dashed / squares)')
    axes[0,0].legend(fontsize=8);fig.suptitle('S414 pressure comparison · nominal reconstruction\nExperiment/MultiFoil CL 0.284 / 1.041; published MSES CL 0.28 / 1.04\nAll physical taps shown; three unsupported centers excluded from scores')
    save(fig,OUT,'full_pressure_comparison')
    fig,axes=plt.subplots(1,3,figsize=(14,4.7),layout='constrained',gridspec_kw={'width_ratios':[1.15,1,1]})
    for element in ('fore','aft'):
        xy=groups['cl028']['nominal']['state']['arrays'][element]['xy']
        axes[0].plot(xy[0],xy[1],label=element)
    axes[0].set(xlabel='x/c',ylabel='z/c',title='Reconstructed modified test shape');axes[0].set_aspect('equal',adjustable='datalim');axes[0].legend()
    for ax,c in zip(axes[1:],CASES):
        methods=('mses',)+SCENARIOS;ref=refs[c]
        values=[ref['mses_minus_experiment_counts']]+[groups[c][s]['target']['error_CD_counts'] for s in SCENARIOS]
        bars=ax.bar(range(5),values,color=[COLORS[s] for s in methods]);ax.bar_label(bars,fmt='%+.2f',padding=3,fontsize=9)
        ax.set_xticks(range(5),['MSES','Nominal','Cubic','Aft −','Aft +'],rotation=25)
        ax.axhline(0,color='black',lw=.8);ax.set(title=f'CL {ref["CL"]:.3f}',ylabel='CD minus experiment (counts)',ylim=(-17,9));ax.grid(axis='y',alpha=.2)
    fig.suptitle('S414 matched-lift drag · unchanged model, Ncrit 9\nGeometry variations are declared sensitivity cases, not uncertainty bounds')
    save(fig,OUT,'geometry_and_drag')
    fig,axes=plt.subplots(1,2,figsize=(12,4.5),layout='constrained')
    for ax,c in zip(axes,CASES):
        for method in ('mses',)+SCENARIOS:
            selected=[r for r in features if r['case']==c and r['coordinate']=='table' and r['method']==method]
            x=[.5*(r['left_x']+r['right_x']) for r in selected]
            ax.plot(x,[r['delta_Cp'] for r in selected],'-o',ms=3,color=COLORS[method],label=LABELS[method])
        ax.plot(x,[r['experiment_delta_Cp'] for r in selected],'ko-',ms=4,label='Experiment')
        ax.axhline(0,color='.5',lw=.7);ax.set(title=f'CL {refs[c]["CL"]:.3f}',xlabel='Midpoint between successive taps (x/c)',ylabel='Cp rise between taps');ax.grid(alpha=.2)
    axes[0].legend(fontsize=8);fig.suptitle('Measured recovery occurs one sampled interval earlier\nInterval midpoints are descriptive; they are not transition or reattachment measurements')
    save(fig,OUT,'pressure_recovery_intervals')


if __name__=='__main__':report()
