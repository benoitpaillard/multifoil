#!/usr/bin/env python3
"""Audit and plot S702 pressure snapshots, saved BL states and published drag.

No new flow solves. Requires audit_s702_saved_pressure_state.py's complete,
hash-checked residual/force re-evaluation receipt.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key]='1'
import argparse
import json
from pathlib import Path
import sys
import platform
import importlib.metadata
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from report_s702_wake_extent import ROOT, Inputs, audit, sha, write, csv_write, require, close, clearance, side_data

STUDY=ROOT/'validation/results/s702_matched_pressure_20261002'
QUAL=ROOT/'validation/results/s702_pressure_qualification_20261002'
OUT=ROOT/'validation/results/s702_matched_pressure_analysis_20261002'
STUDY_HASH='786752736ff78e82005be0fb22009bbdba50015db57d595ab33bc51126bb26c4'
COLORS={9.:'#2364aa',6.25:'#d55e00'}
REGIONS=(('aft_front_recovery',.80,.86),('hinge',.86,.94),('aft_tail',.94,1.1))


def pressure_branch(surface,side):
    """Pressure taps refer to geometric surfaces, not stagnation-based BL branches.

    The saved contour runs upper TE -> nose -> lower TE. Use minimum global x
    for these non-overhanging S702 branches; reject nonmonotone branches rather
    than silently sorting across an overhang. Cf/transition retain flow branches.
    """
    require(side in ('upper','lower'),'Unknown geometric surface')
    nose=int(np.argmin(surface['x_over_cref']))
    branch=surface[:nose+1][::-1] if side=='upper' else surface[nose:]
    require(len(branch)>=2 and np.all(np.diff(branch['x_over_cref'])>0),
            'Nonmonotone geometric branch requires a different tap mapping')
    return branch


def at_taps(branch,experiment):
    rows=[]
    for r in experiment:
        x=float(r['x_over_cref']);duplicate=int(r['endpoint_duplicate'])==1
        inside=branch['x_over_cref'][0]<=x<=branch['x_over_cref'][-1]
        reason='duplicate_endpoint' if duplicate else ('outside_model_support' if not inside else '')
        value=None if reason else float(np.interp(x,branch['x_over_cref'],branch['Cp']))
        cp=float(r['Cp'])
        rows.append(dict(drawing_index=int(r['drawing_index']),x_over_cref=x,experiment_Cp=cp,
            model_Cp=value,error_Cp=None if value is None else value-cp,included=not reason,exclusion=reason))
    return rows


def metric(rows):
    selected=[r for r in rows if r['included']]
    require(bool(selected),'Cannot score a surface without supported taps')
    error=np.array([r['error_Cp'] for r in selected])
    return dict(scored_taps=len(selected),total_plotted_taps=len(rows),
        excluded_taps=len(rows)-len(selected),bias=float(error.mean()),rmse=float(np.sqrt(np.mean(error**2))),
        maximum_absolute_error=float(np.max(abs(error))))


def report(out):
    inputs=Inputs();require(sha(inputs.add(STUDY/'study.json'))==STUDY_HASH,'Unexpected pressure study')
    spec=inputs.json(STUDY/'study.json')
    for name,digest in spec['sha256'].items():
        require(sha(inputs.add(STUDY/'source_snapshot'/name))==digest,'Frozen snapshot differs: '+name)
        require(sha(inputs.add(ROOT/name))==digest,'Current frozen input differs: '+name)
    retrieval=inputs.json(ROOT/'validation/results/s702_matched_pressure_launch_20261002/retrieved_files_sha256.json')
    for name,digest in retrieval.items():
        require(sha(inputs.add(STUDY/name))==digest,'Retrieved output changed: '+name)
    source=inputs.json(QUAL/'source.json');cases={c['name']:c for c in source['cases']}
    context={r['case']:r for r in inputs.csv(QUAL/'snapshot_features_and_drag_context.csv')}
    for c in cases.values():
        for kind in ('experiment','mses'):
            drag=inputs.csv(ROOT/f'validation/data/s702/d+{c["flap"]}/{kind}_drag.csv')
            drag=sorted(drag,key=lambda r:float(r['CL']))
            require(float(drag[0]['CL'])<=c['target_cl']<=float(drag[-1]['CL']),'Drag context outside support')
            expected=float(np.interp(c['target_cl'],[float(r['CL']) for r in drag],[float(r['CD']) for r in drag]))
            close(expected,float(context[c['name']][kind+'_CD']),context='Recomputed published drag context')
    runtime=inputs.json(STUDY/'runtime.json');status=inputs.json(STUDY/'status.json');automatic=inputs.json(STUDY/'summary.json')
    require(runtime['versions']==spec['expected_versions'] and runtime['python']=='3.12.13','Solver runtime changed')
    require(set(runtime['thread_limits'].values())=={'1'},'Thread policy changed')
    require(status['state']=='finished' and status['completed_jobs']==status['total_jobs']==8,'Incomplete execution')
    require(automatic['converged_jobs']==automatic['required_jobs']==8 and automatic['execution_checks_pass'],'Coverage failed')
    recheck=inputs.json(out/'residual_recheck.json')
    require(recheck['complete'] and len(recheck['states'])==8,'Full residual recheck is required')
    for name,digest in recheck['input_sha256'].items():
        require(sha(inputs.add(ROOT/name))==digest,'Recomputed-state receipt invalidated: '+name)
    states={};targets=[];pressures=[];metrics=[];features=[];layers=[];separations=[];histories=[];wakes=[]
    experiment={}
    for name in cases:
        experiment[name]={}
        for e in ('fore','aft'):
            for side in ('upper','lower'):
                experiment[name][e,side]=inputs.array(QUAL/name/(e+'_'+side+'.csv'))
    for j in spec['jobs']:
        s=audit(STUDY/j['name'],j,dict(spec,ncrit=j['ncrit']),runtime,inputs)
        r=s['result'];require(r['feasible'],'Incomplete pressure comparison');states[j['name']]=s
        c=cases[j['pressure_case']];ctx=context[c['name']]
        require(j['target_cl']==c['target_cl'],'Pressure/lift association changed')
        ex_cd=float(ctx['experiment_CD']);mses_cd=float(ctx['mses_CD'])
        auto=next(v for v in automatic['jobs'] if v['job']==j['name'])
        for key in ('CL','CD','alpha_deg'):close(r[key],auto[key],context='Execution summary '+key)
        check=next(v for v in recheck['states'] if v['job']==j['name'])
        require(check['max_residual']<=spec['residual_tolerance'],'Recomputed residual fails tolerance')
        targets.append(dict(job=j['name'],case=c['name'],flap=j['flap'],ncrit=j['ncrit'],target_cl=j['target_cl'],
            CL=r['CL'],CD=r['CD'],CM=r['CM'],experiment_CD=ex_cd,published_mses_n9_CD=mses_cd,
            model_drag_error_counts=(r['CD']-ex_cd)*1e4,mses_drag_error_counts=(mses_cd-ex_cd)*1e4,
            model_alpha_deg=r['alpha_deg'],experiment_alpha_deg=c['alpha'],alpha_error_deg=r['alpha_deg']-c['alpha'],
            lift_error=r['CL']-j['target_cl'],recorded_residual=s['residual'],recomputed_residual=check['max_residual'],
            selected_path=r['selected_path'],requested_evaluations=len(r['records']),
            failed_requests=sum(not v['admissible'] for v in r['records']),
            failed_internal_points=sum(not p['converged'] for d in s['metadata']['diagnostics'] for p in d['continuation']),
            wall_seconds=r['wall_seconds']))
        histories.append(dict(job=j['name'],selected_flow_index=r['selected_flow_index'],records=r['records'],
            attempts=r['attempts'],continuation=[d['continuation'] for d in s['metadata']['diagnostics']]))
        for e in ('fore','aft'):
            other='aft' if e=='fore' else 'fore';gap=clearance(s['arrays'][e]['wake_xy'].T,s['arrays'][other]['xy'].T)
            require(gap>0,'Wake intersects the other body')
            wakes.append(dict(job=j['name'],element=e,other_body_clearance=gap,**r['element_forces'][e]))
            for side in ('upper','lower'):
                branch=pressure_branch(s['surfaces'][e],side)
                tap=at_taps(branch,experiment[c['name']][e,side])
                rows=[dict(job=j['name'],case=c['name'],element=e,side=side,**v) for v in tap]
                pressures.extend(rows)
                metrics.append(dict(job=j['name'],case=c['name'],element=e,side=side,region='all',**metric(rows)))
                if e=='aft' and side=='upper':
                    for region,lo,hi in REGIONS:
                        subset=[v for v in rows if lo<=v['x_over_cref']<hi]
                        metrics.append(dict(job=j['name'],case=c['name'],element=e,side=side,region=region,**metric(subset)))
                    a,b=[float(ctx[key]) for key in ('prehinge_left_x','prehinge_right_x')]
                    val=np.interp([a,b],branch['x_over_cref'],branch['Cp']);difference=float(val[1]-val[0])
                    features.append(dict(job=j['name'],case=c['name'],flap=j['flap'],ncrit=j['ncrit'],
                        left_x=a,right_x=b,experiment_delta_Cp=float(ctx['prehinge_delta_Cp']),
                        model_delta_Cp=difference,error_delta_Cp=difference-float(ctx['prehinge_delta_Cp'])))
                flow=side_data(s['surfaces'][e],side);bl=r['boundary_layers'][e][side];tr=bl['transition']
                layers.append(dict(job=j['name'],case=c['name'],element=e,side=side,
                    transition_x_over_cref=None if tr is None else tr['xy_over_cref'][0],
                    transition_x_over_local_chord=None if tr is None else tr['x_over_local_chord'],
                    separation_regions=len(bl['separation_regions']),min_cf=bl['min_cf'],max_Hk=bl['max_Hk'],
                    TE_theta_over_cref=float(flow['theta'][-1]),TE_delta_star_over_cref=float(flow['delta_star'][-1]),
                    TE_Hk=float(flow['Hk'][-1]),TE_cf=float(flow['cf'][-1])))
                for k,reg in enumerate(bl['separation_regions']):
                    separations.append(dict(job=j['name'],case=c['name'],element=e,side=side,interval=k,
                        start_x=reg['start']['xy_over_cref'][0],end_x=reg['end']['xy_over_cref'][0],
                        start_local_x=reg['start']['x_over_local_chord'],end_local_x=reg['end']['x_over_local_chord'],
                        reaches_TE=reg['reaches_trailing_edge']))
        print(j['name']+': selected-state audit and tap comparison passed',flush=True)
    excluded=[r for r in pressures if not r['included']]
    unsupported=[r for r in excluded if r['exclusion']=='outside_model_support']
    require(len(unsupported)==8 and all(r['element']=='aft' and r['side']=='upper' and .699<r['x_over_cref']<.700 for r in unsupported),
            'Unexpected experimental support exclusions; inspect before reporting')
    from s702_geometry import geometry
    geometry_rows=[]
    for c in cases.values():
        raw=dict(geometry(c['flap']))['aft'];raw_min=float(raw[:,0].min())
        j=next(j for j in spec['jobs'] if j['pressure_case']==c['name'])
        actual=states[j['name']]['arrays']['aft']['xy'];actual_min=float(actual[0].min())
        first=experiment[c['name']]['aft','upper'][0]
        # Table II(b), PDF 26: x/c of the foremost upper tap. Its z coordinate
        # is not tabulated (the table's y is spanwise), so do not invent a point.
        printed=.69941
        geometry_rows.append(dict(case=c['name'],flap=c['flap'],table_ii_aft_upper_first_x=printed,
            plotted_first_tap_x=float(first['x_over_cref']),raw_contour_minimum_x=raw_min,
            realized_solver_minimum_x=actual_min,raw_nose_minus_table_tap_c=raw_min-printed,
            raw_nose_minus_table_tap_mm=(raw_min-printed)*457.48,
            solver_nose_minus_plotted_tap_c=actual_min-float(first['x_over_cref']),
            paneling_change_in_minimum_x=actual_min-raw_min))
    for name,table in (('targets',targets),('pressure_at_taps',pressures),('pressure_metrics',metrics),
                       ('pressure_features',features),('boundary_layers',layers),('separation_intervals',separations),
                       ('excluded_taps',excluded),('element_forces_and_wakes',wakes),('geometry_support',geometry_rows)):
        csv_write(out/(name+'.csv'),table)
    summary=dict(required_targets=8,converged_targets=8,selected_states_audited=8,full_state_audit_pass=True,
        full_residual_recomputation_pass=True,new_flow_solves=0,matched_reynolds=1e6,matched_mach=.1,
        max_absolute_lift_error=max(abs(r['lift_error']) for r in targets),
        max_recorded_residual=max(r['recorded_residual'] for r in targets),
        max_recomputed_residual=max(r['recomputed_residual'] for r in targets),
        requested_evaluations=sum(r['requested_evaluations'] for r in targets),
        failed_requests=sum(r['failed_requests'] for r in targets),failed_internal_points=sum(r['failed_internal_points'] for r in targets),
        selected_paths=sorted({r['selected_path'] for r in targets}),
        scored_tap_comparisons=sum(r['included'] for r in pressures),duplicate_endpoint_exclusions=len(excluded)-len(unsupported),
        outside_support_exclusions=len(unsupported),minimum_wake_body_clearance=min(r['other_body_clearance'] for r in wakes),
        aft_TE_separation=any(r['reaches_TE'] for r in separations if r['element']=='aft'),
        numerical_convergence_established=False,experimental_transition_locations_available=False,
        matched_mses_pressure_available=False,qualification='Pressure/BL mechanism evidence on calibrated S702; no causal proof or independent geometry validation.')
    summary.update(matched_as_built_geometry_established=False,
        raw_aft_nose_minus_printed_foremost_tap_c=geometry_rows[0]['raw_nose_minus_table_tap_c'])
    write(out/'summary.json',summary);write(out/'initialization_histories.json',histories)
    plot(out,spec,cases,states,experiment,targets,features)
    for name in ('validation/report_s702_matched_pressure.py','validation/audit_s702_saved_pressure_state.py',
        'validation/report_s702_wake_extent.py','validation/report_s702_drag_mechanism.py','validation/report_design_initialization.py'):
        inputs.add(ROOT/name)
    write(out/'input_manifest.json',dict(study_sha256=STUDY_HASH,files=inputs.files,
        analysis_runtime=dict(python=platform.python_version(),versions={k:importlib.metadata.version(k) for k in ('numpy','scipy','matplotlib')})))
    return summary


def save(fig,out,name):
    for ext in ('png','svg','pdf'):fig.savefig(out/(name+'.'+ext),dpi=180)
    plt.close(fig)


def plot(out,spec,cases,states,experiment,targets,features):
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    def jobs(c):return [j for j in spec['jobs'] if j['pressure_case']==c['name']]
    def title(c):return f"Flap {c['flap']}° · CL ≈ {c['target_cl']:.3f} · experiment α = {c['alpha']:.2f}°"
    def cp(ax,c,e,side,label=True):
        ex=experiment[c['name']][e,side];ex=ex[ex['endpoint_duplicate']==0]
        ax.plot(ex['x_over_cref'],ex['Cp'],'o',color='black',ms=4,mfc='white',zorder=4,label='Experiment' if label else None)
        for j in jobs(c):
            v=pressure_branch(states[j['name']]['surfaces'][e],side)
            ax.plot(v['x_over_cref'],v['Cp'],color=COLORS[j['ncrit']],lw=1.5,
                ls='-' if side=='upper' else '--',label=f"MultiFoil Ncrit {j['ncrit']:g}" if label else None)
        ax.grid(alpha=.18);ax.set_ylabel('$C_p$')
    fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    for ax,c in zip(axes.flat,cases.values()):
        cp(ax,c,'aft','upper');ax.axvline(.9001,color='.5',ls=':',lw=1)
        ax.axvspan(.86,.94,color='.6',alpha=.08)
        ax.set_xlim(.80,.975);ax.set_ylim(.15,-1.05);ax.set_xlabel('x / reference chord')
        ax.set_title(title(c),fontsize=11);ax.legend(fontsize=8,loc='upper right')
    fig.suptitle('Aft upper pressure recovery at matched lift\nRe = 1,000,000 · M = 0.10 · dots are measured taps; dotted line is the hinge',fontsize=14)
    save(fig,out,'aft_pressure_recovery')
    fig,axes=plt.subplots(4,2,figsize=(12,15),layout='constrained')
    for pair,c in zip(axes,cases.values()):
        for ax,e in zip(pair,('fore','aft')):
            cp(ax,c,e,'upper');cp(ax,c,e,'lower',False)
            ax.invert_yaxis();ax.set_xlabel('x / reference chord');ax.set_title(title(c)+'\n'+e+' element',fontsize=10)
            if e=='aft':ax.axvline(.9001,color='.5',ls=':',lw=1)
            ax.legend(fontsize=8)
    fig.suptitle('Full pressure comparison — geometric upper/lower surfaces\nSolid: upper · dashed: lower · no measured-pressure interpolation between operating points',fontsize=13)
    save(fig,out,'full_pressure')
    fig,axes=plt.subplots(4,3,figsize=(15,13),layout='constrained')
    for trio,c in zip(axes,cases.values()):
        ax=trio[0];cp(ax,c,'aft','upper');ax.set_xlim(.80,.99);ax.set_ylim(.2,-1.05)
        ax.set_title(title(c),fontsize=9);ax.axvline(.9001,color='.5',ls=':',lw=1)
        for j in jobs(c):
            s=states[j['name']];n=j['ncrit'];color=COLORS[n]
            for side,style in (('upper','-'),('lower','--')):
                v=side_data(s['surfaces']['aft'],side);bl=s['result']['boundary_layers']['aft'][side]
                trio[1].plot(v['x_over_cref'],v['cf'],ls=style,color=color,label=f'N {n:g}, {side}')
                if bl['transition']:
                    x=bl['transition']['xy_over_cref'][0]
                    trio[1].plot(x,np.interp(x,v['x_over_cref'],v['cf']),marker='v' if side=='upper' else '^',color=color,ms=5)
                # U[2] changes meaning after transition. Plot amplification only
                # on laminar nodes, normalized by each declared Ncrit.
                lam=v[v['turbulent']==0]
                trio[2].plot(lam['x_over_cref'],lam['amplification_or_shear']/n,ls=style,color=color)
                for reg in bl['separation_regions']:
                    trio[1].axvspan(reg['start']['xy_over_cref'][0],reg['end']['xy_over_cref'][0],color=color,alpha=.06)
        for ax in trio:
            ax.grid(alpha=.18);ax.set_xlabel('x / reference chord')
        trio[1].set(xlim=(.78,1.0),ylim=(-.001,.0045),ylabel='$C_f$ referenced to freestream')
        trio[1].axhline(0,color='.3',lw=.7);trio[1].legend(fontsize=7,ncol=2)
        trio[2].set(xlim=(.78,1.0),ylim=(0,1.05),ylabel='Laminar amplification / Ncrit')
        trio[2].axhline(1,color='.4',ls=':',lw=1)
    fig.suptitle('Calculated aft boundary layers alongside measured pressure\nNegative Cf: separation · triangles: calculated transition · no measured transition locations',fontsize=13)
    save(fig,out,'aft_boundary_layers')
    fig,axes=plt.subplots(2,2,figsize=(12,7.4),layout='constrained')
    from s702_geometry import geometry
    for ax,flap in zip(axes[0],(5,10)):
        for e,xy in geometry(flap):ax.plot(xy[:,0],xy[:,1],label=e)
        ax.plot(.9001,-.0114,'k+',ms=9,label='hinge')
        ax.set_aspect('equal',adjustable='box');ax.set(xlim=(-.02,1.03),ylim=(-.10,.15),xlabel='x / c',ylabel='z / c',title=f'Reconstructed S702 geometry — flap {flap}°')
        ax.legend(fontsize=8);ax.grid(alpha=.2)
    names=list(cases);x=np.arange(4);width=.19
    for offset,(label,color,n) in enumerate((('Experiment','black',None),('Published MSES N9','#009e73','mses'),('MultiFoil N9',COLORS[9.],9.),('MultiFoil N6.25',COLORS[6.25],6.25))):
        vals=[]
        for name in names:
            r=next(t for t in targets if t['case']==name and (not isinstance(n,float) or t['ncrit']==n))
            vals.append(r['experiment_CD'] if n is None else r['published_mses_n9_CD'] if n=='mses' else r['CD'])
        axes[1,0].bar(x+(offset-1.5)*width,np.array(vals)*1e4,width,color=color,label=label)
    axes[1,0].set_ylabel('Drag counts (10⁴ CD)');axes[1,0].set_title('Drag at each snapshot’s own CL');axes[1,0].legend(fontsize=8)
    for offset,(label,color,n) in enumerate((('Experiment','black',None),('MultiFoil N9',COLORS[9.],9.),('MultiFoil N6.25',COLORS[6.25],6.25))):
        vals=[]
        for name in names:
            r=next(t for t in features if t['case']==name and (n is None or t['ncrit']==n))
            vals.append(r['experiment_delta_Cp'] if n is None else r['model_delta_Cp'])
        axes[1,1].bar(x+(offset-1)*.24,vals,.24,color=color,label=label)
    axes[1,1].set_ylabel('Cp(≈0.889c) − Cp(≈0.872c)');axes[1,1].set_title('Measured pre-hinge pressure-rise feature');axes[1,1].legend(fontsize=8)
    for ax in axes[1]:
        ax.set_xticks(x,[f"{cases[name]['flap']}°\nCL {cases[name]['target_cl']:.3f}" for name in names]);ax.grid(axis='y',alpha=.18)
    fig.suptitle('Pressure and drag give different tests of the model\nRe = 1,000,000 · M = 0.10 · unequal-CL snapshots are not a flap-drag difference',fontsize=13)
    save(fig,out,'geometry_pressure_and_drag')
    fig,axes=plt.subplots(1,2,figsize=(10,5),layout='constrained')
    for ax,flap in zip(axes,(5,10)):
        raw=dict(geometry(flap))['aft']
        j=next(j for j in spec['jobs'] if j['flap']==flap)
        xy=states[j['name']]['arrays']['aft']['xy'].T
        ax.plot(raw[:,0],raw[:,1],'-o',ms=3,lw=1,label='Figure-1 contour input')
        ax.plot(xy[:,0],xy[:,1],'.',ms=3,color='#d55e00',label='Realized solver nodes')
        ax.axvline(.69941,color='black',ls='--',lw=1,label='Foremost tap x/c (Table II)')
        ax.set(xlim=(.698,.707),ylim=(-.034,-.011),xlabel='x / c',ylabel='z / c',title=f'Flap {flap}° — aft nose')
        ax.grid(alpha=.2);ax.legend(fontsize=8,loc='upper left')
    fig.suptitle('A pressure tap lies ahead of the reconstructed contour\nOnly the tap’s x is known here; the vertical line does not specify its height',fontsize=12)
    save(fig,out,'aft_nose_geometry_support')


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--out',type=Path,default=OUT)
    args=ap.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    print(json.dumps(report(args.out),indent=2))
