#!/usr/bin/env python3
"""Report the complete 2x2 transition diagnostic against measured S702 pressure."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key]='1'
import argparse
import importlib.metadata
import platform
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from audit_s702_transition_state import (ROOT,STUDY,PARENT,OUT,STUDY_HASH,Inputs,
    load_states,thresholds,require,close,sha,write,side_data)
from report_s702_matched_pressure import pressure_branch,at_taps,metric,REGIONS,save
from report_s702_wake_extent import csv_write,clearance
from s702_geometry import geometry

QUAL=ROOT/'validation/results/s702_pressure_qualification_20261002'
KEYS=((9.,9.),(6.25,9.),(9.,6.25),(6.25,6.25))
LABELS={KEYS[0]:'Both 9',KEYS[1]:'Fore 6.25, aft 9',KEYS[2]:'Fore 9, aft 6.25',KEYS[3]:'Both 6.25'}
COLORS={KEYS[0]:'#2364aa',KEYS[1]:'#8e44ad',KEYS[2]:'#009e73',KEYS[3]:'#d55e00'}
STYLES={KEYS[0]:'-',KEYS[1]:'--',KEYS[2]:'-.',KEYS[3]:'-'}


def factorial(values):
    require(set(values)==set(KEYS),'A complete four-state factorial is required')
    require(all(np.isfinite(v) for v in values.values()),'Nonfinite factorial endpoint')
    q00,q10,q01,q11=[float(values[k]) for k in KEYS]
    return dict(baseline=q00,fore_only=q10,aft_only=q01,both=q11,
        fore_effect=q10-q00,aft_effect=q01-q00,interaction=q11-q10-q01+q00,
        conditional_fore_effect=q11-q01,conditional_aft_effect=q11-q10,combined_change=q11-q00)


def report(out):
    out.mkdir(exist_ok=True);inputs=Inputs()
    spec,entries,states,controls=load_states(inputs)
    receipt=inputs.json(out/'residual_recheck.json')
    require(receipt['complete'] and len(receipt['states'])==18,'Complete residual receipt is required')
    for name,digest in receipt['input_sha256'].items():
        require(sha(inputs.add(ROOT/name))==digest,'Residual receipt invalidated: '+name)
    checks={r['job']:r for r in receipt['states']}
    source=inputs.json(QUAL/'source.json');cases={c['name']:c for c in source['cases']}
    context={r['case']:r for r in inputs.csv(QUAL/'snapshot_features_and_drag_context.csv')}
    experiment={case:{(e,side):inputs.array(QUAL/case/(e+'_'+side+'.csv'))
        for e in ('fore','aft') for side in ('upper','lower')} for case in cases}
    for c in cases.values():
        for kind in ('experiment','mses'):
            drag=sorted(inputs.csv(ROOT/f'validation/data/s702/d+{c["flap"]}/{kind}_drag.csv'),key=lambda r:float(r['CL']))
            require(float(drag[0]['CL'])<=c['target_cl']<=float(drag[-1]['CL']),'Drag context outside support')
            cd=np.interp(c['target_cl'],[float(r['CL']) for r in drag],[float(r['CD']) for r in drag])
            close(cd,float(context[c['name']][kind+'_CD']),context='Published drag interpolation')
    groups={case:{} for case in cases};targets=[];pressures=[];metrics=[];features=[];layers=[];separations=[];elements=[];histories=[]
    for entry in entries:
        j=entry['job'];s=states[j['name']];r=s['result'];n=thresholds(j);key=(n['fore'],n['aft'])
        c=cases[j['pressure_case']];ctx=context[c['name']];control=j.get('control',False)
        require(c['target_cl']==j['target_cl'],'Experimental lift association changed')
        check=checks[j['name']];require(check['ncrit_by_element']==n,'Residual audit used wrong thresholds')
        require(check['max_residual']<=spec['residual_tolerance'],'Residual exceeds tolerance')
        prefix=dict(job=j['name'],case=c['name'],fore_ncrit=n['fore'],aft_ncrit=n['aft'])
        target=dict(**prefix,flap=j['flap'],control=control,new=entry['new'],target_cl=j['target_cl'],CL=r['CL'],CD=r['CD'],
            model_CD_counts=r['CD']*1e4,experiment_CD_counts=float(ctx['experiment_CD'])*1e4,
            published_mses_n9_CD_counts=float(ctx['mses_CD'])*1e4,
            error_CD_counts=(r['CD']-float(ctx['experiment_CD']))*1e4,
            alpha_deg=r['alpha_deg'],experiment_alpha_deg=c['alpha'],alpha_error_deg=r['alpha_deg']-c['alpha'],
            lift_error=r['CL']-j['target_cl'],residual=s['residual'],recomputed_residual=check['max_residual'],
            residual_evidence=check['evidence'],selected_path=r['selected_path'],requested_evaluations=len(r['records']),
            failed_requests=sum(not v['admissible'] for v in r['records']),
            failed_internal_points=sum(not v['converged'] for d in s['metadata']['diagnostics'] for v in d['continuation']),
            wall_seconds=r['wall_seconds'])
        targets.append(target);histories.append(dict(**prefix,records=r['records'],attempts=r['attempts']))
        if control:continue
        require(key not in groups[c['name']],'Duplicate factorial endpoint')
        quantities=dict(CD_counts=r['CD']*1e4,alpha_deg=r['alpha_deg'])
        groups[c['name']][key]=dict(state=s,target=target,quantities=quantities)
        for e in ('fore','aft'):
            gap=clearance(s['arrays'][e]['wake_xy'].T,s['arrays']['aft' if e=='fore' else 'fore']['xy'].T)
            require(gap>0,'Wake intersects the other body')
            elements.append(dict(**prefix,element=e,other_body_clearance=gap,**r['element_forces'][e]))
            quantities[e+'_CD_counts']=r['element_forces'][e]['CD_wake']*1e4
            for side in ('upper','lower'):
                branch=pressure_branch(s['surfaces'][e],side);ex=experiment[c['name']][e,side]
                rows=[dict(**prefix,element=e,side=side,**v) for v in at_taps(branch,ex)];pressures.extend(rows)
                regions=[('all',-np.inf,np.inf)]+(list(REGIONS) if (e,side)==('aft','upper') else [])
                for region,lo,hi in regions:
                    values=metric([v for v in rows if lo<=v['x_over_cref']<hi])
                    metrics.append(dict(**prefix,element=e,side=side,region=region,**values))
                    for m in ('rmse','bias'):quantities[f'Cp_{e}_{side}_{region}_{m}']=values[m]
                bl=r['boundary_layers'][e][side];tr=bl['transition'];flow=side_data(s['surfaces'][e],side)
                layers.append(dict(**prefix,element=e,side=side,transition_x=None if tr is None else tr['xy_over_cref'][0],
                    transition_local_x=None if tr is None else tr['x_over_local_chord'],min_cf=bl['min_cf'],max_Hk=bl['max_Hk'],
                    separation_regions=len(bl['separation_regions']),TE_theta=float(flow['theta'][-1]),TE_cf=float(flow['cf'][-1])))
                for i,reg in enumerate(bl['separation_regions']):
                    separations.append(dict(**prefix,element=e,side=side,interval=i,start_x=reg['start']['xy_over_cref'][0],
                        end_x=reg['end']['xy_over_cref'][0],reaches_TE=reg['reaches_trailing_edge']))
                if (e,side)==('aft','upper'):
                    a,b=[float(ctx[k]) for k in ('prehinge_left_x','prehinge_right_x')]
                    val=np.interp([a,b],branch['x_over_cref'],branch['Cp']);q=float(val[1]-val[0])
                    quantities['prehinge_delta_Cp']=q
                    features.append(dict(**prefix,left_x=a,right_x=b,experiment_delta_Cp=float(ctx['prehinge_delta_Cp']),
                        model_delta_Cp=q,error_delta_Cp=q-float(ctx['prehinge_delta_Cp'])))
    effects=[]
    for case,group in groups.items():
        require(set(group)==set(KEYS),'Incomplete factorial coverage')
        for quantity in group[KEYS[0]]['quantities']:
            f=factorial({k:group[k]['quantities'][quantity] for k in KEYS})
            close(f['fore_effect']+f['aft_effect']+f['interaction'],f['combined_change'],context='Factorial sum')
            effects.append(dict(case=case,quantity=quantity,all_selected_paths_match=len({v['target']['selected_path'] for v in group.values()})==1,**f))
    excluded=[r for r in pressures if not r['included']];unsupported=[r for r in excluded if r['exclusion']=='outside_model_support']
    require(len(unsupported)==16 and all(r['element']=='aft' and r['side']=='upper' and .699<r['x_over_cref']<.700 for r in unsupported),'Experimental support changed')
    require(sum(r['exclusion']=='duplicate_endpoint' for r in excluded)==64,'Duplicate exclusion set changed')
    for name,rows in [('targets',targets),('factorial_effects',effects),('pressure_at_taps',pressures),('pressure_metrics',metrics),
            ('pressure_features',features),('boundary_layers',layers),('separation_intervals',separations),
            ('element_forces_and_wakes',elements),('excluded_taps',excluded)]:csv_write(out/(name+'.csv'),rows)
    new=[r for r in targets if r['new']]
    summary=dict(required_new_targets=10,converged_new_targets=10,selected_states_audited=18,factorial_endpoints=16,
        full_state_audit_pass=True,full_residual_audit_pass=True,fresh_residual_recomputations=10,reused_hash_verified_residual_rechecks=8,
        new_flow_solves=0,controls_pass=True,control_saved_arrays_exact=True,
        maximum_lift_error=max(abs(r['lift_error']) for r in targets),maximum_residual=max(r['recomputed_residual'] for r in targets),
        new_requested_evaluations=sum(r['requested_evaluations'] for r in new),new_failed_requests=sum(r['failed_requests'] for r in new),
        new_failed_internal_points=sum(r['failed_internal_points'] for r in new),selected_paths=sorted({r['selected_path'] for r in targets}),
        all_factorial_paths_match=all(r['all_selected_paths_match'] for r in effects),scored_tap_comparisons=sum(r['included'] for r in pressures),
        duplicate_endpoint_exclusions=64,outside_support_exclusions=16,minimum_wake_body_clearance=min(r['other_body_clearance'] for r in elements),
        aft_TE_separation=any(r['reaches_TE'] for r in separations if r['element']=='aft'),
        numerical_convergence_established=False,exact_tested_geometry_established=False,transition_causality_established=False,
        matched_mses_Cp_available=False,qualification='Element-specific threshold sensitivity on the fixed reconstructed geometry; no Ncrit fit or independent design validation.')
    write(out/'summary.json',summary);write(out/'controls.json',controls);write(out/'initialization_histories.json',histories)
    plot(out,cases,groups,experiment,context,effects)
    for name in ('validation/report_s702_transition_isolation.py','validation/report_s702_matched_pressure.py',
                 'validation/audit_s702_transition_state.py','validation/report_s702_wake_extent.py',
                 'validation/report_s702_drag_mechanism.py','validation/report_design_initialization.py'):
        inputs.add(ROOT/name)
    write(out/'input_manifest.json',dict(study_sha256=STUDY_HASH,files=inputs.files,
        analysis_runtime=dict(python=platform.python_version(),versions={k:importlib.metadata.version(k) for k in ('numpy','scipy','matplotlib')})))
    print(summary)
    return summary


def plot(out,cases,groups,experiment,context,effects):
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    def title(c):return f"Flap {c['flap']}° · CL ≈ {c['target_cl']:.3f}"
    def cp(ax,c,e,side='upper',labels=True):
        ex=experiment[c['name']][e,side];ex=ex[ex['endpoint_duplicate']==0]
        ax.plot(ex['x_over_cref'],ex['Cp'],'o',color='black',mfc='white',ms=4,zorder=5,label='Experiment' if labels else None)
        for k in KEYS:
            branch=pressure_branch(groups[c['name']][k]['state']['surfaces'][e],side)
            ax.plot(branch['x_over_cref'],branch['Cp'],color=COLORS[k],ls=STYLES[k],lw=1.6,label=LABELS[k] if labels else None)
        ax.grid(alpha=.2);ax.set_xlabel('x / reference chord');ax.set_ylabel('$C_p$')
    fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    for ax,c in zip(axes.flat,cases.values()):
        cp(ax,c,'aft');ax.axvline(.9001,color='.5',ls=':',lw=1);ax.axvspan(.86,.94,color='.6',alpha=.07)
        ax.set(xlim=(.80,.975),ylim=(.15,-1.05),title=title(c));ax.legend(fontsize=8)
    fig.suptitle('Separating fore and aft transition effects on aft pressure\nRe = 1,000,000 · M = 0.10 · each condition matched to measured lift')
    save(fig,out,'aft_pressure_factorial')
    fig,axes=plt.subplots(4,3,figsize=(15,13),layout='constrained')
    for trio,c in zip(axes,cases.values()):
        cp(trio[0],c,'aft');trio[0].set(xlim=(.82,.96),ylim=(.15,-1.05),title=title(c))
        for k in KEYS:
            s=groups[c['name']][k]['state'];up=side_data(s['surfaces']['aft'],'upper')
            trio[1].plot(up['x_over_cref'],up['cf'],color=COLORS[k],ls=STYLES[k],label=LABELS[k])
            bl=s['result']['boundary_layers']['aft']['upper'];tr=bl['transition']
            if tr:trio[1].plot(tr['xy_over_cref'][0],0.,'^',color=COLORS[k],ms=6)
            fore=side_data(s['surfaces']['fore'],'upper');lam=fore['turbulent']==0
            trio[2].plot(fore['x_over_cref'][lam],fore['amplification_or_shear'][lam]/k[0],color=COLORS[k],ls=STYLES[k])
        trio[1].axhline(0,color='.4',lw=.8);trio[1].set(xlim=(.82,.99),ylim=(-.002,.006),ylabel='Aft upper $C_f$',xlabel='x / reference chord')
        trio[2].axhline(1,color='.4',ls=':',lw=.8);trio[2].set(xlim=(0.,.82),ylim=(0.,1.08),ylabel='Fore upper amplification / Ncrit',xlabel='x / reference chord')
        for ax in trio:ax.grid(alpha=.15)
    axes[0,0].legend(fontsize=7);axes[0,1].legend(fontsize=7)
    fig.suptitle('Pressure, aft separation and fore amplification\nTriangles: calculated aft transition; no experimental transition endpoints available')
    save(fig,out,'boundary_layer_factorial')
    fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    for ax,c in zip(axes.flat,cases.values()):
        f=next(r for r in effects if r['case']==c['name'] and r['quantity']=='CD_counts')
        values=[f[k] for k in ('fore_effect','aft_effect','interaction','combined_change')]
        bars=ax.bar(['Fore only','Aft only','Interaction','Both'],values,color=[COLORS[KEYS[1]],COLORS[KEYS[2]],'#777777',COLORS[KEYS[3]]])
        ax.bar_label(bars,fmt='%+.2f',padding=4);ax.axhline(0,color='black',lw=.8)
        ax.set(title=title(c),ylabel='Change in total drag (counts)');ax.margins(y=.25);ax.grid(axis='y',alpha=.2)
    fig.suptitle('What drives the drag response when Ncrit changes from 9 to 6.25?\nFore + aft + interaction = both-elements change; each condition uses its own baseline')
    save(fig,out,'drag_factorial_effects')
    fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    for ax,c in zip(axes.flat,cases.values()):
        measured=float(context[c['name']]['prehinge_delta_Cp'])
        values=[groups[c['name']][k]['quantities']['prehinge_delta_Cp'] for k in KEYS]
        bars=ax.bar(['Both 9','Fore only','Aft only','Both 6.25'],values,color=[COLORS[k] for k in KEYS])
        ax.axhline(measured,color='black',ls='--',label=f'Experiment {measured:.4f}');ax.bar_label(bars,fmt='%.4f',padding=3,fontsize=9)
        ax.set_ylim(0,max(max(values),measured)*1.3);ax.set(title=title(c),ylabel='$C_p$(≈0.889c) − $C_p$(≈0.872c)');ax.legend(fontsize=9)
    fig.suptitle('Does the transition change reproduce the pre-hinge pressure feature?\nSame two measured taps; no pressure interpolation between operating conditions')
    save(fig,out,'pressure_feature_factorial')
    fig,axes=plt.subplots(2,1,figsize=(12,8),layout='constrained',gridspec_kw={'height_ratios':[1,2]})
    for flap,ls in ((5,'-'),(10,'--')):
        for e,xy in geometry(flap):
            if e=='fore' and flap==10:continue
            axes[0].plot(*xy.T,ls=ls,label=f'{e}, flap {flap}°' if e=='aft' else 'fore')
    axes[0].plot(.9001,-.0114,'k+',ms=10,label='hinge');axes[0].set_aspect('equal',adjustable='box')
    axes[0].set(xlabel='x / c',ylabel='z / c',title='Fixed reconstructed geometry; exact tested contour remains unqualified');axes[0].legend(ncol=4,fontsize=9)
    ax=axes[1];x=np.arange(len(cases));width=.15
    for i,k in enumerate(KEYS):
        ax.bar(x+(i-2)*width,[groups[c][k]['target']['error_CD_counts'] for c in cases],width,color=COLORS[k],label=LABELS[k])
    ax.bar(x+2*width,[(float(context[c]['mses_CD'])-float(context[c]['experiment_CD']))*1e4 for c in cases],width,color='.5',label='Published MSES N9')
    ax.axhline(0,color='black',lw=1);ax.set_xticks(x,[title(c).replace(' · ','\n') for c in cases.values()]);ax.set_ylabel('Calculated − experimental drag (counts)')
    ax.legend(fontsize=9,ncol=3);ax.grid(axis='y',alpha=.2)
    fig.suptitle('Total drag context at each pressure snapshot’s own lift\nPublished MSES curves only; these unequal-CL snapshots are not a flap-drag difference')
    save(fig,out,'geometry_drag_mses')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--out',type=type(OUT),default=OUT)
    report(parser.parse_args().out)
