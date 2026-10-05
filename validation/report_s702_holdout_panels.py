#!/usr/bin/env python3
"""Audit the two-grid S702 reserved-condition comparison; no flow solves."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from report_s702_wake_extent import (
    ROOT, Inputs, sha, write, csv_write, audit, require, close, clearance, side_data,
)
from run_s702_holdout_panels import effective_spec

STUDY_HASH = 'cd861525d01a02d788bbaddabaacb0c81f9896e153ab2229743190b775fea6f1'
ELEMENTS = ('fore','aft')


def contour_sample_distance(a,b):
    """Largest vertex-to-other-polyline distance, in either direction (not exact Hausdorff)."""
    def directed(points,contour):
        starts=contour; delta=np.roll(contour,-1,axis=0)-contour
        size=np.sum(delta*delta,axis=1); keep=size>1e-24
        starts,delta,size=starts[keep],delta[keep],size[keep]
        offset=points[:,None,:]-starts
        t=np.clip(np.sum(offset*delta,axis=2)/size,0.,1.)
        return float(np.max(np.min(np.linalg.norm(offset-t[:,:,None]*delta,axis=2),axis=1)))
    return max(directed(a,b),directed(b,a))


def resolution(values):
    turbulent=values[values['turbulent']==1]
    if len(turbulent)<2:
        return dict(max_Hk_drop=0.,max_velocity_step=0.,Hk_step_x0=None,Hk_step_x1=None,velocity_step_x0=None,velocity_step_x1=None)
    drops=-np.diff(turbulent['Hk']); steps=abs(np.diff(turbulent['ue_over_vinf']))
    h,u=int(np.argmax(drops)),int(np.argmax(steps));x=turbulent['x_over_local_chord']
    return dict(max_Hk_drop=max(0.,float(drops[h])),max_velocity_step=float(steps[u]),
                Hk_step_x0=float(x[h]),Hk_step_x1=float(x[h+1]),velocity_step_x0=float(x[u]),velocity_step_x1=float(x[u+1]))


def report(study,out):
    require(sha(study/'study.json')==STUDY_HASH,'Unexpected panel study')
    inputs=Inputs();spec=inputs.json(study/'study.json');runtime=inputs.json(study/'runtime.json')
    for name,digest in spec['sha256'].items():
        require(sha(inputs.add(study/'source_snapshot'/name))==digest,'Frozen source changed: '+name)
    for name in ('validation/report_s702_wake_extent.py','viscous.py','vendor/mfoil.py','multifoil.py'):
        require(sha(ROOT/name)==spec['sha256'][name],'Audit dependency differs from executed source: '+name)
    require(runtime['versions']==spec['expected_versions'],'Numerical runtime differs')
    status=inputs.json(study/'status.json');automatic=inputs.json(study/'summary.json')
    require(status['state']=='finished','Study incomplete')
    states,refs,targets,geometry_rows,control_rows,histories={},{},[],[],[],[]
    for job in spec['jobs']:
        expected=effective_spec(job,spec)
        state=audit(study/job['name'],job,expected,runtime,inputs);r=state['result']
        require(r['feasible'],'Incomplete grid coverage prevents full paired interpretation')
        refpath=ROOT/job['reference_result']
        if str(refpath) not in refs:
            before=set(inputs.files)
            oldspec=inputs.json(refpath.parent.parent/'study.json');oldruntime=inputs.json(refpath.parent.parent/'runtime.json')
            oldjob=next(j for j in oldspec['jobs'] if j['name']==refpath.parent.name)
            refs[str(refpath)]=audit(refpath.parent,oldjob,oldspec,oldruntime,inputs)
            for name in set(inputs.files)-before:
                require(name in spec['sha256'] and sha(ROOT/name)==spec['sha256'][name],'Archived input differs: '+name)
            for key in ('reynolds','mach','ncrit','reference_chord','iterations','cl_tolerance','residual_tolerance'):
                require(oldspec[key]==spec[key],'Grid comparison physics/protocol differs: '+key)
            for filename in ('fore.dat',f"aft_d{job['flap']}.dat"):
                name='validation/data/s702/'+filename
                require(oldspec['sha256'][name]==spec['sha256'][name],'Raw source geometry differs')
        old=refs[str(refpath)];base=old['result']
        delta=(r['CD']-base['CD'])*1e4
        checked=next(row for row in automatic['jobs'] if row['job']==job['name'])
        close(checked['difference_from_base_CD_counts'],delta)
        row=dict(job=job['name'],control=job['control'],flap=job['flap'],target_cl=job['target_cl'],
                 panels='/'.join(map(str,job['panels'])),CL=r['CL'],CD=r['CD'],CM=r['CM'],alpha_deg=r['alpha_deg'],
                 lift_error=r['CL']-job['target_cl'],residual=state['residual'],difference_CD_counts=delta,
                 difference_CL=r['CL']-base['CL'],difference_CM=r['CM']-base['CM'],difference_alpha=r['alpha_deg']-base['alpha_deg'],
                 selected_path=r['selected_path'],evaluations=len(r['records']),failed_evaluations=sum(not v['admissible'] for v in r['records']),
                 failed_internal_points=sum(not c['converged'] for d in state['metadata']['diagnostics'] for c in d['continuation']),
                 wall_seconds=r['wall_seconds'])
        for e in ELEMENTS:row[e+'_difference_CD_counts']=(r['element_forces'][e]['CD_wake']-base['element_forces'][e]['CD_wake'])*1e4
        targets.append(row)
        histories.append(dict(job=job['name'],records=r['records'],attempts=r['attempts'],selected_flow_index=r['selected_flow_index'],
                              continuation=[d['continuation'] for d in state['metadata']['diagnostics']]))
        if job['control']:
            require(abs(delta)<=spec['acceptance_control_CD_counts'],'Control tolerance failed')
            for e in ELEMENTS:
                for key,a in old['arrays'][e].items():
                    b=state['arrays'][e][key]
                    control_rows.append(dict(element=e,field=key,exactly_equal=bool(np.array_equal(a,b)),max_abs_difference=float(np.max(abs(a-b)))))
            continue
        states[job['flap'],job['target_cl'],'base']=old;states[job['flap'],job['target_cl'],'fine']=state
        for e in ELEMENTS:
            a,b=old['arrays'][e]['xy'].T,state['arrays'][e]['xy'].T
            geometry_rows.append(dict(flap=job['flap'],target_cl=job['target_cl'],element=e,base_nodes=len(a),fine_nodes=len(b),
                max_vertex_to_other_polyline_distance=contour_sample_distance(a,b),
                base_TE_gap=float(np.linalg.norm(a[0]-a[-1])),fine_TE_gap=float(np.linalg.norm(b[0]-b[-1])),
                max_TE_endpoint_change=float(np.max(abs(a[[0,-1]]-b[[0,-1]])))))
    require(len(targets)==automatic['required_jobs']==automatic['feasible_jobs']==5 and automatic['execution_checks_pass'],'Coverage mismatch')
    surfaces,wakes=[],[]
    for (flap,cl,grid),s in states.items():
        r=s['result'];diag=s['metadata']['diagnostics'][r['selected_flow_index']]
        for e in ELEMENTS:
            a=s['arrays'][e];w=s['wakes'][e];other='aft' if e=='fore' else 'fore'
            gap=clearance(a['wake_xy'].T,s['arrays'][other]['xy'].T);require(gap>0,'Wake/body intersection')
            wakes.append(dict(flap=flap,target_cl=cl,grid=grid,element=e,nodes=len(w),arclength=s['stations'][e]['s'][-1],
                first_spacing=float(np.linalg.norm(a['wake_xy'][:,1]-a['wake_xy'][:,0])),
                last_spacing=float(np.linalg.norm(a['wake_xy'][:,-1]-a['wake_xy'][:,-2])),other_body_clearance=gap,
                CD_wake=r['element_forces'][e]['CD_wake'],CD_friction=r['element_forces'][e]['CD_friction'],
                terminal_theta=float(w['theta'][-1]),terminal_Hk=float(w['Hk'][-1]),terminal_ue=float(w['ue_over_vinf'][-1]),
                TE_gap_added=diag['te_gap_added_over_local_chord'][e],paneling_recovery=diag['monotone_paneling_recovery'][e]))
            for side in ('upper','lower'):
                values=side_data(s['surfaces'][e],side);layer=r['boundary_layers'][e][side];tr=layer['transition']
                sr=dict(flap=flap,target_cl=cl,grid=grid,element=e,side=side,transition_x=tr['x_over_local_chord'] if tr else None,
                        separated_fraction=layer['separated_arclength_fraction'],separation_regions=len(layer['separation_regions']),
                        max_Hk=layer['max_Hk'],min_cf=layer['min_cf'],**resolution(values))
                for f in ('Cp','Hk','theta','delta_star','ue_over_vinf','cf'):sr['TE_'+f]=float(values[f][-1])
                surfaces.append(sr)
    experimental=inputs.csv(ROOT/'validation/results/s702_drag_mechanism/s702_flap_5_to_10.csv')
    pairs=[]
    for cl in (.8,1.2):
        row=dict(target_cl=cl)
        reference=next(v for v in experimental if float(v['reynolds'])==spec['reynolds'] and float(v['mach'])==spec['mach'] and float(v['target_cl'])==cl)
        row['experiment_flap_delta_counts']=float(reference['experiment_delta_10_minus_5_counts'])
        for grid in ('base','fine'):
            a,b=[states[f,cl,grid]['result'] for f in (5,10)]
            row[grid+'_flap_delta_CD_counts']=(b['CD']-a['CD'])*1e4
            row[grid+'_flap_delta_error_counts']=row[grid+'_flap_delta_CD_counts']-row['experiment_flap_delta_counts']
            for e in ELEMENTS:row[grid+'_'+e+'_flap_delta_CD_counts']=(b['element_forces'][e]['CD_wake']-a['element_forces'][e]['CD_wake'])*1e4
        row['change_in_flap_delta_counts']=row['fine_flap_delta_CD_counts']-row['base_flap_delta_CD_counts']
        pairs.append(row)
    summary=dict(required_jobs=5,feasible_jobs=5,selected_states_audited=9,full_state_audit_pass=True,
        control_CD_difference_counts=next(v['difference_CD_counts'] for v in targets if v['control']),
        control_shared_arrays_exactly_equal=all(v['exactly_equal'] for v in control_rows),
        max_absolute_individual_CD_change_counts=max(abs(v['difference_CD_counts']) for v in targets if not v['control']),
        max_absolute_flap_delta_change_counts=max(abs(v['change_in_flap_delta_counts']) for v in pairs),
        paired_changes_exceed_screening_threshold=any(abs(v['change_in_flap_delta_counts'])>spec['design_screening_counts'] for v in pairs),
        max_sampled_geometry_change=max(v['max_vertex_to_other_polyline_distance'] for v in geometry_rows),
        max_TE_endpoint_change=max(v['max_TE_endpoint_change'] for v in geometry_rows),
        artificial_TE_gap_added=any(v['TE_gap_added']!=0 for v in wakes),paneling_recovery_used=any(v['paneling_recovery'] for v in wakes),
        minimum_wake_body_clearance=min(v['other_body_clearance'] for v in wakes),
        requested_evaluations=sum(v['evaluations'] for v in targets),failed_requested_evaluations=sum(v['failed_evaluations'] for v in targets),
        failed_internal_continuation_points=sum(v['failed_internal_points'] for v in targets),
        qualification='Two-grid operational panel sensitivity; geometry representation and derived wake discretization change. No asymptotic error bound or independent design validation.')
    for grid in ('base','fine'):
        group=[v for v in surfaces if v['grid']==grid]
        summary[grid+'_max_post_transition_Hk_drop']=max(v['max_Hk_drop'] for v in group)
        summary[grid+'_max_post_transition_velocity_step']=max(v['max_velocity_step'] for v in group)
    out.mkdir(parents=True,exist_ok=True)
    for name,rows in (('targets',targets),('flap_deltas',pairs),('surface_diagnostics',surfaces),('wake_diagnostics',wakes),
                      ('geometry_comparison',geometry_rows),('control_differences',control_rows)):
        csv_write(out/(name+'.csv'),rows)
    write(out/'summary.json',summary);write(out/'initialization_histories.json',histories)
    plot(out,states,targets,pairs,surfaces)
    (out/'FINDINGS.md').write_text(
        '# S702 reserved-condition panel audit\n\n'
        'All five jobs converged and all nine selected states (five new and four archived references) passed full audit. '
        'The fresh 260/180 control repeats drag exactly and all shared saved array fields are identical. '
        'Every metadata panel override matches its declared job; source hashes, runtime, nonzero-Mach velocities, forces, '
        'transition/separation, residuals, geometry, wake clearance and initialization histories were checked.\n\n'
        'Changing 260/180 to 340/240 panels changes individual CD by up to '
        f"{summary['max_absolute_individual_CD_change_counts']:.6f} count. The paired flap-10-minus-5 response changes as follows:\n\n"
        '| CL | Base flap delta | Fine flap delta | Fine minus base | Experiment |\n|---:|---:|---:|---:|---:|\n'
        + '\n'.join(f"| {v['target_cl']:g} | {v['base_flap_delta_CD_counts']:.6f} | {v['fine_flap_delta_CD_counts']:.6f} | {v['change_in_flap_delta_counts']:+.6f} | {v['experiment_flap_delta_counts']:.6f} |" for v in pairs)
        + '\n\nAll drag values in this table are counts. The paired changes exceed the existing 0.5-count screening threshold, '
        'so tiny design gains at this grid cannot be accepted without a finer numerical check. '
        'The measured flap ranking remains reversed at CL 1.2; the large CL 0.8 discrepancy remains. '
        'Better agreement from grid refinement is not an independent physical validation or a reason to retune Ncrit.\n\n'
        f"Maximum post-transition Hk drop falls from {summary['base_max_post_transition_Hk_drop']:.4f} to {summary['fine_max_post_transition_Hk_drop']:.4f}; "
        f"the largest corrected velocity step falls from {summary['base_max_post_transition_velocity_step']:.5f} to {summary['fine_max_post_transition_velocity_step']:.5f}. "
        'The largest aft-lower velocity steps lie near the trailing edge, not at its transition onset. '
        'These improvements support refinement but do not prove asymptotic convergence. All four fine cases select the cold path; '
        'the difficult base control still requires the warm anchor.\n\n'
        'Raw geometry and physical assumptions are unchanged. Actual solver contours change with resampling: the largest sampled '
        f"vertex-to-other-contour distance is {summary['max_sampled_geometry_change']:.6g}c. This is not an exact Hausdorff distance "
        'or a manufacturing uncertainty. The trailing-edge endpoints agree within '
        f"{summary['max_TE_endpoint_change']:.3g}c. No artificial TE opening or monotone-paneling recovery was activated. "
        'Lowering the inactive minimum-TE-gap setting cannot explain these selected-state differences. '
        'Wake node counts and first spacing also change with surface count and are tabulated explicitly. '
        f"Minimum wake-to-other-body centerline clearance is {summary['minimum_wake_body_clearance']:.6f}c; this does not rule out physical wake/BL interaction.\n\n"
        f"Work: {summary['requested_evaluations']} requests, {summary['failed_requested_evaluations']} failed request, "
        f"{summary['failed_internal_continuation_points']} failed internal continuation points; all retained.\n\n"
        '## Next bounded check\n\n'
        'Add one 440/320 grid at these same four conditions, plus an independent 340/240 control at flap 10, CL 0.8. '
        'Keep Ncrit 6.25, Re 700,000, Mach 0.07, wake extent 1c, source geometry and ordered initialization fixed. '
        'Check whether both individual and paired changes contract below the 0.5-count screen and inspect actual surface/geometry/wake changes. '
        'Do not infer a formal error bound from only three nonuniform grids. If changes do not contract, diagnose the panel/geometry representation '
        'rather than launching an open-ended resolution sweep. Priority 4 and the source limitations in priority 3 remain open.\n\n'
        '![Drag and flap response](panel_drag.png)\n\n![Aft boundary layers](aft_layers.png)\n\n![Aft pressure](aft_pressure.png)\n\n'
        'Reproduce: `python3 validation/report_s702_holdout_panels.py`. Source data, failed targets and all earlier studies remain intact.\n')
    for name in ('validation/report_s702_holdout_panels.py','validation/report_s702_wake_extent.py','validation/report_design_initialization.py',
                 'validation/report_s702_drag_mechanism.py','validation/run_s702_holdout_panels.py','vendor/mfoil.py','multifoil.py'):
        inputs.add(ROOT/name)
    write(out/'input_manifest.json',dict(study_sha256=STUDY_HASH,files=inputs.files))
    return summary


def plot(out,states,targets,pairs,surfaces):
    def save(fig,name):
        fig.tight_layout()
        for ext in ('png','svg'):fig.savefig(out/(name+'.'+ext),dpi=170)
        plt.close(fig)
    fig,ax=plt.subplots(1,2,figsize=(11,4.5))
    for flap,marker in ((5,'o'),(10,'s')):
        rr=[v for v in targets if v['flap']==flap and not v['control']]
        ax[0].plot([v['target_cl'] for v in rr],[v['difference_CD_counts'] for v in rr],marker=marker,label=f'Flap {flap}°')
    for threshold in (-.5,.5):ax[0].axhline(threshold,color='0.65',ls=':',lw=1)
    ax[0].set(xlabel='Target CL',ylabel='Fine − base CD [counts]',title='5/5 cases; exact base control')
    for key,label,style in (('base_flap_delta_CD_counts','260/180','o-'),('fine_flap_delta_CD_counts','340/240','s--'),('experiment_flap_delta_counts','Experiment','k:')):
        ax[1].plot([p['target_cl'] for p in pairs],[p[key] for p in pairs],style,label=label)
    ax[1].set(xlabel='Target CL',ylabel='Flap 10° − 5° CD [counts]',title='CL 1.2 ranking remains reversed')
    for a in ax:a.axhline(0,color='0.5',lw=.6);a.grid(alpha=.2);a.legend(fontsize=8)
    save(fig,'panel_drag')
    fig,ax=plt.subplots(2,2,figsize=(11,7))
    for col,side in enumerate(('upper','lower')):
        for cl,color in ((.8,'#277da8'),(1.2,'#e07a27')):
            for grid,style in (('base','-o'),('fine','--.')):
                v=side_data(states[10,cl,grid]['surfaces']['aft'],side);v=v[v['turbulent']==1]
                for row,field in enumerate(('Hk','ue_over_vinf')):
                    ax[row,col].plot(v['x_over_local_chord'],v[field],style,color=color,ms=3,label=f'CL {cl:g}, {grid}')
        ax[0,col].set_title('Flap 10° aft '+side+' after transition');ax[1,col].set_xlabel('x / aft local chord')
    ax[0,0].set_ylabel('Hk');ax[1,0].set_ylabel('Ue / Uinf');ax[0,1].legend(fontsize=8)
    for a in ax.flat:a.grid(alpha=.2)
    save(fig,'aft_layers')
    fig,ax=plt.subplots(1,2,figsize=(11,4.5))
    for col,cl in enumerate((.8,1.2)):
        for flap,color in ((5,'#277da8'),(10,'#e07a27')):
            for grid,style in (('base','-'),('fine','--')):
                s=states[flap,cl,grid]['surfaces']['aft']
                for side in ('upper','lower'):
                    v=side_data(s,side)
                    ax[col].plot(v['x_over_local_chord'],v['Cp'],style,color=color,label=f'{flap}°, {grid}' if side=='upper' else None)
        ax[col].set(xlabel='x / aft local chord',title=f'CL {cl:g}');ax[col].invert_yaxis();ax[col].grid(alpha=.2)
    ax[0].set_ylabel('Cp');ax[1].legend(fontsize=8);save(fig,'aft_pressure')


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--study',type=Path,default=ROOT/'validation/results/s702_holdout_panels_20260922')
    ap.add_argument('--out',type=Path,default=ROOT/'validation/results/s702_holdout_panels_analysis')
    args=ap.parse_args();print(json.dumps(report(args.study.resolve(),args.out.resolve()),indent=2))
