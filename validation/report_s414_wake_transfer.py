#!/usr/bin/env python3
"""Trace and test fore-wake information channels at archived S414 states.

Only restore operators and evaluate residuals; no nonlinear solver is called.
The restoration is an explicit adaptation of the unchanged S414 state auditor.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key]='1'
import ast
import json
import platform
import importlib.metadata
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from viscous import CoupledSolver,mf,te_info
from s414_geometry import elements
from audit_s414_state import ROOT,STUDY,Inputs,load_states,sha,require,close,array_close
from report_s702_wake_extent import csv_write
from report_s702_matched_pressure import pressure_branch,save
from s702_ncrit_sweep import write

AUDITED=ROOT/'validation/results/s414_reconstruction_analysis_20261003'
OUT=ROOT/'validation/results/s414_wake_transfer_20261003'
VARIABLES=('theta','delta_star','root_shear','edge_velocity')
STEPS=(1e-3,1e-4)
CODE=['validation/report_s414_wake_transfer.py','validation/S414_WAKE_TRANSFER_PROTOCOL.md',
      'validation/report_s702_matched_pressure.py']
SOURCE_PDF=ROOT/'validation/results/s414_geometry_qualification_20261002/sources/s414_report_2010.pdf'
ARTICLE=ROOT/'validation/results/reference_review_20260924/sources/s414_article.pdf'


def restore(job, spec, state):
    raw=elements(job['scenario'],spec['panels'])
    solver=CoupledSolver(raw,reynolds=spec['reynolds'],mach=spec['mach'],ncrit=spec['ncrit'],
        iterations=spec['iterations'],tolerance=spec['residual_tolerance'],wake_length=job['wake_length'])
    result=state['result'];alpha=result['alpha_deg']
    direction=np.array([np.cos(np.deg2rad(alpha)),np.sin(np.deg2rad(alpha))])
    for model in solver.models:
        model.oper.alpha=alpha;model.isol.gam=model.isol.gamref@direction
    solver._build_wakes(alpha)
    for model in solver.models:
        a=state['arrays'][model.geom.name]
        array_close(model.foil.x,a['xy'],'Rebuilt foil geometry',2e-12)
        array_close(model.wake.x,a['wake_xy'],'Rebuilt inviscid wake geometry',2e-11)
        array_close(model.wake.t,a['wake_t'],'Rebuilt wake tangents',2e-11)
        array_close(model.wake.s,a['wake_s'],'Rebuilt wake arclength',2e-11)
        model.glob.U=a['U'].copy();model.glob.Nsys=a['U'].shape[1]
        model.vsol.turb=a['turbulent'].copy()
        model.isol.Istag=[int(a['lower'][0]),int(a['upper'][0])]
        model.isol.sgnue=np.ones(model.foil.N);model.isol.sgnue[a['lower']]=-1
        model.isol.sstag=float(model.foil.s[a['lower'][0]]+a['xi'][a['lower'][0]])
        mf.identify_surfaces(model);mf.stagpoint_move(model)
        array_close(model.isol.xi,a['xi'],'Restored stagnation distances',2e-12)
        for k,label in enumerate(('lower','upper')):
            require(np.array_equal(model.vsol.Is[k],a[label]),'Restored branch indices differ')
        if abs(te_info(model.foil.x)[1])<1e-12:model.vsol.wgap=np.zeros(model.wake.N)
        else:mf.set_wake_gap(model)
    solver._build_source_coupling()
    return solver


def nearest_polyline(points,wake):
    """Project onto open wake segments; never close the wake or extrapolate."""
    points,wake=np.asarray(points,float),np.asarray(wake,float)
    require(points.ndim==wake.ndim==2 and points.shape[1]==wake.shape[1]==2,'Expected planar points')
    segment=np.diff(wake,axis=0);length2=np.sum(segment**2,axis=1)
    require(len(segment)>0 and np.all(length2>0),'Degenerate wake segment')
    fractions=np.clip(np.sum((points[:,None,:]-wake[None,:-1,:])*segment,axis=2)/length2,0,1)
    projected=wake[None,:-1,:]+fractions[:,:,None]*segment
    squared=np.sum((points[:,None,:]-projected)**2,axis=2);index=np.argmin(squared,axis=1)
    row=np.arange(len(points));return np.sqrt(squared[row,index]),index,fractions[row,index],projected[row,index]


def segment_values(values,index,fraction):
    values=np.asarray(values,float)
    return values[index]*(1-fraction)+values[index+1]*fraction


def source_parts(solver):
    """Partition the implemented source columns; this is an algebraic identity."""
    state=np.hstack([m.glob.U for m in solver.models]);sign=np.ones(solver.nsys)
    for m,ss in zip(solver.models,solver.ss):sign[ss.start:ss.start+m.foil.N]=m.isol.sgnue
    parts={}
    for m,ss,src in zip(solver.models,solver.ss,solver.sources):
        sigma=np.asarray(m.vsol.sigma_m@(state[1,ss]*state[3,ss]))
        require(len(sigma)==m.foil.N+m.wake.N-2,'Source partition dimensions differ')
        for label,a,b in [('surface',0,m.foil.N-1),('wake',m.foil.N-1,len(sigma))]:
            cols=slice(src.start+a,src.start+b)
            parts[m.geom.name+'_'+label]=sign*(solver.ue_sigma[:,cols]@sigma[a:b])
    total=sum(parts.values());expected=solver.D@(state[1]*state[3])
    array_close(total,expected,'Source-column sum differs from mass operator',2e-12)
    return parts


def probe(solver):
    require([m.geom.name for m in solver.models]==['fore','aft'],'Unexpected model order')
    saved=[m.glob.U.copy() for m in solver.models]
    flags=[m.vsol.turb.copy() for m in solver.models]
    baseline,J=solver._system();ns=solver.nsys
    norm=max(float(np.linalg.norm(baseline[:3*ns])),float(abs(baseline[3*ns:]).max()))
    require(norm<=1e-7,'Baseline residual is not admissible')
    fore,aft=solver.models;fs,ats=solver.ss
    wake_local=np.arange(fore.foil.N,fore.glob.Nsys);wake_global=fs.start+wake_local
    receiver_bl=slice(3*ats.start,3*ats.stop);receiver_velocity=slice(3*ns+ats.start,3*ns+ats.stop)
    positive_control=slice(3*fs.start,3*fs.stop)
    rows=[];blocks=[]
    for variable,name in enumerate(VARIABLES):
        columns=4*wake_global+variable
        blmax=float(abs(J[receiver_bl][:,columns]).max())
        velocitymax=float(abs(J[receiver_velocity][:,columns]).max())
        require(blmax==0.,'Foreign state unexpectedly enters aft BL equations directly')
        if variable in (0,2):require(velocitymax==0.,'Unexpected direct momentum/shear transfer')
        else:require(velocitymax>0.,'Expected mean-flow channel is missing')
        blocks.append(dict(variable=name,aft_BL_block_max_abs=blmax,aft_velocity_block_max_abs=velocitymax))
        direction=np.zeros((4,ns));direction[variable,wake_global]=saved[0][variable,wake_local]
        expected=J@direction.ravel(order='F')
        for step in STEPS:
            plus=None;minus=None
            for sign in (1,-1):
                fore.glob.U[:]=saved[0]
                fore.glob.U[variable,wake_local]*=1+sign*step
                residual,_=solver._system()
                if sign==1:plus=residual.copy()
                else:minus=residual.copy()
            fore.glob.U[:]=saved[0]
            numerical=(plus-minus)/(2*step)
            # Test the receiving equations; own-BL response is an independent positive control.
            errors=[]
            for selection in (receiver_bl,receiver_velocity):
                scale=max(float(abs(expected[selection]).max()),1e-12)
                errors.append(float(abs(numerical[selection]-expected[selection]).max())/scale)
            require(max(errors)<1e-7,'Cross-element directional derivative check failed')
            blchange=float(abs(plus[receiver_bl]-baseline[receiver_bl]).max())
            velocitychange=float(abs(plus[receiver_velocity]-baseline[receiver_velocity]).max())
            ownchange=float(abs(plus[positive_control]-baseline[positive_control]).max())
            require(ownchange>1e-10,'Perturbation did not change the source element')
            require(blchange==0.,'Aft BL residual changed despite fixed local state')
            if variable in (0,2):require(velocitychange==0.,'Unexpected direct shear/momentum channel')
            else:require(velocitychange>0.,'Missing fore-wake mean-flow response')
            rows.append(dict(variable=name,relative_step=step,
                max_aft_BL_residual_change=blchange,max_aft_velocity_residual_change=velocitychange,
                max_fore_BL_residual_change=ownchange,
                aft_BL_Jacobian_relative_discrepancy=errors[0],aft_velocity_Jacobian_relative_discrepancy=errors[1],
                whole_system_directional_discrepancy=float(abs(numerical-expected).max()),
                diagnostic_only=True))
    final,_=solver._system();array_close(final,baseline,'Baseline residual not restored',0.)
    for m,u,turb in zip(solver.models,saved,flags):
        require(np.array_equal(m.glob.U,u) and np.array_equal(m.vsol.turb,turb),'State/branch flags mutated')
        require(m.param.ncrit==9.,'Ncrit changed')
    parts=source_parts(solver)
    inv=np.concatenate([mf.get_ueinv(m) for m in solver.models]);U=np.hstack(saved)
    array_close(U[3]-inv-sum(parts.values()),baseline[3*ns:],'Source identity differs from residual',2e-12)
    aft_upper_start=aft.foil.N-1-int(np.argmin(aft.foil.x[0,::-1]))
    locations=[]
    for local in range(aft_upper_start,aft.foil.N):
        index=ats.start+local
        locations.append(dict(x_over_c=float(aft.foil.x[0,local]),y_over_c=float(aft.foil.x[1,local]),
            edge_velocity_state=float(U[3,index]),inviscid_edge_velocity=float(inv[index]),
            **{name:float(value[index]) for name,value in parts.items()},
            summed_displacement_velocity=float(sum(value[index] for value in parts.values())),
            coupling_residual=float(baseline[3*ns+index])))
    return dict(baseline_residual=norm,blocks=blocks,finite_differences=rows,source_columns=locations,
        states_restored_exactly=True,ncrit=9.,new_flow_solves=0,
        qualification='Frozen-state residual partial derivatives and source-column decomposition; no equilibrated force or transition response.')


def geometry_context(spec,states):
    proximity=[];contexts=[]
    for job in spec['jobs']:
        if 'repeat_of' in job:continue
        s=states[job['name']];a=s['arrays']['fore'];wake=s['wakes']['fore']
        branch=pressure_branch(s['surfaces']['aft'],'upper')
        branch=branch[(branch['x_over_cref']>=.8)&(branch['x_over_cref']<=.98)]
        points=np.column_stack((branch['x_over_cref'],branch['y_over_cref']))
        distance,index,fraction,projected=nearest_polyline(points,a['wake_xy'].T)
        ds=segment_values(wake['delta_star'],index,fraction)
        theta=segment_values(wake['theta'],index,fraction)
        stress=segment_values(a['U'][2,a['xy'].shape[1]:],index,fraction)
        ue=segment_values(wake['ue_over_vinf'],index,fraction)
        for i,point in enumerate(points):
            proximity.append(dict(job=job['name'],case=job['pressure_case'],scenario=job['scenario'],
                x_over_c=float(point[0]),aft_surface_z=float(point[1]),wake_nearest_x=float(projected[i,0]),
                wake_nearest_z=float(projected[i,1]),centerline_distance_c=float(distance[i]),
                fore_wake_delta_star_c=float(ds[i]),fore_wake_theta_c=float(theta[i]),
                fore_wake_root_shear=float(stress[i]),fore_wake_ue=float(ue[i]),
                aft_delta_star_c=float(branch['delta_star'][i]),
                distance_over_displacement_sum=float(distance[i]/(ds[i]+branch['delta_star'][i]))))
        n=a['xy'].shape[1];result=s['result'];fore_bl=result['boundary_layers']['fore']
        contexts.append(dict(job=job['name'],case=job['pressure_case'],scenario=job['scenario'],
            CL=result['CL'],alpha=result['alpha_deg'],
            fore_upper_transition_x=None if fore_bl['upper']['transition'] is None else fore_bl['upper']['transition']['xy_over_cref'][0],
            fore_lower_transition_x=None if fore_bl['lower']['transition'] is None else fore_bl['lower']['transition']['xy_over_cref'][0],
            fore_upper_TE_turbulent=bool(a['turbulent'][a['upper'][-1]]),
            fore_lower_TE_turbulent=bool(a['turbulent'][a['lower'][-1]]),
            all_fore_wake_nodes_flagged_turbulent=bool(a['turbulent'][n:].all()),
            fore_wake_initial_root_shear=float(a['U'][2,n]),
            fore_wake_initial_delta_star=float(a['U'][1,n]),
            fore_wake_initial_theta=float(a['U'][0,n])))
    return proximity,contexts


def source_trace(inputs):
    entries=[]
    selectors={
        'viscous.py':{'new_model':'Independent per-element Ncrit and integral state',
            'dense_boundary_system':'Local element station/transition/wake equations',
            'CoupledSolver._build_wakes':'Inviscid streamline geometry, separate nonintersecting wakes',
            'CoupledSolver._build_source_coupling':'Cross-element displacement source influence',
            'CoupledSolver._mass_operator':'Mass defect to source-strength map',
            'CoupledSolver._system':'Block-local BL rows; global ue - ue_inv - D(ue delta*) rows'},
        'vendor/mfoil.py':{'wake_sys':'TE-layer sums; laminar TE shear converted with get_cttr',
            'build_param':'Wake closure always turbulent',
            'residual_station':'Local amplification or turbulent shear-lag equation',
            'residual_transition':'Local amplification reaches fixed Ncrit',
            'get_damp':'Local Hk, Re_theta, theta and amplification determine growth',
            'get_cttr':'Root shear initialization, not turbulence intensity'}}
    for file,selected in selectors.items():
        path=inputs.add(ROOT/file);tree=ast.parse(path.read_text());nodes={}
        for node in tree.body:
            if isinstance(node,ast.FunctionDef):nodes[node.name]=node
            elif isinstance(node,ast.ClassDef):
                for method in node.body:
                    if isinstance(method,ast.FunctionDef):nodes[node.name+'.'+method.name]=method
        for name,role in selected.items():
            node=nodes[name];entries.append(dict(file=file,function=name,line=node.lineno,end_line=node.end_lineno,role=role,sha256=sha(path)))
    return entries


def main():
    OUT.mkdir(exist_ok=True);inputs=Inputs();spec,states,_=load_states(inputs)
    verification=inputs.json(AUDITED/'verification.json')
    for name in ('residual_recheck.json','summary.json','diagnosis.json','FINDINGS.md'):
        require(sha(inputs.add(AUDITED/name))==verification['output_sha256'][name],'Prior completed audit changed: '+name)
    residual=inputs.json(AUDITED/'residual_recheck.json')
    require(residual['complete'] and len(residual['states'])==9,'Previous audit incomplete')
    for name,digest in residual['input_sha256'].items():
        require(sha(inputs.add(ROOT/name))==digest,'Prior residual evidence invalidated: '+name)
    for name in CODE:inputs.add(ROOT/name)
    inputs.add(SOURCE_PDF);require(sha(SOURCE_PDF)=='52756797afa2f7089be9e7d2ac807fff024b77769751b84b81ad7d6312fc0c5c','Report changed')
    inputs.add(ARTICLE);require(sha(ARTICLE)=='82b7e6292278e459dc11bfc9d9d64259c1129bc4d99e0f4f4f9870af0d91763b','Article changed')
    trace=source_trace(inputs)
    environment=dict(python=platform.python_version(),versions={k:importlib.metadata.version(k) for k in ('numpy','scipy','matplotlib')})
    receipt=dict(input_sha256=dict(inputs.files),runtime=environment,cases=[],complete=False,new_flow_solves=0)
    path=OUT/'operator_checks.json'
    if path.exists():
        prior=json.loads(path.read_text());require(prior['input_sha256']==receipt['input_sha256'] and prior['runtime']==environment,
            'Diagnostic inputs/runtime changed; preserve previous output before restarting')
        receipt=prior
    for job in spec['jobs']:
        if job['scenario']!='nominal' or 'repeat_of' in job:continue
        if any(row['job']==job['name'] for row in receipt['cases']):continue
        solver=restore(job,spec,states[job['name']]);row=dict(job=job['name'],case=job['pressure_case'],**probe(solver))
        receipt['cases'].append(row);write(path,receipt)
        print(job['name']+': all direct-transfer and finite-difference checks passed',flush=True)
    require(len(receipt['cases'])==2,'Both nominal operator checks are required')
    receipt['complete']=True;write(path,receipt)
    geometry,contexts=geometry_context(spec,states)
    csv_write(OUT/'wake_proximity.csv',geometry);csv_write(OUT/'fore_layer_and_wake_context.csv',contexts)
    csv_write(OUT/'finite_difference_checks.csv',[dict(job=c['job'],**r) for c in receipt['cases'] for r in c['finite_differences']])
    csv_write(OUT/'source_column_contributions.csv',[dict(job=c['job'],**r) for c in receipt['cases'] for r in c['source_columns']])
    write(OUT/'source_trace.json',dict(functions=trace,primary_report=str(SOURCE_PDF.relative_to(ROOT)),article=str(ARTICLE.relative_to(ROOT)),
        notes=['2010 report printed pp.10–11: laminar fore surfaces in bucket; qualitative oil-flow shortening near lower limit.',
               '2014 article section VI: horn interpretation; full-e^N MSES formulation differs from the 2010 curves.',
               'No new MSES run or proof of the original MSES setup internals.']))
    plot(states,geometry,contexts,receipt)
    summary=dict(complete=True,operator_states_checked=len(receipt['cases']),geometry_states_checked=len(contexts),
        central_difference_checks=sum(len(c['finite_differences']) for c in receipt['cases']),
        largest_receiving_Jacobian_relative_discrepancy=max(r['aft_velocity_Jacobian_relative_discrepancy'] for c in receipt['cases'] for r in c['finite_differences']),
        mean_flow_displacement_channel_present=True,direct_wake_shear_to_aft_transition_channel_present=False,
        wake_closure_forced_turbulent=True,all_states_restored_exactly=all(c['states_restored_exactly'] for c in receipt['cases']),
        new_flow_solves=0,ncrit=9.,missing_channel_proven_cause=False,
        qualification='Discrete information-flow audit at fixed states; no new equilibrium, closure validation or turbulence-intensity inference.')
    write(OUT/'summary.json',summary);write(OUT/'input_manifest.json',dict(files=dict(inputs.files),runtime=environment))
    print(json.dumps(summary,indent=2))


def plot(states,geometry,contexts,receipt):
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    for column,case in enumerate(('cl028','cl104')):
        name=case+'_nominal';s=states[name];ax=axes[0,column]
        for element in ('fore','aft'):
            xy=s['arrays'][element]['xy'];ax.plot(xy[0],xy[1],color='.25',label=element+' contour')
        wake=s['arrays']['fore']['wake_xy'];ax.plot(wake[0],wake[1],color='#2166ac',label='Computed fore wake line')
        branch=pressure_branch(s['surfaces']['aft'],'upper');layer=s['result']['boundary_layers']['aft']['upper']
        for region in layer['separation_regions']:
            a,b=region['start']['xy_over_cref'][0],region['end']['xy_over_cref'][0]
            use=(branch['x_over_cref']>=a)&(branch['x_over_cref']<=b)
            ax.plot(branch['x_over_cref'][use],branch['y_over_cref'][use],color='#d6604d',lw=3,label='Calculated separated surface')
        ax.axvspan(.88215,.89313,color='.5',alpha=.12)
        ax.set(xlim=(.78,1.01),ylim=(-.04,.075),xlabel='x/c',ylabel='z/c',title=f'CL {s["result"]["CL"]:.3f}')
        ax.set_aspect('equal',adjustable='box');ax.grid(alpha=.2)
        rows=[r for r in geometry if r['job']==name];ax=axes[1,column];x=[r['x_over_c'] for r in rows]
        for field,label,color in [('centerline_distance_c','Surface to wake-line distance','.2'),
            ('fore_wake_delta_star_c','Fore wake displacement thickness','#2166ac'),
            ('aft_delta_star_c','Aft displacement thickness','#d6604d')]:
            ax.plot(x,[r[field] for r in rows],color=color,label=label)
        ax.set(xlabel='Aft surface x/c',ylabel='Length / reference chord',xlim=(.8,.98));ax.grid(alpha=.2)
    axes[0,0].legend(fontsize=8);axes[1,0].legend(fontsize=8)
    fig.suptitle('S414 fore wake and aft recovery: geometry and integral scales\nGrey band: strongest measured sampled recovery; displacement thickness is not wake width')
    save(fig,OUT,'wake_geometry_and_scales')
    fig,axes=plt.subplots(1,2,figsize=(12,4.6),layout='constrained')
    colors={'fore_surface':'#9970ab','fore_wake':'#2166ac','aft_surface':'#d6604d','aft_wake':'#e6ab02'}
    for ax,c in zip(axes,receipt['cases']):
        rows=c['source_columns'];x=[r['x_over_c'] for r in rows]
        for key,color in colors.items():ax.plot(x,[r[key] for r in rows],color=color,label=key.replace('_',' '))
        ax.plot(x,[r['summed_displacement_velocity'] for r in rows],'k--',lw=1,label='Sum')
        visible=[r[key] for r in rows if .80<=r['x_over_c']<=.96 for key in (*colors,'summed_displacement_velocity')]
        low,high=min(visible),max(visible);padding=.08*max(high-low,1e-3)
        ax.axhline(0,color='.5',lw=.8);ax.set(xlim=(.80,.96),ylim=(low-padding,high+padding),xlabel='Aft upper x/c',ylabel='Contribution to uncorrected edge velocity / V∞',title=c['case']+' nominal');ax.grid(alpha=.2)
    axes[0].legend(fontsize=8);fig.suptitle('The fore wake participates in mean-flow coupling\nDiscrete source-column contributions at the archived solution; not removal experiments')
    save(fig,OUT,'mean_flow_source_channels')


if __name__=='__main__':main()
