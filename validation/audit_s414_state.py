#!/usr/bin/env python3
"""Audit S414 saved states and freshly evaluate residuals; no flow solves.

Reuses the unchanged full-export auditor. The residual recheck is adapted from
 audit_s702_saved_pressure_state.py only to reconstruct each declared S414 input
 and to cross-check all separation intervals. Original receipts remain intact.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key]='1'
import sys
from pathlib import Path
import platform
import importlib.metadata
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from viscous import CoupledSolver,boundary_layer_summary,mf,te_info
from s414_geometry import elements
from run_s414_reconstruction import check_controls
from report_s702_wake_extent import ROOT,Inputs,sha,audit,require,close,array_close
from s702_ncrit_sweep import write

STUDY=ROOT/'validation/results/s414_reconstruction_20261002'
OUT=ROOT/'validation/results/s414_reconstruction_analysis_20261003'
LAUNCH=ROOT/'validation/results/s414_reconstruction_launch_20261002'
STUDY_HASH='56e91de0f2c19c16c2d732c084e3141cd19c765d38bcd29db3d074ae5685d9f9'
AUDIT_CODE=['validation/audit_s414_state.py','validation/report_s702_wake_extent.py',
 'validation/report_s702_drag_mechanism.py','validation/report_design_initialization.py',
 'validation/plot_s702_flap_effect.py','validation/report_s702_ncrit_sweep.py']


def load_states(inputs):
    require(sha(inputs.add(STUDY/'study.json'))==STUDY_HASH,'Unexpected S414 study')
    spec=inputs.json(STUDY/'study.json')
    for name,digest in spec['sha256'].items():
        require(sha(inputs.add(ROOT/name))==digest,'Frozen input changed: '+name)
        require(sha(inputs.add(STUDY/'source_snapshot'/name))==digest,'Snapshot changed: '+name)
    receipt=inputs.json(LAUNCH/'retrieved_files_sha256.json')
    for name,digest in receipt.items():
        require(sha(inputs.add(STUDY/name))==digest,'Retrieved output changed: '+name)
    runtime=inputs.json(STUDY/'runtime.json');status=inputs.json(STUDY/'status.json')
    automatic=inputs.json(STUDY/'summary.json')
    require(status['state']=='finished' and status['completed_jobs']==status['total_jobs']==9,'Incomplete execution')
    require(automatic['converged_jobs']==automatic['required_jobs']==9 and automatic['execution_checks_pass'],'Incomplete coverage')
    require(runtime['python']=='3.12.13' and runtime['versions']==spec['expected_versions'],'Solver environment changed')
    require(set(runtime['thread_limits'].values())=={'1'},'Thread limits changed')
    controls=check_controls(STUDY,spec)
    require(controls['passed'] and controls==automatic['controls']==inputs.json(STUDY/'controls.json'),'Repeat gate differs')
    states={j['name']:audit(STUDY/j['name'],j,spec,runtime,inputs) for j in spec['jobs']}
    require(all(s['result']['feasible'] for s in states.values()),'Incomplete pressure comparison')
    for name in AUDIT_CODE:inputs.add(ROOT/name)
    return spec,states,controls


def recompute(job, spec, state):
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
    residual,_=solver._system()
    blnorm=float(np.linalg.norm(residual[:3*solver.nsys]))
    coupling=float(np.max(abs(residual[3*solver.nsys:])))
    require(max(blnorm,coupling)<=spec['residual_tolerance'],'Restored Newton residual fails tolerance')
    rows=[]
    for model in solver.models:
        e=model.geom.name;a=state['arrays'][e];surface=state['surfaces'][e]
        array_close(model.glob.U,a['U'],'Audit must not change saved solution',0.)
        mf.calc_force(model);mf.get_distributions(model)
        if model.oper.Ma>0:
            for side in range(3):
                param=mf.build_param(model,side)
                for j in model.vsol.Is[side]:
                    rho,_=mf.get_rho(model.glob.U[:,j],param)
                    model.post.cf[j]*=rho/model.oper.rho
        array_close(model.post.cp[:model.foil.N][::-1],surface['Cp'],'Recomputed Cp',2e-12)
        array_close(model.post.cf[:model.foil.N][::-1],surface['cf'],'Recomputed Cf',2e-12)
        forces=result['element_forces'][e]
        for name,field in (('CL','cl'),('CM','cm'),('CD_wake','cd'),('CD_friction','cdf')):
            close(getattr(model.post,field),forces[name],context='Recomputed '+e+' '+name)
        layers=boundary_layer_summary(model)
        for side in ('upper','lower'):
            expected=result['boundary_layers'][e][side]
            tr=layers[side]['transition'];old=expected['transition']
            require((tr is None)==(old is None),'Restored transition presence differs')
            if tr:close(tr['s_over_cref'],old['s_over_cref'],context='Recomputed transition')
            for key in ('min_cf','max_Hk','separated_arclength_fraction'):
                close(layers[side][key],expected[key],context='Recomputed BL '+key)
            actual_regions=layers[side]['separation_regions'];expected_regions=expected['separation_regions']
            require(len(actual_regions)==len(expected_regions),'Separation interval count differs')
            for actual_region,expected_region in zip(actual_regions,expected_regions):
                require(actual_region['reaches_trailing_edge']==expected_region['reaches_trailing_edge'],'TE separation flag differs')
                for endpoint in ('start','end'):
                    close(actual_region[endpoint]['s_over_cref'],expected_region[endpoint]['s_over_cref'],context='Recomputed separation')
        rows.append(dict(element=e,recomputed_CL=float(model.post.cl),recomputed_CD=float(model.post.cd),
            recomputed_CM=float(model.post.cm),recomputed_friction_CD=float(model.post.cdf),
            maximum_Cp_difference=float(np.max(abs(model.post.cp[:model.foil.N][::-1]-surface['Cp'])))))
    return dict(job=job['name'],BL_residual_l2=blnorm,coupling_residual_max=coupling,
        recorded_residual=state['residual'],max_residual=max(blnorm,coupling),
        new_flow_solves=0,elements=rows)


def main():
    OUT.mkdir(exist_ok=True);inputs=Inputs();spec,states,controls=load_states(inputs)
    environment=dict(python=platform.python_version(),
        versions={k:importlib.metadata.version(k) for k in ('numpy','scipy','matplotlib')})
    receipt=dict(input_sha256=inputs.files,audit_runtime=environment,states=[],complete=False,new_flow_solves=0,
        method='Nine full export/state audits and fresh coupled residual/force/pressure/BL evaluations; no nonlinear solve',controls=controls)
    path=OUT/'residual_recheck.json'
    if path.exists():
        import json
        old=json.loads(path.read_text())
        require(old['input_sha256']==receipt['input_sha256'] and old['audit_runtime']==environment,
            'Audit inputs/runtime changed; preserve previous receipt')
        receipt=old
    for job in spec['jobs']:
        if any(r['job']==job['name'] for r in receipt['states']):continue
        row=recompute(job,spec,states[job['name']]);receipt['states'].append(row);write(path,receipt)
        print(f"{len(receipt['states'])}/9 {job['name']}: residual {row['max_residual']:.3g}",flush=True)
    receipt['complete']=True;write(path,receipt)


if __name__=='__main__':main()
