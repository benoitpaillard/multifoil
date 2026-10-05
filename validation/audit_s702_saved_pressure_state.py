"""Re-evaluate saved Newton states without marching or solving a flow.

Rebuild geometry, inviscid/wake influences and residual operators, restore U
and branch flags, then evaluate once. Never call solve, _initialize or _newton.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key]='1'
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from multifoil import repanel
from viscous import CoupledSolver, boundary_layer_summary, mf, te_info
from s702_geometry import geometry
from report_s702_wake_extent import require, close, array_close


def recompute(job, spec, state):
    raw=[(n,repanel(x,p)) for (n,x),p in zip(geometry(job['flap']),spec['panels'])]
    solver=CoupledSolver(raw,reynolds=spec['reynolds'],mach=spec['mach'],ncrit=job['ncrit'],
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
        rows.append(dict(element=e,recomputed_CL=float(model.post.cl),recomputed_CD=float(model.post.cd),
            recomputed_CM=float(model.post.cm),recomputed_friction_CD=float(model.post.cdf),
            maximum_Cp_difference=float(np.max(abs(model.post.cp[:model.foil.N][::-1]-surface['Cp'])))))
    return dict(job=job['name'],BL_residual_l2=blnorm,coupling_residual_max=coupling,
        recorded_residual=state['residual'],max_residual=max(blnorm,coupling),
        new_flow_solves=0,elements=rows)


def main():
    import json
    import platform
    import importlib.metadata
    from report_s702_wake_extent import ROOT, Inputs, audit, sha
    from s702_ncrit_sweep import write
    study=ROOT/'validation/results/s702_matched_pressure_20261002'
    out=ROOT/'validation/results/s702_matched_pressure_analysis_20261002'
    out.mkdir(exist_ok=True)
    inputs=Inputs();spec=inputs.json(study/'study.json');runtime=inputs.json(study/'runtime.json')
    require(sha(study/'study.json')=='786752736ff78e82005be0fb22009bbdba50015db57d595ab33bc51126bb26c4',
            'Unexpected matched-pressure study')
    for name,digest in spec['sha256'].items():
        require(sha(inputs.add(ROOT/name))==digest,'Frozen input changed: '+name)
        require(sha(inputs.add(study/'source_snapshot'/name))==digest,'Frozen snapshot changed: '+name)
    for name in ('validation/audit_s702_saved_pressure_state.py','validation/report_s702_wake_extent.py',
                 'validation/report_s702_drag_mechanism.py','validation/report_design_initialization.py'):
        inputs.add(ROOT/name)
    states={j['name']:audit(study/j['name'],j,dict(spec,ncrit=j['ncrit']),runtime,inputs) for j in spec['jobs']}
    require(all(s['result']['feasible'] for s in states.values()),'All eight states required')
    environment=dict(python=platform.python_version(),
        versions={k:importlib.metadata.version(k) for k in ('numpy','scipy','matplotlib')})
    receipt=dict(input_sha256=inputs.files,audit_runtime=environment,states=[],complete=False,
        new_flow_solves=0,method='Restore saved U and evaluate the coupled residual once; no nonlinear solve')
    path=out/'residual_recheck.json'
    if path.exists():
        old=json.loads(path.read_text())
        require(old['input_sha256']==receipt['input_sha256'] and old['audit_runtime']==environment,
                'Recheck inputs/runtime changed; retain the previous receipt')
        receipt=old
    for job in spec['jobs']:
        if any(r['job']==job['name'] for r in receipt['states']):continue
        row=recompute(job,spec,states[job['name']]);receipt['states'].append(row);write(path,receipt)
        print(f"{len(receipt['states'])}/8 {job['name']}: recomputed residual {row['max_residual']:.3g}",flush=True)
    receipt['complete']=True;write(path,receipt)


if __name__=='__main__':main()
