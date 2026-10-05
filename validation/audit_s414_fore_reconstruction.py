#!/usr/bin/env python3
"""Audit three selected states and the failed fourth target without filling its forces."""
import os
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
import platform
import importlib.metadata
import numpy as np
from audit_s414_state import ROOT, Inputs, sha, audit, require, close, array_close, elements, mf, te_info, boundary_layer_summary
from amplification_jacobian import AmplificationJacobianSolver
from run_s414_fore_reconstruction import controls, geometry_elements as elements
from s702_ncrit_sweep import write

STUDY = ROOT/'validation/results/s414_fore_reconstruction_20261004'
OUT = ROOT/'validation/results/s414_fore_reconstruction_analysis_20261004'
LAUNCH = ROOT/'validation/results/s414_fore_reconstruction_launch_20261004'
STUDY_HASH = 'b512f46a3704a474daf54e256c62188aa7ed90450c2de1003ec242dddfbbec74'


def load_states(inputs):
    require(sha(inputs.add(STUDY/'study.json')) == STUDY_HASH, 'Study identity changed')
    spec = inputs.json(STUDY/'study.json')
    for name, digest in spec['sha256'].items():
        require(sha(inputs.add(ROOT/name)) == digest, 'Frozen input changed: '+name)
        require(sha(inputs.add(STUDY/'source_snapshot'/name)) == digest, 'Snapshot changed: '+name)
    retrieved = inputs.json(LAUNCH/'retrieved_files_sha256.json')
    for name, digest in retrieved.items():
        require(sha(inputs.add(STUDY/name)) == digest, 'Retrieved file changed: '+name)
    receipt = inputs.json(LAUNCH/'retrieval.json')
    require(receipt['all_files_match'] and receipt['files_verified'] == len(retrieved) == 430, 'Retrieval incomplete')
    runtime = inputs.json(STUDY/'runtime.json'); status = inputs.json(STUDY/'status.json')
    automatic = inputs.json(STUDY/'summary.json')
    require(status['state'] == 'finished' and status['completed_jobs'] == status['total_jobs'] == 4, 'Run incomplete')
    require(not automatic['execution_checks_pass'] and automatic['required_jobs'] == 4 and automatic['converged_jobs'] == 3, 'Recorded coverage changed')
    require(runtime['python'] == '3.12.13' and runtime['versions'] == spec['expected_versions'], 'Runtime changed')
    require(set(runtime['thread_limits'].values()) == {'1'}, 'Thread settings changed')
    control = controls(STUDY, spec)
    require(control['passed'] and control == automatic['controls'] == inputs.json(STUDY/'controls.json'), 'Control differs')
    states = {}
    for job in spec['jobs']:
        state = audit(STUDY/job['name'], job, spec, runtime, inputs)
        require(state['result']['scenario'] == job['scenario'], 'Result geometry scenario differs')
        require(state['result']['target_cl'] == job['target_cl'], 'Result lift target differs')
        for element, raw in elements(job['scenario'], spec['panels']):
            exported = np.loadtxt(inputs.add(STUDY/job['name']/'flow'/('geometry_'+element+'.dat')), skiprows=1)
            array_close(exported, raw, 'Exported input geometry differs', 5e-12)
        mode = job['amplification_jacobian']
        require(state['result']['amplification_jacobian'] == mode, 'Result mode differs')
        require(all(d['amplification_jacobian'] == mode for d in state['metadata']['diagnostics']), 'Evaluation mode differs')
        states[job['name']] = state
    for name in ('validation/audit_s414_fore_reconstruction.py', 'validation/audit_s414_state.py', 'amplification_jacobian.py'):
        inputs.add(ROOT/name)
    require(sum(s['result']['feasible'] for s in states.values()) == 3, 'Coverage inconsistent')
    for name, state in states.items():
        if state['result']['feasible']: continue
        for element in ('fore', 'aft'):
            require(not (STUDY/name/(element+'_state.npz')).exists(), 'Failed target has selected arrays')
    return spec, states, control


def main():
    OUT.mkdir(exist_ok=True); inputs = Inputs(); spec, states, control = load_states(inputs)
    receipt = dict(complete=False, input_sha256=inputs.files, states=[], controls=control,
        new_nonlinear_flow_solves=0, runtime=dict(python=platform.python_version(),
            versions={k: importlib.metadata.version(k) for k in ('numpy', 'scipy', 'matplotlib')}))
    path = OUT/'residual_recheck.json'
    if path.exists():
        old = inputs.json(path)
        require(old['input_sha256'] == receipt['input_sha256'] and old['runtime'] == receipt['runtime'], 'Preserve changed audit')
        receipt = old
    for job in spec['jobs']:
        if not states[job['name']]['result']['feasible']: continue
        if any(r['job'] == job['name'] for r in receipt['states']): continue
        row = recompute(job, spec, states[job['name']]); row['amplification_jacobian'] = job['amplification_jacobian']
        receipt['states'].append(row); write(path, receipt)
        print(job['name'], 'fresh residual', row['max_residual'], flush=True)
    receipt['complete'] = True; receipt['target_coverage_complete'] = False
    receipt['selected_states_audited'] = len(receipt['states'])
    receipt['failed_targets'] = [j['name'] for j in spec['jobs'] if not states[j['name']]['result']['feasible']]
    write(path, receipt)


# Fresh recomputation below is adapted from the unchanged prior auditor.

def recompute(job, spec, state):
    raw=elements(job['scenario'],spec['panels'])
    solver=AmplificationJacobianSolver(raw,amplification_jacobian=job['amplification_jacobian'],reynolds=spec['reynolds'],mach=spec['mach'],ncrit=spec['ncrit'],
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

if __name__ == '__main__':
    main()
