#!/usr/bin/env python3
"""Audit element-specific Ncrit states without nonlinear flow solves.

The export and residual routines are explicitly adapted copies of the frozen
wake-extent auditor and saved-pressure-state recheck. Those original scripts
and receipts remain unchanged. Only threshold provenance/assignment differs.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key]='1'
import json
import platform
import importlib.metadata
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from multifoil import repanel
from viscous import CoupledSolver,boundary_layer_summary,mf,te_info
from s702_geometry import geometry
from report_s702_wake_extent import (ROOT,Inputs,sha,write,require,close,array_close,
    thermo,ELEMENTS,selected_index,audit_surface,side_data)
from run_s702_transition_isolation import check_controls

STUDY=ROOT/'validation/results/s702_transition_isolation_20261002'
PARENT=ROOT/'validation/results/s702_matched_pressure_20261002'
PREVIOUS=ROOT/'validation/results/s702_matched_pressure_analysis_20261002'
OUT=ROOT/'validation/results/s702_transition_isolation_analysis_20261002'
STUDY_HASH='08f05e608a161798a1d38472ad46db6f7121e412d98d3790cb19bd24bbfed9ea'


def thresholds(job):
    values=job.get('ncrit_by_element')
    if values is None:values=dict(fore=job['ncrit'],aft=job['ncrit'])
    require(set(values)=={'fore','aft'} and all(np.isfinite(v) and v>0 for v in values.values()),'Invalid Ncrit map')
    return dict(values)


def check_threshold_provenance(diag,expected):
    close(diag['ncrit'],expected['fore'],context='Legacy scalar must denote fore Ncrit')
    require(diag['ncrit_by_element']==expected,'Per-element Ncrit differs')


def load_states(inputs):
    require(sha(inputs.add(STUDY/'study.json'))==STUDY_HASH,'Unexpected transition study')
    spec=inputs.json(STUDY/'study.json')
    for name,digest in spec['sha256'].items():
        require(sha(inputs.add(ROOT/name))==digest,'Frozen input changed: '+name)
        require(sha(inputs.add(STUDY/'source_snapshot'/name))==digest,'Frozen snapshot changed: '+name)
    receipt=inputs.json(ROOT/'validation/results/s702_transition_isolation_launch_20261002/retrieved_files_sha256.json')
    for name,digest in receipt.items():
        require(sha(inputs.add(STUDY/name))==digest,'Retrieved file changed: '+name)
    runtime=inputs.json(STUDY/'runtime.json');status=inputs.json(STUDY/'status.json')
    automatic=inputs.json(STUDY/'summary.json')
    require(status['state']=='finished' and automatic['converged_jobs']==automatic['required_jobs']==10
            and automatic['execution_checks_pass'],'Incomplete execution')
    require(runtime['python']=='3.12.13' and runtime['versions']==spec['expected_versions'],'Runtime changed')
    require(set(runtime['thread_limits'].values())=={'1'},'Thread limits changed')
    controls=check_controls(STUDY,spec);require(controls['passed'],'Fresh controls failed')
    require(controls==automatic['controls']==inputs.json(STUDY/'controls.json'),'Control checks differ')
    oldspec=inputs.json(PARENT/'study.json');oldruntime=inputs.json(PARENT/'runtime.json')
    require(oldruntime['versions']==runtime['versions'] and oldruntime['python']==runtime['python'],'Baseline runtime differs')
    entries=[]
    for j in spec['jobs']:entries.append(dict(job=j,spec=spec,runtime=runtime,folder=STUDY/j['name'],new=True))
    for j in oldspec['jobs']:
        entries.append(dict(job=j,spec=dict(oldspec,ncrit=j['ncrit']),runtime=oldruntime,folder=PARENT/j['name'],new=False))
    states={e['job']['name']:audit(e['folder'],e['job'],e['spec'],e['runtime'],inputs) for e in entries}
    require(len(states)==18 and all(s['result']['feasible'] for s in states.values()),'All 18 audited states required')
    # Identical fixed geometry is essential to the factorial contrasts.
    for e in entries:
        case=e['job']['pressure_case'];s=states[e['job']['name']]
        ref=next(v for v in entries if not v['new'] and v['job']['pressure_case']==case)
        for element in ELEMENTS:
            array_close(s['arrays'][element]['xy'],states[ref['job']['name']]['arrays'][element]['xy'],'Factorial geometry changed',0.)
    return spec,entries,states,controls


def audit(folder, job, spec, runtime, inputs):
    """Check all exports plus saved arrays, including nonzero-Mach velocity conversion."""
    expected = thresholds(job)
    r = inputs.json(folder/'result.json')
    if 'ncrit_by_element' in job:
        require(r['ncrit_by_element']==expected, 'Result Ncrit map differs')
        if r['feasible']:require(r['selected_solver_ncrit_by_element']==expected, 'Selected solver Ncrit differs')
    meta = inputs.json(folder/'flow/metadata.json'); polar = inputs.csv(folder/'flow/polar.csv')
    require(meta['case']['study'] == spec and meta['case']['job'] == job, 'Export settings differ from manifest')
    require(meta['versions'] == dict(python=runtime['python'], **runtime['versions']), 'Export runtime differs')
    require(len(r['records']) == len(polar) == len(meta['diagnostics']), 'Evaluation counts differ')
    for rec, p, diag in zip(r['records'], polar, meta['diagnostics']):
        require(rec['converged'] == bool(int(p['converged'])) == (diag['status']=='converged'), 'Convergence differs')
        for field in ('alpha_deg', 'CL', 'CD', 'CM', 'residual'):
            if rec[field] is None: require(not np.isfinite(float(p[field])), 'Null evaluation has finite export')
            else: close(rec[field], float(p[field]), context='record/export '+field)
        for key in ('mach', 'reynolds'): close(diag[key], spec[key], context=key)
        check_threshold_provenance(diag, expected)
    if not r['feasible']:
        require(all(r[k] is None for k in ('CL','CD','CM','alpha_deg')) and 'selected_flow_index' not in r,
                'Failed target has selected force/state')
        return dict(result=r, metadata=meta, arrays={}, surfaces={}, wakes={})
    i = selected_index(r, polar); diag = meta['diagnostics'][i]
    require(r['records'][i]['admissible'] and not diag['failure_reason'], 'Selected record is inadmissible')
    residual = float(polar[i]['residual'])
    require(residual <= spec['residual_tolerance'] and abs(r['CL']-job['target_cl']) <= spec['cl_tolerance'], 'Tolerance failure')
    close(residual, max(diag['history'][-1]['bl_residual_l2'], diag['history'][-1]['coupling_residual_max']))
    for key in ('element_forces', 'boundary_layers'): require(r[key] == diag[key], 'Selected diagnostics differ')
    for field, term in (('CL','CL'), ('CD','CD_wake'), ('CM','CM')):
        close(r[field], sum(f[term] for f in r['element_forces'].values()), context='element sum '+field)
    require(r['selected_path'] == r['records'][i]['initialization']['path'], 'Selected path differs')
    expected_paths = [dict(name='cold_local',start_alpha=job['center'],seed_alphas=[]),
                      dict(name='warm_anchor',start_alpha=spec['warm_start_alpha'],seed_alphas=[])]
    if 'entry_offset' in spec:
        expected_paths.append(dict(name='warm_entry',start_alpha=spec['warm_start_alpha'],
                                   seed_alphas=[job['center']+spec['entry_offset']]))
    require(r['initialization_paths'] == expected_paths, 'Declared initialization differs')
    require(r['records'][i]['initialization']['phase'] == 'fixed_lift', 'Seed used as target')
    require(r['attempts'][-1]['feasible'] and sum(a['feasible'] for a in r['attempts']) == 1, 'First success rule differs')
    order = list(dict.fromkeys(a['path'] for a in r['attempts']))
    declared = [p['name'] for p in r['initialization_paths']]
    require(order == declared[:len(order)] and order[-1] == r['selected_path'], 'Initialization reordered')
    for row in r['transition_audit']:
        require(row['inside_cell'] and abs(row['amplification_residual']) <= 1e-7 and row['station_residual_max'] <= 1e-7,
                'Subcell transition audit failed')
    param = thermo(spec); arrays, surfaces, wakes, stations = {}, {}, {}, {}
    for element, panels in zip(ELEMENTS, spec['panels']):
        inputs.add(folder/f'flow/geometry_{element}.dat')
        surface = inputs.array(folder/f'flow/surface_{i:03d}_{element}.csv')
        wake = inputs.array(folder/f'flow/wake_{i:03d}_{element}.csv')
        with np.load(inputs.add(folder/f'{element}_state.npz'), allow_pickle=False) as saved:
            a = {k:saved[k].copy() for k in saved.files}
        n = len(surface)
        require(n == panels+1 and a['U'].shape == (4, n+len(wake)), 'State dimensions differ')
        require(all(np.isfinite(v).all() for v in a.values()), 'Nonfinite Newton state')
        for k, field in ((0,'theta'), (1,'delta_star'), (2,'amplification_or_shear')):
            array_close(a['U'][k,:n][::-1], surface[field], 'Newton/surface '+field)
        for k, field in ((0,'theta'), (1,'delta_star')):
            array_close(a['U'][k,n:], wake[field], 'Newton/wake '+field)
        # U[3] is uncorrected speed; exports contain Karman–Tsien corrected speed.
        uk = mf.get_uk(a['U'][3], param)[0]
        array_close(uk[:n][::-1], surface['ue_over_vinf'], 'Newton/surface corrected velocity')
        array_close(uk[n:], wake['ue_over_vinf'], 'Newton/wake corrected velocity')
        hk = np.array([mf.get_Hk(a['U'][:,j], param)[0] for j in range(n+len(wake))])
        array_close(hk[:n][::-1], surface['Hk'], 'Newton/surface Hk')
        array_close(hk[n:], wake['Hk'], 'Newton/wake Hk')
        require(np.array_equal(a['turbulent'][:n][::-1], surface['turbulent']), 'Turbulence mismatch')
        labels = np.full(n, 'upper', dtype='<U5'); labels[a['lower']] = 'lower'
        require(np.array_equal(labels[::-1], surface['side']), 'Surface indices differ')
        require(set(a['lower']) | set(a['upper']) == set(range(n)) and not set(a['lower']) & set(a['upper']), 'Surface index partition differs')
        for k, field in enumerate(('x_over_cref','y_over_cref')):
            array_close(a['xy'][k][::-1], surface[field], 'Surface geometry differs')
            array_close(a['wake_xy'][k], wake[field], 'Wake geometry differs')
        require(np.all(wake['theta'] > 0) and np.all(wake['ue_over_vinf'] > 0), 'Invalid wake')
        for side in ('upper','lower'): audit_surface(side_data(surface,side), r['boundary_layers'][element][side])
        ds = np.linalg.norm(np.diff(a['wake_xy'], axis=1), axis=0)
        require(np.all(ds > 0), 'Nonpositive wake spacing')
        arclength = np.r_[0., np.cumsum(ds)]
        if 'wake_s' in a: array_close(a['wake_s']-a['wake_s'][0], arclength, 'Wake arclength differs')
        if 'wake_t' in a: array_close(np.linalg.norm(a['wake_t'], axis=0), np.ones(len(wake)), 'Wake tangent normalization')
        if 'wake_discretization' in r:
            info = r['wake_discretization'][element]
            require(info['nodes'] == len(wake), 'Wake node count differs')
            for k, value in (('arclength',arclength[-1]),('first_spacing',ds[0]),('last_spacing',ds[-1])):
                close(info[k],value,context=k)
        # Same Squire–Young estimator evaluated downstream; bookkeeping, not independent drag.
        h = a['U'][1,n:]/a['U'][0,n:]
        sy = 2* a['U'][0,n:] * (uk[n:]/param.Vinf)**((5+h)/2)
        close(sy[-1], r['element_forces'][element]['CD_wake'], context='terminal Squire–Young drag')
        stations[element] = dict(s=arclength, sy=sy)
        arrays[element], surfaces[element], wakes[element] = a, surface, wake
    return dict(result=r, metadata=meta, residual=residual, arrays=arrays, surfaces=surfaces, wakes=wakes, stations=stations)


def recompute(job, spec, state):
    expected_ncrit=thresholds(job)
    raw=[(n,repanel(x,p)) for (n,x),p in zip(geometry(job['flap']),spec['panels'])]
    solver=CoupledSolver(raw,reynolds=spec['reynolds'],mach=spec['mach'],ncrit=expected_ncrit['fore'],
        iterations=spec['iterations'],tolerance=spec['residual_tolerance'],wake_length=job['wake_length'])
    result=state['result'];alpha=result['alpha_deg']
    direction=np.array([np.cos(np.deg2rad(alpha)),np.sin(np.deg2rad(alpha))])
    for model in solver.models:
        model.param.ncrit=expected_ncrit[model.geom.name]
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
        new_flow_solves=0,ncrit_by_element=expected_ncrit,elements=rows)


def main():
    OUT.mkdir(exist_ok=True);inputs=Inputs()
    spec,entries,states,controls=load_states(inputs)
    original=inputs.json(PREVIOUS/'residual_recheck.json')
    require(original['complete'] and len(original['states'])==8,'Previous residual audit incomplete')
    for name,digest in original['input_sha256'].items():
        require(sha(inputs.add(ROOT/name))==digest,'Previous residual evidence changed: '+name)
    for name in ('validation/audit_s702_transition_state.py','validation/report_s702_wake_extent.py',
                 'validation/report_s702_drag_mechanism.py','validation/report_design_initialization.py'):
        inputs.add(ROOT/name)
    environment=dict(python=platform.python_version(),
        versions={k:importlib.metadata.version(k) for k in ('numpy','scipy','matplotlib')})
    receipt=dict(input_sha256=inputs.files,audit_runtime=environment,states=[],complete=False,new_flow_solves=0,
        method='Full export audit of 18 states; fresh residual/force evaluation of 10 new states with actual fore/aft Ncrit; hash-verified reuse of 8 original residual rechecks',
        reused_audit_runtime=original['audit_runtime'],controls=controls)
    path=OUT/'residual_recheck.json'
    if path.exists():
        prior=json.loads(path.read_text())
        require(prior['input_sha256']==receipt['input_sha256'] and prior['audit_runtime']==environment,
                'Audit inputs/runtime changed; preserve prior receipt')
        receipt=prior
    for e in entries:
        j=e['job']
        if any(row['job']==j['name'] for row in receipt['states']):continue
        if e['new']:
            row=recompute(j,e['spec'],states[j['name']]);row['evidence']='fresh_residual_recomputation'
        else:
            row=dict(next(r for r in original['states'] if r['job']==j['name']),
                ncrit_by_element=thresholds(j),evidence='reused_hash_verified_original_residual_recomputation')
        require(row['max_residual']<=spec['residual_tolerance'],'Recomputed residual fails tolerance')
        require(row['ncrit_by_element']==thresholds(j),'Residual receipt threshold provenance differs')
        close(row['recorded_residual'],states[j['name']]['residual'],context='Recorded residual mismatch')
        receipt['states'].append(row);write(path,receipt)
        print(f"{len(receipt['states'])}/18 {j['name']}: residual {row['max_residual']:.3g} ({row['evidence']})",flush=True)
    receipt['complete']=True;write(path,receipt)


if __name__=='__main__':main()
