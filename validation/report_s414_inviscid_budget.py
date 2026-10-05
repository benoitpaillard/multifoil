#!/usr/bin/env python3
"""Same-geometry, same-angle inviscid/viscous pressure split; no nonlinear solve."""
import os
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
from datetime import datetime, timezone
from pathlib import Path
import platform
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np
import scipy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from viscous import CoupledSolver, mf
from s414_geometry import elements
from report_s414_pressure_budget import (Inputs, require, sha, write, csv_write, branch,
    interp, integrate, loads, lift, physical_taps, SURFACES, REGIONS)

PRIOR = ROOT/'validation/results/s414_pressure_budget_20261003'
PRIOR_HASH = '1f23027fa45433222e5cf5417559d0dab5bc6a03221334202fa5facf19529c5b'
INCIDENCE = ROOT/'validation/results/s414_incidence_analysis_20261003'
PREFLIGHT = ROOT/'validation/results/s414_incidence_preflight_20261003'
STUDY = ROOT/'validation/results/s414_incidence_20261003'
OUT = ROOT/'validation/results/s414_inviscid_budget_20261003'
JOBS = ('cl0056_measured_alpha', 'cl028_measured_alpha')


def split_error(inviscid, viscous, experiment):
    a, b, c = (np.asarray(x, float) for x in (inviscid, viscous, experiment))
    require(a.shape == b.shape == c.shape and all(np.isfinite(x).all() for x in (a,b,c)),
            'Inconsistent or nonfinite comparison')
    return a-c, b-a, b-c


def cp_from_raw_speed(speed, mach):
    """Independent Karman–Tsien expression, using uncorrected speed / Vinf."""
    u = np.asarray(speed, float)
    require(np.isfinite(u).all() and np.isfinite(mach) and 0 <= mach < .7, 'Invalid speed/Mach')
    beta = np.sqrt(1-mach*mach)
    cp0 = 1-u*u
    denominator = beta + mach*mach*cp0/(2*(1+beta))
    require(np.all(denominator > 0), 'Singular pressure correction')
    return cp0/denominator


def admit_angle(target):
    require(target['comparison_basis'] == 'fixed_alpha', 'Matched-lift state cannot replace measured-angle case')
    require(abs(target['alpha_deg']-target['reference_caption_alpha']) <= 1e-12,
            'Incidence mismatch')
    require(target['matched_lift_CD_error_counts'] is None, 'Fixed-angle drag silently scored at matched lift')


def freeze_inputs(inputs):
    require(sha(inputs.add(PRIOR/'verification.json')) == PRIOR_HASH, 'Previous budget receipt changed')
    receipt = inputs.json(PRIOR/'verification.json')
    require(receipt['complete'], 'Previous budget incomplete')
    for name, digest in receipt['output_sha256'].items():
        require(sha(inputs.add(PRIOR/name)) == digest, 'Previous budget output changed: '+name)
    for name, digest in inputs.json(PRIOR/'input_manifest.json')['files'].items():
        require(sha(inputs.add(ROOT/name)) == digest, 'Frozen input changed: '+name)
    for name, digest in inputs.json(PREFLIGHT/'input_manifest.json')['files'].items():
        require(sha(inputs.add(ROOT/name)) == digest, 'Preflight input changed: '+name)
    inputs.add(PREFLIGHT/'pressure_at_taps.csv')
    return receipt


def main():
    started = time.monotonic(); OUT.mkdir(exist_ok=True)
    inputs = Inputs(); freeze_inputs(inputs)
    spec = inputs.json(STUDY/'study.json')
    targets = {r['job']: r for r in inputs.json(INCIDENCE/'summary.json')['targets']}
    source_taps = inputs.csv(INCIDENCE/'pressure_at_taps.csv')
    prior_surfaces = inputs.csv(PRIOR/'surface_budgets.csv')
    prior_budgets = inputs.csv(PRIOR/'budgets.csv')
    prior_intervals = inputs.csv(PRIOR/'interval_budgets.csv')
    preflight = inputs.csv(PREFLIGHT/'pressure_at_taps.csv')
    print('Frozen prior inputs verified; building one unchanged inviscid operator.', flush=True)
    solver = CoupledSolver(elements('nominal', spec['panels']), reynolds=spec['reynolds'],
        mach=spec['mach'], ncrit=spec['ncrit'], wake_length=1.)
    # Independently solve the assembled linear equations, including streamfunction constants.
    reference = np.linalg.solve(solver.aic, solver.rhs)
    matrix_residual = float(np.max(np.abs(solver.aic@reference-solver.rhs)))
    basis_difference = float(np.max(np.abs(reference[:solver.nfoil]-solver.gref)))
    require(matrix_residual < 1e-12 and basis_difference < 2e-11, 'Linear operator check failed')
    np.savez_compressed(OUT/'inviscid_operator.npz', matrix=solver.aic, rhs=solver.rhs,
        reference=reference, production_gref=solver.gref,
        fore_xy=solver.models[0].foil.x, aft_xy=solver.models[1].foil.x)
    nodes, taps, surfaces_out, regions, dense, totals, pressure_metrics, checks = [], [], [], [], [], [], [], []
    plot_surfaces = {}
    for job in JOBS:
        t = targets[job]; admit_angle(t); case = t['case']; alpha = t['alpha_deg']
        folder = STUDY/job; result = inputs.json(folder/'result.json'); index = result['selected_flow_index']
        require(result['feasible'] and result['comparison_basis'] == 'fixed_alpha', 'Invalid archived state')
        direction = np.array([np.cos(np.deg2rad(alpha)), np.sin(np.deg2rad(alpha))])
        state_surfaces = {}; full = {'inviscid': np.zeros(2), 'viscous': np.zeros(2)}
        for model in solver.models:
            e = model.geom.name; n = model.foil.N
            s = np.genfromtxt(inputs.add(folder/f'flow/surface_{index:03d}_{e}.csv'), delimiter=',',
                names=True, dtype=None, encoding='utf8')
            with np.load(inputs.add(folder/f'{e}_state.npz'), allow_pickle=False) as saved:
                xy = saved['xy'].copy(); raw_speed = saved['U'][3,:n].copy()
            geom_error = float(np.max(np.abs(model.foil.x-xy)))
            require(geom_error <= 2e-12 and np.array_equal(xy[0,::-1],s['x_over_cref'])
                    and np.array_equal(xy[1,::-1],s['y_over_cref']), 'Rebuilt/archived geometry mismatch')
            gamma = model.isol.gamref@direction
            cp_inv = mf.get_cp(gamma, model.param)[0][::-1]
            cp_error = float(np.max(np.abs(cp_inv-cp_from_raw_speed(gamma, spec['mach'])[::-1])))
            cp_restore = float(np.max(np.abs(s['Cp']-cp_from_raw_speed(raw_speed, spec['mach'])[::-1])))
            require(cp_error < 1e-12 and cp_restore < 1e-12, 'Pressure convention differs')
            fi = loads(s['x_over_cref'],s['y_over_cref'],cp_inv,closed=True)
            fv = loads(s['x_over_cref'],s['y_over_cref'],s['Cp'],closed=True)
            require(abs(lift(fv,alpha)-result['element_forces'][e]['CL']) < 1e-12,
                    'Dense viscous lift differs from saved state')
            full['inviscid'] += fi; full['viscous'] += fv
            dense.append(dict(job=job,case=case,element=e,alpha_deg=alpha,
                inviscid_CN=float(fi[0]),viscous_CN=float(fv[0]),viscous_change_CN=float(fv[0]-fi[0]),
                inviscid_CC_pressure=float(fi[1]),viscous_CC_pressure=float(fv[1]),
                inviscid_CL_pressure=lift(fi,alpha),viscous_CL_pressure=lift(fv,alpha),
                viscous_change_CL_pressure=lift(fv-fi,alpha)))
            checks.append(dict(job=job,element=e,geometry_max_error=geom_error,
                independent_inviscid_Cp_error=cp_error,restored_viscous_Cp_error=cp_restore,
                inviscid_TE_Kutta_error=float(abs(gamma[0]+gamma[-1]))))
            require(checks[-1]['inviscid_TE_Kutta_error'] < 1e-12, 'Kutta condition failed')
            for k, (x,y,ci,cv,g,u) in enumerate(zip(s['x_over_cref'],s['y_over_cref'],cp_inv,s['Cp'],gamma[::-1],raw_speed[::-1])):
                nodes.append(dict(job=job,case=case,element=e,node=k,alpha_deg=alpha,x=float(x),y=float(y),
                    inviscid_signed_speed=float(g),viscous_raw_speed=float(u),inviscid_Cp=float(ci),
                    viscous_Cp=float(cv),viscous_change_Cp=float(cv-ci)))
            inv_surface = s.copy(); inv_surface['Cp'] = cp_inv
            state_surfaces[e] = (s,inv_surface)
        plot_surfaces[job] = state_surfaces
        for coord in ('table','plotted'):
            base = dict(job=job,case=case,coordinate=coord,alpha_deg=alpha)
            partial = {'experiment':np.zeros(2),'inviscid':np.zeros(2),'viscous':np.zeros(2)}
            dense_supported = {'inviscid':np.zeros(2),'viscous':np.zeros(2)}
            for surface in SURFACES:
                e,side = surface.split('_'); s,iv = state_surfaces[e]
                b,ib = branch(s,side),branch(iv,side); gx,gy = b['x_over_cref'],b['y_over_cref']
                raw = [r for r in source_taps if r['job']==job and r['coordinate']==coord and r['surface']==surface]
                supported,tx,ex = physical_taps(raw)
                ci,cv = interp(gx,ib['Cp'],tx),interp(gx,b['Cp'],tx)
                require(np.max(np.abs(cv-[float(r['model_Cp']) for r in supported])) < 1e-13, 'Prior viscous taps differ')
                ie,change,ve = split_error(ci,cv,ex)
                for k,r in enumerate(supported):
                    taps.append(dict(**base,surface=surface,table_index=int(r['table_index']),x=float(tx[k]),
                        experiment_Cp=float(ex[k]),inviscid_Cp=float(ci[k]),viscous_Cp=float(cv[k]),
                        inviscid_error_Cp=float(ie[k]),viscous_change_Cp=float(change[k]),viscous_error_Cp=float(ve[k])))
                    if coord=='table':
                        old = next(q for q in preflight if q['case']==case and q['incidence']=='measured_caption'
                            and q['surface']==surface and q['table_index']==r['table_index'])
                        require(abs(float(old['alpha_deg'])-alpha)<1e-12 and float(old['x'])==tx[k]
                            and abs(float(old['inviscid_Cp'])-ci[k])<2e-12, 'Archived inviscid preflight differs')
                forces = {name:integrate(gx,gy,tx,cp,tx[0],tx[-1],side)
                          for name,cp in [('experiment',ex),('inviscid',ci),('viscous',cv)]}
                for name,f in forces.items(): partial[name] += f
                for name,cp in [('inviscid',ib['Cp']),('viscous',b['Cp'])]:
                    dense_supported[name] += integrate(gx,gy,gx,cp,tx[0],tx[-1],side)
                ie,change,ve = split_error(forces['inviscid'],forces['viscous'],forces['experiment'])
                old = next(r for r in prior_surfaces if r['job']==job and r['coordinate']==coord and r['surface']==surface)
                require(abs(ve[0]-float(old['delta_CN']))<1e-12, 'Prior surface budget changed')
                surfaces_out.append(dict(**base,surface=surface,taps=len(tx),x_first=float(tx[0]),x_last=float(tx[-1]),
                    experiment_CN=float(forces['experiment'][0]),inviscid_CN=float(forces['inviscid'][0]),
                    viscous_CN=float(forces['viscous'][0]),inviscid_error_CN=float(ie[0]),
                    viscous_change_CN=float(change[0]),viscous_error_CN=float(ve[0]),
                    inviscid_error_CL_pressure=lift(ie,alpha),viscous_change_CL_pressure=lift(change,alpha),
                    viscous_error_CL_pressure=lift(ve,alpha)))
                for region,lo,hi in REGIONS[e]:
                    left,right = max(lo,tx[0]),min(hi,tx[-1])
                    if right <= left: continue
                    f = {name:integrate(gx,gy,tx,cp,left,right,side)
                         for name,cp in [('experiment',ex),('inviscid',ci),('viscous',cv)]}
                    ie,change,ve = split_error(f['inviscid'],f['viscous'],f['experiment'])
                    prior_value = sum(float(r['delta_CN']) for r in prior_intervals if
                        r['job']==job and r['coordinate']==coord and r['surface']==surface and r['region']==region)
                    require(abs(ve[0]-prior_value)<1e-12, 'Prior region masks differ')
                    regions.append(dict(**base,surface=surface,region=region,left_x=float(left),right_x=float(right),
                        inviscid_error_CN=float(ie[0]),viscous_change_CN=float(change[0]),viscous_error_CN=float(ve[0])))
            old = next(r for r in prior_budgets if r['job']==job and r['coordinate']==coord)
            ie,change,ve = split_error(partial['inviscid'],partial['viscous'],partial['experiment'])
            exp_cn = float(old['reported_appendix_CN'])
            require(abs(ve[0]-float(old['delta_partial_CN']))<1e-12, 'Prior total budget differs')
            totals.append(dict(**base,taps=87,reported_appendix_CN=exp_cn,
                inviscid_CL=lift(full['inviscid'],alpha),viscous_CL=lift(full['viscous'],alpha),
                experiment_partial_CN=float(partial['experiment'][0]),inviscid_partial_CN=float(partial['inviscid'][0]),
                viscous_partial_CN=float(partial['viscous'][0]),inviscid_error_partial_CN=float(ie[0]),
                viscous_change_partial_CN=float(change[0]),viscous_error_partial_CN=float(ve[0]),
                inviscid_full_CN=float(full['inviscid'][0]),viscous_full_CN=float(full['viscous'][0]),
                inviscid_error_full_CN=float(full['inviscid'][0]-exp_cn),
                viscous_change_full_CN=float(full['viscous'][0]-full['inviscid'][0]),
                viscous_error_full_CN=float(full['viscous'][0]-exp_cn),
                inviscid_sampling_correction_CN=float(dense_supported['inviscid'][0]-partial['inviscid'][0]),
                viscous_sampling_correction_CN=float(dense_supported['viscous'][0]-partial['viscous'][0]),
                inviscid_uncovered_CN=float(full['inviscid'][0]-dense_supported['inviscid'][0]),
                viscous_uncovered_CN=float(full['viscous'][0]-dense_supported['viscous'][0]),
                experimental_unresolved_CN=float(old['experimental_unresolved_CN'])))
            for surface in (*SURFACES,'all_surfaces'):
                selected = [r for r in taps if r['job']==job and r['coordinate']==coord and
                            (surface=='all_surfaces' or r['surface']==surface)]
                for method in ('inviscid','viscous'):
                    err = np.array([r[method+'_error_Cp'] for r in selected])
                    pressure_metrics.append(dict(**base,surface=surface,method=method,taps=len(selected),
                        rmse_Cp=float(np.sqrt(np.mean(err*err))),mean_error_Cp=float(np.mean(err))))
        print(job, 'inviscid / viscous CL:',lift(full['inviscid'],alpha),lift(full['viscous'],alpha),flush=True)
    for name,rows in [('nodes',nodes),('pressure_at_taps',taps),('surface_budgets',surfaces_out),
                      ('region_budgets',regions),('dense_forces',dense),('totals',totals),('pressure_metrics',pressure_metrics),
                      ('state_checks',checks)]:
        csv_write(OUT/(name+'.csv'),rows)
    write(OUT/'summary.json',dict(complete=True,comparison='Same geometry, measured caption angle, Mach and pressure convention',
        cases=2,operator_rebuilds=1,inviscid_angle_evaluations=2,new_nonlinear_flow_solves=0,
        archived_preflight_taps_repeated=174,table_and_plotted_rows=len(taps),matrix_residual=matrix_residual,
        independent_basis_difference=basis_difference,geometry_max_error=max(r['geometry_max_error'] for r in checks),
        maximum_restored_Cp_error=max(r['restored_viscous_Cp_error'] for r in checks),
        exact_tested_geometry_established=False,numerical_convergence_established=False,
        pure_inviscid_experimental_reference=False,isolated_transition_or_wake_intervention=False,
        new_mses_runs=0,totals=totals,wall_seconds=time.monotonic()-started))
    plot(targets,plot_surfaces,taps,surfaces_out,totals)
    for path in (Path(__file__),ROOT/'test_s414_inviscid_budget.py',ROOT/'validation/report_s414_pressure_budget.py',
                 ROOT/'viscous.py',ROOT/'vendor/mfoil.py',ROOT/'multifoil.py',ROOT/'paneling.py',ROOT/'validation/s414_geometry.py'):
        inputs.add(path)
    write(OUT/'input_manifest.json',dict(created_utc=datetime.now(timezone.utc).isoformat(),files=inputs.files,
        runtime=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,matplotlib=matplotlib.__version__),
        thread_limits={k:os.environ[k] for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')}))


def save(fig,name):
    for suffix in ('png','pdf','svg'): fig.savefig(OUT/(name+'.'+suffix),dpi=170,bbox_inches='tight')
    plt.close(fig)


def plot(targets,states,taps,surfaces,totals):
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    for detail in (False,True):
        names=('fore_lower','aft_upper') if detail else SURFACES
        fig,axes=plt.subplots(len(names),2,figsize=(12,8 if detail else 13),layout='constrained')
        for col,job in enumerate(JOBS):
            t=targets[job]; total=next(r for r in totals if r['job']==job and r['coordinate']=='table')
            for row,surface in enumerate(names):
                ax=axes[row,col]; e,side=surface.split('_'); s,iv=states[job][e]
                for data,color,style,label in [(iv,'#b45a00','--',f'Inviscid, CL {total["inviscid_CL"]:.4f}'),
                                               (s,'#2364aa','-',f'Viscous, CL {total["viscous_CL"]:.4f}')]:
                    b=branch(data,side);ax.plot(b['x_over_cref'],b['Cp'],style,color=color,label=label)
                raw=[r for r in taps if r['job']==job and r['coordinate']=='table' and r['surface']==surface]
                ax.plot([r['x'] for r in raw],[r['experiment_Cp'] for r in raw],'ko',ms=3.5,label=f'Experiment, CL {t["reference_caption_CL"]:.3f}')
                if detail: ax.set_xlim((0,.06) if row==0 else (.80,1.005))
                ax.invert_yaxis();ax.grid(alpha=.2);ax.legend(fontsize=7)
                ax.set(title=f'α {t["alpha_deg"]:.2f}° · {surface.replace("_"," ")}',xlabel='x/c',ylabel='Cp')
        fig.suptitle('Same geometry and incidence: inviscid baseline versus saved viscous solution\nMach 0.10 pressure correction is identical; experimental data are viscous; no matched-angle MSES curve is available')
        save(fig,'local_pressure_response' if detail else 'full_pressure_split')
    fig,axes=plt.subplots(1,2,figsize=(12,5.3),layout='constrained')
    for ax,job in zip(axes,JOBS):
        group=[r for r in surfaces if r['job']==job and r['coordinate']=='table']
        for offset,field,label,color in [(-.25,'inviscid_error_CN','Inviscid − experiment','#b45a00'),
                                        (0,'viscous_change_CN','Viscous − inviscid','#777777'),
                                        (.25,'viscous_error_CN','Viscous − experiment','#2364aa')]:
            values=[next(r[field] for r in group if r['surface']==s) for s in SURFACES]
            bars=ax.bar(np.arange(4)+offset,values,width=.23,color=color,label=label)
            ax.bar_label(bars,fmt='%+.3f',fontsize=7,padding=3)
        ax.axhline(0,color='k',lw=.8);ax.grid(axis='y',alpha=.2);ax.legend(fontsize=8)
        ax.set_xticks(range(4),[s.replace('_',' ') for s in SURFACES],rotation=18)
        ax.set(title=f'α {targets[job]["alpha_deg"]:.2f}°',ylabel='Partial normal-force coefficient difference',ylim=(-.046,.093))
    fig.suptitle('Viscous interaction unloads both elements but leaves a net excess\nIdentical 87 supported taps per case; brown + gray = blue; these are load differences, not drag')
    save(fig,'load_decomposition')


if __name__=='__main__': main()
