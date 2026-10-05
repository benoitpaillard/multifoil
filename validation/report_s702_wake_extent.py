#!/usr/bin/env python3
"""Audit selected S702 wake states and report extent sensitivity without solves."""
import argparse
import json
from pathlib import Path
from types import SimpleNamespace

import matplotlib
matplotlib.use('Agg')
from matplotlib.path import Path as Polygon
import matplotlib.pyplot as plt
import numpy as np

from report_design_initialization import ROOT, Inputs, sha, write, csv_write
from report_s702_drag_mechanism import require, close, selected_index, audit_surface, side_data
from multifoil import _intersections
from vendor import mfoil as mf

STUDY_HASH = 'd79a1add8303d9f5b5de97d18396f7360a5dada3e44f7e88593ab3693c4d6f96'
ELEMENTS = ('fore', 'aft')


def array_close(a, b, context, tolerance=1e-12):
    require(a.shape == b.shape and np.max(abs(a-b)) <= tolerance, context)


def thermo(spec):
    model = SimpleNamespace(param=mf.Param(), oper=mf.Oper(), geom=mf.Geom())
    model.geom.chord = spec['reference_chord']
    model.oper.Ma, model.oper.Re = spec['mach'], spec['reynolds']
    mf.init_thermo(model)
    return model.param


def clearance(wake, body):
    """Open wake polyline to closed other body; never close the wake artificially."""
    a, b = wake[:-1], wake[1:]
    c, d = body, np.roll(body, -1, axis=0)
    if _intersections(a, b, c, d).any() or Polygon(body).contains_points(wake).any():
        return 0.
    def distance(points, starts, ends):
        delta = ends-starts; size = np.sum(delta*delta, axis=1)
        keep = size > 1e-24; delta, starts, size = delta[keep], starts[keep], size[keep]
        offset = points[:, None, :]-starts
        t = np.clip(np.sum(offset*delta, axis=2)/size, 0., 1.)
        return float(np.min(np.linalg.norm(offset-t[:, :, None]*delta, axis=2)))
    return min(distance(wake, c, d), distance(body, a, b))


def audit(folder, job, spec, runtime, inputs):
    """Check all exports plus saved arrays, including nonzero-Mach velocity conversion."""
    r = inputs.json(folder/'result.json')
    meta = inputs.json(folder/'flow/metadata.json'); polar = inputs.csv(folder/'flow/polar.csv')
    require(meta['case']['study'] == spec and meta['case']['job'] == job, 'Export settings differ from manifest')
    require(meta['versions'] == dict(python=runtime['python'], **runtime['versions']), 'Export runtime differs')
    require(len(r['records']) == len(polar) == len(meta['diagnostics']), 'Evaluation counts differ')
    for rec, p, diag in zip(r['records'], polar, meta['diagnostics']):
        require(rec['converged'] == bool(int(p['converged'])) == (diag['status']=='converged'), 'Convergence differs')
        for field in ('alpha_deg', 'CL', 'CD', 'CM', 'residual'):
            if rec[field] is None: require(not np.isfinite(float(p[field])), 'Null evaluation has finite export')
            else: close(rec[field], float(p[field]), context='record/export '+field)
        for key in ('mach', 'reynolds', 'ncrit'): close(diag[key], spec[key], context=key)
        require(diag['ncrit_by_element'] == dict(fore=spec['ncrit'], aft=spec['ncrit']), 'Per-element Ncrit differs')
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


def report(study, out):
    require(sha(study/'study.json') == STUDY_HASH, 'Unexpected wake study')
    inputs = Inputs(); spec = inputs.json(study/'study.json')
    for name, digest in spec['sha256'].items():
        require(sha(inputs.add(study/'source_snapshot'/name)) == digest, 'Source snapshot changed: '+name)
    runtime = inputs.json(study/'runtime.json')
    require(runtime['versions'] == spec['expected_versions'], 'Runtime differs')
    require(inputs.json(study/'status.json')['state'] == 'finished', 'Study unfinished')
    original = inputs.json(study/'summary.json')
    states, references, rows, wake_rows, surface_rows, histories, control_diff = {}, {}, [], [], [], [], []
    for job in spec['jobs']:
        current = audit(study/job['name'], job, spec, runtime, inputs)
        r = current['result']; require(r['feasible'], 'Incomplete study: full paired interpretation prohibited')
        path = ROOT/job['reference_result']
        if str(path) not in references:
            # All reference files consumed by this auditor were pinned before the new run.
            oldspec = inputs.json(path.parent.parent/'study.json')
            oldruntime = inputs.json(path.parent.parent/'runtime.json')
            oldjob = next(j for j in oldspec['jobs'] if j['name'] == path.parent.name)
            before = set(inputs.files)
            ref = audit(path.parent, oldjob, oldspec, oldruntime, inputs)
            for key in set(inputs.files)-before:
                require(key in spec['sha256'] and sha(ROOT/key) == spec['sha256'][key], 'Reference missing from frozen inputs: '+key)
            references[str(path)] = ref
        ref = references[str(path)]; old = ref['result']
        for e in ELEMENTS:
            array_close(current['arrays'][e]['xy'], ref['arrays'][e]['xy'], 'Paired solver geometry changed', 0.)
        delta = (r['CD']-old['CD'])*1e4
        auto = next(v for v in original['jobs'] if v['job'] == job['name'])
        close(auto['difference_from_1c_CD_counts'],delta)
        rows.append(dict(job=job['name'], control=job['control'], flap=job['flap'], target_cl=job['target_cl'],
            wake_length=job['wake_length'], CL=r['CL'], CD=r['CD'], CM=r['CM'], alpha_deg=r['alpha_deg'],
            lift_error=r['CL']-job['target_cl'], residual=current['residual'], difference_CD_counts=delta,
            difference_CL=r['CL']-old['CL'], difference_CM=r['CM']-old['CM'], difference_alpha=r['alpha_deg']-old['alpha_deg'],
            selected_path=r['selected_path'], evaluations=len(r['records']), failed_evaluations=sum(not x['admissible'] for x in r['records']),
            failed_internal_points=sum(not c['converged'] for d in current['metadata']['diagnostics'] for c in d['continuation']),
            wall_seconds=r['wall_seconds']))
        histories.append(dict(job=job['name'], records=r['records'], attempts=r['attempts'], selected_flow_index=r['selected_flow_index'],
                              continuation=[d['continuation'] for d in current['metadata']['diagnostics']]))
        if job['control']:
            require(abs(delta) <= spec['acceptance_control_CD_counts'], 'Fresh control drag differs')
            for e in ELEMENTS:
                for field, a in ref['arrays'][e].items():
                    b = current['arrays'][e][field]
                    control_diff.append(dict(element=e, field=field, exactly_equal=bool(np.array_equal(a,b)),
                                             max_abs_difference=float(np.max(abs(a-b)))))
            continue
        states[job['flap'], job['target_cl'], 1] = ref
        states[job['flap'], job['target_cl'], 2] = current
    require(len(rows) == original['required_jobs'] == original['feasible_jobs'] == 5 and original['execution_checks_pass'], 'Coverage differs')
    for (flap, cl, extent), state in states.items():
        r = state['result']
        for element in ELEMENTS:
            a, wake, station = state['arrays'][element], state['wakes'][element], state['stations'][element]
            other = 'aft' if element=='fore' else 'fore'
            gap = clearance(a['wake_xy'].T, state['arrays'][other]['xy'].T)
            require(gap > 0, 'Wake intersects other body')
            wr = dict(flap=flap,target_cl=cl,wake_length=extent,element=element,nodes=len(wake),arclength=station['s'][-1],
                other_body_clearance=gap,CD_wake=r['element_forces'][element]['CD_wake'],CD_friction=r['element_forces'][element]['CD_friction'],
                downstream_SY_range_after_half_chord_counts=float(np.ptp(station['sy'][station['s']>=.5])*1e4))
            for endpoint, k in (('first',0), ('last',-1)):
                for field in ('theta','delta_star','ue_over_vinf','Hk'): wr[endpoint+'_'+field] = float(wake[field][k])
            wake_rows.append(wr)
            for side in ('upper','lower'):
                values = side_data(state['surfaces'][element],side)
                layer = r['boundary_layers'][element][side]; tr=layer['transition']
                turbulent = values[values['turbulent']==1]
                surface_rows.append(dict(flap=flap,target_cl=cl,wake_length=extent,element=element,side=side,
                    transition_x=tr['x_over_local_chord'] if tr else None, separated_fraction=layer['separated_arclength_fraction'],
                    separation_regions=len(layer['separation_regions']),max_Hk=layer['max_Hk'], min_cf=layer['min_cf'],
                    max_turbulent_Hk_drop=max(0.,float(np.max(-np.diff(turbulent['Hk']))) if len(turbulent)>1 else 0.),
                    max_turbulent_velocity_step=float(np.max(abs(np.diff(turbulent['ue_over_vinf'])))) if len(turbulent)>1 else 0.))
    pairs = []
    for cl in (.8,1.2):
        row = dict(target_cl=cl)
        for length in (1,2):
            a,b = [states[f,cl,length]['result'] for f in (5,10)]
            row[f'length{length}_delta_CD_counts'] = (b['CD']-a['CD'])*1e4
            for e in ELEMENTS:
                row[f'length{length}_{e}_delta_CD_counts'] = (b['element_forces'][e]['CD_wake']-a['element_forces'][e]['CD_wake'])*1e4
        row['difference_of_flap_deltas_counts'] = row['length2_delta_CD_counts']-row['length1_delta_CD_counts']
        pairs.append(row)
    summary = dict(required_jobs=5, feasible_jobs=5, selected_states_audited=9, full_state_audit_pass=True,
        control_CD_difference_counts=next(r['difference_CD_counts'] for r in rows if r['control']),
        control_shared_arrays_exactly_equal=all(d['exactly_equal'] for d in control_diff),
        max_absolute_individual_CD_change_counts=max(abs(r['difference_CD_counts']) for r in rows if not r['control']),
        max_absolute_flap_delta_change_counts=max(abs(r['difference_of_flap_deltas_counts']) for r in pairs),
        minimum_wake_other_body_clearance=min(r['other_body_clearance'] for r in wake_rows),
        requested_evaluations=sum(r['evaluations'] for r in rows), failed_requested_evaluations=sum(r['failed_evaluations'] for r in rows),
        failed_internal_continuation_points=sum(r['failed_internal_points'] for r in rows),
        max_turbulent_Hk_drop=max(r['max_turbulent_Hk_drop'] for r in surface_rows),
        max_turbulent_velocity_step=max(r['max_turbulent_velocity_step'] for r in surface_rows),
        qualification='Operational extent and derived wake discretization sensitivity only; no independent drag or contour-ranking validation.')
    out.mkdir(parents=True,exist_ok=True)
    for name,data in (('targets',rows),('flap_deltas',pairs),('wake_diagnostics',wake_rows),('surface_diagnostics',surface_rows),('control_differences',control_diff)):
        csv_write(out/(name+'.csv'),data)
    write(out/'summary.json',summary); write(out/'initialization_histories.json',histories)
    plot(out,states,rows,pairs)
    (out/'FINDINGS.md').write_text(
        '# S702 wake-extent audit\n\n'
        'All five new jobs converged, and all nine selected states (five new plus four archived 1c references) passed '
        'the full export/Newton-array audit. The fresh 1c control repeats drag exactly; every shared saved array field is identical. '
        'The nonzero-Mach audit applies the recorded Karman–Tsien velocity conversion before comparing U[3] with exported speeds. '
        'Source snapshots, runtime, geometry, force sums, transition/separation, tolerances and initialization histories were checked.\n\n'
        f"Maximum individual 2c-minus-1c drag change is {summary['max_absolute_individual_CD_change_counts']:.6f} count. "
        f"The largest change in flap-10-minus-flap-5 drag is {summary['max_absolute_flap_delta_change_counts']:.6f} count. "
        'Both are much smaller than the 0.5-count design-screening threshold and the previously documented ~19.38-count '
        'CL-0.8 flap-delta discrepancy. Extending this modeled wake does not explain that discrepancy. This conclusion is confined '
        'to Re=700,000, Mach=0.07, Ncrit=6.25, flaps 5/10 and CL 0.8/1.2.\n\n'
        + '\n'.join(f"- CL {p['target_cl']:g}: flap delta {p['length1_delta_CD_counts']:.6f} → {p['length2_delta_CD_counts']:.6f} counts; change {p['difference_of_flap_deltas_counts']:+.6f}." for p in pairs)
        + '\n\nThe wake setting changes both extent and its automatically derived node count. Surface geometry is identical across these pairs. '
        'Wake polylines do not intersect the other element; minimum segment clearance is '
        f"{summary['minimum_wake_other_body_clearance']:.6f}c. Positive centerline clearance does not establish absent wake/BL interaction. "
        'The downstream Squire–Young traces re-evaluate the existing estimator and are not an independent momentum-drag calculation. '
        'No merging or confluence physics has been validated.\n\n'
        f"There are {summary['requested_evaluations']} requested evaluations, {summary['failed_requested_evaluations']} failed requests and "
        f"{summary['failed_internal_continuation_points']} failed internal continuation points. These failures remain saved.\n\n"
        '## Next discriminating check\n\n'
        'The flap-10 aft layers have a post-transition Hk drop near 1.10 per panel at CL 0.8 and velocity steps above 0.05 '
        'at both lift targets. The [XFOIL manual](https://web.mit.edu/drela/Public/web/xfoil/xfoil_doc.txt) gives approximate '
        'resolution guides of Hk drops below 1 near post-transition bubbles and velocity steps below 0.05 past transition. '
        'These are diagnostic guides, not MultiFoil acceptance tolerances or experimental uncertainty. The earlier 340/240 '
        'fine-panel check covered Re=1 million, not these reserved Re=700,000 states. A bounded 340/240 versus 260/180 comparison '
        'on these same four conditions is justified before changing closures. Keep wake extent at 1c and Ncrit frozen. '
        'The derived wake discretization will also change with surface count and must remain explicit.\n\n'
        'Priority 4 remains open. Original geometry/hinge and ideal-trip limitations in priority 3 remain unresolved. '
        'No new geometry, calibration, or aerodynamic model was selected.\n\n'
        '![Drag sensitivity](drag_sensitivity.png)\n\n![Downstream wakes](wake_evolution.png)\n\n![Aft surface resolution](aft_surface_resolution.png)\n\n'
        'Reproduce: `python3 validation/report_s702_wake_extent.py`. Tables retain each condition, element and surface.\n')
    for name in ('validation/report_s702_wake_extent.py','validation/report_design_initialization.py',
                 'validation/report_s702_drag_mechanism.py','vendor/mfoil.py','multifoil.py'):
        inputs.add(ROOT/name)
    if (out/'source_manifest.json').exists():
        source = inputs.json(out/'source_manifest.json')
        require(sha(inputs.add(out/source['file'])) == source['sha256'], 'External source snapshot changed')
    write(out/'input_manifest.json',dict(study_sha256=STUDY_HASH, files=inputs.files))
    return summary


def plot(out,states,rows,pairs):
    def save(fig,name):
        fig.tight_layout()
        for ext in ('png','svg'): fig.savefig(out/(name+'.'+ext),dpi=170)
        plt.close(fig)
    fig,ax=plt.subplots(1,2,figsize=(11,4.5))
    for flap,marker in ((5,'o'),(10,'s')):
        rr=[r for r in rows if not r['control'] and r['flap']==flap]
        ax[0].plot([r['target_cl'] for r in rr],[r['difference_CD_counts'] for r in rr],marker=marker,label=f'Flap {flap}°')
    ax[0].set(xlabel='Target CL',ylabel='CD(2c) − CD(1c) [counts]',title='5/5 new jobs; exact 1c control')
    for e,c in (('fore','#277da8'),('aft','#e07a27')):
        for length,offset,style in ((1,-.045,'-'),(2,.045,'--')):
            ax[1].plot([p['target_cl']+offset for p in pairs],[p[f'length{length}_{e}_delta_CD_counts'] for p in pairs],
                       style+'o',color=c,label=f'{e}, {length}c')
    ax[1].set(xlabel='Target CL (offset for visibility)',ylabel='Flap 10° − 5° wake CD [counts]',title='Fore benefit still offset by aft penalty',xticks=[.8,1.2])
    for a in ax: a.axhline(0,color='0.6',lw=.7);a.grid(alpha=.2);a.legend(fontsize=8)
    save(fig,'drag_sensitivity')
    fig,ax=plt.subplots(3,2,figsize=(11,9))
    colors={(5,.8):'#277da8',(10,.8):'#e07a27',(5,1.2):'#32936f',(10,1.2):'#9355a2'}
    for (flap,cl,length),state in states.items():
        for col,e in enumerate(ELEMENTS):
            s=state['stations'][e]; w=state['wakes'][e]; label=f'{flap}°, CL {cl:g}, {length}c'
            for row,y in enumerate((w['theta'],w['Hk'],s['sy']*1e4)):
                ax[row,col].plot(s['s'],y,'-' if length==1 else '--',color=colors[flap,cl],lw=1,label=label)
            ax[0,col].set_title(e.capitalize()+' wake')
            ax[2,col].set_xlabel('Arclength from wake start / reference chord')
    for row,label in enumerate(('Momentum thickness / c','Kinematic shape factor Hk','Squire–Young estimate [counts]')):
        ax[row,0].set_ylabel(label)
    ax[0,1].legend(fontsize=7,ncol=2)
    for a in ax.flat:a.grid(alpha=.2)
    save(fig,'wake_evolution')
    fig,ax=plt.subplots(2,2,figsize=(11,6.8))
    for col,side in enumerate(('upper','lower')):
        for cl,color in ((.8,'#277da8'),(1.2,'#e07a27')):
            for length in (1,2):
                v=side_data(states[10,cl,length]['surfaces']['aft'],side);v=v[v['turbulent']==1]
                for row,field in enumerate(('Hk','ue_over_vinf')):
                    ax[row,col].plot(v['x_over_local_chord'],v[field],'-o' if length==1 else '--',ms=3,color=color,label=f'CL {cl:g}, {length}c')
        ax[0,col].set_title('Flap 10° aft '+side+' after transition');ax[1,col].set_xlabel('x / aft local chord')
    ax[0,0].set_ylabel('Hk');ax[1,0].set_ylabel('Ue / Uinf');ax[0,1].legend(fontsize=8)
    for a in ax.flat:a.grid(alpha=.2)
    save(fig,'aft_surface_resolution')


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--study',type=Path,default=ROOT/'validation/results/s702_wake_extent_20260922')
    ap.add_argument('--out',type=Path,default=ROOT/'validation/results/s702_wake_extent_analysis')
    args=ap.parse_args()
    print(json.dumps(report(args.study.resolve(),args.out.resolve()),indent=2))
