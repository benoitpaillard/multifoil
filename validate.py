#!/usr/bin/env python3
"""Run published validation cases or compare archived results, without fitting.

Run from any directory. See validation/README.md for provenance and limitations.
Exit 0: every requested solve and reference comparison available; 2: incomplete;
1: invalid input. Exit 0 does NOT mean the aerodynamic model is validated.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import sys

os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np

ROOT = Path(__file__).resolve().parent
DATA = ROOT / 'validation'
CASES = {
    'naca0012_ladson': ('naca0012', 'Ladson, NASA TM-4074, Table I, p.6'),
    's809_somers': ('s809', 'Somers, NREL/SR-440-6918, Figs.17(c)(i), 15(c)'),
    'nlr7301': ('nlr7301', 'van den Berg/Oskam; Ranneberg digitization'),
    'nlr7301_full_polar': ('nlr7301', 'Optional full NLR7301 polar including lift loss'),
    'nlr7301_mses': ('nlr7301_mses', 'Published MSES: Neto et al., ICAS2006, Figure6, NLR7301 2.6% gap'),
    'naca0012_mses_drag': ('naca0012_mses_drag', 'Published MSES: Lima e Silva et al., COBEM2005, Figure3, NACA0012 drag'),
}
DEFAULT_CASES = ['naca0012_ladson','s809_somers','nlr7301','nlr7301_mses','naca0012_mses_drag']


def number(value):
    try:
        value = float(value)
        return value if math.isfinite(value) else None
    except (ValueError, TypeError):
        return None


def reference_info(cfg):
    info=dict(cfg.get('benchmark',{}))
    info.setdefault('reference_kind','experiment')
    if info['reference_kind'] not in ('experiment','published_solver'):
        raise ValueError('Unknown reference_kind; use experiment or published_solver')
    if info['reference_kind']=='published_solver' and not info.get('reference_solver'):
        raise ValueError('A published-solver reference must identify reference_solver')
    info.setdefault('reference_label',f"Published {info['reference_solver']}"
                    if info['reference_kind']=='published_solver' else 'Experiment')
    return info


def read_csv(path):
    with Path(path).open(newline='') as f:
        return list(csv.DictReader(f))


def write_csv(path, rows, fields):
    with Path(path).open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def verify_data():
    """Reject silent edits to the distributed reference data and geometries."""
    manifest = json.loads((DATA / 'sources.json').read_text())
    for filename, expected in manifest['data_sha256'].items():
        actual = hashlib.sha256((DATA / filename).read_bytes()).hexdigest()
        if actual != expected:
            raise ValueError(f'Reference-data hash mismatch: {filename}')
    return manifest


def at_angle(rows, alpha):
    """Exact reference-angle match (floating-point tolerance only), no interpolation."""
    matches = [r for r in rows if abs(float(r['alpha_deg']) - alpha) < 1e-7]
    if len(matches) > 1:
        raise ValueError(f'Duplicate result at alpha={alpha}')
    return matches[0] if matches else None


def error_stats(rows):
    values = [r['error'] for r in rows if r['error'] is not None]
    return {
        'reference_values_requested': len(rows),
        'values_compared': len(values),
        'unavailable': len(rows) - len(values),
        'bias': float(np.mean(values)) if values else None,
        'mae': float(np.mean(np.abs(values))) if values else None,
        'rmse': float(np.sqrt(np.mean(np.square(values)))) if values else None,
        'max_abs': float(np.max(np.abs(values))) if values else None,
    }


def compare_polar(reference, predictions, requested):
    comparisons = []
    for ref in reference:
        alpha = float(ref['alpha_deg'])
        if not any(abs(alpha-a) < 1e-7 for a in requested):
            continue
        pred = at_angle(predictions, alpha)
        converged = pred is not None and pred.get('converged') == '1'
        for quantity in ('CL', 'CD', 'CM'):
            observed = number(ref.get(quantity))
            if observed is None:
                continue  # A blank experimental CD is not a zero-drag observation.
            computed = number(pred.get(quantity)) if converged else None
            low, high = number(ref.get(quantity+'_low')), number(ref.get(quantity+'_high'))
            error = computed-observed if computed is not None else None
            outside = (max(low-computed, computed-high, 0.)
                       if computed is not None and low is not None and high is not None else None)
            comparisons.append(dict(alpha_deg=alpha, quantity=quantity,
                reference=observed, predicted=computed, error=error,
                reference_low=low, reference_high=high, distance_outside_interval=outside,
                status='compared' if computed is not None else 'unavailable'))
    return comparisons


def compare_transition(reference, predictions, requested, diagnostics=None):
    rows = []
    for ref in reference:
        alpha = float(ref['alpha_deg'])
        if not any(abs(alpha-a) < 1e-7 for a in requested):
            continue
        subset = [r for r in predictions if r['element'] == ref['element'] and r['side'] == ref['side']]
        match = at_angle(subset, alpha)
        computed = number(match.get('transition_x_over_local_chord')) if match else None
        onset = computed
        observable = ref.get('observable','transition')
        if computed is not None and observable == 'attached_turbulence' and diagnostics is not None:
            diagnostic = next((d for a,d in diagnostics if abs(a-alpha)<1e-7), {})
            bl = diagnostic.get('boundary_layers',{}).get(ref['element'],{}).get(ref['side'],{})
            tr = bl.get('transition')
            if tr:
                st = tr['s_over_cref']
                # A microphone sees attached turbulent flow, not e^N onset inside a bubble.
                for region in bl.get('separation_regions',[]):
                    if region['start']['s_over_cref'] <= st <= region['end']['s_over_cref']:
                        computed = (None if region['reaches_trailing_edge'] else
                                    region['end']['x_over_local_chord'])
                        break
        low, high = float(ref['x_low']), float(ref['x_high'])
        allowance = float(ref['reading_allowance'])
        outside = max(low-computed, computed-high, 0.) if computed is not None else None
        expanded = max(low-allowance-computed, computed-high-allowance, 0.) if computed is not None else None
        rows.append(dict(alpha_deg=alpha, element=ref['element'], side=ref['side'],
            alpha_tunnel_deg=number(ref.get('alpha_tunnel_deg')), observable=observable,
            model_transition_onset=onset,
            reference_low=low, reference_high=high, reading_allowance=allowance,
            predicted=computed, distance_outside_interval=outside,
            distance_outside_with_reading_allowance=expanded,
            status='compared' if computed is not None else 'no_transition_result'))
    return rows


def sampled_peak(reference, predictions, requested, peak_alpha=None):
    """Compare a measured peak only if its adjacent reference angles were attempted."""
    refs = sorted([r for r in reference if number(r.get('CL')) is not None],
                  key=lambda r: float(r['alpha_deg']))
    if not refs:
        return {'reference_peak_selection': 'unavailable: reference contains no CL',
                'reference_sampled_maximum': None, 'reference_peak_bracket_covered': False,
                'predicted_CL_at_reference_peak': None, 'error_at_reference_peak': None,
                'predicted_drop_across_reference_peak': None,
                'reference_drop_across_reference_peak': None}
    i = (max(range(len(refs)), key=lambda k: float(refs[k]['CL'])) if peak_alpha is None else
         min(range(len(refs)), key=lambda k: abs(float(refs[k]['alpha_deg'])-peak_alpha)))
    peak = refs[i]
    output = {'reference_peak_selection': 'global sampled maximum' if peak_alpha is None else
                'author-defined usable maximum; later separated-flow values are suspect',
              'reference_sampled_maximum': {'alpha_deg': float(peak['alpha_deg']),
                'CL': float(peak['CL'])}, 'reference_peak_bracket_covered': False,
                'predicted_CL_at_reference_peak': None, 'error_at_reference_peak': None,
                'predicted_drop_across_reference_peak': None,
                'reference_drop_across_reference_peak': None}
    if 0 < i < len(refs)-1:
        bracket = refs[i-1:i+2]
        output['reference_peak_bracket_deg'] = [float(r['alpha_deg']) for r in bracket]
        covered = all(any(abs(float(r['alpha_deg'])-a) < 1e-7 for a in requested) for r in bracket)
        output['reference_peak_bracket_covered'] = covered
        if covered:
            preds = [at_angle(predictions, float(r['alpha_deg'])) for r in bracket]
            if all(r is not None and r['converged'] == '1' and number(r.get('CL')) is not None for r in preds):
                output['predicted_CL_at_reference_peak'] = float(preds[1]['CL'])
                output['error_at_reference_peak'] = float(preds[1]['CL'])-float(peak['CL'])
                output['predicted_drop_across_reference_peak'] = float(preds[1]['CL'])-float(preds[2]['CL'])
                output['reference_drop_across_reference_peak'] = float(peak['CL'])-float(bracket[2]['CL'])
    return output


def check_metadata(folder, cfg):
    meta = json.loads((folder/'metadata.json').read_text())
    if meta.get('backend') != 'viscous':
        raise ValueError(f'{folder}: validation requires the viscous backend')
    actual = meta['case']
    for key in ('reynolds', 'mach', 'ncrit', 'iterations', 'alpha_deg', 'elements', 'reference_chord'):
        if actual.get(key) != cfg.get(key):
            raise ValueError(f'{folder}: metadata mismatch for {key}; use matching case/settings')
    return meta


def plot_comparison(folder, reference, predictions, transition, title, info=None, context=()):
    info=info or reference_info({})
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    quantities = [q for q in ('CL', 'CD', 'CM') if any(number(r.get(q)) is not None for r in reference)]
    drag_only = quantities == ['CD']
    if drag_only:
        quantities = ['CD', 'CD']  # Full range plus a readable low-incidence detail.
    fig, axes = plt.subplots(1, len(quantities), figsize=(5*len(quantities), 4), squeeze=False)
    for ax, q in zip(axes[0], quantities):
        rr = [r for r in reference if number(r.get(q)) is not None]
        x, y = [float(r['alpha_deg']) for r in rr], [float(r[q]) for r in rr]
        regular=[r for r in rr if not r.get('interpretation','').startswith('extensive')]
        stress=[r for r in rr if r.get('interpretation','').startswith('extensive')]
        ax.plot([float(r['alpha_deg']) for r in regular],[float(r[q]) for r in regular],
                's' if info['reference_kind']=='published_solver' else 'o', ms=4, mfc='none',
                color='#25813b' if info['reference_kind']=='published_solver' else 'black',
                label=info['reference_label'])
        contextual=[r for r in context if number(r.get(q)) is not None]
        if contextual:
            ax.plot([float(r['alpha_deg']) for r in contextual],[float(r[q]) for r in contextual],
                    'x',ms=4,color='gray',label='Experiment (context only)')
        if stress:
            ax.plot([float(r['alpha_deg']) for r in stress],[float(r[q]) for r in stress],
                    '+', ms=6, color='#9b653d', label='Separated-flow data (qualified)')
        if all(number(r.get(q+'_low')) is not None for r in rr):
            ax.errorbar(x, y, yerr=[[float(r[q])-float(r[q+'_low']) for r in rr],
                                   [float(r[q+'_high'])-float(r[q]) for r in rr]],
                        fmt='none', ecolor='gray', lw=1,
                        label='Graphical reading allowance' if info['reference_kind']=='published_solver' else 'Reading / branch range')
        ax.plot([float(r['alpha_deg']) for r in predictions],
                [number(r.get(q)) if r['converged']=='1' and number(r.get(q)) is not None else np.nan for r in predictions],
                '.-', color='#1774af', label='MultiFoil')
        flagged=[r for r in predictions if r['converged']=='1' and number(r.get(q)) is not None and
                 any(abs(float(r['alpha_deg'])-a)<1e-7 for a in info.get('model_flagged_angles_deg',[]))]
        if flagged:
            ax.plot([float(r['alpha_deg']) for r in flagged],[float(r[q]) for r in flagged],
                    'x',color='#b22b35',ms=7,label='Model: local supersonic flag')
        for r in predictions:
            if r['converged'] != '1':
                ax.plot(float(r['alpha_deg']), .03, 'x', color='red', transform=ax.get_xaxis_transform())
        ax.set(xlabel='Angle of attack [deg]', ylabel=q)
        if drag_only and ax is axes[0,1]:
            low_refs=[r for r in rr if abs(float(r['alpha_deg']))<=10]
            low_preds=[r for r in predictions if abs(float(r['alpha_deg']))<=10 and
                       r['converged']=='1' and number(r.get(q)) is not None]
            low_x=[float(r['alpha_deg']) for r in low_refs+low_preds]
            low_y=[float(r.get(q+'_high') or r[q]) for r in low_refs]+[float(r[q]) for r in low_preds]
            if low_x and low_y:
                ax.set(xlim=(min(low_x)-.5,max(low_x)+.5),ylim=(0,1.15*max(low_y)),
                       title='Detail: |alpha| <= 10 degrees')
        ax.grid(alpha=.2)
        ax.legend(fontsize=8)
    suffix=('MultiFoil versus '+info['reference_label'] if info['reference_kind']=='published_solver'
            else 'measured coefficients versus model')
    fig.suptitle(title+' — '+suffix, fontsize=11)
    fig.tight_layout()
    fig.savefig(folder/'comparison.png', dpi=160, bbox_inches='tight')
    plt.close(fig)
    if transition:
        fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
        for ax, side in zip(axes, ('upper', 'lower')):
            rr = [r for r in transition if r['side'] == side]
            for r in rr:
                ax.plot([r['reference_low']-r['reading_allowance'], r['reference_high']+r['reading_allowance']],
                        [r['alpha_deg']]*2, color='lightgray', lw=5)
                ax.plot([r['reference_low'], r['reference_high']], [r['alpha_deg']]*2, color='black', lw=2)
            ax.plot([r['predicted'] if r['predicted'] is not None else np.nan for r in rr],
                    [r['alpha_deg'] for r in rr], 'o-', ms=4, label='MultiFoil')
            missing=[r for r in rr if r['predicted'] is None]
            if missing:
                ax.plot([.98]*len(missing),[r['alpha_deg'] for r in missing], 'x', color='red',
                        transform=ax.get_yaxis_transform(),label='Unavailable')
            ax.plot([], [], color='black', lw=2, label='Measured range')
            ax.plot([], [], color='lightgray', lw=5, label='Reading allowance')
            ax.set(title=side.capitalize(), xlabel='Transition x/c', xlim=(-.03, 1.03))
            ax.grid(alpha=.2)
            ax.legend(fontsize=8)
        axes[0].set_ylabel('Angle of attack [deg]')
        fig.suptitle(title+' — transition / attached turbulence; inferred angle correction', fontsize=11)
        fig.tight_layout()
        fig.savefig(folder/'transition_comparison.png', dpi=160, bbox_inches='tight')
        plt.close(fig)


def compare_case(name, folder, cfg):
    info=reference_info(cfg)
    metadata=check_metadata(folder, cfg)
    predictions = read_csv(folder/'polar.csv')
    requested = cfg['alpha_deg']
    if any(not any(abs(float(r['alpha_deg'])-a)<1e-7 for a in requested) for r in predictions):
        raise ValueError(f'{folder}: unexpected angle in results')
    directory = DATA/'data'/CASES[name][0]
    reference = read_csv(directory/'polar.csv')
    if any(r.get('reference_kind',info['reference_kind'])!=info['reference_kind'] for r in reference):
        raise ValueError('Reference CSV kind does not match benchmark classification')
    rows = compare_polar(reference, predictions, requested)
    for row in rows:
        row['reference_kind']=info['reference_kind']
        row['reference_label']=info['reference_label']
    write_csv(folder/'force_errors.csv', rows, ['alpha_deg','reference_kind','reference_label','quantity','reference','predicted','error',
        'reference_low','reference_high','distance_outside_interval','status'])
    metrics = {}
    for quantity in ('CL','CD','CM'):
        qr = [r for r in rows if r['quantity'] == quantity]
        if not qr:
            continue
        metrics[quantity] = {region: error_stats([r for r in qr if select(r)]) for region,select in [
            ('all',lambda r:True),('abs_alpha_le_10',lambda r:abs(r['alpha_deg'])<=10),
            ('abs_alpha_gt_10',lambda r:abs(r['alpha_deg'])>10)]}
    transition = []
    if (directory/'transition.csv').exists():
        transition = compare_transition(read_csv(directory/'transition.csv'),
            read_csv(folder/'transition_separation.csv'), requested,
            list(zip([float(r['alpha_deg']) for r in predictions],metadata['diagnostics'])))
        write_csv(folder/'transition_errors.csv', transition, ['alpha_deg','alpha_tunnel_deg','element','side',
            'observable','model_transition_onset','reference_low',
            'reference_high','reading_allowance','predicted','distance_outside_interval',
            'distance_outside_with_reading_allowance','status'])
    conv = [a for a in requested if (at_angle(predictions,a) or {}).get('converged')=='1']
    result = {'case':name,'reference':CASES[name][1],
        'reference_kind':info['reference_kind'],'reference_label':info['reference_label'],
        'reference_qualifications':info,
        'conditions':{k:cfg[k] for k in ('reynolds','mach','ncrit')},
        'requested_flow_points':len(requested),'converged_flow_points':len(conv),
        'failed_or_missing_angles_deg':[a for a in requested if a not in conv],
        'local_supersonic_angles_deg':[float(p['alpha_deg']) for p,d in zip(predictions,metadata['diagnostics'])
            if any(d.get('local_supersonic_detected',{}).values())],
        'metrics':metrics, 'transition_values_requested':len(transition),
        'transition_values_compared':sum(r['predicted'] is not None for r in transition),
        'transition_outside_with_reading_allowance':sum(
            r['distance_outside_with_reading_allowance'] is not None and
            r['distance_outside_with_reading_allowance']>0 for r in transition),
        'stall_comparison':sampled_peak(reference,predictions,requested,9.22 if name=='s809_somers' else None),
        'solver_stall_diagnostic':json.loads((folder/'stall.json').read_text()),
        'status':'comparison_only_no_accuracy_pass_claim',
        'notes':['Errors are prediction minus reference; no interpolation of MultiFoil solutions or fitting.',
                 'Statistics use available values; explicit coverage counts include missing solves.',
                 'Low/high angle bins are numerical bins, not declarations of attached/separated flow.',
                 'No experimental uncertainty was invented for tabulated data. Digitized ranges are reading allowances.']}
    result['notes'].extend(info.get('notes',[]))
    (folder/'comparison.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    context=read_csv(directory/info['context_data']) if info.get('context_data') else []
    plot_info={**info,'model_flagged_angles_deg':result['local_supersonic_angles_deg']}
    plot_comparison(folder,reference,predictions,transition,name,plot_info,context)
    return result


def report(results, out):
    lines = ['# Benchmark comparison results','',
        'Experimental references assess measurement agreement; published solver references assess',
        'code-to-code agreement. Convergence is not an accuracy pass. See',
        '`validation/README.md` and `sources.json` for conditions and data qualifications.','',
        '| Case | Reference | Converged / requested | Quantity | Compared / reference values | Bias | RMSE | Max absolute error |',
        '|---|---|---:|---|---:|---:|---:|---:|']
    flat = []
    fmt = lambda x:'unavailable' if x is None else f'{x:.6g}'
    for result in results:
        for q, stat in result['metrics'].items():
            s=stat['all']
            lines.append(f"| {result['case']} | {result['reference_label']} | {result['converged_flow_points']} / {result['requested_flow_points']} | {q} | {s['values_compared']} / {s['reference_values_requested']} | {fmt(s['bias'])} | {fmt(s['rmse'])} | {fmt(s['max_abs'])} |")
            flat.append(dict(case=result['case'],reference_kind=result['reference_kind'],quantity=q,**s))
    lines += ['', 'One drag count = 0.0001 CD. CL and CM errors are absolute coefficients.',
        'S809 CL errors use digitized interval centers; these are less precise than the NACA table.',
        'MSES errors use samples from the published graphical curve, not original solver files.',
        'An empty reference CD entry is omitted, never replaced by zero.','']
    for r in results:
        lines += [f"## {r['case']}",'',f"Reference: {r['reference_label']} ({r['reference_kind']}).",
                  f"Conditions: Re={r['conditions']['reynolds']:g}, Mach={r['conditions']['mach']:g}, Ncrit={r['conditions']['ncrit']:g}.",
                  f"Failed or missing angles: {r['failed_or_missing_angles_deg']}.",
                  f"Sampled lift-peak diagnostic: **{r['solver_stall_diagnostic']['status']}**."]
        if r['local_supersonic_angles_deg']:
            lines += [f"Local supersonic flow flagged at {r['local_supersonic_angles_deg']} degrees: the compressibility correction cannot resolve shocks. Convergence does not establish model validity there."]
        st = r['stall_comparison']
        if st['reference_sampled_maximum'] is None:
            lines += ['No reference CL was supplied; no reference lift-peak comparison is available.']
        else:
            lines += [f"Reference peak: {st['reference_sampled_maximum']} ({st['reference_peak_selection']}); its three-point bracket was {'covered' if st['reference_peak_bracket_covered'] else 'not covered'} by the requested angles."]
        if 'CD' in r['metrics']:
            drag=r['metrics']['CD']['all']
            counts=lambda v:fmt(v*10000 if v is not None else None)
            lines += [f"Drag error in counts: bias {counts(drag['bias'])}, RMSE {counts(drag['rmse'])}, maximum absolute {counts(drag['max_abs'])}."]
        if st['predicted_CL_at_reference_peak'] is not None:
            lines += [f"Predicted CL at that angle: {fmt(st['predicted_CL_at_reference_peak'])}; error {fmt(st['error_at_reference_peak'])}.",
                      f"Lift drop to next reference angle: reference {fmt(st['reference_drop_across_reference_peak'])}, predicted {fmt(st['predicted_drop_across_reference_peak'])} (negative means lift still rises)."]
        if r['reference_kind']=='published_solver':
            lines += ['',*r['reference_qualifications'].get('notes',[]),'']
        if r['transition_values_requested']:
            lines += [f"Transition: {r['transition_values_compared']} / {r['transition_values_requested']} values compared; {r['transition_outside_with_reading_allowance']} outside the measured range plus digitization allowance.",
                      f"![Transition comparison]({r['case']}/transition_comparison.png)"]
        lines += ['',f"![Force comparison]({r['case']}/comparison.png)",'']
    (out/'REPORT.md').write_text('\n'.join(lines)+'\n')
    (out/'summary.json').write_text(json.dumps(results,indent=2,allow_nan=False)+'\n')
    write_csv(out/'summary.csv',flat,['case','reference_kind','quantity','reference_values_requested','values_compared',
                                     'unavailable','bias','mae','rmse','max_abs'])


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', choices=list(CASES), nargs='+', default=DEFAULT_CASES)
    parser.add_argument('--list', action='store_true')
    parser.add_argument('--out', type=Path, default=Path('validation_run'))
    parser.add_argument('--compare-only', action='store_true', help='Compare existing outputs with matching metadata')
    parser.add_argument('--ncrit', type=float, help='Explicit sensitivity study; recorded in outputs')
    parser.add_argument('--panel-scale', type=float, default=1., help='Scale every element panel count for refinement')
    args=parser.parse_args(argv)
    if args.list:
        for name, (_,description) in CASES.items(): print(name+': '+description)
        return 0
    try:
        verify_data()
        if not math.isfinite(args.panel_scale) or args.panel_scale <= 0:
            raise ValueError('--panel-scale must be finite and positive')
        out=args.out.resolve()
        results=[]
        for name in dict.fromkeys(args.case):
            from multifoil import load_case, save_results, plot_results
            cfg,elements,ref=load_case(DATA/'cases'/f'{name}.json')
            defaults={'reynolds':1e6,'mach':0.,'ncrit':9.,'iterations':80,
                'amplification_jacobian':'corrected_af',
                'tolerance':1e-7,'wake_length':1.,'minimum_te_gap':1e-5,
                'continuation_step':1.,'min_step':.25,'start_alpha':None}
            for key,value in defaults.items():cfg.setdefault(key,value)
            if args.ncrit is not None: cfg['ncrit']=args.ncrit
            if args.panel_scale != 1.:
                # Rebuild from the original coordinates, not an already coarse mesh.
                import tempfile
                for e,(_,x) in zip(cfg['elements'],elements):
                    e['panels']=max(40,int(round((len(x)-1)*args.panel_scale)))
                absolute=json.loads(json.dumps(cfg))
                for e in absolute['elements']:
                    if 'file' in e:e['file']=str((DATA/'cases'/e['file']).resolve())
                with tempfile.TemporaryDirectory(prefix='multifoil_validation_') as tmp:
                    temporary=Path(tmp)/'case.json'
                    temporary.write_text(json.dumps(absolute))
                    _,elements,ref=load_case(temporary)
            folder=out/name
            if not args.compare_only:
                if folder.exists() and any(folder.iterdir()):
                    raise ValueError(f'{folder} is not empty; choose another --out or --compare-only')
                from viscous import CoupledSolver
                solver=CoupledSolver(elements,ref,**{k:cfg[k] for k in defaults})
                solutions=[]
                for alpha in cfg['alpha_deg']:
                    r=solver.solve(alpha);solutions.append(r)
                    print(f'{name}: alpha={alpha:.6g}, CL={r.cl:.6g}, CD={r.cd:.6g}, converged={r.converged}',flush=True)
                    save_results(folder,solutions,solver.elements,cfg,'viscous')
                plot_results(folder,solutions,solver.elements,solver)
            result=compare_case(name,folder,cfg)
            results.append(result)
            report(results,out)
        complete=all(r['converged_flow_points']==r['requested_flow_points'] and
                     r['transition_values_compared']==r['transition_values_requested'] and
                     all(s['all']['unavailable']==0 for s in r['metrics'].values()) for r in results)
        print(f'Saved {out / "REPORT.md"}; comparison complete={complete}; no accuracy pass implied.')
        return 0 if complete else 2
    except (ValueError, OSError, KeyError) as exc:
        print(f'Validation input error: {exc}',file=sys.stderr)
        return 1


if __name__=='__main__':
    raise SystemExit(main())
