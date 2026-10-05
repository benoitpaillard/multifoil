#!/usr/bin/env python3
"""V3: inviscid cross-formulation check on identical S414 polygons; linear analysis only.

Production streamfunction/linear-vortex operator against the existing Hess–Smith
source/vortex backend, after qualifying both on controls. No viscous solve.
Protocol: validation/S414_FORMULATION_CHECK_PROTOCOL.md.
"""
import os
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
import argparse
import importlib.metadata
from pathlib import Path
import platform
import sys
import time
from types import SimpleNamespace
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from report_design_initialization import Inputs, sha, write, csv_write
from report_s702_drag_mechanism import require
from multifoil import HessSmithSolver, clean_contour, signed_area, transform
from viscous import CoupledSolver, new_model, te_info, vortex_velocity_matrix
from paneling import curvature_paneling
import run_s414_transition_boundary as deck

OUT = ROOT/'validation/results/s414_formulation_check_20261004'
PROTOCOL = ROOT/'validation/S414_FORMULATION_CHECK_PROTOCOL.md'
LOWER = ROOT/'validation/results/s414_lower_bucket_20261003'
LOWER_RECEIPT = ROOT/'validation/results/s414_lower_bucket_launch_20261003/retrieved_files_sha256.json'
FORE = ROOT/deck.ARCHIVE
TAPS = ROOT/'validation/results/s414_lower_bucket_analysis_20261003/pressure_at_taps.csv'
SCHEDULE = ((340, 240), (480, 340), (680, 480))  # fixed before execution
MACH, NOSE, LINEAR = .1, .06, 1e-8
CONTOUR = dict(center=(.5, 0.), radius=3., points=720)
SURFACES = ('fore_upper', 'fore_lower', 'aft_upper', 'aft_lower')
KT_PANELS = (120, 240, 480)


def cp_kt(speed, mach=MACH):
    """Karman–Tsien pressure from uncorrected speed / Vinf; identical for both methods."""
    cp0 = 1-np.asarray(speed, float)**2; beta = np.sqrt(1-mach*mach)
    return cp0/(beta+mach*mach*cp0/(2*(1+beta)))


def production_operator(polygons):
    """Unchanged production inviscid assembly on prescribed clockwise nodes; no repaneling."""
    solver = object.__new__(CoupledSolver); solver.models = []; solver.fs = []; start = 0
    for name, x in polygons:
        x = np.ascontiguousarray(x, dtype=float)
        solver.models.append(SimpleNamespace(geom=SimpleNamespace(name=name), isol=SimpleNamespace(),
                                             foil=SimpleNamespace(x=x, N=x.shape[1])))
        solver.fs.append(slice(start, start+x.shape[1])); start += x.shape[1]
    solver.nfoil = start; solver._build_inviscid()
    return solver


def frame(xy):
    """Chord frame of the solver input, exactly as viscous.new_model builds it."""
    X = np.ascontiguousarray(clean_contour(xy)[::-1].T)
    te = .5*(X[:, 0]+X[:, -1]); ile = int(np.argmax(np.linalg.norm(X-te[:, None], axis=0)))
    le = X[:, ile].copy(); chord = np.linalg.norm(te-le); axis = (te-le)/chord
    require(te_info(X)[1] >= 1e-5*chord, 'Finite trailing edge expected; no gap regularization may act here')
    basis = np.column_stack((axis, [-axis[1], axis[0]]))
    return le, chord, basis, basis.T@(X-le[:, None])/chord


def sample(xy, count):
    """`count` production-spaced nodes on the spline of the unchanged baseline input."""
    le, chord, basis, local = frame(xy)
    nodes, _, _, recovered = curvature_paneling(local, count, 2., .1)
    return le[:, None]+chord*(basis@nodes), bool(recovered)


def unit(alpha):
    return np.array([np.cos(np.deg2rad(alpha)), np.sin(np.deg2rad(alpha))])


def contour_points():
    t = 2*np.pi*np.arange(CONTOUR['points'])/CONTOUR['points']
    circle = np.array([np.cos(t), np.sin(t)])
    return np.asarray(CONTOUR['center'])[:, None]+CONTOUR['radius']*circle, circle, 2*np.pi*CONTOUR['radius']/CONTOUR['points']


class Pair:
    """Both formulations on one identical set of clockwise polygons."""
    def __init__(self, polygons, far=True):
        self.polygons = [(n, np.ascontiguousarray(x, dtype=float)) for n, x in polygons]
        ccw = [(n, x.T[::-1].copy()) for n, x in self.polygons]
        require(all(signed_area(xy) > 0 for _, xy in ccw), 'Hess–Smith input must be counter-clockwise')
        self.production = production_operator(self.polygons); self.hess = HessSmithSolver(ccw)
        self.sharp = [abs(te_info(x)[1]) < 1e-12 for _, x in self.polygons]
        self.far = None
        if far:
            points, normal, ds = contour_points()
            self.far = (points, normal, ds, [vortex_velocity_matrix(x, points) for _, x in self.polygons])

    def solve(self, alpha):
        u = unit(alpha); p = self.production; h = self.hess.solve(alpha); out = dict(alpha=alpha, elements={})
        solution = np.linalg.solve(p.aic, p.rhs@u)
        out['production_linear_residual'] = float(np.max(abs(p.aic@solution-p.rhs@u)))
        out['hess_linear_residual'] = float(h.residual)
        out['hess_kutta_residual'] = float(h.diagnostics['kutta_residual'])
        out['hess_source_flux'] = dict(zip([n for n, _ in self.polygons], h.diagnostics['source_flux_by_element']))
        out['hess_CL_circulation'] = float(h.diagnostics['cl_circulation_total'])
        kutta = 0.; totals = dict(production=np.zeros(2), hess=np.zeros(2), production_nodal=np.zeros(2))
        for k, ((name, x), sl) in enumerate(zip(self.polygons, p.fs)):
            gamma = p.gref[sl]@u; n = x.shape[1]; kutta = max(kutta, abs(gamma[0]+gamma[-1]))
            require(float(np.max(abs(gamma-solution[sl]))) < 1e-9, 'Production basis differs from a direct solve')
            vt = h.surfaces[k]['ue_signed']
            # Hess–Smith panels run counter-clockwise; reverse order and sign to the clockwise convention.
            hs = -vt[:n-1][::-1]; mid = .5*(gamma[:-1]+gamma[1:])
            d = np.diff(x, axis=1); ds = np.linalg.norm(d, axis=0); s = np.r_[0., np.cumsum(ds)]
            base = dict(production=.5*(2-gamma[0]**2-gamma[-1]**2), hess=None if self.sharp[k] else float(1-vt[n-1]**2))
            gap = x[:, 0]-x[:, -1]  # clockwise closing segment: upper TE -> lower TE
            force = {}
            for label, speed in (('production', mid), ('hess', hs)):
                cp = 1-speed**2; f = np.array([-np.sum(cp*d[0]), np.sum(cp*d[1])])
                if not self.sharp[k]: f += base[label]*np.array([-gap[0], gap[1]])
                force[label] = f; totals[label] += f
            cpn = 1-gamma**2; closed = np.c_[x, x[:, :1]]; cpc = np.r_[cpn, cpn[0]]; dc = np.diff(closed, axis=1)
            totals['production_nodal'] += np.array([-np.sum(.5*(cpc[:-1]+cpc[1:])*dc[0]), np.sum(.5*(cpc[:-1]+cpc[1:])*dc[1])])
            out['elements'][name] = dict(x=x, s_nodes=s, s_mid=.5*(s[:-1]+s[1:]), xm=.5*(x[:, :-1]+x[:, 1:]),
                d=d, gamma=gamma, production=mid, hess=hs, force=force, base_cp=base)
        out['production_kutta_residual'] = float(kutta)
        lift = lambda f: float(f[0]*u[0]-f[1]*u[1])
        out['CL_pressure'] = {k: lift(v) for k, v in totals.items()}
        if self.far:
            points, normal, ds, matrices = self.far; tangent = np.array([-normal[1], normal[0]])
            vp = np.tile(u, (points.shape[1], 1))
            for matrix, sl in zip(matrices, p.fs): vp += np.einsum('ijk,k->ij', matrix, p.gref[sl]@u)
            vh = self.hess.velocity(points.T, alpha)
            for label, v in (('production', vp), ('hess', vh)):
                out[label+'_far_CL'] = float(-2*np.sum(v[:, 0]*tangent[0]+v[:, 1]*tangent[1])*ds)
                out[label+'_far_flux'] = float(np.sum(v[:, 0]*normal[0]+v[:, 1]*normal[1])*ds)
        return out


def crossing(s, u):
    """Stagnation arclength: sign change of clockwise speed nearest the leading edge."""
    k = np.flatnonzero((u[:-1] < 0) & (u[1:] >= 0))
    require(len(k) >= 1, 'No stagnation crossing')
    k = k[np.argmin(abs(k-len(u)/2))]
    return float(s[k]-u[k]*(s[k+1]-s[k])/(u[k+1]-u[k]))


def surface_metrics(element, label):
    """Loads, nose minimum and stagnation from midpoint samples of one method."""
    e = element; x = e['x']; u = e[label]; cp = cp_kt(u); n = x.shape[1]
    ile = int(np.argmin(x[0])); xm = e['xm']
    rows = {}
    for side, panels in (('lower', slice(0, ile)), ('upper', slice(ile, n-1))):
        cpi = 1-u[panels]**2
        nose = np.flatnonzero(xm[0, panels] < x[0, ile]+NOSE)
        k = nose[np.argmin(cp[panels][nose])]
        rows[side] = dict(CN_incompressible=float(-np.sum(cpi*e['d'][0, panels])),
            nose_min_Cp=float(cp[panels][k]), nose_min_Cp_x=float(xm[0, panels][k]),
            nose_max_speed=float(abs(u[panels][k])))
    s0 = crossing(e['s_mid'], u)
    rows['stagnation'] = dict(s=s0, x=float(np.interp(s0, e['s_nodes'], x[0])), y=float(np.interp(s0, e['s_nodes'], x[1])))
    return rows


def at_positions(element, label, fractions):
    return np.interp(fractions*element['s_nodes'][-1], element['s_mid'], element[label])


def tap_values(solution, label, taps):
    """Corrected Cp of one method at supported tap x, interpolated on its own surface branch."""
    values = []
    for surface, x in taps:
        e = solution['elements'][surface.split('_')[0]]; ile = int(np.argmin(e['x'][0]))
        panels = slice(0, ile) if surface.endswith('lower') else slice(ile, e['x'].shape[1]-1)
        xs = e['xm'][0, panels]; cp = cp_kt(e[label][panels]); order = np.argsort(xs)
        # No extrapolation: a tap beyond the first or last panel midpoint is counted, not scored.
        values.append(float(np.interp(x, xs[order], cp[order])) if xs[order][0] <= x <= xs[order][-1] else np.nan)
    return np.array(values)


def karman_trefftz(panels, tau_deg=10., eps=.1, kappa=.1):
    """Counter-clockwise sharp-TE contour (upper TE first), exact speed and lift slope data."""
    n = 2-np.deg2rad(tau_deg)/np.pi; mu = complex(-eps, kappa); a = abs(1-mu); beta = np.arcsin(kappa/a)
    def mapped(theta):
        zeta = mu+a*np.exp(1j*theta); r = (zeta-1)/(zeta+1); rn = r**n
        return n*(1+rn)/(1-rn), zeta, r, rn
    def speed(theta, alpha):
        _, zeta, r, rn = mapped(theta); al = np.deg2rad(alpha)
        w = np.exp(-1j*al)-a*a*np.exp(1j*al)/(zeta-mu)**2+2j*a*np.sin(al+beta)/(zeta-mu)
        return abs(w)/abs(4*n*n*r**(n-1)/((1-rn)**2*(zeta+1)**2))
    theta = -beta+2*np.pi*np.arange(panels+1)/panels
    z = mapped(theta)[0]; z[-1] = z[0]; lo, hi = z.real.min(), z.real.max(); chord = hi-lo
    xy = np.column_stack(((z.real-lo)/chord, z.imag/chord))
    return xy, theta, speed, lambda alpha: 8*np.pi*a*np.sin(np.deg2rad(alpha)+beta)/chord


def naca0012_open(panels):
    x = .5*(1-np.cos(np.linspace(0, np.pi, panels//2+1)))
    y = .6*(.2969*np.sqrt(x)-.126*x-.3516*x**2+.2843*x**3-.1015*x**4)
    return np.concatenate((np.column_stack((x, y))[::-1], np.column_stack((x, -y))[1:]))


def clockwise(xy):
    return np.ascontiguousarray(xy[::-1].T)


def controls():
    rows, ok = [], True
    previous = {}
    for panels in KT_PANELS:
        xy, theta, speed, lift = karman_trefftz(panels); alpha = 4.
        s = Pair([('kt', clockwise(xy))]).solve(alpha); e = s['elements']['kt']
        mid_theta = (.5*(theta[:-1]+theta[1:]))[::-1]; node_theta = theta[::-1]
        keep = e['xm'][0] < .98; nodes = (e['x'][0] < .98)
        errors = dict(hess=abs(e['hess'])-speed(mid_theta, alpha),
                      production=abs(e['gamma'][nodes])-speed(node_theta[nodes], alpha))
        for label in ('production', 'hess'):
            err = errors[label][keep] if label == 'hess' else errors[label]
            rms = float(np.sqrt(np.mean(err**2))); cl = s[label+'_far_CL']
            row = dict(control='karman_trefftz', method=label, panels=panels, alpha_deg=alpha,
                rms_speed_error=rms, max_speed_error=float(np.max(abs(err))), CL_exact=lift(alpha),
                CL_circulation=cl, CL_relative_error=abs(cl/lift(alpha)-1), CL_pressure=s['CL_pressure'][label],
                far_flux=s[label+'_far_flux'], linear_residual=s[label+'_linear_residual'],
                kutta_residual=s[label+'_kutta_residual'])
            if label in previous:
                row['monotone'] = bool(rms < previous[label][0] and row['CL_relative_error'] < previous[label][1])
                ok = ok and row['monotone']
            previous[label] = (rms, row['CL_relative_error']); rows.append(row)
            ok = ok and row['linear_residual'] < LINEAR and row['kutta_residual'] < LINEAR
            if panels == KT_PANELS[-1]: ok = ok and row['CL_relative_error'] < .005 and rms < .005
    foil = clockwise(naca0012_open(200)); pair = Pair([('naca', foil)])
    zero, plus, minus = (pair.solve(a) for a in (0., 4., -4.))
    for label in ('production', 'hess'):
        mirror = float(np.max(abs(zero['elements']['naca'][label]+zero['elements']['naca'][label][::-1])))
        anti = float(np.max(abs(plus['elements']['naca'][label]+minus['elements']['naca'][label][::-1])))
        lift_sum = abs(plus['CL_pressure'][label]+minus['CL_pressure'][label])
        row = dict(control='symmetric_finite_te', method=label, panels=200, zero_incidence_CL=zero['CL_pressure'][label],
            mirror_speed_difference=mirror, antisymmetry_speed_difference=anti, antisymmetry_CL_sum=lift_sum,
            CL_pressure=plus['CL_pressure'][label], CL_circulation=plus[label+'_far_CL'], far_flux=plus[label+'_far_flux'],
            linear_residual=plus[label+'_linear_residual'], kutta_residual=plus[label+'_kutta_residual'])
        ok = ok and max(abs(row['zero_incidence_CL']), mirror, anti, lift_sum, row['linear_residual'], row['kutta_residual']) < LINEAR
        rows.append(row)
    flap = clockwise(transform(naca0012_open(120), .3, -20., (0., 0.), (1.02, -.05)))
    main = clockwise(naca0012_open(160))
    a = Pair([('main', main), ('flap', flap)], far=False).solve(4.); b = Pair([('flap', flap), ('main', main)], far=False).solve(4.)
    for label in ('production', 'hess'):
        order = max(float(np.max(abs(a['elements'][n][label]-b['elements'][n][label]))) for n in ('main', 'flap'))
        rows.append(dict(control='element_order', method=label, panels=280, ordering_speed_difference=order,
            CL_pressure=a['CL_pressure'][label], linear_residual=a[label+'_linear_residual'], kutta_residual=a[label+'_kutta_residual']))
        ok = ok and order < LINEAR and rows[-1]['linear_residual'] < LINEAR and rows[-1]['kutta_residual'] < LINEAR
    return rows, ok


def uniform(rows):
    keys = list(dict.fromkeys(k for r in rows for k in r))
    return [{k: r.get(k) for k in keys} for r in rows]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=OUT)
    args = parser.parse_args(); out = args.out.resolve(); out.mkdir(parents=True, exist_ok=True)
    require(not (out/'summary.json').exists(), 'Preserve completed check; choose a new output directory')
    started = time.monotonic(); inputs = Inputs(); inputs.add(PROTOCOL)
    control_rows, qualified = controls(); csv_write(out/'controls.csv', uniform(control_rows))
    print('Controls qualified:', qualified, flush=True)
    receipt = inputs.json(LOWER_RECEIPT); fore_receipt = inputs.json(ROOT/deck.RECEIPT)
    def archived(study, record, job, name):
        require(sha(inputs.add(study/job/name)) == record[job+'/'+name], 'Archived file changed: '+job+'/'+name)
        return study/job/name
    angle = lambda study, record, job: inputs.json(archived(study, record, job, 'result.json'))['alpha_deg']
    spec = inputs.json(LOWER/'study.json')
    measured = {j['name']: j['experiment_alpha'] for j in spec['jobs']}
    (lo, hi, _), _ = deck.archived_endpoints()
    incidences = {'nominal': [('measured_cl0056', measured['cl0056_nominal']), ('measured_cl0284', measured['cl028_control']),
            ('matched_cl0056', angle(LOWER, receipt, 'cl0056_nominal')), ('matched_cl0175', angle(LOWER, receipt, 'cl0175_nominal')),
            ('matched_cl0284', angle(LOWER, receipt, 'cl028_control')), ('bracket_lower', lo), ('bracket_upper', hi)],
        'fore_source_cubic': [('matched_cl0056_cubic', angle(FORE, fore_receipt, 'cl0056_fore_source_cubic')),
            ('bracket_lower', lo), ('bracket_upper', hi)]}
    saved = {'nominal': (LOWER, receipt, 'cl0056_nominal'), 'fore_source_cubic': (FORE, fore_receipt, 'cl0056_fore_source_cubic')}
    taps = {c: [(r['surface'], float(r['x_over_c'])) for r in inputs.csv(TAPS)
                if r['case'] == 'cl0056' and r['method'] == 'multifoil' and r['coordinate'] == c and r['included'] == 'True']
            for c in ('table', 'plotted')}
    cases, differences, refinement, profiles, geometry = [], [], [], [], []
    for scenario, angles in incidences.items():
        raw = deck.geometry_elements(scenario, list(SCHEDULE[0])); solutions = {}; rule = {}
        for panels in SCHEDULE:
            sampled = [(name, *sample(xy, count+1)) for (name, xy), count in zip(raw, panels)]
            polygons = [(name, x) for name, x, _ in sampled]
            for name, _, recovered in sampled:
                # The production monotone-spacing recovery must act identically on every grid,
                # otherwise the refined nodes would sample a different rule, not the same contour.
                require(rule.setdefault(name, recovered) == recovered, 'Sampling rule changes with resolution: '+name)
            if panels == SCHEDULE[0]:
                for (name, xy), (_, x) in zip(raw, polygons):
                    model = new_model(name, xy, 1e6, MACH, 9., (.25, 0.), 1e-5)
                    require(np.array_equal(model.foil.x, x), 'Baseline sampling differs from production nodes')
                    study, record, job = saved[scenario]
                    with np.load(archived(study, record, job, name+'_state.npz')) as state:
                        error = float(np.max(abs(state['xy']-x)))
                    require(error <= 2e-12, 'Polygon differs from archived saved-state nodes')
                    require(bool(model.paneling_recovery) == rule[name], 'Recovery flag differs from production')
                    geometry.append(dict(scenario=scenario, element=name, nodes=x.shape[1], archived_state=job,
                        monotone_spacing_recovery_all_grids=rule[name],
                        maximum_difference_from_archived_nodes=error, trailing_edge_gap=float(np.linalg.norm(x[:, 0]-x[:, -1]))))
            pair = Pair(polygons); print(scenario, panels, 'operators built', round(time.monotonic()-started), 's', flush=True)
            for label, alpha in angles:
                s = solutions[(panels, label)] = pair.solve(alpha)
                metrics = {m: {e: surface_metrics(s['elements'][e], m) for e in ('fore', 'aft')} for m in ('production', 'hess')}
                s['metrics'] = metrics
                for m in ('production', 'hess'):
                    row = dict(scenario=scenario, panels_fore=panels[0], panels_aft=panels[1], incidence=label, alpha_deg=alpha, method=m,
                        CL_pressure=s['CL_pressure'][m], CL_far_circulation=s[m+'_far_CL'], far_flux=s[m+'_far_flux'],
                        linear_residual=s[m+'_linear_residual'], kutta_residual=s[m+'_kutta_residual'])
                    if m == 'production': row['CL_pressure_nodal_closed'] = s['CL_pressure']['production_nodal']
                    else: row.update(CL_vortex_strength=s['hess_CL_circulation'],
                        source_flux_fore=s['hess_source_flux']['fore'], source_flux_aft=s['hess_source_flux']['aft'])
                    for e in ('fore', 'aft'):
                        for side in ('upper', 'lower'):
                            row[f'{e}_{side}_CN'] = metrics[m][e][side]['CN_incompressible']
                        row[f'{e}_stagnation_s'] = metrics[m][e]['stagnation']['s']
                        row[f'{e}_stagnation_x'] = metrics[m][e]['stagnation']['x']
                    for side in ('upper', 'lower'):
                        row[f'fore_{side}_nose_min_Cp'] = metrics[m]['fore'][side]['nose_min_Cp']
                        row[f'fore_{side}_nose_min_Cp_x'] = metrics[m]['fore'][side]['nose_min_Cp_x']
                    cases.append(row)
                    require(row['linear_residual'] < LINEAR and row['kutta_residual'] < LINEAR, 'Linear criterion failed')
                d = dict(scenario=scenario, panels_fore=panels[0], panels_aft=panels[1], incidence=label, alpha_deg=alpha)
                for e in ('fore', 'aft'):
                    el = s['elements'][e]; du = el['hess']-el['production']; dcp = cp_kt(el['hess'])-cp_kt(el['production'])
                    d[f'{e}_max_speed_difference'] = float(np.max(abs(du))); d[f'{e}_rms_speed_difference'] = float(np.sqrt(np.mean(du**2)))
                    d[f'{e}_max_Cp_difference'] = float(np.max(abs(dcp))); d[f'{e}_rms_Cp_difference'] = float(np.sqrt(np.mean(dcp**2)))
                    d[f'{e}_max_Cp_difference_x'] = float(el['xm'][0][np.argmax(abs(dcp))])
                    d[f'{e}_stagnation_s_difference'] = metrics['hess'][e]['stagnation']['s']-metrics['production'][e]['stagnation']['s']
                    for side in ('upper', 'lower'):
                        d[f'{e}_{side}_CN_difference'] = metrics['hess'][e][side]['CN_incompressible']-metrics['production'][e][side]['CN_incompressible']
                fore = s['elements']['fore']; nose = fore['xm'][0] < fore['x'][0].min()+NOSE
                d['fore_nose_max_Cp_difference'] = float(np.max(abs(cp_kt(fore['hess'][nose])-cp_kt(fore['production'][nose]))))
                for side in ('upper', 'lower'):
                    d[f'fore_{side}_nose_min_Cp_difference'] = metrics['hess']['fore'][side]['nose_min_Cp']-metrics['production']['fore'][side]['nose_min_Cp']
                d['CL_pressure_difference'] = s['CL_pressure']['hess']-s['CL_pressure']['production']
                d['CL_circulation_difference'] = s['hess_far_CL']-s['production_far_CL']
                for coordinate, positions in taps.items():
                    dt = tap_values(s, 'hess', positions)-tap_values(s, 'production', positions)
                    d[f'taps_{coordinate}_outside_midpoint_range'] = int(np.sum(~np.isfinite(dt))); dt = dt[np.isfinite(dt)]
                    d[f'taps_{coordinate}_count'] = len(dt); d[f'taps_{coordinate}_max_Cp_difference'] = float(np.max(abs(dt)))
                    d[f'taps_{coordinate}_rms_Cp_difference'] = float(np.sqrt(np.mean(dt**2)))
                differences.append(d)
                if panels in (SCHEDULE[0], SCHEDULE[-1]):
                    keep = np.flatnonzero(fore['xm'][0] < fore['x'][0].min()+2.5*NOSE)
                    for k in keep:
                        profiles.append(dict(scenario=scenario, panels_fore=panels[0], incidence=label, alpha_deg=alpha,
                            panel=int(k), s=float(fore['s_mid'][k]), x=float(fore['xm'][0, k]), y=float(fore['xm'][1, k]),
                            production_speed=float(fore['production'][k]), hess_speed=float(fore['hess'][k]),
                            production_Cp=float(cp_kt(fore['production'][k])), hess_Cp=float(cp_kt(fore['hess'][k]))))
        for label, alpha in angles:
            base, fine = solutions[(SCHEDULE[0], label)], solutions[(SCHEDULE[-1], label)]
            fore = base['elements']['fore']; fractions = fore['s_mid']/fore['s_nodes'][-1]
            nose = fore['xm'][0] < fore['x'][0].min()+NOSE
            row = dict(scenario=scenario, incidence=label, alpha_deg=alpha)
            between = {}
            for panels in SCHEDULE:
                s = solutions[(panels, label)]; e = s['elements']['fore']
                cp = {m: cp_kt(at_positions(e, m, fractions))[nose] for m in ('production', 'hess')}
                between[panels] = float(np.max(abs(cp['hess']-cp['production'])))
                row[f'between_methods_nose_Cp_{panels[0]}'] = between[panels]
                for m in ('production', 'hess'):
                    reference = cp_kt(at_positions(fine['elements']['fore'], m, fractions))[nose]
                    row[f'{m}_nose_Cp_change_to_finest_{panels[0]}'] = float(np.max(abs(cp[m]-reference)))
                    row[f'{m}_CL_pressure_{panels[0]}'] = s['CL_pressure'][m]
                    row[f'{m}_far_flux_{panels[0]}'] = s[m+'_far_flux']
                    row[f'{m}_fore_lower_nose_min_Cp_{panels[0]}'] = s['metrics'][m]['fore']['lower']['nose_min_Cp']
            spread = max(row[f'{m}_nose_Cp_change_to_finest_{SCHEDULE[0][0]}'] for m in ('production', 'hess'))
            values = [between[p] for p in SCHEDULE]
            row['resolved_discretization_spread_baseline'] = spread
            row['between_methods_decreasing'] = bool(values[0] > values[1] > values[2])
            row['agree_within_resolved_spread'] = bool(row['between_methods_decreasing'] and values[0] <= spread)
            refinement.append(row)
    for name, rows in (('s414_cases.csv', cases), ('s414_differences.csv', differences), ('refinement.csv', refinement),
                       ('fore_nose_profiles.csv', profiles), ('geometry_checks.csv', geometry)):
        csv_write(out/name, uniform(rows))
    for name in ('multifoil.py', 'viscous.py', 'paneling.py', 'vendor/mfoil.py', 'validation/s414_geometry.py',
                 'validation/run_s414_transition_boundary.py', 'validation/report_s414_formulation_check.py',
                 deck.QUAL+'/report_test_fore.dat', deck.QUAL+'/report_test_aft.dat'):
        inputs.add(ROOT/name)
    write(out/'input_manifest.json', dict(files=inputs.files))
    base = [d for d in differences if d['panels_fore'] == SCHEDULE[0][0]]; fine = [d for d in differences if d['panels_fore'] == SCHEDULE[-1][0]]
    summary = dict(complete=True, controls_qualified=bool(qualified), schedule=[list(p) for p in SCHEDULE],
        contour='Spline of the unchanged baseline solver input, production curvature spacing; baseline equals production nodes exactly',
        geometries=list(incidences), incidences={k: [dict(name=n, alpha_deg=a) for n, a in v] for k, v in incidences.items()},
        solves=len(cases)//2, maximum_linear_residual=max(r['linear_residual'] for r in cases),
        maximum_kutta_residual=max(r['kutta_residual'] for r in cases),
        maximum_polygon_difference_from_archived_nodes=max(g['maximum_difference_from_archived_nodes'] for g in geometry),
        baseline=dict(max_fore_nose_Cp_difference=max(d['fore_nose_max_Cp_difference'] for d in base),
            max_Cp_difference_any_panel=max(max(d['fore_max_Cp_difference'], d['aft_max_Cp_difference']) for d in base),
            max_tap_Cp_difference_table=max(d['taps_table_max_Cp_difference'] for d in base),
            max_tap_Cp_difference_plotted=max(d['taps_plotted_max_Cp_difference'] for d in base),
            max_abs_CL_pressure_difference=max(abs(d['CL_pressure_difference']) for d in base),
            max_abs_CL_circulation_difference=max(abs(d['CL_circulation_difference']) for d in base)),
        finest=dict(max_fore_nose_Cp_difference=max(d['fore_nose_max_Cp_difference'] for d in fine),
            max_Cp_difference_any_panel=max(max(d['fore_max_Cp_difference'], d['aft_max_Cp_difference']) for d in fine),
            max_tap_Cp_difference_table=max(d['taps_table_max_Cp_difference'] for d in fine),
            max_abs_CL_pressure_difference=max(abs(d['CL_pressure_difference']) for d in fine),
            max_abs_CL_circulation_difference=max(abs(d['CL_circulation_difference']) for d in fine)),
        cases_agreeing_within_resolved_spread=sum(r['agree_within_resolved_spread'] for r in refinement),
        cases_with_decreasing_difference=sum(r['between_methods_decreasing'] for r in refinement), cases=len(refinement),
        supported_taps=dict(table=len(taps['table']), plotted=len(taps['plotted'])),
        new_nonlinear_flow_solves=0, wall_seconds=time.monotonic()-started,
        runtime=dict(python=platform.python_version(), versions={k: importlib.metadata.version(k) for k in ('numpy', 'scipy')}),
        qualification='Cross-formulation numerical comparison on identical polygons; not validation against the viscous experiment')
    write(out/'summary.json', summary)
    print('complete; controls', qualified, '; agreeing cases', summary['cases_agreeing_within_resolved_spread'], '/', len(refinement))


if __name__ == '__main__':
    main()
