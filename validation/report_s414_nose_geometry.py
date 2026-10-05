#!/usr/bin/env python3
"""V4: trace the S414 fore-nose turning reversal and test its inviscid pressure effect.

Geometry trace, admission of the declared alternatives, then a linear pressure
comparison on final solver nodes. No viscous or nonlinear solve.
Protocol: validation/S414_NOSE_GEOMETRY_PROTOCOL.md.
"""
import os
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
import argparse
import importlib.metadata
from pathlib import Path
import platform
import sys
import numpy as np
from scipy.interpolate import CubicSpline

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from report_design_initialization import Inputs, sha, write, csv_write
from report_s702_drag_mechanism import require
from viscous import new_model
from s414_geometry import cubic_repanel, elements as declared_elements, QUAL
import s414_nose_geometry as nose
import report_s414_formulation_check as check

OUT = ROOT/'validation/results/s414_nose_geometry_20261004'
PROTOCOL = ROOT/'validation/S414_NOSE_GEOMETRY_PROTOCOL.md'
LOWER = ROOT/'validation/results/s414_lower_bucket_20261003'
T2 = ROOT/'validation/results/s414_transition_boundary_t2_20261004'
TAPS = ROOT/'validation/results/s414_lower_bucket_analysis_20261003/pressure_at_taps.csv'
CHORD_MM, MATERIAL, GRID_SENSITIVITY = 457.2, .05, .013
RECOVERY = (.008, .06)
STAGES = ('nominal', 'archived_cubic', 'single_curve', 'convex_nose')


def to_polyline(points, polyline, chunk=256):
    """Distance from each point to an open polyline, evaluated in chunks to bound memory."""
    a, b = polyline[:-1], polyline[1:]; d = b-a; size = np.sum(d*d, axis=1); out = np.empty(len(points))
    for start in range(0, len(points), chunk):
        offset = points[start:start+chunk, None, :]-a[None]; t = np.clip(np.sum(offset*d[None], axis=2)/size, 0, 1)
        out[start:start+chunk] = np.min(np.linalg.norm(offset-t[:, :, None]*d[None], axis=2), axis=1)
    return out


def concave_intervals(points):
    """Runs of negative turning on a dense polyline: x-range, integrated turn (deg), most negative curvature."""
    turn = nose.turns(points); length = np.linalg.norm(np.diff(points, axis=0), axis=1)
    curvature = np.radians(turn)/(.5*(length[:-1]+length[1:])); negative = turn < 0
    edges = np.flatnonzero(np.diff(np.r_[0, negative.astype(int), 0])); out = []
    for a, b in zip(edges[::2], edges[1::2]):
        out.append(dict(first_index=int(a+1), x_start=float(points[a+1, 0]), x_end=float(points[b, 0]),
            integrated_turn_degrees=float(turn[a:b].sum()), minimum_curvature=float(curvature[a:b].min())))
    return out, curvature


def curve_function(stage, xy, le, s):
    if stage == 'archived_cubic':
        # The archived fore-only variant: natural cubic per branch (s414_geometry.cubic_repanel).
        upper = CubicSpline(s[:le+1], xy[:le+1], axis=0, bc_type='natural'); lower = CubicSpline(s[le:]-s[le], xy[le:], axis=0, bc_type='natural')
        return lambda v: upper(v) if v <= s[le] else lower(v-s[le])
    curve = nose.Reconstruction(stage)
    return lambda v: curve.point(float(v))


def dense_curve(stage, xy, le, s, i0, i1, count=40001):
    q = np.linspace(s[i0], s[i1], count); point = curve_function(stage, xy, le, s)
    return q, np.array([point(v) for v in q]), point


def final_nodes(name, points):
    return new_model(name, points, 1e6, check.MACH, 9., (.25, 0.), 1e-5).foil.x


def nose_turns(polygon):
    """Negative turns of a counter-clockwise polygon inside the nose support, and the smallest turn there."""
    turn = nose.turns(polygon); inside = polygon[1:-1, 0] <= nose.NOSE
    rows = [dict(x=float(polygon[k+1, 0]), y=float(polygon[k+1, 1]), turn_degrees=float(turn[k]),
                 surface='upper' if k+1 < int(np.argmin(polygon[:, 0])) else 'lower')
            for k in np.flatnonzero(inside & (turn < 0))]
    return rows, float(turn[inside].min())


def surface_samples(element, label):
    """Midpoint x, arclength, speed and Cp on the lower and upper surface of one clockwise polygon."""
    x = element['x']; ile = int(np.argmin(x[0])); n = x.shape[1]; out = {}
    for side, panels in (('lower', slice(0, ile)), ('upper', slice(ile, n-1))):
        out[side] = dict(x=element['xm'][0, panels], s=element['s_mid'][panels], speed=abs(element[label][panels]),
                         cp=check.cp_kt(element[label][panels]))
    return out


def at_x(samples, query):
    order = np.argsort(samples['x'])
    return np.interp(query, samples['x'][order], samples['cp'][order])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=OUT)
    out = parser.parse_args().out.resolve(); out.mkdir(parents=True, exist_ok=True)
    require(not (out/'summary.json').exists(), 'Preserve completed study; choose a new output directory')
    inputs = Inputs(); inputs.add(PROTOCOL); inputs.add(ROOT/QUAL/'report_test_fore.dat'); inputs.add(ROOT/QUAL/'report_test_aft.dat')
    xy, le = nose.source(); s = nose.chord_parameter(xy); i0, i1 = nose.support(xy, le)
    # ---- 1. reported polygon
    turn = nose.turns(xy); length = np.linalg.norm(np.diff(xy, axis=0), axis=1)
    vertices = [dict(vertex=i, x=float(xy[i, 0]), y=float(xy[i, 1]), surface='upper' if i < le else 'lower' if i > le else 'split vertex',
        incoming_length=float(length[i-1]), outgoing_length=float(length[i]), turn_degrees=float(turn[i-1]),
        discrete_curvature=float(np.radians(turn[i-1])/(.5*(length[i-1]+length[i]))), in_nose_support=bool(i0 <= i <= i1))
        for i in range(i0-2, i1+3)]
    csv_write(out/'reported_nose_vertices.csv', vertices)
    require(all(v['turn_degrees'] > 0 for v in vertices if v['in_nose_support']), 'Reported polygon is not convex in the nose support')
    # ---- 2. continuous curves, input samples and final solver nodes
    inputs_by_stage = {'nominal': nose.fore_input('nominal'), 'archived_cubic': cubic_repanel(xy, 340),
                       'single_curve': nose.fore_input('single_curve'), 'convex_nose': nose.fore_input('convex_nose')}
    require(np.array_equal(inputs_by_stage['nominal'], declared_elements('nominal', [340, 240])[0][1]), 'Nominal input differs from production')
    aft_input = declared_elements('nominal', [340, 240])[1][1]
    with np.load(inputs.add(LOWER/'cl0056_nominal/fore_state.npz')) as saved: archived_fore = saved['xy'].copy()
    with np.load(inputs.add(LOWER/'cl0056_nominal/aft_state.npz')) as saved: archived_aft = saved['xy'].copy()
    nodes = {stage: final_nodes('fore', points) for stage, points in inputs_by_stage.items()}; aft_nodes = final_nodes('aft', aft_input)
    require(float(np.max(abs(nodes['nominal']-archived_fore))) <= 2e-12 and float(np.max(abs(aft_nodes-archived_aft))) <= 2e-12,
            'Rebuilt nominal nodes differ from the archived saved state')
    q, nominal_dense, _ = dense_curve('nominal', xy, le, s, i0, i1); nominal_polygon = nodes['nominal'].T[::-1]
    concavity, polygons, departures, curvature_rows, node_rows = [], [], [], [], []
    for stage in STAGES:
        _, dense, point = dense_curve(stage, xy, le, s, i0, i1); intervals, curvature = concave_intervals(dense)
        for row in intervals:
            if abs(row['integrated_turn_degrees']) >= 1e-3:
                concavity.append(dict(curve=stage, surface='upper' if q[row['first_index']] < s[le] else 'lower',
                    **{k: v for k, v in row.items() if k != 'first_index'}))
        lower = q >= s[le]; keep = np.flatnonzero(lower[1:-1])[::40]
        curvature_rows += [dict(curve=stage, x=float(dense[k+1, 0]), curvature=float(curvature[k])) for k in keep]
        points = inputs_by_stage[stage]; polygon = nodes[stage].T[::-1]
        for label, poly in (('input_samples', points), ('final_solver_nodes', polygon)):
            negative, smallest = nose_turns(poly)
            polygons.append(dict(variant=stage, stage=label, vertices=len(poly), negative_turns_in_nose_support=len(negative),
                smallest_turn_degrees=smallest, negative_turn_positions='; '.join(f"{r['surface']} x={r['x']:.5f} {r['turn_degrees']:.3f}°" for r in negative)))
        turns_final = nose.turns(polygon); spacing = np.linalg.norm(np.diff(polygon, axis=0), axis=1); ile = int(np.argmin(polygon[:, 0]))
        node_rows += [dict(variant=stage, node=k+1, surface='upper' if k+1 < ile else 'lower', x=float(polygon[k+1, 0]), y=float(polygon[k+1, 1]),
                           turn_degrees=float(turns_final[k]), left_spacing=float(spacing[k]), right_spacing=float(spacing[k+1]))
                      for k in range(len(turns_final)) if polygon[k+1, 0] <= 2*nose.NOSE]
        moved = np.linalg.norm(points-inputs_by_stage['nominal'], axis=1); inside = polygon[:, 0] <= nose.NOSE
        offset = to_polyline(polygon, nominal_polygon)
        departures.append(dict(variant=stage, input_samples_changed=int(np.sum(moved > 0)),
            maximum_input_sample_displacement_mm=float(moved.max()*CHORD_MM),
            changed_samples_x_minimum=float(points[moved > 0, 0].min()) if (moved > 0).any() else None,
            changed_samples_x_maximum=float(points[moved > 0, 0].max()) if (moved > 0).any() else None,
            maximum_curve_distance_from_nominal_curve_mm=float(to_polyline(dense[::8], nominal_dense[::4]).max()*CHORD_MM),
            maximum_curve_distance_from_reported_polygon_mm=float(to_polyline(dense[::8], xy[i0:i1+1]).max()*CHORD_MM),
            maximum_reported_vertex_miss_mm=float(max(np.linalg.norm(point(s[i])-xy[i]) for i in range(i0, i1+1))*CHORD_MM),
            final_node_offset_from_nominal_polygon_inside_support_mm=float(offset[inside].max()*CHORD_MM),
            final_node_offset_from_nominal_polygon_outside_support_mm=float(offset[~inside].max()*CHORD_MM),
            trailing_edge_gap=float(np.linalg.norm(polygon[0]-polygon[-1]))))
    for name, rows in (('curve_concavity.csv', concavity), ('polygon_turns.csv', polygons), ('departures.csv', departures),
                       ('lower_nose_curvature.csv', curvature_rows), ('final_nodes_nose.csv', node_rows)):
        csv_write(out/name, check.uniform(rows))
    convex = nose.Reconstruction('convex_nose')
    csv_write(out/'convex_nose_arcs.csv', [dict(vertex=i, x=float(xy[i, 0]), y=float(xy[i, 1]),
        tangent_angle_degrees=float(np.degrees(convex.tangent_angles[i])),
        first_arc_curvature=convex.arc_curvatures[i][0] if i < i1 else None, second_arc_curvature=convex.arc_curvatures[i][1] if i < i1 else None)
        for i in range(i0, i1+1)])
    # ---- 3. admission, decided on geometry alone
    final = {r['variant']: r for r in polygons if r['stage'] == 'final_solver_nodes'}; moved = {r['variant']: r for r in departures}
    admission = {}
    for variant in ('single_curve', 'convex_nose'):
        replaced = nose.Reconstruction(variant).changed(nose.cosine_parameters(s, le, 340))
        admission[variant] = dict(final_solver_negative_turns_in_nose_support=final[variant]['negative_turns_in_nose_support'],
            reported_vertices_reproduced=bool(moved[variant]['maximum_reported_vertex_miss_mm'] < 1e-9),
            outside_support_samples_identical=bool(np.array_equal(inputs_by_stage[variant][~replaced], inputs_by_stage['nominal'][~replaced])
                                                   and inputs_by_stage[variant][replaced][:, 0].max() <= nose.NOSE),
            trailing_edge_unchanged=bool(abs(moved[variant]['trailing_edge_gap']-moved['nominal']['trailing_edge_gap']) < 1e-15))
        admission[variant]['admitted'] = bool(admission[variant]['final_solver_negative_turns_in_nose_support'] == 0
            and all(v for k, v in admission[variant].items() if k != 'final_solver_negative_turns_in_nose_support'))
    for variant in ('nominal', 'single_curve', 'convex_nose'):
        np.savetxt(out/f'fore_input_{variant}.dat', inputs_by_stage[variant], fmt='%.17g',
                   header=f'S414 fore solver input, {variant}; 341 points, upper TE first; V4 reconstruction, not the tested contour', comments='')
    # ---- 4. inviscid pressure on final nodes
    spec = inputs.json(LOWER/'study.json'); t2 = inputs.json(T2/'study.json')
    angle = lambda job: inputs.json(LOWER/job/'result.json')['alpha_deg']
    measured = next(j['experiment_alpha'] for j in spec['jobs'] if j['name'] == 'cl0056_nominal')
    incidences = [('measured_cl0056', measured), ('matched_cl0056', angle('cl0056_nominal')),
        ('transition_interval_lower', t2['nominal_new_incidences_deg'][-1]), ('transition_interval_upper', t2['tranche1_angles_deg'][0]),
        ('matched_cl0175', angle('cl0175_nominal')), ('matched_cl0284', angle('cl028_control'))]
    taps = sorted({(r['surface'].split('_')[1], float(r['x_over_c'])) for r in inputs.csv(TAPS) if r['case'] == 'cl0056' and r['method'] == 'multifoil'
                   and r['coordinate'] == 'table' and r['included'] == 'True' and r['surface'].startswith('fore') and float(r['x_over_c']) <= nose.NOSE})
    query = dict(lower=np.linspace(*RECOVERY, 261), upper=np.linspace(.003, nose.NOSE, 286))
    solutions, cases, differences, profiles = {}, [], [], []
    for variant in ('nominal', 'single_curve', 'convex_nose'):
        pair = check.Pair([('fore', nodes[variant]), ('aft', aft_nodes)], far=False)
        for label, alpha in incidences:
            solved = pair.solve(alpha)
            require(max(solved['production_linear_residual'], solved['hess_linear_residual'], solved['production_kutta_residual'],
                        solved['hess_kutta_residual']) < check.LINEAR, 'Linear criterion failed')
            for method in ('production', 'hess'):
                surfaces = surface_samples(solved['elements']['fore'], method); low = surfaces['lower']
                window = low['x'] <= nose.NOSE; k = np.flatnonzero(window)[np.argmin(low['cp'][window])]
                recovery = (low['x'] >= RECOVERY[0]) & (low['x'] <= RECOVERY[1])
                order = np.argsort(low['s'][recovery]); gradient = np.diff(low['speed'][recovery][order])/np.diff(low['s'][recovery][order])
                slope = np.diff(low['cp'][recovery][order])
                row = dict(variant=variant, incidence=label, alpha_deg=alpha, method=method, CL_pressure=solved['CL_pressure'][method],
                    fore_lower_minimum_Cp=float(low['cp'][k]), fore_lower_minimum_Cp_x=float(low['x'][k]),
                    # Flow on the lower surface runs against the clockwise arclength, hence the sign.
                    strongest_deceleration_on_recovery=float(-gradient.max()),
                    recovery_Cp_slope_sign_changes=int(np.sum(np.diff(np.sign(slope)) != 0)))
                for side, x in taps: row[f'tap_{side}_{x:.5f}_Cp'] = float(at_x(surfaces[side], x))
                cases.append(row); solutions[(variant, label, method)] = (surfaces, row)
                if method == 'production':
                    for side in ('lower', 'upper'):
                        keep = surfaces[side]['x'] <= 2*nose.NOSE
                        profiles += [dict(variant=variant, incidence=label, alpha_deg=alpha, surface=side, x=float(a), s=float(b), Cp=float(c))
                                     for a, b, c in zip(surfaces[side]['x'][keep], surfaces[side]['s'][keep], surfaces[side]['cp'][keep])]
    for variant in ('single_curve', 'convex_nose'):
        for label, alpha in incidences:
            for method in ('production', 'hess'):
                (a, ra), (b, rb) = solutions[(variant, label, method)], solutions[('nominal', label, method)]
                delta = {side: at_x(a[side], query[side])-at_x(b[side], query[side]) for side in ('lower', 'upper')}
                row = dict(variant=variant, incidence=label, alpha_deg=alpha, method=method,
                    lower_recovery_max_abs_delta_Cp=float(np.max(abs(delta['lower']))),
                    lower_recovery_max_abs_delta_Cp_x=float(query['lower'][np.argmax(abs(delta['lower']))]),
                    upper_nose_max_abs_delta_Cp=float(np.max(abs(delta['upper']))),
                    upper_nose_max_abs_delta_Cp_x=float(query['upper'][np.argmax(abs(delta['upper']))]),
                    delta_minimum_Cp=ra['fore_lower_minimum_Cp']-rb['fore_lower_minimum_Cp'],
                    delta_minimum_Cp_x=ra['fore_lower_minimum_Cp_x']-rb['fore_lower_minimum_Cp_x'],
                    delta_strongest_deceleration=ra['strongest_deceleration_on_recovery']-rb['strongest_deceleration_on_recovery'],
                    delta_CL_pressure=ra['CL_pressure']-rb['CL_pressure'])
                for side, x in taps: row[f'tap_{side}_{x:.5f}_delta_Cp'] = ra[f'tap_{side}_{x:.5f}_Cp']-rb[f'tap_{side}_{x:.5f}_Cp']
                row['material'] = bool(row['lower_recovery_max_abs_delta_Cp'] > MATERIAL or abs(row['delta_minimum_Cp']) > MATERIAL)
                differences.append(row)
    for name, rows in (('pressure_cases.csv', cases), ('pressure_differences.csv', differences), ('nose_pressure_profiles.csv', profiles)):
        csv_write(out/name, check.uniform(rows))
    for name in ('multifoil.py', 'viscous.py', 'paneling.py', 'vendor/mfoil.py', 'validation/s414_geometry.py', 'validation/s414_nose_geometry.py',
                 'validation/report_s414_nose_geometry.py', 'validation/report_s414_formulation_check.py'):
        inputs.add(ROOT/name)
    write(out/'input_manifest.json', dict(files=inputs.files))
    def worst(variant, method, key):
        return max(abs(r[key]) for r in differences if r['variant'] == variant and r['method'] == method)
    decision = {variant: dict(admitted=admission[variant]['admitted'],
        production=dict(lower_recovery_max_abs_delta_Cp=worst(variant, 'production', 'lower_recovery_max_abs_delta_Cp'),
            max_abs_delta_minimum_Cp=worst(variant, 'production', 'delta_minimum_Cp'),
            upper_nose_max_abs_delta_Cp=worst(variant, 'production', 'upper_nose_max_abs_delta_Cp')),
        hess=dict(lower_recovery_max_abs_delta_Cp=worst(variant, 'hess', 'lower_recovery_max_abs_delta_Cp'),
            max_abs_delta_minimum_Cp=worst(variant, 'hess', 'delta_minimum_Cp')),
        material=any(r['material'] for r in differences if r['variant'] == variant and r['method'] == 'production'))
        for variant in ('single_curve', 'convex_nose')}
    summary = dict(complete=True, nose_support_vertices=[i0, i1], nose_support_x_limit=nose.NOSE,
        reported_vertices=len(xy), reported_nose_vertices=i1-i0+1,
        smallest_reported_turn=min((v for v in vertices if v['in_nose_support']), key=lambda v: v['turn_degrees']),
        polygon_turns=polygons, admission=admission, fairing=convex.fairing, incidences=[dict(name=n, alpha_deg=a) for n, a in incidences],
        material_threshold_Cp=MATERIAL, production_grid_sensitivity_Cp=GRID_SENSITIVITY, recovery_window_x=list(RECOVERY), decision=decision,
        frozen_inputs={f'fore_input_{v}.dat': sha(out/f'fore_input_{v}.dat') for v in ('nominal', 'single_curve', 'convex_nose')},
        new_nonlinear_flow_solves=0, runtime=dict(python=platform.python_version(),
            versions={k: importlib.metadata.version(k) for k in ('numpy', 'scipy')}),
        qualification='Reconstruction sensitivity on a figure-derived polygon; neither alternative is the tested contour. '
                      'Inviscid comparison only; no viscous effect is predicted')
    write(out/'summary.json', summary)
    print('admission', {k: v['admitted'] for k, v in admission.items()}); print('decision', decision)


if __name__ == '__main__':
    main()
