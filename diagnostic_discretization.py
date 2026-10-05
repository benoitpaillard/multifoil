"""Opt-in discretization experiment with the production viscous core unchanged.

The source contour, surface node count and wake controls can be varied
independently. Wake tracing below retains CoupledSolver's predictor/corrector
and source-system layout. This module changes no boundary-layer closure.
Controls are dictionaries keyed by element name. Lengths use assembly reference
chords; wake_first_spacing is the requested tracing step, not its curved chord.
"""
import numpy as np
from matplotlib.path import Path as Polygon

from multifoil import clean_contour, validate_geometry, _closed_segments, _intersections
from paneling import curvature_paneling
from viscous import CoupledSolver, te_info, mf


class DiscretizedSolver(CoupledSolver):
    def __init__(self, elements, *, surface_nodes=None, wake_nodes=None,
                 wake_first_spacing=None, **kwargs):
        names = [name for name, _ in elements]
        if len(names) != len(set(names)):
            raise ValueError('Discretization controls require unique element names')
        self.discretization = {}
        for label, values in (('surface_nodes', surface_nodes), ('wake_nodes', wake_nodes),
                              ('wake_first_spacing', wake_first_spacing)):
            values = {} if values is None else dict(values)
            if set(values)-set(names):
                raise ValueError('Unknown element in '+label)
            for value in values.values():
                if label.endswith('nodes'):
                    minimum = 41 if label == 'surface_nodes' else 16
                    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)) or value < minimum:
                        raise ValueError(label+' requires integer counts >= '+str(minimum))
                elif isinstance(value, (bool, np.bool_)) or not np.isfinite(value) or value <= 0:
                    raise ValueError('wake_first_spacing must be finite and positive')
            self.discretization[label] = values
        if sum(self.discretization['surface_nodes'].get(n, len(x)) for n, x in elements) > 1200:
            raise ValueError('Viscous limit is 1200 surface nodes')
        if sum(self.discretization['wake_nodes'].values()) > 1800:
            raise ValueError('Wake node limit exceeded')
        # Keep the original, pre-regularization source for the identical frame.
        self.contour_sources = [(n, np.array(x, copy=True)) for n, x in elements]
        super().__init__(elements, **kwargs)

    def _build_inviscid(self):
        controls = self.discretization['surface_nodes']
        for model, (name, xy) in zip(self.models, self.contour_sources):
            count = controls.get(name, model.foil.N)
            if count == model.foil.N:
                continue
            source = np.ascontiguousarray(clean_contour(xy)[::-1].T)
            te = .5*(source[:, 0]+source[:, -1])
            le = model.chord_origin
            chord = np.linalg.norm(te-le)
            axis = model.chord_axis
            frame = np.column_stack((axis, [-axis[1], axis[0]]))
            # The source remains the same regularized contour; only sampling changes.
            local = frame.T@(model.geom.xpoint-le[:, None])/chord
            y, s, tangent, recovered = curvature_paneling(local, count, 2., .1)
            model.foil.N = count
            model.foil.x = le[:, None]+chord*(frame@y)
            model.foil.s = chord*s
            model.foil.t = frame@tangent
            model.paneling_recovery = recovered
        self.elements = [(m.geom.name, m.foil.x.T[::-1].copy()) for m in self.models]
        validate_geometry(self.elements)
        self.fs = []; start = 0
        for model in self.models:
            self.fs.append(slice(start, start+model.foil.N)); start += model.foil.N
        self.nfoil = start
        super()._build_inviscid()

    def _build_wakes(self, alpha):
        nodes = self.discretization['wake_nodes']
        steps = self.discretization['wake_first_spacing']
        if not nodes and not steps:
            return super()._build_wakes(alpha)
        stream = np.array([np.cos(np.deg2rad(alpha)), np.sin(np.deg2rad(alpha))])
        te_points = [.5*(m.foil.x[:, 0]+m.foil.x[:, -1]) for m in self.models]
        end_plane = max(stream@pt for pt in te_points)+self.wake_length
        lengths = [max(self.wake_length, end_plane-stream@pt) for pt in te_points]
        counts = [nodes.get(m.geom.name, max(16, int(np.ceil(m.foil.N/10+10+20*np.log(
            max(1., (self.wake_length+max(np.linalg.norm(pt-other) for other in te_points))/m.local_chord))))))
            for m, pt in zip(self.models, te_points)]
        if self.nfoil+sum(counts) > 1800:
            raise ValueError('Surface plus wake node limit (1800) exceeded')
        for idx, (model, length, count) in enumerate(zip(self.models, lengths, counts)):
            auto_step = .5*((model.foil.s[1]-model.foil.s[0])+(model.foil.s[-1]-model.foil.s[-2]))
            step = steps.get(model.geom.name, auto_step)
            if step*(count-1) >= length:
                raise ValueError('Explicit wake controls require increasing geometric spacing')
            distances = mf.space_geom(step, length, count)
            if not np.isfinite(distances).all() or np.any(np.diff(distances) <= 0):
                raise ValueError('Invalid wake station distribution')
            t = te_info(model.foil.x)[0]
            x = np.zeros((2, count)); tangent = np.zeros_like(x)
            x[:, 0] = .5*(model.foil.x[:, 0]+model.foil.x[:, -1])+t*(1e-6*model.local_chord)
            for j in range(count-1):
                velocity = self._velocity(x[:, j:j+1], alpha)[:, 0]
                if np.linalg.norm(velocity) < 1e-5:
                    raise RuntimeError('Wake encountered stagnation')
                t0 = velocity/np.linalg.norm(velocity); ds = distances[j+1]-distances[j]
                pred = x[:, j]+ds*t0
                v2 = self._velocity(pred[:, None], alpha)[:, 0]; t1 = v2/np.linalg.norm(v2)
                x[:, j+1] = x[:, j]+.5*ds*(t0+t1); tangent[:, j] = t0; tangent[:, j+1] = t1
            for k, other in enumerate(self.models):
                if k != idx:
                    a, b = _closed_segments(other.foil.x.T)
                    if (Polygon(other.foil.x.T).contains_points(x.T).any()
                            or _intersections(x[:, :-1].T, x[:, 1:].T, a, b).any()):
                        raise RuntimeError(f'Wake from {model.geom.name} intersects {other.geom.name}; this topology is unsupported')
            model.wake.N = count; model.wake.x = x; model.wake.t = tangent
            model.wake.s = model.foil.s[-1]+np.r_[0, np.cumsum(np.linalg.norm(np.diff(x, axis=1), axis=0))]
            refs = []
            for angle in (0., 90.):
                refs.append(np.sum(self._velocity(x, angle)*tangent, axis=0))
            model.isol.uewiref = np.column_stack(refs)
        self.ss = []; self.sources = []; isys = isrc = 0
        for model in self.models:
            ns = model.foil.N+model.wake.N; nc = ns-2
            self.ss.append(slice(isys, isys+ns)); self.sources.append(slice(isrc, isrc+nc))
            isys += ns; isrc += nc
        self.nsys, self.nsource = isys, isrc

    def _solve_point(self, alpha):
        result = super()._solve_point(alpha)
        result.diagnostics['discretization_controls'] = {k: dict(v) for k, v in self.discretization.items()}
        return result
