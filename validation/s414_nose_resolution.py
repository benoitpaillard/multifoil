"""Diagnostic subdivision of the frozen S414 fore-lower solver polygon.

No source interpolant, closure or production module is changed. Original
station coordinates/arclength/tangents are retained exactly. Inserted points
are collinear with their parent panel; arclength parameters interpolate the
archived values, so this does not silently replace the legacy metric.
"""
from pathlib import Path
import hashlib
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from diagnostic_discretization import DiscretizedSolver
from multifoil import validate_geometry
from viscous import CoupledSolver

FACTORS = (1, 2, 4)
NOSE_LIMIT = 0.04


def build_elements():
    from s414_geometry import elements
    return elements('nominal', [340, 240])


def digest(array):
    return hashlib.sha256(np.ascontiguousarray(array).tobytes()).hexdigest()


def subdivide(xy, arclength, tangent, factor, nose_limit=NOSE_LIMIT):
    """Split lower-TE-to-LE panels fully inside the declared total-chord window."""
    if type(factor) is not int or factor not in FACTORS:
        raise ValueError('subdivision_factor must be 1, 2 or 4')
    xy, arclength, tangent = map(np.asarray, (xy, arclength, tangent))
    n = xy.shape[1]
    if (xy.shape != (2, n) or tangent.shape != (2, n) or
            arclength.shape != (n,) or not np.isfinite(xy).all() or
            not np.isfinite(arclength).all() or not np.isfinite(tangent).all() or
            np.any(np.diff(arclength) <= 0)):
        raise ValueError('Invalid ordered panel data')
    te = .5*(xy[:, 0]+xy[:, -1])
    le = int(np.argmax(np.linalg.norm(xy-te[:, None], axis=0)))
    selected = [j for j in range(le) if max(xy[0, j], xy[0, j+1]) <= nose_limit]
    if not selected:
        raise ValueError('No fore-lower nose panels in the declared window')
    points, distances, tangents, old_to_new = [], [], [], []
    for j in range(n):
        old_to_new.append(len(points))
        points.append(xy[:, j].copy())
        distances.append(arclength[j])
        tangents.append(tangent[:, j].copy())
        if j not in selected or factor == 1:
            continue
        dx = xy[:, j+1]-xy[:, j]
        ds = arclength[j+1]-arclength[j]
        for k in range(1, factor):
            f = k/factor
            points.append(xy[:, j]+f*dx)
            distances.append(arclength[j]+f*ds)
            tangents.append(dx/ds)
    new_xy, new_s, new_t = np.asarray(points).T, np.asarray(distances), np.asarray(tangents).T
    mapping = np.asarray(old_to_new, dtype=int)
    if not (np.array_equal(new_xy[:, mapping], xy) and
            np.array_equal(new_s[mapping], arclength) and
            np.array_equal(new_t[:, mapping], tangent)):
        raise ValueError('Subdivision moved an original station')
    metadata = dict(
        subdivision_factor=factor, original_nodes=n, refined_nodes=len(new_s),
        nose_limit_total_chord=nose_limit, leading_edge_index=le,
        selected_parent_panels=selected, original_to_refined=mapping.tolist(),
        original_xy_sha256=digest(xy), refined_xy_sha256=digest(new_xy),
        original_s_sha256=digest(arclength), refined_s_sha256=digest(new_s),
        metric='Original foil.s retained; inserted s linear on each parent interval',
        tangent='Original tangents retained; inserted tangent dx/ds of parent panel',
        maximum_parent_arc_minus_chord=float(np.max(
            np.diff(arclength)[selected]-np.linalg.norm(np.diff(xy, axis=1)[:, selected], axis=0))),
    )
    return new_xy, new_s, new_t, metadata


class Solver(DiscretizedSolver):
    """Native factor-1 path; refined panels with original wake spacing/counts."""
    def __init__(self, elements, *, subdivision_factor=1, **kwargs):
        if type(subdivision_factor) is not int or subdivision_factor not in FACTORS:
            raise ValueError('subdivision_factor must be 1, 2 or 4')
        if [name for name, _ in elements] != ['fore', 'aft']:
            raise ValueError('Require declared fore/aft geometry in its original order')
        if any(k in kwargs for k in ('surface_nodes', 'wake_nodes', 'wake_first_spacing')):
            raise ValueError('Wake and surface controls are frozen by the diagnostic')
        self.subdivision_factor = subdivision_factor
        self.resolution_controls = {}
        self.baseline_geometry = {}
        super().__init__(elements, **kwargs)

    def _build_inviscid(self):
        te_points = [.5*(m.foil.x[:, 0]+m.foil.x[:, -1]) for m in self.models]
        for m, pt in zip(self.models, te_points):
            name = m.geom.name
            x, s, t = m.foil.x.copy(), m.foil.s.copy(), m.foil.t.copy()
            self.baseline_geometry[name] = dict(xy=x, s=s, t=t)
            wake_count = max(16, int(np.ceil(m.foil.N/10+10+20*np.log(max(
                1., (self.wake_length+max(np.linalg.norm(pt-other) for other in te_points))/m.local_chord)))))
            first_spacing = float(.5*((s[1]-s[0])+(s[-1]-s[-2])))
            if name == 'fore':
                new_x, new_s, new_t, metadata = subdivide(x, s, t, self.subdivision_factor)
                if self.subdivision_factor != 1:
                    m.foil.x, m.foil.s, m.foil.t = new_x, new_s, new_t
                    m.foil.N = new_s.size
            else:
                metadata = dict(subdivision_factor=1, original_nodes=m.foil.N,
                    refined_nodes=m.foil.N, original_xy_sha256=digest(x),
                    refined_xy_sha256=digest(x), original_to_refined=list(range(m.foil.N)))
            metadata.update(wake_nodes=wake_count, wake_first_spacing=first_spacing,
                            wake_length=self.wake_length)
            self.resolution_controls[name] = metadata
            if self.subdivision_factor != 1:
                self.discretization['wake_nodes'][name] = wake_count
                self.discretization['wake_first_spacing'][name] = first_spacing
        self.elements = [(m.geom.name, m.foil.x.T[::-1].copy()) for m in self.models]
        validate_geometry(self.elements)
        self.fs = []
        start = 0
        for m in self.models:
            self.fs.append(slice(start, start+m.foil.N))
            start += m.foil.N
        self.nfoil = start
        if self.nfoil > 1200:
            raise ValueError('Refinement exceeds 1200 surface nodes')
        CoupledSolver._build_inviscid(self)

    def _solve_point(self, alpha):
        result = super()._solve_point(alpha)
        result.diagnostics['nose_resolution'] = self.resolution_controls
        return result
