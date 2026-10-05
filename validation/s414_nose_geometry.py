"""S414 fore-nose reconstruction alternatives for stage V4; no aerodynamic fitting.

Every curve passes through the same 83 reported fore vertices. Only the
interpolant between them changes, and only inside the declared nose support.
The aft element, registration, reference chord and trailing edges are untouched.
"""
from pathlib import Path
import sys
import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from multifoil import clean_contour, read_dat, repanel, validate_geometry
from s414_geometry import QUAL, elements as declared_elements

NOSE = .06  # declared nose support: reported vertices with x / total chord at or below this value
MARGIN = .05  # each vertex tangent stays at least this fraction of the local turn inside its two chords
VARIANTS = ('nominal', 'single_curve', 'convex_nose')


def source():
    """Reported fore vertices in solver input order (upper TE first) and the branch-split vertex."""
    xy = clean_contour(read_dat(ROOT/QUAL/'report_test_fore.dat'))
    le = int(np.argmax(np.linalg.norm(xy-.5*(xy[0]+xy[-1]), axis=1)))
    return xy, le


def chord_parameter(xy):
    return np.r_[0., np.cumsum(np.linalg.norm(np.diff(xy, axis=0), axis=1))]


def cosine_parameters(s, le, panels):
    """The production sampling of repanel(): cosine spacing in chord arclength on each branch."""
    upper = .5*(1-np.cos(np.linspace(0, np.pi, panels//2+1)))*s[le]
    lower = s[le]+.5*(1-np.cos(np.linspace(0, np.pi, panels-panels//2+1)))*(s[-1]-s[le])
    return np.r_[upper, lower[1:]]


def turns(polygon):
    """Signed exterior turn at each interior vertex, degrees; positive is the convex sense of this contour."""
    d = np.diff(np.asarray(polygon, float), axis=0); angle = np.arctan2(d[:, 1], d[:, 0])
    return np.degrees((np.diff(angle)+np.pi) % (2*np.pi)-np.pi)


def support(xy, le):
    """First and last reported vertex of the nose support (x <= NOSE on both surfaces)."""
    inside = np.flatnonzero(xy[:, 0] <= NOSE)
    if not (inside[0] < le < inside[-1] and np.array_equal(inside, np.arange(inside[0], inside[-1]+1))):
        raise ValueError('Nose support is not one run of reported vertices around the leading edge')
    return int(inside[0]), int(inside[-1])


def circle_tangent_angle(a, b, c, incoming):
    """Tangent angle at b of the circle through three consecutive vertices, continuous with the incoming chord angle."""
    u, v = b-a, c-b; lu, lv = np.linalg.norm(u), np.linalg.norm(v)
    turn = np.arctan2(u[0]*v[1]-u[1]*v[0], u@v)
    return incoming+turn*lv/(lu+lv)


def bending_energy(angles, chords, lengths):
    """Sum of curvature squared times length over all arcs of the convex biarc chain."""
    a = chords-angles[:-1]; b = angles[1:]-chords
    denominator = np.sin(a)*(1+np.cos(b))+np.sin(b)*(1+np.cos(a))
    d1 = lengths*np.sin(b)/denominator; d2 = lengths*np.sin(a)/denominator
    return float(np.sum(a*np.tan(a/2)/d1+b*np.tan(b/2)/d2))


def fair_convex_tangents(xy, i0, i1, first, last):
    """Tangent angles at support vertices: fixed at both ends, elsewhere the minimum-bending-energy
    choice among all convex-consistent ones. Geometry only; no aerodynamic quantity enters."""
    d = np.diff(xy[i0-1:i1+2], axis=0); lengths = np.linalg.norm(d, axis=1)
    chord = np.unwrap(np.arctan2(d[:, 1], d[:, 0]))  # chord[k] joins vertex i0-1+k to i0+k
    turn = np.diff(chord)                             # turn[k] at vertex i0+k
    if not np.all(turn > 0):
        raise ValueError('Reported polygon is not convex in the nose support; convexity is not justified')
    ends = [np.arctan2(t[1], t[0]) for t in (first, last)]
    ends = [chord[k]+(e-chord[k]+np.pi) % (2*np.pi)-np.pi for e, k in zip(ends, (0, len(chord)-2))]
    lower = chord[:-1]+MARGIN*turn; upper = chord[1:]-MARGIN*turn
    if not (lower[0] <= ends[0] <= upper[0] and lower[-1] <= ends[1] <= upper[-1]):
        raise ValueError('Join tangents are not convex-consistent with the reported polygon')
    start = np.array([circle_tangent_angle(xy[i-1], xy[i], xy[i+1], chord[i-i0]) for i in range(i0, i1+1)])
    start = np.clip(start, lower, upper); start[0], start[-1] = ends
    inner = slice(1, -1)
    def objective(free):
        angles = start.copy(); angles[inner] = free
        return bending_energy(angles, chord[1:-1], lengths[1:-1])
    best = minimize(objective, start[inner], method='L-BFGS-B', bounds=list(zip(lower[inner], upper[inner])),
                    options=dict(ftol=1e-15, gtol=1e-12, maxiter=2000))
    angles = start.copy(); angles[inner] = best.x
    return angles, dict(energy_three_point_start=objective(start[inner]), energy_minimum=float(best.fun),
                        iterations=int(best.nit), at_bound=int(np.sum((best.x <= lower[inner]+1e-12) | (best.x >= upper[inner]-1e-12))))


def biarc(p0, t0, p1, t1):
    """Convex pair of circular arcs joining two points with given unit tangents.

    The junction tangent is parallel to the chord, so the first arc turns by the
    angle between the start tangent and the chord and the second by the angle
    between the chord and the end tangent. Both must be positive: the data must
    be convex in the sense of this contour, otherwise the construction refuses.
    Returns a sampler by arclength fraction, the two arc curvatures and the length.
    """
    v = p1-p0; length = np.linalg.norm(v); e = v/length
    a = np.arctan2(t0[0]*e[1]-t0[1]*e[0], t0@e); b = np.arctan2(e[0]*t1[1]-e[1]*t1[0], e@t1)
    if not (a > 0 and b > 0):
        raise ValueError('Vertex tangents are not convex-consistent; convexity is not justified here')
    denominator = np.sin(a)*(1+np.cos(b))+np.sin(b)*(1+np.cos(a))
    d1 = length*np.sin(b)/denominator
    junction = p0+d1*t0+d1*e; arcs = []
    for start_point, tangent, end_point in ((p0, t0, junction), (p1, -t1, junction)):
        chord = end_point-start_point; n = np.array([-tangent[1], tangent[0]]); h = chord@n
        radius = (chord@chord)/(2*h); center = start_point+radius*n
        start = np.arctan2(start_point[1]-center[1], start_point[0]-center[0])
        end = np.arctan2(end_point[1]-center[1], end_point[0]-center[0])
        sweep = (end-start+np.pi) % (2*np.pi)-np.pi
        arcs.append(dict(center=center, radius=abs(radius), start=start, sweep=sweep,
                         length=float(abs(radius*sweep)), curvature=float(np.sign(sweep)/abs(radius))))
    first, second = arcs; total = first['length']+second['length']
    def point(arc, fraction):
        angle = arc['start']+fraction*arc['sweep']
        return arc['center']+arc['radius']*np.array([np.cos(angle), np.sin(angle)])
    def sample(fraction):
        length = fraction*total
        if length <= first['length']: return point(first, length/first['length'])
        # The second arc was built from p1 backwards, so run it in reverse.
        return point(second, 1-(length-first['length'])/second['length'])
    return sample, (first['curvature'], -second['curvature']), total


class Reconstruction:
    """One continuous fore contour through the reported vertices, sampled like the production input."""
    def __init__(self, variant):
        if variant not in VARIANTS:
            raise ValueError('Unknown nose reconstruction')
        self.variant = variant; self.xy, self.le = source(); self.s = chord_parameter(self.xy)
        xy, le, s = self.xy, self.le, self.s
        if variant == 'nominal':
            # Production rule: the two branches are interpolated separately, each in its own
            # arclength from zero, and meet at the split vertex with no tangent condition.
            self.branches = (PchipInterpolator(s[:le+1], xy[:le+1], axis=0), PchipInterpolator(s[le:]-s[le], xy[le:], axis=0))
        else:
            self.whole = PchipInterpolator(s, xy, axis=0)
        if variant == 'convex_nose':
            i0, i1 = support(xy, le); self.support = (i0, i1); self.segments = {}; self.arc_curvatures = {}
            unit = lambda i: self.whole(s[i], 1)/np.linalg.norm(self.whole(s[i], 1))
            # The join tangents are those of the production curve, so the contour outside is untouched.
            angles, self.fairing = fair_convex_tangents(xy, i0, i1, unit(i0), unit(i1))
            self.tangent_angles = dict(zip(range(i0, i1+1), angles))
            tangents = {i: np.array([np.cos(a), np.sin(a)]) for i, a in self.tangent_angles.items()}
            for i in range(i0, i1):
                self.segments[i], self.arc_curvatures[i], _ = biarc(xy[i], tangents[i], xy[i+1], tangents[i+1])

    def point(self, q):
        s, le = self.s, self.le
        if self.variant == 'nominal':
            return self.branches[0](q) if q <= s[le] else self.branches[1](q-s[le])
        if self.variant == 'convex_nose':
            i0, i1 = self.support
            if s[i0] < q < s[i1]:
                i = min(int(np.searchsorted(s, q, side='right'))-1, i1-1)
                return self.segments[i]((q-s[i])/(s[i+1]-s[i]))
        return self.whole(q)

    def changed(self, parameters):
        """Sample parameters whose point this variant replaces; everything else stays the production sample."""
        s = self.s
        if self.variant == 'nominal': return np.zeros(len(parameters), dtype=bool)
        lo, hi = (s[self.le-1], s[self.le+1]) if self.variant == 'single_curve' else (s[self.support[0]], s[self.support[1]])
        return (parameters > lo) & (parameters < hi) & ~np.isin(parameters, s)

    def input(self, panels):
        """The array handed to the solver in place of repanel(source, panels)."""
        points = repanel(self.xy, panels)  # the production function itself, bit for bit
        parameters = cosine_parameters(self.s, self.le, panels)
        for k in np.flatnonzero(self.changed(parameters)):
            points[k] = self.point(float(parameters[k]))
        return points


def fore_input(variant, panels=340):
    return Reconstruction(variant).input(panels)


def elements(variant, panels=(340, 240)):
    """Fore alternative with the unchanged nominal aft input."""
    nominal = declared_elements('nominal', list(panels))
    if variant == 'nominal':
        return nominal
    result = [('fore', fore_input(variant, panels[0])), nominal[1]]
    validate_geometry(result)
    return result
