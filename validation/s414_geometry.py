"""Declared S414 reconstruction sensitivities; no aerodynamic fitting."""
from pathlib import Path
import sys
import numpy as np
from scipy.interpolate import CubicSpline

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from multifoil import clean_contour, read_dat, repanel, validate_geometry

QUAL = 'validation/results/s414_geometry_qualification_20261002'
SHIFT = .13 / 457.2
SCENARIOS = ('nominal', 'source_cubic', 'aft_z_minus', 'aft_z_plus')


def cubic_repanel(xy, panels):
    """Only replace the source interpolator; preserve the production sampling."""
    xy = clean_contour(xy)
    if not isinstance(panels, int) or not 40 <= panels <= 2000:
        raise ValueError('Invalid panel count')
    le = int(np.argmax(np.linalg.norm(xy - .5*(xy[0]+xy[-1]), axis=1)))
    if not 3 <= le <= len(xy)-4:
        raise ValueError('Unsupported leading-edge split')
    def branch(pts, count):
        s = np.r_[0., np.cumsum(np.linalg.norm(np.diff(pts, axis=0), axis=1))]
        q = .5*(1-np.cos(np.linspace(0, np.pi, count+1))) * s[-1]
        sampled = CubicSpline(s, pts, axis=0, bc_type='natural')(q)
        sampled[0], sampled[-1] = pts[0], pts[-1]
        return sampled
    return np.concatenate((branch(xy[:le+1], panels//2), branch(xy[le:], panels-panels//2)[1:]))


def elements(scenario, panels=(340, 240)):
    if scenario not in SCENARIOS or len(panels) != 2:
        raise ValueError('Unknown S414 geometry scenario')
    result = []
    for name, count in zip(('fore', 'aft'), panels):
        source = read_dat(ROOT/QUAL/f'report_test_{name}.dat')
        xy = cubic_repanel(source, count) if scenario == 'source_cubic' else repanel(source, count)
        if name == 'aft' and scenario in ('aft_z_minus', 'aft_z_plus'):
            xy[:, 1] += SHIFT * (-1 if scenario == 'aft_z_minus' else 1)
        result.append((name, xy))
    validate_geometry(result)
    return result
