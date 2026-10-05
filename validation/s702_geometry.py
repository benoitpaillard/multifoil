"""Published S702 vector-figure reconstructions in the original common frame."""
from pathlib import Path
import sys
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from multifoil import read_dat, validate_geometry

DATA = ROOT / 'validation/data/s702'
FILES = {-2: 'aft_dm2.dat', 0: 'aft_d0.dat', 5: 'aft_d5.dat', 10: 'aft_d10.dat'}


def geometry(deflection=0):
    """Return fore and aft outlines; only the small aft-element flap moves.

    No element-wise normalization, incidence adjustment or drag fitting is done.
    Only the four published configurations are accepted.
    """
    if deflection not in FILES:
        raise ValueError('Available published flap deflections: -2, 0, 5, 10 degrees')
    elements = [('fore', read_dat(DATA/'fore.dat')),
                ('aft', read_dat(DATA/FILES[deflection]))]
    validate_geometry(elements)
    return elements
