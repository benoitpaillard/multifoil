#!/usr/bin/env python3
"""Post-hoc V3 diagnostic: does the between-method lift offset follow the trailing-edge base?

Added after the S414 comparison showed a resolution-independent circulation offset
and its largest local differences on the trailing-edge panels. Symmetric NACA 0012
control only; the two methods model a blunt base differently (solid base panel in
Hess–Smith, source/vortex wake panel in the production operator). Linear solves only.
"""
import os
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
from pathlib import Path
import numpy as np
import report_s414_formulation_check as check

OUT = check.OUT
FACTORS, PANELS, ALPHA = (0., .5, 1., 2., 4.), (200, 400), 4.


def naca0012(panels, factor):
    """factor 0 closes the trailing edge; 1 is the standard open section; larger widens the base."""
    x = .5*(1-np.cos(np.linspace(0, np.pi, panels//2+1)))
    y = .6*(.2969*np.sqrt(x)-.126*x-.3516*x**2+.2843*x**3-(.1036-.0021*factor)*x**4)
    if factor == 0: y[-1] = 0.
    return np.concatenate((np.column_stack((x, y))[::-1], np.column_stack((x, -y))[1:]))


def main():
    check.require(not (OUT/'te_treatment_diagnostic.csv').exists(), 'Preserve completed diagnostic')
    rows = []
    for factor in FACTORS:
        for panels in PANELS:
            xy = naca0012(panels, factor); s = check.Pair([('naca', check.clockwise(xy))]).solve(ALPHA)
            e = s['elements']['naca']; inner = e['xm'][0] < .9
            rows.append(dict(te_factor=factor, te_gap_over_chord=float(np.linalg.norm(xy[0]-xy[-1])), panels=panels,
                alpha_deg=ALPHA, production_CL_circulation=s['production_far_CL'], hess_CL_circulation=s['hess_far_CL'],
                CL_circulation_difference=s['hess_far_CL']-s['production_far_CL'],
                relative_CL_difference=s['hess_far_CL']/s['production_far_CL']-1,
                production_far_flux=s['production_far_flux'], hess_far_flux=s['hess_far_flux'],
                max_speed_difference_ahead_of_90_percent=float(np.max(abs(e['hess'][inner]-e['production'][inner]))),
                max_speed_difference_last_panels=float(np.max(abs(e['hess'][~inner]-e['production'][~inner]))),
                linear_residual=max(s['production_linear_residual'], s['hess_linear_residual'])))
            print(rows[-1], flush=True)
    check.csv_write(OUT/'te_treatment_diagnostic.csv', rows)


if __name__ == '__main__':
    main()
