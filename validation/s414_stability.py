#!/usr/bin/env python3
"""Linear stability of a parallel boundary-layer profile (Orr-Sommerfeld), for stage V5.

Independent of the inherited envelope correlation: Chebyshev collocation of the
Orr-Sommerfeld equation on a mapped semi-infinite domain, a temporal global
solve to identify the unstable wave, and Newton iteration for the spatial
eigenvalue at a real frequency. Lengths and speeds are those of the supplied
profile; the Reynolds number must use the same scales. Non-parallel effects are
neglected, as in every e^N method.
"""
import numpy as np
from scipy.linalg import eig, lu_factor, lu_solve


def chebyshev(n):
    """Gauss-Lobatto nodes s_j = cos(pi j/n) and the first-derivative matrix."""
    j = np.arange(n+1); s = np.cos(np.pi*j/n)
    c = np.where((j == 0) | (j == n), 2., 1.)*(-1.)**j
    ds = s[:, None]-s[None, :]
    d = (c[:, None]/c[None, :])/(ds+np.eye(n+1))
    d -= np.diag(d.sum(axis=1))
    return s, d


class Grid:
    """Collocation points on 0 <= y <= y_max, half of them below y = cluster."""
    def __init__(self, count=140, y_max=60., cluster=3.):
        s, d = chebyshev(count)
        b = 1+2*cluster/y_max
        self.y = cluster*(1+s)/(b-s)
        d1 = d/(cluster*(b+1)/(b-s)**2)[:, None]
        self.d1, self.d2 = d1, d1@d1
        self.d4 = self.d2@self.d2
        self.count, self.y_max, self.cluster = count, y_max, cluster


def operator(alpha, omega, re, velocity, curvature, grid):
    """Orr-Sommerfeld operator L(alpha, omega) and dL/dalpha, with clamped rows at both ends.

    L = (D^2 - a^2)^2 - i Re [(a U - w)(D^2 - a^2) - a U'']
    """
    identity = np.eye(grid.count+1)
    lap = grid.d2-alpha**2*identity
    u, upp = np.diag(velocity), np.diag(curvature)
    op = grid.d4-2*alpha**2*grid.d2+alpha**4*identity-1j*re*((alpha*u-omega*identity)@lap-alpha*upp)
    op_alpha = -4*alpha*lap-1j*re*(u@lap-2*alpha**2*u+2*alpha*omega*identity-upp)
    for row, source in ((0, identity[0]), (1, grid.d1[0]), (-1, identity[-1]), (-2, grid.d1[-1])):
        op[row], op_alpha[row] = source, 0.
    return op, op_alpha


def temporal(alpha, re, velocity, curvature, grid):
    """Most unstable wave for a real wavenumber: complex frequency and eigenfunction.

    Only waves with phase speed between 0 and 1 are considered; returns
    (None, None) if there is none.
    """
    identity = np.eye(grid.count+1)
    lap = grid.d2-alpha**2*identity
    a = grid.d4-2*alpha**2*grid.d2+alpha**4*identity-1j*re*alpha*(np.diag(velocity)@lap-np.diag(curvature))
    b = -1j*re*lap
    for row, source in ((0, identity[0]), (1, grid.d1[0]), (-1, identity[-1]), (-2, grid.d1[-1])):
        a[row], b[row] = source, 0.
    values, vectors = eig(a, b)
    ok = np.isfinite(values) & (np.abs(values) < 50*max(alpha, 1.))
    speed = np.where(ok, values.real/alpha, -1.)
    ok &= (speed > 0.02) & (speed < .98)
    if not np.any(ok):
        return None, None
    k = np.where(ok)[0][np.argmax(values.imag[ok])]
    return values[k], vectors[:, k]


def spatial(omega, re, velocity, curvature, grid, alpha, shape=None, iterations=30, tolerance=1e-10):
    """Newton iteration for the complex wavenumber at a real frequency.

    Returns (alpha, eigenfunction, converged). The eigenfunction is normalized
    to one at the point where the starting shape is largest.
    """
    alpha = complex(alpha)
    if shape is None:
        frequency, shape = temporal(alpha.real, re, velocity, curvature, grid)
        if shape is None:
            return alpha, None, False
    k = int(np.argmax(np.abs(shape))); phi = shape/shape[k]
    n = len(phi)
    for _ in range(iterations):
        op, op_alpha = operator(alpha, omega, re, velocity, curvature, grid)
        matrix = np.zeros((n+1, n+1), dtype=complex)
        matrix[:n, :n], matrix[:n, n] = op, op_alpha@phi
        matrix[n, k] = 1.
        right = np.r_[-op@phi, 0.]
        try:
            step = lu_solve(lu_factor(matrix), right)
        except (ValueError, np.linalg.LinAlgError):
            return alpha, phi, False
        if not np.all(np.isfinite(step)):
            return alpha, phi, False
        phi = phi+step[:n]; alpha = alpha+step[n]
        if abs(step[n]) < tolerance*max(1., abs(alpha)):
            return alpha, phi, True
    return alpha, phi, False


def scan(re, velocity, curvature, grid, wavenumbers):
    """Temporal solutions over a band of real wavenumbers: list of (alpha, omega, eigenfunction)."""
    found = [(a,)+temporal(a, re, velocity, curvature, grid) for a in wavenumbers]
    return [(a, w, v) for a, w, v in found if w is not None]


def seed(omega, found):
    """Starting wavenumber for a real frequency from a temporal scan (Gaster's relation).

    Returns (alpha_guess, shape) or (None, None) when no travelling wave with
    that real frequency exists in the scanned band. Where several wavenumbers
    match, the one with the largest temporal growth is taken.
    """
    best = None
    for (a0, w0, v0), (a1, w1, v1) in zip(found[:-1], found[1:]):
        if (w0.real-omega)*(w1.real-omega) <= 0 and w1.real != w0.real:
            t = (omega-w0.real)/(w1.real-w0.real)
            group = (w1.real-w0.real)/(a1-a0)
            growth = w0.imag+t*(w1.imag-w0.imag)
            if group > 0 and (best is None or growth > best[0]):
                best = (growth, complex(a0+t*(a1-a0), -growth/group), v0 if t < .5 else v1)
    return (None, None) if best is None else best[1:]


def on_grid(grid, y, velocity, curvature):
    """Profile and curvature on the collocation points; uniform flow beyond the data."""
    inside = grid.y <= y[-1]
    u = np.where(inside, np.interp(grid.y, y, velocity), 1.)
    upp = np.where(inside, np.interp(grid.y, y, curvature), 0.)
    return u, upp
