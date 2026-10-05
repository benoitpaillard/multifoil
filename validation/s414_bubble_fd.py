#!/usr/bin/env python3
"""Finite-difference laminar boundary layer through separation, for V5 tranche 2.

Extends the tranche-1 solver (which is left untouched) with an inverse mode:
the displacement thickness is prescribed and the edge velocity is the unknown,
which is well posed through a separation bubble. Reversed flow is marched with
the FLARE approximation (streamwise convection u du/dx dropped where u < 0).

Variables: eta = y sqrt(u_ref/(nu xi)), u/u_ref = f', with a known reference
speed u_ref(xi) and the unknown ratio w = ue/u_ref. With m = dln(u_ref)/dln(xi),
k = (m-1)/2 and C = 1 (0 in reversed flow under FLARE):

    f''' + w dw/dln(xi) + m w^2 - C f' [df'/dln(xi) + m f' + k eta f'']
         + [df/dln(xi) + (m+1)/2 f + k eta f'] f'' = 0

Direct mode fixes w = 1; inverse mode adds eta_e - f_e/w = scaled displacement
thickness. Streamwise derivatives are second-order backward differences.
"""
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu
import s414_laminar_fd as fd


def eta_grid(eta_max=30., count=1201):
    return np.linspace(0., eta_max, count)


def _newton(eta, residual_and_jacobian, f, u, v, extra, iterations=60, tolerance=1e-10):
    """Shared Newton loop for the bordered box system; extra is the scalar unknown."""
    for iteration in range(iterations):
        residual, jacobian = residual_and_jacobian(f, u, v, extra)
        norm = float(np.max(np.abs(residual)))
        if not np.isfinite(norm):
            return None
        if norm < tolerance:
            return f, u, v, extra, norm, iteration
        try:
            step = splu(jacobian).solve(-residual)
        except RuntimeError:
            return None
        if not np.all(np.isfinite(step)):
            return None
        # Limit the change of the scalar unknown; the profile follows.
        limit = min(1., .2/max(abs(step[-1]), 1e-30)) if abs(step[-1]) > .2 else 1.
        f = f+limit*step[0:-1:3]; u = u+limit*step[1:-1:3]; v = v+limit*step[2:-1:3]; extra = extra+limit*step[-1]
    return None


def _assemble(eta, r3, d_f, d_u, d_v1, d_v0, d_extra, f, u, v, edge_row, closure_row):
    """Bordered box system: 3n profile unknowns and one scalar."""
    n = len(eta); h = np.diff(eta)
    j = np.arange(1, n)
    rows1, rows2, rows3 = 3*j-1, 3*j, 3*j+1
    cf0, cu0, cv0, cf1, cu1, cv1 = 3*(j-1), 3*(j-1)+1, 3*(j-1)+2, 3*j, 3*j+1, 3*j+2
    one = np.ones(n-1); last = 3*n
    residual = np.empty(3*n+1)
    residual[0], residual[1] = f[0], u[0]
    residual[rows1] = f[1:]-f[:-1]-.5*h*(u[1:]+u[:-1])
    residual[rows2] = u[1:]-u[:-1]-.5*h*(v[1:]+v[:-1])
    residual[rows3] = r3
    residual[3*n-1], edge_columns, edge_values = edge_row
    residual[last], closure_columns, closure_values = closure_row
    rows = np.concatenate(([0, 1], rows1, rows1, rows1, rows1, rows2, rows2, rows2, rows2,
                           rows3, rows3, rows3, rows3, rows3, rows3, rows3,
                           np.full(len(edge_columns), 3*n-1), np.full(len(closure_columns), last)))
    cols = np.concatenate(([0, 1], cf1, cf0, cu1, cu0, cu1, cu0, cv1, cv0, cf1, cf0, cu1, cu0, cv1, cv0,
                           np.full(n-1, last), edge_columns, closure_columns))
    data = np.concatenate(([1., 1.], one, -one, -.5*h, -.5*h, one, -one, -.5*h, -.5*h,
                           d_f, d_f, d_u, d_u, d_v1, d_v0, d_extra, edge_values, closure_values))
    return residual, coo_matrix((data, (rows, cols)), shape=(3*n+1, 3*n+1)).tocsc()


def similar_profile(displacement, guess=None, exponent=0., eta=None, flare=False):
    """Falkner-Skan profile with a prescribed scaled displacement thickness.

    The exponent m of ue ~ xi**m is the unknown, which follows the family
    through separation onto the reversed-flow branch. Returns (m, f, u, v).
    With flare=True the reversed part obeys the FLARE-modified equation; that
    quantifies the approximation and is not a Falkner-Skan solution.
    """
    eta = eta_grid() if eta is None else eta
    n = len(eta); h = np.diff(eta); mid = .5*(eta[1:]+eta[:-1])
    if guess is None:
        guess = fd.similarity(exponent, eta)
    f, u, v = (a.copy() for a in guess)

    def system(f, u, v, m):
        fb, ub, vb = .5*(f[1:]+f[:-1]), .5*(u[1:]+u[:-1]), .5*(v[1:]+v[:-1])
        a1, k = .5*(m+1), .5*(m-1)
        forward = np.ones(n-1) if not flare else (ub >= 0).astype(float)
        dv = (v[1:]-v[:-1])/h
        r3 = dv+m+np.where(forward > 0, a1*fb*vb-m*ub**2, (a1*fb+k*mid*ub)*vb)
        coefficient = np.where(forward > 0, a1*fb, a1*fb+k*mid*ub)
        d_f = .5*a1*vb
        d_u = np.where(forward > 0, -m*ub, .5*k*mid*vb)
        d_m = 1+np.where(forward > 0, .5*fb*vb-ub**2, .5*fb*vb+.5*mid*ub*vb)
        return _assemble(eta, r3, d_f, d_u, 1/h+.5*coefficient, -1/h+.5*coefficient, d_m, f, u, v,
                         (u[-1]-1, np.array([3*n-2]), np.array([1.])),
                         (eta[-1]-f[-1]-displacement, np.array([3*n-3]), np.array([-1.])))
    result = _newton(eta, system, f, u, v, exponent)
    if result is None:
        raise ValueError('Similar profile did not converge')
    f, u, v, m = result[:4]
    return float(m), f, u, v


def family(displacements, eta=None, flare=False):
    """Continuation of the similar family in the scaled displacement thickness, starting from Blasius."""
    eta = eta_grid() if eta is None else eta
    profile = fd.similarity(0., eta); m = 0.; out = []
    start = fd.measures(eta, *profile)[0]
    wanted = np.sort(np.asarray(displacements, dtype=float))
    for value in np.unique(np.r_[np.arange(start, wanted[-1], .05)[1:], wanted]):
        m, *profile = similar_profile(value, profile, m, eta, flare)
        if np.any(wanted == value):
            out.append((float(value), m, tuple(a.copy() for a in profile)))
    return out


def measures(eta, f, u, v, w=1.):
    """Scaled displacement and momentum integrals, wall shear and strongest reversed flow."""
    ratio = u/w
    d1 = float(eta[-1]-f[-1]/w)
    integrand = ratio*(1-ratio)
    d2 = float(np.sum(.5*(integrand[1:]+integrand[:-1])*np.diff(eta)))
    return d1, d2, float(v[0]), float(ratio.min())


def march(xi_nodes, reference, nu, switch=None, displacement=None, eta=None, flare=True, maximum_depth=10,
          keep_profiles=False):
    """March from a stagnation similarity start; direct up to xi = switch, inverse beyond.

    reference(xi) is the edge velocity in direct mode and the scaling speed in
    inverse mode; displacement(xi) is the prescribed displacement thickness.
    If displacement has a join(xi, value) method it is called once, with the
    solver's own displacement thickness at the last direct station, before the
    first inverse step, so that the prescription can start from that value.
    A failed step is halved; if no step succeeds the march stops with status
    'failed' (inverse) or 'separated' (direct, as in tranche 1).
    """
    eta = eta_grid() if eta is None else eta
    n = len(eta); h = np.diff(eta); mid = .5*(eta[1:]+eta[:-1])
    xi_nodes = np.asarray(xi_nodes, dtype=float)
    f, u, v = fd.similarity(1., eta)
    xi = float(xi_nodes[0]); w = 1.
    records, profiles = [], []

    def store(x, f, u, v, w, mode, newton):
        speed = float(reference(x)); scale = np.sqrt(nu*x/speed)
        d1, d2, shear, reverse = measures(eta, f, u, v, w)
        records.append((x, w*speed, w, d1*scale, d2*scale, d1/d2, 2*shear/(w**2*np.sqrt(speed*x/nu)), 2*shear*d2/w, reverse,
                        abs(float(v[-1])), 1. if mode == 'inverse' else 0., newton))
        if keep_profiles:
            profiles.append(dict(xi=x, eta=eta.copy(), ratio=u/w, curvature=np.gradient(v, eta)/w, scaled_displacement=d1))

    store(xi, f, u, v, w, 'direct', 0)
    older, older_step = (f, u, w), None
    status, reason, joined = 'completed', '', [False]

    def attempt(trial):
        speed, back = float(reference(trial)), trial*(1-1e-6)
        m = np.log(speed/float(reference(back)))/np.log(trial/back)
        a1, k = .5*(m+1), .5*(m-1)
        step_log = np.log(trial/xi)
        ratio = 1. if older_step is None else step_log/older_step
        c0, c1, c2 = ((1+2*ratio)/((1+ratio)*step_log), -(1+ratio)/step_log, ratio**2/((1+ratio)*step_log)) \
            if ratio <= 2 else (1/step_log, -1/step_log, 0.)
        hf, hu, hw = c1*f+c2*older[0], c1*u+c2*older[1], c1*w+c2*older[2]
        hfb, hub = .5*(hf[1:]+hf[:-1]), .5*(hu[1:]+hu[:-1])
        inverse = switch is not None and trial > switch
        if inverse and not joined[0]:
            if hasattr(displacement, 'join'):
                displacement.join(xi, records[-1][3])
            joined[0] = True
        target = float(displacement(trial))*np.sqrt(speed/(nu*trial)) if inverse else None

        def system(f, u, v, w):
            fb, ub, vb = .5*(f[1:]+f[:-1]), .5*(u[1:]+u[:-1]), .5*(v[1:]+v[:-1])
            df, du, dw = c0*fb+hfb, c0*ub+hub, c0*w+hw
            forward = (ub >= 0) if flare else np.ones(n-1, dtype=bool)
            r3 = (v[1:]-v[:-1])/h + w*dw + m*w**2 + np.where(
                forward, a1*fb*vb-m*ub**2-ub*du+vb*df, (df+a1*fb+k*mid*ub)*vb)
            coefficient = np.where(forward, a1*fb+df, df+a1*fb+k*mid*ub)
            d_f = .5*(a1+c0)*vb
            d_u = np.where(forward, -m*ub-.5*(du+c0*ub), .5*k*mid*vb)
            d_w = np.full(n-1, dw+c0*w+2*m*w)
            closure = (f[-1]-w*(eta[-1]-target), np.array([3*n-3, 3*n]), np.array([1., -(eta[-1]-target)])) if inverse \
                else (w-1., np.array([3*n]), np.array([1.]))
            return _assemble(eta, r3, d_f, d_u, 1/h+.5*coefficient, -1/h+.5*coefficient, d_w, f, u, v,
                             (u[-1]-w, np.array([3*n-2, 3*n]), np.array([1., -1.])), closure)
        result = _newton(eta, system, f, u, v, w)
        if result is None:
            return None, 'Newton failure', inverse
        if not inverse and result[2][0] <= 0:
            return None, 'non-positive wall shear', inverse
        if result[3] <= 0:
            return None, 'non-positive edge velocity', inverse
        return result, '', inverse

    for target_node in xi_nodes[1:]:
        smallest, failed = (target_node-xi)*2.**(-maximum_depth), None
        while xi < target_node:
            final = failed is not None and failed-xi < smallest
            trial = target_node if failed is None else (failed if final else .5*(xi+failed))
            result, message, inverse = attempt(trial)
            if result is None:
                reason = message
                if final:
                    status = 'failed' if inverse else 'separated'; break
                failed = trial
                continue
            older, older_step = (f, u, w), float(np.log(trial/xi))
            f, u, v, w = result[:4]; xi = float(trial)
            store(xi, f, u, v, w, 'inverse' if inverse else 'direct', result[5])
            if final:
                failed = None
        if status != 'completed':
            break
    table = np.array(records)
    out = dict(status=status, stop_reason=reason if status != 'completed' else '', xi=table[:, 0], ue=table[:, 1],
               ratio=table[:, 2], delta_star=table[:, 3], theta=table[:, 4], H=table[:, 5], cf=table[:, 6],
               cf_re_theta=table[:, 7], reverse=table[:, 8], edge_gradient=table[:, 9], inverse=table[:, 10] > 0,
               newton=table[:, 11].astype(int), eta_max=float(eta[-1]), eta_count=n)
    if keep_profiles:
        out['profiles'] = profiles
    return out
