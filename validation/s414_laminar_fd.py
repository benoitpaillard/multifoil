#!/usr/bin/env python3
"""Independent laminar boundary-layer calculations for stage V5.

Two methods that share no formula with the vendored integral closure:

* a finite-difference solution of the incompressible laminar boundary-layer
  equations in Falkner-Skan variables (box discretization across the layer,
  second-order backward differences along it), marched from a similarity start
  to the end of the domain or to laminar separation (zero wall shear). The
  centred streamwise differencing of the original Keller box is kept as
  scheme='box': it is undamped and gives a sawtooth wall shear when the edge
  velocity has a kink, so it is not used for results;
* Thwaites' one-parameter integral method with the Cebeci-Bradshaw fits.

Both take the edge velocity as the only flow input. Direct marching stops at
separation: nothing downstream of zero wall shear is computed or implied.
"""
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu

SEPARATION_LAMBDA = -0.09  # Thwaites' separation value of theta^2 ue'/nu


def eta_grid(eta_max=14., count=561):
    return np.linspace(0., eta_max, count)


def _solve(eta, m, alpha, previous, guess, iterations=40, tolerance=1e-11, backward=None):
    """Newton solution of one station.

    previous is None for a similarity station (no streamwise term). With
    backward=None the momentum equation is centred between the previous station
    and this one (Keller box), with exponent m = dln(ue)/dln(xi) and
    alpha = 1/dln(xi) of the step. With backward=(c0, history_f, history_u) the
    equation is collocated at this station and d/dln(xi) of f and u is
    c0*value + history (backward differences); m is then the local exponent.
    """
    n = len(eta); h = np.diff(eta)
    f, u, v = (a.copy() for a in guess)
    a1 = .5*(m+1)
    if backward is not None:
        c0, history_f, history_u = backward
        hfb, hub = .5*(history_f[1:]+history_f[:-1]), .5*(history_u[1:]+history_u[:-1])
    elif previous is not None:
        fp, up, vp = previous
        fbp, ubp, vbp = .5*(fp[1:]+fp[:-1]), .5*(up[1:]+up[:-1]), .5*(vp[1:]+vp[:-1])
        dvp = (vp[1:]-vp[:-1])/h
    j = np.arange(1, n)
    rows1, rows2, rows3 = 3*j-1, 3*j, 3*j+1
    cf0, cu0, cv0 = 3*(j-1), 3*(j-1)+1, 3*(j-1)+2
    cf1, cu1, cv1 = 3*j, 3*j+1, 3*j+2
    one = np.ones(n-1)
    for iteration in range(iterations):
        fb, ub, vb = .5*(f[1:]+f[:-1]), .5*(u[1:]+u[:-1]), .5*(v[1:]+v[:-1])
        r1 = f[1:]-f[:-1]-h*ub
        r2 = u[1:]-u[:-1]-h*vb
        if backward is not None:
            df, du = c0*fb+hfb, c0*ub+hub
            r3 = (v[1:]-v[:-1])/h + a1*fb*vb + m*(1-ub**2) - ub*du + vb*df
            d_v1, d_v0 = 1/h+.5*a1*fb+.5*df, -1/h+.5*a1*fb+.5*df
            d_f = .5*a1*vb + .5*c0*vb
            d_u = -m*ub - .5*(du+c0*ub)
        elif previous is None:
            r3 = (v[1:]-v[:-1])/h + a1*fb*vb + m*(1-ub**2)
            d_v1, d_v0 = 1/h+.5*a1*fb, -1/h+.5*a1*fb
            d_f = .5*a1*vb
            d_u = -m*ub
        else:
            r3 = (.5*((v[1:]-v[:-1])/h+dvp) + .5*a1*(fb*vb+fbp*vbp) + m*(1-.5*(ub**2+ubp**2))
                  - alpha*(.5*(ub**2-ubp**2) - .5*(vb+vbp)*(fb-fbp)))
            common = .25*a1*fb + .25*alpha*(fb-fbp)
            d_v1, d_v0 = .5/h+common, -.5/h+common
            d_f = .25*a1*vb + .25*alpha*(vb+vbp)
            d_u = -.5*m*ub - .5*alpha*ub
        residual = np.empty(3*n)
        residual[0], residual[1], residual[-1] = f[0], u[0], u[-1]-1
        residual[rows1], residual[rows2], residual[rows3] = r1, r2, r3
        norm = float(np.max(np.abs(residual)))
        if not np.isfinite(norm):
            return None
        if norm < tolerance:
            return f, u, v, norm, iteration
        rows = np.concatenate(([0, 1, 3*n-1], rows1, rows1, rows1, rows1, rows2, rows2, rows2, rows2,
                               rows3, rows3, rows3, rows3, rows3, rows3))
        cols = np.concatenate(([0, 1, 3*n-2], cf1, cf0, cu1, cu0, cu1, cu0, cv1, cv0,
                               cf1, cf0, cu1, cu0, cv1, cv0))
        data = np.concatenate(([1., 1., 1.], one, -one, -.5*h, -.5*h, one, -one, -.5*h, -.5*h,
                               d_f, d_f, d_u, d_u, d_v1, d_v0))
        jacobian = coo_matrix((data, (rows, cols)), shape=(3*n, 3*n)).tocsc()
        try:
            step = splu(jacobian).solve(-residual)
        except RuntimeError:
            return None
        if not np.all(np.isfinite(step)):
            return None
        f += step[0::3]; u += step[1::3]; v += step[2::3]
    return None


def similarity(m, eta, guess=None):
    """Falkner-Skan similarity profile for ue ~ xi**m (attached branch)."""
    if guess is None:
        u = 1-np.exp(-eta)
        guess = (eta-u, u, np.exp(-eta))
    result = _solve(eta, m, 0., None, guess)
    if result is None:
        raise ValueError('Similarity profile did not converge')
    return result[:3]


def measures(eta, f, u, v):
    """Similarity-scaled displacement and momentum integrals and wall shear."""
    d1 = float(eta[-1]-f[-1])
    integrand = u*(1-u)
    d2 = float(np.sum(.5*(integrand[1:]+integrand[:-1])*np.diff(eta)))
    return d1, d2, float(v[0])


def march(xi_nodes, ue_function, nu, m_start=1., eta=None, maximum_depth=12, edge_tolerance=1e-6,
          keep_profiles=False, scheme='bdf2'):
    """March from a similarity start at xi_nodes[0] through the requested nodes.

    scheme='bdf2' (default): second-order backward differences in ln(xi), first
    order where the step grows by more than a factor two; the pressure-gradient
    exponent is the upstream-sided local value, so a kink at a node is taken
    from the interval just marched. scheme='box': centred Keller box.

    A step that fails, or that returns non-positive wall shear, is halved up to
    maximum_depth times. If no admissible step remains the march stops and the
    separation position is estimated from the square-root (Goldstein) behaviour
    of the wall shear over the last accepted stations. Returns a dictionary;
    status is 'completed' or 'separated'.
    """
    eta = eta_grid() if eta is None else eta
    xi_nodes = np.asarray(xi_nodes, dtype=float)
    assert np.all(np.diff(xi_nodes) > 0) and xi_nodes[0] > 0
    f, u, v = similarity(m_start, eta)
    xi = float(xi_nodes[0]); ue = float(ue_function(xi))
    records = []; profiles = []

    def store(x, speed, f, u, v, newton):
        d1, d2, shear = measures(eta, f, u, v)
        scale = np.sqrt(nu*x/speed)
        records.append((x, speed, d1*scale, d2*scale, d1/d2, 2*shear/np.sqrt(speed*x/nu), shear, d2,
                        abs(float(v[-1])), newton))
        if keep_profiles:
            profiles.append(u.copy())

    store(xi, ue, f, u, v, 0)
    status, reason = 'completed', ''
    # Upstream of the first node the solution is the similarity profile, constant in ln(xi).
    older, older_step = (f, u), None

    def attempt(trial):
        speed = float(ue_function(trial))
        if speed <= 0:
            return None, speed, 'non-positive edge velocity'
        step_log = np.log(trial/xi)
        if scheme == 'box':
            result = _solve(eta, np.log(speed/ue)/step_log, 1/step_log, (f, u, v), (f, u, v))
        else:
            back = trial*(1-1e-6)
            exponent = np.log(speed/float(ue_function(back)))/np.log(trial/back)
            ratio = 1. if older_step is None else step_log/older_step
            if ratio <= 2:
                c0, c1, c2 = (1+2*ratio)/((1+ratio)*step_log), -(1+ratio)/step_log, ratio**2/((1+ratio)*step_log)
            else:
                c0, c1, c2 = 1/step_log, -1/step_log, 0.
            result = _solve(eta, exponent, 0., (f, u, v), (f, u, v),
                            backward=(c0, c1*f+c2*older[0], c1*u+c2*older[1]))
        if result is None:
            return None, speed, 'Newton failure'
        if result[2][0] <= 0:
            return None, speed, 'non-positive wall shear'
        return result, speed, ''

    for target in xi_nodes[1:]:
        smallest, failed = (target-xi)*2.**(-maximum_depth), None
        while xi < target:
            # Bisect towards the nearest failed position; a failure reached by a
            # vanishing step is the end of the attached solution.
            final = failed is not None and failed-xi < smallest
            trial = target if failed is None else (failed if final else .5*(xi+failed))
            result, speed, message = attempt(trial)
            if result is None:
                reason = message
                if final:
                    status = 'separated'; break
                failed = trial
                continue
            older, older_step = (f, u), float(np.log(trial/xi))
            f, u, v = result[:3]; xi, ue = float(trial), speed
            store(xi, ue, f, u, v, result[4])
            if final:
                failed = None
        if status == 'separated':
            break
    table = np.array(records)
    out = dict(status=status, stop_reason=reason if status == 'separated' else '', xi=table[:, 0], ue=table[:, 1],
               delta_star=table[:, 2], theta=table[:, 3], H=table[:, 4], cf=table[:, 5], wall_shear=table[:, 6],
               theta_scaled=table[:, 7], edge_gradient=table[:, 8], newton=table[:, 9].astype(int),
               eta_max=float(eta[-1]), eta_count=len(eta), xi_separation=None)
    out['edge_resolved'] = bool(np.max(table[:, 8]) < edge_tolerance)
    if status == 'separated':
        out['xi_separation'] = goldstein(table[:, 0], table[:, 6])
    if keep_profiles:
        out['profiles'] = np.array(profiles); out['eta'] = eta
    return out


def oscillation(values):
    """Sawtooth index: median |second difference| over the median level; ~0 for a smooth record."""
    values = np.asarray(values, dtype=float)
    if len(values) < 5:
        return 0.
    return float(np.median(np.abs(values[2:]-2*values[1:-1]+values[:-2]))/np.median(np.abs(values)))


def goldstein(xi, shear, points=4):
    """Zero of a straight line through shear**2 over the last accepted stations."""
    xi, shear = np.asarray(xi)[-points:], np.asarray(shear)[-points:]
    if len(xi) < 2:
        return float(xi[-1])
    slope, intercept = np.polyfit(xi, shear**2, 1)
    if slope >= 0:
        return float(xi[-1])
    estimate = -intercept/slope
    # The estimate may not precede the last attached station.
    return float(max(estimate, xi[-1]))


def thwaites(xi, ue, due, nu, k_stagnation=None):
    """Thwaites' method on the given nodes; returns theta, lambda, H, l and separation.

    xi[0] is treated as the end of a stagnation region ue = K xi when
    k_stagnation is given, otherwise the momentum thickness starts from zero.
    """
    xi, ue, due = (np.asarray(a, dtype=float) for a in (xi, ue, due))
    integrand = ue**5
    integral = np.r_[0., np.cumsum(.5*(integrand[1:]+integrand[:-1])*np.diff(xi))]
    if k_stagnation is not None:
        integral += ue[0]**5*xi[0]/6
    theta = np.sqrt(.45*nu*integral/ue**6)
    if k_stagnation is None:
        theta[0] = 0.
    lam = theta**2*due/nu
    clipped = np.clip(lam, -.1, .1)
    H = np.where(clipped >= 0, 2.61-3.75*clipped+5.24*clipped**2, 2.088+.0731/(clipped+.14))
    l = np.where(clipped >= 0, .22+1.57*clipped-1.8*clipped**2, .22+1.402*clipped+.018*clipped/(clipped+.107))
    separation = None
    below = np.where(lam <= SEPARATION_LAMBDA)[0]
    if len(below) and below[0] > 0:
        i = below[0]
        w = (SEPARATION_LAMBDA-lam[i-1])/(lam[i]-lam[i-1])
        separation = float(xi[i-1]+w*(xi[i]-xi[i-1]))
    elif len(below):
        separation = float(xi[0])
    return dict(theta=theta, lam=lam, H=H, l=l, xi_separation=separation)
