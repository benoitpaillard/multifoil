#!/usr/bin/env python3
"""e^N amplification from Orr-Sommerfeld growth rates along a sequence of profiles (V5 tranche 2).

Each station supplies a velocity profile in units of its own displacement
thickness and edge velocity. For every physical frequency the spatial growth
rate is followed from station to station and integrated from its neutral
point; N is the envelope over frequencies. Independent of the inherited
envelope correlation.
"""
import numpy as np
import s414_stability as stability

WAVENUMBERS = np.geomspace(.03, 8., 44)  # scan band, in units of 1/displacement thickness


def station(xi, delta_star, ue, nu, eta, ratio, curvature, scaled_displacement, grid):
    """Collocation data of one profile from the finite-difference variables."""
    y = eta/scaled_displacement
    u, upp = stability.on_grid(grid, y, ratio, curvature*scaled_displacement**2)
    return dict(xi=float(xi), delta_star=float(delta_star), ue=float(ue), re=float(ue*delta_star/nu), U=u, Upp=upp)


def _accept(alpha, reduced):
    return alpha.real > 0 and .02 < reduced/alpha.real < .98


def growth(stations, omegas, grid, quiet_limit=3, scans=6):
    """Spatial growth rate -alpha_i per unit xi, array [frequency, station]; nan where no wave is found.

    The unstable wave is identified by temporal scans at a few equally spaced
    stations, shared by all frequencies. Each frequency starts at the scanned
    station where it grows fastest and is followed upstream and downstream by
    Newton continuation. Upstream tracking stops when the wave is lost or after
    quiet_limit consecutive damped stations (below its neutral point);
    downstream it continues until the wave is lost. Untracked stations stay nan.
    """
    count = len(stations)
    rate = np.full((len(omegas), count), np.nan)
    wavenumber = np.full((len(omegas), count), np.nan, dtype=complex)
    cache = {}

    def scanned(i):
        if i not in cache:
            item = stations[i]
            cache[i] = stability.scan(item['re'], item['U'], item['Upp'], grid, WAVENUMBERS)
        return cache[i]

    def from_scan(i, omega):
        item = stations[i]; reduced = omega*item['delta_star']/item['ue']
        start, shape = stability.seed(reduced, scanned(i))
        if start is None:
            return None, None
        alpha, phi, ok = stability.spatial(reduced, item['re'], item['U'], item['Upp'], grid, start, shape)
        return (alpha, phi) if ok and _accept(alpha, reduced) else (None, None)

    def continued(i, omega, guess, shape):
        item = stations[i]; reduced = omega*item['delta_star']/item['ue']
        alpha, phi, ok = stability.spatial(reduced, item['re'], item['U'], item['Upp'], grid, guess, shape)
        return (alpha, phi) if ok and _accept(alpha, reduced) and abs(alpha-guess) < .5*abs(guess) else (None, None)

    seeds = sorted({count-1-round(j*(count-1)/max(scans-1, 1)) for j in range(scans)}, reverse=True)
    for k, omega in enumerate(omegas):
        seeded = None
        for i in seeds:
            alpha, phi = from_scan(i, omega)
            local = None if alpha is None else -alpha.imag/stations[i]['delta_star']
            if alpha is not None and (seeded is None or local > seeded[3]):
                seeded = (i, alpha, phi, local)
        if seeded is None or seeded[3] <= 0:
            continue
        i0, alpha0, phi0 = seeded[:3]
        rate[k, i0] = -alpha0.imag/stations[i0]['delta_star']; wavenumber[k, i0] = alpha0
        for direction in (-1, 1):
            alpha, phi, quiet, i = alpha0, phi0, 0, i0+direction
            while 0 <= i < count:
                guess = alpha*stations[i]['delta_star']/stations[i-direction]['delta_star']
                new, shape = continued(i, omega, guess, phi)
                if new is None:
                    break
                alpha, phi = new, shape
                rate[k, i] = -alpha.imag/stations[i]['delta_star']; wavenumber[k, i] = alpha
                quiet = quiet+1 if alpha.imag > 0 else 0
                if direction < 0 and quiet >= quiet_limit:
                    break
                i += direction
    return rate, wavenumber


def envelope(xi, rate):
    """N of each frequency, integrated from its first amplified station, and the envelope over frequencies."""
    xi = np.asarray(xi, dtype=float)
    each = np.zeros_like(rate)
    for k in range(rate.shape[0]):
        value, started = 0., False
        for i in range(1, len(xi)):
            a, b = rate[k, i-1], rate[k, i]
            if not started:
                if np.isfinite(b) and b > 0:
                    started = True
                    # Neutral point by linear interpolation where the upstream station is damped.
                    if np.isfinite(a) and a < 0:
                        x0 = xi[i-1]+(xi[i]-xi[i-1])*(-a)/(b-a)
                        value = .5*b*(xi[i]-x0)
                    elif np.isfinite(a):
                        value = .5*(a+b)*(xi[i]-xi[i-1])
                    else:
                        value = 0.
            elif np.isfinite(a) and np.isfinite(b):
                value += .5*(a+b)*(xi[i]-xi[i-1])
            else:
                value = np.nan  # the wave was lost downstream of its start: excluded from here on
            each[k, i] = value
    return each, np.nanmax(np.maximum(each, 0.), axis=0)
