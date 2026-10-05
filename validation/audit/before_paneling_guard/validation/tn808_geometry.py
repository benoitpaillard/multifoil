"""Reconstruct NACA TN808 slotted flap 1-b from Table I and Figure 1.

Coordinates use the original wing chord. Interpolation and datum interpretations
are explicit; this is a reconstruction, not a claim to the original CAD.
"""
from pathlib import Path
import sys
import numpy as np
from scipy.interpolate import CubicSpline, PchipInterpolator

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from multifoil import validate_geometry
from optimization import minimum_clearance

# NACA TN808 Table I, printed p.10, PDF11. Percent of original wing chord.
MAIN_TABLE = np.array([
    [1.25,2.67,-1.23],[2.5,3.61,-1.71],[5,4.91,-2.26],
    [7.5,5.80,-2.61],[10,6.43,-2.92],[15,7.19,-3.50],
    [20,7.50,-3.97],[25,7.60,-4.28],[30,7.55,-4.46],
    [40,7.14,-4.48],[50,6.41,-4.17],[60,5.47,-3.67],
    [70,4.36,-3.00],[80,3.08,-2.16],[90,1.68,-1.23],
    [95,.92,-.70],[100,.13,-.13]]) / 100
FLAP_TABLE = np.array([
    [0,1.05,1.05],[.375,1.64,.58],[.75,1.96,.44],
    [1.50,2.39,.28],[2.25,2.69,.19],[3.,2.92,.13],
    [4.50,3.27,.05],[6.,3.48,.01],[9.,3.63,0],
    [12.,3.60,0],[15.,3.40,0],[18.,3.05,0],[20.,2.76,0],
    [21.,2.57,0],[24.,1.86,0],[27.,1.08,0],[28.50,.70,0],
    [30.,.25,0]]) / 100
SLOT_TABLE = np.array([
    [69.96,-3.10],[69.96,-1.88],[70.83,-.90],[72.92,.23],
    [77.08,1.25],[80.,1.58],[81.25,1.67],[83.33,1.79],
    [85.42,1.88],[86.88,1.90],[90.,1.59]]) / 100


def naca23012(n=4001):
    """Standard five-digit mean line plus normal four-digit thickness offsets.

    Comparing with Table I is part of geometry_audit(). Keep the finite TE.
    The normal thickness offset recovers the small negative-x nose overhang
    absent from the tabulation; do not close two y(x) curves at the chord origin.
    """
    x = .5*(1-np.cos(np.linspace(0,np.pi,n)))
    m, k = .2025, 15.957
    yc = k/6*np.where(x<m, x**3-3*m*x*x+m*m*(3-m)*x, m**3*(1-x))
    slope = k/6*np.where(x<m, 3*x*x-6*m*x+m*m*(3-m), -m**3)
    theta = np.arctan(slope)
    yt = .6*(.2969*np.sqrt(x)-.126*x-.3516*x*x+.2843*x**3-.1015*x**4)
    upper = np.column_stack([x-yt*np.sin(theta), yc+yt*np.cos(theta)])
    lower = np.column_stack([x+yt*np.sin(theta), yc-yt*np.cos(theta)])
    return upper, lower


def table_corrected_main(n=801):
    upper,lower=naca23012(n)
    dense=naca23012()
    for col,curve,reference in ((1,upper,dense[0]),(2,lower,dense[1])):
        reference=reference[np.argmin(reference[:,0]):]
        correction=MAIN_TABLE[:,col]-np.interp(MAIN_TABLE[:,0],reference[:,0],reference[:,1])
        offset=CubicSpline(np.r_[0,MAIN_TABLE[:,0]],np.r_[0,correction],
                           bc_type=((1,0.),(2,0.)))
        curve[:,1]+=offset(np.maximum(curve[:,0],0.))
    return upper,lower


def main_element(entry_adjustment=0.):
    upper, lower = table_corrected_main(401)
    # Table ordinates anchor the surfaces. Lower surface joins the slot table's
    # slightly displaced first point via an explicitly defined Hermite segment.
    lip = np.array([.9, .0168])
    u = upper[upper[:,0] < .9].copy()
    u = np.vstack([u,lip])
    start_x = .65
    j = np.searchsorted(lower[:,0],start_x)
    l = lower[:j].copy()
    y0 = np.interp(start_x,lower[:,0],lower[:,1])
    slope0 = (lower[j+1,1]-lower[j-1,1])/(lower[j+1,0]-lower[j-1,0])
    slot = SLOT_TABLE.copy()
    slot[0,1] += entry_adjustment
    xend,yend = slot[0]
    # Horizontal tangent at the lower corner; Figure1 has an unrounded step.
    bridge = CubicSpline([start_x,xend],[y0,yend],bc_type=((1,slope0),(1,0.)))
    bx = np.linspace(start_x,xend,60)
    bridge_xy = np.column_stack([bx,bridge(bx)])
    # Monotone coordinate interpolation in arc length preserves every slot row
    # including its vertical entrance without x overshoot.
    s = np.r_[0,np.cumsum(np.linalg.norm(np.diff(slot,axis=0),axis=1))]
    ss = np.linspace(0,s[-1],180)
    inside = PchipInterpolator(s,slot,axis=0)(ss)
    contour = np.vstack([u[::-1],l[1:],bridge_xy,inside[1:]])
    return contour


def flap_element(n=301):
    # Interpolate in sqrt(x), imposing Table I's radius .0045*c at the nose.
    # A simple y(x) spline would create a pointed leading edge.
    length,radius = .30,.0045
    t = np.sqrt(FLAP_TABLE[:,0]/length)
    sample = np.linspace(0,1,n)
    curves=[]
    for col,sign in ((1,1),(2,-1)):
        spline=CubicSpline(t,FLAP_TABLE[:,col],
            bc_type=((1,sign*np.sqrt(2*radius*length)),(2,0.)))
        curves.append(np.column_stack([length*sample**2,spline(sample)]))
    return np.vstack([curves[0][::-1],curves[1][1:]])


def geometry(deflection, gap=.02, datum="nose_to_te", entry_adjustment=0.):
    """Vertical gap at lip x=.9 to the flap upper surface; nose station x=.8975.

    Main lip underside is (.9,.0159). Positive deflection is clockwise.
    Nominal flap chord joins the tabulated nose to the midpoint of its finite
    trailing edge. `table_axis` is a separate datum-sensitivity interpretation.
    Neither interpretation is silently fitted to the experimental polar.
    """
    if datum not in ("nose_to_te","table_axis") or not 0<gap<.1:
        raise ValueError("Invalid datum or gap")
    main = main_element(entry_adjustment)
    flap = flap_element()
    nose=FLAP_TABLE[0,[0,1]]
    te=.5*(flap[0]+flap[-1])
    original_angle=np.arctan2(*(te-nose)[::-1])
    rotation=-np.radians(float(deflection))
    if datum=="nose_to_te":
        rotation-=original_angle
    c,s=np.cos(rotation),np.sin(rotation)
    matrix=np.array([[c,-s],[s,c]])
    flap=(flap-nose)@matrix.T + np.array([.8975,0.])
    # Vertical projection onto upper surface, before any y translation.
    upper=flap[:301][::-1]
    if np.any(np.diff(upper[:,0])<=0):
        raise ValueError("Upper flap x is not monotone for the requested deflection")
    at_lip=float(np.interp(.9,upper[:,0],upper[:,1]))
    flap[:,1]+=.0159-gap-at_lip
    raw=[("main",main),("flap",flap)]
    validate_geometry(raw)
    actual_nose=flap[300].copy()
    actual_te=.5*(flap[0]+flap[-1])
    actual_angle=np.degrees(np.arctan2(*(actual_te-actual_nose)[::-1]))
    return raw,dict(deflection_deg=float(deflection),gap_vertical=float(gap),
        datum=datum,entry_adjustment=float(entry_adjustment),
        nose=actual_nose.tolist(),te=actual_te.tolist(),
        physical_chord_angle_deg=float(actual_angle),
        table_axis_to_chord_deg=float(np.degrees(original_angle)),
        minimum_clearance=minimum_clearance(main,flap))


def geometry_audit():
    upper,lower=naca23012()
    # Use the monotone upper branch away from its nose overhang.
    upper=upper[np.argmin(upper[:,0]):]
    errors=[]
    for col,curve in ((1,upper),(2,lower)):
        errors.extend(np.interp(MAIN_TABLE[:,0],curve[:,0],curve[:,1])-MAIN_TABLE[:,col])
    errors=np.asarray(errors)
    corrected=table_corrected_main(4001)
    remaining=[]
    for col,curve in ((1,corrected[0]),(2,corrected[1])):
        curve=curve[np.argmin(curve[:,0]):]
        remaining.extend(np.interp(MAIN_TABLE[:,0],curve[:,0],curve[:,1])-MAIN_TABLE[:,col])
    flap=flap_element()
    min_thickness=float(np.min(flap[:301][::-1,1]-flap[300:,1]))
    return dict(main_table_max_abs_ordinate_error=float(np.max(np.abs(errors))),
                main_table_rms_ordinate_error=float(np.sqrt(np.mean(errors**2))),
                corrected_main_table_max_abs_error=float(np.max(np.abs(remaining))),
                flap_minimum_input_thickness=min_thickness,
                main_lip=[.9,.0168,.0159],
                note="Analytic NACA nose with smooth table-ordinate corrections; "
                     "table-constrained slot; prescribed .65-.6996 "
                     "lower-surface bridge; sqrt(x) flap interpolation. "
                     "Finite original trailing-edge thicknesses retained.")
