"""Recover monotone panel locations when mfoil's inverse cubic spline folds.

The fallback reuses the curvature measure from MIT-licensed vendor/mfoil.py
(Fidkowski). Only the inverse cumulative-spacing interpolator is changed.
Valid upstream results are returned unchanged, including their tangent arrays.
"""
import numpy as np
from scipy.interpolate import PchipInterpolator
from vendor import mfoil as mf


def curvature_paneling(points, count, uniformity=2., te_factor=.1):
    result=mf.spline_curvature(points,count,uniformity,te_factor,None)
    if all(np.isfinite(a).all() for a in result) and np.all(np.diff(result[1])>0):
        return (*result,False)
    # A positive curvature measure has a monotone cumulative integral. Cubic
    # interpolation of its inverse need not be monotone near sharp curvature
    # changes. Reconstruct the same measure and invert with shape preservation.
    coarse=mf.spline2d(points)
    dense_s=np.linspace(0,coarse['X'].x[-1],501)
    dense=mf.splineval(coarse,dense_s)
    spline=mf.spline2d(dense)
    s=spline['X'].x
    weight=np.zeros(len(s))
    q,w=mf.quadseg()
    xmin,xmax=np.min(points[0]),np.max(points[0])
    for i in range(len(s)-1):
        ds=s[i+1]-s[i]
        t=q*ds
        px,py=spline['X'].c[:,i],spline['Y'].c[:,i]
        xss=6*px[0]*t+2*px[1]
        yss=6*py[0]*t+2*py[1]
        integral=.01*uniformity+.5*np.dot(w,np.sqrt(xss*xss+yss*yss))*ds
        x=(.5*(dense[0,i]+dense[0,i+1])-xmin)/(xmax-xmin)
        integral+=te_factor*.5*np.exp(-100*(1-x))
        weight[i+1]=weight[i]+integral
    weight+=2*np.sum(weight)/len(weight)
    locations=PchipInterpolator(weight,s)(np.linspace(weight[0],weight[-1],count))
    xy,tangent=mf.splineval(spline,locations),mf.splinetan(spline,locations)
    if not (np.isfinite(xy).all() and np.isfinite(tangent).all()
            and np.all(np.diff(locations)>0)):
        raise ValueError('Curvature paneling failed to produce finite, ordered points')
    return xy,locations,tangent,True
