#!/usr/bin/env python3
"""Multi-element 2-D airfoil analysis with coupled viscous boundary layers.

The default backend includes natural transition and limited separation.
See README.md for model limitations, validation and inviscid/Viiflow options.
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import json
import math
import os
from pathlib import Path
import re
import sys
import time
import warnings

# These small dense systems generally run faster without BLAS oversubscription.
# Users can override these before launching Python.
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.linalg import lu_factor, lu_solve


def signed_area(xy):
    return 0.5 * np.sum(xy[:, 0] * np.roll(xy[:, 1], -1)
                        - np.roll(xy[:, 0], -1) * xy[:, 1])


def naca4(code="0012", nodes=201):
    """Cosine-spaced NACA 4-digit coordinates, closed TE, upper TE first."""
    code = str(code)
    if not re.fullmatch(r"\d{4}", code) or int(code[2:]) == 0:
        raise ValueError("naca must be a four-digit string with nonzero thickness")
    m, p, t = int(code[0]) / 100, int(code[1]) / 10, int(code[2:]) / 100
    if m and not 0 < p < 1:
        raise ValueError("A cambered NACA requires a nonzero camber position")
    x = (1 - np.cos(np.linspace(0, np.pi, max(21, nodes // 2 + 1)))) / 2
    yt = 5*t*(.2969*np.sqrt(x) - .126*x - .3516*x**2 + .2843*x**3 - .1036*x**4)
    yc, dy = np.zeros_like(x), np.zeros_like(x)
    if m:
        fore = x < p
        yc[fore] = m/p**2*(2*p*x[fore]-x[fore]**2)
        yc[~fore] = m/(1-p)**2*((1-2*p)+2*p*x[~fore]-x[~fore]**2)
        dy[fore] = 2*m/p**2*(p-x[fore])
        dy[~fore] = 2*m/(1-p)**2*(p-x[~fore])
    th = np.arctan(dy)
    upper = np.column_stack((x-yt*np.sin(th), yc+yt*np.cos(th)))
    lower = np.column_stack((x+yt*np.sin(th), yc-yt*np.cos(th)))
    upper[-1] = lower[-1] = [1., 0.]
    return np.concatenate((upper[::-1], lower[1:]))


def read_dat(path):
    """Read a Selig contour; reject Lednicer count headers explicitly."""
    rows, started = [], False
    for no, line in enumerate(Path(path).read_text().splitlines(), 1):
        line = line.split("#")[0].split("!")[0].strip()
        if not line:
            continue
        fields = line.replace(",", " ").split()
        try:
            if len(fields) != 2:
                raise ValueError
            row = [float(v.replace("D", "E").replace("d", "e")) for v in fields]
        except ValueError:
            if started:
                raise ValueError(f"{path}:{no}: expected two coordinates") from None
            continue
        started = True
        rows.append(row)
    if (rows and all(v >= 3 and v.is_integer() for v in rows[0])
            and int(sum(rows[0])) == len(rows)-1):
        raise ValueError(f"{path}: Lednicer count header: convert to Selig contour order")
    return clean_contour(np.asarray(rows, dtype=float))


def clean_contour(xy):
    xy = np.asarray(xy, dtype=float)
    if xy.ndim != 2 or xy.shape[1] != 2 or len(xy) < 9 or not np.isfinite(xy).all():
        raise ValueError("An element needs at least 9 finite (x,y) coordinates")
    scale = np.linalg.norm(np.ptp(xy, axis=0))
    if scale <= 0:
        raise ValueError("Zero-size element")
    xy = xy[np.r_[True, np.linalg.norm(np.diff(xy, axis=0), axis=1) > scale*1e-12]]
    if signed_area(xy) < 0:
        xy = xy[::-1].copy()
    if abs(signed_area(xy)) < scale**2*1e-8:
        raise ValueError("Degenerate element contour")
    # Closing a blunt TE by repeating the first node destroys the Kutta indices.
    if np.linalg.norm(xy[0]-xy[-1]) < scale*1e-10:
        first = xy[1]-xy[0]
        last = xy[-1]-xy[-2]
        if np.dot(first, last)/(np.linalg.norm(first)*np.linalg.norm(last)) > -.3:
            raise ValueError("Remove the closing base segment from this blunt-TE .dat file")
    if np.linalg.norm(xy[0]-xy[-1]) > .05*scale:
        raise ValueError("Contour must start and end on the two trailing-edge sides")
    return xy


def repanel(xy, panels):
    """Interpolate upper/lower branches in arclength with cosine spacing.

    The LE is the node farthest from the TE midpoint. This assumes conventional
    airfoils. Splines preserve the input shape; coarse source data stay coarse.
    """
    xy = clean_contour(xy)
    if not isinstance(panels, int) or panels < 40 or panels > 2000:
        raise ValueError("panels must be an integer between 40 and 2000 per element")
    le = int(np.argmax(np.linalg.norm(xy - .5*(xy[0]+xy[-1]), axis=1)))
    if le < 3 or le > len(xy)-4:
        raise ValueError("Cannot locate a leading edge between the trailing-edge ends")
    def branch(pts, n):
        s = np.r_[0, np.cumsum(np.linalg.norm(np.diff(pts, axis=0), axis=1))]
        q = .5*(1-np.cos(np.linspace(0, np.pi, n+1))) * s[-1]
        return PchipInterpolator(s, pts, axis=0)(q)
    return np.concatenate((branch(xy[:le+1], panels//2),
                           branch(xy[le:], panels-panels//2)[1:]))


def transform(xy, scale=1., angle_deg=0., pivot=(0., 0.), translation=(0., 0.)):
    """Scale, rotate CCW about scaled pivot, then translate; row coordinates."""
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError("scale must be positive")
    a = np.deg2rad(angle_deg)
    r = np.array([[np.cos(a), -np.sin(a)], [np.sin(a), np.cos(a)]])
    pivot = np.asarray(pivot)*scale
    return (xy*scale-pivot) @ r.T + pivot + np.asarray(translation)


def _cross(a, b):
    return a[..., 0]*b[..., 1]-a[..., 1]*b[..., 0]


def _closed_segments(xy):
    if np.linalg.norm(xy[0]-xy[-1]) < 1e-12:
        return xy[:-1], xy[1:]
    return xy, np.roll(xy, -1, axis=0)


def _intersections(a, b, c, d):
    """All-pairs segment touch/cross test, including collinear overlaps."""
    ab, cd = b-a, d-c
    o1 = _cross(ab[:, None], c[None]-a[:, None])
    o2 = _cross(ab[:, None], d[None]-a[:, None])
    o3 = _cross(cd[None], a[:, None]-c[None])
    o4 = _cross(cd[None], b[:, None]-c[None])
    eps = 1e-13
    bbox = np.all(np.maximum(np.minimum(a,b)[:,None], np.minimum(c,d)[None])
                  <= np.minimum(np.maximum(a,b)[:,None], np.maximum(c,d)[None])+eps, axis=2)
    return bbox & (o1*o2 <= eps**2) & (o3*o4 <= eps**2)


def validate_geometry(elements):
    """Reject self-intersections, intersecting elements, and nested bodies."""
    from matplotlib.path import Path as MPath
    for i, (name, xy) in enumerate(elements):
        a,b = _closed_segments(xy)
        hit = _intersections(a,b,a,b)
        n = len(a)
        for k in range(n):
            hit[k,k] = hit[k,(k+1)%n] = hit[k,(k-1)%n] = False
        if hit.any():
            raise ValueError(f"Self-intersection in {name}")
        for other, uv in elements[:i]:
            c,d = _closed_segments(uv)
            if (_intersections(a,b,c,d).any() or MPath(xy).contains_point(uv[0])
                    or MPath(uv).contains_point(xy[0])):
                raise ValueError(f"Elements {name} and {other} intersect or contain each other")


def load_case(path):
    path = Path(path).resolve()
    cfg = json.loads(path.read_text())
    cref = float(cfg.get("reference_chord", 1.))
    if not np.isfinite(cref) or cref <= 0:
        raise ValueError("reference_chord must be positive")
    ref = np.asarray(cfg.get("moment_reference", [.25*cref, 0.]), dtype=float)
    if ref.shape != (2,) or not np.isfinite(ref).all():
        raise ValueError("moment_reference must be two finite coordinates")
    elements = []
    for i, e in enumerate(cfg["elements"]):
        n = e.get("panels", cfg.get("panels", 180))
        if ("naca" in e) == ("file" in e):
            raise ValueError("Each element needs exactly one of naca or file")
        xy = naca4(e["naca"], 801) if "naca" in e else read_dat(path.parent/e["file"])
        xy = repanel(xy, n)
        xy = transform(xy, e.get("scale",1.), e.get("angle_deg",0.),
                       e.get("pivot",[0.,0.]), e.get("translation",[0.,0.])) / cref
        if not np.isfinite(xy).all():
            raise ValueError("Geometry transforms must be finite")
        name = e.get("name", f"element_{i}")
        if not re.fullmatch(r"[A-Za-z0-9_-]+", name) or name in [a for a,b in elements]:
            raise ValueError("Element names must be unique ASCII letters/digits/hyphens/underscores")
        elements.append((name,xy))
    if not elements:
        raise ValueError("At least one element is required")
    if sum(len(xy) for _,xy in elements) > 3000:
        raise ValueError("Total panel budget exceeded (3000 nodes); reduce resolution")
    validate_geometry(elements)
    return cfg, elements, ref/cref


@dataclass
class Solution:
    alpha: float
    cl: float
    cm: float
    cd: float | None
    cd_pressure_numerical: float | None
    converged: bool
    residual: float
    surfaces: list
    diagnostics: dict


class HessSmithSolver:
    """Hess-Smith: N constant sources + one constant vortex sheet per body.

    Enforces no penetration and one Kutta condition per element in one dense
    system. Unit freestream and coordinates divided by reference chord.
    Inviscid incompressible potential flow: no Reynolds-number dependence.
    """
    def __init__(self, elements, moment_reference=(.25,0.)):
        self.elements = elements
        self.ref = np.asarray(moment_reference)
        starts, ends, owners, self.slices, self.kutta = [], [], [], [], []
        offset = 0
        for k, (_, xy) in enumerate(elements):
            a,b = _closed_segments(xy)
            starts.append(a); ends.append(b)
            owners.extend([k]*len(a))
            self.slices.append(slice(offset, offset+len(a)))
            self.kutta.append((offset, offset+len(xy)-2))
            offset += len(a)
        self.a, self.b = np.concatenate(starts), np.concatenate(ends)
        self.owner = np.asarray(owners)
        self.length = np.linalg.norm(self.b-self.a, axis=1)
        if np.min(self.length) < 1e-12:
            raise ValueError("Zero-length panel")
        self.t = (self.b-self.a)/self.length[:,None]
        self.n = np.column_stack((self.t[:,1], -self.t[:,0]))
        self.mid = .5*(self.a+self.b)
        ns, ne = len(self.a), len(elements)
        vs = self._source_velocity(self.mid)
        # Exact exterior self limits avoid branch ambiguity at a surface.
        vs[np.arange(ns), np.arange(ns)] = .5*self.n
        vv = np.stack((-vs[:,:,1],vs[:,:,0]),axis=2)
        membership = np.eye(ne)[self.owner]
        normal_s = np.einsum("ijk,ik->ij", vs, self.n)
        tangent_s = np.einsum("ijk,ik->ij", vs, self.t)
        normal_v = np.einsum("ijk,ik->ij", vv, self.n) @ membership
        tangent_v = np.einsum("ijk,ik->ij", vv, self.t) @ membership
        matrix = np.zeros((ns+ne, ns+ne))
        matrix[:ns,:ns], matrix[:ns,ns:] = normal_s, normal_v
        rhs = np.zeros((ns+ne,2))
        rhs[:ns] = -self.n
        for k,(i,j) in enumerate(self.kutta):
            matrix[ns+k,:ns] = tangent_s[i]+tangent_s[j]
            matrix[ns+k,ns:] = tangent_v[i]+tangent_v[j]
            rhs[ns+k] = -(self.t[i]+self.t[j])
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            lu = lu_factor(matrix)
        self.basis = lu_solve(lu, rhs)
        self.matrix, self.rhs = matrix, rhs
        self.tangent_operator = np.hstack((tangent_s,tangent_v))
        self.ne, self.ns = ne, ns

    def _source_velocity(self, points):
        d = np.asarray(points)[:,None,:]-self.a[None,:,:]
        x = np.einsum("ijk,jk->ij",d,self.t)
        y = np.einsum("ijk,jk->ij",d,self.n)
        r1 = x*x+y*y
        r2 = (x-self.length)**2+y*y
        tiny = np.finfo(float).tiny
        u = np.log(np.maximum(r1,tiny)/np.maximum(r2,tiny))/(4*np.pi)
        angle = np.arctan2(y,x-self.length)-np.arctan2(y,x)
        angle = (angle+np.pi) % (2*np.pi)-np.pi
        v = angle/(2*np.pi)
        return u[:,:,None]*self.t[None]+v[:,:,None]*self.n[None]

    def solve(self, alpha):
        if not np.isfinite(alpha):
            raise ValueError("alpha must be finite")
        angle = np.deg2rad(alpha)
        u = np.array([np.cos(angle), np.sin(angle)])
        z = self.basis @ u
        vt = self.t @ u + self.tangent_operator @ z
        cp = 1-vt*vt
        force = -cp[:,None]*self.n*self.length[:,None]
        total = force.sum(axis=0)
        lift = np.array([-u[1],u[0]])
        # Conventional airfoil pitching moment: nose-up positive.
        cm = -np.sum(_cross(self.mid-self.ref,force))
        surfaces, source_flux, circulation_cl = [], [], []
        for k, ((name,_), ix) in enumerate(zip(self.elements,self.slices)):
            surfaces.append({"name":name,"xy":self.mid[ix].copy(),"cp":cp[ix].copy(),
                             "ue_signed":vt[ix].copy(),"cl_pressure":float(force[ix].sum(axis=0)@lift)})
            source_flux.append(float(z[:self.ns][ix]@self.length[ix]))
            circulation_cl.append(float(-2*z[self.ns+k]*self.length[ix].sum()))
        residual = float(np.max(np.abs(self.matrix@z-self.rhs@u)))
        return Solution(float(alpha),float(total@lift),float(cm),None,float(total@u),
                        bool(np.isfinite(cp).all() and residual < 1e-8),residual,surfaces,
                        {"source_flux_by_element":source_flux,
                         "cl_circulation_total":sum(circulation_cl),
                         "kutta_residual":float(max(abs(vt[i]+vt[j]) for i,j in self.kutta)),
                         "model":"inviscid Hess-Smith (constant source/vortex)",
                         "note":"Pressure drag is discretization error, not viscous drag."})

    def velocity(self, points, alpha, chunk_size=256):
        angle = np.deg2rad(alpha)
        u = np.array([np.cos(angle),np.sin(angle)])
        z = self.basis @ u
        gamma = z[self.ns:][self.owner]
        out = np.empty_like(np.asarray(points,dtype=float))
        for start in range(0,len(out),chunk_size):
            ix = slice(start,start+chunk_size)
            vs = self._source_velocity(np.asarray(points)[ix])
            vv = np.stack((-vs[:,:,1],vs[:,:,0]),axis=2)
            out[ix] = u + np.einsum("ijk,j->ik",vs,z[:self.ns])+np.einsum("ijk,j->ik",vv,gamma)
        return out


class LinearVortexSolver(HessSmithSolver):
    """Linear vortex sheets for sharp trailing edges; all bodies coupled.

    Unknowns are vortex strengths at each surface node. There is a separate
    pair of coincident upper/lower TE nodes on every body. The extra equation
    gamma_upper_TE + gamma_lower_TE = 0 selects finite trailing-edge flow.
    Total lift is obtained from far-field circulation (Kutta-Joukowski).
    """
    def __init__(self,elements,moment_reference=(.25,0.)):
        self.elements, self.ref = elements,np.asarray(moment_reference)
        self.slices, self.node_slices = [],[]
        starts,ends,first,last = [],[],[],[]
        ip,inode = 0,0
        for _,xy in elements:
            if np.linalg.norm(xy[0]-xy[-1]) > 1e-10:
                raise ValueError("Linear vortex backend requires sharp trailing edges")
            n = len(xy)-1
            starts.append(xy[:-1]); ends.append(xy[1:])
            self.slices.append(slice(ip,ip+n))
            self.node_slices.append(slice(inode,inode+n+1))
            first.extend(range(inode,inode+n)); last.extend(range(inode+1,inode+n+1))
            ip += n; inode += n+1
        self.a,self.b = np.concatenate(starts),np.concatenate(ends)
        self.first,self.last = np.asarray(first),np.asarray(last)
        self.length = np.linalg.norm(self.b-self.a,axis=1)
        self.t = (self.b-self.a)/self.length[:,None]
        self.n = np.column_stack((self.t[:,1],-self.t[:,0]))
        self.mid = .5*(self.a+self.b)
        self.ns,self.nn = ip,inode
        left,right = self._linear_velocity(self.mid,on_surface=True)
        velocity = np.zeros((ip,inode,2))
        # A node appears in adjacent panels; accumulate both contributions.
        for ps,ns in zip(self.slices,self.node_slices):
            velocity[:,ns.start:ns.stop-1] += left[:,ps]
            velocity[:,ns.start+1:ns.stop] += right[:,ps]
        matrix = np.zeros((inode,inode))
        matrix[:ip] = np.einsum("ijk,ik->ij",velocity,self.n)
        rhs = np.zeros((inode,2)); rhs[:ip] = -self.n
        for k,ns in enumerate(self.node_slices):
            matrix[ip+k,[ns.start,ns.stop-1]] = 1
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            lu = lu_factor(matrix)
        self.basis = lu_solve(lu,rhs)
        self.matrix,self.rhs = matrix,rhs
        self.tangent_operator = np.einsum("ijk,ik->ij",velocity,self.t)

    def _linear_velocity(self,points,on_surface=False):
        d = np.asarray(points)[:,None,:]-self.a[None]
        x = np.einsum("ijk,jk->ij",d,self.t)
        y = np.einsum("ijk,jk->ij",d,self.n)
        tiny = np.finfo(float).tiny
        u = np.log(np.maximum(x*x+y*y,tiny)/np.maximum((x-self.length)**2+y*y,tiny))/(4*np.pi)
        angle = np.arctan2(y,x-self.length)-np.arctan2(y,x)
        v = ((angle+np.pi)%(2*np.pi)-np.pi)/(2*np.pi)
        if on_surface:
            np.fill_diagonal(u,0.); np.fill_diagonal(v,.5)
        # Analytical first moment of the panel kernel, weighted by s/L.
        u1 = x/self.length*u - 1/(2*np.pi) + y/self.length*v
        v1 = x/self.length*v - y/self.length*u
        uniform = v[:,:,None]*self.t[None]-u[:,:,None]*self.n[None]
        right = v1[:,:,None]*self.t[None]-u1[:,:,None]*self.n[None]
        return uniform-right,right

    def solve(self,alpha):
        if not np.isfinite(alpha):
            raise ValueError("alpha must be finite")
        angle = np.deg2rad(alpha)
        u = np.array([np.cos(angle),np.sin(angle)])
        gamma = self.basis @ u
        vt = self.t@u+self.tangent_operator@gamma
        cp = 1-vt**2
        force = -cp[:,None]*self.n*self.length[:,None]
        lift = np.array([-u[1],u[0]])
        circ = .5*(gamma[self.first]+gamma[self.last])*self.length
        cl = -2*np.sum(circ)
        cm = -np.sum(_cross(self.mid-self.ref,force))
        surfaces = []
        for (name,_),ix in zip(self.elements,self.slices):
            surfaces.append({"name":name,"xy":self.mid[ix].copy(),"cp":cp[ix].copy(),
                             "ue_signed":vt[ix].copy(),"cl_pressure":float(force[ix].sum(axis=0)@lift)})
        residual = float(np.max(abs(self.matrix@gamma-self.rhs@u)))
        return Solution(float(alpha),float(cl),float(cm),None,float(force.sum(axis=0)@u),
                        bool(np.isfinite(cp).all() and residual<1e-8),residual,surfaces,
                        {"model":"inviscid linear vortex panels",
                         "cl_pressure_total":float(force.sum(axis=0)@lift),
                         "cl_circulation_total":float(cl),
                         "kutta_residual":float(max(abs(gamma[s.start]+gamma[s.stop-1]) for s in self.node_slices)),
                         "note":"CL uses total circulation; pressure drag is discretization error, not viscous drag."})

    def velocity(self,points,alpha,chunk_size=256):
        angle = np.deg2rad(alpha)
        u = np.array([np.cos(angle),np.sin(angle)])
        gamma = self.basis@u
        out = np.empty_like(np.asarray(points,dtype=float))
        for start in range(0,len(out),chunk_size):
            ix = slice(start,start+chunk_size)
            left,right = self._linear_velocity(np.asarray(points)[ix])
            out[ix] = (u + np.einsum("ijk,j->ik",left,gamma[self.first])
                       + np.einsum("ijk,j->ik",right,gamma[self.last]))
        return out


def PanelSolver(elements,moment_reference=(.25,0.)):
    """Choose linear vortex panels for sharp TEs, Hess-Smith for blunt TEs."""
    if not elements or sum(len(xy) for _,xy in elements)>3000:
        raise ValueError("Use at least one element and at most 3000 total nodes")
    elements = [(name,clean_contour(xy)) for name,xy in elements]
    validate_geometry(elements)
    sharp = all(np.linalg.norm(xy[0]-xy[-1])<1e-10 for _,xy in elements)
    solver_type = LinearVortexSolver if sharp else HessSmithSolver
    return solver_type(elements,moment_reference)


class ViiflowSolver:
    """Public API adapter. Multi-element use requires a suitable Pro license.

    No fallback is performed when a license, import, or convergence fails.
    """
    def __init__(self,elements,moment_reference=(.25,0.),reynolds=1e6,mach=0.,ncrit=9.,
                 iterations=150,iterate_wakes=False):
        try:
            import viiflow as vf
        except ImportError as exc:
            raise RuntimeError("Viiflow is not installed. See README.md; multi-element needs Pro.") from exc
        if not np.isfinite(reynolds) or reynolds <= 0:
            raise ValueError("Reynolds number must be positive")
        if not np.isfinite(mach) or not 0 <= mach < 1:
            raise ValueError("Mach number must satisfy 0 <= Ma < 1")
        if not np.isfinite(ncrit) or ncrit <= 0 or iterations < 1:
            raise ValueError("Ncrit and iteration count must be positive")
        self.vf, self.elements = vf, elements
        self.settings = vf.setup(Re=reynolds, Ma=mach, Ncrit=ncrit, Silent=True,
                                 Itermax=iterations, IterateWakes=iterate_wakes,
                                 CmRefX=float(moment_reference[0]), CmRefY=float(moment_reference[1]))
        self.state = None

    def _init(self):
        try:
            self.state = self.vf.init([np.ascontiguousarray(xy.T) for _,xy in self.elements],self.settings)
        except Exception as exc:
            self.state = None
            raise RuntimeError(f"Viiflow initialization failed: {exc}. "
                               "The Basic 1.9.4 binary only permits one element; multi-element needs Pro.") from exc

    def solve(self,alpha):
        if not np.isfinite(alpha):
            raise ValueError("alpha must be finite")
        self.settings.Alpha = float(alpha)
        was_warm = self.state is not None
        if self.state is None:
            self._init()
        p,bl,x = self.state
        out = self.vf.iter(x,bl,p,self.settings)
        # A failed continuation gets one fresh initialization at this alpha.
        if not out[1] and was_warm:
            self._init()
            p,bl,x = self.state
            out = self.vf.iter(x,bl,p,self.settings)
        x, flag, res = out[:3]
        cd = float(sum(float(b.CD) for b in bl))
        cl,cm = float(p.CL),float(p.CM)
        flag = bool(flag and np.isfinite([cl,cm,cd]).all())
        self.state = (p,bl,x) if flag else None
        surfaces, offset = [],0
        if flag:
            cp = np.asarray(p.cp).reshape(-1)
            for k,(name,_) in enumerate(self.elements):
                xy = np.asarray(p.foils[k].X).T.copy()
                n = len(xy)
                data = {"name":name,"xy":xy,"cp":cp[offset:offset+n].copy(),"cd_wake":float(bl[k].CD)}
                nodes = bl[k].bl_fl.nodes
                for key, attr in [("cf","CF"),("theta","theta"),("delta_star","delta"),("H","H")]:
                    values = np.asarray(nodes[attr]).reshape(-1)
                    if len(values) == n:
                        data[key] = values.copy()
                surfaces.append(data)
                offset += n
        return Solution(float(alpha),cl if flag else math.nan,cm if flag else math.nan,
                        cd if flag else math.nan,None,flag,float(np.linalg.norm(res)),surfaces,
                        {"model":"Viiflow viscous-inviscid", "cd_definition":"sum of element wake drag coefficients"})


def _save_csv(path, rows, fields):
    with Path(path).open("w",newline="") as f:
        w = csv.DictWriter(f,fieldnames=fields)
        w.writeheader(); w.writerows(rows)


def save_results(folder, solutions, elements, cfg, backend):
    folder = Path(folder)
    folder.mkdir(parents=True,exist_ok=True)
    fields = ["alpha_deg","converged","CL","CD","CM","CD_pressure_numerical","residual"]
    rows = []
    for i,r in enumerate(solutions):
        rows.append(dict(zip(fields,[r.alpha,int(r.converged),r.cl,r.cd,r.cm,r.cd_pressure_numerical,r.residual])))
        for s in r.surfaces:
            arrays = {"x_over_cref":s["xy"][:,0],"y_over_cref":s["xy"][:,1],"Cp":s["cp"]}
            arrays.update({k:s[k] for k in ["x_over_local_chord","side","ue_signed","ue_over_vinf", "cf", "theta",
                "delta_star","H","amplification_or_shear","turbulent"] if k in s})
            if backend=='viscous' and 'H' in arrays:arrays['Hk']=arrays.pop('H')
            _save_csv(folder/f"surface_{i:03d}_{s['name']}.csv",
                      [dict(zip(arrays,vals)) for vals in zip(*arrays.values())],list(arrays))
            if "wake" in s:
                arrays=s['wake']
                _save_csv(folder/f"wake_{i:03d}_{s['name']}.csv",
                          [dict(zip(arrays,vals)) for vals in zip(*arrays.values())],list(arrays))
    _save_csv(folder/"polar.csv",rows,fields)
    for name,xy in elements:
        np.savetxt(folder/f"geometry_{name}.dat",xy,header=name,comments="",fmt="%.12g")
    import importlib.metadata as metadata
    versions = {"python":sys.version.split()[0]}
    for name in ["numpy","scipy","matplotlib"]+(["viiflow"] if backend=="viiflow" else []):
        try:
            versions[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            versions[name] = "unknown"
    meta = {"backend":backend,"versions":versions,"case":cfg,"coordinate_units":"x/reference_chord, y/reference_chord",
            "diagnostics":[r.diagnostics for r in solutions]}
    # Use strict JSON: failed solves have NaN in CSV/Python, null in JSON.
    def json_safe(x):
        if isinstance(x,(float,np.floating)) and not np.isfinite(x):return None
        if isinstance(x,dict):return {k:json_safe(v) for k,v in x.items()}
        if isinstance(x,(list,tuple)):return [json_safe(v) for v in x]
        return x
    if backend=='panel':meta['note']='Inviscid: no viscous CD, Reynolds effects, transition or stall.'
    (folder/"metadata.json").write_text(json.dumps(json_safe(meta),indent=2,allow_nan=False))
    if backend=='viscous':
        from viscous import analyze_stall
        (folder/'stall.json').write_text(json.dumps(json_safe(analyze_stall(solutions)),indent=2,allow_nan=False))
        fields=['alpha_deg','element','side','transition_x_over_local_chord',
                'separated_arclength_fraction','min_cf','max_Hk']
        rows=[]
        for r in solutions:
            for name,bl in r.diagnostics.get('boundary_layers',{}).items():
                for side,data in bl.items():
                    tr=data['transition']
                    rows.append(dict(zip(fields,[r.alpha,name,side,None if tr is None else tr['x_over_local_chord'],
                        data['separated_arclength_fraction'],data['min_cf'],data['max_Hk']])))
        _save_csv(folder/'transition_separation.csv',rows,fields)


def plot_results(folder, solutions, elements, solver=None, flow=False):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.path import Path as MPath
    plt.rcParams.update({"font.size":10,"axes.spines.top":False,"axes.spines.right":False})
    folder = Path(folder)
    ok = [r for r in solutions if r.converged]
    if not ok:
        return
    picked = ok[len(ok)//2]
    fig,(ax,geo) = plt.subplots(2,1,figsize=(10,7),gridspec_kw={"height_ratios":[2,1]},layout="constrained")
    for s,(name,xy) in zip(picked.surfaces,elements):
        ax.plot(s["xy"][:,0],s["cp"],label=name)
        geo.plot(xy[:,0],xy[:,1],label=name)
    ax.invert_yaxis(); ax.grid(alpha=.2); ax.legend(); ax.set_ylabel("Pressure coefficient $C_p$")
    ax.set_title(f"{picked.diagnostics['model']} | alpha={picked.alpha:g} deg | CL={picked.cl:.4f}")
    geo.set_aspect("equal"); geo.set_xlabel("x / reference chord"); geo.set_ylabel("y / reference chord")
    fig.savefig(folder/"pressure.png",dpi=160,bbox_inches='tight'); plt.close(fig)
    fig,axes = plt.subplots(1,2,figsize=(10,4),layout="constrained")
    axes[0].plot([r.alpha for r in solutions],[r.cl if r.converged else np.nan for r in solutions],"o-")
    axes[0].set_xlabel("Angle of attack [deg]"); axes[0].set_ylabel("$C_L$")
    failed=[r.alpha for r in solutions if not r.converged]
    if failed:
        axes[0].plot(failed,[.025]*len(failed),'x',color='#777777',
                     transform=axes[0].get_xaxis_transform(),label='Not converged')
        axes[0].legend(loc='lower left')
    if picked.cd is not None:
        axes[1].plot([r.cd if r.converged else np.nan for r in solutions],
                     [r.cl if r.converged else np.nan for r in solutions],"o-")
        axes[1].set_xlabel("$C_D$ (viscous)"); axes[1].set_ylabel("$C_L$")
    else:
        axes[1].plot([r.alpha for r in solutions],[r.cm for r in solutions],"o-")
        axes[1].set_xlabel("Angle of attack [deg]"); axes[1].set_ylabel("$C_M$ (nose-up positive)")
    for ax in axes: ax.grid(alpha=.2)
    fig.suptitle(picked.diagnostics["model"],fontsize=11)
    fig.savefig(folder/"polar.png",dpi=160,bbox_inches='tight'); plt.close(fig)
    if any('side' in s for s in picked.surfaces):
        separated=max(ok,key=lambda r:r.alpha)
        fig,axes=plt.subplots(2,1,figsize=(10,7),sharex=True,layout='constrained')
        for s in separated.surfaces:
            for side in ['upper','lower']:
                ids=np.flatnonzero(s['side']==side)
                ids=ids[np.argsort(s['x_over_local_chord'][ids])]
                label=s['name']+' '+side
                axes[0].plot(s['x_over_local_chord'][ids],s['cf'][ids],label=label)
                axes[1].plot(s['x_over_local_chord'][ids],s['H'][ids],label=label)
        axes[0].axhline(0,color='black',lw=.7);axes[0].legend(ncol=2)
        axes[0].set_ylabel('$C_f$ (freestream based)');axes[1].set_ylabel('Kinematic shape factor $H_k$')
        axes[1].set_xlabel('x / local element chord')
        for ax in axes:ax.grid(alpha=.2)
        fig.suptitle(f'Boundary layers | alpha={separated.alpha:g} deg | negative Cf indicates separation')
        fig.savefig(folder/'boundary_layers.png',dpi=160,bbox_inches='tight');plt.close(fig)
    if flow and isinstance(solver,HessSmithSolver):
        all_xy = np.concatenate([xy for _,xy in elements])
        lo,hi = all_xy.min(axis=0),all_xy.max(axis=0)
        extent = hi[0]-lo[0]
        xx = np.linspace(lo[0]-.3*extent,hi[0]+.45*extent,220)
        yy = np.linspace(lo[1]-.5*extent,hi[1]+.5*extent,140)
        x,y = np.meshgrid(xx,yy)
        pts = np.column_stack((x.ravel(),y.ravel()))
        mask = np.zeros(len(pts),dtype=bool)
        for _,xy in elements: mask |= MPath(xy).contains_points(pts)
        vel = solver.velocity(pts,picked.alpha)
        vel[mask] = np.nan
        u,v = vel[:,0].reshape(x.shape),vel[:,1].reshape(x.shape)
        fig,ax = plt.subplots(figsize=(11,5),layout="constrained")
        speed = np.sqrt(u*u+v*v)
        im = ax.pcolormesh(x,y,speed,cmap="viridis",vmin=0,vmax=3,shading="auto")
        ax.streamplot(xx,yy,u,v,color="white",density=1.5,linewidth=.55,arrowsize=.65)
        for _,xy in elements: ax.fill(xy[:,0],xy[:,1],color="#18222c",zorder=5)
        ax.set_aspect("equal"); ax.set_xlabel("x / reference chord"); ax.set_ylabel("y / reference chord")
        ax.set_title(f"Inviscid flow | alpha={picked.alpha:g} deg")
        fig.colorbar(im,ax=ax,label="Speed / freestream speed")
        fig.savefig(folder/"flow.png",dpi=160); plt.close(fig)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case",type=Path,help="Case JSON; .dat paths are relative to this file")
    parser.add_argument("--backend",choices=["viscous","panel","viiflow"],default="viscous")
    parser.add_argument("--alpha",type=float,nargs="+",help="Angles in degrees, e.g. --alpha -2 0 4 8")
    parser.add_argument("--re",type=float,dest="reynolds",help="Reynolds number based on reference_chord")
    parser.add_argument("--mach",type=float)
    parser.add_argument("--ncrit",type=float)
    parser.add_argument("--amplification-jacobian",choices=["corrected_af","legacy"],
                        help="Viscous amplification derivative (default corrected_af)")
    parser.add_argument("--iterations",type=int,help="Maximum Newton iterations at each continuation point")
    parser.add_argument("--continuation-step",type=float,help="Maximum internal alpha step, degrees (default 1)")
    parser.add_argument("--start-alpha",type=float,help="Initial angle for viscous continuation")
    parser.add_argument("--wake-length",type=float,help="Wake extension past the assembly, in reference chords")
    parser.add_argument("--verbose",action="store_true",help="Print Newton convergence history")
    parser.add_argument("--out",type=Path,default=Path("results"))
    parser.add_argument("--flow",action="store_true",help="Also draw inviscid streamlines (panel backend)")
    args = parser.parse_args(argv)
    try:
        cfg,elements,ref = load_case(args.case)
        if args.amplification_jacobian is not None and args.backend != "viscous":
            raise ValueError("--amplification-jacobian requires --backend viscous")
        angles = args.alpha if args.alpha is not None else cfg.get("alpha_deg",[-2,0,2,4,6,8])
        if not isinstance(angles,list) or not angles or not np.isfinite(angles).all():
            raise ValueError("alpha_deg must be a nonempty list of finite angles")
        cfg["alpha_deg"] = angles
        for key,value in [("reynolds",args.reynolds),("mach",args.mach),("ncrit",args.ncrit),
                          ("amplification_jacobian",args.amplification_jacobian),
                          ("iterations",args.iterations),("continuation_step",args.continuation_step),
                          ("wake_length",args.wake_length)]:
            if value is not None:
                cfg[key] = value
        if args.start_alpha is not None:cfg['start_alpha']=args.start_alpha
        if args.backend == "panel":
            if args.mach not in (None,0.) or args.reynolds is not None or args.ncrit is not None:
                raise ValueError("--re/--ncrit and nonzero --mach require a viscous backend")
            print("INVISCID: CD, transition and stall are not predicted; case Re/Mach/Ncrit are unused.")
            solver = PanelSolver(elements,ref)
        elif args.backend=='viscous':
            if args.flow:raise ValueError("--flow is available only with --backend panel")
            from viscous import CoupledSolver
            for key,value in [('reynolds',1e6),('mach',0.),('ncrit',9.),('iterations',80),
                ('amplification_jacobian','corrected_af'),
                ('tolerance',1e-7),('wake_length',1.),('minimum_te_gap',1e-5),
                ('continuation_step',1.),('min_step',.25),('start_alpha',None)]:cfg.setdefault(key,value)
            solver=CoupledSolver(elements,ref,**{key:cfg[key] for key in ['reynolds','mach','ncrit','iterations',
                'amplification_jacobian',
                'tolerance','wake_length','minimum_te_gap','continuation_step','min_step','start_alpha']},verbose=args.verbose)
            elements=solver.elements
            print('VISCOUS: coupled mfoil boundary layers, natural transition and limited separation; independent wakes.')
        else:
            if args.flow:
                raise ValueError("--flow is currently available only with --backend panel")
            def setting(arg,key,default):
                return arg if arg is not None else cfg.get(key,default)
            solver = ViiflowSolver(elements,ref,setting(args.reynolds,"reynolds",1e6),
                                   setting(args.mach,"mach",0.),setting(args.ncrit,"ncrit",9.),
                                   cfg.get("iterations",150),cfg.get("iterate_wakes",False))
        if args.out.exists() and any(args.out.iterdir()):
            raise ValueError(f"Output directory {args.out} is not empty; choose a new --out")
        start = time.perf_counter()
        solutions = []
        for alpha in angles:
            r = solver.solve(alpha)
            solutions.append(r)
            cd_text = "unavailable (inviscid)" if r.cd is None else f"{r.cd:.6f}"
            print(f"alpha={r.alpha:7.3f} CL={r.cl:9.5f} CD={cd_text} CM={r.cm:9.5f} converged={r.converged}",flush=True)
            if not r.converged:print('  '+r.diagnostics.get('failure_reason','Solver did not converge'),flush=True)
            if any(r.diagnostics.get('local_supersonic_detected',{}).values()):
                print('  Local supersonic flow detected: compressibility correction is outside its reliable range.',flush=True)
        save_results(args.out,solutions,elements,cfg,args.backend)
        plot_results(args.out,solutions,elements,solver,args.flow)
        if args.backend=='viscous':
            from viscous import analyze_stall
            stall=analyze_stall(solutions)
            print('Stall diagnostic: '+stall['status']+' (see stall.json; failure is not stall)')
        print(f"Saved {args.out.resolve()} in {time.perf_counter()-start:.2f} s")
        return 0 if all(r.converged for r in solutions) else 2
    except (ValueError,RuntimeError,OSError,np.linalg.LinAlgError) as exc:
        print(f"ERROR: {exc}",file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
