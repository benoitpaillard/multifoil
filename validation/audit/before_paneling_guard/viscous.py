"""Multi-element, Newton-coupled viscous/inviscid solver using mfoil closures.

Boundary-layer physics: Krzysztof J. Fidkowski's MIT-licensed mfoil (vendor/).
This module adds a global multibody potential flow and viscous source coupling.
It predicts transition and limited separation, not fully separated unsteady flow.
"""
from __future__ import annotations
import copy
import math
import os
import time
from types import SimpleNamespace
import warnings

os.environ.setdefault("OPENBLAS_NUM_THREADS","1")
os.environ.setdefault("OMP_NUM_THREADS","1")
import numpy as np
from scipy import linalg

from vendor import mfoil as mf
from multifoil import Solution, clean_contour, validate_geometry, _closed_segments, _intersections


def te_info(X):
    """mfoil's TE geometry, with a well-defined zero-gap limit."""
    t1,t2 = X[:,0]-X[:,1],X[:,-1]-X[:,-2]
    t1,t2 = t1/np.linalg.norm(t1),t2/np.linalg.norm(t2)
    t=(t1+t2)/np.linalg.norm(t1+t2)
    s=X[:,-1]-X[:,0]
    gap=np.linalg.norm(s)
    hte=-s[0]*t[1]+s[1]*t[0]
    dtdx=t1[0]*t2[1]-t2[0]*t1[1]
    p=s/gap if gap>1e-12 else np.array([-t[1],t[0]])
    return t,hte,dtdx,abs(t[0]*p[1]-t[1]*p[0]),float(t@p)


# The original TE_info divides by zero for a sharp TE. This numerical limit
# leaves all finite-gap values and all boundary-layer closures unchanged.
mf.TE_info=te_info


def new_model(name,xy,reynolds,mach,ncrit,ref,minimum_te_gap):
    M=SimpleNamespace()
    for key,typ in [("geom",mf.Geom),("foil",mf.Panel),("wake",mf.Panel),
                    ("oper",mf.Oper),("isol",mf.Isol),("vsol",mf.Vsol),
                    ("glob",mf.Glob),("post",mf.Post),("param",mf.Param)]:
        setattr(M,key,typ())
    # mfoil uses clockwise ordering: lower TE -> LE -> upper TE.
    X=np.ascontiguousarray(clean_contour(xy)[::-1].T)
    te=.5*(X[:,0]+X[:,-1]);ile=int(np.argmax(np.linalg.norm(X-te[:,None],axis=0)))
    le=X[:,ile].copy();chord=np.linalg.norm(te-le);axis=(te-le)/chord
    M.regularization=0.
    current=te_info(X)[1]
    required=minimum_te_gap*chord
    if current<required:
        # Smooth, tiny opening avoids the singular sharp-TE viscous source limit.
        normal=np.array([-axis[1],axis[0]])
        weight=np.clip(axis@(X-le[:,None])/chord,0,1)**12
        signs=np.where(np.arange(X.shape[1])<=ile,-1.,1.)
        increment=(required-current)/max(abs(np.dot(normal,np.array([-te_info(X)[0][1],te_info(X)[0][0]]))),.1)
        X+=normal[:,None]*(.5*increment*weight*signs)[None,:]
        M.regularization=float(increment/chord)
    M.chord_origin=le;M.chord_axis=axis
    M.foil.N=X.shape[1]
    # Curvature spacing resolves the leading-edge bubble and avoids vanishingly
    # short cosine panels at a nearly closed TE. Perform it in local coordinates
    # so rotating an element does not change its node distribution.
    frame=np.column_stack((axis,[-axis[1],axis[0]]))
    local=frame.T@(X-le[:,None])/chord
    Y,S,T=mf.spline_curvature(local,M.foil.N,2.,.1,None)
    M.foil.x=le[:,None]+chord*(frame@Y)
    M.foil.s=chord*S;M.foil.t=frame@T
    M.geom.name=name; M.geom.xpoint=X; M.geom.npoint=len(xy)
    M.geom.chord=1. # Coordinates and Reynolds number use the assembly reference.
    M.geom.xref=np.asarray(ref,dtype=float)
    M.local_chord=float(np.max(np.linalg.norm(X-.5*(X[:,0]+X[:,-1])[:,None],axis=0)))
    M.oper.Re=reynolds; M.oper.Ma=mach; M.oper.viscous=True
    M.param.ncrit=ncrit; M.param.verb=0; M.param.doplot=False
    mf.init_thermo(M)
    return M


def vortex_stream_matrix(X,targets,cross_body=False):
    n=X.shape[1]
    out=np.zeros((targets.shape[1],n))
    _,h,_,tcp,tdp=te_info(X)
    if abs(h)>1e-12:
        base=X[:,[-1,0]]
        base_psi=np.array([mf.panel_constsource_stream(base,pt) for pt in targets.T])
        if cross_body:base_psi=np.unwrap(base_psi,period=np.linalg.norm(base[:,1]-base[:,0]))
    for i,pt in enumerate(targets.T):
        for j in range(n-1):
            a,b=mf.panel_linvortex_stream(X[:,j:j+2],pt)
            out[i,j]+=a; out[i,j+1]+=b
        if abs(h)>1e-12:
            panel=X[:,[-1,0]]
            src=base_psi[i]
            a,b=mf.panel_linvortex_stream(panel,pt)
            val=-.5*tcp*src+.5*tdp*(a+b)
            out[i,0]+=val; out[i,-1]-=val
    return out


def vortex_velocity_matrix(X,targets):
    n=X.shape[1]
    out=np.zeros((targets.shape[1],2,n))
    _,h,_,tcp,tdp=te_info(X)
    for i,pt in enumerate(targets.T):
        for j in range(n-1):
            a,b=mf.panel_linvortex_velocity(X[:,j:j+2],pt,None,False)
            out[i,:,j]+=a; out[i,:,j+1]+=b
        if abs(h)>1e-12:
            panel=X[:,[-1,0]]
            src=mf.panel_constsource_velocity(panel,pt,None)
            a,b=mf.panel_linvortex_velocity(panel,pt,None,False)
            val=-.5*tcp*src+.5*tdp*(a+b)
            out[i,:,0]+=val; out[i,:,-1]-=val
    return out


def source_stream_matrix(M,targets,cross_body=False):
    """Surface constant sources and continuous piecewise-linear wake sources.

    Discretization follows mfoil.calc_ue_m, extended to arbitrary target bodies.
    """
    N,Nw=M.foil.N,M.wake.N
    out=np.zeros((targets.shape[1],N+Nw-2))
    # Native source kernels choose branch cuts suited to their own airfoil.
    # Those cuts can cross another body. Unwrap EACH elementary kernel around
    # the target contour before assembly; a constant shift is absorbed by that
    # body's independent streamfunction constant. A source outside a closed
    # non-intersecting target has a continuous streamfunction on that target.
    def lin(panel):
        values=np.array([mf.panel_linsource_stream(panel,pt) for pt in targets.T])
        if cross_body:values=np.unwrap(values,period=.5*np.linalg.norm(panel[:,1]-panel[:,0]),axis=0)
        return values.T
    for j in range(N-1):
        panel=M.foil.x[:,j:j+2]
        values=np.array([mf.panel_constsource_stream(panel,pt) for pt in targets.T])
        if cross_body:values=np.unwrap(values,period=np.linalg.norm(panel[:,1]-panel[:,0]))
        out[:,j]=values
    for j in range(Nw-1):
        ends=M.wake.x[:,j:j+2]
        X=np.column_stack((ends[:,0],ends.mean(axis=1),ends[:,1]))
        if j==Nw-2:X[:,2]=2*X[:,2]-X[:,1]
        a,b=lin(X[:,:2]);out[:,N-1+j]+=b
        if j>0:
            out[:,N-1+j]+=.5*a;out[:,N-2+j]+=.5*a
        elif cross_body:
            # Continuity of the source at the TE: the first wake-node source
            # equals the sum of the upper/lower terminal surface sources.
            # For the body's own streamfunction retain mfoil's near-TE
            # approximation (this very short subpanel contribution is omitted).
            out[:,0]+=a;out[:,N-2]+=a
        a,b=lin(X[:,1:]);out[:,N-1+j]+=a+.5*b
        out[:,N+j if j<Nw-2 else N-1+j]+=.5*b
    return out


def source_velocity_matrix(M,targets,tangents,self_wake=False):
    """Wake-tangential velocity response; native analytical self limits."""
    N,Nw=M.foil.N,M.wake.N
    out=np.zeros((targets.shape[1],N+Nw-2))
    for i,(pt,ti) in enumerate(zip(targets.T,tangents.T)):
        at_te=self_wake and i==0
        for j in range(int(at_te),N-1-int(at_te)):
            out[i,j]=mf.panel_constsource_velocity(M.foil.x[:,j:j+2],pt,ti)
        for j in range(Nw):
            X=M.wake.x[:,[max(j-1,0),j,min(j+1,Nw-1)]].copy()
            X[:,0]=.5*(X[:,0]+X[:,1]); X[:,2]=.5*(X[:,1]+X[:,2])
            if j==Nw-1: X[:,2]=2*X[:,1]-X[:,0]
            d1,d2=np.linalg.norm(X[:,1]-X[:,0]),np.linalg.norm(X[:,2]-X[:,1])
            if self_wake and i==j:
                if j==0:
                    dl=np.linalg.norm(M.foil.x[:,1]-M.foil.x[:,0])
                    du=np.linalg.norm(M.foil.x[:,-1]-M.foil.x[:,-2])
                    out[i,0]+=(np.log(dl/d2)+1)/(2*np.pi)
                    out[i,N-2]+=(np.log(du/d2)+1)/(2*np.pi)
                    out[i,N-1]-=1/(2*np.pi)
                elif j<Nw-1:
                    a=np.log(d1/d2)/(4*np.pi)
                    out[i,N-2+j]+=a+1/(2*np.pi)
                    out[i,N-1+j]+=a-1/(2*np.pi)
            elif j==0:
                a,b=mf.panel_linsource_velocity(X[:,1:],pt,ti)
                out[i,N-1]+=b; out[i,0]+=a; out[i,N-2]+=a
            elif j==Nw-1:
                out[i,N+Nw-3]+=mf.panel_constsource_velocity(X[:,[0,2]],pt,ti)
            else:
                a1,b1=mf.panel_linsource_velocity(X[:,:2],pt,ti)
                a2,b2=mf.panel_linsource_velocity(X[:,1:],pt,ti)
                # Shared node strength is (sigma_left+sigma_right)/2 on
                # BOTH half panels. The 2023 Python source omits two terms.
                out[i,N-2+j]+=a1+.5*(b1+a2)
                out[i,N-1+j]+=b2+.5*(b1+a2)
    return out


def dense_boundary_system(M):
    """mfoil's BL residual/Jacobian with dense, modern-NumPy-compatible assembly.

    Same station, transition, wake and stagnation equations as the MIT source;
    avoids the original sparse-shape comparison bug and sparse slice overhead.
    """
    ns=M.glob.Nsys
    R=np.zeros(3*ns); RU=np.zeros((3*ns,4*ns)); RX=np.zeros((3*ns,ns))
    for side in range(3):
        Is=np.asarray(list(M.vsol.Is[side]),dtype=int)
        xi=M.isol.xi[Is]; U=M.glob.U[:,Is]; aux=np.zeros(len(Is))
        param=mf.build_param(M,side)
        if side==2: aux[:]=M.vsol.wgap
        i0=int(side<2 and xi[0]<1e-8*xi[-1])
        if side<2:
            ip=[i0,i0+1]
            ust,ust_U,ust_x,xst=mf.stagnation_state(U[:,ip],xi[ip])
            param.turb=False; param.simi=True
            r,ru,rx=mf.residual_station(param,np.array([xst,xst]),np.column_stack((ust,ust)),aux[[i0,i0]])
            param.simi=False
            rust=ru[:,:4]+ru[:,4:]
            ru=rust@ust_U; rx=rust@ust_x; js=Is[ip]
            if i0:
                sl=slice(3*Is[0],3*Is[0]+3)
                R[sl]=U[:3,0]-ust[:3]
                RU[sl,4*Is[0]:4*Is[0]+4]+=np.eye(3,4)
                for k,j in enumerate(js):
                    RU[sl,4*j:4*j+4]-=ust_U[:3,4*k:4*k+4]
                    RX[sl,j]-=ust_x[:3,k]
        else:
            r,ru,js=mf.wake_sys(M,param); rx=np.zeros((3,0))
            param.turb=True;param.wake=True
        sl=slice(3*Is[i0],3*Is[i0]+3); R[sl]=r
        for k,j in enumerate(js):
            RU[sl,4*j:4*j+4]+=ru[:,4*k:4*k+4]
            if len(rx) and rx.shape[1]: RX[sl,j]+=rx[:,k]
        for i in range(i0+1,len(Is)):
            ip=[i-1,i]; js=Is[ip]
            tr=bool(M.vsol.turb[js[0]] ^ M.vsol.turb[js[1]])
            if tr:
                r,ru,rx=mf.residual_transition(M,param,xi[ip],U[:,ip],aux[ip])
                mf.store_transition(M,side,i)
            else:
                r,ru,rx=mf.residual_station(param,xi[ip],U[:,ip],aux[ip])
            sl=slice(3*Is[i],3*Is[i]+3);R[sl]+=r
            for k,j in enumerate(js):
                RU[sl,4*j:4*j+4]+=ru[:,4*k:4*k+4]
                RX[sl,j]+=rx[:,k]
            if tr: param.turb=True
    rst=RX@np.r_[-M.isol.sgnue,-np.ones(M.wake.N)]
    for j,deriv in zip(M.isol.Istag,M.isol.sstag_ue): RU[:,4*j+3]+=rst*deriv
    return R,RU


def boundary_layer_summary(M):
    """Natural transition and Cf < 0 intervals, in local chord coordinates."""
    result={}
    for side,label in enumerate(("lower","upper")):
        ids=np.asarray(M.vsol.Is[side],dtype=int)
        xi=M.isol.xi[ids];cf=M.post.cf[ids]
        xyz=M.foil.x[:,ids]
        def location(s):
            point=np.array([np.interp(s,xi,row) for row in xyz])
            return {"s_over_cref":float(s),"xy_over_cref":point.tolist(),
                    "x_over_local_chord":float(M.chord_axis@(point-M.chord_origin)/M.local_chord)}
        turbulent=np.asarray(M.vsol.turb[ids],dtype=bool)
        tr=None
        if turbulent.any():
            first=int(np.flatnonzero(turbulent)[0])
            st=float(M.vsol.Xt[side,0])
            if first and not xi[first-1]-1e-9<=st<=xi[first]+1e-9:
                st=float(xi[first])
            tr=location(st)
        # Insert linearly interpolated zero crossings before finding intervals.
        breaks=[float(xi[0])]
        for j in range(len(ids)-1):
            if cf[j]*cf[j+1]<0:
                breaks.append(float(xi[j]-cf[j]*(xi[j+1]-xi[j])/(cf[j+1]-cf[j])))
        breaks.append(float(xi[-1]));regions=[];length=0.
        for a,b in zip(breaks[:-1],breaks[1:]):
            if np.interp(.5*(a+b),xi,cf)<-1e-9:
                regions.append({"start":location(a),"end":location(b),
                                "reaches_trailing_edge":bool(abs(b-xi[-1])<1e-9)})
                length+=b-a
        result[label]={"transition":tr,"separation_regions":regions,
                       "separated_arclength_fraction":float(length/(xi[-1]-xi[0])),
                       "min_cf":float(min(cf)),"max_Hk":float(max(M.post.Hk[ids]))}
    return result


def analyze_stall(solutions):
    """Bracket sampled positive-lift peaks; never interpret failure as stall.

    This is a polar diagnostic, not a post-stall extension or a CL limiter.
    A bracket still needs angle/grid refinement and physical validation.
    """
    ordered=sorted(solutions,key=lambda r:r.alpha)
    good=[r for r in ordered if r.converged and np.isfinite(r.cl)]
    best=max(good,key=lambda r:r.cl) if good else None
    peaks=[]
    for left,mid,right in zip(ordered,ordered[1:],ordered[2:]):
        if not all(r.converged and np.isfinite(r.cl) for r in (left,mid,right)):continue
        if not left.alpha<mid.alpha<right.alpha:continue
        if mid.cl>left.cl+1e-3 and mid.cl>right.cl+1e-3:
            peaks.append({"alpha_sample_deg":mid.alpha,"CL_sample":mid.cl,
                          "alpha_bracket_deg":[left.alpha,right.alpha]})
    return {"status":"peak_bracketed" if peaks else "not_bracketed",
            "sampled_maximum":None if best is None else {"alpha_deg":best.alpha,"CL":best.cl},
            "peak_candidates":peaks,
            "failed_angles_deg":[r.alpha for r in ordered if not r.converged],
            "interpretation":"A sampled lift peak is a stall candidate, not a validated stall angle. "
              "Refine panels and alpha spacing. Numerical failure is not stall. "
              "The integral boundary-layer model does not resolve massive separation or unsteady post-stall flow."}


class CoupledSolver:
    """One global Newton system for all foil boundary layers and their wakes."""
    def __init__(self,elements,moment_reference=(.25,0.),reynolds=1e6,mach=0.,ncrit=9.,
                 iterations=80,tolerance=1e-7,wake_length=1.,verbose=False,minimum_te_gap=1e-5,
                 continuation_step=1.,min_step=.25,start_alpha=None):
        if not np.isfinite(reynolds) or reynolds<=0: raise ValueError("Re must be positive")
        if not np.isfinite(mach) or not 0<=mach<.7: raise ValueError("Use 0 <= Mach < 0.7; shocks are unsupported")
        if not np.isfinite(ncrit) or ncrit<=0: raise ValueError("Ncrit must be positive")
        if iterations<1 or tolerance<=0 or wake_length<=0: raise ValueError("Invalid solver limits")
        if not np.isfinite([tolerance,wake_length,continuation_step,min_step]).all() or not 0<min_step<=continuation_step:
            raise ValueError("Require finite limits and 0 < min_step <= continuation_step")
        if start_alpha is not None and not np.isfinite(start_alpha):raise ValueError('start_alpha must be finite')
        if sum(len(x) for _,x in elements)>1200: raise ValueError("Viscous limit is 1200 surface nodes; reduce panels")
        validate_geometry(elements)
        self.elements=elements;self.ref=np.asarray(moment_reference)
        if not 1e-7<=minimum_te_gap<=1e-3: raise ValueError("minimum_te_gap must be between 1e-7 and 1e-3 local chords")
        self.models=[new_model(n,x,reynolds,mach,ncrit,self.ref,minimum_te_gap) for n,x in elements]
        self.elements=[(M.geom.name,M.foil.x.T[::-1].copy()) for M in self.models]
        validate_geometry(self.elements)
        self.iterations=iterations;self.tolerance=tolerance;self.wake_length=wake_length
        self.verbose=verbose;self.history=[];self.initialized=False;self.last_alpha=None
        self.continuation_step=continuation_step;self.min_step=min_step;self.last_solution=None
        self.start_alpha=start_alpha
        self.fs=[];start=0
        for M in self.models:
            self.fs.append(slice(start,start+M.foil.N));start+=M.foil.N
        self.nfoil=start
        self._build_inviscid()

    def _build_inviscid(self):
        nf,ne=self.nfoil,len(self.models)
        A=np.zeros((nf+ne,nf+ne));rhs=np.zeros((nf+ne,2))
        self.sharp_rows=[]
        for i,(M,sl) in enumerate(zip(self.models,self.fs)):
            rhs[sl]=np.column_stack((-M.foil.x[1],M.foil.x[0]))
            for j,(P,col) in enumerate(zip(self.models,self.fs)):
                A[sl,col]=vortex_stream_matrix(P.foil.x,M.foil.x,i!=j)
            A[sl,nf+i]=-1
            if abs(te_info(M.foil.x)[1])<1e-12:
                row=sl.stop-1; A[row]=0;rhs[row]=0
                ids=np.array([0,1,2,M.foil.N-3,M.foil.N-2,M.foil.N-1])+sl.start
                A[row,ids]=[1,-2,1,-1,2,-1]
                self.sharp_rows.append(row)
            A[nf+i,[sl.start,sl.stop-1]]=1
        self.lu=linalg.lu_factor(A)
        self.gref=linalg.lu_solve(self.lu,rhs)[:nf]
        self.aic=A;self.rhs=rhs
        for M,sl in zip(self.models,self.fs):M.isol.gamref=self.gref[sl].copy()

    def _velocity(self,points,alpha):
        a=np.deg2rad(alpha); u=np.array([np.cos(a),np.sin(a)])
        v=np.tile(u,(points.shape[1],1))
        for M,sl in zip(self.models,self.fs):
            v+=np.einsum('ijk,k->ij',vortex_velocity_matrix(M.foil.x,points),self.gref[sl]@u)
        return v.T

    def _build_wakes(self,alpha):
        from matplotlib.path import Path as Polygon
        stream=np.array([np.cos(np.deg2rad(alpha)),np.sin(np.deg2rad(alpha))])
        te_points=[.5*(M.foil.x[:,0]+M.foil.x[:,-1]) for M in self.models]
        end_plane=max(stream@pt for pt in te_points)+self.wake_length
        lengths=[max(self.wake_length,end_plane-stream@pt) for pt in te_points]
        # Keep node counts fixed across alpha continuation. Use a geometric
        # bound on TE separation; otherwise a ceiling change would remap U.
        counts=[max(16,int(np.ceil(M.foil.N/10+10+20*np.log(max(1.,(self.wake_length+
                    max(np.linalg.norm(pt-other) for other in te_points))/M.local_chord)))))
                for M,pt in zip(self.models,te_points)]
        if self.nfoil+sum(counts)>1800:
            raise ValueError('Surface plus wake node limit (1800) exceeded; reduce panels or wake_length')
        for idx,(M,wlength,nw) in enumerate(zip(self.models,lengths,counts)):
            N=M.foil.N
            ds0=.5*((M.foil.s[1]-M.foil.s[0])+(M.foil.s[-1]-M.foil.s[-2]))
            distances=mf.space_geom(ds0,wlength,nw)
            t=te_info(M.foil.x)[0]
            X=np.zeros((2,nw));T=np.zeros_like(X)
            X[:,0]=.5*(M.foil.x[:,0]+M.foil.x[:,-1])+t*(1e-6*M.local_chord)
            for j in range(nw-1):
                v=self._velocity(X[:,j:j+1],alpha)[:,0]
                if np.linalg.norm(v)<1e-5: raise RuntimeError("Wake encountered stagnation")
                t0=v/np.linalg.norm(v); ds=distances[j+1]-distances[j]
                pred=X[:,j]+ds*t0
                v2=self._velocity(pred[:,None],alpha)[:,0];t1=v2/np.linalg.norm(v2)
                X[:,j+1]=X[:,j]+.5*ds*(t0+t1);T[:,j]=t0;T[:,j+1]=t1
            for k,other in enumerate(self.models):
                if k!=idx:
                    a,b=_closed_segments(other.foil.x.T)
                    if (Polygon(other.foil.x.T).contains_points(X.T).any()
                            or _intersections(X[:,:-1].T,X[:,1:].T,a,b).any()):
                        raise RuntimeError(f"Wake from {M.geom.name} intersects {other.geom.name}; this topology is unsupported")
            M.wake.N=nw;M.wake.x=X;M.wake.t=T
            M.wake.s=M.foil.s[-1]+np.r_[0,np.cumsum(np.linalg.norm(np.diff(X,axis=1),axis=0))]
            refs=[]
            for a in (0.,90.):refs.append(np.sum(self._velocity(X,a)*T,axis=0))
            M.isol.uewiref=np.column_stack(refs)
        self.ss=[];self.sources=[];isys=isrc=0
        for M in self.models:
            ns=M.foil.N+M.wake.N;nc=ns-2
            self.ss.append(slice(isys,isys+ns));self.sources.append(slice(isrc,isrc+nc))
            isys+=ns;isrc+=nc
        self.nsys,self.nsource=isys,isrc

    def _build_source_coupling(self):
        nf,ne=self.nfoil,len(self.models)
        B=np.zeros((nf+ne,self.nsource))
        for i,(M,rows) in enumerate(zip(self.models,self.fs)):
            for j,(P,cols) in enumerate(zip(self.models,self.sources)):
                B[rows,cols]=source_stream_matrix(P,M.foil.x,i!=j)
        B[self.sharp_rows]=0
        response=-linalg.lu_solve(self.lu,B)[:nf]
        D=np.zeros((self.nsys,self.nsource))
        for i,(M,foil,sys) in enumerate(zip(self.models,self.fs,self.ss)):
            D[sys.start:sys.start+M.foil.N]=response[foil]
            C=np.zeros((M.wake.N,nf));S=np.zeros((M.wake.N,self.nsource))
            for j,(P,fsl,src) in enumerate(zip(self.models,self.fs,self.sources)):
                vel=vortex_velocity_matrix(P.foil.x,M.wake.x)
                C[:,fsl]=np.einsum('ijk,ji->ik',vel,M.wake.t)
                S[:,src]=source_velocity_matrix(P,M.wake.x,M.wake.t,i==j)
            wake=C@response+S
            wake[0]=response[foil.stop-1]
            D[sys.start+M.foil.N:sys.stop]=wake
        self.ue_sigma=D
        for M,ss,sc in zip(self.models,self.ss,self.sources):
            M.vsol.ue_sigma=D[ss,sc].copy()
        self._mass_operator()

    def _mass_operator(self):
        dp=np.zeros((self.nsource,self.nsys));sgn=np.ones(self.nsys)
        for M,ss,sc in zip(self.models,self.ss,self.sources):
            mf.rebuild_ue_m(M) # also supplies the local matrix needed by stagpoint_move
            dp[sc,ss]=M.vsol.sigma_m.toarray()
            sgn[ss.start:ss.start+M.foil.N]=M.isol.sgnue
        self.D=sgn[:,None]*(self.ue_sigma@dp)

    def _initialize(self,alpha,cold):
        u=np.array([np.cos(np.deg2rad(alpha)),np.sin(np.deg2rad(alpha))])
        previous=[]
        for M in self.models:
            previous.append(None if cold else mf.get_ueinv(M).copy())
            M.oper.alpha=alpha;M.isol.gam=M.isol.gamref@u
        self._build_wakes(alpha)
        for M in self.models:
            if cold:
                mf.stagpoint_find(M);mf.identify_surfaces(M)
            else:
                M.isol.xi=np.r_[abs(M.foil.s-M.isol.sstag),M.wake.s-M.isol.sstag]
            if abs(te_info(M.foil.x)[1])<1e-12:M.vsol.wgap=np.zeros(M.wake.N)
            else:mf.set_wake_gap(M)
        self._build_source_coupling()
        for M,oldinv in zip(self.models,previous):
            if cold:
                M.oper.initbl=True
                M.initialization='inviscid marching'
                try:
                    mf.init_boundary_layer(M)
                except (FloatingPointError,ValueError,AssertionError):
                    # A very adverse inviscid gradient can break the one-way
                    # marching initializer. Seed with a gentle accelerating
                    # profile, then restore the exact global inviscid field.
                    # Only the initial guess changes; the converged equations
                    # and residual tolerances remain identical.
                    gref=M.isol.gamref.copy();wref=M.isol.uewiref.copy()
                    xi=M.isol.xi[:M.foil.N]
                    seed=np.sqrt(xi/(xi+.005*M.local_chord))
                    M.isol.gamref=np.outer(M.isol.sgnue*seed,u)
                    M.isol.uewiref=np.outer(np.ones(M.wake.N),u)
                    try:mf.init_boundary_layer(M)
                    finally:M.isol.gamref=gref;M.isol.uewiref=wref
                    M.initialization='gentle-profile seed; exact global equations restored'
            else:
                M.glob.U[3]+=mf.get_ueinv(M)-oldinv
            mf.stagpoint_move(M)
        self._mass_operator()

    def _system(self):
        ns=self.nsys
        R=np.zeros(4*ns);J=np.zeros((4*ns,4*ns))
        state=np.hstack([M.glob.U for M in self.models])
        ue,ds=state[3],state[1]
        inv=np.concatenate([mf.get_ueinv(M) for M in self.models])
        for M,sl in zip(self.models,self.ss):
            r,j=dense_boundary_system(M)
            M.glob.R=r;M.glob.R_U=j
            R[3*sl.start:3*sl.stop]=r
            J[3*sl.start:3*sl.stop,4*sl.start:4*sl.stop]=j
        R[3*ns:]=ue-inv-self.D@(ue*ds)
        J[3*ns:,3::4]=np.eye(ns)-self.D*ds[None,:]
        J[3*ns:,1::4]=-self.D*ue[None,:]
        return R,J

    def _newton(self):
        self.history=[]
        for it in range(self.iterations):
            R,J=self._system()
            blnorm=float(np.linalg.norm(R[:3*self.nsys]));venorm=float(np.max(abs(R[3*self.nsys:])))
            self.history.append({"iteration":it,"bl_residual_l2":blnorm,"coupling_residual_max":venorm})
            if self.verbose:print(f"Newton {it:02d}: BL={blnorm:.3e}, interaction={venorm:.3e}",flush=True)
            if blnorm<self.tolerance and venorm<self.tolerance:return True
            if not np.isfinite(R).all() or not np.isfinite(J).all():return False
            # Thickness, amplification and velocity have different scales.
            # Equilibration improves the numerical linear solve without
            # changing the nonlinear equations or their convergence criteria.
            U=np.hstack([M.glob.U for M in self.models])
            scale=np.maximum(abs(U),np.array([1e-5,1e-5,.05,1.])[:,None]).ravel(order='F')
            J*=scale[None,:]
            rowscale=np.maximum(np.max(abs(J),axis=1),1e-30)
            J/=rowscale[:,None]
            with warnings.catch_warnings():
                warnings.simplefilter('error',linalg.LinAlgWarning)
                d=-linalg.solve(J,R/rowscale,check_finite=False)*scale
            # Native physical step limits, with ONE common relaxation for all bodies.
            old=[];omegas=[]
            for M,sl in zip(self.models,self.ss):
                old.append(M.glob.U.copy())
                M.glob.dU=d[4*sl.start:4*sl.stop].reshape(4,-1,order='F');M.glob.dalpha=0.
                mf.update_state(M)
                du=M.glob.dU;good=abs(du)>1e-10
                ratio=(M.glob.U-old[-1])[good]/du[good]
                omegas.append(float(np.clip(np.median(ratio),0,1)) if ratio.size else 1.)
            om=min(omegas)
            for M,prev in zip(self.models,old):
                M.glob.U=prev
                M.glob.dU*=om
                mf.update_state(M) # apply common step plus the native physical repairs
                mf.stagpoint_move(M);mf.update_transition(M)
            self._mass_operator()
        return False

    def _solve_point(self,alpha):
        if not np.isfinite(alpha):raise ValueError("alpha must be finite")
        start=time.perf_counter();reason="";ok=False;self.history=[]
        cold=not self.initialized
        try:
            with np.errstate(divide='raise',invalid='raise',over='raise'):
                self._initialize(float(alpha),cold)
                ok=self._newton()
        except (ValueError,AssertionError,FloatingPointError,linalg.LinAlgError,linalg.LinAlgWarning,RuntimeError) as exc:
            reason=f"{type(exc).__name__}: {exc}"
            if self.verbose:print(reason,flush=True)
        self.initialized=ok;self.last_alpha=float(alpha) if ok else None
        if not ok and not reason:reason="Newton iteration limit reached"
        residual=max(self.history[-1]['bl_residual_l2'],self.history[-1]['coupling_residual_max']) if self.history else math.inf
        diag={"model":"Newton-coupled mfoil boundary layers, multi-element potential flow",
              "iterations":len(self.history),"history":self.history.copy(),"failure_reason":reason,
              "wall_seconds":time.perf_counter()-start,"ncrit":self.models[0].param.ncrit,
              "reynolds":self.models[0].oper.Re,"mach":self.models[0].oper.Ma,
              "te_gap_added_over_local_chord":{M.geom.name:M.regularization for M in self.models},
              "initialization":{M.geom.name:getattr(M,'initialization','not reached') for M in self.models},
              "cd_definition":"Sum of Squire-Young wake momentum drag estimates; independent, non-merging wakes",
              "wake_length_definition":"Approximate extension beyond the aftmost TE, measured along freestream in reference chords",
              "status":"converged" if ok else "nonconverged"}
        if not ok:return Solution(float(alpha),math.nan,math.nan,math.nan,None,False,residual,[],diag)
        surfaces=[];cl=cm=cd=0.;diag['boundary_layers']={};diag['element_forces']={}
        diag['local_supersonic_detected']={}
        for M in self.models:
            mf.calc_force(M);mf.get_distributions(M)
            # Convert the upstream velocity-normalized Cf to tau_w/q_infinity
            # at nonzero Mach, including the local density ratio.
            if M.oper.Ma>0:
                for si in range(3):
                    param=mf.build_param(M,si)
                    for j in M.vsol.Is[si]:
                        rho,_=mf.get_rho(M.glob.U[:,j],param)
                        M.post.cf[j]*=rho/M.oper.rho
            cl+=M.post.cl;cm+=M.post.cm;cd+=M.post.cd
            N=M.foil.N
            diag['local_supersonic_detected'][M.geom.name]=bool(M.oper.Ma>0 and min(M.post.cp[:N])<M.param.cps)
            diag['boundary_layers'][M.geom.name]=boundary_layer_summary(M)
            diag['element_forces'][M.geom.name]={"CL":float(M.post.cl),"CM":float(M.post.cm),
                "CD_wake":float(M.post.cd),"CD_friction":float(M.post.cdf)}
            side=np.full(N,'upper',dtype='<U5');side[M.vsol.Is[0]]='lower'
            surfaces.append({"name":M.geom.name,"xy":M.foil.x.T[::-1].copy(),
                "x_over_local_chord":(M.chord_axis@(M.foil.x-M.chord_origin[:,None])/M.local_chord)[::-1].copy(),
                "side":side[::-1].copy(),"ue_over_vinf":np.asarray(M.post.ue[:N])[::-1].copy(),
                "cp":np.asarray(M.post.cp[:N])[::-1].copy(),
                "cf":np.asarray(M.post.cf[:N])[::-1].copy(),
                "theta":M.glob.U[0,:N][::-1].copy(),"delta_star":M.glob.U[1,:N][::-1].copy(),
                "H":np.asarray(M.post.Hk[:N])[::-1].copy(),
                "amplification_or_shear":M.glob.U[2,:N][::-1].copy(),
                "turbulent":np.asarray(M.vsol.turb[:N])[::-1].copy(),"cd_wake":float(M.post.cd),
                "wake":{"x_over_cref":M.wake.x[0].copy(),"y_over_cref":M.wake.x[1].copy(),
                        "theta":M.glob.U[0,N:].copy(),"delta_star":M.glob.U[1,N:].copy(),
                        "ue_over_vinf":np.asarray(M.post.ue[N:]).copy(),"Hk":M.post.Hk[N:].copy()}})
        return Solution(float(alpha),float(cl),float(cm),float(cd),None,True,residual,surfaces,diag)

    def solve(self,alpha):
        """Reach alpha by continuation, retaining the last converged state.

        Internal continuation points are recorded in diagnostics; returned
        forces always refer to the requested angle, or NaN on failure.
        """
        alpha=float(alpha)
        if not np.isfinite(alpha):raise ValueError("alpha must be finite")
        if self.initialized and abs(alpha-self.last_alpha)<1e-12:
            return copy.deepcopy(self.last_solution)
        started=time.perf_counter();trail=[];step=self.continuation_step
        if not self.initialized:
            initial=self.start_alpha if self.start_alpha is not None else (0. if abs(alpha)>self.continuation_step else alpha)
            result=self._solve_point(initial)
            trail.append({"alpha_deg":initial,"converged":result.converged,"CL":result.cl})
            if result.converged:self.last_solution=result
            elif initial!=alpha:
                result=self._solve_point(alpha)
                trail.append({"alpha_deg":alpha,"converged":result.converged,"CL":result.cl})
                if result.converged:self.last_solution=result
            if not result.converged or abs(self.last_alpha-alpha)<1e-12:
                result.alpha=alpha;result.diagnostics['continuation']=trail
                result.diagnostics['wall_seconds']=time.perf_counter()-started
                return result
        while abs(alpha-self.last_alpha)>1e-12:
            base=self.last_alpha
            target=base+np.sign(alpha-base)*min(step,abs(alpha-base))
            checkpoint=copy.deepcopy(self.models)
            result=self._solve_point(target)
            trail.append({"alpha_deg":target,"converged":result.converged,"CL":result.cl})
            if result.converged:
                self.last_solution=result
                step=min(self.continuation_step,step*2)
            else:
                self.models=checkpoint;self.initialized=True;self.last_alpha=base
                attempted=abs(target-base)
                if attempted<=self.min_step+1e-12:
                    result.alpha=alpha
                    result.diagnostics['failure_reason']=f"Continuation failed at {target:g} deg: "+result.diagnostics['failure_reason']
                    break
                step=max(self.min_step,.5*attempted)
        result.diagnostics['continuation']=trail
        result.diagnostics['wall_seconds']=time.perf_counter()-started
        return result
