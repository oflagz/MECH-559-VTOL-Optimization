"""Numba-accelerated implementation of exactly the same signed BEMT closure."""
from functools import lru_cache
import math
import numpy as np
from numba import njit
from .bemt import (SectionPolar,RotorGeometry,Performance,get_reference_geometry,
                   RHO,MU,SPEED_SOUND)

@lru_cache(maxsize=20)
def mesh(n):
    edges=np.linspace(.155,.995,n+1)
    mid=(edges[1:]+edges[:-1])/2
    from .bemt import get_reference_geometry
    ref,cc,tt=get_reference_geometry()
    shape=np.interp(mid,ref,cc);shape/=np.mean(shape)
    return mid.astype(np.float64),np.diff(edges).astype(np.float64),shape.astype(np.float64),np.interp(mid,ref,tt).astype(np.float64)

@njit(cache=True)
def _core(R,c,om,tw,B,theta75,V,pv,mids,widths,shape,refpitch,reference,limit,tol):
    n=len(mids)
    vi=np.full(n,4. if V==0 else .6)
    sw=np.zeros(n)
    dt=np.zeros(n);dq=np.zeros(n)
    reason=0;converged=False;err=1e9
    for iteration in range(1,limit+1):
        vi_next=np.zeros(n);sw_next=np.zeros(n);err=0.
        for k in range(n):
            mu=mids[k];r=R*mu
            chord=R*shape[k]*.181 if reference else c*shape[k]
            theta=(refpitch[k] if reference else theta75+tw*(mu-.75))*math.pi/180
            tang=om*r-sw[k];va=V+vi[k]
            if tang<=1e-6 or va<=1e-8:
                reason=1;break
            w=math.hypot(va,tang);phi=math.atan2(va,tang)
            a=theta-phi; re=max(RHO*w*chord/MU,8000.)
            lr=min(1.3,max(.75,(re/60000)**pv[6]));aa=a-pv[1]*math.pi/180
            cl=pv[2]*math.tanh(pv[0]*lr*aa/pv[2]);cd=pv[3]*min(1.5,max(.65,(60000/re)**pv[7]))+pv[4]*cl*cl
            cd+=pv[5]*max(0.,abs(aa)-18*math.pi/180)**2
            L=.5*RHO*w*w*B*chord*cl;D=.5*RHO*w*w*B*chord*cd
            dt[k]=L*math.cos(phi)-D*math.sin(phi)
            dq[k]=r*(L*math.sin(phi)+D*math.cos(phi))
            sf=max(abs(math.sin(phi)),1e-4)
            tip=(B/2)*(R-r)/(r*sf);hub=(B/2)*(r-.15*R)/(.15*R*sf)
            Ft=(2/math.pi)*math.acos(math.exp(-max(tip,0)))
            Fh=(2/math.pi)*math.acos(math.exp(-max(hub,0)))
            F=max(Ft*Fh,.06)
            disc=V*V+dt[k]/(math.pi*r*RHO*F)
            if disc<0:
                reason=2;break
            tiv=(-V+math.sqrt(disc))/2
            flow=V+tiv
            if flow<=1e-8:
                reason=3;break
            tsi=dq[k]/(4*math.pi*RHO*F*r*r*flow)
            if abs(tsi)>.8*om*r:
                reason=4;break
            vi_next[k]=tiv;sw_next[k]=tsi
            err=max(err,abs(tiv-vi[k]),abs(tsi-sw[k]))
        if reason>0:break
        for k in range(n):
            vi[k]=.65*vi[k]+.35*vi_next[k]
            sw[k]=.65*sw[k]+.35*sw_next[k]
        if err<tol:
            converged=True;break
                                                                                 
    T=0.;Q=0.;min_alpha=1e9;max_alpha=-1e9;min_re=1e9;max_re=0.;negative=0
    for k in range(n):
        mu=mids[k];r=R*mu
        chord=R*shape[k]*.181 if reference else c*shape[k]
        theta=(refpitch[k] if reference else theta75+tw*(mu-.75))*math.pi/180
        tang=max(om*r-sw[k],1e-8);va=V+vi[k]
        w=math.hypot(va,tang);phi=math.atan2(va,tang);a=theta-phi
        re=max(RHO*w*chord/MU,8000.);lr=min(1.3,max(.75,(re/60000)**pv[6]));aa=a-pv[1]*math.pi/180
        cl=pv[2]*math.tanh(pv[0]*lr*aa/pv[2]);cd=pv[3]*min(1.5,max(.65,(60000/re)**pv[7]))+pv[4]*cl*cl
        cd+=pv[5]*max(0.,abs(aa)-18*math.pi/180)**2
        L=.5*RHO*w*w*B*chord*cl;D=.5*RHO*w*w*B*chord*cd
        dT=L*math.cos(phi)-D*math.sin(phi)
        dQ=r*(L*math.sin(phi)+D*math.cos(phi))
        T+=dT*R*widths[k];Q+=dQ*R*widths[k]
        min_alpha=min(min_alpha,a);max_alpha=max(max_alpha,a)
        min_re=min(min_re,re);max_re=max(max_re,re)
        if dT<0:negative+=1
    P=Q*om;nrev=om/(2*math.pi);D=2*R
    ct=T/(RHO*nrev*nrev*D**4);cp=P/(RHO*nrev**3*D**5)
    fm=T**1.5/math.sqrt(2*RHO*math.pi*R*R)/P if T>0 and P>0 else math.nan
    eta=T*V/P if P>0 else math.nan
    return (1. if converged else 0.,float(reason),T,Q,P,ct,cp,fm,eta,
            min_alpha*180/math.pi,max_alpha*180/math.pi,min_re,max_re,err,float(iteration),negative/n)

def fast_performance(g,velocity,polar=SectionPolar(),n_sections=32,reference=False):
    mids,widths,shape,refpitch=mesh(n_sections)
    pv=np.array([polar.cl_slope,polar.alpha_zero_deg,polar.cl_max,polar.cd0,polar.drag_k,
                 polar.stall_drag,polar.re_slope,polar.re_cd],dtype=np.float64)
    z=_core(g.radius,g.mean_chord,g.omega,g.twist_deg,g.blade_count,g.theta75_deg,
            velocity,pv,mids,widths,shape,refpitch,reference,120,2e-5)
    reason={0:'iteration limit',1:'nonpositive disk throughflow',2:'negative axial momentum discriminant',
            3:'nonpositive momentum throughflow',4:'swirl outside low-induction regime'}
    msg='converged' if bool(z[0]) else reason[int(z[1])]
    return Performance(bool(z[0]),msg,float(z[2]),float(z[3]),float(z[4]),float(z[5]),float(z[6]),
                       float(z[7]),float(z[8]),float(z[9]),float(z[10]),float(z[11]),float(z[12]),
                       float(z[13]),int(z[14]),float(z[15]))
