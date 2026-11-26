"""Independent signed-induction blade-element/momentum-theory rotor model."""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import numpy as np

RHO=1.225
MU=1.81e-5
SPEED_SOUND=340.
ROOT=Path(__file__).resolve().parents[1]

@dataclass(frozen=True)
class SectionPolar:
    cl_slope:float=5.7
    alpha_zero_deg:float=-1.0
    cl_max:float=1.25
    cd0:float=0.024
    drag_k:float=0.035
    stall_drag:float=0.8
                                                              
    re_slope:float=0.08
    re_cd:float=0.12
    def lift_drag(self,alpha,re):
        re=np.maximum(re,8000.)
        lr=np.clip((re/60000.)**self.re_slope,0.75,1.3)
        a=alpha-np.deg2rad(self.alpha_zero_deg)
        cl=self.cl_max*np.tanh(self.cl_slope*lr*a/self.cl_max)
        cd=self.cd0*np.clip((60000./re)**self.re_cd,0.65,1.5)+self.drag_k*cl**2
        cd+=self.stall_drag*np.maximum(np.abs(a)-np.deg2rad(18.),0.)**2
        return cl,cd

@dataclass(frozen=True)
class RotorGeometry:
    radius:float
    mean_chord:float
    omega:float
    twist_deg:float
    blade_count:int
    theta75_deg:float
    rpm:float|None=None

@dataclass
class Performance:
    valid:bool
    reason:str
    thrust:float
    torque:float
    shaft_power:float
    ct:float
    cp:float
    fm:float
    eta:float
    min_alpha_deg:float
    max_alpha_deg:float
    min_re:float
    max_re:float
    max_induction_residual:float
    iterations:int
    negative_load_fraction:float


def get_reference_geometry():
    path=ROOT/'data/uiuc/apcsf_10x7_geometry.csv'
    arr=np.genfromtxt(path,delimiter=',',names=True)
    return arr['mu'],arr['chord_over_R'],arr['pitch_deg']


def chord_and_pitch(g:RotorGeometry, mu, reference=False):
    if reference:
        rr,cc,tt=get_reference_geometry()
        return (np.interp(mu,rr,cc)*g.radius,
                np.deg2rad(np.interp(mu,rr,tt)))
                                                                                
    ref_mu,ref_chord,_=get_reference_geometry()
    shape=np.interp(mu,ref_mu,ref_chord)
    shape/=np.mean(shape)
    return g.mean_chord*shape,np.deg2rad(g.theta75_deg+g.twist_deg*(mu-.75))


def performance(g:RotorGeometry, axial_mps:float,polar:SectionPolar=SectionPolar(),n_sections=48,reference=False,
                max_iter=120,tol=2e-5):
    if (g.radius<=0 or g.mean_chord<=0 or g.omega<=0 or g.blade_count not in (2,3,4) or axial_mps<0):
        raise ValueError('invalid physical arguments')
    edges=np.linspace(.155, .995, n_sections+1)
    mu=.5*(edges[1:]+edges[:-1]); r=g.radius*mu
    dr=g.radius*np.diff(edges)
    chord,theta=chord_and_pitch(g,mu,reference)
                                                                              
    vi=np.full_like(r,4.0 if axial_mps==0 else .6)
    swirl=np.zeros_like(r)
    omega=g.omega;rho=RHO;B=g.blade_count
    reason='iteration limit';valid=False
    for iteration in range(1,max_iter+1):
        tangential=omega*r-swirl
        if np.min(tangential)<=1e-6:
            reason='swirl exceeds blade tangential velocity';break
        axial=axial_mps+vi
        if np.min(axial)<=1e-8:
            reason='nonpositive disk throughflow';break
        W=np.hypot(axial,tangential)
        phi=np.arctan2(axial,tangential)
        alpha=theta-phi
        re=rho*W*chord/MU
        cl,cd=polar.lift_drag(alpha,re)
        L=.5*rho*W*W*B*chord*cl
        D=.5*rho*W*W*B*chord*cd
        dt=L*np.cos(phi)-D*np.sin(phi)
        dq=r*(L*np.sin(phi)+D*np.cos(phi))
        fs=np.maximum(np.abs(np.sin(phi)),1e-4)
                                                                           
        ftip=(B/2)*(g.radius-r)/(r*fs)
        fhub=(B/2)*(r-.15*g.radius)/(.15*g.radius*fs)
        Ft=(2/np.pi)*np.arccos(np.exp(-np.maximum(ftip,0)))
        Fh=(2/np.pi)*np.arccos(np.exp(-np.maximum(fhub,0)))
        F=np.maximum(Ft*Fh,.06)
        discriminant=axial_mps**2+dt/(np.pi*r*rho*F)
        if np.min(discriminant)<0:
            reason='negative axial momentum discriminant';break
        vi_new=(-axial_mps+np.sqrt(discriminant))/2
        disk_flow=axial_mps+vi_new
        if np.min(disk_flow)<=1e-8:
            reason='nonpositive momentum throughflow';break
        swirl_new=dq/(4*np.pi*rho*F*r*r*disk_flow)
        if np.max(np.abs(swirl_new)/(omega*r))>.8:
            reason='swirl outside low-induction regime';break
        residual=max(np.max(np.abs(vi_new-vi)),np.max(np.abs(swirl_new-swirl)))
        vi=.65*vi+.35*vi_new
        swirl=.65*swirl+.35*swirl_new
        if residual<tol:
            valid=True;reason='converged';break
                                                                                      
    tangential=np.maximum(omega*r-swirl,1e-8);axial=axial_mps+vi
    W=np.hypot(axial,tangential);phi=np.arctan2(axial,tangential)
    alpha=theta-phi;re=rho*W*chord/MU
    cl,cd=polar.lift_drag(alpha,re)
    L=.5*rho*W*W*B*chord*cl
    D=.5*rho*W*W*B*chord*cd
    dt=L*np.cos(phi)-D*np.sin(phi)
    dq=r*(L*np.sin(phi)+D*np.cos(phi))
    T=float(np.dot(dt,dr));Q=float(np.dot(dq,dr));P=Q*omega
    n=omega/(2*np.pi);diam=2*g.radius
    ct=T/(rho*n*n*diam**4);cp=P/(rho*n**3*diam**5)
    fm=T**1.5/(np.sqrt(2*rho*np.pi*g.radius**2)*P) if T>0 and P>0 else float('nan')
    eta=T*axial_mps/P if P>0 else float('nan')
    return Performance(valid,reason,T,Q,P,ct,cp,fm,eta,float(np.rad2deg(alpha.min())),
                       float(np.rad2deg(alpha.max())),float(re.min()),float(re.max()),
                       float(residual) if 'residual' in locals() else float('inf'),iteration,float(np.mean(dt<0)))


def reference_performance(rpm,j,polar=SectionPolar(),n_sections=48):
    radius=.127;omega=rpm*2*np.pi/60
    g=RotorGeometry(radius,.022,omega,-10,2,14.38)
    axial=j*(rpm/60)*2*radius
    return performance(g,axial,polar,n_sections=n_sections,reference=True)
