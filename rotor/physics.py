"""Non-clipping annular momentum closure for *diagnostic* BEMT calculations.

Not calibrated for any APC blade. It rejects regimes where the simple
momentum closure has no admissible real solution rather than hiding
negative blade-element loading through zero-clipping. This corrects the
specific numerical inconsistency in the archived 2026 solver; it does NOT
constitute experimental validation of a physics-based propeller model.
"""
from dataclasses import dataclass
import numpy as np

@dataclass(frozen=True)
class SectionResult:
    converged:bool
    reason:str
    thrust_n:float
    torque_nm:float
    shaft_power_w:float
    negative_thrust_sections:int
    negative_torque_sections:int
    iterations:int


def induction_from_annulus(dTdr,dQdr,r,V,F,rho,omega):
    """Return signed disk-plane axial and tangential induced velocities.

    Uses dT/dr=4*pi*r*rho*F*(V+vi)*vi and
         dQ/dr=4*pi*r*r*rho*F*(V+vi)*vtheta.
    This is the propeller branch with V+vi>0. Negative annular thrust is
    allowed in axial flight only if that branch remains real and admissible.
    """
    r=np.asarray(r,dtype=float)
    dTdr=np.asarray(dTdr,dtype=float)
    dQdr=np.asarray(dQdr,dtype=float)
    F=np.asarray(F,dtype=float)
    disc=V**2+dTdr/(np.pi*r*rho*F)
    if np.any(disc<0):raise ValueError('no real axial momentum solution')
    vi=(-V+np.sqrt(disc))/2
    flow=V+vi
    if np.any(flow<1e-6):raise ValueError('nonpositive axial disk flow')
    swirl=dQdr/(4*np.pi*r*r*rho*F*flow)
    if np.any(np.abs(swirl)>0.8*omega*r):
        raise ValueError('swirl outside low-induction model domain')
    return vi,swirl


def bemt_diagnostic(radius,chord,omega,pitch_deg,blades,axial_mps,
                    polar,n=128,rho=1.225,maxiter=500,tol=1e-6):
    """Constant chord/pitch, user-supplied polar(alpha rad, Re)->(CL,CD).

    For broad negative-loaded/windmill regimes return a model limitation
    rather than falsely claiming steady-state momentum convergence.
    """
    if min(radius,chord,omega,blades,n)<0 or radius==0 or axial_mps<0:
        raise ValueError('positive dimensions and nonnegative axial velocity required')
    edge=np.linspace(.20*radius,radius,n+1);r=(edge[1:]+edge[:-1])/2;dr=np.diff(edge)
    vi=np.full_like(r,.8 if axial_mps>0 else 4.)
    swirl=np.zeros_like(r)
    theta=np.deg2rad(pitch_deg)
    def forces(va,vt):
        phi=np.arctan2(va,vt);w=np.hypot(va,vt)
        CL,CD=polar(theta-phi,rho*w*chord/1.81e-5)
        L=.5*rho*w*w*chord*blades*CL*dr
        D=.5*rho*w*w*chord*blades*CD*dr
        return L*np.cos(phi)-D*np.sin(phi), r*(L*np.sin(phi)+D*np.cos(phi)),phi
    reason='iteration limit';converged=False
    for it in range(1,maxiter+1):
        dt,dq,phi=forces(axial_mps+vi,np.maximum(omega*r-swirl,1e-8))
        sinphi=np.maximum(np.abs(np.sin(phi)),1e-5)
        ftip=.5*blades*(radius-r)/(r*sinphi)
        fhub=.5*blades*(r-.2*radius)/(.2*radius*sinphi)
        F=np.maximum(.05,(2/np.pi*np.arccos(np.exp(-ftip)))*
                           (2/np.pi*np.arccos(np.exp(-fhub))))
        try:
            tiv,tis=induction_from_annulus(dt/dr,dq/dr,r,axial_mps,F,rho,omega)
        except ValueError as ex:
            reason=str(ex)
            break
        next_vi=.75*vi+.25*tiv
        next_swirl=.75*swirl+.25*tis
        change=max(np.max(abs(next_vi-vi)),np.max(abs(next_swirl-swirl)))
        vi,swirl=next_vi,next_swirl
        if change<tol:
            converged=True;reason='converged';break
    dt,dq,_=forces(axial_mps+vi,np.maximum(omega*r-swirl,1e-8))
    thrust=float(sum(dt));torque=float(sum(dq))
    return SectionResult(converged,reason,thrust,torque,torque*omega,
                         int(sum(dt<0)),int(sum(dq<0)),it)
