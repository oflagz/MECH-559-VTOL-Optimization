"""UIUC APC 10x7 reference-propeller BEMT diagnostic.

Only the geometry and C_T/C_P benchmark come from measurements.
Blade-section aerodynamic POLARS remain uncalibrated synthetic estimates.
Model is not suitable for windmilling or high-induction reverse-flow regimes.
"""
from pathlib import Path
import numpy as np
from .model import polar

DATA=Path(__file__).resolve().parents[1]/'data'/'uiuc'
RADIUS=0.127  # APC 10-inch diameter / 2, in metres
D=2*RADIUS
NBLADES=2
RHO=1.225

def load_geometry(path=DATA/'apcsf_10x7_geometry.csv'):
    a=np.genfromtxt(path,delimiter=',',names=True)
    assert np.all(np.diff(a['r_over_R'])>0)
    return a

def predict_reference(rpm, J=0., n_elements=120, rho=RHO, use_tip_hub_loss=True):
    """Blade geometry from UIUC; illustrative surrogate airfoil only.

    Returns nondimensional CT, CP and efficiency, or convergence diagnostics.
    No free parameters are fitted against UIUC thrust or power observations.
    """
    geom=load_geometry()
    n=rpm/60.
    omega=2*np.pi*n
    V=J*n*D
    rr=np.linspace(float(geom['r_over_R'][0]),1.,n_elements+1)
    r=RADIUS*(rr[:-1]+rr[1:])/2
    dr=RADIUS*np.diff(rr)
    x=r/RADIUS
    chord=RADIUS*np.interp(x,geom['r_over_R'],geom['c_over_R'])
    pitch=np.deg2rad(np.interp(x,geom['r_over_R'],geom['pitch_deg']))
    vi=np.full_like(r,1.0)
    swirl=np.zeros_like(r)
    convergence=False
    relax=.27

    def loads(vi,swirl):
        va=V+vi
        vt=np.maximum(omega*r-swirl,1e-8)
        phi=np.arctan2(va,vt)
        vrel=np.hypot(va,vt)
        re=rho*vrel*chord/1.81e-5
        cl,cd=polar(pitch-phi,re)
        dL=.5*rho*vrel**2*chord*NBLADES*cl*dr
        dD=.5*rho*vrel**2*chord*NBLADES*cd*dr
        dT=dL*np.cos(phi)-dD*np.sin(phi)
        dQ=r*(dL*np.sin(phi)+dD*np.cos(phi))
        return dT,dQ,phi,re

    for it in range(1,351):
        dT,dQ,phi,re=loads(vi,swirl)
        if use_tip_hub_loss:
            sp=np.maximum(np.abs(np.sin(phi)),.02)
            ft=NBLADES/2*(RADIUS-r)/(r*sp)
            fh=NBLADES/2*(r-rr[0]*RADIUS)/(rr[0]*RADIUS*sp)
            F=np.maximum(.05,2/np.pi*np.arccos(np.exp(-ft)) *
                              2/np.pi*np.arccos(np.exp(-fh)))
        else:F=np.ones_like(r)
        dtdr=np.maximum(dT/dr,0)
        target_vi=(np.sqrt(V*V+dtdr/(np.pi*rho*F*r))-V)/2
        target_swirl=np.maximum(dQ/dr,0)/(4*np.pi*rho*F*r*r*np.maximum(V+target_vi,.05))
        target_swirl=np.minimum(target_swirl,.6*omega*r)
        vi_new=(1-relax)*vi+relax*target_vi
        swirl_new=(1-relax)*swirl+relax*target_swirl
        err=max(np.max(np.abs(vi_new-vi)),np.max(np.abs(swirl_new-swirl)))
        vi,swirl=vi_new,swirl_new
        if err<5e-6:
            convergence=True
            break
    dT,dQ,phi,re=loads(vi,swirl)
    thrust=float(np.sum(dT)); torque=float(np.sum(dQ)); shaft=torque*omega
    ct=thrust/(rho*n*n*D**4)
    cp=shaft/(rho*n**3*D**5)
    eta=J*ct/cp if cp>0 else np.nan
    return {'CT':ct,'CP':cp,'eta':eta,'thrust_n':thrust,'shaft_w':shaft,
            'converged':convergence,'iterations':it,'re_min':float(min(re)),
            're_max':float(max(re))}
