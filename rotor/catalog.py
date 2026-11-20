"""Source-specific propeller interpolation. Never extrapolate the experimental envelope."""
from dataclasses import dataclass
from pathlib import Path
import csv
import numpy as np
from scipy.interpolate import PchipInterpolator

DATA=Path(__file__).resolve().parents[1]/'data'/'rotor'
@dataclass(frozen=True)
class Propeller:
    name:str
    diameter_m:float
    mass_kg:float
    source:str
    @property
    def radius(self):return self.diameter_m/2

# Published APC product weights; exactly 2 blades in the measured UIUC set.
CATALOG={
    '16x12':Propeller('16x12',16*.0254,1.83*0.028349523125,'https://www.apcprop.com/product/16x12e/'),
    '18x10':Propeller('18x10',18*.0254,2.57*0.028349523125,'https://www.apcprop.com/product/18x10e/'),
    '18x12':Propeller('18x12',18*.0254,2.61*0.028349523125,'https://www.apcprop.com/product/18x12e/'),
    '20x10':Propeller('20x10',20*.0254,3.39*0.028349523125,'https://www.apcprop.com/product/20x10e/'),
}

def load_rows(filename):
    with (DATA/filename).open(newline='') as f:
        return [{k:float(v) if k not in ('model',) else v for k,v in r.items()} for r in csv.DictReader(f)]

STATIC=load_rows('static.csv')
FORWARD=load_rows('forward.csv')

class MeasuredPropeller:
    """Measured, zero-extrapolation CT/CP lookup. rpm in revolutions/minute."""
    def __init__(self, name, rho=1.225):
        self.prop=CATALOG[name]
        self.rho=rho
        sta=sorted([r for r in STATIC if r['model']==name],key=lambda r:r['rpm'])
        self.static_rpm=np.array([r['rpm'] for r in sta])
        self.static_ct=PchipInterpolator(self.static_rpm,[r['CT'] for r in sta],extrapolate=False)
        self.static_cp=PchipInterpolator(self.static_rpm,[r['CP'] for r in sta],extrapolate=False)
        self.forward={}
        for rpm in sorted({r['rpm'] for r in FORWARD if r['model']==name}):
            rows=sorted([r for r in FORWARD if r['model']==name and r['rpm']==rpm],key=lambda r:r['J'])
            J=np.array([r['J'] for r in rows]);CT=np.array([r['CT'] for r in rows]);CP=np.array([r['CP'] for r in rows])
            assert np.all(np.diff(J)>0)
            self.forward[rpm]=(J,PchipInterpolator(J,CT,extrapolate=False),PchipInterpolator(J,CP,extrapolate=False))
        self.forward_rpms=np.array(sorted(self.forward))

    def _dimension(self,rpm,ct,cp,velocity):
        D=self.prop.diameter_m;n=rpm/60
        thrust=float(self.rho*n*n*D**4*ct)
        shaft=float(self.rho*n**3*D**5*cp)
        torque=shaft/(2*np.pi*n)
        return dict(thrust_n=thrust,shaft_w=shaft,torque_nm=torque,rpm=float(rpm),
                    ct=float(ct),cp=float(cp),advance_ratio=velocity/(n*D) if n>0 else 0,
                    efficiency=thrust*velocity/shaft if shaft>0 else float('nan'),
                    tip_mach=np.hypot(np.pi*D*n,velocity)/340)

    def static(self,rpm):
        if not self.static_rpm[0]<=rpm<=self.static_rpm[-1]:return None
        return self._dimension(rpm,self.static_ct(rpm),self.static_cp(rpm),0.0)

    def axial(self,rpm,velocity):
        if rpm<self.forward_rpms[0] or rpm>self.forward_rpms[-1] or velocity<=0:return None
        hi=np.searchsorted(self.forward_rpms,rpm,side='left')
        if hi==len(self.forward_rpms):return None
        if abs(rpm-self.forward_rpms[hi])<=1e-8:lo=hi
        elif hi==0:return None
        else:lo=hi-1
        D=self.prop.diameter_m
        J=velocity/(rpm/60*D)
        def ev(i):
            A,ct,cp=self.forward[self.forward_rpms[i]]
            if J<A[0] or J>A[-1]:return None
            return float(ct(J)),float(cp(J))
        left=ev(lo);right=ev(hi)
        if left is None or right is None:return None
        a=0 if lo==hi else (rpm-self.forward_rpms[lo])/(self.forward_rpms[hi]-self.forward_rpms[lo])
        ct=(1-a)*left[0]+a*right[0]
        cp=(1-a)*left[1]+a*right[1]
        if ct<=0 or cp<=0:return None
        return self._dimension(rpm,ct,cp,velocity)

    def coverage(self,velocity):
        grid=np.linspace(self.forward_rpms[0],self.forward_rpms[-1],1801)
        valid=[float(x) for x in grid if self.axial(float(x),velocity) is not None]
        return (min(valid),max(valid)) if valid else None
