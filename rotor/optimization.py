"""Discrete measured-hardware selection; exact hover/cruise trim; SLSQP cross-check.

No undocumented variable interpolation: radius/pitch/blade count are catalogue
choices. Within each candidate, hover and cruise RPM are continuous variables.
"""
from dataclasses import dataclass, replace, asdict
import numpy as np
from scipy.optimize import brentq, minimize
from .catalog import MeasuredPropeller,CATALOG

@dataclass(frozen=True)
class Mission:
    base_aircraft_mass_kg:float=9.5 # 5-kg payload included; excludes 4 propeller blade masses.
    payload_mass_kg:float=5.0
    n_rotors:int=4
    rho:float=1.225
    gravity:float=9.81
    download_fraction:float=.05
    reserve_fraction:float=.15
    cruise_drag_n:float=7.0
    cruise_mps:float=20.
    distance_m:float=40000.
    hover_s:float=300.
    eta_drive:float=.85  # illustrative; replace using Damian's map
    motor_torque_max_nm:float=1.5
    total_electric_peak_limit_w:float=2400.
    mount_spacing_m:float=.65
    clearance_m:float=.05
    tip_mach_limit:float=.65
    def __post_init__(self):
        if self.base_aircraft_mass_kg<self.payload_mass_kg:raise ValueError('payload exceeds takeoff mass')
        if self.cruise_mps<=0:raise ValueError('positive cruise speed required')


def find_roots(f,target,a,b,n=800):
    grid=np.linspace(a,b,n)
    vals=[]
    for g in grid:
        d=f(float(g));vals.append(np.nan if d is None else d['thrust_n']-target)
    roots=[]
    for i in range(n-1):
        p,q=vals[i],vals[i+1]
        if not(np.isfinite(p) and np.isfinite(q)):continue
        if p==0:roots.append(grid[i])
        elif p*q<0:
            roots.append(brentq(lambda x:f(x)['thrust_n']-target,grid[i],grid[i+1],xtol=1e-8))
    return sorted(set(round(float(x),8) for x in roots))


def solve_catalog(name, m=Mission()):
    prop=MeasuredPropeller(name,m.rho)
    D=prop.prop.diameter_m
    result=dict(propeller=name,diameter_m=D,blade_count=2,prop_mass_kg=prop.prop.mass_kg,
                aircraft_mass_kg=m.base_aircraft_mass_kg+m.n_rotors*prop.prop.mass_kg,
                cruise_mps=m.cruise_mps,reason='',feasible=False)
    if D+m.clearance_m>m.mount_spacing_m:
        result['reason']='mount clearance';return result
    th=(1+m.download_fraction)*result['aircraft_mass_kg']*m.gravity/m.n_rotors
    tc=m.cruise_drag_n/m.n_rotors
    result.update(required_hover_n=th,required_cruise_n=tc)
    roots_h=find_roots(prop.static,th,prop.static_rpm[0],prop.static_rpm[-1])
    heads=find_roots(prop.static,(1+m.reserve_fraction)*th,prop.static_rpm[0],prop.static_rpm[-1])
    if not roots_h or not heads:
        result['reason']='no measured hover trim/headroom';return result
    roots_c=find_roots(lambda x:prop.axial(x,m.cruise_mps),tc,
                       prop.forward_rpms[0],prop.forward_rpms[-1],n=1800)
    if not roots_c:
        result['reason']='no measured-domain cruise trim';return result
    candidates=[]
    for h in roots_h:
        H=prop.static(h)
        for c in roots_c:
            C=prop.axial(c,m.cruise_mps)
            if C is None:continue
            Hmax=prop.static(heads[0])
            assert Hmax is not None
            if any(x['torque_nm']>m.motor_torque_max_nm+1e-9 or x['tip_mach']>m.tip_mach_limit+1e-9
                   for x in (H,C,Hmax)):continue
            if any(m.n_rotors*x['shaft_w']/m.eta_drive>m.total_electric_peak_limit_w+1e-9 for x in (H,C,Hmax)):continue
            e=m.n_rotors/m.eta_drive*(H['shaft_w']*m.hover_s+C['shaft_w']*m.distance_m/m.cruise_mps)/3600
            candidates.append((e,h,c,H,C,Hmax))
    if not candidates:
        result['reason']='motor/power/tip limit';return result
    e,h,c,H,C,Hmax=min(candidates)
    # independent SLSQP verification at the calculated exact-trim roots
    x0=np.array([h,c])
    lo=np.array([prop.static_rpm[0],prop.forward_rpms[0]],float)
    hi=np.array([prop.static_rpm[-1],prop.forward_rpms[-1]],float)
    def evaluate(z):
        speed=lo+z*(hi-lo)
        a=prop.static(float(speed[0]));b=prop.axial(float(speed[1]),m.cruise_mps)
        return a,b
    def f(z):
        a,b=evaluate(z)
        if a is None or b is None:return 1e8
        return m.n_rotors/m.eta_drive*(a['shaft_w']*m.hover_s+b['shaft_w']*m.distance_m/m.cruise_mps)/3600/100
    def eq(z):
        a,b=evaluate(z)
        if a is None or b is None:return np.array([-1e3,-1e3])
        return np.array([a['thrust_n']/th-1,b['thrust_n']/tc-1])
    opt=minimize(f,(x0-lo)/(hi-lo),method='SLSQP',bounds=[(0,1),(0,1)],
                 constraints=[{'type':'eq','fun':eq}],options={'ftol':1e-11,'maxiter':120})
    # Jac full rank is a 2 x 2 equality constraint Jacobian: local regularity check.
    dh=0.05
    Ah=(prop.static(h+dh)['thrust_n']-prop.static(h-dh)['thrust_n'])/(2*dh*th)
    Ac=(prop.axial(c+dh,m.cruise_mps)['thrust_n']-prop.axial(c-dh,m.cruise_mps)['thrust_n'])/(2*dh*tc)
    jac=np.diag([Ah,Ac])
    result.update(feasible=True,reason='measured coverage and screening constraints satisfied',
                  hover_rpm=h,cruise_rpm=c,hover_shaft_w=H['shaft_w'],cruise_shaft_w=C['shaft_w'],
                  hover_electrical_total_w=m.n_rotors*H['shaft_w']/m.eta_drive,
                  cruise_electrical_total_w=m.n_rotors*C['shaft_w']/m.eta_drive,
                  energy_wh=e,hover_torque_nm=H['torque_nm'],cruise_torque_nm=C['torque_nm'],
                  peak_reserve_rpm=heads[0],reserve_shaft_w=Hmax['shaft_w'],
                  hover_fm=th**1.5/np.sqrt(2*m.rho*np.pi*(D/2)**2)/H['shaft_w'],
                  cruise_efficiency=C['efficiency'],
                  tip_mach_max=max(H['tip_mach'],C['tip_mach'],Hmax['tip_mach']),
                  jacobian_rank=int(np.linalg.matrix_rank(jac)),jacobian_condition=float(np.linalg.cond(jac)),
                  slsqp_success=bool(opt.success),slsqp_iterations=int(opt.nit),
                  slsqp_nfev=int(opt.nfev),slsqp_message=str(opt.message),
                  slsqp_max_trim_error=float(np.max(abs(eq(opt.x)))),
                  slsqp_energy_wh=float(f(opt.x)*100),
                  root_vs_slsqp_wh=float(abs(f(opt.x)*100-e)))
    return result


def evaluate_catalog(m=Mission()):
    return [solve_catalog(name,m) for name in sorted(CATALOG)]
