"""Constrained mission-energy optimization with physically matched thrust.

Decision vector is [radius, chord, hover omega, cruise omega, theta75].
Two equality constraints impose hover/cruise steady-state thrust trim.
An additional hover RPM headroom check tests the proposal's thrust margin.
Power/torque are always evaluated at the SAME thrust and RPM.
"""
from functools import lru_cache
from dataclasses import asdict
import numpy as np
from scipy.optimize import minimize, brentq
from .model import Mission,predict

NBLADES=(2,3,4)
# Bounds are exploratory geometric/RPM/pitch ranges, NOT manufacturer certifications.
LOW=np.array([.18,.018,170.,170.,10.])
HIGH=np.array([.30,.052,650.,750.,33.])
NAMES=['radius_m','chord_m','omega_hover_rad_s','omega_cruise_rad_s','theta75_deg']

@lru_cache(maxsize=12000)
def _cache(key, n_blades, mission):
    R,c,oh,oc,theta=key
    h=predict(R,c,oh,theta,n_blades,0.,mission)
    cruise=predict(R,c,oc,theta,n_blades,mission.cruise_speed_mps,mission)
    reserve=predict(R,c,1.1*oh,theta,n_blades,0.,mission)
    return h,cruise,reserve

def analyze(x, n_blades, mission=Mission()):
    x=np.asarray(x,float)
    h,c,r=_cache(tuple(float(v) for v in x),n_blades,mission)
    R,ch,oh,oc,theta=x
    m=mission
    energy_wh=m.n_rotors*(h['electrical_power_w']*m.hover_seconds+
                            c['electrical_power_w']*m.cruise_seconds)/3600
    # Scaled constraint margins; >= 0 means feasible.
    margins={
      'hover_power':1-m.n_rotors*h['electrical_power_w']/m.total_max_electrical_w,
      'cruise_power':1-m.n_rotors*c['electrical_power_w']/m.total_max_electrical_w,
      'motor_hover_torque':1-h['torque_nm']/m.motor_max_torque_nm,
      'motor_cruise_torque':1-c['torque_nm']/m.motor_max_torque_nm,
      'hover_tip_mach':1-h['tip_mach']/m.max_tip_mach,
      'cruise_tip_mach':1-c['tip_mach']/m.max_tip_mach,
      'mount_clearance':(m.mount_spacing_m-2*R-m.clearance_m)/m.mount_spacing_m,
      'hover_figure_of_merit':h['figure_of_merit']-m.min_hover_fm,
      'hover_thrust_margin':r['thrust_n']/m.hover_thrust_n-(1+m.thrust_margin),
      'reserve_power':1-m.n_rotors*r['electrical_power_w']/m.total_max_electrical_w,
      'reserve_torque':1-r['torque_nm']/m.motor_max_torque_nm,
      'reserve_tip_mach':1-r['tip_mach']/m.max_tip_mach,
      'bemt_convergence':1. if all(z['converged'] for z in (h,c,r)) else -1.,
      'cruise_positive_power':c['shaft_power_w']/100.,
    }
    trim=np.array([h['thrust_n']/m.hover_thrust_n-1,
                   c['thrust_n']/m.cruise_thrust_n-1])
    return {'energy_wh':float(energy_wh),'hover':h,'cruise':c,'reserve':r,
            'trim_residual':trim,'margins':margins,
            'geometry':dict(zip(NAMES,[float(v) for v in x])),
            'n_blades':n_blades,'disk_area_total_m2':m.n_rotors*np.pi*R*R,
            'disk_loading_n_m2':m.mass_kg*m.gravity/(m.n_rotors*np.pi*R*R)}

def _initial_speed(R,c,theta,n_blades,mission,axial_speed,target):
    # Scan to bracket first positive-thrust trim crossing. This is not an
    # optimization shortcut: it generates candidates satisfying equalities.
    def residual(omega):
        z=predict(R,c,omega,theta,n_blades,axial_speed,mission)
        return z['thrust_n']-target if z['converged'] else np.nan
    speeds=np.array([LOW[2], HIGH[3]])
    last_o,last_f=speeds[0],residual(speeds[0])
    for o in speeds[1:]:
        f=residual(o)
        if np.isfinite(f) and np.isfinite(last_f) and f*last_f<=0:
            try:
                return float(brentq(residual,last_o,o,xtol=1e-4))
            except (ValueError,RuntimeError):
                pass
        last_o,last_f=o,f
    return None

def initial_guesses(mission,n_blades):
    for R,c,theta in [(.25,.030,20),(.27,.036,22),(.22,.041,24),(.29,.025,18),(.20,.046,28)]:
        oh=_initial_speed(R,c,theta,n_blades,mission,0.,mission.hover_thrust_n)
        oc=_initial_speed(R,c,theta,n_blades,mission,mission.cruise_speed_mps,
                          mission.cruise_thrust_n)
        if oh is not None and oc is not None and LOW[2]<=oh<=HIGH[2] and LOW[3]<=oc<=HIGH[3]:
            yield np.array([R,c,oh,oc,theta])

def solve(mission=Mission(),n_blades=3,max_starts=3):
    starts=list(initial_guesses(mission,n_blades))[:max_starts]
    if not starts:
        return None,[],[]
    results=[]
    for x0 in starts:
        def full(y):return analyze(LOW+y*(HIGH-LOW),n_blades,mission)
        def objective(y):return full(y)['energy_wh']/100
        def eq(y):return full(y)['trim_residual']
        def ineq(y):return np.array(list(full(y)['margins'].values()))
        res=minimize(objective,(x0-LOW)/(HIGH-LOW),method='SLSQP',
                     bounds=[(0.,1.)]*5,
                     constraints=[{'type':'eq','fun':eq},{'type':'ineq','fun':ineq}],
                     options={'maxiter':110,'ftol':1e-7,'disp':False})
        x=LOW+res.x*(HIGH-LOW)
        a=analyze(x,n_blades,mission)
        feasible=(np.max(abs(a['trim_residual']))<1e-3 and
                  min(a['margins'].values())>=-1e-3 and
                  np.all(x>=LOW-1e-6) and np.all(x<=HIGH+1e-6))
        results.append({'success':bool(res.success),'feasible':bool(feasible),
                        'message':str(res.message),'iterations':int(res.nit),
                        'x':x.tolist(),'analysis':a})
    feasible=[r for r in results if r['feasible']]
    return (min(feasible,key=lambda r:r['analysis']['energy_wh'])
            if feasible else None),results,starts
