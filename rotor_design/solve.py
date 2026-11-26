"""Constrained SLSQP optimize continuous radius, mean chord and twist with thrust-trim RPMs."""
from dataclasses import dataclass,replace
from functools import lru_cache
import numpy as np
from scipy.optimize import brentq,minimize
from .bemt import RotorGeometry,SectionPolar
from .fast_bemt import fast_performance as performance
from .mission import Mission,evaluate_design,takeoff_mass,cruise_drag

BOUNDS=[(.185,.297),(.009,.034),(-60.,-12.)]                 
SPEED_BOUNDS_H=(180.,730.)
SPEED_BOUNDS_C=(200.,780.)

def optimal_trim(geo,blades,theta75_deg,polar,m=Mission(),n=32):
    R,c,tw=geo;mass=takeoff_mass(R,c,blades,m)
    targets=((1+m.download_fraction)*mass*9.81/4,cruise_drag(mass,m)/4)
    speeds=[]
    for target,V,(lo,hi) in zip(targets,(0.,m.cruise_mps),(SPEED_BOUNDS_H,SPEED_BOUNDS_C)):
                                                                                       
        nodes=np.linspace(lo,hi,16)
        vv=[]
        for omega in nodes:
            g=RotorGeometry(R,c,omega,tw,blades,theta75_deg)
            p=performance(g,V,polar,n_sections=n)
            vv.append(p.thrust-target if p.valid else np.nan)
        roots=[]
        for k in range(len(nodes)-1):
            if not (np.isfinite(vv[k]) and np.isfinite(vv[k+1])):continue
            if vv[k]==0:roots.append(float(nodes[k]))
            if vv[k]*vv[k+1]<0:
                def fn(omega):
                    p=performance(RotorGeometry(R,c,omega,tw,blades,theta75_deg),V,polar,n_sections=n)
                    if not p.valid:raise ValueError('intermediate invalid')
                    return p.thrust-target
                try: roots.append(float(brentq(fn,float(nodes[k]),float(nodes[k+1]),xtol=3e-5)))
                except ValueError:continue
        if not roots:return None
        speeds.append(min(roots))
    x=[R,c,speeds[0],speeds[1],tw]
    q=evaluate_design(x,blades,theta75_deg,m,polar,n)
    return q if q['valid'] else None


def optimize_one(blades,theta75_deg,polar,m=Mission(),starts=8,n=32,seed=3401,
                 enforce_empirical_domain=True,enforce_battery=False,maxiter=75):
    rng=np.random.default_rng(seed+blades)
    scales=np.array([b[1]-b[0] for b in BOUNDS]);low=np.array([b[0] for b in BOUNDS])
    def decode(z):return low+scales*np.array(z)
    @lru_cache(maxsize=3200)
    def trial(*vals):
        geo=decode(vals)
                                                                                            
        return optimal_trim(geo,blades,theta75_deg,polar,m,n)
    def get(z):return trial(*tuple(np.round(z,9)))
    def fun(z):
        q=get(z)
        return q['energy_wh']/400 if q else 1000.
    def constraints(z):
        q=get(z)
        if q is None:return np.full((13 if enforce_battery else 12)+(3 if enforce_empirical_domain else 0),-5.)
        g=q['g'] if enforce_battery else q['g'][:-1]
        if enforce_empirical_domain:
            j=m.cruise_mps/((q['design'][3]/(2*np.pi))*2*q['design'][0])
            gh=max(abs(q['hover'].min_alpha_deg),abs(q['hover'].max_alpha_deg))/25.-1
            gc=max(abs(q['cruise'].min_alpha_deg),abs(q['cruise'].max_alpha_deg))/25.-1
            g=np.r_[g,gh,gc,j/.96-1]
        return -g
    starts_set=[[.45,.5,.5],[.7,.3,.25],[.2,.65,.65]]
    for _ in range(max(0,starts-len(starts_set))):starts_set.append(rng.uniform(.10,.90,3))
    attempts=[];candidates=[]
    for k,z0 in enumerate(starts_set):
                                                                             
        opt=minimize(fun,np.array(z0,dtype=float),bounds=[(0,1)]*3,method='SLSQP',
           constraints=[{'type':'ineq','fun':constraints}],
           options={'maxiter':maxiter,'ftol':3e-7,'disp':False})
        q=get(opt.x)
        feasible=bool(q is not None and np.min(constraints(opt.x))>=-1e-4 and np.max(np.abs(q['eq']))<=1e-4)
        attempts.append({'blade_count':blades,'theta75_deg':theta75_deg,'start':k,
                         'x0':list(map(float,decode(z0))),'x_final':list(map(float,decode(opt.x))),
                         'success':bool(opt.success),'feasible':feasible,'nfev':int(opt.nfev),
                         'iterations':int(opt.nit),'message':opt.message,
                         'energy_wh':float(q['energy_wh']) if q else None,
                         'max_g':float(-np.min(constraints(opt.x))) if q else None})
        if feasible:candidates.append(q)
    best=min(candidates,key=lambda x:x['energy_wh']) if candidates else None
    return best,attempts
