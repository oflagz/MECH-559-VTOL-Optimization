"""Numerical agreement diagnostics (measurement interpolation, NOT BEMT validation)."""
from .catalog import MeasuredPropeller,CATALOG
import numpy as np
from scipy.interpolate import PchipInterpolator

def holdout_quality():
    results=[]
    for name in CATALOG:
        p=MeasuredPropeller(name)
        for mode,curves in [('static',[('all',p.static_rpm,np.array([float(p.static_ct(x)) for x in p.static_rpm]),np.array([float(p.static_cp(x)) for x in p.static_rpm]))]),
                            ('forward',[(str(rpm),z[0],z[1](z[0]),z[2](z[0])) for rpm,z in p.forward.items()])]:
            for key,x,ct,cp in curves:
                x=np.array(x);n=len(x)
                if n<6:continue
                mask=np.zeros(n,dtype=bool);mask[::2]=True;mask[-1]=True
                valid=(~mask)&(x>=min(x[mask]))&(x<=max(x[mask]))
                for what,y in [('CT',ct),('CP',cp)]:
                    pred=PchipInterpolator(x[mask],y[mask],extrapolate=False)(x[valid])
                    error=pred-y[valid]
                    mae=float(np.mean(np.abs(error)))
                    norm=mae/max(1e-6,float(np.max(y)-np.min(y)))
                    results.append(dict(propeller=name,condition=mode,run=key,quantity=what,
                                        train_count=int(sum(mask)),holdout_count=int(sum(valid)),
                                        mae_coefficient=mae,normalized_mae_pct=100*norm))
    return results
