"""Finite-difference active-set and approximate KKT diagnostics (L8-L9).

Assumes objective f=E/100, equality trim h=0, inequality margins g>=0,
normalized x: y=(x-LOW)/(HIGH-LOW), 0<=y<=1.
KKT: grad f + J_h.T lambda - J_g.T mu = 0, mu>=0.

Finite-difference KKT residuals are diagnostics, not proofs of global optima.
"""
import numpy as np
from scipy.optimize import lsq_linear
from .optimize import LOW, HIGH, NAMES, analyze
from .model import Mission

def central_jac(func,y, step=3e-5):
    y=np.asarray(y,float)
    result=[]
    for i in range(len(y)):
        eps=np.zeros_like(y);eps[i]=step
        result.append((np.asarray(func(y+eps),float)-
                       np.asarray(func(y-eps),float))/(2*step))
    return np.stack(result,axis=-1)

def kkt_audit(x,nb=3,m=Mission(),active_tol=2e-3):
    y=(np.asarray(x,float)-LOW)/(HIGH-LOW)
    def values(y_):
        z=analyze(LOW+y_*(HIGH-LOW),nb,m)
        return z['energy_wh']/100.,np.asarray(z['trim_residual']),z['margins']
    f,h,gdict=values(y)
    keys=list(gdict)
    g=np.array([gdict[k] for k in keys])
    df=central_jac(lambda z:np.array([values(z)[0]]),y).ravel()
    Jh=central_jac(lambda z:values(z)[1],y)
    Jg=central_jac(lambda z:np.array(list(values(z)[2].values())),y)
    active_idx=np.flatnonzero(g<=active_tol)
    active_names=[keys[j] for j in active_idx]
    Jactive=Jg[active_idx] if len(active_idx) else np.zeros((0,len(y)))
    for j in range(len(y)):
        if y[j] <= active_tol:
            vec=np.eye(len(y))[j]
            Jactive=np.vstack([Jactive,vec])
            active_names.append(NAMES[j]+':lower')
        if 1-y[j] <= active_tol:
            vec=-np.eye(len(y))[j]
            Jactive=np.vstack([Jactive,vec])
            active_names.append(NAMES[j]+':upper')
    A=np.column_stack([Jh.T,-Jactive.T])
    lb=np.r_[np.full(len(h),-np.inf),np.zeros(len(active_names))]
    ub=np.full(A.shape[1],np.inf)
    solution=lsq_linear(A,-df,bounds=(lb,ub),tol=1e-12,max_iter=1000)
    stat=df+A@solution.x
    Jtot=np.vstack([Jh,Jactive])
    rank=int(np.linalg.matrix_rank(Jtot,tol=1e-6))
    return dict(
        energy_wh=float(f*100),
        eq_residual_max=float(np.max(abs(h))),
        inequality_min_margin=float(np.min(g)),
        active_constraints=active_names,
        equality_multipliers=solution.x[:len(h)].tolist(),
        active_inequality_multipliers=dict(zip(active_names,solution.x[len(h):].tolist())),
        kkt_stationarity_l2=float(np.linalg.norm(stat)),
        kkt_stationarity_inf=float(np.max(abs(stat))),
        gradient_normalized=df.tolist(),
        jacobian_rows=int(len(Jtot)),jacobian_rank=rank,
        licq_full_row_rank=bool(rank==len(Jtot)),
        condition_singular_values=np.linalg.svd(Jtot,compute_uv=False).tolist(),
        normalized_x=y.tolist(),
        status='first-order numerical diagnostic only; no guarantee of global optimality')
