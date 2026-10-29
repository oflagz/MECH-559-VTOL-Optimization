"""Two separate aerodynamic options: legacy illustrative polar, and real
UIUC Clark-Y(B) section test data used as an APC section *proxy*.

For measured data, C_l is interpolated from the independent CLARKYB.LFT
measurement at the exact angles appearing in CLARKYB.DRG. Both raw lift
and drag inputs originate from UIUC LSAT Vol.3. Between Reynolds groups:
log-linear interpolation; beyond measured Re / alpha: constant endpoint
extension *with reported coverage flags*, never silent extrapolation.
This is not a polar measurement of the proprietary APC E63/Clark-Y blend.
"""
from functools import lru_cache
from pathlib import Path
import numpy as np

DATA = Path(__file__).resolve().parents[1]/'data'/'airfoil_polars'/'clarky_uiuc_experimental.csv'

@lru_cache(maxsize=1)
def _table():
    table=np.genfromtxt(DATA,delimiter=',',names=True,dtype=None,encoding='utf8')
    grid={}
    for re in np.unique(table['reynolds']):
        rows=table[table['reynolds']==re]
        order=np.argsort(rows['alpha_deg'])
        a=rows['alpha_deg'][order]
        assert len(a)>=8 and np.all(np.diff(a)>0)
        grid[int(re)]=(a,rows['cl_interpolated'][order],rows['cd_measured'][order])
    return grid

def coefficients(alpha_rad,reynolds,model='surrogate'):
    """Return C_l, C_d, and flags showing whether the source data cover inputs.

    alpha_rad and reynolds may be scalar/arrays of any broadcastable shape.
    """
    alpha, re=np.broadcast_arrays(np.asarray(alpha_rad,dtype=float),
                                  np.asarray(reynolds,dtype=float))
    if model=='surrogate':
        cl=1.18*np.tanh(2*np.pi*alpha/1.18)
        cd=.015+.027*cl**2+.12*np.maximum(abs(alpha)-.23,0)**2
        flags=np.zeros_like(alpha,dtype=bool)  # no measured-data validity for synthetic model
        return cl,cd,flags,flags
    if model!='clarky_uiuc':
        raise ValueError('Unknown polar model: '+str(model))
    grid=_table()
    keys=np.array(sorted(grid))
    low=np.searchsorted(keys,re,side='right')-1
    low=np.clip(low,0,len(keys)-1)
    high=np.clip(low+1,0,len(keys)-1)
    lnre=np.log(np.maximum(re,1.))
    denom=np.log(keys[high])-np.log(keys[low])
    factor=np.divide(lnre-np.log(keys[low]),denom,
                     out=np.zeros_like(re,dtype=float),where=denom>0)
    factor=np.clip(factor,0,1)
    degrees=np.rad2deg(alpha)
    shape=alpha.shape
    cl0=np.empty(shape);cd0=np.empty(shape);cl1=np.empty(shape);cd1=np.empty(shape)
    avalid=np.ones(shape,dtype=bool)
    # Per-bucket interpolation; alpha clamps in out-of-data regions, flagged.
    for k in np.unique(np.concatenate([low.ravel(),high.ravel()])):
        a,c,d=grid[int(keys[k])]
        for sel,outcl,outcd in ((low==k,cl0,cd0),(high==k,cl1,cd1)):
            if np.any(sel):
                outcl[sel]=np.interp(degrees[sel],a,c)
                outcd[sel]=np.interp(degrees[sel],a,d)
                avalid[sel]&=(degrees[sel]>=a[0])&(degrees[sel]<=a[-1])
    return (cl0+(cl1-cl0)*factor,cd0+(cd1-cd0)*factor,
            avalid,(re>=keys[0])&(re<=keys[-1]))
