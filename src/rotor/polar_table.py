"""Reynolds-aware interpolation for tabulated airfoil polars.

CSV columns: reynolds, alpha_deg, cl, cd. Each Reynolds group must contain
identical alpha grids. Out-of-range queries are clipped, not extrapolated.
"""
import csv
import numpy as np

class PolarTable:
    def __init__(self, csv_path):
        with open(csv_path, newline="") as handle:
            rows = list(csv.DictReader(handle))
        if not rows:
            raise ValueError("Empty polar table")
        data = np.array([[float(row[k]) for k in
                          ("reynolds","alpha_deg","cl","cd")] for row in rows])
        self.re = np.unique(data[:,0])
        self.alpha = np.unique(data[:,1])
        self.cl = np.empty((len(self.re),len(self.alpha)))
        self.cd = np.empty_like(self.cl)
        for i,re in enumerate(self.re):
            group = data[data[:,0]==re]
            if len(group)!=len(self.alpha) or not np.array_equal(
                np.sort(group[:,1]),self.alpha):
                raise ValueError("Every Reynolds group must share alpha grid")
            group = group[np.argsort(group[:,1])]
            self.cl[i] = group[:,2]
            self.cd[i] = group[:,3]
        if np.any(self.cd < 0) or np.any(self.re <= 0):
            raise ValueError("Invalid polar values")

    def __call__(self, alpha, reynolds, p=None):
        a, re = np.broadcast_arrays(np.rad2deg(alpha),reynolds)
        a = np.clip(a,self.alpha[0],self.alpha[-1])
        re = np.clip(re,self.re[0],self.re[-1])
        def interpolate(table):
            low = np.searchsorted(self.re,re,side="right")-1
            low = np.clip(low,0,len(self.re)-1)
            high = np.clip(low+1,0,len(self.re)-1)
            denom = np.log(self.re[high])-np.log(self.re[low])
            t = np.divide(np.log(re)-np.log(self.re[low]), denom,
                          out=np.zeros_like(re,dtype=float),
                          where=denom!=0)
            v0 = np.array([np.interp(ai,self.alpha,table[i])
                           for ai,i in zip(a.flat,low.flat)]).reshape(a.shape)
            v1 = np.array([np.interp(ai,self.alpha,table[i])
                           for ai,i in zip(a.flat,high.flat)]).reshape(a.shape)
            return v0+(v1-v0)*t
        return interpolate(self.cl),interpolate(self.cd)
