"""Independent physics-model fitting to UIUC SF 10x7, with reserved measurement groups."""
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from .bemt import SectionPolar,reference_performance,ROOT

SOURCE=ROOT/'data/uiuc'

def training_data():
    rows=[]
    st=np.genfromtxt(SOURCE/'apcsf_10x7_static.csv',delimiter=',',names=True)
    for v in st[::2]: rows.append(('static',float(v['rpm']),0.,float(v['ct']),float(v['cp'])))
    for speed in (3008,4011):
        arr=np.genfromtxt(SOURCE/f'apcsf_10x7_{speed}.csv',delimiter=',',names=True)
        for row in arr[::3]:rows.append(('dynamic',float(speed),float(row['j']),float(row['ct']),float(row['cp'])))
    return rows

def testing_data():
    rows=[]
    for speed in (5003,5006,6014):
        arr=np.genfromtxt(SOURCE/f'apcsf_10x7_{speed}.csv',delimiter=',',names=True)
        for row in arr[::2]:rows.append(('dynamic',float(speed),float(row['j']),float(row['ct']),float(row['cp'])))
    st=np.genfromtxt(SOURCE/'apcsf_10x7_static.csv',delimiter=',',names=True)
    for row in st[1::2]:rows.append(('static',float(row['rpm']),0.,float(row['ct']),float(row['cp'])))
    return rows

def params_to_polar(z):
    return SectionPolar(cl_slope=z[0],alpha_zero_deg=z[1],cl_max=z[2],cd0=z[3],drag_k=z[4],stall_drag=z[5],re_slope=.08,re_cd=.12)

def residuals(z,data,ns=32):
    polar=params_to_polar(z);res=[]
    for _,rpm,j,ct,cp in data:
        p=reference_performance(rpm,j,polar,ns)
        if not p.valid or p.ct<-.1 or p.cp<-.1:
            res.extend([1.,1.]);continue
        res.extend([(p.ct-ct)/.05,(p.cp-cp)/.03])
    return res

def calibrate():
    bounds=([3.5,-5.,.7,.008,.002,0],[8.,6.,1.9,.07,.11,2.])
    init=np.array([5.7,-1.,1.3,.024,.035,.8])
    fit=least_squares(residuals,init,bounds=bounds,args=(training_data(),),
                      max_nfev=90,ftol=1e-5,xtol=1e-5,gtol=1e-5)
    return fit,params_to_polar(fit.x)

if __name__=='__main__':
    import time
    t=time.time();o,p=calibrate()
    print('calibration',o.success,o.nfev,o.cost,o.x,'seconds',time.time()-t)
    for label,data in [('train',training_data()),('holdout',testing_data())]:
        errors=[];failed=0
        for _,rpm,j,ct,cp in data:
            v=reference_performance(rpm,j,p)
            if v.valid:errors.append((abs(v.ct-ct),abs(v.cp-cp)))
            else:failed+=1
        print(label,'n',len(errors),'invalid',failed,'MAE',np.mean(errors,axis=0) if errors else None)
