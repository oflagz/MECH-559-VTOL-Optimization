from dataclasses import asdict
import numpy as np
from scipy.optimize import least_squares
from .bemt import RotorGeometry,SectionPolar
from .fast_bemt import fast_performance
from .calibrate import training_data,testing_data

BOUNDS=([3.,-8.,.8,.006,.002,.0,0.,0.],[9.,5.,1.7,.08,.15,1.,.35,.25])
def make(z):return SectionPolar(*z)
def evaluate(z,data):
    polar=make(z);ans=[]
    for _,rpm,j,ct,cp in data:
        g=RotorGeometry(.127,.022,rpm*2*np.pi/60,-10.,2,14.38)
        q=fast_performance(g,j*rpm/60*.254,polar,n_sections=48,reference=True)
        if q.valid: ans.extend([(q.ct-ct)/.05,(q.cp-cp)/.035])
        else:ans.extend([.5,.5])
    return np.array(ans)
def calibrate_fast():
    init=np.array([6.,-4.,1.5,.03,.03,.5,.1,.1])
    f=least_squares(evaluate,init,bounds=BOUNDS,args=(training_data(),),max_nfev=500,
                   ftol=1e-9,xtol=1e-9,gtol=1e-9)
    return f,make(f.x)
def quality(p,data):
    a=[];invalid=[]
    for typ,rpm,j,ct,cp in data:
        g=RotorGeometry(.127,.022,rpm*2*np.pi/60,-10.,2,14.38)
        q=fast_performance(g,j*rpm/60*.254,p,n_sections=64,reference=True)
        if q.valid:a.append([ct,cp,q.ct,q.cp,rpm,j])
        else:invalid.append([typ,rpm,j,q.reason])
    a=np.array(a)
    return {'valid_n':len(a),'invalid_n':len(invalid),'ct_mae':float(np.mean(abs(a[:,0]-a[:,2]))),
            'cp_mae':float(np.mean(abs(a[:,1]-a[:,3]))),
            'ct_median_relative_pct':float(np.median(abs((a[:,0]-a[:,2])/np.maximum(abs(a[:,0]),.02)))*100),
            'cp_median_relative_pct':float(np.median(abs((a[:,1]-a[:,3])/np.maximum(abs(a[:,1]),.02)))*100),
            'invalid_points':invalid, 'observations':a.tolist()}
if __name__=='__main__':
    import time
    t=time.time();f,p=calibrate_fast()
    print('fit time',round(time.time()-t,2),'cost',f.cost,'params',f.x)
    print('train', {k:v for k,v in quality(p,training_data()).items() if k not in ('observations','invalid_points')})
    print('test',{k:v for k,v in quality(p,testing_data()).items() if k not in ('observations','invalid_points')})
