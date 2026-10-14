import numpy as np
from rotoropt.model import Mission
from rotoropt.optimize import solve,analyze,LOW,HIGH

def test_three_blade_constrained_solution():
    best,runs,starts=solve(Mission(),3,max_starts=2)
    assert best is not None
    assert all(r['feasible'] for r in runs)
    x=np.array(best['x'])
    a=best['analysis']
    assert np.all(x>=LOW-1e-5) and np.all(x<=HIGH+1e-5)
    assert max(abs(a['trim_residual']))<1e-3
    assert min(a['margins'].values())>-1e-3
    assert a['energy_wh']>0
    assert len(starts)>0
