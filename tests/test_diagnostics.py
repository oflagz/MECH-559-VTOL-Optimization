import numpy as np
from rotoropt.diagnostics import central_jac,kkt_audit
from rotoropt.model import Mission

def test_central_finite_difference():
    J=central_jac(lambda x:np.array([x[0]**2+x[1],x[0]*x[1]]),np.array([2.,3.]))
    assert np.allclose(J,[[4,1],[3,2]],atol=1e-7)

def test_kkt_returns_feasibility_and_multiplier_metadata():
    x=[.2541,.018,254.68,197.89,33.0]
    a=kkt_audit(x,3,Mission())
    assert len(a['gradient_normalized'])==5
    assert a['jacobian_rank']<=5
    assert a['kkt_stationarity_l2']>=0
