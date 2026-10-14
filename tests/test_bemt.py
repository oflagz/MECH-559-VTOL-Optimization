import numpy as np
from src.rotor.parameters import Params
from src.rotor.bemt import bemt_axial

def test_bemt_returns_finite_outputs():
    result = bemt_axial(0.25,0.03,300.,-10.,3,0.,Params())
    for key in ("T","Q","P_sh","vi","Re_min","Re_max"):
        assert np.isfinite(result[key])
    assert isinstance(result["converged"],bool)
    assert result["iterations"] <= Params().bemt_maxiter
