"""Numerical/reference-data sanity checks; not accuracy certification."""
import numpy as np
from rotoropt.experimental import load_geometry,predict_reference,RADIUS,NBLADES

def test_measured_propeller_geometry_is_complete():
    g=load_geometry()
    assert len(g)==18
    assert np.all(np.diff(g['r_over_R'])>0)
    assert np.all(g['c_over_R']>0)
    assert np.all(g['pitch_deg']>0)
    assert NBLADES==2
    assert np.isclose(RADIUS,.127)

def test_bemt_reference_static_runs_and_outputs_finite_coefficients():
    out=predict_reference(3029,0.,120)
    assert out['converged']
    assert out['CT']>0 and out['CP']>0
    assert np.isfinite(out['eta'])
    assert out['re_min']>1000

def test_bemt_reference_forward_runs_without_fitting_to_measurements():
    a=predict_reference(3008,.383,120)
    assert a['converged']
    assert a['CT']>0 and a['CP']>0
    assert 0<a['eta']<1.1
