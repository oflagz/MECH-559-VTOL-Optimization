import numpy as np
from rotoropt.experimental import predict_reference

def test_surrogate_static_ct_is_almost_rpm_independent():
    a=predict_reference(2283,polar_model='surrogate')
    b=predict_reference(5987,polar_model='surrogate')
    assert abs(a['CT']-b['CT'])<1e-5

def test_measured_proxy_changes_with_rpm_but_is_not_fully_covered():
    a=predict_reference(2283,polar_model='clarky_uiuc')
    b=predict_reference(5987,polar_model='clarky_uiuc')
    assert b['CT']>a['CT']
    assert a['measured_polar_coverage_weighted']==0
    assert b['measured_polar_coverage_weighted']<.9

def test_induction_corrections_are_not_noops():
    a=predict_reference(5987,polar_model='clarky_uiuc',
                        use_tip_hub_loss=True,include_swirl=True)
    b=predict_reference(5987,polar_model='clarky_uiuc',
                        use_tip_hub_loss=False,include_swirl=False)
    assert a['converged'] and b['converged']
    assert not np.isclose(a['CT'],b['CT'],rtol=1e-4)
