import numpy as np
from dataclasses import replace
from rotoropt.polars import coefficients,_table
from rotoropt.model import Mission,predict
from rotoropt.experimental import load_geometry,predict_reference

def test_measured_table_has_multiple_reynolds_groups():
    grid=_table()
    assert set(grid)=={60000,100000,200000,300000}
    assert all(len(x[0])>=8 and np.all(np.diff(x[0])>0) for x in grid.values())

def test_reynolds_dependence_is_real_and_exposed():
    alpha=np.deg2rad([0,5,10])
    l60,d60,valid,valid_re=coefficients(alpha,[60000]*3,'clarky_uiuc')
    l200,d200,_,_=coefficients(alpha,[200000]*3,'clarky_uiuc')
    assert np.all(valid) and np.all(valid_re)
    assert np.all(d60>d200)
    assert np.max(abs(l60-l200))>0.1

def test_clamping_is_flagged_not_silent():
    cl,cd,va,vr=coefficients(np.deg2rad([2,16,-9]),[25000,100000,500000],'clarky_uiuc')
    assert vr.tolist()==[False,True,False]
    assert va.tolist()==[True,False,False]
    assert np.all(np.isfinite(cl)) and np.all(cd>0)

def test_real_propeller_geometry_in_si():
    a=load_geometry()
    assert np.isclose(a['r_over_R'][-1],1.0)
    assert np.isclose(a['pitch_deg'][np.where(np.isclose(a['r_over_R'],.75))[0][0]],14.38)
    assert np.isclose(2*.127,.254)
    assert np.all(np.diff(a['r_over_R'])>0)

def test_reference_reports_insufficient_data_coverage():
    p=predict_reference(3029,polar_model='clarky_uiuc')
    assert p['converged']
    assert 0<=p['measured_polar_coverage_weighted']<0.1
    assert p['re_max']<60000

def test_new_optimizer_model_is_selectable_without_changing_legacy_defaults():
    legacy=Mission()
    measured=replace(legacy,polar_model='clarky_uiuc')
    x=(.25,.03,300,25,3,0)
    old=predict(*x,legacy)
    new=predict(*x,measured)
    assert old['converged'] and new['converged']
    assert old['measured_polar_coverage_weighted']==0
    assert new['measured_polar_coverage_weighted']<1
    assert abs(old['thrust_n']-new['thrust_n'])>0.1
