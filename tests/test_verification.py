import numpy as np
from dataclasses import replace
from rotor.verification import holdout_quality
from rotor.catalog import MeasuredPropeller
from rotor.optimization import Mission,solve_catalog
from rotor.physics import induction_from_annulus,bemt_diagnostic

def test_measured_interpolation_holdout_well_behaved():
    data=holdout_quality()
    assert len(data)>10
    assert all(item['mae_coefficient']>=0 for item in data)
    dynamic=[item for item in data if item['condition']=='forward']
    assert np.median([x['normalized_mae_pct'] for x in dynamic])<5

def test_nominal_propeller_mass_feedback():
    cases=[]
    for mass in (8.5,9.5,10.5):
        q=solve_catalog('18x12',Mission(base_aircraft_mass_kg=mass))
        assert q['feasible']
        cases.append(q)
    assert cases[0]['energy_wh']<cases[1]['energy_wh']<cases[2]['energy_wh']
    assert cases[0]['hover_rpm']<cases[1]['hover_rpm']<cases[2]['hover_rpm']

def test_cruise_drag_sensitivity():
    cases=[]
    for drag in (5,7,9):
        q=solve_catalog('18x12',Mission(cruise_drag_n=drag))
        assert q['feasible']
        cases.append(q['energy_wh'])
    assert cases[0]<cases[1]<cases[2]

def test_cruise_efficiency_limited_by_one():
    p=MeasuredPropeller('18x12')
    for rpm in (3032,4017,4991):
        for J in np.linspace(.48,.60,10):
            V=J*(rpm/60)*p.prop.diameter_m
            q=p.axial(rpm,V)
            if q:
                assert 0<q['efficiency']<1

def test_signed_negative_blade_loading_diagnostic():
    def surrogate(a,re):
        cl=1.18*np.tanh(2*np.pi*a/1.18)
        cd=.015+.027*cl*cl+.12*np.maximum(abs(a)-.23,0)**2
        return cl,cd
    old=bemt_diagnostic(.254102,.018,197.8851,33,3,20,surrogate)
    assert not old.converged
    assert old.reason!='converged'

def test_no_unverified_advance_ratio_extrapolation():
    prop=MeasuredPropeller('18x12')
    at=prop.axial(3500.,20)
    assert at is not None
    assert prop.axial(3500.,33.333333) is None
