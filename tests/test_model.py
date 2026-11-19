import numpy as np
import pytest
from rotor.catalog import CATALOG,MeasuredPropeller
from rotor.optimization import Mission,solve_catalog,evaluate_catalog
from rotor.physics import induction_from_annulus

def test_zero_extrapolation():
    p=MeasuredPropeller('18x12')
    assert p.static(100.) is None
    assert p.axial(1000.,20.) is None
    assert p.axial(3200.,33.3333) is None

def test_coefficients_units_and_efficiency():
    p=MeasuredPropeller('18x12')
    a=p.static(3500)
    assert a['shaft_w']>0 and a['thrust_n']>0
    assert abs(a['torque_nm']*3500*2*np.pi/60-a['shaft_w'])<1e-8
    c=p.axial(4017.,20)
    assert c is not None and 0<c['efficiency']<1
    assert abs(c['thrust_n']*20/c['shaft_w']-c['efficiency'])<1e-8

def test_catalog_result_physical_limits():
    m=Mission()
    results=evaluate_catalog(m)
    best=[r for r in results if r['feasible']]
    assert len(best)>0
    for r in best:
        assert r['slsqp_success']
        assert r['jacobian_rank']==2
        assert r['slsqp_max_trim_error']<1e-6
        assert r['root_vs_slsqp_wh']<1e-5
        assert r['hover_fm']<=1+1e-8
        assert r['tip_mach_max']<=m.tip_mach_limit
        assert r['hover_torque_nm']<=m.motor_torque_max_nm
        assert r['required_hover_n']>m.payload_mass_kg*m.gravity/4

def test_official_mission_not_claimed_valid():
    m=Mission(cruise_mps=40000/1200)
    assert not any(r['feasible'] for r in evaluate_catalog(m))

def test_negative_axial_load_signed_instead_of_clipped():
    r=np.array([.18,.20]);F=np.ones(2)
    vi,swirl=induction_from_annulus(np.array([-2.,2.]),np.array([-.01,.02]),r,20.,F,1.225,300.)
    assert vi[0]<0<vi[1]
    assert swirl[0]<0<swirl[1]

def test_invalid_hover_negative_load_rejected():
    with pytest.raises(ValueError):
        induction_from_annulus(np.array([-1.]),np.array([0.]),np.array([.2]),0.,np.array([1.]),1.225,300.)

def test_mass_increase_increases_hover_thrust():
    x=solve_catalog('18x12',Mission(base_aircraft_mass_kg=9.0))
    y=solve_catalog('18x12',Mission(base_aircraft_mass_kg=10.0))
    assert y['required_hover_n']>x['required_hover_n']
