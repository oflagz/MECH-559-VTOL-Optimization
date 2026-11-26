"""Physics, mission, fitting, constraints and repeatability regression tests."""
import json
from dataclasses import replace
from pathlib import Path
import numpy as np
import pytest
from rotor_design.bemt import RotorGeometry,performance
from rotor_design.fast_bemt import fast_performance
from rotor_design.calibration_check import calibrate_fast,quality
from rotor_design.calibrate import testing_data as uiuc_test_data
from rotor_design.mission import Mission,takeoff_mass,cruise_drag,evaluate_design
from rotor_design.solve import optimal_trim

ROOT=Path(__file__).resolve().parents[1]
DATA=json.loads((ROOT/'results/tables/study_summary.json').read_text())

@pytest.fixture(scope='module')
def polar():return calibrate_fast()[1]

def test_proposal_mission_not_replaced():
    m=Mission()
    assert m.distance_m==40000 and m.hover_seconds==300 and m.cruise_seconds==1200
    assert m.cruise_mps==pytest.approx(100/3)
    assert m.payload_kg==5 and m.rotors==4

def test_fast_solver_matches_independent_reference(polar):
    g=RotorGeometry(.24,.018,440,-21,2,20)
    a=performance(g,0.,polar,n_sections=48)
    b=fast_performance(g,0.,polar,n_sections=48)
    assert a.valid and b.valid
    assert abs(a.thrust-b.thrust)<1e-7
    assert abs(a.torque-b.torque)<1e-7

def test_mass_ledger_feedback():
    m=Mission()
    assert takeoff_mass(.27,.02,3,m)>takeoff_mass(.23,.02,3,m)
    assert takeoff_mass(.23,.03,3,m)>takeoff_mass(.23,.02,3,m)
    assert cruise_drag(10.,m)>cruise_drag(8.,m)
    assert takeoff_mass(.23,.02,3,m)>m.payload_kg+m.propulsion_mass_kg+m.wing_mass_kg

def test_heldout_data_not_presented_as_complete_validation(polar):
    v=quality(polar,uiuc_test_data())
    assert v['valid_n']+v['invalid_n']==len(uiuc_test_data())
    assert v['invalid_n']>0                               
    assert 0<v['ct_mae']<.03
    assert 0<v['cp_mae']<.03

def test_all_five_original_design_variables_present():
    best=DATA['optimum'];assert len(best['design'])==5
    assert set(q['blades'] for q in DATA['blade_optima'])=={2,3,4}

def test_selected_design_exact_trim_and_numerical_limits(polar):
    best=DATA['optimum'];q=evaluate_design(best['design'],best['blade_count'],20.,polar=polar,sections=64)
    assert q['valid']
    assert max(abs(q['eq']))<1e-5
    assert np.max(q['g'])<1e-5
    assert abs(q['energy_wh']-best['energy_wh'])<.01
    assert q['hover'].fm>=.58-1e-5
    assert q['energy_wh']<Mission().usable_battery_wh

def test_no_clipping_negative_local_cruise_loading(polar):
    best=DATA['optimum']
    q=evaluate_design(best['design'],best['blade_count'],20.,polar=polar,sections=64)
    assert q['cruise'].negative_load_fraction>0.10
    assert q['cruise'].thrust>0

def test_screening_does_not_imply_existing_motor_voltage_feasible():
    best=DATA['optimum']
    Q=best['hover_torque_nm'];omega=best['design'][2];kt=.0957
    Rm=.069*(kt/.028)**2
    i0=1.1*.028/kt
    voltage=kt*omega+(Q/kt+i0)*Rm
    assert voltage>42.                                            

def test_optimization_continuity_under_refinement(polar):
    b=DATA['optimum'];x=b['design'];B=b['blade_count']
    q64=evaluate_design(x,B,20,polar=polar,sections=64)
    q128=evaluate_design(x,B,20,polar=polar,sections=128)
    assert q64['valid'] and q128['valid']
    assert abs(q64['energy_wh']-q128['energy_wh'])/q128['energy_wh']<.005

def test_decomposition_retrims_both_speed_variables(polar):
    b=DATA['optimum'];x=b['design']
    q=optimal_trim([x[0],x[1],x[4]],b['blade_count'],20.,polar,n=64)
    assert q is not None and np.max(abs(q['eq']))<1e-5
    assert abs(q['energy_wh']-b['energy_wh'])<.03
