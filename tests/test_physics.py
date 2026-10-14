import math
from dataclasses import replace
import numpy as np
from rotoropt.model import Mission,predict,hover_momentum_check
from rotoropt.optimize import analyze, LOW, HIGH

def test_proposal_analytical_hover_baseline():
    a=hover_momentum_check(.25,.03,300,3)
    assert abs(a['thrust_n']-12.88)<.02
    assert abs(a['electrical_total_w']-443)<2
    assert abs(a['figure_of_merit']-.708)<.002

def test_axial_bemt_converges_and_has_plausible_signs():
    m=Mission()
    h=predict(.25,.03,230,22,3,0,m)
    c=predict(.25,.03,350,22,3,m.cruise_speed_mps,m)
    assert h['converged'] and c['converged']
    assert h['thrust_n']>0 and h['torque_nm']>0 and h['shaft_power_w']>0
    assert c['thrust_n']>0 and c['shaft_power_w']>0
    assert 0<h['figure_of_merit']<1
    assert np.isclose(h['electrical_power_w'],h['shaft_power_w']/m.eta_motor)

def test_thrust_increases_with_hover_speed_near_design_point():
    x=[predict(.25,.03,w,22,3,0)['thrust_n'] for w in [220,260,300]]
    assert x[0]<x[1]<x[2]

def test_energy_duration_and_mass_are_coupled():
    m=Mission()
    x=[.25,.03,230,350,22]
    base=analyze(x,3,m)
    longer=analyze(x,3,replace(m,cruise_distance_m=2*m.cruise_distance_m))
    assert longer['energy_wh']>base['energy_wh']
    assert np.isclose(longer['energy_wh']-base['energy_wh'],
           4*base['cruise']['electrical_power_w']*m.cruise_seconds/3600)
    assert replace(m,mass_kg=6).hover_thrust_n>m.hover_thrust_n

def test_screening_not_confused_with_payload_mission():
    m=Mission()
    assert m.cruise_seconds==1800
    assert m.mass_kg==5.0
    assert m.hover_seconds==120
