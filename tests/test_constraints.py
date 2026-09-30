from src.rotor.parameters import Params
from src.rotor.constraints import constraint_margins

def test_constraint_names_match_proposal_g1_to_g10():
    margins = constraint_margins([0.25,0.030,300.,400.,-10.],3,Params())
    assert set(margins) == {
        "hover_thrust","cruise_thrust","hover_tip_mach","cruise_tip_mach",
        "hover_power","cruise_power","figure_of_merit","hover_torque",
        "cruise_torque","clearance",
    }
