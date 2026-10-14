from src.rotor.parameters import Params
from src.rotor.constraints import constraint_margins

def test_constraint_names_match_proposal_and_solver_guards():
    margins=constraint_margins([0.25,0.030,300.,400.,20.],3,Params())
    assert set(margins)=={
        "hover_thrust","cruise_thrust","hover_tip_mach","cruise_tip_mach",
        "hover_power","cruise_power","figure_of_merit","hover_torque",
        "cruise_torque","clearance","hover_bemt_converged",
        "cruise_bemt_converged"}
