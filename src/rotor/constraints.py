import numpy as np
from .hover import hover_model
from .bemt import hover_capability, cruise_model

def constraint_margins(x, Nb, p):
    """Return positive feasible margins for proposal constraints g1-g10."""
    R, c, omega_h, omega_c, theta_tw = x
    hover = hover_model(R, c, omega_h, Nb, p)
    hover_cap = hover_capability(R, c, omega_h, theta_tw, Nb, p)
    cruise = cruise_model(R, c, omega_c, theta_tw, Nb, p)
    cruise_el = cruise["P_sh"] / p.eta_c
    weight = p.mTO * p.g

    return {
        "hover_thrust": p.Nr * hover_cap["T"] - (1+p.MT)*(1+p.delta_dl)*weight,
        "cruise_thrust": p.Nr * cruise["T"] - p.Dc,
        "hover_tip_mach": p.Mtip_max - omega_h * R / p.a,
        "cruise_tip_mach": p.Mtip_max - np.sqrt(p.Vc**2+(omega_c*R)**2)/p.a,
        "hover_power": p.Pavail_h - p.Nr * hover["P_el"],
        "cruise_power": p.Pavail_c - p.Nr * cruise_el,
        "figure_of_merit": hover["FM"] - p.FM_min,
        "hover_torque": p.Qmotor_max - hover["Q"],
        "cruise_torque": p.Qmotor_max - cruise["Q"],
        "clearance": p.s_mount - (2*R + p.delta_clear),
    }

def scipy_constraints(Nb, p):
    names = tuple(constraint_margins([0.25,0.03,300.,400.,-10.], Nb, p))
    return [
        {"type": "ineq", "fun": lambda x, name=name: constraint_margins(x,Nb,p)[name]}
        for name in names
    ]
