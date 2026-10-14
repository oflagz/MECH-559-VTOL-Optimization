import numpy as np
from .hover import hover_model
from .bemt import hover_capability, cruise_model

def constraint_margins(x, Nb, p):
    """Positive values are feasible. theta75 is optimized; twist fixed at zero."""
    R,c,omega_h,omega_c,theta75 = x
    hover = hover_model(R,c,omega_h,Nb,p)
    cap = hover_capability(R,c,omega_h,0.0,Nb,p,theta75_deg=theta75)
    cruise = cruise_model(R,c,omega_c,0.0,Nb,p,theta75_deg=theta75)
    weight = p.mTO*p.g
    return {
        "hover_thrust":p.Nr*cap["T"]-(1+p.MT)*(1+p.delta_dl)*weight,
        "cruise_thrust":p.Nr*cruise["T"]-p.Dc,
        "hover_tip_mach":p.Mtip_max-omega_h*R/p.a,
        "cruise_tip_mach":p.Mtip_max-np.hypot(p.Vc,omega_c*R)/p.a,
        "hover_power":p.Pavail_h-p.Nr*hover["P_el"],
        "cruise_power":p.Pavail_c-p.Nr*cruise["P_sh"]/p.eta_c,
        "figure_of_merit":hover["FM"]-p.FM_min,
        "hover_torque":p.Qmotor_max-hover["Q"],
        "cruise_torque":p.Qmotor_max-cruise["Q"],
        "clearance":p.s_mount-(2*R+p.delta_clear),
        "hover_bemt_converged":1.0 if cap["converged"] else -1.0,
        "cruise_bemt_converged":1.0 if cruise["converged"] else -1.0,
    }

def scipy_constraints(Nb,p):
    names=tuple(constraint_margins([0.25,0.03,300.,400.,20.],Nb,p))
    return [{"type":"ineq","fun":lambda x,name=name:
             constraint_margins(x,Nb,p)[name]} for name in names]
