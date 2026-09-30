from .hover import hover_model
from .bemt import cruise_model

def mission_energy_wh(x, Nb, p):
    """Rotor electrical mission energy, matching the proposal rotor objective."""
    R, c, omega_h, omega_c, theta_tw = x
    hover = hover_model(R, c, omega_h, Nb, p)
    cruise = cruise_model(R, c, omega_c, theta_tw, Nb, p)
    cruise_electrical = cruise["P_sh"] / p.eta_c
    cruise_time = p.Lc / p.Vc
    energy_j = p.Nr * (
        hover["P_el"] * p.th + cruise_electrical * cruise_time
    )
    return energy_j / 3600.0
