from .hover import hover_model
from .bemt import cruise_model

def mission_energy_wh(x, Nb, p):
    """Rotor electrical mission energy, with theta75 as the fifth variable."""
    R, c, omega_h, omega_c, theta75 = x
    hover = hover_model(R, c, omega_h, Nb, p)
    cruise = cruise_model(R, c, omega_c, 0.0, Nb, p, theta75_deg=theta75)
    return p.Nr*(hover["P_el"]*p.th +
                 cruise["P_sh"]/p.eta_c*p.Lc/p.Vc)/3600.0
