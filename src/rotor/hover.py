import numpy as np
from .geometry import rotor_geometry

def hover_required_thrust(p):
    return (1.0 + p.delta_dl) * p.mTO * p.g / p.Nr

def hover_model(R, c, omega_h, Nb, p):
    """Proposal analytical hover screening model."""
    area, solidity = rotor_geometry(R, c, Nb)
    thrust = hover_required_thrust(p)
    vi = np.sqrt(thrust / (2.0 * p.rho * area))
    induced_power = thrust * vi
    profile_power = (
        solidity * p.Cd0_profile / 8.0
        * p.rho * area * (omega_h * R) ** 3
    )
    shaft_power = p.kappa * induced_power + profile_power
    return {
        "A": area,
        "sigma": solidity,
        "T_req": thrust,
        "vi": vi,
        "Pi": induced_power,
        "P0": profile_power,
        "P_sh": shaft_power,
        "P_el": shaft_power / p.eta_h,
        "FM": induced_power / shaft_power,
        "Q": shaft_power / omega_h,
    }
