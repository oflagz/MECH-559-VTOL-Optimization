import numpy as np
from .airfoil import surrogate_polar
from .geometry import pitch_distribution
from .hover import hover_required_thrust

def bemt_axial(R, c, omega, theta_tw_deg, Nb, V_axial, p):
    """Preliminary axial blade-element/momentum model.

    Uses uniform induced velocity and no swirl. This makes the optimization
    executable now, but it is not the final BEMT promised by the proposal.
    """
    if R <= 0 or c <= 0 or omega <= 0:
        return {"T": -1e9, "P_sh": 1e12, "Q": 1e12, "vi": 0.0}

    edges = np.linspace(p.hub_fraction * R, R, p.n_elements + 1)
    r = 0.5 * (edges[:-1] + edges[1:])
    dr = np.diff(edges)
    theta = pitch_distribution(r, R, theta_tw_deg, p.theta75_deg)
    area = np.pi * R**2
    vi = max(0.1, np.sqrt(hover_required_thrust(p) / (2 * p.rho * area)))

    for _ in range(p.bemt_maxiter):
        va = V_axial + vi
        vt = omega * r
        phi = np.arctan2(va, vt)
        wrel = np.sqrt(va**2 + vt**2)
        cl, cd = surrogate_polar(theta - phi, p)
        q = 0.5 * p.rho * wrel**2
        dL = q * c * cl * dr * Nb
        dD = q * c * cd * dr * Nb
        thrust = np.sum(dL * np.cos(phi) - dD * np.sin(phi))
        torque = np.sum((dL * np.sin(phi) + dD * np.cos(phi)) * r)

        if thrust <= 0:
            target = 0.0
        else:
            target = 0.5 * (
                -V_axial + np.sqrt(V_axial**2 + 2 * thrust / (p.rho * area))
            )
        updated = (1 - p.bemt_relax) * vi + p.bemt_relax * target
        if abs(updated - vi) < p.bemt_tol:
            vi = updated
            break
        vi = updated

    return {"T": thrust, "P_sh": max(torque * omega, 0.0), "Q": torque, "vi": vi}

def hover_capability(R, c, omega_h, theta_tw, Nb, p):
    return bemt_axial(R, c, omega_h, theta_tw, Nb, 0.0, p)

def cruise_model(R, c, omega_c, theta_tw, Nb, p):
    return bemt_axial(R, c, omega_c, theta_tw, Nb, p.Vc, p)
