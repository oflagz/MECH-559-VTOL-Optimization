"""Annular axial/swirl BEMT prototype.

Positive thrust and torque correspond to propeller operation. Momentum closure
uses annular hover/axial flow equations without Prandtl tip/hub corrections.
Outputs include convergence status; use validated polar data for design work.
"""
import numpy as np
from .airfoil import surrogate_polar
from .geometry import pitch_distribution
from .hover import hover_required_thrust

def bemt_axial(R, c, omega, theta_tw_deg, Nb, V_axial, p, polar=None, theta75_deg=None):
    if R <= 0 or c <= 0 or omega <= 0:
        return {"T": -1e9, "P_sh": 1e12, "Q": 1e12, "vi": 0.0,
                "converged": False, "iterations": 0}
    polar = polar or surrogate_polar
    edges = np.linspace(p.hub_fraction * R, R, p.n_elements + 1)
    r = 0.5 * (edges[:-1] + edges[1:])
    dr = np.diff(edges)
    theta = pitch_distribution(r, R, theta_tw_deg, p.theta75_deg if theta75_deg is None else theta75_deg)
    vi = np.full_like(r, max(0.1, np.sqrt(hover_required_thrust(p) /
                                         (2*p.rho*np.pi*R**2))))
    swirl = np.zeros_like(r)
    converged = False
    for iteration in range(1, p.bemt_maxiter + 1):
        va = V_axial + vi
        vt = omega*r - swirl
        phi = np.arctan2(va, vt)
        speed = np.hypot(va, vt)
        alpha = theta - phi
        reynolds = p.rho * speed * c / 1.81e-5
        try:
            cl, cd = polar(alpha, reynolds, p)
        except TypeError:
            cl, cd = polar(alpha, p)
        q = 0.5*p.rho*speed**2
        dL = q*c*cl*dr*Nb
        dD = q*c*cd*dr*Nb
        dT = dL*np.cos(phi) - dD*np.sin(phi)
        dQ = (dL*np.sin(phi) + dD*np.cos(phi))*r
        # dT = 4*pi*r*rho*(V+vi)*vi*dr; solve positive induced velocity.
        loading = np.maximum(dT/dr, 0.0)
        vi_target = 0.5*(-V_axial + np.sqrt(V_axial**2 +
                         loading/(np.pi*p.rho*r)))
        # dQ = 4*pi*r^2*rho*(V+vi)*vtheta*dr.
        swirl_target = np.maximum(dQ/dr, 0.0) / (
            4*np.pi*p.rho*r**2*np.maximum(V_axial+vi_target, 1e-5))
        swirl_target = np.minimum(swirl_target, 0.8*omega*r)
        vi_new = (1-p.bemt_relax)*vi + p.bemt_relax*vi_target
        swirl_new = (1-p.bemt_relax)*swirl + p.bemt_relax*swirl_target
        err = max(float(np.max(abs(vi_new-vi))),
                  float(np.max(abs(swirl_new-swirl))))
        vi, swirl = vi_new, swirl_new
        if err < p.bemt_tol:
            converged = True
            break
    va = V_axial+vi
    vt = omega*r-swirl
    phi = np.arctan2(va, vt)
    speed = np.hypot(va, vt)
    reynolds = p.rho*speed*c/1.81e-5
    try:
        cl, cd = polar(theta-phi, reynolds, p)
    except TypeError:
        cl, cd = polar(theta-phi, p)
    q = 0.5*p.rho*speed**2
    dL = q*c*cl*dr*Nb
    dD = q*c*cd*dr*Nb
    thrust = float(np.sum(dL*np.cos(phi)-dD*np.sin(phi)))
    torque = float(np.sum((dL*np.sin(phi)+dD*np.cos(phi))*r))
    return {"T": thrust, "P_sh": max(torque*omega, 0.0), "Q": torque,
            "vi": float(np.mean(vi)), "converged": converged,
            "iterations": iteration, "Re_min": float(np.min(reynolds)),
            "Re_max": float(np.max(reynolds))}

def hover_capability(R, c, omega_h, theta_tw, Nb, p, theta75_deg=None):
    return bemt_axial(R, c, omega_h, theta_tw, Nb, 0.0, p, theta75_deg=theta75_deg)

def cruise_model(R, c, omega_c, theta_tw, Nb, p, theta75_deg=None):
    return bemt_axial(R, c, omega_c, theta_tw, Nb, p.Vc, p, theta75_deg=theta75_deg)
