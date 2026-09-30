import numpy as np

def rotor_geometry(R: float, c: float, Nb: int):
    area = np.pi * R**2
    solidity = Nb * c / (np.pi * R)
    return area, solidity

def pitch_distribution(r, R, theta_tw_deg, theta75_deg):
    """theta(r)=theta75+theta_tw*(r/R-0.75), as defined in the proposal."""
    return np.deg2rad(theta75_deg + theta_tw_deg * (r / R - 0.75))
