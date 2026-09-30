from dataclasses import dataclass
import numpy as np

@dataclass(frozen=True)
class Params:
    """Rotor-model inputs.

    Mission values reproduce the provisional rotor screening case in the
    proposal. Values explicitly marked provisional must be reconciled with the
    final system-level mission/candidate hardware before design conclusions.
    """
    Nr: int = 4
    mTO: float = 5.0
    rho: float = 1.225
    g: float = 9.81
    a: float = 340.0
    delta_dl: float = 0.05
    MT: float = 0.15
    Vc: float = 20.0
    Dc: float = 7.0
    th: float = 120.0
    Lc: float = 36000.0
    kappa: float = 1.15
    Cd0_profile: float = 0.012
    eta_h: float = 0.85
    eta_c: float = 0.85

    # Provisional until a propeller/airfoil family is selected.
    theta75_deg: float = 20.0
    Mtip_max: float = 0.65
    FM_min: float = 0.60
    Pavail_h: float = 1000.0
    Pavail_c: float = 1000.0
    Qmotor_max: float = 1.5
    delta_clear: float = 0.05
    s_mount: float = 0.65

    # Preliminary smooth polar surrogate; replace with measured/interpolated data.
    cl_alpha: float = 2 * np.pi
    cl_max: float = 1.20
    cd_min: float = 0.012
    k_cd: float = 0.020
    hub_fraction: float = 0.20
    n_elements: int = 40
    bemt_relax: float = 0.35
    bemt_tol: float = 1e-5
    bemt_maxiter: int = 150
