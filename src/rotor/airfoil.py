import numpy as np

def surrogate_polar(alpha, p):
    """Smooth preliminary Cl/Cd model.

    This is intentionally a placeholder. Replace it with Reynolds-number-aware
    interpolation of selected airfoil/propeller polar data before final design.
    """
    cl_linear = p.cl_alpha * alpha
    cl = p.cl_max * np.tanh(cl_linear / p.cl_max)
    cd = p.cd_min + p.k_cd * cl**2
    return cl, cd
