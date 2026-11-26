"""Same search as sprint-3 (every Ns, Np, SLSQP on m_M and K_T) but the
inputs from the other subsystems are a dict so they can be swept.
"""

import numpy as np
from scipy.optimize import minimize

import propulsion_model as pm

BASE = {
    "lift_to_drag": 7.01,
    "fm": 0.708,
    "radius": 0.25,
    "omega_h": 300.0,
    "omega_c": 400.0,
    "eta_prop": 0.60,
    "k_download": 1.05,
    "airframe": 0.0,
    "r_ref": 0.0,
    "i0_scaled": False,
    "mm_max": 0.50,
    "kt_max": 0.060,
    "np_max": 8,
    "esc_volts": False,
    "datasheet_wind": False,
}

# final model: MN5212 resistance and no-load current, rewound with K_T.
# no mass term. m_M, K_T and Np upper bounds only close the search
FINAL = dict(BASE, datasheet_wind=True, kt_max=0.3, mm_max=3.0, np_max=20, esc_volts=True)

# T-Motor ALPHA 60A 12S: 18-50.4 V input. Molicel P45B: 4.2 V per cell at full charge
V_ESC_MAX = 50.4
V_CELL_CHARGE = 4.2

LO = np.array([0.05, 0.015])
NS = range(4, 17)
NAMES = ["g1", "g2", "g3", "g4", "g5", "g6", "g7"]

# rough size of each constraint so SLSQP sees them on the same scale
SCALE = np.array([300.0, 1000.0, 0.5, 0.5, 10.0, 10.0, 10.0])

# winding resistance of the reference motor is given at this K_T
KT_REF = 0.03
# MN5212 KV340: 69 mOhm, 1.1 A idle, Kv 340 -> K_T = 60/(2*pi*340)
KT_DS = 0.028
R_DS = 0.069
I0_DS = 1.1


def r_winding(m_M, k_t, r_ref, datasheet=False):
    if datasheet:
        return R_DS * (k_t / KT_DS) ** 2
    # R goes with K_T^2 for the same copper, and with 1/L^5 = m^(-5/3) for size
    return r_ref * (k_t / KT_REF) ** 2 * (pm.M_REF / m_M) ** (5 / 3)


def run(ns, np_, m_M, k_t, p, mission):
    design = {"Ns": ns, "Np": np_, "m_M": m_M, "K_T": k_t}
    design["R_m"] = r_winding(m_M, k_t, p["r_ref"], p["datasheet_wind"])
    # i0 = 2 A at the reference motor; scaled means no-load power does not depend on K_T
    pm.I_0 = 2.0
    if p["datasheet_wind"]:
        pm.I_0 = I0_DS * KT_DS / k_t
    elif p["i0_scaled"]:
        pm.I_0 = 2.0 * KT_REF / k_t
    craft = pm.default_aircraft(pm.PAYLOAD_KG)
    for k in ["lift_to_drag", "fm", "omega_h", "omega_c", "eta_prop", "k_download"]:
        craft[k] = p[k]
    craft["a_rot"] = np.pi * p["radius"] ** 2
    try:
        out = pm.close_mass(design, craft, mission, pm.PAYLOAD_KG, p["airframe"])
    except (RuntimeError, OverflowError):
        return None
    if p["esc_volts"]:
        out["g"]["g7 ESC voltage V"] = ns * V_CELL_CHARGE - V_ESC_MAX
    return out


def n_con(p):
    if p["esc_volts"]:
        return 7
    return 6


def to_x(z, p):
    hi = np.array([p["mm_max"], p["kt_max"]])
    return LO + np.clip(z, 0, 1) * (hi - LO)


def solve_pack(ns, np_, p, mission, starts=([0.5, 0.5],), method="SLSQP", objective="f"):
    cache = {}

    def out(z):
        key = (round(z[0], 10), round(z[1], 10))
        if key not in cache:
            x = to_x(z, p)
            cache[key] = run(ns, np_, x[0], x[1], p, mission)
        return cache[key]

    def f(z):
        o = out(z)
        if o is None:
            return 1e3
        return o[objective]

    n = n_con(p)

    def con(z):
        o = out(z)
        if o is None:
            return -np.ones(n)
        return -np.array(list(o["g"].values())) / SCALE[:n]

    # g7 only depends on Ns, so it is checked after and not given to the optimizer
    def con6(z):
        return con(z)[:6]

    runs = []
    for z0 in starts:
        if method == "SLSQP":
            r = minimize(f, z0, method="SLSQP", bounds=[(0, 1), (0, 1)],
                         constraints={"type": "ineq", "fun": con6},
                         options={"maxiter": 200, "ftol": 1e-10})
        else:
            r = minimize(f, z0, method="COBYLA", bounds=[(0, 1), (0, 1)],
                         constraints={"type": "ineq", "fun": con6},
                         options={"maxiter": 2000, "rhobeg": 0.1, "tol": 1e-8})
        o = out(r.x)
        ok = o is not None and min(con(r.x)) > -1e-6
        x = to_x(r.x, p)
        bounds = []
        if r.x[0] < 1e-4 or r.x[0] > 1 - 1e-4:
            bounds.append("m_M bound")
        if r.x[1] < 1e-4 or r.x[1] > 1 - 1e-4:
            bounds.append("K_T bound")
        runs.append({"Ns": ns, "Np": np_, "m_M": x[0], "K_T": x[1], "z0": z0,
                     "ok": ok, "out": o, "nit": r.get("nit"), "nfev": r.nfev,
                     "msg": r.message, "slack": con(r.x), "bounds": bounds})
    return runs


def search(p, mission, starts=([0.5, 0.5],), objective="f"):
    best = None
    table = {}
    for ns in NS:
        if p["esc_volts"] and ns * V_CELL_CHARGE > V_ESC_MAX + 1e-9:
            continue
        for np_ in range(1, p["np_max"] + 1):
            runs = [r for r in solve_pack(ns, np_, p, mission, starts, objective=objective) if r["ok"]]
            if not runs:
                continue
            r = min(runs, key=lambda r: r["out"][objective])
            table[(ns, np_)] = r
            if best is None or r["out"][objective] < best["out"][objective]:
                best = r
    return best, table


def active(r, tol=1e-3):
    return [n for n, s in zip(NAMES, r["slack"]) if s < tol] + r["bounds"]
