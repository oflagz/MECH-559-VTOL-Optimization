"""Tries every cell count and lets SLSQP pick the motor.

Rotor and wing numbers are still the screening set in INTERFACE.
Swap them when Ashish and Antoine actually optimize, not from one feasible point.

Run: python3 optimize.py
"""

import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, str(Path(__file__).resolve().parent))
import propulsion_model as pm

# Still the section 2.2 screening numbers. The mission itself is already
# the one in the introduction. Don't drop in a single feasible-point L/D.
INTERFACE = {
    "k_download": 1.05,
    "fm": 0.708,
    "omega_h": 300.0,
    "omega_c": 400.0,
    "lift_to_drag": 7.01,
    "eta_prop": 0.60,
    "airframe_kg": pm.AIRFRAME_KG,
}

LO = np.array([pm.BOUNDS["m_M"][0], pm.BOUNDS["K_T"][0]])
HI = np.array([pm.BOUNDS["m_M"][1], pm.BOUNDS["K_T"][1]])
NS_RANGE = range(int(pm.BOUNDS["Ns"][0]), int(pm.BOUNDS["Ns"][1]) + 1)
NP_RANGE = range(int(pm.BOUNDS["Np"][0]), int(pm.BOUNDS["Np"][1]) + 1)
STARTS = (
    ("middle", np.array([0.5, 0.5])),
    (
        "hand check",
        (np.array([pm.FEASIBLE_DESIGN["m_M"], pm.FEASIBLE_DESIGN["K_T"]]) - LO) / (HI - LO),
    ),
)


def aircraft(m_to):
    craft = pm.default_aircraft(m_to)
    craft.update({k: INTERFACE[k] for k in INTERFACE if k != "airframe_kg"})
    return craft


def evaluate(ns, np_, m_M, k_t, mission):
    design = {
        "Ns": float(ns),
        "Np": float(np_),
        "m_M": float(m_M),
        "K_T": float(k_t),
        "R_m": 0.0,
    }
    try:
        return pm.close_mass(
            design, aircraft(pm.PAYLOAD_KG), mission, pm.PAYLOAD_KG, INTERFACE["airframe_kg"]
        )
    except RuntimeError:
        return None


def solve_pair(ns, np_, mission):
    found = []

    def physical(z):
        return LO + np.clip(z, 0.0, 1.0) * (HI - LO)

    cache = {}

    def out_of(z):
        key = (round(float(z[0]), 8), round(float(z[1]), 8))
        if key not in cache:
            x = physical(z)
            cache[key] = evaluate(ns, np_, x[0], x[1], mission)
        return cache[key]

    def objective(z):
        out = out_of(z)
        return 1e6 if out is None else out["f"]

    def constraint(name):
        def fun(z):
            out = out_of(z)
            if out is None:
                return -1.0
            return -out["g"][name]
        return fun

    probe = evaluate(ns, np_, 0.22, 0.020, mission)
    if probe is None:
        return None
    names = list(probe["g"])
    cons = [{"type": "ineq", "fun": constraint(name)} for name in names]

    for name, z0 in STARTS:
        result = minimize(
            objective,
            z0,
            method="SLSQP",
            bounds=[(0.0, 1.0), (0.0, 1.0)],
            constraints=cons,
            options={"maxiter": 80, "ftol": 1e-9, "disp": False},
        )
        if not result.success:
            continue
        x = physical(result.x)
        out = evaluate(ns, np_, x[0], x[1], mission)
        if out is None or not pm.constraints_feasible(out, tol=1e-4):
            continue
        out = dict(out)
        out["start"] = name
        out["Ns"] = ns
        out["Np"] = np_
        out["m_M"] = float(x[0])
        out["K_T"] = float(x[1])
        found.append(out)
    return found


def on_bound(out):
    hits = []
    if out["m_M"] - LO[0] < 1e-3 or HI[0] - out["m_M"] < 1e-3:
        hits.append("m_M")
    if out["K_T"] - LO[1] < 1e-4 or HI[1] - out["K_T"] < 1e-4:
        hits.append("K_T")
    return hits


def sensitivity(best, mission):
    base = best["f"]
    rows = []
    cases = (
        ("L/D", "lift_to_drag", 7.71, 6.31),
        ("FM", "fm", 0.779, 0.637),
        ("hover_s", "hover_s", 330.0, 270.0),
    )
    for label, key, up, down in cases:
        deltas = []
        for value in (up, down):
            if key in INTERFACE:
                saved = INTERFACE[key]
                INTERFACE[key] = value
            else:
                saved = mission[key]
                mission[key] = value
            out = evaluate(best["Ns"], best["Np"], best["m_M"], best["K_T"], mission)
            if key in INTERFACE:
                INTERFACE[key] = saved
            else:
                mission[key] = saved
            deltas.append(None if out is None else (out["f"] - base, out["e_mission"] - best["e_mission"]))
        rows.append((label, deltas[0], deltas[1]))
    return rows


def main():
    mission = pm.default_mission()
    best = None
    best_starts = []
    n_ok = 0
    for ns in NS_RANGE:
        for np_ in NP_RANGE:
            found = solve_pair(ns, np_, mission)
            if not found:
                continue
            n_ok += 1
            out = min(found, key=lambda item: item["f"])
            if best is None or out["f"] < best["f"]:
                best = out
                best_starts = found

    if best is None:
        print("no feasible point")
        return 1

    usable = best["e_pack"] * pm.DOD
    print(
        "best Ns=%d Np=%d m_M=%.3f kg K_T=%.4f" % (best["Ns"], best["Np"], best["m_M"], best["K_T"])
    )
    print("starts on that pack:")
    seen = {item["start"] for item in best_starts}
    for item in best_starts:
        print(
            "  %-12s m_M=%.3f kg  K_T=%.4f  f=%.3f kg"
            % (item["start"], item["m_M"], item["K_T"], item["f"])
        )
    for name, _ in STARTS:
        if name not in seen:
            print("  %-12s no feasible motor" % name)
    print(
        "f=%.3f kg  m_to=%.2f kg  mission=%.0f Wh  usable=%.0f Wh"
        % (best["f"], best["m_to"], best["e_mission"], usable)
    )
    n_pairs = (NS_RANGE.stop - NS_RANGE.start) * (NP_RANGE.stop - NP_RANGE.start)
    print("feasible cell pairs: %d / %d" % (n_ok, n_pairs))
    print("on bound:", ", ".join(on_bound(best)) or "none")
    for name, residual in best["g"].items():
        flag = " active" if residual > -1e-2 else ""
        print("  %-28s %+10.3f%s" % (name, residual, flag))

    print("sensitivity at the best pack (up / down): df kg, dE Wh")
    for key, up, down in sensitivity(best, mission):
        print(
            "  %-8s  df %+7.3f / %+7.3f   dE %+6.1f / %+6.1f"
            % (key, up[0], down[0], up[1], down[1])
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
