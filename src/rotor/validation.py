"""Comparison against independent propeller test measurements.

Input CSV: rpm, velocity_mps, thrust_n, power_w, radius_m, chord_m,
twist_deg, blades, density_kg_m3. Optional: theta75_deg.
All geometry columns describe the tested propeller, not an inferred fit.
"""
import csv
from dataclasses import replace
import numpy as np
from .bemt import bemt_axial
from .parameters import Params

def validate_csv(path, polar, params=None):
    params = params or Params()
    output = []
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle):
            get = lambda key: float(row[key])
            p = replace(params, rho=get("density_kg_m3"),
                        theta75_deg=float(row.get("theta75_deg") or params.theta75_deg))
            pred = bemt_axial(get("radius_m"),get("chord_m"),
                              get("rpm")*2*np.pi/60,get("twist_deg"),
                              int(row["blades"]),get("velocity_mps"),p,polar=polar)
            thrust = get("thrust_n")
            power = get("power_w")
            output.append({
                "rpm": get("rpm"),"velocity_mps":get("velocity_mps"),
                "measured_thrust_n":thrust,"predicted_thrust_n":pred["T"],
                "thrust_error_n":pred["T"]-thrust,
                "measured_power_w":power,"predicted_power_w":pred["P_sh"],
                "power_error_w":pred["P_sh"]-power,
                "converged":pred["converged"],
            })
    return output

def summary(rows):
    if not rows:
        raise ValueError("No validation observations")
    return {
        "n":len(rows),
        "converged":sum(r["converged"] for r in rows),
        "thrust_mae_n":float(np.mean([abs(r["thrust_error_n"]) for r in rows])),
        "power_mae_w":float(np.mean([abs(r["power_error_w"]) for r in rows])),
    }
