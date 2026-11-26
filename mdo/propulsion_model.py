"""Propulsion model from section 2.3.

Design variables are Ns, Np, motor mass, and KT. Rotor and wing numbers
come in as arguments. Running this file only checks the 9s3p hand point.

Run: python3 propulsion_model.py
"""

import math

G = 9.81
RHO = 1.225

# Middle of the ranges in the proposal table. Not a cell I actually picked.
V_CELL = 0.5 * (3.6 + 3.7)
Q_CELL = 0.5 * (4.0 + 5.0)
M_CELL = 0.5 * (0.060 + 0.070)
R_CELL = 0.5 * (0.010 + 0.020)
C_MAX = 0.5 * (8.0 + 15.0)
DOD = 0.5 * (0.70 + 0.85)
GCTP = 0.5 * (0.65 + 0.80)
I_0 = 0.5 * (1.0 + 3.0)
M_REF = 0.5 * (0.2 + 0.4)
Q_NOM_REF = 0.5 * (1.0 + 1.1)
Q_MAX_REF = 0.5 * (1.5 + 1.6)
ETA_ESC = 0.5 * (0.90 + 0.95)
K_PEAK = 0.5 * (1.3 + 1.5)
I_ESC_MAX = 60.0  # the writeup says 60 A, not the middle of 40-80

# Bounds from the proposal update.
BOUNDS = {
    "Ns": (4.0, 16.0),
    "Np": (1.0, 8.0),
    "m_M": (0.05, 0.50),
    "K_T": (0.015, 0.060),
}

# The 9s3p hand check. Airframe is 0 because wing and rotor mass are not in yet.
PAYLOAD_KG = 5.0
AIRFRAME_KG = 0.0
ROTOR_RADIUS_M = 0.25
FEASIBLE_DESIGN = {
    "Ns": 9.0,
    "Np": 3.0,
    "m_M": 0.22,
    "K_T": 0.020,
    "R_m": 0.0,
}
CLOSED_MASS_KG = 8.47


def default_mission():
    """5 min hover and 40 km in 20 min. Reserve and transition are still zero."""
    cruise_mps = 40.0e3 / (20.0 * 60.0)
    return {
        "hover_s": 300.0,
        "reserve_s": 0.0,
        "cruise_s": 20.0 * 60.0,
        "cruise_mps": cruise_mps,
        "transition_s": 0.0,
    }


def default_aircraft(m_to):
    """Screening rotor from section 2.2. L/D is 7.01 from that check, not Antoine's wing."""
    area = math.pi * ROTOR_RADIUS_M ** 2
    return {
        "m_to": m_to,
        "n_rot": 4,
        "a_rot": area,
        "k_download": 1.05,
        "fm": 0.708,
        "omega_h": 300.0,
        "omega_c": 400.0,
        "lift_to_drag": 7.01,
        "eta_prop": 0.60,
    }


def hover_powers(m_to, k_download, n_rot, a_rot, fm):
    """Hover power from momentum theory. FM divides the ideal power.

    The proposal puts FM under the square root. That is not what figure of
    merit means, so I left it out of the square root.
    """
    thrust = k_download * m_to * G
    ideal = thrust ** 1.5 / math.sqrt(2.0 * RHO * n_rot * a_rot)
    return ideal, ideal / fm


def motor_torque(m_M):
    """Torque grows with motor mass the way the proposal says."""
    scale = (m_M / M_REF) ** (3.5 / 3.0)
    return Q_NOM_REF * scale, Q_MAX_REF * scale


def evaluate(design, aircraft, mission):
    """One pass at a fixed takeoff mass. A constraint is ok when g is below zero."""
    ns = design["Ns"]
    np_ = design["Np"]
    m_M = design["m_M"]
    k_t = design["K_T"]
    r_m = design["R_m"]
    m_to = aircraft["m_to"]
    n_rot = aircraft["n_rot"]

    v_bus = ns * V_CELL
    q_pack = np_ * Q_CELL
    e_pack = v_bus * q_pack
    r_pack = (ns / np_) * R_CELL
    m_batt = ns * np_ * M_CELL / GCTP

    ideal_h, p_shaft_h = hover_powers(
        m_to, aircraft["k_download"], n_rot, aircraft["a_rot"], aircraft["fm"]
    )
    q_h = p_shaft_h / (n_rot * aircraft["omega_h"])
    i_h = q_h / k_t + I_0
    v_h = k_t * aircraft["omega_h"] + i_h * r_m
    p_h = n_rot * v_h * i_h / ETA_ESC
    i_pack_h = p_h / v_bus

    p_shaft_c = (
        m_to * G * mission["cruise_mps"]
        / (aircraft["lift_to_drag"] * aircraft["eta_prop"])
    )
    q_c = p_shaft_c / (n_rot * aircraft["omega_c"])
    i_c = q_c / k_t + I_0
    v_c = k_t * aircraft["omega_c"] + i_c * r_m
    p_c = n_rot * v_c * i_c / ETA_ESC

    e_mission = (
        p_h * (mission["hover_s"] + mission["reserve_s"])
        + p_c * mission["cruise_s"]
        + (p_h + p_c) * mission["transition_s"]
    ) / 3600.0

    m_esc = 8.42e-4 * K_PEAK * i_h
    f = m_batt + n_rot * (m_M + m_esc)
    q_nom, q_max = motor_torque(m_M)

    g = {
        "g1 energy Wh": e_mission - e_pack * DOD,
        "g2 hover power W": p_h - C_MAX * q_pack * v_bus,
        "g3 continuous torque Nm": q_h - q_nom,
        "g4 peak torque Nm": K_PEAK * q_h - q_max,
        "g5 bus voltage V": v_h - (v_bus - i_pack_h * r_pack),
        "g6 ESC current A": i_h - I_ESC_MAX,
    }
    return {
        "f": f,
        "g": g,
        "m_batt": m_batt,
        "m_esc": m_esc,
        "v_bus": v_bus,
        "e_pack": e_pack,
        "e_mission": e_mission,
        "ideal_h": ideal_h,
        "p_shaft_h": p_shaft_h,
        "p_shaft_c": p_shaft_c,
        "p_h": p_h,
        "p_c": p_c,
        "q_h": q_h,
        "q_c": q_c,
        "q_nom": q_nom,
        "q_max": q_max,
        "i_h": i_h,
        "v_h": v_h,
        "i_pack_h": i_pack_h,
        "r_pack": r_pack,
    }


def close_mass(design, aircraft, mission, payload_kg, airframe_kg):
    """Takeoff mass has to equal payload plus airframe plus this subsystem."""
    m_to = payload_kg + airframe_kg + 1.0
    for _ in range(40):
        aircraft["m_to"] = m_to
        out = evaluate(design, aircraft, mission)
        updated = payload_kg + airframe_kg + out["f"]
        if abs(updated - m_to) < 1e-9:
            out["m_to"] = updated
            out["h_mass"] = updated - (payload_kg + airframe_kg + out["f"])
            return out
        m_to = updated
    raise RuntimeError("takeoff mass did not converge")


def constraints_feasible(out, tol=1e-6):
    return all(residual <= tol for residual in out["g"].values())


def _fail(failures, message):
    failures.append(message)


def _near(actual, expected, rel_tol, failures, message):
    if expected == 0.0:
        if abs(actual) > rel_tol:
            _fail(failures, message)
        return
    if abs(actual - expected) / abs(expected) > rel_tol:
        _fail(failures, message)


def _verify_rotor_screening_algebra(failures):
    """Check the hover numbers printed in the rotor screening table."""
    radius = ROTOR_RADIUS_M
    area = math.pi * radius ** 2
    weight = 5.0 * G
    thrust_one = 1.05 * weight / 4.0
    ideal_one = thrust_one ** 1.5 / math.sqrt(2.0 * RHO * area)
    _near(thrust_one, 12.88, 0.002, failures, "hover thrust does not match the printed 12.88 N")
    _near(ideal_one, 66.7, 0.01, failures, "ideal hover power does not match the printed 66.7 W")

    profile = 17.5
    kappa = 1.15
    shaft_one = kappa * ideal_one + profile
    fm = ideal_one / shaft_one
    ideal, shaft = hover_powers(5.0, 1.05, 4, area, fm)
    if abs(ideal / 4.0 - ideal_one) > 1e-6:
        _fail(failures, "total ideal power is not four times the one-rotor value")
    if abs(shaft / 4.0 - shaft_one) > 1e-6:
        _fail(failures, "shaft power does not match ideal power divided by FM")

    wrong = (1.05 * weight) ** 1.5 / math.sqrt(2.0 * RHO * 4.0 * area * fm)
    if abs(shaft - wrong) / shaft < 0.02:
        _fail(failures, "FM-under-the-square-root form was not actually different")


def _verify_pack_formulas(failures):
    m_batt = 9 * 3 * M_CELL / GCTP
    e_pack = (9 * V_CELL) * (3 * Q_CELL)
    _near(m_batt, 2.42, 0.01, failures, "9s3p battery mass does not match 2.42 kg")
    _near(e_pack, 443.0, 0.01, failures, "9s3p pack energy does not match 443 Wh")

    scale = (0.22 / M_REF) ** (3.5 / 3.0)
    q_nom, q_max = motor_torque(0.22)
    if abs(q_nom - Q_NOM_REF * scale) > 1e-12:
        _fail(failures, "continuous torque scaling does not match the proposal")
    if abs(q_max - Q_MAX_REF * scale) > 1e-12:
        _fail(failures, "peak torque scaling does not match the proposal")


def _verify_feasible_point(failures):
    mission = default_mission()
    design = dict(FEASIBLE_DESIGN)
    aircraft = default_aircraft(CLOSED_MASS_KG)
    closed = close_mass(design, aircraft, mission, PAYLOAD_KG, AIRFRAME_KG)

    _near(closed["m_to"], CLOSED_MASS_KG, 0.002, failures, "closed mass is not 8.47 kg")
    if abs(closed["h_mass"]) > 1e-8:
        _fail(failures, "mass closure residual is not zero")
    if not constraints_feasible(closed):
        _fail(failures, "feasible point violates a constraint")

    _near(closed["m_batt"], 2.42, 0.01, failures, "battery mass at closed point")
    _near(closed["e_pack"], 443.0, 0.01, failures, "pack energy at closed point")
    _near(closed["ideal_h"], 588.0, 0.02, failures, "ideal hover power at closed point")
    _near(closed["p_shaft_h"], 830.0, 0.02, failures, "hover shaft power at closed point")
    _near(closed["p_shaft_c"], 659.0, 0.02, failures, "cruise shaft power at closed point")
    _near(closed["p_h"], 950.0, 0.02, failures, "hover electrical power at closed point")
    _near(closed["p_c"], 782.0, 0.02, failures, "cruise electrical power at closed point")
    _near(closed["e_mission"], 340.0, 0.02, failures, "mission energy at closed point")
    _near(closed["q_h"], 0.692, 0.02, failures, "hover torque at closed point")
    _near(closed["i_h"], 36.6, 0.02, failures, "hover current at closed point")
    _near(closed["f"], 3.47, 0.02, failures, "propulsion mass at closed point")

    if abs(closed["i_pack_h"] - closed["p_h"] / closed["v_bus"]) > 1e-9:
        _fail(failures, "pack current is not electrical power divided by bus voltage")
    if design["R_m"] == 0.0 and abs(closed["v_h"] - design["K_T"] * 300.0) > 1e-9:
        _fail(failures, "with R_m = 0 the motor voltage should be only the back-EMF")

    usable = closed["e_pack"] * DOD
    if closed["e_mission"] >= usable:
        _fail(failures, "mission energy is not below usable pack energy")

    return closed


def verify():
    failures = []
    _verify_rotor_screening_algebra(failures)
    _verify_pack_formulas(failures)
    closed = _verify_feasible_point(failures)

    print("feasible point Ns=9 Np=3 m_M=0.22 K_T=0.020 (Sprint 1 / section 2.3)")
    print(
        "closed mass %.2f kg  propulsion mass %.3f kg  battery %.3f kg"
        % (closed["m_to"], closed["f"], closed["m_batt"])
    )
    print(
        "pack %.0f Wh  usable %.0f Wh  mission %.0f Wh"
        % (closed["e_pack"], closed["e_pack"] * DOD, closed["e_mission"])
    )
    print(
        "hover shaft %.0f W  cruise shaft %.0f W  hover I %.1f A"
        % (closed["p_shaft_h"], closed["p_shaft_c"], closed["i_h"])
    )
    for name, residual in closed["g"].items():
        print("  %-28s %+10.3f" % (name, residual))

    if failures:
        print("VERIFICATION FAILED")
        for item in failures:
            print("  " + item)
        return 1
    print("checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(verify())
