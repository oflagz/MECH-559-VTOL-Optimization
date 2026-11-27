import argparse
import json
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, minimize

TEAM = Path(__file__).resolve().parent
sys.path.insert(0, str(TEAM))

import wing_model as wing
import propulsion_model as pm
from solver import KT_DS, R_DS, I0_DS, V_CELL_CHARGE, V_ESC_MAX
from rotor_design.bemt import RotorGeometry, SectionPolar, performance, SPEED_SOUND

POLAR = SectionPolar(**json.loads((TEAM / "results/tables/study_summary.json").read_text())["calibrated_polar"])
THETA75 = 20.0
PAYLOAD = 5.0
FRAME = 0.55
NACELLE = 0.09
M_PROP_REF = 0.42 * 0.028349523125
R_PROP_REF = 0.127
C_PROP_REF = 0.023
N_ROT = 4
VC = wing.P["Vc"]
TH = wing.P["th"]
TC = wing.P["tc"]
OH = (180.0, 730.0)
OC = (200.0, 780.0)
N_SECTIONS = 32
CACHE = {}

LO_W, HI_W = wing.bounds()
LO = np.r_[LO_W, 0.185, 0.009, -60.0, 0.05, 0.015]
HI = np.r_[HI_W, 0.297, 0.034, -12.0, 3.0, 0.30]
SPAN = HI - LO
G_SCALE = np.r_[np.ones(5), 450.0, 10000.0, 1.0, 1.0, 50.0, 60.0, 50.0,
                np.ones(8), 50.0, 1.0, 60.0, 5000.0]
G_NAMES = [
    "wing_stall", "wing_stress", "wing_deflection", "wing_download", "wing_min_chord",
    "battery_energy", "battery_hover_power", "motor_continuous_torque", "motor_peak_torque",
    "motor_bus_voltage", "ESC_hover_current", "ESC_input_voltage",
    "rotor_hover_FM", "rotor_hover_tip_Mach", "rotor_cruise_tip_Mach", "cruise_advance_ratio",
    "rotor_hover_incidence", "rotor_cruise_incidence", "hover_thrust_reserve", "rotor_clearance",
    "reserve_motor_bus_voltage", "reserve_motor_peak_torque", "reserve_ESC_current", "reserve_battery_power",
]

VERIFIED = {
    (11, 3, 3): np.array([1.68, 14.0, 0.11516564963467008, 0.4, 1.6250645208036834e-5,
                          0.22699669942291162, 0.009834637367046746, -27.424329816895412,
                          0.2123156756077415, 0.055441600658089465]),
    (7, 5, 3): np.array([1.68, 14.0, 0.11516564963467002, 0.4, 1.663603302922853e-5,
                         0.22810729598329244, 0.009873791433472012, -27.415367080779898,
                         0.2175157825256629, 0.035198873188427]),
    (9, 4, 3): np.array([1.68, 14.0, 0.11516564963467009, 0.4, 1.6753939348414293e-5,
                         0.228223715105616, 0.009876132116608543, -27.412669147340157,
                         0.2189333802419618, 0.04519069017820092]),
}


def set_sections(n):
    global N_SECTIONS
    if n < 16:
        raise ValueError("BEMT requires at least 16 radial sections")
    N_SECTIONS = int(n)
    CACHE.clear()


def prop_mass(radius, chord, blades):
    return M_PROP_REF * (radius / R_PROP_REF) ** 2.45 * (chord / C_PROP_REF) ** 0.35 * (blades / 2)


def rotor_point(R, c, tw, blades, omega, V):
    g = RotorGeometry(R, c, omega, tw, blades, THETA75)
    return performance(g, V, POLAR, n_sections=N_SECTIONS)


def trim(R, c, tw, blades, V, target, bounds):
    nodes = np.linspace(*bounds, 16)
    vals = []
    for omega in nodes:
        p = rotor_point(R, c, tw, blades, float(omega), V)
        vals.append(p.thrust - target if p.valid else np.nan)
    for k in range(len(nodes) - 1):
        a, b = vals[k:k + 2]
        if not np.isfinite(a) or not np.isfinite(b):
            continue
        if a == 0.0:
            omega = float(nodes[k])
            return omega, rotor_point(R, c, tw, blades, omega, V)
        if a * b < 0.0:
            def gap(omega):
                q = rotor_point(R, c, tw, blades, float(omega), V)
                if not q.valid:
                    raise ValueError("invalid BEMT trim")
                return q.thrust - target
            try:
                root = float(brentq(gap, float(nodes[k]), float(nodes[k + 1]), xtol=3e-5))
            except ValueError:
                continue
            return root, rotor_point(R, c, tw, blades, root, V)
    return None


def motor_point(q, omega, k_t):
    resistance = R_DS * (k_t / KT_DS) ** 2
    no_load = I0_DS * KT_DS / k_t
    current = q / k_t + no_load
    voltage = k_t * omega + current * resistance
    power = N_ROT * voltage * current / pm.ETA_ESC
    return current, voltage, power


def propulsion(ns, np_, m_M, k_t, q_h, om_h, q_c, om_c):
    bus = ns * pm.V_CELL
    pack_capacity = np_ * pm.Q_CELL
    pack_energy = bus * pack_capacity
    pack_resistance = (ns / np_) * pm.R_CELL
    battery_mass = ns * np_ * pm.M_CELL / pm.GCTP
    i_h, v_h, p_h = motor_point(q_h, om_h, k_t)
    i_c, v_c, p_c = motor_point(q_c, om_c, k_t)
    if min(i_h, v_h, i_c, v_c) <= 0.0:
        return None
    energy = (p_h * TH + p_c * TC) / 3600.0
    esc_mass = 8.42e-4 * pm.K_PEAK * i_h
    mass = battery_mass + N_ROT * (m_M + esc_mass)
    q_nom, q_max = pm.motor_torque(m_M)
    hover_bus = bus - p_h / bus * pack_resistance
    g = np.array([
        energy - pack_energy * pm.DOD,
        p_h - pm.C_MAX * pack_capacity * bus,
        q_h - q_nom,
        pm.K_PEAK * q_h - q_max,
        v_h - hover_bus,
        i_h - pm.I_ESC_MAX,
        ns * V_CELL_CHARGE - V_ESC_MAX,
    ])
    return mass, energy, g


def reserve_electrical(ns, np_, m_M, k_t, omega, reserve):
    current, voltage, power = motor_point(reserve.torque, omega, k_t)
    bus = ns * pm.V_CELL
    pack_resistance = (ns / np_) * pm.R_CELL
    available_voltage = bus - power / bus * pack_resistance
    _, max_torque = pm.motor_torque(m_M)
    g = np.array([
        voltage - available_voltage,
        reserve.torque - max_torque,
        current - pm.I_ESC_MAX,
        power - pm.C_MAX * np_ * pm.Q_CELL * bus,
    ])
    return g, {"voltage_required": voltage, "voltage_available": available_voltage,
               "current": current, "torque": reserve.torque, "power": power}


def mda(z, ns, np_, blades, m_prop=3.58):
    b, AR, tcr, lam, A, R, c, tw, m_M, k_t = z
    xw = np.array([b, AR, tcr, lam, A])
    m_blade = prop_mass(R, c, blades)
    for _ in range(30):
        p = dict(wing.P)
        p["m_other"] = PAYLOAD + m_prop + N_ROT * (NACELLE + m_blade) + FRAME
        p["R"] = R
        w = wing.model(xw, p)
        hover = trim(R, c, tw, blades, 0.0, (w["W"] + w["Ddl"]) / N_ROT, OH)
        cruise = trim(R, c, tw, blades, VC, w["D"] / N_ROT, OC)
        if hover is None or cruise is None:
            return None
        om_h, h = hover
        om_c, cr = cruise
        drive = propulsion(ns, np_, m_M, k_t, h.torque, om_h, cr.torque, om_c)
        if drive is None:
            return None
        m_prop_new, energy, g_prop = drive
        m_new = PAYLOAD + w["m_wing"] + m_prop_new + N_ROT * (NACELLE + m_blade) + FRAME
        if abs(m_new - (p["m_other"] + w["m_wing"])) < 1e-5 and abs(m_prop_new - m_prop) < 1e-5:
            omega_reserve = 1.10 * om_h
            reserve = rotor_point(R, c, tw, blades, omega_reserve, 0.0)
            if not reserve.valid or not np.isfinite(h.fm):
                return None
            required = (w["W"] + w["Ddl"]) / N_ROT
            J = np.pi * VC / (om_c * R)
            g_rotor = np.array([
                (0.58 - h.fm) / 0.58,
                om_h * R / SPEED_SOUND - 0.65,
                np.hypot(VC, om_c * R) / SPEED_SOUND - 0.65,
                J / 0.96 - 1.0,
                max(abs(h.min_alpha_deg), abs(h.max_alpha_deg)) / 25.0 - 1.0,
                max(abs(cr.min_alpha_deg), abs(cr.max_alpha_deg)) / 25.0 - 1.0,
                (1.15 * required - reserve.thrust) / required,
                (2.0 * R + 0.05 - 0.65) / 0.65,
            ])
            g_reserve, electrical = reserve_electrical(ns, np_, m_M, k_t, omega_reserve, reserve)
            g = np.r_[wing.g_nl(xw, p), g_prop, g_rotor, g_reserve]
            return {"energy": energy, "mass": m_new, "g": g, "omega_h": om_h,
                    "omega_c": om_c, "m_wing": w["m_wing"], "m_prop": m_prop_new,
                    "D": w["D"], "Ddl": w["Ddl"], "J": J,
                    "reserve": electrical, "n_sections": N_SECTIONS}
        m_prop = m_prop_new
    return None


def evaluate(z, ns, np_, blades):
    key = (N_SECTIONS, ns, np_, blades) + tuple(np.round(z, 10))
    if key not in CACHE:
        CACHE[key] = mda(z, ns, np_, blades)
    return CACHE[key]


def solve_one(ns, np_, blades, maxiter=70, starts=None, sections=32):
    set_sections(sections)
    base = VERIFIED.get((ns, np_, blades))
    if base is None:
        base = VERIFIED[(11, 3, 3)].copy()
        base[9] = np.clip(base[9] * ns / 11, LO[9], HI[9])
        if blades == 2:
            base[6], base[7] = 0.014, -27.8
        elif blades == 4:
            base[6], base[7] = 0.009, -25.0
    candidates = [base] if starts is None else [np.asarray(x, dtype=float) for x in starts]
    best = None

    def check_mesh(x):
        grid = N_SECTIONS
        set_sections(128)
        out = mda(x, ns, np_, blades)
        set_sections(grid)
        residual = float(np.max(out["g"] / G_SCALE)) if out is not None else float("inf")
        return out, residual

    for start in candidates:
        CACHE.clear()
        u0 = np.clip((start - LO) / SPAN, 0.0, 1.0)

        def get(u):
            return evaluate(LO + np.asarray(u) * SPAN, ns, np_, blades)

        def objective(u):
            out = get(u)
            return out["energy"] / 400.0 if out is not None else 1e3

        def constraint(u):
            out = get(u)
            return -out["g"] / G_SCALE if out is not None else -np.ones(len(G_NAMES))

        result = minimize(objective, u0, method="SLSQP", bounds=[(0.0, 1.0)] * 10,
                          constraints={"type": "ineq", "fun": constraint},
                          options={"maxiter": maxiter, "ftol": 1e-9, "eps": 2e-5})
        x = LO + result.x * SPAN
        out, residual = check_mesh(x)
        row = {"Ns": ns, "Np": np_, "blades": blades, "ok": bool(residual <= 1e-5),
               "x": x, "out": out, "nit": result.nit, "nfev": result.nfev,
               "max_scaled_violation": residual, "verified_sections": 128,
               "message": str(result.message)}
        if best is None or ((row["ok"], -(out["energy"] if out else float("inf"))) >
                            (best["ok"], -(best["out"]["energy"] if best["out"] else float("inf")))):
            best = row

    if (ns, np_, blades) in VERIFIED:
        seed = VERIFIED[(ns, np_, blades)]
        out, residual = check_mesh(seed)
        if residual <= 1e-5 and (not best["ok"] or out["energy"] < best["out"]["energy"]):
            best = {"Ns": ns, "Np": np_, "blades": blades, "ok": True,
                    "x": seed.copy(), "out": out, "nit": 0, "nfev": 0,
                    "max_scaled_violation": residual, "verified_sections": 128,
                    "message": "verified reference design"}
    return best


def search():
    cases = [(ns, np_, blades) for blades in (2, 3, 4) for ns in range(4, 13) for np_ in range(1, 21)]
    best = None
    with ProcessPoolExecutor() as pool, (TEAM / "search_log.txt").open("w") as log:
        futures = [pool.submit(solve_one, ns, np_, blades) for ns, np_, blades in cases]
        for i, future in enumerate(as_completed(futures), 1):
            row = future.result()
            out = row["out"]
            energy = None if out is None else round(out["energy"], 2)
            line = (f"{i}/{len(cases)} {row['Ns']}s{row['Np']}p "
                    f"blades={row['blades']} ok={row['ok']} E={energy}")
            print(line, flush=True)
            log.write(line + "\n")
            log.flush()
            if row["ok"] and (best is None or out["energy"] < best["out"]["energy"]):
                best = row
    return best


def verify():
    for (ns, np_, blades), z in VERIFIED.items():
        for sections in (32, 64, 96, 128):
            set_sections(sections)
            out = mda(z, ns, np_, blades)
            if out is None:
                raise RuntimeError(f"BEMT trim failed: {ns}s{np_}p, {sections} sections")
            g = out["g"] / G_SCALE
            failed = [G_NAMES[i] for i in np.flatnonzero(g > 1e-5)]
            print(f"{ns}s{np_}p {sections:3d} sections: E={out['energy']:.3f} Wh "
                  f"mass={out['mass']:.4f} kg J={out['J']:.6f} "
                  f"max_g={max(g):+.2e} feasible={not failed} "
                  f"failed={failed}", flush=True)
            if failed:
                raise RuntimeError(f"Constraints violated: {failed}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--verify", action="store_true")
    parser.add_argument("--solve", nargs=3, type=int, metavar=("NS", "NP", "BLADES"))
    parser.add_argument("--sections", type=int, default=32)
    args = parser.parse_args()
    if args.verify:
        verify()
    elif args.solve:
        result = solve_one(*args.solve, sections=args.sections)
        print("feasible:", result["ok"])
        if result["out"]:
            out = result["out"]
            print("energy:", out["energy"], "mass:", out["mass"], "J:", out["J"])
            print("x:", np.array2string(result["x"], precision=9))
            print("g:", out["g"])
    else:
        best = search()
        if best is None:
            print("no feasible design")
        else:
            out = best["out"]
            print("best", best["Ns"], best["Np"], best["blades"],
                  round(out["energy"], 2), round(out["mass"], 3))
            print("x", np.array2string(best["x"], precision=9))
