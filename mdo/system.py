import json
import sys
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


def prop_mass(radius, chord, blades):
    return M_PROP_REF * (radius / R_PROP_REF) ** 2.45 * (chord / C_PROP_REF) ** 0.35 * (blades / 2)
N_ROT = 4
VC = wing.P["Vc"]
TH = wing.P["th"]
TC = wing.P["tc"]
OH = (180.0, 730.0)
OC = (200.0, 780.0)
CACHE = {}


def rotor_point(R, c, tw, blades, omega, V):
    return performance(RotorGeometry(R, c, omega, tw, blades, THETA75), V, POLAR, n_sections=32)


def trim(R, c, tw, blades, V, target, bounds):
    lo, hi = bounds
    nodes = np.linspace(lo, hi, 16)
    vals = []
    for omega in nodes:
        p = rotor_point(R, c, tw, blades, float(omega), V)
        vals.append(p.thrust - target if p.valid else np.nan)
    for k in range(len(nodes) - 1):
        a, b = vals[k], vals[k + 1]
        if not (np.isfinite(a) and np.isfinite(b)):
            continue
        if a == 0.0:
            p = rotor_point(R, c, tw, blades, float(nodes[k]), V)
            return float(nodes[k]), p
        if a * b < 0.0:
            def gap(omega):
                q = rotor_point(R, c, tw, blades, float(omega), V)
                if not q.valid:
                    raise ValueError
                return q.thrust - target
            try:
                root = float(brentq(gap, float(nodes[k]), float(nodes[k + 1]), xtol=3e-5))
            except ValueError:
                continue
            return root, rotor_point(R, c, tw, blades, root, V)
    return None


def propulsion(ns, np_, m_M, k_t, q_h, om_h, q_c, om_c):
    r_m = R_DS * (k_t / KT_DS) ** 2
    i0 = I0_DS * KT_DS / k_t
    v_bus = ns * pm.V_CELL
    q_pack = np_ * pm.Q_CELL
    e_pack = v_bus * q_pack
    r_pack = (ns / np_) * pm.R_CELL
    m_batt = ns * np_ * pm.M_CELL / pm.GCTP

    def phase(q, om):
        i = q / k_t + i0
        v = k_t * om + i * r_m
        if i <= 0.0 or v <= 0.0:
            return None
        p_el = N_ROT * v * i / pm.ETA_ESC
        return i, v, p_el

    h = phase(q_h, om_h)
    c = phase(q_c, om_c)
    if h is None or c is None:
        return None
    i_h, v_h, p_h = h
    i_c, v_c, p_c = c
    e = (p_h * TH + p_c * TC) / 3600.0
    m_esc = 8.42e-4 * pm.K_PEAK * i_h
    mass = m_batt + N_ROT * (m_M + m_esc)
    q_nom, q_max = pm.motor_torque(m_M)
    i_pack = p_h / v_bus
    g = np.array([
        e - e_pack * pm.DOD,
        p_h - pm.C_MAX * q_pack * v_bus,
        q_h - q_nom,
        pm.K_PEAK * q_h - q_max,
        v_h - (v_bus - i_pack * r_pack),
        i_h - pm.I_ESC_MAX,
        ns * V_CELL_CHARGE - V_ESC_MAX,
    ])
    return mass, e, g


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
            reserve = rotor_point(R, c, tw, blades, 1.10 * om_h, 0.0)
            if not reserve.valid or not np.isfinite(h.fm):
                return None
            required = (w["W"] + w["Ddl"]) / N_ROT
            J = VC / ((om_c / (2.0 * np.pi)) * 2.0 * R)
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
            g = np.concatenate([wing.g_nl(xw, p), g_prop, g_rotor])
            return {
                "energy": energy,
                "mass": m_new,
                "g": g,
                "omega_h": om_h,
                "omega_c": om_c,
                "m_wing": w["m_wing"],
                "m_prop": m_prop_new,
                "D": w["D"],
                "Ddl": w["Ddl"],
            }
        m_prop = m_prop_new
    return None


def evaluate(z, ns, np_, blades):
    key = (ns, np_, blades) + tuple(np.round(z, 10))
    if key not in CACHE:
        CACHE[key] = mda(z, ns, np_, blades)
    return CACHE[key]


def solve_one(ns, np_, blades, maxiter=40):
    lo_w, hi_w = wing.bounds()
    lo = np.r_[lo_w, 0.185, 0.009, -60.0, 0.05, 0.015]
    hi = np.r_[hi_w, 0.297, 0.034, -12.0, 3.0, 0.30]
    x0 = np.array([1.507, 12.56, 0.10, 0.40, 12.2e-6, 0.2306, 0.010, -27.49, 0.213, 0.028])
    x0 = np.minimum(np.maximum(x0, lo), hi)
    n_fail = 20

    def fun(z):
        out = evaluate(z, ns, np_, blades)
        if out is None:
            return 1e6
        return out["energy"]

    def con(z):
        out = evaluate(z, ns, np_, blades)
        if out is None:
            return -np.ones(n_fail)
        return -out["g"]

    res = minimize(fun, x0, method="SLSQP", bounds=list(zip(lo, hi)),
                   constraints={"type": "ineq", "fun": con},
                   options={"maxiter": maxiter, "ftol": 1e-8})
    out = evaluate(res.x, ns, np_, blades)
    ok = out is not None and np.min(con(res.x)) >= -1e-4
    return {"Ns": ns, "Np": np_, "blades": blades, "ok": ok, "x": res.x, "out": out, "nit": res.nit, "nfev": res.nfev}


def search():
    from concurrent.futures import ProcessPoolExecutor, as_completed
    cases = [(ns, np_, blades) for blades in (2, 3, 4) for ns in range(4, 13) for np_ in range(1, 21)]
    best = None
    done = 0
    log = Path(__file__).resolve().parent / "search_log.txt"
    with ProcessPoolExecutor() as pool:
        futs = [pool.submit(solve_one, ns, np_, blades) for ns, np_, blades in cases]
        for fut in as_completed(futs):
            row = fut.result()
            done += 1
            energy = None if row["out"] is None else round(row["out"]["energy"], 2)
            line = f"{done}/{len(cases)} {row['Ns']}s{row['Np']}p blades={row['blades']} ok={row['ok']} E={energy} nit={row['nit']} nfev={row['nfev']}"
            print(line, flush=True)
            with log.open("a") as f:
                f.write(line + "\n")
            if row["ok"] and (best is None or row["out"]["energy"] < best["out"]["energy"]):
                best = row
    return best


if __name__ == "__main__":
    best = search()
    if best is None:
        print("no feasible design")
    else:
        o = best["out"]
        print("best", best["Ns"], best["Np"], best["blades"], round(o["energy"], 2), round(o["mass"], 3))
        print("x", np.array2string(best["x"], precision=5))
