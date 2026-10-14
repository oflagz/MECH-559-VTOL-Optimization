# MECH 559 | VTOL Tiltrotor Rotor Optimization (Rebuilt)

**Status:** fully reproducible preliminary *mathematical* optimization; **NOT** a verified real propeller design. Model inputs are exploratory and section polars are synthetic/illustrative. No measured thrust/power validation or rotor-mass coupling has yet been completed.

This repository was restarted to remove the previous mixed analytical/BEMT objective, which could assign shaft power independently of thrust capability. It now evaluates thrust and power at **the same rotor geometry, pitch and RPM**, and enforces both hover and cruise trim as equalities. The Git history retains the earlier prototype.

## Run locally (Python 3.10+)

```bash
python -m pip install -r requirements.txt
python -m pytest -q
python -m scripts.week1_study
```

All figures are regenerated in `results/figures/`, with full tables and optimization diagnostics in `results/tables/`.

## Design and mathematics

- Continuous decision variables: `x=[R,c,omega_hover,omega_cruise,theta75]`, with `theta75` explicitly optimized following the professor's feedback.
- Fixed rotor assemblies: 4. Enumerated blades per rotor: 2, 3, 4.
- Blade twist fixed to 0 during the preliminary screen (constant pitch), not optimized in place of theta75.
- SLSQP minimizes `4*(P_hover,electrical*t_hover+P_cruise,electrical*t_cruise)`; all terms use one blade-element and momentum model.
- Exact trim equality constraints `T_hover=T_required` and `T_cruise=D_cruise/4` prevent overthrust from masquerading as useful steady-flight energy.
- Inequalities enforce power, motor torque, tip Mach, mounting clearance, hover figure of merit and a 10% RPM reserve used to test >=15% hover thrust headroom. Solver outputs are independently re-evaluated and feasibility audited.
- Multistart optimizer, constrained blade-count comparison, feasible radius/pitch map, wing mass/cruise drag input sensitivity and blade-element convergence plots.
- **Mass penalty not yet included in the mission objective:** empirical rotor/motor mass data and the team's mass closure are not currently available. The local fixed-mass solution is not the global aircraft optimum.

## Important distinction: official mission vs provisional screening

**Project abstract:** 5 kg **payload**, 40 km sortie, 5 minutes hover, stated 20 minutes cruise. **Provisional rotor screening:** 5 kg **total takeoff mass**, 36 km at 20 m/s = 30 minutes cruise, 2 minutes hover. The project report does not reconcile these. Do not silently interchange them. See `docs/ASSUMPTIONS.md`.

## Repository map

- `rotoropt/model.py`: consistent BEMT performance, Prandtl loss corrections, independent analytic hover check
- `rotoropt/optimize.py`: exact-thrust constraints, geometry/power limits, scaled SLSQP multistart
- `scripts/week1_study.py`: reproducible tables and six sprint figures
- `tests/`: physics and solver regression tests
- `docs/ROADMAP.md`: four-week schedule and outputs
- `docs/INTERFACES.md`: input priorities, outputs for Antoine and Damian
- `docs/FEEDBACK_AUDIT.md`: professor feedback mapped to implementations
- `docs/SPRINT2_ROTORS.md`: report-ready progress with explicit limitations

## Validation pathway

For measured small-UAV propeller thrust/power data use the [UIUC Propeller Database](https://m-selig.ae.illinois.edu/props/propDB.html). Match a tested propeller with manufacturer **blade chord/twist distribution and airfoil polar data** before claiming BEMT validation; a propeller's nominal diameter and pitch alone are insufficient. See Week 2 plan.
