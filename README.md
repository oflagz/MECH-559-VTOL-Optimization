# MECH 559 — VTOL Tiltrotor Hybrid UAV Optimization

Python models for the MECH 559 multidisciplinary design project: **energy optimization of a VTOL tiltrotor hybrid UAV**.

The system is decomposed into **VTOL rotors**, **wing**, and **electric propulsion / battery** subsystems. This repository currently focuses on the rotor model and establishes interfaces for later coupled system optimization.

## Current rotor formulation

Continuous rotor design vector:

```text
x = [R, c, omega_h, omega_c, theta_tw]
```

Blade count `Nb` is enumerated over `{2, 3, 4}`. The current optimizer uses SciPy SLSQP with multiple starting points.

Implemented:
- analytical hover momentum/profile-power model from the proposal;
- preliminary axial blade-element/momentum model for thrust and cruise power;
- rotor mission-energy objective;
- proposal constraints g1-g10;
- blade-count enumeration and multistart SLSQP;
- regression test for the proposal hover screening calculation.

The BEMT implementation is deliberately **preliminary**. It currently uses a smooth surrogate airfoil polar, uniform induced velocity and no swirl. It must be upgraded and validated before optimized dimensions are treated as final engineering results.

## Repository structure

```text
src/
  rotor/          rotor geometry, hover, BEMT, objective and constraints
  wing/           wing subsystem interface (placeholder)
  propulsion/     propulsion/battery interface (placeholder)
  system/         coupled system model (placeholder)
optimization/     SLSQP optimization routines
scripts/          executable validation/optimization entry points
tests/            regression and model tests
docs/             assumptions and model-status documentation
data/             validation/polar/mass data (added as selected)
results/          generated figures/tables (not source data)
```

## Setup

```bash
python -m venv .venv
# Windows
.venv\Scripts\activate
# macOS/Linux
source .venv/bin/activate

pip install -r requirements.txt
```

## Validate the proposal screening point

```bash
python scripts/validate_screening_case.py
```

Expected values are approximately 12.88 N required thrust per rotor, 0.115 solidity, 94.2 W shaft power per rotor, figure of merit 0.708, and 443 W total hover electrical power.

## Run the preliminary rotor optimization

```bash
python scripts/run_rotor_optimization.py
```

## Run tests

```bash
pytest
```

## Important model-status note

The rotor section of the proposal contains provisional screening mission values that do not yet fully match the project-level mission statement. See [docs/assumptions.md](docs/assumptions.md) before interpreting optimization results.

## Next development priorities

1. Add selected airfoil/propeller polar data with Reynolds-number interpolation.
2. Implement annular BEMT with axial and swirl induction.
3. Validate BEMT against published/manufacturer propeller data.
4. Build rotor/motor mass surrogate from candidate hardware.
5. Connect wing drag/mount spacing and battery/power constraints.
6. Close the takeoff-mass/energy sizing loop and add system-level constraints.
