# MECH 559 · VTOL Tiltrotor Optimization

**McGill University · Engineering Systems Optimization · Fall 2026**

Collaborative project exploring mission-energy minimization for a four-rotor tiltrotor UAV. The aircraft is decomposed into **wing**, **rotor**, and **electric propulsion** subsystems. This repository hosts reusable engineering models, source data, numerical studies, and interfaces for multidisciplinary optimization.

## Project scope

| Subsystem | Scope | Current implementation |
| --- | --- | --- |
| Wing | Aerodynamic drag, structural mass, stall/deflection, rotor–wing download | `wing_model.py`, `wing_sprint3.py` |
| VTOL rotors | Hover/cruise thrust and shaft power, propeller selection, operating RPM, sensitivity | `rotor/`, `scripts/run_rotor_study.py` |
| Electric propulsion | Motor and battery sizing, electrical constraints, mass closure | `propulsion/` |

> **Mission distinction:** The original proposal calls for a **5 kg payload**, a **40 km sortie**, **5 minutes of hover**, and **20 minutes of cruise** (33.33 m/s). The rotor catalogue also includes an explicitly *alternative* 20 m/s cruise case to illustrate performance within measured propeller data. The alternative requires 33.33 minutes to cover 40 km and **does not satisfy the original 20-minute cruise requirement**.

## Quick start — rotors

Python 3.10+ is recommended. Run commands from the repository root.

```bash
python -m pip install -r requirements.txt
python -m pytest -q
python -m scripts.run_rotor_study
```

This regenerates `results/rotor/tables/` and `results/rotor/figures/`. Source-table transcription (usually unnecessary because the measured CSVs are already committed):

```bash
python data/rotor/make_source_tables.py
```

### Rotor model status

- The **active rotor sizing method** enumerates published APC propeller configurations and optimizes/solves hover and cruise RPM within **measured UIUC static and forward-flight performance ranges**, with no extrapolation.
- A signed-loading BEMT diagnostic is retained in `rotor/physics.py` to document why an older continuous-geometry surrogate optimum was rejected. It **is not** the active design selector.
- A provisional alternative case (base aircraft mass 9.5 kg including payload and excluding propellers; 7 N drag; 20 m/s cruise) selects an APC 18×12E, with approximately **278 Wh rotor electrical mission energy**. This is not a complete aircraft energy budget or a design for the original 20-minute cruise mission.
- **No candidate among the measured catalogue was supported for 33.33 m/s cruise**, within the provided test domains. This is a *data-envelope result*, not proof that an aircraft cannot fly that mission.
- Rotor, motor, and battery masses and efficiencies require common inputs. Nominal catalogue diameter versus actual test diameter remains an explicit source of uncertainty.

See [Rotor modeling and verification](docs/ROTOR_MODEL.md) and [experimental data provenance](data/rotor/SOURCES.md).

## Repository layout

```text
.github/workflows/     Automated checks for the rotor subsystem
rotor/                Rotor aerodynamic lookup, optimization, and verification
propulsion/           Electric propulsion component models
wing_model.py         Wing component model (team-owned)
wing_sprint3.py       Wing optimization studies (team-owned)
data/rotor/           Rotor source data and transcription script
results/rotor/        Generated rotor figures and numerical tables
scripts/              Reproducible analysis entry points
tests/                Rotor numerical/physics regression checks
docs/                 Shared technical documentation (no individual reports)
```

## Subsystem interfaces

| Quantity | Provider | Consumer |
| --- | --- | --- |
| Takeoff mass and payload allocation | System/propulsion | Wing and rotor thrust requirements |
| Cruise drag and rotor installation spacing | Wing | Rotor cruise trim and radius feasibility |
| Hover download/wing interference | Wing and rotor | Hover thrust sizing |
| Rotor thrust, RPM, torque, and shaft power | Rotor | Motor/ESC and battery sizing |
| Propeller mass and rotor disk diameter | Rotor | Aircraft mass ledger and wing installation |
| Motor limits, efficiencies, and usable battery energy | Propulsion | Rotor feasibility and system closure |

Detailed rotor source data and assumptions are documented in `docs/ROTOR_MODEL.md`. **Component-level solutions are not interchangeable with a converged coupled aircraft optimum.**

## Collaboration

- Please preserve existing teammates' implementations when working on one subsystem.
- Commit reproducible code, original/source-derived data with provenance, tests, and interpretable engineering results.
- Keep individual assignment reports and generated LaTeX/PDF artifacts **outside this shared project repository**.
- Run `python -m pytest -q` before merging changes.

See [CONTRIBUTING.md](CONTRIBUTING.md) for conventions. Licensed as indicated in [LICENSE](LICENSE).
