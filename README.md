# MECH 559 — VTOL Tiltrotor Optimization

McGill University · Group 5 · Fall 2026

## Subsystems

- Wing: `wing_model.py`, `wing_sprint3.py` (Antoine)
- Rotor: `rotor_design/`, `analysis/` (Ashish)
- Electric propulsion: `propulsion/` (Damian)

Mission: 5 kg payload, four rotors, 40 km cruise in 1,200 s (33.33 m/s), 300 s hover.

## Rotor analysis

The individual-report rotor model uses signed-induction BEMT calibrated against UIUC APC Slow Flyer 10×7 data. It optimizes radius, chord and twist with thrust-trimmed hover/cruise speeds using SLSQP, enumerates blade counts, and evaluates COBYLA, Pareto, KKT, mesh-convergence and interval-sensitivity checks.

```bash
python -m pip install -r requirements.txt
python -m pytest -q
python -m analysis.run_all
```

Outputs: `results/figures/` and `results/tables/`. Experiment provenance: `data/uiuc/SOURCES.md`. Engineering assumptions: `docs/ROTOR_ASSUMPTIONS.md`.

**Scope:** Component-level rotor optimization only. The model is not a validated propeller design or a converged wing–rotor–propulsion MDO solution. The selected motor/battery screening configuration is not compatible with the rotor's hover-voltage requirement.
