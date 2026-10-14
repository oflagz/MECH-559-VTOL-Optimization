# Sprint 2 — Rotor subsystem only (Ashish)

## Sprint goal
Rebuild the VTOL rotor optimizer with a consistent hover/cruise performance model, implement SLSQP for the agreed five-variable formulation and obtain defensible preliminary **model** results with diagnostic plots and numerical checks.

## Planned/completed tasks

- [x] Replace the earlier independently mixed power and thrust calculations with one BEMT thrust/torque/power model evaluated at matched geometry, pitch and RPM.
- [x] Treat blade pitch at 75% radius, theta75, as an explicit continuous variable; enumerate N_blades=2/3/4 with N_rotors=4 fixed.
- [x] Implement SciPy SLSQP with normalized bounds, two exact thrust-trim equalities, power/torque/Mach/clearance/FM and reserve thrust inequalities, and three initial points per blade count.
- [x] Reproduce the proposal's *analytical screening* calculation independently.
- [x] Generate and test numerical optimum tables, mission-energy plots, feasible radius/pitch sweep, trim curve, wing mass/drag input sensitivity and discretization study.
- [x] Define the required numerical interfaces with wing and propulsion models, prioritizing takeoff mass, cruise drag, wing download, mount spacing, motor torque/power/efficiency and battery information.
- [ ] Validate real propeller performance with measured data (W2; not completed).
- [ ] Add calibrated rotor/nacelle/motor mass and aircraft mass closure (W2-3; not completed).

## Preliminary numerical results (surrogate airfoil; 96 blade elements; screening mission only)

| Blades | Energy, Wh | Hover electrical W total | Cruise electrical W total | R, m | c, m | Pitch75, deg | omega_h, rad/s | omega_c, rad/s |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 2 | 152.76 | 750.65 | 255.48 | 0.1800 | 0.02392 | 27.38 | 466.84 | 363.79 |
| 3 | **138.10** | 531.74 | 240.75 | 0.2541 | 0.01800 | 33.00 | 254.68 | 197.89 |
| 4 | 146.48 | 516.72 | 258.52 | 0.2488 | 0.01800 | 33.00 | 229.40 | 198.45 |

At the 3-blade screening optimum the expected thrust is 12.876 N/rotor in hover and 1.750 N/rotor in cruise, matching the scenario's required thrust; hover electrical energy 17.72 Wh and cruise 120.38 Wh. The cruise portion contributes approximately 87.2% of total rotor electrical mission energy because the scenario spends 30 minutes cruising but only 2 minutes hovering. These are provisional numerical model outputs, NOT experimental/design certification.

## Physical interpretation

1. **Cruise dominates the objective** despite hover using more instantaneous power, because cruise lasts 15 times longer. Thus a lowest-hover-power-only design is not mission-optimal.
2. **3 blades won in this surrogate screening**, not because 3 blades are universally best but because the trade-off of lower hover demand and cruise blade/profile losses worked out that way with the chosen geometry and loss assumptions. Change the airfoil, bounds, mission or mass model and the ranking can change.
3. **Several bounds are active**: the 3-blade optimum uses the minimum chord 0.018 m and maximum pitch 33 degrees; its hover figure-of-merit exactly reaches the assumed 0.58 minimum. This signals strong sensitivity to guessed aerodynamic polars and variable limits. It is not evidence that such a high-pitch, narrow blade is a recommended manufactured rotor.
4. **Analytical hover cross-check** at original proposal point recovers 12.876 N/rotor, 442.57 W total electrical input and FM 0.708; this is a separate idealized calculation, not the 3-blade BEMT optimum. Two models cannot be substituted for one another without calibration.
5. **Numerical verification:** all nine SLSQP starts converged to constraint-feasible solutions with energy repeatability within ~1e-4 Wh per blade count. Separate 96-to-192-element BEMT evaluations assess mesh sensitivity. Because the optimizer moves when changing discretization, the numerical geometry should not yet be called converged to the accuracy of the printed decimal places.
6. **Coupling sensitivity** is tabulated by changing supplied aircraft mass and cruise drag and re-trimming RPM at unchanged rotor geometry. Higher mass increases hover demand; higher wing drag increases cruise energy. Scenarios without a trim inside the RPM bounds are correctly marked infeasible, not assigned zero power.

## Impediments/lessons learned

The previous hover-only parametric plots were insufficient to establish an optimal rotor because neither required thrust nor cruise trim was enforced. The revised formulation avoids that error. The dominant blocker is the absence of measured low-Reynolds-number blade polars, actual blade geometry, realistic rotor/motor mass, and an agreed mission and team-provided drag/motor parameters. Numerical feasibility here is not physical validation.

## Next sprint

Select one matched propeller dataset from UIUC/manufacturer sources; implement actual blade chord/twist sections and airfoil polar interpolation; quantify thrust/power errors against independent measurements; fit rotor/motor/nacelle assembly mass and evaluate numerical uncertainty. Confirm team-level takeoff mass, drag, power limits and mission prior to mass-energy coupling.
