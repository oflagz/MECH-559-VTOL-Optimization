# Explicit engineering assumptions and limitations

## Different missions in the team's proposal

- Project abstract: 5 kg **payload**, 40 km route, 300 seconds total hover, 20 minutes wing cruise.
- Rotor mathematical screening case: 5 kg **gross takeoff mass**, 36 km at 20 m/s = 1800 seconds (30 minutes) cruise, 120 seconds total hover.
- The abstract's 40 km / 20 minutes implies 33.33 m/s if all distance is in cruise, not 20 m/s. A future common-mission decision is necessary. The screening case is NOT a substitute for the official mission.

## Unvalidated assumptions

- Four tiltrotors, constant chord, fixed zero blade twist, continuously optimized pitch at 75% radius; blade count enumerated 2, 3, 4. Pitch does not change between hover and cruise, matching a fixed-pitch propeller concept.
- Illustrative smooth Cl(alpha), Cd(alpha) surrogate with no experimental airfoil Reynolds dependence. Re is reported but unused in surrogate coefficients. Prandtl corrections and angular induction are implemented. Do not describe this as validated BEMT.
- Radius 0.18–0.30 m; chord 0.018–0.052 m; hover omega 170–650 rad/s; cruise omega 170–750 rad/s; theta75 10–33 deg. Bounds chosen for numerical screening only.
- Download 5%, thrust margin 15%, 0.65 m rotor mount spacing, 0.05 m clearance, 0.58 min hover figure of merit, 1 kW total electrical power budget, 1.5 Nm per motor torque, max tip Mach .65, motor+electrical efficiency .85; all provisional and NOT certified hardware limits.
- Extra thrust capability evaluated at 1.10 times normal hover RPM (>=15% additional required thrust), subject to same torque/power/tip limits at this reserve speed. This is an assumed screening operating headroom model, not actuator dynamics.
- Steady level cruise: each rotor exactly balances Dc/4; hover each rotor exactly balances (1+download)*mTO*g/4. Analytic hover momentum is an **independent screening check**, not the optimization objective.
- All mission energy is rotor/motor electrical use; transitions, avionics, reserve and mass-feedback absent. Usable battery energy constraint is deliberately deferred to multidisciplinary system optimization.
- No rotor blade/nacelle/motor assembly mass curve has been calibrated. Fixed-mass local energy optimum must not be called a total-aircraft mass-energy optimum.

## Interpretation rules

- Boundary optima in pitch/chord/FM are **evidence that the answer is sensitive to provisional bounds and polar assumptions**. Do not select hardware from them.
- A feasible SLSQP result is feasible **only within the stated screening model**.
- No measured experimental comparison has been made in W1. Week 2 must establish the predictive error using matched geometry, aerodynamic polars and published test data.
- Changing blade discretization from 24 to 96 elements changes the optimal 3-blade radius materially. The study therefore optimizes at 96 elements; compare to 192-element evaluations in the convergence table. The full convergence of the *optimized design* is not certified, only spot-checks at its operating point.
