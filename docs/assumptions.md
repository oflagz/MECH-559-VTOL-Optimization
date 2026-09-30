# Assumptions and model status

## Source of truth
The rotor formulation follows the MECH 559 project proposal. The repository keeps
proposal-defined equations separate from provisional numerical assumptions.

## Mission reconciliation required
The rotor screening section uses a provisional 5 kg takeoff mass, 120 s hover,
36 km cruise at 20 m/s, 7 N cruise drag, and 1 kW hover electrical allowance.
The project-level mission description uses a 5 kg payload, 40 km sortie and
5 minutes total hover. These values must be reconciled by the team before the
coupled system optimization is treated as final.

## Implemented
- Hover momentum theory and profile-power model.
- Proposal rotor energy objective.
- Constraints g1-g10 in positive-margin form for SciPy.
- Enumeration of 2, 3 and 4 blades with SLSQP multistart.
- A preliminary axial blade-element/momentum calculation.

## Provisional / not final
- theta75 = 20 deg.
- Tip-Mach, torque, cruise-power, clearance and mount-spacing numerical limits.
- Smooth analytical Cl/Cd surrogate.
- Uniform induced velocity and no swirl in BEMT.
- Continuous design-variable bounds.

## Still required
1. Select airfoil/propeller family and import Reynolds-aware Cl/Cd polar data.
2. Upgrade BEMT to annular axial/swirl induction with convergence checks.
3. Validate thrust/power against published or manufacturer propeller data.
4. Fit rotor/motor mass models from candidate hardware.
5. Couple takeoff mass, wing drag/mount spacing and battery usable energy.
6. Add system reserve/transition/avionics energy and g11/g12 where applicable.
