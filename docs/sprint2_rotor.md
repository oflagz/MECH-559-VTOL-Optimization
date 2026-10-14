# Sprint 2: VTOL rotors (Ashish)

## Goal
Assign and implement the rotor subsystem numerical analysis and constrained
optimization algorithm. Align the design variables with Sprint 1's agreed
theta75 pitch variable.

## Completed code
- Python/SciPy SLSQP multistart, blade-count enumeration and positive-margin
  constraints.
- Corrected design vector to R, c, omega_hover, omega_cruise, theta75.
- Annular axial/swirl induction solver with explicit convergence flags.
- Reynolds-aware airfoil polar CSV interpolation and measurement comparison
  script (MAE for thrust and power).
- CI workflow for pytest and analytical hover screening reproduction.

## Important qualifications
The blade-element solver still uses surrogate airfoil data by default.
No experimental propeller dataset has yet been committed or evaluated.
Numerical optimization results must not be presented as validated physical
designs. Constant-chord/zero-twist blade geometry is provisional. Mission
parameters still require reconciliation with the team.

## Next technical milestone
Choose a real propeller with blade geometry, polar grids and independent
thrust/power measurements; run matched-condition validation and convergence
sweeps. Integrate the validated rotor model into the system mass/power loop.
