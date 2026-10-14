# Next sprint: aerodynamic model upgrade

## Implemented in this iteration
- Annular blade-element calculations with per-annulus induced axial velocity.
- Swirl velocity/torque coupling and relaxed fixed-point iterations.
- Numerical convergence flag, iteration count and local Reynolds range.
- CSV-based Cl/Cd lookup with angle-of-attack and log-Reynolds interpolation.
- Basic numerical smoke test for BEMT.

## Not yet validated
The existing optimizer still defaults to the analytical surrogate polar.
The annular model has no Prandtl tip/hub loss, high-induction corrections,
compressibility model, reverse-flow handling, or empirical calibration.
Its output is **not validated design evidence**.

## Reference-data validation checklist
1. Select a specific propeller and obtain geometry (radius, chord and twist).
2. Obtain airfoil polar grids at representative Reynolds numbers.
3. Collect published thrust and torque/power curves versus RPM and advance ratio.
4. Evaluate matched operating conditions with the same geometry and density.
5. Report signed and absolute errors and convergence; reject extrapolated points.
6. Only then rerun constrained optimization and assess sensitivity.

The user must supply or approve an actual propeller reference dataset before
a defensible quantitative validation can be reported.
