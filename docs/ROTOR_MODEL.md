# VTOL rotor subsystem — modeling notes

## Scope

The rotor model supports a four-rotor tiltrotor concept by estimating the shaft torque, shaft power, hover thrust capability, forward-flight thrust, and mission energy of discrete commercial propellers. It evaluates hover and cruise operating RPM for a fixed candidate propeller and checks availability of experimentally measured operation points.

## Governing relationships

For propeller diameter \(D\), rotation rate \(n\) (revolutions per second), density \(\rho\), and forward velocity \(V\), the lookup uses the source-defined nondimensional quantities:

\[
J=\frac{V}{nD},\quad T=C_T\rho n^2D^4,\quad
P_{shaft}=C_P\rho n^3D^5,\quad
Q=\frac{P_{shaft}}{2\pi n}.
\]

Hover and cruise thrust must be trimmed to the mission-required values; electrical energy over the rotor flight phases is

\[
E_r=\frac{N_r}{\eta_{drive}}\frac{P_{h,shaft}t_h+P_{c,shaft}(L/V_c)}{3600}\quad [\mathrm{Wh}].
\]

The model assumes one constant drive efficiency, treated as provisional until motor/ESC maps are provided by propulsion.

## Optimization and numerical checks

- Discrete search over published APC propeller catalogue candidates.
- For each candidate, bracketed thrust roots determine hover and cruise RPM inside the published data region.
- SLSQP independently checks the same two-variable RPM trim as an algorithmic cross-check.
- Limits include mount clearance, tip Mach, nominal torque, peak electrical power, and 15% hover thrust capability headroom.
- Source measurements are interpolated by shape-preserving cubic polynomials within observed RPM or advance-ratio bounds; no extrapolation is permitted.
- Held-out interpolation tests quantify numerical interpolation behavior; they are **not validation for untested blade geometries**.

## Engineering boundaries

1. **Official vs alternative mission:** The official requirement is 40 km in 20 minutes (33.33 m/s) and 5 minutes of hover with a 5 kg payload. The lower-speed demonstration is 40 km at 20 m/s and does not meet the official cruise-time requirement.
2. **Measured data envelope:** Failure to identify a measured-domain trim does not establish an aeronautical impossibility. It means that the current experimental catalogue is insufficient to substantiate that operating point.
3. **Continuous geometry:** A prior BEMT-based continuous-geometry optimum was withdrawn after signed-loading consistency issues. The active method selects from a discrete catalogue. It does not optimize independent chord, radius, and pitch in continuous space.
4. **Mass and power:** Propeller masses are included in the rotor catalogue; teammate-approved complete motor/nacelle/battery mass feedback is not closed. Drive efficiency, rotor mounting spacing, drag, download, and other limits are provisional.
5. **Test diameters:** Published coefficients use actual test propeller diameter, whereas the initial dimensional estimates use nominal catalogue diameter; the effect is investigated in `results/rotor/tables/nominal_diameter_sensitivity.csv`.

## Study outputs

Run `python -m scripts.run_rotor_study` to regenerate the following:

- `results/rotor/tables/slower_cruise_catalog.csv`: alternatives in the slower-cruise screening mission.
- `results/rotor/tables/official_mission_catalog.csv`: status of each candidate at the official speed.
- `results/rotor/tables/mass_drag_sensitivity.csv`: parameter changes and re-trimmed energy.
- `results/rotor/tables/heldout_interpolation.csv`: interpolation checks.
- `results/rotor/tables/study_summary.json`: machine-readable scenario assumptions and results.
- `results/rotor/figures/`: plots generated from these studies.

Data provenance: [`data/rotor/SOURCES.md`](../data/rotor/SOURCES.md).
