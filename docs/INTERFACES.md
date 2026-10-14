# Interface requirements: rotor first, enabling other two subsystems

## Inputs required *now* (owner / units / why)

| Parameter | Owner | Units | Priority | Effect |
|---|---|---|---|---|
| Agreed gross takeoff mass including payload | Team/Damian | kg | P0 | Hover thrust, rotor radius and total mass closure |
| Agreed cruise speed/distance/hover duration | Team | m/s, m, s | P0 | Cruise operating point, mission energy |
| Cruise drag D_c at selected cruise speed | Antoine | N | P0 | Exact cruise propeller thrust, motor electrical energy |
| Wing hover-download fraction and rotor/wing interference | Antoine | dimensionless | P0 | Hover thrust and wing/rotor couplings |
| Rotor-to-rotor/nacelle mount spacing | Antoine | m | P0 | Maximum rotor diameter and clearance |
| Motor continuous/peak torque vs speed | Damian | Nm vs rad/s | P0 | Power/torque feasibility and reserve headroom |
| Motor/ESC efficiency vs torque/RPM, total electrical limit | Damian | efficiency, W | P0 | Electrical mission energy and feasibility |
| Usable battery energy, pack mass, reserve policy | Damian | Wh, kg | P1 | Actual system constraint, mass loop |
| Blade r/R chord, pitch/twist distributions and airfoil identity | Rotor | m, deg | P1 | Actual BEMT geometry |
| Experimental propeller CT/CP curves with matched operating conditions | Rotor | nondimensional | P1 | Predictive validation against measurement |
| Rotor+hub+motor+nacelle/tilt mechanism masses | Rotor/Damian | kg | P1 | Radius/motor mass penalty in system objective |

## Outputs generated for the other subsystems

- `results/tables/optimization_by_blade_count.csv` includes chosen blade count, geometry, pitch, RPM at hover/cruise, total electrical power per phase, per-rotor shaft torque and thrust, total rotor disk area, disk loading, efficiencies, and energy in Wh.
- `results/tables/wing_mass_drag_interface.csv` shows the sensitivity of energy/torque to changed takeoff mass and wing drag at fixed rotor geometry, re-trimming RPM to matched thrust; it is not a new fully coupled optimum.
- `results/tables/mesh_convergence.csv` captures numerical BEMT sensitivity.
- `results/tables/multistart_audit.csv` records SLSQP repeatability and feasibility.
- `results/tables/study_manifest.json` includes all screening inputs and optimization results.

**Important:** Rotor assembly mass is UNKNOWN, not silently assumed zero. Rotor mass and battery closure remain required interface outputs for Weeks 2–3.

## Reference sources

UIUC Propeller Data Site: https://m-selig.ae.illinois.edu/props/propDB.html. Matches the Brandt and Selig 2011 low-Reynolds-number propeller study cited in the group's proposal. The UIUC database supplies measured thrust/power coefficient performance; matched blade-geometry and polar data must also be obtained to validate a specific propeller.
