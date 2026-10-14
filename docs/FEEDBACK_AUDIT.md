# Professor feedback and Sprint 1 correction audit

The earlier professor comments shared with the team requested explicit literature/prior work, common mission/units, declared variables/outputs and subsystem coupling, rotor/wing/propulsion interface owners, physical trade-offs (including rotor mass), integer decision handling, and numerical verification. Sprint 1 records the more specific correction to optimize theta75 instead of blade twist.

| Requirement | This week's implementation | Remaining obligation |
|---|---|---|
| Fixed four-rotor architecture | Mission.n_rotors=4 | None for topology |
| Use theta75 as continuous design variable | SLSQP x[4]=theta75; constant twist is fixed at 0 | Verify selected propeller family physically supports stated pitch |
| Blade-count integer handling | Enumerated 2,3,4 independently | Catalogue-driven options in W2 |
| Optimize mission energy, not hover alone | Exact hover+cruise shaft/electrical power objective | System total energy with mass feedback in W3 |
| Thrust capability at same operating point as power | BEMT thrust equalities for both phases; RPM reserve capability inequality | Cross-check with independent measurements in W2 |
| Physical/geometric limits | Scaled power, torque, Mach, clearance, FM/headroom inequalities | Agreed motor and wing component limits |
| Reproducible feasible-point calculation | Analytical hover regression; model trim, constraints, output CSV/plots | Validate reality of polar and geometry |
| Multistart/active constraints/sensitivity | 3 starts per blade-count; binding constraints, radius-pitch sweep, mass/drag interface, mesh convergence | Independent optimizer verification W4 |
| Rotor mass feedback | Interface specified, mass model omitted until fit | Empirical rotor+motor+nacelle fit W2, aircraft mass closure W3 |
| Common mission and data ownership | Explicit inconsistent mission cases, ownership table | Team to confirm common mission |

This audit refers to the comments previously supplied by the user and to Sprint 1; it is not a new professor evaluation.
