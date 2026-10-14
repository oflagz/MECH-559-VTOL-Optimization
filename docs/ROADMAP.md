# Four-week VTOL rotor roadmap (2026-10-13 to 2026-11-09)

| Week | Dates | Primary engineering goal | Concrete, testable deliverable | Dependencies |
|---|---|---|---|---|
| W1 | Oct 13-19 | Rebuild mathematically consistent rotor model and solver, with preliminary data | Tested axial BEMT; 5-variable SLSQP with hover/cruise thrust trim equalities; 2/3/4 blades; multistart; >=5 figures; scenario interface tables; reproducible repo | Provisional screening assumptions, documented distinctly from official mission |
| W2 | Oct 20-26 | Model verification against independent measured aerodynamic data | Geometry/polar ingestion; UIUC propeller measured thrust/torque/power comparison; CT/CP vs J; quantified model error; empirical assembly mass dataset and fit | Tested propeller geometry, airfoil polars and matching RPM/wind-speed data |
| W3 | Oct 27-Nov 2 | System coupling and physical mass closure | Iterated mTO = fixed + wing + battery + 4*rotor mass; validated rotor torque/power transferred to powertrain; download/drag/spacing update; minimize system energy; compare sequential vs coupled | Agreed system mission, wing and propulsion numerical interface values, motor data |
| W4 | Nov 3-9 | Credible final comparison and report | Active constraints, finite-difference gradient verification, perturbation/sensitivity analysis, derivative-free cross-check, uncertainty bands, consolidated final figures and conclusion | W2 validation + W3 interface closure |

## Priority / critical path

1. **Now:** reconcile official 5 kg payload and provisional 5 kg total mass; reconcile 40 km/20 minutes with adopted cruise airspeed; request available motor torque/power curves; request wing cruise drag, rotor spacing and hover download.
2. **Before W2 ends:** select *one* propeller family with measured data and geometry; identify blade airfoil and aerodynamic polars; obtain prop/motor/nacelle masses.
3. **Before W3 begins:** freeze dated revision of shared mission and subsystem-interface CSV; include data ownership and provenance.
4. **W4:** freeze solver versions/figures with test and validation output; explicitly distinguish *mathematical feasibility* from physical evidence.

## Acceptance criteria

**W1 complete:** tests pass; exact matched thrust in both phases; every reported design satisfies all implemented inequalities; multistart repeatability; plots regenerate without manual edits. **W2 complete only when independent measurements quantitatively compare to a matched real propeller.** **W3 complete only with mass and battery loop closed using agreed data.** **W4 complete with sensitivity/uncertainty and constraints independently checked.
