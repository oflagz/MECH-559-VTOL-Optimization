# Sprint 3: rotor aerodynamic benchmark and optimization diagnostics

## What's NEW since Sprint 2

1. **Measured data introduced:** APC Slow Flyer 10x7 physical blade geometry and wind-tunnel thrust/power coefficients from the UIUC database, from the Brandt & Selig reference cited in the proposal. Sixteen static and sixteen forward-flow data points are stored with source URLs. This 2-bladed *benchmark* is **not** our proposed 3-blade VTOL rotor.
2. **Quantified (failed) aerodynamic predictive comparison:** The existing smooth surrogate section lift/drag model gave static mean absolute percentage errors of **37.11%** in CT and **52.25%** in CP on the published APC 10x7 test. In the positive-thrust/power forward-flow subset where the model and measurement both operate as propellers (9 data points), the corresponding errors are **56.43%** and **58.54%**. At higher J the current model enters negative power/thrust before the measured data do. **Do not present as successful experimental validation.** This is evidence that blade-section polars and induced-flow treatment require improvement.
3. **MECH 559 Lecture 8 KKT and regularity / Lecture 9 active sets and SQP:** For the *surrogate screening* 3-blade optimum, the binding restrictions are minimum chord 0.018 m, maximum blade pitch 33 deg, and minimum assumed hover FM 0.58, alongside two thrust-trim equalities. The 5x5 matrix of those active/equality constraint gradients is full rank (LICQ at this numerical point), and all three estimated inequality multipliers are nonnegative. The finite-difference stationarity residual is ~1.7e-14 in scaled coordinates. Because all five degrees of freedom are constrained, the local stationarity fit is algebraically determined; **this does not certify global optimality** and does not establish the aerodynamic model is physically correct.
4. **True reoptimization after system input changes:** Unlike the previous sprint's unchanged-geometry re-trims, SLSQP was rerun at each of 9 mass/drag scenarios (5/5.5/6 kg takeoff mass, 5/7/9 N cruise drag). Rotor mission-energy predictions ranged **109.85--172.26 Wh**. At 5 kg and 7 N, the baseline is **138.10 Wh**. These are provisional *scenario* values, **not** teammate-supplied values or a closed aircraft mass/battery loop.
5. **Mesh-dependent reoptimization:** Reoptimizing at 192 instead of 96 radial elements moved the model optimum from **138.1006 to 138.2390 Wh** (+0.10%), with radius changing 0.25410 to 0.24993 m (~-1.64%). Critically, simply evaluating the old 96-element design at 192 elements did **not** satisfy exact thrust trim: hover residual -0.067% and cruise -0.235%. Optimize again when fidelity changes.

## Course concepts employed, with scope

- L8: equality and inequality KKT feasibility, complementary slackness, active-constraint regularity, Lagrange multipliers, local vs global optimality.
- L9: scaled SLSQP, active-set constraints, finite-difference Jacobians, multistart optimization.
- Sensitivity: recompute a new feasible SLSQP solution when coupled mass or drag changes; do not mistake fixed-design evaluation for coupled optimization.
- Mesh verification: distinguish numerical (grid) convergence from comparison with independent measurements.

## Next technical dependencies

- **Required from Antoine:** measured/estimated mission cruise drag, wing-download fraction and rotor mount spacing.
- **Required from Damian/system:** agreed takeoff mass including payload, efficiency and motor torque-vs-speed envelope, usable battery energy; avoid motor/prop mass double-counting.
- **Required from rotor literature/manufacturer:** real blade section Cl(alpha,Re), Cd(alpha,Re) and suitably documented motor/prop mass data before a credible numerical optimum can be called a real aircraft design.
- Reconcile proposal's 5 kg payload and 40 km mission with provisional 5 kg gross mass / 36 km screening before coupling.

## Run and reproduce

```bash
python -m pip install -r requirements.txt
python -m pytest -q
python -m scripts.sprint3_study
```

Output: `results/sprint3/` CSV, JSON, and five figures. Source data: `data/uiuc/SOURCES.md`.
