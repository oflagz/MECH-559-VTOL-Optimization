# NEXT SPRINT development — aerodynamic accuracy audit (NOT Sprint 3 report)

## Status as of 28 October 2026

This is preparatory work for next week's rotor effort. **No Sprint 3 report files were changed.** The previous 138.10 Wh SLSQP figure remains only a legacy screening result and is not experimentally validated.

## 1. Replace Re-independent airfoil assumptions — COMPLETED AS A LIMITED PROXY

Implemented `rotoropt/polars.py`: a selectable `Mission(polar_model='clarky_uiuc')` uses 45 derived pairs from UIUC **measured** Clark-Y(B) airfoil lift/drag tests, organized into nominal Re = 60k, 100k, 200k, 300k groups, with linear angle interpolation and log-Re interpolation. Outside the measured domain the values are endpoint-clamped for continuity but the model reports loading-weighted coverage. This is *not* a true APC blade polar; the inboard Eppler E63 section and its blend remain proprietary/unresolved.

## 2. Verify blade geometry and conventions — COMPLETED FOR REFERENCE PROPELLER

Published APC Slow Flyer **10×7** benchmark from UIUC uses D=0.254 m (R=0.127 m), two blades, r/R=0.15..1, chord c/R, and blade geometric pitch in degrees. At r/R=0.75 the tabulated pitch is 14.38°. `tests/test_measured_polars.py` verifies these values and ascending radial stations. The benchmark geometry is distinct from the optimization problem's uniform-chord, constant-pitch abstract rotor. Do not substitute one geometry for the other.

## 3. Check induction/loss contributions — NUMERICAL SENSITIVITY COMPLETED

`results/next_sprint/induction_and_section_sensitivity.csv` compares tip/hub-loss on/off, angular swirl induction on/off, and a deliberately uncertain inboard section substitution at both 3029 and 5987 RPM. **These are numerical model variations**, not experimental proof that a particular loss model is correct. No hardware calibration was performed.

## 4. Recompare against independent propeller data — COMPLETED, PHYSICAL VALIDATION FAILED

Full 16 static RPM and 16 forward-flight J cases were run, without fitting any section coefficient to those measurements. Reported MAPE percentages are absolute percentage differences, not measurement uncertainty. For forward flight, errors on points where both model and experiment predict positive thrust/power are reported separately; the false-sign count is always reported, so unsupported points are not hidden.

| Comparison | Old surrogate | Measured Clark-Y section proxy | Interpretation |
|---|---:|---:|---|
| Static C_T MAPE | 37.11% | 25.67% | Improved numerically, **not validated** |
| Static C_P MAPE | 52.25% | 39.67% | Improved numerically, **not validated** |
| Forward C_T MAPE, both positive | 56.43% (9 cases) | 50.24% (10 cases) | Different sample counts; not a controlled pooled comparison |
| Forward C_P MAPE, both positive | 58.54% (9 cases) | 45.71% (10 cases) | Different sample counts; not directly paired |

**Major limitation:** at approximately 3029 RPM, the APC 10×7 blade-section Re range is about 5.1k–50.4k, while the source Clark-Y data start at 60k. Therefore the **weighted polar coverage is 0% at this condition**, and even at 5987 RPM it reaches only about 61%. A lower numerical error under extensive clipping is **not** meaningful experimental validation.

## 5. Rerun SLSQP — COMPLETED, ALL NEW CANDIDATES UNQUALIFIED

The same five-variable SLSQP + 2/3/4-blade enumeration was rerun using the selectable Clark-Y proxy. All three cases satisfy their *programmed* trim and inequality constraints, but measured airfoil coverage is insufficient in hover and cruise, so **no design passes physical-coverage qualification**.

| Blades | Legacy surrogate (Wh) | Clark-Y proxy (Wh) | Verified optimum? |
|---:|---:|---:|---|
| 2 | 152.76 | 122.47 | **No** |
| 3 | 138.10 | 127.84 | **No** |
| 4 | 146.48 | 133.18 | **No** |

The ordering changes, illustrating sensitivity to the airfoil assumption. **Do not report the proxy energy numbers as optimized real hardware or conclude that 2 blades are better.** All three proxy optima use a 170 rad/s cruise RPM lower bound and a 33° pitch upper bound. Additionally, no rotor mass or full aircraft mass closure is included. Candidate screening is not completion of the project.

## What is genuinely needed next

- Obtain low-Re polars **covering roughly 5k–100k** for an APC-like E63/Clark-Y blend, or select a benchmark propeller operating inside available measured polar limits. If using XFOIL/NeuralFoil, label results computationally generated, not experimentally measured.
- Prefer actual manufactured blade section geometry and measured static/dynamic performance, with validation points reserved from fitting.
- Implement a documented data-domain constraint or uncertainty penalty so SLSQP cannot exploit clipped stall/low-Re polars.
- Add rotor/motor assembly mass and teammates' confirmed total mass, wing drag, download, spacing, torque limits and efficiencies. Reconcile the proposal's 5 kg payload and provisional 5 kg gross takeoff mass before final optimization.
- Repeat optimization after the aerodynamic model covers operating conditions and evaluate KKT, multistart, and mesh sensitivity, with a separate experimental accuracy report.

## Reproduction

```bash
python -m pip install -r requirements.txt
python data/airfoil_polars/build_clarky_data.py
python -m pytest -q
python -m scripts.next_sprint_aerodynamics
```

Results (CSV, plots, manifest) are saved in `results/next_sprint`. UIUC measured airfoil data license and citations are documented under `data/airfoil_polars`.
