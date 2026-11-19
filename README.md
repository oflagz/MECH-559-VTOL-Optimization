# MECH 559 — Clean rotor subsystem (Group 5)

This directory contains the **2026 individual rotor engineering study** rebuilt from scratch after auditing the earlier BEMT optimization. It is based on published experimental UIUC propeller performance and has explicit boundaries on extrapolation. **Do not reuse the retired 138.10 Wh surrogate-model optimum.**

## Run everything

```bash
python -m pip install -r requirements.txt
python -m pytest -q
python data/make_source_tables.py
python -m scripts.run_study
cd . && pdflatex -interaction=nonstopmode -halt-on-error -output-directory=docs docs/individual_report.tex
```

`data/make_source_tables.py` reproduces the source-transcribed static/dynamic CSVs.
The primary engineering results are in `results/tables/study_summary.json`; the five plots and derived tables are in `results/figures` and `results/tables`. The PDF is `docs/individual_report.pdf`.

## Folder layout

- `rotor/catalog.py` — experimental propeller lookup; PCHIP-in-J and bracketed-RPM interpolation; no extrapolation.
- `rotor/optimization.py` — measured-data catalogue enumeration, thrust trim, constraints, SLSQP and independent brentq check.
- `rotor/physics.py` — corrected signed annular momentum equations; explicitly flags unsupported BEMT regimes. **Diagnostic, not physical design selector.**
- `rotor/verification.py` — held-out measured-curve interpolation checks.
- `data/` — referenced public experiments, APC propeller masses, and source transcription script.
- `scripts/run_study.py` — figure/table generation and all sensitivity sweeps.
- `docs/individual_report.tex` — exact five-section individual course report; compiled to PDF.
- `tests/` — reproducibility, sign handling, measured-domain constraints and sanity checks.

## Technical status and limits

- 13 tests pass locally.
- Main alternative (not official mission): 5 kg payload included in an assumed **9.5 kg base aircraft excluding propellers**, 40 km at 20 m/s, 5 minutes hover: APC 18x12E (two-blade) 3890 RPM hover / 3500 RPM cruise, **278.3 Wh rotor-only electrical**. The 40-km cruise lasts **33.3 minutes**, not the official 20 minutes.
- Official 40 km in 20 min requires **33.3 m/s**: none of the four measured catalogue propellers has a verified matching forward-flight trim in the covered RPM/J domain. This is not proof of infeasibility of every possible hardware design.
- Catalogue optimization has discrete propeller selection and two continuous RPMs. Original five-variable geometry BEMT model retired due inconsistent negative loading. No independent claim of optimal diameter/chord/pitch in continuous space is made.
- The manufacturer's 2026 propeller mass specifications and historic UIUC test specimens may not be identical. Source CT/CP use true measured test diameter, but available screening uses nominal diameter: +/-1% diameter sensitivity recorded.
- Assumed power/torque and cruise-drag interfaces must be approved by teammates; no full coupled MDO or aircraft certification is claimed.

## Sources

See `data/SOURCES.md` for each measured UIUC curve and each APC product weight. UIUC Database Volume 4: https://m-selig.ae.illinois.edu/props/volume-4/propDB-volume-4.html .
