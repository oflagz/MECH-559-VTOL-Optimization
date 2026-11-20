# Contributing

This repository is shared by the wing, rotor, and electric-propulsion project contributors.

## Before a change

1. Create a short-lived feature branch for a subsystem-specific change.
2. Keep input data provenance and assumptions with the model that uses them.
3. Prefer small, reviewable commits; do not rewrite another contributor's model without coordinating with its owner.
4. Run the relevant tests and document changes to shared quantities (mass, power, velocity, drag, energy).

## Reproducibility

- Python dependencies: `requirements.txt`.
- Rotor tests: `python -m pytest -q`.
- Rotor figures/tables: `python -m scripts.run_rotor_study`.
- Technical results belong under `results/<subsystem>/`, with source datasets under `data/<subsystem>/`.
- Record units and distinguish measured inputs from estimated/assumed parameters.

## Files that do not belong in the shared repository

Individual project reports, compiled PDFs/LaTeX documents, temporary notebooks, Python caches, rendered logs, and redundant output copies. Keep individual submissions separately (e.g., in a private course-submission folder). Do not add any personally identifying student records.
