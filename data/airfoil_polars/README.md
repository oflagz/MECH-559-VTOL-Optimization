# UIUC LSAT Clark-Y(B) experimental section data — NEXT sprint development

## Source and licensing

Experimental measurements by the University of Illinois at Urbana-Champaign Low-Speed Airfoil Tests (UIUC LSATs), **Summary of Low-Speed Airfoil Data, Volume 3**. The original files are:

- Lift: <https://m-selig.web.engr.illinois.edu/pd/pub/lsat/vol3/CLARKYB.LFT>
- Drag: <https://m-selig.web.engr.illinois.edu/pd/pub/lsat/vol3/CLARKYB.DRG>
- Overview and redistribution conditions: <https://m-selig.web.engr.illinois.edu/pd.html>
- Original GPL-2.0 license: <https://m-selig.web.engr.illinois.edu/pd/pub/lsat/GPL.TXT>
- UIUC LSATs Manifesto: <https://m-selig.web.engr.illinois.edu/pd/pub/lsat/MANIFEST.html>

**Notice:** Derived CSV from the UIUC LSAT wind-tunnel tests; redistributed for education and reproducible modeling under the UIUC data's GNU GPL v2 conditions, independently of the main project's software license. See `LICENSE-GPL-2.0.txt`, `COPYRIGHT-NOTICE.txt`, and the official manifesto link. Data are provided as-is, without aerodynamic accuracy warranties. Selection and interpolation of the data are modifications made October 2026.

## Processing and limitations

Run `python data/airfoil_polars/build_clarky_data.py` to reproduce the included 45-row `clarky_uiuc_experimental.csv`. Source values are separately transcribed from the ascending (non-hysteresis) clean lift runs in `.LFT` and the drag runs in `.DRG`, not invented section polars. The lift measurement at each **drag-data angle** is calculated by linear interpolation between lift measurements at the corresponding nominal Reynolds group. Values were drawn from nominal 60,000, 100,000, 200,000, 300,000 section Reynolds groups; small differences between original lift- and drag-test Re are approximated as matched. Between groups, log(Re) interpolation is used.

Input outside the measured angle/Re data is **endpoint-clamped and recorded as out-of-coverage**. Clipping is for solver continuity *only*; predictions there are not experimentally supported. Most APC SF 10×7 points at 3000 RPM fall below the minimum Re, and many high-pitch optimization cases lie outside the measured angle ranges. No attempt is made to invent E63 airfoil polars.

**Very important:** APC says slow-fly propellers blend **Eppler E63-like inboard sections with Clark-Y-like outboard sections**; the Clark-Y(B) wind-tunnel section is *not* a direct measurement of the actual APC blade. The proxy is used here for sensitivity/testing, **not validation of the optimized aircraft**.

## Reuse

This dataset comes from the UIUC Low-Speed Airfoil Tests program and is subject to the original conditions. Cite the UIUC LSAT Volume 3 report and retain the source notices and GPL when redistributing the derived measurements. This data directory is separate from the project code and its license.
