# Data provenance and units

UIUC Propeller Database, Volume 4, Dantsker, Caccamo, Deters & Selig (2022):
https://m-selig.ae.illinois.edu/props/volume-4/propDB-volume-4.html

The included measured samples are from the following public data files (with a subset of original J rows for manageable data storage; no synthetic intermediate samples):

Static:
- 13x10: https://m-selig.ae.illinois.edu/props/volume-4/data/apce_13x10_static_0559od.txt
- 16x12: https://m-selig.ae.illinois.edu/props/volume-4/data/apce_16x12_static_2167od.txt
- 18x12: https://m-selig.ae.illinois.edu/props/volume-4/data/apce_18x12_static_2201od.txt
- 20x10: https://m-selig.ae.illinois.edu/props/volume-4/data/apce_20x10_static_2210od.txt

Forward:
- 13x10 5999 and 6991 RPM: apce_13x10_0567od_5999.txt and apce_13x10_0569od_6991.txt
- 16x12 4011 and 5000 RPM: apce_16x12_2172od_4011.txt and apce_16x12_2174od_5000.txt
- 18x12 3032, 4017 and 4991 RPM: apce_18x12_2195od_3032.txt, apce_18x12_2198od_4017.txt, apce_18x12_2200od_4991.txt
- 20x10 4025 and 4506 RPM: apce_20x10_2206od_4025.txt and apce_20x10_2208od_4506.txt
All forward files under https://m-selig.ae.illinois.edu/props/volume-4/data/ .

CT=T/(rho*n^2*D^4), CP=Pshaft/(rho*n^3*D^5), advance ratio J=V/(n*D).
n is revolutions per second; RPM=60n. All propellers here are 2-bladed APC Thin Electric.

**Important qualification:** we use catalogue nominal diameters; UIUC normalizes its CT and CP using actual measured diameters, which may differ. This introduces roughly 4%-5% force/power uncertainty per 1% diameter mismatch. Values must be treated as approximate pending diameter-specific corrections. The 13x10 static curve is sourced for supplementary screening but is not included in the product-mass-known primary candidate set.

Propeller weights (not motor/nacelle weights) from manufacturer:
- 16x12E, 1.83 oz: https://www.apcprop.com/product/16x12e/
- 18x12E, 2.61 oz: https://www.apcprop.com/product/18x12e/
- 20x10E, 3.39 oz: https://www.apcprop.com/product/20x10e/
Product series revision differences may exist between historic test specimens and current listings.
