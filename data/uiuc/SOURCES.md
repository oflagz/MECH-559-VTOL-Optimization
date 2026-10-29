# UIUC Propeller Database measured benchmark (not optimized aircraft)

Selected reference: **APC Slow Flyer 10 x 7** (two-blade propeller), from
UIUC Propeller Database Volume 1 (Brandt and Selig 2011 collection).

Downloaded/transcribed factual measurement tables, retaining the source row values:

- Geometry (radial r/R, c/R, pitch degrees):
  https://m-selig.ae.illinois.edu/props/volume-1/data/apcsf_10x7_geom.txt
- Static RPM sweep (CT, CP):
  https://m-selig.ae.illinois.edu/props/volume-1/data/apcsf_10x7_static_kt0827.txt
- Forward-flight data (J, CT, CP, eta), nominal 3008 RPM:
  https://m-selig.ae.illinois.edu/props/volume-1/data/apcsf_10x7_kt0828_3008.txt
- Database description:
  https://m-selig.ae.illinois.edu/props/volume-1/propDB-volume-1.html

D = 10 inches, R = 0.127 m, N_blades = 2. The blade's **radial geometry is
measured/published**. The blade section airfoil lift/drag polars are NOT
provided. Consequently the BEMT evaluation uses the repository's **illustrative
surrogate Cl/Cd**; this is a model-to-measurement comparison/diagnostic, NOT a
validated predictive model. Different RPM/flow test runs should not be pooled
without noting their experimental conditions.

Definitions: CT = T/(rho*n**2*D**4), CP = P_shaft/(rho*n**3*D**5),
J = V/(n*D), eta = J*CT/CP, n in revolutions/second.

Observed 3008 RPM forward-flight data include negative-thrust windmilling
points at J>=0.862; we do not evaluate those points with the present
propeller-only positive-induction solver.
