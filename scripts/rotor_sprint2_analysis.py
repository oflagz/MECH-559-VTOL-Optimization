"""Preliminary analytical hover screening and sensitivity plots.

Run from repository root: python -m scripts.rotor_sprint2_analysis
Outputs CSV/PNG under results/. Not a constrained design optimum.
"""
from pathlib import Path
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from src.rotor.parameters import Params
from src.rotor.hover import hover_model

OUT=Path("results")/"sprint2"
OUT.mkdir(parents=True,exist_ok=True)
p=Params()
radii=np.linspace(0.15,0.30,121)
omegas=np.linspace(150.,650.,121)
def model(r,omega=300.):
    h=hover_model(r,0.03,omega,3,p)
    return h["Pi"],h["P0"],h["P_sh"],h["FM"],p.Nr*h["P_el"]
rows=[]
for r in radii:
    pi,p0,shaft,fm,power=model(r)
    rows.append([r,pi,p0,shaft,fm,power,power*p.th/3600])
with (OUT/"radius_sweep.csv").open("w",newline="") as f:
    w=csv.writer(f)
    w.writerow(["radius_m","induced_w_per_rotor","profile_w_per_rotor",
                "shaft_w_per_rotor","figure_of_merit","hover_electric_w",
                "hover_energy_wh"])
    w.writerows(rows)
def finish(filename):
    plt.grid(alpha=.25);plt.tight_layout()
    plt.savefig(OUT/filename,dpi=180);plt.close()
a=np.array(rows)
plt.figure(figsize=(7,4.5))
plt.plot(radii,a[:,5],label="Total hover electrical power")
plt.axvline(.25,ls="--",color="gray",label="Proposal screening radius")
plt.xlabel("Rotor radius (m)");plt.ylabel("Electrical power (W)")
plt.title("Hover power vs radius");plt.legend();finish("01_hover_power_vs_radius.png")
plt.figure(figsize=(7,4.5))
for col,label in [(1,"Ideal induced"),(2,"Profile"),(3,"Shaft")]:
    plt.plot(radii,a[:,col],label=label)
plt.xlabel("Rotor radius (m)");plt.ylabel("Power per rotor (W)")
plt.title("Hover power components");plt.legend();finish("02_hover_power_components.png")
plt.figure(figsize=(7,4.5))
plt.plot(radii,a[:,4])
plt.xlabel("Rotor radius (m)");plt.ylabel("Figure of merit")
plt.title("Hover figure of merit");finish("03_figure_of_merit.png")
grid=np.array([[model(r,o)[4] for r in radii] for o in omegas])
with (OUT/"radius_rpm_grid.csv").open("w",newline="") as f:
    w=csv.writer(f);w.writerow(["omega_rad_s","radius_m","hover_electric_w"])
    w.writerows((o,r,grid[i,j]) for i,o in enumerate(omegas)
                for j,r in enumerate(radii))
plt.figure(figsize=(7,4.5))
cs=plt.contourf(radii,omegas,grid,levels=25)
plt.colorbar(cs,label="Electrical power (W)")
plt.scatter([.25],[300],color="red",label="Proposal screening")
plt.xlabel("Rotor radius (m)");plt.ylabel("Angular speed (rad/s)")
plt.title("Hover power sensitivity");plt.legend();finish("04_hover_power_contour.png")
baseline=model(.25)
best=a[np.argmin(a[:,5])]
print(f"Screening: thrust/rotor={(1+p.delta_dl)*p.mTO*p.g/p.Nr:.3f} N, "
      f"shaft/rotor={baseline[2]:.2f} W, FM={baseline[3]:.3f}, "
      f"hover electrical={baseline[4]:.2f} W")
print(f"Fixed-RPM radius sweep minimum: R={best[0]:.3f} m, "
      f"power={best[5]:.2f} W (not a constrained optimum)")
