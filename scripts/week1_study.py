"""Run from project root: python -m scripts.week1_study

Writes numerical results, plots, a JSON audit, and a report-ready summary.
All findings are MODEL RESULTS with non-validated surrogate airfoil input.
"""
from __future__ import annotations
import csv
import json
from dataclasses import asdict, replace
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from rotoropt.model import Mission, predict, hover_momentum_check
from rotoropt.optimize import analyze, solve, _initial_speed, LOW, HIGH

ROOT=Path(__file__).resolve().parents[1]
TABLE=ROOT/'results'/'tables'
FIG=ROOT/'results'/'figures'
TABLE.mkdir(parents=True,exist_ok=True)
FIG.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'figure.figsize':(7.2,4.6),'figure.dpi':115,
                     'axes.spines.top':False,'axes.spines.right':False,
                     'font.size':10,'savefig.dpi':170})

def write_csv(name,rows):
    if not rows:return
    with (TABLE/name).open('w',newline='',encoding='utf-8') as file:
        w=csv.DictWriter(file,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)

def savefig(name):
    plt.tight_layout();plt.savefig(FIG/name,facecolor='white');plt.close()

def trim_geometry(R,c,theta,nb,m):
    oh=_initial_speed(R,c,theta,nb,m,0.,m.hover_thrust_n)
    oc=_initial_speed(R,c,theta,nb,m,m.cruise_speed_mps,m.cruise_thrust_n)
    if oh is None or oc is None or not(LOW[2]<=oh<=HIGH[2] and LOW[3]<=oc<=HIGH[3]):
        return None
    return analyze([R,c,oh,oc,theta],nb,m)

def main():
    m=Mission()
    analytical=hover_momentum_check(.25,.030,300.,3,m)
    with (TABLE/'screening_analytical.json').open('w') as f:
        json.dump(analytical,f,indent=2)
    results=[];optim_runs=[]
    for nb in (2,3,4):
        best,runs,starts=solve(m,n_blades=nb,max_starts=3)
        for i,run in enumerate(runs):
            optim_runs.append(dict(blades=nb,start=i+1,success=run['success'],
                                  feasible=run['feasible'],energy_wh=run['analysis']['energy_wh'],
                                  iterations=run['iterations'],message=run['message']))
        if best:
            a=best['analysis'];x=a['geometry']
            results.append(dict(blades=nb,energy_wh=a['energy_wh'],
                hover_energy_wh=m.n_rotors*a['hover']['electrical_power_w']*m.hover_seconds/3600,
                cruise_energy_wh=m.n_rotors*a['cruise']['electrical_power_w']*m.cruise_seconds/3600,
                hover_power_total_w=m.n_rotors*a['hover']['electrical_power_w'],
                cruise_power_total_w=m.n_rotors*a['cruise']['electrical_power_w'],
                hover_thrust_each_n=a['hover']['thrust_n'],
                cruise_thrust_each_n=a['cruise']['thrust_n'],
                hover_torque_each_nm=a['hover']['torque_nm'],
                cruise_torque_each_nm=a['cruise']['torque_nm'],
                hover_fm=a['hover']['figure_of_merit'],
                cruise_prop_efficiency=a['cruise']['prop_efficiency'],
                disk_area_total_m2=a['disk_area_total_m2'],
                disk_loading_n_m2=a['disk_loading_n_m2'],
                min_margin=min(a['margins'].values()),
                max_abs_trim_error=float(max(abs(a['trim_residual']))),
                all_converged=a['margins']['bemt_convergence']>0,
                **x))
    write_csv('optimization_by_blade_count.csv',results)
    write_csv('multistart_audit.csv',optim_runs)
    if not results: raise RuntimeError('No feasible numerical model solutions')
    winner=min(results,key=lambda r:r['energy_wh'])
    main3=next(r for r in results if r['blades']==3)
    assert all(r['max_abs_trim_error']<.001 and r['min_margin']>=-1e-3 for r in results)

    # Figure 1: Optimized energy by discrete blade count, phase split.
    plt.figure()
    labels=[f"{z['blades']} blades" for z in results]
    x=np.arange(len(results))
    ho=np.array([z['hover_energy_wh'] for z in results]);cr=np.array([z['cruise_energy_wh'] for z in results])
    plt.bar(x,ho,label='Hover (2 min)',color='#367CC7')
    plt.bar(x,cr,bottom=ho,label='Cruise (30 min)',color='#D38D45')
    for i,z in enumerate(results):plt.text(i,z['energy_wh']+2,f"{z['energy_wh']:.1f}",ha='center',fontsize=10)
    plt.xticks(x,labels);plt.ylabel('Electrical mission energy (Wh)')
    plt.ylim(0,max(z['energy_wh'] for z in results)*1.18)
    plt.title('Constrained, thrust-trimmed rotor solutions (surrogate polar)')
    plt.legend(frameon=False);savefig('01_optimized_energy_by_blades.png')

    # Figure 2: operating speed and energy, each option.
    plt.figure()
    for z in results:
        plt.scatter(z['omega_hover_rad_s'],z['hover_power_total_w'],s=95,marker='o',label=f"{z['blades']} blades")
        plt.scatter(z['omega_cruise_rad_s'],z['cruise_power_total_w'],s=95,marker='x',color='black')
    plt.ylabel('Total electrical input (W)');plt.xlabel('Rotational speed (rad/s)')
    plt.title('Trimmed hover (circle) and cruise (cross) operating points')
    plt.legend(frameon=False);savefig('02_trim_operating_points.png')

    # Figure 3: 3-bladed fixed chord radius-pitch energy map. Re-trim both RPMs.
    rr=np.linspace(.18,.30,12);tt=np.linspace(12,33,12)
    E=np.full((len(tt),len(rr)),np.nan);V=np.zeros_like(E,dtype=bool)
    maps=[]
    for i,theta in enumerate(tt):
        for j,R in enumerate(rr):
            z=trim_geometry(float(R),.025,float(theta),3,m)
            if z:
                margins=z['margins']
                feasible=(min(margins.values())>=-1e-3 and
                          max(abs(z['trim_residual']))<1e-3)
                E[i,j]=z['energy_wh'];V[i,j]=feasible
                maps.append(dict(radius_m=R,pitch75_deg=theta,chord_m=.025,
                                 energy_wh=z['energy_wh'],feasible=int(feasible),
                                 omega_hover=z['geometry']['omega_hover_rad_s'],
                                 omega_cruise=z['geometry']['omega_cruise_rad_s'],
                                 hover_fm=z['hover']['figure_of_merit']))
    write_csv('radius_pitch_trim_sweep.csv',maps)
    plt.figure(figsize=(8,5))
    masked=np.ma.masked_where(~V,E)
    if np.any(V):
        ctr=plt.contourf(rr,tt,masked,levels=17,cmap='viridis')
        plt.colorbar(ctr,label='Mission electrical energy (Wh)')
    plt.contourf(rr,tt,~V,levels=[.5,1.5],colors=['#DDDDDD'],alpha=.65)
    plt.xlabel('Radius (m)');plt.ylabel(r'Pitch at $0.75R$ (deg)')
    plt.title('3-blade radius/pitch sweep, chord fixed at 0.025 m\nGrey = fails bounds, trim, or physical constraints')
    savefig('03_feasible_radius_pitch_map.png')

    # Figure 4: thrust against RPM, unlike earlier plots demonstrates feasibility.
    plt.figure()
    for axial,label,omega,p_req in [(0.,'Hover',main3['omega_hover_rad_s'],m.hover_thrust_n),
                                    (m.cruise_speed_mps,'Cruise',main3['omega_cruise_rad_s'],m.cruise_thrust_n)]:
        speeds=np.linspace(160,650,45)
        thrust=[predict(main3['radius_m'],main3['chord_m'],o,
                        main3['theta75_deg'],3,axial,m)['thrust_n'] for o in speeds]
        plt.plot(speeds,thrust,label=f'{label} predicted thrust')
        plt.scatter([omega],[p_req],s=75,zorder=5)
        plt.axhline(p_req,linestyle='--',alpha=.35)
    plt.axhline(0,color='black',linewidth=.5)
    plt.xlabel('Rotor angular speed (rad/s)');plt.ylabel('Thrust per rotor (N)')
    plt.title('Thrust capability and exact mission trim — 3-blade option')
    plt.legend(frameon=False);savefig('04_thrust_trim_curves.png')

    # Figure 5: important *input* sensitivity, hold optimized geometry fixed,
    # retrim the RPM for each varying mass and wing cruise drag input.
    masses=[4.5,5.,5.5,6.];drags=[5.,7.,9.,11.]
    interface=[];sens=np.full((len(masses),len(drags)),np.nan)
    for i,mass in enumerate(masses):
        for j,drag in enumerate(drags):
            scenario=replace(m,mass_kg=mass,cruise_drag_n=drag)
            z=trim_geometry(main3['radius_m'],main3['chord_m'],main3['theta75_deg'],3,scenario)
            valid=bool(z and min(z['margins'].values())>=-1e-3)
            if z:
                interface.append(dict(mass_kg=mass,wing_drag_n=drag,
                    energy_wh=z['energy_wh'],feasible=int(valid),
                    hover_power_w=scenario.n_rotors*z['hover']['electrical_power_w'],
                    cruise_power_w=scenario.n_rotors*z['cruise']['electrical_power_w'],
                    motor_torque_hover_nm=z['hover']['torque_nm'],
                    motor_torque_cruise_nm=z['cruise']['torque_nm']))
            if valid:sens[i,j]=z['energy_wh']
    write_csv('wing_mass_drag_interface.csv',interface)
    plt.figure(figsize=(7.3,4.6))
    image=plt.imshow(sens,origin='lower',cmap='YlGnBu',aspect='auto',
                     interpolation='nearest')
    plt.colorbar(image,label='Energy (Wh), fixed rotor geometry / retrimmed speeds')
    plt.xticks(range(len(drags)),drags);plt.yticks(range(len(masses)),masses)
    plt.xlabel('Wing-provided cruise drag (N)')
    plt.ylabel('System-provided takeoff mass (kg)')
    for i in range(len(masses)):
        for j in range(len(drags)):
            plt.text(j,i,f'{sens[i,j]:.0f}' if np.isfinite(sens[i,j]) else 'X',
                     ha='center',va='center',fontsize=9,color='black')
    plt.title('Sensitivity to teammates’ shared inputs\nX = no feasible trim/constraint compliance')
    savefig('05_mass_drag_coupling.png')

    # Figure 6: BEMT discretization convergence with geometry/RPM held fixed.
    mesh=[]
    for elements in (16,24,32,48,64,96,128,192):
        mm=replace(m,n_elements=elements)
        h=predict(main3['radius_m'],main3['chord_m'],
                  main3['omega_hover_rad_s'],main3['theta75_deg'],3,0.,mm)
        c=predict(main3['radius_m'],main3['chord_m'],
                  main3['omega_cruise_rad_s'],main3['theta75_deg'],3,mm.cruise_speed_mps,mm)
        mesh.append(dict(elements=elements,hover_thrust_n=h['thrust_n'],
                         hover_power_w=h['electrical_power_w'],
                         cruise_thrust_n=c['thrust_n'],
                         cruise_power_w=c['electrical_power_w'],
                         converged=int(h['converged'] and c['converged'])))
    write_csv('mesh_convergence.csv',mesh)
    plt.figure()
    plt.plot([z['elements'] for z in mesh],[z['hover_thrust_n'] for z in mesh],
             marker='o',label='Hover thrust (N)')
    plt.plot([z['elements'] for z in mesh],[z['cruise_thrust_n'] for z in mesh],
             marker='s',label='Cruise thrust (N)')
    plt.xlabel('Blade elements');plt.ylabel('Predicted thrust per rotor (N)')
    plt.title('BEMT radial discretization convergence')
    plt.legend(frameon=False);savefig('06_bemt_mesh_convergence.png')

    winners={str(r['blades']):r for r in results}
    with (TABLE/'study_manifest.json').open('w',encoding='utf-8') as file:
        json.dump({'model_status':'illustrative surrogate polar, no experimental validation',
                   'mission':asdict(m),'by_blades':winners,
                   'analytical_screening':analytical,'optimization_multistart':optim_runs,
                   'default_solution':'3-blade case used for figures 2-6; energy optimum across blade counts listed separately'},
                   file,indent=2)
    print('Analytical hover cross-check:',analytical)
    print('Best by blade count:')
    for r in results:print(r)
    print('Feasible sweep points:',int(V.sum()),'/',V.size)
    print('Figures:',sorted(f.name for f in FIG.glob('*.png')))

if __name__=='__main__':main()
