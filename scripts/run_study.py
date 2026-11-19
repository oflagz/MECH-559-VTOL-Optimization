"""Run full MECH559 rotor analysis; single reproducible entry point.

python -m scripts.run_study
Outputs report tables, plots and the JSON summary. No extrapolation.
"""
from pathlib import Path
from dataclasses import asdict, replace
import csv,json
import numpy as np
from scipy.optimize import minimize_scalar
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from rotor.catalog import CATALOG,MeasuredPropeller
from rotor.optimization import Mission,evaluate_catalog,solve_catalog
from rotor.verification import holdout_quality
from rotor.physics import bemt_diagnostic

ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/'results';(OUT/'figures').mkdir(parents=True,exist_ok=True)
(OUT/'tables').mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.size':10,'figure.figsize':(7.5,4.8),'savefig.dpi':170,'axes.spines.top':False,'axes.spines.right':False})

def csvout(name,rows):
    if not rows:return
    keys=list(dict.fromkeys(k for row in rows for k in row))
    with (OUT/'tables'/name).open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=keys);writer.writeheader();writer.writerows(rows)

def save(name):
    plt.grid(alpha=.15)
    plt.tight_layout()
    plt.savefig(OUT/'figures'/name)
    plt.close()

def main():
    # Scenarios are explicitly separated: official speed 40km/1200s vs 40km at 20m/s.
    slower=Mission(); official=replace(slower,cruise_mps=40000/1200.)
    catalog_results=evaluate_catalog(slower);official_results=evaluate_catalog(official)
    assert all(r['slsqp_success'] for r in catalog_results if r['feasible'])
    csvout('slower_cruise_catalog.csv',catalog_results)
    csvout('official_mission_catalog.csv',official_results)
    good=sorted([r for r in catalog_results if r['feasible']],key=lambda r:r['energy_wh'])
    best=good[0] if good else None
    # Purely numerical interpolation checks on withheld UIUC rows.
    holdout=holdout_quality();csvout('heldout_interpolation.csv',holdout)
    # Static thrust measured test curves & requirement.
    plt.figure()
    for name in CATALOG:
        p=MeasuredPropeller(name)
        n=np.linspace(p.static_rpm[0],p.static_rpm[-1],140)
        t=[p.static(r)['thrust_n'] for r in n]
        plt.plot(n,t,label=f'APC {name} (2 blades)')
        plt.scatter(p.static_rpm,[p.static(r)['thrust_n'] for r in p.static_rpm],s=12)
    plt.axhline(best['required_hover_n'] if best else 25.,ls='--',color='black',label='Selected-case hover load')
    plt.xlabel('RPM');plt.ylabel('Measured static thrust, N per rotor')
    plt.title('Static hover capability: actual UIUC wind-tunnel data')
    plt.legend();save('01_measured_hover_thrust.png')
    # Static hover power vs diameter at different mass choices.
    rows=[]
    for base in (8.5,9.5,10.5):
        for name in CATALOG:
            x=solve_catalog(name,replace(slower,base_aircraft_mass_kg=base))
            rows.append({'base_aircraft_mass_kg':base,'model':name,'diameter_m':x['diameter_m'],
                         'hover_rpm':x.get('hover_rpm'), 'feasible_mission':x['feasible'],
                         'reason':x['reason'],'hover_electrical_w':x.get('hover_electrical_total_w')})
    csvout('mass_case_candidates.csv',rows)
    # 18x12 measured curves and trim point.
    if best:
        p=MeasuredPropeller(best['propeller']);plt.figure()
        for rpm,(J,ct,cp) in p.forward.items():
            plt.plot(J,ct(J),label=f'{rpm:.0f} RPM measured')
        j=slower.cruise_mps/(best['cruise_rpm']/60*best['diameter_m'])
        ct=p.axial(best['cruise_rpm'],slower.cruise_mps)['ct']
        plt.scatter([j],[ct],color='black',s=80,zorder=5,label='Selected cruise trim')
        plt.xlabel('Advance ratio $J=V/(nD)$');plt.ylabel('Thrust coefficient $C_T$')
        plt.title('Measured APC 18x12E forward-flight curves and interpolated trim')
        plt.legend();save('02_measured_cruise_ct.png')
        # Forward power with measured CP.
        plt.figure()
        for rpm,(J,ct,cp) in p.forward.items():
            plt.plot(J,cp(J),label=f'{rpm:.0f} RPM measured')
        cp=p.axial(best['cruise_rpm'],slower.cruise_mps)['cp']
        plt.scatter([j],[cp],color='black',s=80,zorder=5,label='Selected cruise trim')
        plt.xlabel('Advance ratio $J=V/(nD)$');plt.ylabel('Power coefficient $C_P$')
        plt.title('Measured APC 18x12E shaft-power coefficient')
        plt.legend();save('03_measured_cruise_cp.png')

    # Vary only supplied parameters and re-trim; infeasible bins are blank, not fake zero.
    mgrid=(8.5,9.0,9.5,10.0,10.5)
    draggrid=(5,6,7,8,9)
    scenario=[];E=np.full((len(mgrid),len(draggrid)),np.nan)
    for i,base in enumerate(mgrid):
        for k,drag in enumerate(draggrid):
            m=replace(slower,base_aircraft_mass_kg=base,cruise_drag_n=drag)
            res=evaluate_catalog(m)
            feasible=[r for r in res if r['feasible']]
            win=min(feasible,key=lambda r:r['energy_wh']) if feasible else None
            if win:E[i,k]=win['energy_wh']
            scenario.append({'aircraft_mass_excluding_props_kg':base,'cruise_drag_n':drag,
                             'has_data_qualified_solution':bool(win),
                             'propeller':win['propeller'] if win else '',
                             'energy_wh':win['energy_wh'] if win else None})
    csvout('mass_drag_sensitivity.csv',scenario)
    plt.figure(figsize=(7.3,4.7));masked=np.ma.masked_invalid(E)
    img=plt.imshow(masked,origin='lower',aspect='auto',interpolation='nearest',cmap='cividis')
    plt.xticks(range(len(draggrid)),draggrid);plt.yticks(range(len(mgrid)),mgrid)
    plt.xlabel('Cruise drag (N)');plt.ylabel('Base aircraft mass excluding propellers (kg)')
    plt.title('Re-trimmed rotor energy: 40 km at 20 m/s, 5-min hover')
    plt.colorbar(img,label='Electrical rotor energy (Wh)')
    for i in range(len(mgrid)):
        for k in range(len(draggrid)):
            txt=f'{E[i,k]:.0f}' if np.isfinite(E[i,k]) else 'N/A'
            plt.text(k,i,txt,ha='center',va='center',fontsize=8,
                     color='white' if np.isfinite(E[i,k]) and E[i,k] < np.nanmedian(E) else 'black')
    save('04_mass_drag_sensitivity.png')

    # Cruise speed sweep (includes strict measured envelope; speed vs flown time tradeoff).
    vv=np.arange(15.,33.334,.25)
    speeds=[]
    for v in vv:
        m=replace(slower,cruise_mps=float(v))
        r=evaluate_catalog(m)
        viable=[x for x in r if x['feasible']]
        winner=min(viable,key=lambda x:x['energy_wh']) if viable else None
        speeds.append({'cruise_mps':v,'duration_minutes':m.distance_m/m.cruise_mps/60,
                       'feasible_measured_domain':bool(winner),
                       'energy_wh':winner['energy_wh'] if winner else None,
                       'propeller':winner['propeller'] if winner else ''})
    csvout('cruise_speed_envelope.csv',speeds)
    plt.figure()
    valid=[r for r in speeds if r['feasible_measured_domain']]
    plt.plot([r['cruise_mps'] for r in valid],[r['energy_wh'] for r in valid],marker='o',ms=3,
             label='Feasible measured-domain candidates')
    plt.axvline(20,color='gray',ls='--',label='Prior screening cruise speed')
    plt.axvline(40000/1200,color='black',ls=':',label='Official 20-min cruise requirement')
    plt.xlabel('Cruise speed (m/s)');plt.ylabel('Rotor electrical mission energy (Wh)')
    plt.title('Available measured cruise-trim envelope (not full aircraft feasibility)')
    plt.legend(fontsize=8);save('05_speed_energy_envelope.png')

    # Evidence of model audit: old 'optimum' fails non-clipping momentum regime.
    def surrogate(a,re):
        cl=1.18*np.tanh(2*np.pi*a/1.18)
        cd=.015+.027*cl*cl+.12*np.maximum(abs(a)-.23,0)**2
        return cl,cd
    from rotor.physics import bemt_diagnostic
    oldh=bemt_diagnostic(.254102,.018,254.6827,33,3,0.,surrogate)
    oldc=bemt_diagnostic(.254102,.018,197.8851,33,3,20.,surrogate)
    # Additional engineering sensitivities: uncertainty in actual diameter affects both CT and CP.
    diams=[]
    if best:
        for ratio in (.99,1.,1.01):
            original=CATALOG[best['propeller']]
            from rotor import catalog
            catalog.CATALOG[best['propeller']]=replace(original,diameter_m=original.diameter_m*ratio)
            x=solve_catalog(best['propeller'],slower)
            catalog.CATALOG[best['propeller']]=original
            diams.append(dict(diameter_ratio=ratio,feasible=x['feasible'],
                              energy_wh=x.get('energy_wh'),hover_rpm=x.get('hover_rpm'),cruise_rpm=x.get('cruise_rpm')))
    csvout('nominal_diameter_sensitivity.csv',diams)
    out={'model_type':'strict experimental-interpolation catalogue; original continuous BEMT not used for design selection',
         'official_mission':asdict(official),'slower_cruise_alternative':asdict(slower),
         'catalogue_results_alternative':catalog_results,'catalogue_results_official':official_results,
         'best_alternative':best,'interpolation_holdout':holdout,
         'BEMT_old_operating_point':{'hover':asdict(oldh),'cruise':asdict(oldc)},
         'speed_qualified_range_mps':[min(x['cruise_mps'] for x in valid),max(x['cruise_mps'] for x in valid)] if valid else None,
         'mass_drag_energy_range_wh':[float(np.nanmin(E)),float(np.nanmax(E))],
         'nominal_diameter_sensitivity':diams,
         'warning':'No catalogued test-supported rotor can be declared feasible for official 40km/20min mission. Nominal diameter, mass and motor parameter assumptions remain.'}
    (OUT/'tables'/'study_summary.json').write_text(json.dumps(out,indent=2,default=float))
    print('BEST',best)
    print('OFFICIAL',[(q['propeller'],q['reason']) for q in official_results])
    print('SPEED_RANGE',out['speed_qualified_range_mps'])
    print('DIAM',diams)
    print('OLD_BEMT',oldh,oldc)
    print('HOLDOUT_MEDIAN',float(np.median([q['normalized_mae_pct'] for q in holdout if q['condition']=='forward'])))

if __name__=='__main__':main()
