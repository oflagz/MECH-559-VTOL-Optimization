"""Sprint 3: real-data BEMT comparison and optimization verification.

Run from project root: python -m scripts.sprint3_study
Does NOT calibrate the illustrative section polar or assert hardware validity.
"""
from pathlib import Path
import csv,json
from dataclasses import replace,asdict
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from rotoropt.experimental import predict_reference,DATA
from rotoropt.optimize import solve,analyze
from rotoropt.diagnostics import kkt_audit
from rotoropt.model import Mission

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results'/'sprint3';OUT.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.size':10,'figure.figsize':(7.4,4.7),'savefig.dpi':180,
                     'axes.spines.top':False,'axes.spines.right':False})

def savefig(name):
    plt.grid(alpha=.22)
    plt.tight_layout()
    plt.savefig(OUT/name,bbox_inches='tight');plt.close()

def read_csv(path):
    with open(path,newline='') as f:return list(csv.DictReader(f))

def write_csv(path,rows):
    if rows:
        with open(OUT/path,'w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)

def mean_abs_pct(rows,col):
    x=[100*abs(r['pred_'+col]-r['meas_'+col])/abs(r['meas_'+col])
       for r in rows if r['meas_'+col]!=0]
    return float(np.mean(x)) if x else None

def measured_benchmark():
    static=[]
    for r in read_csv(DATA/'apcsf_10x7_static.csv'):
        rpm=float(r['rpm']);a=predict_reference(rpm)
        static.append({'rpm':rpm,'meas_CT':float(r['CT']),'pred_CT':a['CT'],
                       'meas_CP':float(r['CP']),'pred_CP':a['CP'],
                       'converged':a['converged']})
    dynamic=[]
    for r in read_csv(DATA/'apcsf_10x7_3008rpm.csv'):
        J=float(r['J']);a=predict_reference(3008,J)
        dynamic.append({'J':J,'meas_CT':float(r['CT']),'pred_CT':a['CT'],
                        'meas_CP':float(r['CP']),'pred_CP':a['CP'],
                        'meas_eta':float(r['eta']), 'pred_eta':a['eta'],
                        'converged':a['converged'],
                        'propulsive_valid':bool(a['CT']>0 and a['CP']>0)})
    write_csv('uiuc_static_comparison.csv',static)
    write_csv('uiuc_forward_comparison.csv',dynamic)
    x=[r['rpm'] for r in static]
    for coef,label in [('CT',r'Thrust coefficient $C_T$'),('CP',r'Power coefficient $C_P$')]:
        plt.figure()
        plt.plot(x,[r['meas_'+coef] for r in static],marker='o',label='UIUC measured')
        plt.plot(x,[r['pred_'+coef] for r in static],marker='s',ls='--',label='BEMT surrogate')
        plt.xlabel('Propeller RPM');plt.ylabel(label)
        plt.title('APC Slow Flyer 10×7: static performance')
        plt.legend();savefig('01_uiuc_static_'+coef.lower()+'.png')
    for coef,label in [('CT',r'Thrust coefficient $C_T$'),('CP',r'Power coefficient $C_P$')]:
        plt.figure()
        plt.plot([r['J'] for r in dynamic],[r['meas_'+coef] for r in dynamic],marker='o',label='UIUC measured')
        plt.plot([r['J'] for r in dynamic],[r['pred_'+coef] for r in dynamic],marker='s',ls='--',label='BEMT surrogate')
        plt.axhline(0,color='grey',lw=.8)
        plt.xlabel(r'Advance ratio $J=V/(nD)$');plt.ylabel(label)
        plt.title('APC Slow Flyer 10×7: 3008 RPM, forward flight')
        plt.legend();savefig('02_uiuc_forward_'+coef.lower()+'.png')
    # MAPE only where both observed and predicted are operating as positive-thrust/power propeller;
    # excludes windmilling / negative-power regime unsupported by the model.
    positive=[r for r in dynamic if r['propulsive_valid'] and r['meas_CT']>0 and r['meas_CP']>0]
    summary={
      'data_source':'UIUC Propeller Database Volume 1 / Brandt and Selig; APC SF 10x7',
      'geometry_status':'published radial chord/pitch distributions',
      'polar_status':'illustrative synthetic polar; NOT validated / NOT calibrated',
      'static_count':len(static),'dynamic_count':len(dynamic),
      'static_CT_mape_pct':mean_abs_pct(static,'CT'),
      'static_CP_mape_pct':mean_abs_pct(static,'CP'),
      'forward_positive_count':len(positive),
      'forward_CT_mape_pct':mean_abs_pct(positive,'CT'),
      'forward_CP_mape_pct':mean_abs_pct(positive,'CP'),
      'note':'Large prediction errors and wrong sign near windmilling show model is not experimentally predictive. No parameter calibration was performed.'
    }
    return summary

def optimization_checks():
    baseline=Mission()
    base,runs,starts=solve(baseline,n_blades=3,max_starts=3)
    assert base is not None
    a=base['analysis']; audit=kkt_audit(base['x'],3,baseline)
    rows=[]
    for mass in (5.,5.5,6.):
        for drag in (5.,7.,9.):
            m=replace(baseline,mass_kg=mass,cruise_drag_n=drag)
            best,runs2,_=solve(m,3,max_starts=2)
            row=dict(takeoff_mass_kg=mass,cruise_drag_n=drag,feasible=best is not None,
                     nominal_case=int(mass==5. and drag==7.))
            if best:
                z=best['analysis'];row.update(energy_wh=z['energy_wh'],
                     radius_m=best['x'][0],chord_m=best['x'][1],
                     omega_hover=best['x'][2],omega_cruise=best['x'][3],
                     pitch75_deg=best['x'][4],
                     max_abs_trim_residual=float(max(abs(z['trim_residual']))),
                     min_margin=float(min(z['margins'].values())))
            else:
                row.update(energy_wh=np.nan,radius_m=np.nan,chord_m=np.nan,
                     omega_hover=np.nan,omega_cruise=np.nan,pitch75_deg=np.nan,
                     max_abs_trim_residual=np.nan,min_margin=np.nan)
            rows.append(row)
    write_csv('reoptimized_mass_drag_sweep.csv',rows)
    arr=np.array([r['energy_wh'] for r in rows]).reshape(3,3)
    plt.figure(figsize=(7.4,5))
    masked=np.ma.masked_invalid(arr)
    plt.imshow(masked,origin='lower',aspect='auto',cmap='viridis')
    plt.colorbar(label='Reoptimized rotor mission energy (Wh)')
    plt.xticks(range(3),[5,7,9]);plt.yticks(range(3),[5,5.5,6])
    plt.xlabel('Cruise drag from wing, N');plt.ylabel('Gross takeoff mass, kg')
    plt.title('SLSQP reoptimization after changing coupled inputs')
    for i in range(3):
        for j in range(3):
            plt.text(j,i,f'{arr[i,j]:.1f}' if np.isfinite(arr[i,j]) else 'infeasible',
                     ha='center',va='center',color='white' if np.isfinite(arr[i,j]) else 'black')
    savefig('03_reoptimized_mass_drag.png')
    # Mesh resolution: optimize geometry again, not just re-evaluate a fixed solution.
    high=replace(baseline,n_elements=192)
    finer, runs3, _=solve(high,3,max_starts=3)
    mesh_comp={
      'n96':{'energy_wh':a['energy_wh'],'x':base['x'],
             'min_margin':float(min(a['margins'].values()))},
      'n192':None if finer is None else {
             'energy_wh':finer['analysis']['energy_wh'],'x':finer['x'],
             'min_margin':float(min(finer['analysis']['margins'].values()))},
      'n96_solution_evaluated_at_n192': None,
    }
    z=analyze(base['x'],3,high)
    mesh_comp['n96_solution_evaluated_at_n192']={
        'energy_wh':z['energy_wh'],
        'hover_trim_error_rel':float(z['trim_residual'][0]),
        'cruise_trim_error_rel':float(z['trim_residual'][1]),
        'min_margin':float(min(z['margins'].values()))}
    return audit,rows,mesh_comp

def main():
    bench=measured_benchmark()
    audit,scenarios,mesh=optimization_checks()
    result={'benchmark':bench,'kkt_3blade':audit,'resolution_reoptimization':mesh,
            'screening_only':True,'team_mass_drag':'scenario parameters, not teammate-approved values',
            'course_topics':['L8 KKT and regularity','L9 SQP and active set','sensitivity and local-vs-global optimality']}
    (OUT/'sprint3_manifest.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))
    print('SCENARIOS')
    for r in scenarios:print(r['takeoff_mass_kg'],r['cruise_drag_n'],r['feasible'],r['energy_wh'])

if __name__=='__main__':main()
