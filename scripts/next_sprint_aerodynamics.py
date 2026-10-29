"""Development study AFTER Sprint 3. Not part of the previous sprint report.

Run: python -m scripts.next_sprint_aerodynamics
Produces source-audited airfoil/geometry checks, model-vs-measurement
comparisons, loss/swirl sensitivity, and SLSQP reoptimization.
WARNING: Clark-Y measured polar is a PROXY for APC's E63/Clark-Y blend
and out-of-domain polars are endpoint CLAMPED WITH COVERAGE WARNING.
"""
from dataclasses import replace
from pathlib import Path
import csv,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from rotoropt.experimental import DATA as APC, load_geometry, predict_reference
from rotoropt.model import Mission
from rotoropt.optimize import solve
from rotoropt.polars import coefficients

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results'/'next_sprint'
OUT.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'figure.figsize':(7,4.6),'figure.dpi':110,
    'savefig.dpi':180,'axes.spines.top':False,'axes.spines.right':False})

def read_csv(path):
    with path.open(newline='') as f: return list(csv.DictReader(f))

def save_csv(name,rows):
    if not rows:return
    with (OUT/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]))
        w.writeheader();w.writerows(rows)

def fig(name):
    plt.grid(alpha=.20);plt.tight_layout()
    plt.savefig(OUT/name,bbox_inches='tight');plt.close()

def mape(rows, observed, forecast):
    return float(np.mean([100*abs(r[forecast]-r[observed])/abs(r[observed])
              for r in rows if r[observed]!=0]))

def evaluate():
    geometry=load_geometry()
    assert np.isclose(geometry['r_over_R'][0],.15)
    assert np.isclose(geometry['r_over_R'][-1],1.)
    assert np.isclose(geometry['pitch_deg'][np.argmin(abs(geometry['r_over_R']-.75))],14.38)
    assert np.all(geometry['c_over_R']>0)

    static=[];dynamic=[]
    static_measured=read_csv(APC/'apcsf_10x7_static.csv')
    dynamic_measured=read_csv(APC/'apcsf_10x7_3008rpm.csv')
    cases={'legacy':'surrogate','clarky_proxy':'clarky_uiuc'}
    for row in static_measured:
        rpm=float(row['rpm'])
        for label,mode in cases.items():
            p=predict_reference(rpm,polar_model=mode)
            static.append({'rpm':rpm,'model':label,'measured_CT':float(row['CT']),
             'predicted_CT':p['CT'],'measured_CP':float(row['CP']),
             'predicted_CP':p['CP'],'error_CT_pct':100*(p['CT']/float(row['CT'])-1),
             'error_CP_pct':100*(p['CP']/float(row['CP'])-1),
             're_min':p['re_min'],'re_max':p['re_max'],
             'alpha_min_deg':p['alpha_min_deg'],'alpha_max_deg':p['alpha_max_deg'],
             'measured_polar_coverage_weighted':p['measured_polar_coverage_weighted'],
             'converged':p['converged']})
    for row in dynamic_measured:
        j=float(row['J'])
        for label,mode in cases.items():
            p=predict_reference(3008,J=j,polar_model=mode)
            dynamic.append({'J':j,'model':label,'measured_CT':float(row['CT']),
             'predicted_CT':p['CT'],'measured_CP':float(row['CP']),
             'predicted_CP':p['CP'],
             'measured_positive_thrust':float(row['CT'])>0,
             'predicted_positive_thrust':p['CT']>0,
             'measured_polar_coverage_weighted':p['measured_polar_coverage_weighted'],
             'converged':p['converged']})
    save_csv('static_model_comparison.csv',static)
    save_csv('forward_model_comparison.csv',dynamic)
    values={}
    for label in cases:
        st=[r for r in static if r['model']==label]
        fw=[r for r in dynamic if r['model']==label]
        forward_positive_both=[r for r in fw if r['measured_CT']>0 and r['predicted_CT']>0 and r['measured_CP']>0 and r['predicted_CP']>0]
        values[label]={
         'static_CT_mape_pct':mape(st,'measured_CT','predicted_CT'),
         'static_CP_mape_pct':mape(st,'measured_CP','predicted_CP'),
         'forward_positive_both_CT_mape_pct':mape(forward_positive_both,'measured_CT','predicted_CT'),
         'forward_positive_both_CP_mape_pct':mape(forward_positive_both,'measured_CP','predicted_CP'),
         'forward_positive_both_count':len(forward_positive_both),
         'forward_measured_positive_count':sum(r['measured_CT']>0 for r in fw),
         'forward_wrong_sign_count':sum(r['measured_CT']>0 and r['predicted_CT']<=0 for r in fw),
         'static_coverage_min':min(r['measured_polar_coverage_weighted'] for r in st),
         'static_coverage_max':max(r['measured_polar_coverage_weighted'] for r in st),
         'all_bemt_converged':all(r['converged'] for r in st+fw),
        }

    for coeff in ['CT','CP']:
        plt.figure()
        x=[float(row['rpm']) for row in static_measured]
        plt.plot(x,[float(row[coeff]) for row in static_measured],marker='o',label='UIUC measured')
        for label in cases:
            v=[r for r in static if r['model']==label]
            plt.plot(x,[r['predicted_'+coeff] for r in v],marker='.',label=label)
        plt.xlabel('RPM');plt.ylabel(coeff)
        plt.title('APC SF 10x7 static — measured vs two polar models')
        plt.legend();fig('01_static_'+coeff.lower()+'.png')
    for coeff in ['CT','CP']:
        plt.figure()
        x=[float(row['J']) for row in dynamic_measured]
        plt.plot(x,[float(row[coeff]) for row in dynamic_measured],marker='o',label='UIUC measured')
        for label in cases:
            v=[r for r in dynamic if r['model']==label]
            plt.plot(x,[r['predicted_'+coeff] for r in v],marker='.',label=label)
        plt.axhline(0,color='grey',linewidth=.8)
        plt.xlabel('Advance ratio J');plt.ylabel(coeff)
        plt.title('APC SF 10x7 forward — measured vs models, 3008 RPM')
        plt.legend();fig('02_forward_'+coeff.lower()+'.png')
    plt.figure()
    x=[float(row['rpm']) for row in static_measured]
    v=[r for r in static if r['model']=='clarky_proxy']
    plt.plot(x,[100*r['measured_polar_coverage_weighted'] for r in v],marker='o')
    plt.axhline(90,color='grey',ls='--',label='90% reference coverage')
    plt.ylim(0,103);plt.xlabel('RPM')
    plt.ylabel('Weighted blade-section coverage (%)')
    plt.title('Where UIUC Clark-Y measured polar data apply')
    plt.legend();fig('03_reynolds_coverage.png')

    # Loss, swirl, and inboard airfoil uncertainty; no tuning to measurements.
    checks=[]
    for rpm in [3029.,5987.]:
        for polar_mode in ['surrogate','clarky_uiuc']:
            for loss in (True,False):
                for swirl in (True,False):
                    if not loss and not swirl:pass
                    p=predict_reference(rpm,polar_model=polar_mode,
                            use_tip_hub_loss=loss,include_swirl=swirl)
                    checks.append(dict(rpm=rpm,polar_model=polar_mode,
                        tip_hub_loss=loss,swirl=swirl,CT=p['CT'],CP=p['CP'],
                        coverage=p['measured_polar_coverage_weighted']))
    for rpm in [3029.,5987.]:
        p=predict_reference(rpm,polar_model='clarky_uiuc',
                      inner_section='surrogate_inboard')
        checks.append(dict(rpm=rpm,polar_model='clarky_with_surrogate_inboard',
                     tip_hub_loss=True,swirl=True,CT=p['CT'],CP=p['CP'],
                     coverage=p['measured_polar_coverage_weighted']))
    save_csv('induction_and_section_sensitivity.csv',checks)

    # Rerun SLSQP on the new polar, but formally reject a 'validated optimum'
    # when measured section data do not cover the operating range.
    optimal=[]
    for name in ['surrogate','clarky_uiuc']:
        scenario=replace(Mission(),polar_model=name)
        for blades in (2,3,4):
            best,runs,_=solve(scenario,blades,max_starts=3)
            if best is None:
                optimal.append(dict(model=name,blades=blades,solver_feasible=False,
                                    physically_supported=False));continue
            a=best['analysis']
            ch=a['hover']['measured_polar_coverage_weighted']
            cc=a['cruise']['measured_polar_coverage_weighted']
            physically_supported=(name=='clarky_uiuc' and ch>=.9 and cc>=.9)
            optimal.append(dict(model=name,blades=blades,solver_feasible=True,
              physically_supported=physically_supported,
              energy_wh=a['energy_wh'],radius_m=best['x'][0],
              chord_m=best['x'][1],omega_hover_rad_s=best['x'][2],
              omega_cruise_rad_s=best['x'][3],theta75_deg=best['x'][4],
              hover_coverage=ch,cruise_coverage=cc,
              all_converged=all(a[stage]['converged'] for stage in ['hover','cruise','reserve']),
              max_trim_residual=max(abs(a['trim_residual']))))
    save_csv('reoptimized_with_clarky_proxy.csv',optimal)
    plt.figure()
    for i,model in enumerate(['surrogate','clarky_uiuc']):
        rows=[r for r in optimal if r['model']==model and r['solver_feasible']]
        plt.bar(np.arange(len(rows)) + (i-.5)*.32,
                [r['energy_wh'] for r in rows],width=.32,
                label=model)
    plt.xticks(range(3),['2 blades','3 blades','4 blades'])
    plt.ylabel('Predicted energy, Wh')
    plt.title('SLSQP results: conditional screening ONLY\nMeasured-section coverage is insufficient for design approval')
    plt.legend();fig('04_slsqp_comparison.png')

    # Contrast the actual published 3D geometric section angle with the
    # 33 degree full-scale optimizer pitch assumption.
    plt.figure()
    plt.plot(geometry['r_over_R'],geometry['c_over_R'],marker='o',label='c / R')
    plt.xlabel('r / R');plt.ylabel('Chord / radius')
    plt.title('Measured APC SF 10x7 blade chord distribution')
    fig('05_geometry_chord.png')
    plt.figure()
    plt.plot(geometry['r_over_R'],geometry['pitch_deg'],marker='o')
    plt.axhline(14.38,ls='--',label='Measured pitch at 75% radius')
    plt.xlabel('r / R');plt.ylabel('Pitch angle, degrees')
    plt.title('Measured APC SF 10x7 blade pitch distribution')
    plt.legend();fig('06_geometry_pitch.png')

    summary={'source_blade':'UIUC APC SF 10x7, proprietary E63/Clark-Y blend',
      'tested_section':'UIUC measured Clark-Y(B) proxy, NOT exact APC section',
      'section_data_range':{'re_min':60000,'re_max':300000,
         'limitations':'drag alpha coverage varies by Reynolds; endpoint clipping flagged'},
      'benchmark':values,
      'optimization':optimal,
      'invalid_design_count':sum(not r.get('physically_supported',False) for r in optimal if r.get('model')=='clarky_uiuc'),
      'experimental_validation_success':False,
      'notes':['Do not count out-of-range Clark-Y endpoints as validated airfoil behavior.',
               'APC inboard section is an E63-like section, not measured Clark-Y.',
               'Static UIUC points do not directly provide motor/ESC electrical efficiency.',
               'No 5kg-payload system mass closure or actual actuator/motor data.']}
    (OUT/'diagnostic_summary.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2))
    return summary

if __name__=='__main__': evaluate()
