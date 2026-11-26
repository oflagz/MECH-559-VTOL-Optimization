"""Independent verification using MECH 559 Lectures 10, 12, and 14."""
from dataclasses import replace
from itertools import product
from pathlib import Path
import csv
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import minimize

from rotor_design.bemt import SectionPolar
from rotor_design.mission import Mission
from rotor_design.solve import BOUNDS, optimal_trim

ROOT=Path(__file__).resolve().parents[1]
TAB=ROOT/'results/tables'; FIG=ROOT/'results/figures'; REP=ROOT/'report'
LOW=np.array([bounds[0] for bounds in BOUNDS]);WIDTH=np.array([b[1]-b[0] for b in BOUNDS]);
PITCH=20.0


def writecsv(filename,rows):
    if not rows: raise ValueError('No rows to save: '+filename)
    keys=list(dict.fromkeys(k for row in rows for k in row))
    with (TAB/filename).open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)


def energies(q,m):
    factor=m.rotors/(3600*m.eta_motor_esc)
    return factor*q['hover'].shaft_power*m.hover_seconds, factor*q['cruise'].shaft_power*m.cruise_seconds


def domain_g(q,m):
    """All physical constraints, plus proposal-specific polar/J-domain constraints."""
    hover_angle=max(abs(q['hover'].min_alpha_deg),abs(q['hover'].max_alpha_deg))
    cruise_angle=max(abs(q['cruise'].min_alpha_deg),abs(q['cruise'].max_alpha_deg))
    J=m.cruise_mps/(q['design'][3]/(2*np.pi)*2*q['design'][0])
    return np.r_[q['g'],hover_angle/25.-1,cruise_angle/25.-1,J/.96-1]


def qualified(q,m,tol=1.e-5):
    return q is not None and q['valid'] and np.max(np.abs(q['eq']))<2e-4 and np.max(domain_g(q,m))<=tol


def dominance_mask(values,tol=1e-9):
    """Strict Pareto dominance for minimization; returns mask over INPUT rows.

    Works with ties and is not predicated on convexity or smoothness.
    """
    v=np.asarray(values,dtype=float)
    if v.ndim!=2 or v.shape[1]!=2: raise ValueError('need (n,2) objectives')
    if len(v)==0:return np.array([],dtype=bool)
    if not np.isfinite(v).all():raise ValueError('nonfinite objectives')
                                                                               
    return np.array([not np.any(np.all(v<=p+tol,axis=1)&np.any(v<p-tol,axis=1)) for p in v])


class Evaluator:
    def __init__(self,blades,polar,mission=Mission(),n=64):
        self.blades=blades;self.polar=polar;self.m=mission;self.n=n
        self.cache={};self.calls=0

    def get(self,z):
        zz=np.asarray(z,dtype=float)
        if np.any(zz<-.015) or np.any(zz>1.015):return None
        key=tuple(np.round(zz,9))
        if key not in self.cache:
            self.calls+=1
            self.cache[key]=optimal_trim(list(LOW+zz*WIDTH),self.blades,PITCH,self.polar,self.m,n=self.n)
        return self.cache[key]

    def constraints(self,z):
        q=self.get(z)
        return -domain_g(q,self.m) if q is not None else np.full(16,-10.)

    def objective(self,z,part='total'):
        q=self.get(z)
        if q is None:return 50.                                                       
        eh,ec=energies(q,self.m)
        return {'total':eh+ec,'hover':eh,'cruise':ec}[part]/400.


def dfo_comparison(summary,p,m):
    """Derivative-free inequality-constrained solution, compared at three blade counts."""
    results=[]
    for item in summary['blade_optima']:
        if 'energy_wh' not in item:continue
        B=int(item['blades'])
        g0=np.array([item['radius_m'],item['mean_chord_m'],item['twist_deg']])
        z0=(g0-LOW)/WIDTH
        ev=Evaluator(B,p,m,n=64)
                                                                                 
        for label,start in [('perturbed',np.clip(z0+np.array([.028,.016,-.018]),.005,.995)),('slsqp_seed',z0)]:
            o=minimize(lambda z:ev.objective(z),start,method='COBYLA',
                   bounds=[(0,1)]*3,
                   constraints=[{'type':'ineq','fun':ev.constraints}],
                   options={'rhobeg':.055,'tol':3e-6,'catol':3.e-6,'maxiter':800})
            q=ev.get(o.x)
            ok=qualified(q,m,tol=1.e-5)
            results.append({'blades':B,'start':label,'solver':'COBYLA','success':bool(o.success),
                       'feasible':bool(ok),'energy_wh':float(q['energy_wh']) if q else None,
                       'slsqp_energy_wh':float(item['energy_wh']),
                       'difference_wh':float(q['energy_wh']-item['energy_wh']) if q else None,
                       'max_ineq_residual':float(np.max(domain_g(q,m))) if q else None,
                       'max_trim_residual':float(np.max(abs(q['eq']))) if q else None,
                       'radius_m':float(q['design'][0]) if q else None,
                       'mean_chord_m':float(q['design'][1]) if q else None,
                       'twist_deg':float(q['design'][4]) if q else None,
                       'nfev':int(o.nfev),'message':str(o.message),'model_calls_total':ev.calls})
    writecsv('derivative_free_comparison.csv',results)
    return results


def pareto_epsilon_study(summary,p,m):
    """Biobjective epsilon-constraint optimizations; filter by exact Pareto dominance.

    All potential points are re-evaluated with the same physical constraints.
    No global completeness guarantee is made for this nonlinear constrained problem.
    """
    best=summary['optimum'];B=int(best['blade_count'])
    sample=[]
    with (TAB/'pareto_sample.csv').open(newline='') as f:
        for v in csv.DictReader(f):sample.append(v)
    candidates=[]
    def add(q,source,epsilon=None):
        if not qualified(q,m):return
        eh,ec=energies(q,m)
        candidates.append({'hover_wh':eh,'cruise_wh':ec,'total_wh':eh+ec,
             'radius_m':q['design'][0],'mean_chord_m':q['design'][1],
             'twist_deg':q['design'][4], 'source':source,'epsilon_wh':epsilon})
    ev=Evaluator(B,p,m,n=64)
    zopt=(np.array([best['design'][0],best['design'][1],best['design'][4]])-LOW)/WIDTH
    add(ev.get(zopt),'SLSQP minimum total')
    for v in sample:
        z=(np.array([float(v['radius_m']),float(v['chord_m']),float(v['twist_deg'])])-LOW)/WIDTH
        add(ev.get(z),'feasible sampling')
    if not candidates:raise RuntimeError('No feasible candidate points for Pareto analysis')
    hover_start=min(row['hover_wh'] for row in candidates)
    hover_end=max(row['hover_wh'] for row in candidates)
    epsilon_targets=np.linspace(hover_start+0.10,hover_end,11)
    rows=[]
    for eps in epsilon_targets:
        seeds=[r for r in candidates if r['hover_wh']<=eps+1e-5]
        if not seeds:continue
                                                                                   
        seeds.sort(key=lambda r:r['cruise_wh'])
        seed=seeds[0]
        z0=(np.array([seed['radius_m'],seed['mean_chord_m'],seed['twist_deg']])-LOW)/WIDTH
        def epsilon_constraint(z):
            q=ev.get(z)
            if q is None:return -10.
            return (eps-energies(q,m)[0])/400.
        opt=minimize(lambda z:ev.objective(z,part='cruise'),z0,method='SLSQP',
                     bounds=[(0,1)]*3,
                     constraints=[{'type':'ineq','fun':ev.constraints},
                                  {'type':'ineq','fun':epsilon_constraint}],
                     options={'ftol':1e-9,'maxiter':160})
        q=ev.get(opt.x)
        passed=bool(qualified(q,m) and q is not None and energies(q,m)[0]<=eps+0.025)
        rows.append({'epsilon_wh':float(eps),'success':bool(opt.success),
              'feasible':passed,'nfev':int(opt.nfev),'nit':int(opt.nit),
              'hover_wh':float(energies(q,m)[0]) if q else None,
              'cruise_wh':float(energies(q,m)[1]) if q else None,
              'max_g':float(np.max(domain_g(q,m))) if q else None})
        if passed:add(q,'epsilon-constrained SLSQP',float(eps))
    writecsv('epsilon_constraint_solves.csv',rows)
                                                                                              
    uni={}
    for v in candidates:
        key=tuple(round(v[k],5) for k in ('hover_wh','cruise_wh'))
        if key not in uni or v['source']=='epsilon-constrained SLSQP':uni[key]=v
    candidates=list(uni.values())
    mask=dominance_mask([[v['hover_wh'],v['cruise_wh']] for v in candidates])
    for row,nd in zip(candidates,mask):row['nondominated']=bool(nd)
    candidates.sort(key=lambda v:v['hover_wh'])
    writecsv('pareto_qualified_archive.csv',candidates)
    nd=sorted([v for v in candidates if v['nondominated']],key=lambda v:v['hover_wh'])
    fig,ax=plt.subplots(figsize=(7.7,4.9))
    other=[r for r in candidates if not r['nondominated']]
    ax.scatter([v['hover_wh'] for v in other],[v['cruise_wh'] for v in other],
               s=19,alpha=.28,label='dominated, qualified designs')
    ax.plot([v['hover_wh'] for v in nd],[v['cruise_wh'] for v in nd],'-o',
            markersize=4,lw=2,label='nondominated archive (approximation)')
    stars=[v for v in candidates if v['source']=='SLSQP minimum total']
    if stars:
        s=stars[0]
        ax.scatter([s['hover_wh']],[s['cruise_wh']],marker='*',s=170,c='black',zorder=8,label='min total energy')
        ax.annotate(f"{s['total_wh']:.2f} Wh total",(s['hover_wh'],s['cruise_wh']),xytext=(8,9),textcoords='offset points')
    ax.set(xlabel='Hover energy (Wh)',ylabel='Cruise energy (Wh)',
           title='Feasible two-objective archive / epsilon-constraint study')
    ax.grid(alpha=.2);ax.legend(fontsize=8);fig.tight_layout();fig.savefig(FIG/'07_hover_cruise_pareto.png',dpi=175);plt.close(fig)
    return rows,candidates,nd


def interval_sensitivity(summary,p,m):
    """Deterministic interval corners; ranges are scenarios, NOT estimated confidence bounds."""
    d=summary['optimum']['design']
    B=int(summary['optimum']['blade_count']);geometry=[d[0],d[1],d[4]]
                                                                                  
    intervals={
       'nonprop_airframe_kg':(0.40,0.70),                           
       'wing_parasite_drag_n':(m.wing_drag_parasite_n*.90,m.wing_drag_parasite_n*1.10),
       'download_fraction':(0.022,0.042),
       'drive_efficiency':(0.80,0.90),
    }
    param_keys=list(intervals)
    scenarios=[('nominal',{},'nominal')]
                                                                                          
    for bits in product((0,1),repeat=len(param_keys)):
        vals={name:intervals[name][bit] for name,bit in zip(param_keys,bits)}
        scenarios.append(('corner-'+''.join(map(str,bits)),vals,''.join(map(str,bits))))
    rows=[]
    for name,values,tag in scenarios:
        mi=replace(m,
            fuselage_tail_avionics_kg=values.get('nonprop_airframe_kg',m.fuselage_tail_avionics_kg),
            wing_drag_parasite_n=values.get('wing_parasite_drag_n',m.wing_drag_parasite_n),
            download_fraction=values.get('download_fraction',m.download_fraction),
            eta_motor_esc=values.get('drive_efficiency',m.eta_motor_esc))
        q=optimal_trim(geometry,B,PITCH,p,mi,n=64)
        g=domain_g(q,mi) if q else None
        gkeys=('Hover thrust reserve','Cruise thrust reserve','Hover tip Mach','Cruise tip Mach',
               'Hover electric power','Cruise electric power','Figure of merit',
               'Hover torque','Cruise torque','Rotor clearance','Reserve torque',
               'Hover reserve power','Battery energy','Hover incidence','Cruise incidence','Advance ratio')
        residuals={f'g_{i:02d}':float(v) if g is not None else None for i,v in enumerate(g if g is not None else [None]*16)}
        violated=[label for label,v in zip(gkeys,g) if v>1e-5] if g is not None else []
        rows.append({'scenario':name,'corner':tag,
           **{k:values.get(k,{'nonprop_airframe_kg':m.fuselage_tail_avionics_kg,
                             'wing_parasite_drag_n':m.wing_drag_parasite_n,
                             'download_fraction':m.download_fraction,
                             'drive_efficiency':m.eta_motor_esc}[k]) for k in param_keys},
           'trim_found':bool(q is not None),'all_constraints_satisfied':qualified(q,mi),
           'energy_wh':float(q['energy_wh']) if q else None,
           'aircraft_mass_kg':float(q['mass_kg']) if q else None,
           'drag_n':float(q['drag_n']) if q else None,
           'max_constraint_residual':float(np.max(g)) if g is not None else None,
           'battery_margin_wh':float(mi.usable_battery_wh-q['energy_wh']) if q else None,
           'battery_violation': bool(g[12]>1e-5) if g is not None else None,
           'violated_constraints':'; '.join(violated), **residuals})
    writecsv('interval_uncertainty_corners.csv',rows)
    corners=[r for r in rows if r['corner']!='nominal']
    finite=[r for r in rows if r['energy_wh'] is not None]
    nominal=next(r for r in rows if r['scenario']=='nominal')
    cats=np.arange(len(corners))
    fig,ax=plt.subplots(figsize=(11.3,5.2))
    for title,fil,style in (
        ('All constraints satisfied',lambda r:r['all_constraints_satisfied'],dict(color='#2563eb',marker='o')),
        ('Battery energy exceeded',lambda r:not r['all_constraints_satisfied'] and r['battery_violation'],dict(color='#dc2626',marker='x')),
        ('Other constraints exceeded',lambda r:not r['all_constraints_satisfied'] and not r['battery_violation'],dict(color='#d97706',marker='x'))):
        ids=[i for i,r in enumerate(corners) if fil(r)]
        ax.scatter(ids,[corners[i]['energy_wh'] for i in ids],label=title,s=72,zorder=5,**style)
    ax.axhline(m.usable_battery_wh,color='gray',ls='--',lw=1.3,label='382 Wh provisional battery limit')
    ax.scatter([],[],marker='*',c='black',s=140,label=f"Nominal: {nominal['energy_wh']:.2f} Wh")
    ax.set_xticks(cats,[r['corner'] for r in corners],rotation=45,ha='right')
    ax.set(xlabel='Interval endpoint code (mass, parasite drag, download, drive efficiency)',
           ylabel='Electrical rotor energy (Wh)',title='Sixteen re-trimmed uncertainty scenarios')
    ax.grid(axis='y',alpha=.2);ax.legend(ncol=2,fontsize=8)
    fig.tight_layout();fig.savefig(FIG/'10_interval_corner_energy.png',dpi=175);plt.close(fig)
    labels=['Figure of merit','Advance ratio','Hover reserve power','Battery energy']
    idx=[6,15,11,12]
    values=np.array([[max(0,float(r[f'g_{j:02d}']))*100 for j in idx] for r in corners])
    fig,ax=plt.subplots(figsize=(9.7,7.0))
    import matplotlib.colors as mcolors
    masked=np.ma.masked_less_equal(values,1e-3)
    hm=ax.imshow(masked,aspect='auto',cmap='YlOrRd',norm=mcolors.LogNorm(vmin=1e-3,vmax=max(1,values.max())))
    ax.set_facecolor('#f4f4f4')
    ax.set_xticks(range(4),labels)
    ax.set_yticks(range(16),[r['corner'] for r in corners])
    ax.set(xlabel='Constraint',ylabel='Interval endpoint code',title='Constraint exceedance (% of normalized limit)')
    for i in range(16):
        for j in range(4):
            v=values[i,j]
            if v>1e-3:
                ax.text(j,i,f'{v:.3g}%',ha='center',va='center',fontsize=7,color='black' if v<3 else 'white')
            else:
                ax.text(j,i,'–',ha='center',va='center',fontsize=7,color='gray')
    fig.colorbar(hm,ax=ax,label='Positive normalized exceedance (%)',fraction=.045,pad=.025)
    fig.tight_layout();fig.savefig(FIG/'11_interval_constraint_heatmap.png',dpi=175);plt.close(fig)
    return intervals,rows,{'corner_count':len(rows)-1,
      'trimmed_count':sum(x['trim_found'] for x in rows if x['scenario']!='nominal'),
      'qualified_count':sum(x['all_constraints_satisfied'] for x in rows if x['scenario']!='nominal'),
      'energy_min_trimmed_wh':float(min(r['energy_wh'] for r in finite)),
      'energy_max_trimmed_wh':float(max(r['energy_wh'] for r in finite)),
      'nominal_energy_wh':nominal['energy_wh'],
      'nominal_feasible':nominal['all_constraints_satisfied']}


def main():
    summary=json.loads((TAB/'study_summary.json').read_text())
    p=SectionPolar(**summary['calibrated_polar']);m=Mission()
    dfo=dfo_comparison(summary,p,m)
    eps,archive,nd=pareto_epsilon_study(summary,p,m)
    intervals,corners,unc=interval_sensitivity(summary,p,m)
    accepted=[r for r in dfo if r['feasible']]
    comparable=[r for r in accepted if int(r['blades'])==int(summary['optimum']['blade_count'])]
    outcome={'methods':{'primary':'SLSQP with exact thrust trim','derivative_free':'COBYLA inequality-constrained, no derivatives',
        'multiobjective':'epsilon constraint on hover energy; minimize cruise energy',
        'uncertainty':'nonprobabilistic interval-corner scenario screening'},
        'dfo_all':dfo,'dfo_qualified':len(accepted),'dfo_best_for_selected_blade_wh':min((r['energy_wh'] for r in comparable),default=None),
        'pareto_archive_count':len(archive),'pareto_nondominated_count':len(nd),
        'epsilon_solve_count':len(eps),'epsilon_qualified_count':sum(r['feasible'] for r in eps),
        'intervals':intervals,'uncertainty':unc,
        'interpretation':{'derivative_free':'Independent search is not a global-optimality proof; report all failures and objective differences.',
        'pareto':'Nondominated only within the feasible computed archive; not a proven continuous Pareto front.',
        'uncertainty':'Illustrative bounded scenarios, NOT confidence intervals or real-world distributions. No robust feasibility is established unless every scenario satisfies all constraints.',
        'hardware':'UIUC measurement transfer from 2-bladed reference to optimized 3-bladed rotor is not experimentally validated.',
        'coupling':'No closed wing/rotor/motor/battery MDA; mission voltage remains inconsistent with Damian motor.'}}
                                                                                  
    main_refs={int(v['blades']):v for v in summary['blade_optima'] if 'energy_wh' in v}
    report_rows=[]
    for B in (2,3,4):
        if B not in main_refs:continue
        ref_energy=float(main_refs[B]['energy_wh'])
        passed=[r for r in dfo if r['blades']==B and r['feasible']]
        if passed:
            benchmark=min(r['energy_wh'] for r in passed)
            report_rows.append(f'{B} & {ref_energy:.3f} & {benchmark:.3f} & Agreement within solver tolerance '+r'\\')
        else:
            report_rows.append(f'{B} & {ref_energy:.3f} & --- & No qualified COBYLA point '+r'\\')
    (REP/'dfo_table_rows.tex').write_text('\n'.join(report_rows)+'\n\\bottomrule\n')
    (TAB/'advanced_validation_summary.json').write_text(json.dumps(outcome,indent=2,default=lambda v: v.item() if hasattr(v,'item') else str(v)))
    macro=lambda name,value:'\\newcommand{\\'+name+'}{'+str(value)+'}\n'
    macros=macro('DFOCount',len(dfo))+macro('DFOQualified',len(accepted))
    macros+=macro('DFOBestSelected',f"{outcome['dfo_best_for_selected_blade_wh']:.2f}" if comparable else 'N/A')
    macros+=macro('ParetoNondominated',len(nd))+macro('ParetoArchive',len(archive))
    macros+=macro('ParetoEpsilonQualified',outcome['epsilon_qualified_count'])+macro('ParetoEpsilonCount',len(eps))
    macros+=macro('IntervalQualified',unc['qualified_count'])+macro('IntervalCornerCount',unc['corner_count'])
    macros+=macro('IntervalEnergyMin',f"{unc['energy_min_trimmed_wh']:.1f}")+macro('IntervalEnergyMax',f"{unc['energy_max_trimmed_wh']:.1f}")
    macros+=macro('IntervalFMViolations',sum(r['g_06']>1e-5 for r in corners if r['corner']!='nominal'))
    macros+=macro('IntervalJViolations',sum(r['g_15']>1e-5 for r in corners if r['corner']!='nominal'))
    macros+=macro('IntervalReservePViolations',sum(r['g_11']>1e-5 for r in corners if r['corner']!='nominal'))
    macros+=macro('IntervalBatteryViolations',sum(r['g_12']>1e-5 for r in corners if r['corner']!='nominal'))
    (REP/'advanced_results.tex').write_text(macros)
    print(json.dumps({'selected_energy_wh':summary['optimum']['energy_wh'],
       'qualified_dfo':len(accepted),'selected_dfo':outcome['dfo_best_for_selected_blade_wh'],
       'epsilon_qualified':outcome['epsilon_qualified_count'],'pareto_nondominated':len(nd),
       'uncertainty':unc},indent=2,default=lambda v: v.item() if hasattr(v,'item') else str(v)))

if __name__=='__main__':main()
