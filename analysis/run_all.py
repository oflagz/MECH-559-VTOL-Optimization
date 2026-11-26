"""Single reproducible numerical study: official mission, no data from former rotor repo."""
from dataclasses import asdict
from pathlib import Path
import json,csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import nnls
from rotor_design.bemt import RotorGeometry
from rotor_design.fast_bemt import fast_performance
from rotor_design.calibrate import training_data,testing_data
from rotor_design.calibration_check import calibrate_fast,quality
from rotor_design.solve import optimal_trim,optimize_one,BOUNDS
from rotor_design.mission import Mission,evaluate_design,prop_mass,takeoff_mass,cruise_drag

ROOT=Path(__file__).resolve().parents[1]
FIG=ROOT/'results/figures';TAB=ROOT/'results/tables';REP=ROOT/'report'
for folder in (FIG,TAB,REP):folder.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.size':10,'figure.figsize':(7.3,4.5),'savefig.dpi':180,
                     'axes.spines.top':False,'axes.spines.right':False,'figure.autolayout':False})
COLORS=['#245273','#ca763e','#4b8361']

def save(n):
    plt.grid(alpha=.19)
    plt.tight_layout(pad=1.1)
    plt.savefig(FIG/n)
    plt.close()

def csvwrite(path,rows):
    if not rows:return
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with (TAB/path).open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=keys);writer.writeheader();writer.writerows(rows)

def coeff_validation(p):
    tests=[]
    for cohort,rows in [('fit',training_data()),('heldout',testing_data())]:
        for typ,rpm,j,ct,cp in rows:
            g=RotorGeometry(.127,.022,rpm*2*np.pi/60.,-10.,2,14.38)
            v=fast_performance(g,j*rpm/60*.254,p,n_sections=64,reference=True)
            tests.append({'cohort':cohort,'condition':typ,'rpm':rpm,'J':j,
                'CT_measured':ct,'CP_measured':cp,'CT_BEMT':v.ct if v.valid else None,
                'CP_BEMT':v.cp if v.valid else None,'admissible':v.valid,'reason':v.reason})
    csvwrite('benchmark_measured_vs_bemt.csv',tests)
    plt.figure(figsize=(8,4.7))
    for i,(key,y_lab) in enumerate([('CT','Thrust coefficient'),('CP','Power coefficient')]):
        plt.subplot(1,2,i+1)
        for cohort,marker,color in [('fit','o',COLORS[0]),('heldout','s',COLORS[1])]:
            subset=[r for r in tests if r['cohort']==cohort and r['admissible']]
            plt.scatter([r[key+'_measured'] for r in subset],[r[key+'_BEMT'] for r in subset],
                        label=cohort,alpha=.75,s=22,marker=marker,color=color)
        both=[float(r[key+'_measured']) for r in tests if r['admissible']]+[float(r[key+'_BEMT']) for r in tests if r['admissible']]
        lo,hi=min(both),max(both)
        plt.plot([lo,hi],[lo,hi],'k--',lw=1,label='perfect match')
        plt.xlabel('Measured '+key);plt.ylabel('BEMT '+key)
        plt.title(y_lab)
        if i==0:plt.legend(fontsize=8)
    save('01_uiuc_benchmark.png')
    return tests

def reliability(q,m):
    R,c,oh,oc,tw=q['design']
    j=m.cruise_mps/((oc/(2*np.pi))*2*R)
    ah=max(abs(q['hover'].min_alpha_deg),abs(q['hover'].max_alpha_deg))
    ac=max(abs(q['cruise'].min_alpha_deg),abs(q['cruise'].max_alpha_deg))
    return {'J_cruise':j,'max_hover_abs_alpha_deg':ah,'max_cruise_abs_alpha_deg':ac,
            'J_domain_margin':.96-j,'alpha_hover_margin_deg':25-ah,'alpha_cruise_margin_deg':25-ac,
            'battery_margin_wh':m.usable_battery_wh-q['energy_wh']-m.extra_energy_wh}

def main():
    m=Mission()
    fit,p=calibrate_fast()
    train=quality(p,training_data());held=quality(p,testing_data())
    source=coeff_validation(p)
                                                                                       
                                                                                                
    family_pitch=20.
    batch=[];records=[]
    for B in (2,3,4):
        q,a=optimize_one(B,family_pitch,p,m,starts=(7 if B!=4 else 3),n=64,enforce_battery=True,maxiter=(75 if B!=4 else 35))
        if q is None:
                                                                                                      
            q0,a0=optimize_one(B,family_pitch,p,m,starts=4,n=64,enforce_battery=False,maxiter=55)
            q=q0
            a.extend(a0)
        batch.append(q)
        records+=a
    csvwrite('multistart.csv',[{'blade_count':r['blade_count'],'start_id':r['start'],
        'slsqp_success':r['success'],'feasible':r['feasible'],'iterations':r['iterations'],
        'nfev':r['nfev'],'energy_wh':r['energy_wh'],'largest_normalized_violation':r['max_g']}
        for r in records])
    allc=[]
    for B,q in zip((2,3,4),batch):
        if q is None:allc.append({'blades':B,'status':'no converged trim'});continue
        rel=reliability(q,m)
        allc.append({'blades':B,'status':'battery feasible' if q['g'][-1]<=1e-5 else 'battery deficit',
            'energy_wh':q['energy_wh'],'radius_m':q['design'][0],'mean_chord_m':q['design'][1],
            'hover_rad_s':q['design'][2],'cruise_rad_s':q['design'][3],
            'twist_deg':q['design'][4],'mass_kg':q['mass_kg'],**rel})
    csvwrite('blade_count_optima.csv',allc)
    valid=[q for q in batch if q is not None and q['g'][-1]<=1e-5 and reliability(q,m)['J_domain_margin']>=-1e-4]
    best=min(valid,key=lambda x:x['energy_wh']) if valid else None
    if not best:
                                                                                          
        best=min((q for q in batch if q is not None),key=lambda x:x['energy_wh'],default=None)
    if best is None:raise RuntimeError('No valid BEMT trim found; cannot fabricate a design result')
    B=best['blades']; x=best['design']; R,c,oh,oc,tw=x
                                                                                                  
    meshes=[]
    for n in (24,32,48,64,96,128,192):
        q=evaluate_design(x,B,family_pitch,m,p,sections=n)
        meshes.append({'annuli':n,'energy_wh':q['energy_wh'],
            'hover_thrust_n':q['hover'].thrust,'cruise_thrust_n':q['cruise'].thrust,
            'hover_trim_residual':q['eq'][0], 'cruise_trim_residual':q['eq'][1],
            'all_induction_converged':q['valid']})
    csvwrite('mesh_convergence.csv',meshes)
    plt.figure()
    plt.plot([v['annuli'] for v in meshes],[v['energy_wh'] for v in meshes],'-o',color=COLORS[0])
    plt.axhline(meshes[-1]['energy_wh'],color='black',ls='--',lw=1,label='192-annulus reference')
    plt.xlabel('Radial blade elements');plt.ylabel('Computed energy (Wh)')
    plt.title('BEMT mesh-refinement audit at fixed design');plt.legend()
    save('09_mesh_convergence.png')

                                                                                    
    plt.figure()
    xs=[];yy=[]
    for row in allc:
        if row.get('energy_wh') is not None:xs.append(row['blades']);yy.append(row['energy_wh'])
    colors=[COLORS[1] if v!=B else COLORS[2] for v in xs]
    bars=plt.bar([str(b) for b in xs],yy,color=colors,edgecolor='white')
    for bar,v in zip(bars,yy):plt.text(bar.get_x()+bar.get_width()/2,v+2,f'{v:.1f} Wh',ha='center',fontsize=9)
    plt.axhline(m.usable_battery_wh,color='black',ls=':',label='Existing pack usable: 382 Wh')
    plt.ylim(0,max(yy)*1.16);plt.xlabel('Number of blades per rotor');plt.ylabel('Rotor electrical energy (Wh)')
    plt.title('Enumerated blade-count optima; official mission');plt.legend()
    save('08_blade_count_comparison.png')

                                                                                               
    sweeps=[]
    for idx,tag,scale,units in [(0,'radius',.18,'m'),(1,'mean_chord',.30,'m'),(4,'blade_twist',.0,'deg')]:
        low,high=({0:(R-.027,R+.027),1:(max(.009,c-.004),c+.004),4:(tw-12,tw+12)}[idx])
        xx=np.linspace(low,high,33)
        data=[]
        for v in xx:
            geo=[R,c,tw]
            geo[{0:0,1:1,4:2}[idx]]=float(v)
            q=optimal_trim(geo,B,family_pitch,p,m,n=64)
            if q:
                rel=reliability(q,m)
                cond=(np.max(q['g'])<=1e-3 and rel['max_hover_abs_alpha_deg']<=25.001 and
                      rel['max_cruise_abs_alpha_deg']<=25.001 and rel['J_cruise']<=.9601)
                data.append((float(v),q['energy_wh'],bool(cond)))
                sweeps.append({'variable':tag,'value':float(v),'energy_wh':q['energy_wh'],
                    'constraints_satisfied':bool(cond),'hover_speed':q['design'][2],
                    'cruise_speed':q['design'][3]})
            else:data.append((float(v),np.nan,False))
        plt.figure()
        ok=[a for a in data if a[2]];no=[a for a in data if not a[2] and np.isfinite(a[1])]
        if ok:plt.plot([z[0] for z in ok],[z[1] for z in ok],'-o',markersize=3,color=COLORS[0],label='fully feasible')
        if no:plt.scatter([z[0] for z in no],[z[1] for z in no],s=18,color=COLORS[1],alpha=.7,label='violates constraint')
        plt.scatter([x[idx]],[best['energy_wh']],s=105,marker='*',color='black',zorder=5,label='selected optimum')
        plt.annotate(f'{x[idx]:.4g} {units}',(x[idx],best['energy_wh']),xytext=(7,-19),textcoords='offset points',fontsize=9)
        plt.xlabel(tag.replace('_',' ').title()+f' ({units})');plt.ylabel('Mission energy (Wh)')
        plt.title(f'One-variable re-trimmed sensitivity: {tag.replace("_"," ")}')
        plt.legend(fontsize=8);save(f'{["02","03","04"][[0,1,4].index(idx)]}_{tag}_sweep.png')
    csvwrite('design_variable_sweeps.csv',sweeps)

                                                                                                    
    rpmrows=[]
    for V,idx,name,number in [(0.,2,'hover',5),(m.cruise_mps,3,'cruise',6)]:
        opt=x[idx];omegas=np.linspace(opt*.80,opt*1.20,75)
        thrust=[];power=[];status=[]
        for omega in omegas:
            geom=RotorGeometry(R,c,float(omega),tw,B,family_pitch)
            q=fast_performance(geom,V,p,n_sections=64)
            thrust.append(q.thrust if q.valid else np.nan)
            power.append(q.shaft_power if q.valid else np.nan)
            status.append(q.valid)
            rpmrows.append({'phase':name,'omega_rad_s':float(omega),'thrust_n':q.thrust if q.valid else None,
                            'shaft_power_w':q.shaft_power if q.valid else None,'valid':q.valid})
        fig,ax=plt.subplots(figsize=(7.5,4.6));ax.plot(omegas,thrust,color=COLORS[0],label='BEMT thrust')
        target=best['hover_required_n'] if name=='hover' else best['cruise_required_n']
        ax.axhline(target,color='black',ls='--',lw=1.0,label='required thrust')
        ax.scatter([opt],[target],marker='*',s=125,color='black',zorder=5)
        ax.annotate(f'Selected {opt:.1f} rad/s',xy=(opt,target),xytext=(9,10),textcoords='offset points',fontsize=9)
        ax.set_xlabel('Angular speed (rad/s)');ax.set_ylabel('Thrust per rotor (N)')
        ax2=ax.twinx();ax2.plot(omegas,power,color=COLORS[1],ls=':',lw=2,label='shaft power')
        ax2.set_ylabel('Shaft power per rotor (W)');ax.grid(alpha=.2)
        ax.set_title(f'{name.title()} speed optimization: thrust balance and power')
        lines=ax.get_lines()+ax2.get_lines()
        ax.legend(lines,[v.get_label() for v in lines],loc='upper left',fontsize=8)
        fig.tight_layout();fig.savefig(FIG/f'{number:02d}_{name}_rpm_trim.png');plt.close(fig)
    csvwrite('rpm_operating_sweeps.csv',rpmrows)

                                                                                                  
    samples=[];rng=np.random.default_rng(701)
    for k in range(700):
        geo=[np.clip(R+rng.normal(0,.009),BOUNDS[0][0],BOUNDS[0][1]),
             np.clip(c+rng.normal(0,.002),BOUNDS[1][0],BOUNDS[1][1]),
             np.clip(tw+rng.normal(0,5),BOUNDS[2][0],BOUNDS[2][1])]
        q=optimal_trim(geo,B,family_pitch,p,m,n=32)
        if not q:continue
        rel=reliability(q,m)
        if not (np.max(q['g'])<=.002 and rel['max_hover_abs_alpha_deg']<=25.01 and
                rel['max_cruise_abs_alpha_deg']<=25.01 and rel['J_cruise']<=.961):continue
        e_h=4*q['hover'].shaft_power*m.hover_seconds/(3600*m.eta_motor_esc)
        e_c=4*q['cruise'].shaft_power*m.cruise_seconds/(3600*m.eta_motor_esc)
        samples.append({'hover_wh':e_h,'cruise_wh':e_c,'total_wh':e_h+e_c,
                        'radius_m':q['design'][0],'chord_m':q['design'][1],
                        'twist_deg':q['design'][4]})
    csvwrite('pareto_sample.csv',samples)
    eh=4*best['hover'].shaft_power*m.hover_seconds/(3600*m.eta_motor_esc)
    ec=4*best['cruise'].shaft_power*m.cruise_seconds/(3600*m.eta_motor_esc)
    if samples:
        plt.figure()
        for r in samples:plt.scatter(r['hover_wh'],r['cruise_wh'],color=COLORS[0],alpha=.25,s=9)
        s=sorted(samples,key=lambda v:v['hover_wh'])
        frontier=[];mc=float('inf')
        for r in s:
            if r['cruise_wh']<mc:
                frontier.append(r);mc=r['cruise_wh']
        plt.plot([v['hover_wh'] for v in frontier],[v['cruise_wh'] for v in frontier],'-',
                 color=COLORS[1],linewidth=2,label='non-dominated sample')
        plt.scatter([eh],[ec],marker='*',s=170,color='black',zorder=6,label='minimum-energy design')
        plt.annotate(f'({eh:.1f}, {ec:.1f}) Wh',(eh,ec),xytext=(8,9),textcoords='offset points',fontsize=9)
        plt.xlabel('Hover electrical energy (Wh)');plt.ylabel('Cruise electrical energy (Wh)')
        plt.title('Feasible hover/cruise energy trade-off (sampling, not proof of Pareto completeness)')
        plt.legend();save('07_hover_cruise_pareto.png')
    else:
                                                                  
        plt.figure();plt.scatter([eh],[ec],s=150,color='black',marker='*',label='selected design')
        plt.title('No additional feasible Pareto samples in sampled neighborhood');plt.legend()
        plt.xlabel('Hover electrical energy (Wh)');plt.ylabel('Cruise electrical energy (Wh)')
        save('07_hover_cruise_pareto.png')

                                                                                                          
    low=np.array([b[0] for b in BOUNDS]);width=np.array([b[1]-b[0] for b in BOUNDS]);z=(np.array([R,c,tw])-low)/width
    def metrics(z):
        q=optimal_trim(list(low+np.asarray(z)*width),B,family_pitch,p,m,n=64)
        if q is None:return None
        rel=reliability(q,m)
        return q['energy_wh']/400.,np.r_[q['g'],
            rel['max_hover_abs_alpha_deg']/25-1,rel['max_cruise_abs_alpha_deg']/25-1,
            rel['J_cruise']/.96-1]
    center=metrics(z);jac=np.zeros((len(center[1]),3));grad=np.zeros(3)
    eps=2e-5
    for j in range(3):
        u=z.copy();v=z.copy();u[j]+=eps;v[j]-=eps
        aa=metrics(u);bb=metrics(v)
        if aa is not None and bb is not None:
            grad[j]=(aa[0]-bb[0])/(2*eps);jac[:,j]=(aa[1]-bb[1])/(2*eps)
        else:
            grad[j]=np.nan;jac[:,j]=np.nan
    constraints=center[1]
    names=['hover thrust reserve','cruise thrust reserve','hover tip Mach','cruise tip Mach',
           'hover electrical power','cruise electrical power','minimum figure of merit',
           'hover shaft torque','cruise shaft torque','mounting clearance','hover reserve torque',
           'hover reserve power','usable battery energy','hover polar angle','cruise polar angle','advance ratio domain']
    active=np.where(constraints>-5e-4)[0].tolist()
    bound_active=[]
    for k,(name,lo,hi) in enumerate(zip(['R','c','twist'],z,np.ones(3)-z)):
        if lo<5e-4:bound_active.append((name+' lower',k,-1.))
        if hi<5e-4:bound_active.append((name+' upper',k,1.))
    J=np.array([jac[k] for k in active]) if active else np.zeros((0,3))
    if bound_active:
        ext=np.zeros((len(bound_active),3))
        for i,(_,k,sgn) in enumerate(bound_active):ext[i,k]=sgn
        J=np.vstack([J,ext])
    if J.size and np.all(np.isfinite(J)) and np.all(np.isfinite(grad)):
        multipliers,norm=nnls(J.T,-grad)
        licq=int(np.linalg.matrix_rank(J))
        stationarity=float(np.linalg.norm(grad+J.T@multipliers))
    else:
        multipliers=[];licq=0;stationarity=float('nan')
    kkt={'active_constraints':[names[i] for i in active]+[q[0] for q in bound_active],
         'jacobian_rank':licq,'active_constraint_count':len(J),
         'gradient_norm':float(np.linalg.norm(grad)),'stationarity_residual':stationarity,
         'nonnegative_multipliers':[float(x) for x in multipliers],
         'note':'Reduced KKT diagnostic on three geometry variables; two thrust equalities eliminated by exact bracketed trim.'}

    out={'mission':asdict(m),'cruise_speed_mps':m.cruise_mps,'fixed_theta75_deg':family_pitch,
         'calibrated_polar':asdict(p),'fit_success':bool(fit.success),'fit_cost':float(fit.cost),
         'benchmark_training':{k:v for k,v in train.items() if k not in ('observations','invalid_points')},
         'benchmark_holdout':{k:v for k,v in held.items() if k not in ('observations','invalid_points')},
         'benchmark_invalid_holdout':held['invalid_points'],
         'blade_optima':allc,'optimum':{'blade_count':B,'design':x,'energy_wh':best['energy_wh'],
              'hover_shaft_power_w':best['hover'].shaft_power,'cruise_shaft_power_w':best['cruise'].shaft_power,
              'hover_torque_nm':best['hover'].torque,'cruise_torque_nm':best['cruise'].torque,
              'hover_figure_of_merit':best['hover'].fm,'cruise_propulsive_efficiency':best['cruise'].eta,
              'hover_thrust_n':best['hover'].thrust,'cruise_thrust_n':best['cruise'].thrust,
              'required_hover_n':best['hover_required_n'],'required_cruise_n':best['cruise_required_n'],
              'mass_kg':best['mass_kg'],'cruise_drag_n':best['drag_n'],
              'propeller_mass_each_kg':best['prop_mass_each_kg'],
              'hover_energy_wh':eh,'cruise_energy_wh':ec,
              'negative_cruise_loading_fraction':best['cruise'].negative_load_fraction,
              'max_constraint_violation':float(np.max(best['g'])),'trim_residuals':best['eq'].tolist(),
              **reliability(best,m)},'kkt':kkt,'pareto_sample_count':len(samples),
         'caveats':['Calibration is ONLY APC 10x7 Slow Flyer with two blades. Three- and four-blade geometry and optimized larger propellers have not been experimentally validated.',
                    'Some high-J measured points have no admissible low-induction momentum solution; coverage counts are reported.','No full mass-power MDO closure; fuselage and nacelle masses are explicit assumptions.',
                    'Usable pack energy excludes transition, avionics, weather, ageing and mandatory reserves; a small positive margin cannot imply full-system mission feasibility.',
                    'New high-pitch rotor family theta75=20 degrees is an assumed family parameter, not a measured manufacturer specification.']}
    (TAB/'study_summary.json').write_text(json.dumps(out,indent=2,default=float))
                                                                                                
    macros=r'''\newcommand{\OptBlade}{%d}
\newcommand{\OptEnergy}{%.2f}
\newcommand{\OptRadius}{%.4f}
\newcommand{\OptChord}{%.4f}
\newcommand{\OptTwist}{%.2f}
\newcommand{\OptHoverOmega}{%.2f}
\newcommand{\OptCruiseOmega}{%.2f}
\newcommand{\OptHoverRPM}{%.0f}
\newcommand{\OptCruiseRPM}{%.0f}
\newcommand{\OptMass}{%.3f}
\newcommand{\OptDrag}{%.2f}
\newcommand{\OptHoverThrust}{%.3f}
\newcommand{\OptCruiseThrust}{%.3f}
\newcommand{\OptHoverPower}{%.1f}
\newcommand{\OptCruisePower}{%.1f}
\newcommand{\OptHoverTorque}{%.3f}
\newcommand{\OptCruiseTorque}{%.3f}
\newcommand{\OptFM}{%.3f}
\newcommand{\OptEta}{%.3f}
\newcommand{\OptEH}{%.1f}
\newcommand{\OptEC}{%.1f}
\newcommand{\OptJ}{%.3f}
\newcommand{\OptPackMargin}{%.1f}
\newcommand{\BenchCT}{%.2f}
\newcommand{\BenchCP}{%.2f}
\newcommand{\BenchCount}{%d}
\newcommand{\BenchTotal}{%d}
\newcommand{\MeshDelta}{%.3f}
\newcommand{\KKTResidual}{%.2g}
\newcommand{\ParetoCount}{%d}
'''%(B,best['energy_wh'],R,c,tw,oh,oc,oh*60/(2*np.pi),oc*60/(2*np.pi),
           best['mass_kg'],best['drag_n'],best['hover'].thrust,best['cruise'].thrust,
           best['hover'].shaft_power,best['cruise'].shaft_power,best['hover'].torque,best['cruise'].torque,
           best['hover'].fm,best['cruise'].eta,eh,ec,reliability(best,m)['J_cruise'],
           reliability(best,m)['battery_margin_wh'],held['ct_median_relative_pct'],
           held['cp_median_relative_pct'],held['valid_n'],held['valid_n']+held['invalid_n'],
           abs(meshes[-1]['energy_wh']-meshes[3]['energy_wh'])/meshes[-1]['energy_wh']*100,
           stationarity,len(samples))
    (REP/'results.tex').write_text(macros)
    print(json.dumps({'best_blades':B,'E_wh':round(best['energy_wh'],3),'design':x,
      'battery_margin_wh':round(reliability(best,m)['battery_margin_wh'],3),
      'heldout_ct_pct':round(held['ct_median_relative_pct'],2),
      'heldout_cp_pct':round(held['cp_median_relative_pct'],2),
      'heldout_coverage':f"{held['valid_n']}/{held['valid_n']+held['invalid_n']}",
      'kkt':kkt,'pareto_samples':len(samples)},indent=2))

if __name__=='__main__':
    main()
                                                                              
                                                                        
    from .advanced_validation import main as advanced_validation_main
    advanced_validation_main()

