"""Mission and coupled input ledger transcribed independently from submitted proposal and teammates' PDFs."""
from dataclasses import dataclass,asdict
import numpy as np
from .bemt import RotorGeometry,SectionPolar,RHO,SPEED_SOUND
from .fast_bemt import fast_performance as performance

@dataclass(frozen=True)
class Mission:
                                                                              
    distance_m:float=40000.
    hover_seconds:float=300.
    cruise_seconds:float=1200.
    payload_kg:float=5.
    rotors:int=4
                                                                         
    wing_mass_kg:float=.314
    wing_drag_parasite_n:float=10.88
    wing_drag_induced_ref_n:float=1.29
    wing_reference_mass_kg:float=7.31
    download_fraction:float=.032
    installation_spacing_m:float=.65                                                     
    installation_clearance_m:float=.05
                                                                              
    propulsion_mass_kg:float=3.58
    usable_battery_wh:float=382.
                                                                              
    fuselage_tail_avionics_kg:float=.55
    nacelle_tilt_mass_each_kg:float=.09                                          
    eta_motor_esc:float=.85
    peak_shaft_torque_nm:float=1.5
    available_electrical_power_w:float=2400.
                                                                  
    max_tip_mach:float=.65
    min_fm:float=.58
    thrust_reserve:float=.15
                                                                                             
    reference_prop_mass_kg:float=.42*.028349523125
    reference_prop_radius_m:float=.127
    mass_radius_power:float=2.45
    reference_mean_chord_m:float=.023
                                                                              
    extra_energy_wh:float=0.
    @property
    def cruise_mps(self):return self.distance_m/self.cruise_seconds

def prop_mass(radius,chord,blades,m=Mission()):
    return m.reference_prop_mass_kg*(radius/m.reference_prop_radius_m)**m.mass_radius_power*(chord/m.reference_mean_chord_m)**.35*(blades/2)

def takeoff_mass(radius,chord,blades,m=Mission()):
    return (m.payload_kg+m.wing_mass_kg+m.propulsion_mass_kg+
            m.fuselage_tail_avionics_kg+m.rotors*(m.nacelle_tilt_mass_each_kg+prop_mass(radius,chord,blades,m)))

def cruise_drag(mass,m=Mission()):
                                                                          
                                                                               
    return m.wing_drag_parasite_n+m.wing_drag_induced_ref_n*(mass/m.wing_reference_mass_kg)**2

def evaluate_design(x,blades,theta75_deg,m=Mission(),polar=SectionPolar(),sections=48):
    R,c,omega_h,omega_c,twist=x
    gh=RotorGeometry(R,c,omega_h,twist,blades,theta75_deg)
    gc=RotorGeometry(R,c,omega_c,twist,blades,theta75_deg)
    mass=takeoff_mass(R,c,blades,m)
    W=mass*9.81;drag=cruise_drag(mass,m)
    required_hover=(1+m.download_fraction)*W/m.rotors
    required_cruise=drag/m.rotors
    h=performance(gh,0.,polar,n_sections=sections)
    cr=performance(gc,m.cruise_mps,polar,n_sections=sections)
                                                                                
    ghreserve=RotorGeometry(R,c,omega_h*1.10,twist,blades,theta75_deg)
    hr=performance(ghreserve,0.,polar,n_sections=sections)
    gcreserve=RotorGeometry(R,c,omega_c*1.10,twist,blades,theta75_deg)
    crreserve=performance(gcreserve,m.cruise_mps,polar,n_sections=sections)
    valid=h.valid and cr.valid and hr.valid and crreserve.valid and h.shaft_power>0 and cr.shaft_power>0
                                                                          
    g=np.array([
       ((1+m.thrust_reserve)*required_hover-hr.thrust)/required_hover,
       (required_cruise-crreserve.thrust)/required_cruise,
       omega_h*R/SPEED_SOUND-m.max_tip_mach,
       np.hypot(m.cruise_mps,omega_c*R)/SPEED_SOUND-m.max_tip_mach,
       m.rotors*h.shaft_power/(m.eta_motor_esc*m.available_electrical_power_w)-1,
       m.rotors*cr.shaft_power/(m.eta_motor_esc*m.available_electrical_power_w)-1,
       (m.min_fm-h.fm)/m.min_fm if np.isfinite(h.fm) else 9,
       (h.torque-m.peak_shaft_torque_nm)/m.peak_shaft_torque_nm,
       (cr.torque-m.peak_shaft_torque_nm)/m.peak_shaft_torque_nm,
       (2*R+m.installation_clearance_m-m.installation_spacing_m)/m.installation_spacing_m,
       (hr.torque-m.peak_shaft_torque_nm)/m.peak_shaft_torque_nm,
       m.rotors*hr.shaft_power/(m.eta_motor_esc*m.available_electrical_power_w)-1
    ])
    eq=np.array([h.thrust/required_hover-1,cr.thrust/required_cruise-1])
    e=(m.rotors/m.eta_motor_esc)*(h.shaft_power*m.hover_seconds+cr.shaft_power*m.cruise_seconds)/3600
    energy_residual=(e+m.extra_energy_wh-m.usable_battery_wh)/m.usable_battery_wh
    g=np.r_[g,energy_residual]
    if not valid or np.any(~np.isfinite(g)) or not np.isfinite(e):
        g=np.full(13,10.)
        eq=np.array([10.,10.]);e=1e6
    return {'valid':valid,'hover':h,'cruise':cr,'reserve':hr,'cruise_reserve':crreserve,'mass_kg':mass,'drag_n':drag,
            'hover_required_n':required_hover,'cruise_required_n':required_cruise,'energy_wh':e,
            'g':g,'eq':eq,'prop_mass_each_kg':prop_mass(R,c,blades,m),
            'design':list(map(float,x)),'blades':blades,'theta75_deg':theta75_deg}
