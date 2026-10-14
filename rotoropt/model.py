"""Consistent thrust, torque and shaft-power prediction for axial propeller flow.

NOT VALIDATED: section polars are smooth *illustrative* surrogates, not
experimental data. Use results to test optimization/code logic, not to
select hardware. All quantities SI except blade pitch (degrees).
"""
from dataclasses import dataclass, replace
import numpy as np

@dataclass(frozen=True)
class Mission:
    # The 5 kg is provisional GROSS TAKEOFF mass (not the proposal's 5 kg payload).
    mass_kg: float = 5.0
    rho: float = 1.225
    gravity: float = 9.81
    n_rotors: int = 4
    download_fraction: float = .05
    thrust_margin: float = .15
    cruise_drag_n: float = 7.0
    cruise_speed_mps: float = 20.0
    hover_seconds: float = 120.0
    cruise_distance_m: float = 36000.0
    eta_motor: float = .85
    speed_of_sound: float = 340.0
    max_tip_mach: float = .65
    min_hover_fm: float = .58  # screening assumption, not manufacturer constraint
    total_max_electrical_w: float = 1000.
    motor_max_torque_nm: float = 1.5
    mount_spacing_m: float = .65
    clearance_m: float = .05
    hub_fraction: float = .20
    n_elements: int = 96
    max_iter: int = 90
    relax: float = .32
    tolerance: float = 5e-5
    # Intentionally *not* a calibrated airfoil or rotor mass model.

    @property
    def hover_thrust_n(self):
        return self.mass_kg*self.gravity*(1+self.download_fraction)/self.n_rotors

    @property
    def cruise_thrust_n(self):
        return self.cruise_drag_n/self.n_rotors

    @property
    def cruise_seconds(self):
        return self.cruise_distance_m/self.cruise_speed_mps


def polar(alpha_radians, reynolds):
    """Illustrative smooth section polar; deliberately NOT Reynolds-calibrated."""
    cl = 1.18 * np.tanh(2*np.pi*alpha_radians/1.18)
    cd = .015 + .027*cl**2 + .12*np.maximum(abs(alpha_radians)-.23,0)**2
    return cl, cd


def predict(radius_m, chord_m, omega_rad_s, pitch75_deg, n_blades,
            axial_speed_mps, mission=Mission(), loss=True):
    """Annular BEMT, propeller sign convention; Prandtl tip/hub factors.

    Returns diagnostics. Assumes constant chord and constant pitch; theta75
    is an optimized pitch variable, twist is fixed to zero at this stage.
    Does not model windmilling, reverse flow, compressibility or elastic blades.
    """
    if min(radius_m, chord_m, omega_rad_s) <= 0 or axial_speed_mps < 0:
        raise ValueError('positive geometry, speed and axial flow required')
    m=mission
    edges=np.linspace(m.hub_fraction*radius_m,radius_m,m.n_elements+1)
    r=(edges[:-1]+edges[1:])/2
    dr=np.diff(edges)
    vi=np.full_like(r,max(.05,np.sqrt(m.hover_thrust_n/(2*m.rho*np.pi*radius_m**2))))
    swirl=np.zeros_like(r)
    pitch=np.deg2rad(pitch75_deg)
    converged=False
    def forces(vi,swirl):
        va=axial_speed_mps+vi
        vt=np.maximum(omega_rad_s*r-swirl,1.e-8)
        phi=np.arctan2(va,vt)
        vrel=np.hypot(va,vt)
        reynolds=m.rho*vrel*chord_m/(1.81e-5)
        cl,cd=polar(pitch-phi,reynolds)
        lift=.5*m.rho*vrel**2*chord_m*n_blades*cl*dr
        drag=.5*m.rho*vrel**2*chord_m*n_blades*cd*dr
        thrust=lift*np.cos(phi)-drag*np.sin(phi)
        torque=r*(lift*np.sin(phi)+drag*np.cos(phi))
        return thrust,torque,phi,reynolds
    for iteration in range(1,m.max_iter+1):
        dt,dq,phi,reynolds=forces(vi,swirl)
        if loss:
            sinphi=np.maximum(np.abs(np.sin(phi)),.02)
            ftip=n_blades/2*(radius_m-r)/(r*sinphi)
            fhub=n_blades/2*(r-m.hub_fraction*radius_m)/(m.hub_fraction*radius_m*sinphi)
            F=np.maximum(.05,(2/np.pi*np.arccos(np.exp(-ftip))) *
                           (2/np.pi*np.arccos(np.exp(-fhub))))
        else:
            F=np.ones_like(r)
        # Momentum annulus: dT/dr = 4*pi*r*rho*F*(V+vi)*vi.
        load=np.maximum(dt/dr,0)
        target_vi=(np.sqrt(axial_speed_mps**2+load/(np.pi*m.rho*F*r))-
                   axial_speed_mps)/2
        # Angular momentum annulus, using tangential induced speed at disk.
        target_swirl=np.maximum(dq/dr,0)/(4*np.pi*m.rho*F*r*r*
                                              np.maximum(axial_speed_mps+target_vi,0.05))
        target_swirl=np.minimum(target_swirl,.6*omega_rad_s*r)
        vi_new=(1-m.relax)*vi+m.relax*target_vi
        swirl_new=(1-m.relax)*swirl+m.relax*target_swirl
        difference=max(np.max(abs(vi_new-vi)),np.max(abs(swirl_new-swirl)))
        vi,swirl=vi_new,swirl_new
        if difference<m.tolerance:
            converged=True
            break
    dt,dq,phi,reynolds=forces(vi,swirl)
    thrust=float(sum(dt))
    torque=float(sum(dq))
    shaft_power=torque*omega_rad_s
    ideal_hover_power=(max(thrust,0)**1.5 /
                       np.sqrt(2*m.rho*np.pi*radius_m**2))
    fm=ideal_hover_power/shaft_power if axial_speed_mps==0 and shaft_power>0 else np.nan
    return {'thrust_n':thrust,'torque_nm':torque,'shaft_power_w':shaft_power,
            'electrical_power_w':shaft_power/m.eta_motor,
            'figure_of_merit':fm,'converged':converged,'iterations':iteration,
            'reynolds_min':float(min(reynolds)),'reynolds_max':float(max(reynolds)),
            'tip_mach':np.hypot(omega_rad_s*radius_m,axial_speed_mps)/m.speed_of_sound,
            'prop_efficiency':thrust*axial_speed_mps/shaft_power if shaft_power>0 else np.nan}


def hover_momentum_check(radius_m,chord_m,omega,n_blades,mission=Mission()):
    """Independent analytical check at REQUIRED thrust, not BEMT capacity."""
    t=mission.hover_thrust_n
    area=np.pi*radius_m**2
    vi=np.sqrt(t/(2*mission.rho*area))
    solidity=n_blades*chord_m/(np.pi*radius_m)
    pi=t*vi
    p0=(solidity*.012/8)*mission.rho*area*(omega*radius_m)**3
    shaft=1.15*pi+p0
    return {'thrust_n':t,'ideal_induced_w':pi,'profile_w':p0,
            'shaft_w':shaft,'figure_of_merit':pi/shaft,
            'electrical_total_w':mission.n_rotors*shaft/mission.eta_motor}
