"""MECH 559 - Wing subsystem (Antoine). Design vector x = [b, AR, t/c, lam, A_cap].
Sweep removed (no effect on drag in this model). t/c enters CD0 (form factor) AND the
structure (spar depth), wing mass enters W and hover energy. Units SI unless noted.
PLACEHOLDER inputs are flagged: replace with common-mission / teammates' values."""
import numpy as np
from scipy.optimize import minimize
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------- parameters ----------------
P = dict(
    rho=1.225, g=9.81, mu=1.81e-5,
    Vc=40000/1200, tc=1200.0, th=300.0,      # common mission of the report (40 km, 20 min cruise, 5 min hover)
    m_other=7.0,                           # PLACEHOLDER kg: everything except wing (system model)
    CLmax=1.2, nult=4.5,                   # 3.0 g limit x 1.5 safety
    E_mod=60e9, sig_allow=300e6,           # spar-cap modulus / allowable (incl. knockdown) [assumed]
    rho_cap=1600.0, k_struct=2.5,          # CFRP density; web/ribs/joints factor on cap mass [assumed]
    rho_skin=0.45,                         # kg/m2 wetted area (skin+core+ribs) [assumed, calibrate: Kahvecioglu]
    h_frac=0.85,                           # spar depth / local thickness
    f_other=0.0134,                        # m2 drag area of fuselage, nacelles, tail (calibrated so D~7 N)
    eta_prop=0.60, eta_c=0.85, eta_h=0.85, FM=0.708, Nr=4, R=0.25,   # from proposal screening values
    CD_dl=0.4, ov=0.25, dl_max=0.05,       # download model (placeholder for Ashish's k_download), 5% of W
    b_max=3.0, AR_rng=(6, 14), tcr_rng=(0.10, 0.16), lam_rng=(0.4, 1.0), c_min=0.12,
    Acap_rng=(1e-6, 2e-4), defl_frac=0.15,
)
SC = np.array([1.0, 10.0, 0.1, 1.0, 1e-5])   # scaling: z = x / SC

def oswald(lam, AR, CD0):
    d = 0.0524*lam**4 - 0.15*lam**3 + 0.1659*lam**2 - 0.0706*lam + 0.0119  # taper term (Nita & Scholz)
    return 1.0/((1.0 + d) + 0.38*CD0*np.pi*AR)                              # Obert/Kroo form 1/(Q+P*pi*AR)

def beam(x, W, p=P, n=101):
    """Numerical Euler-Bernoulli half-wing, elliptical load, constant cap area, depth ~ local chord."""
    b, AR, tcr, lam, A = x
    S = b*b/AR; L = b/2; croot = 2*S/(b*(1+lam))
    y = np.linspace(0, L, n); s = y/L
    c = croot*(1-(1-lam)*s); h = p['h_frac']*tcr*c; I = A*h**2/2
    Ph = p['nult']*W/2
    w = 4*Ph/(np.pi*L)*np.sqrt(np.maximum(1-s**2, 0))
    def rcum(f):  # integral from y to L
        out = np.zeros_like(f); out[:-1] = np.cumsum(((f[:-1]+f[1:])/2*np.diff(y))[::-1])[::-1]; return out
    def cum(f):
        out = np.zeros_like(f); out[1:] = np.cumsum((f[:-1]+f[1:])/2*np.diff(y)); return out
    V = rcum(w); M = rcum(V)
    sig = M/(A*h)
    delta = cum(cum(M/(p['E_mod']*I)))
    return sig.max(), delta[-1], c[0]

def model(x, p=P):
    b, AR, tcr, lam, A = x
    S = b*b/AR; cbar = S/b
    Swet = 2*S*(1+0.25*tcr)
    m_cap = p['rho_cap']*A*2*b
    m_wing = p['rho_skin']*Swet + p['k_struct']*m_cap
    mtot = p['m_other'] + m_wing; W = mtot*p['g']
    q = 0.5*p['rho']*p['Vc']**2
    Re = p['rho']*p['Vc']*cbar/p['mu']
    Cf = 0.455/np.log10(Re)**2.58
    FF = 1 + 2*tcr + 100*tcr**4                      # Raymer, (x/c)m=0.3, unswept -> t/c affects drag
    CD0 = Cf*FF*Swet/S
    e = oswald(lam, AR, CD0); K = 1/(np.pi*e*AR)
    D = q*(S*CD0 + p['f_other']) + K*W**2/(q*S)
    Pcr = p['Vc']*D
    Vstall = np.sqrt(2*W/(p['rho']*S*p['CLmax']))
    # hover download (rotor wake on exposed wing) and hover energy with mass feedback
    T_rot = 1.05*W/p['Nr']; A_rot = np.pi*p['R']**2
    vi = np.sqrt(T_rot/(2*p['rho']*A_rot))
    Sexp = p['Nr']*(2*p['R'])*cbar*p['ov']
    Ddl = p['CD_dl']*0.5*p['rho']*(2*vi)**2*Sexp
    Ttot = W + Ddl
    Ph_el = Ttot**1.5/(np.sqrt(2*p['rho']*p['Nr']*A_rot)*p['FM']*p['eta_h'])
    E = (Pcr/(p['eta_prop']*p['eta_c'])*p['tc'] + Ph_el*p['th'])/3600   # Wh
    sig, dtip, croot = beam(x, W, p)
    return dict(S=S, cbar=cbar, croot=croot, m_wing=m_wing, mtot=mtot, W=W, CD0=CD0, e=e, K=K, D=D, Pcr=Pcr,
                LD=W/D, CL=W/(q*S), Vstall=Vstall, Ddl=Ddl, dl_frac=Ddl/W, Ph_el=Ph_el, E=E, sig=sig, dtip=dtip)

def g_nl(x, p=P):  # nonlinear constraints, g<=0, normalised
    r = model(x, p)
    return np.array([
        (1.2*r['Vstall'] - p['Vc'])/p['Vc'],                    # g1 stall margin
        (r['sig'] - p['sig_allow'])/p['sig_allow'],             # g2 root/spar stress
        (r['dtip'] - p['defl_frac']*x[0]/2)/(x[0]/2),           # g3 tip deflection
        (r['Ddl'] - p['dl_max']*r['W'])/(p['dl_max']*r['W']),   # g4 hover download
        (p['c_min'] - x[0]/x[1])/p['c_min'],                    # g14 min chord
    ])
NM = ['b', 'AR', 't/c', 'lambda', 'A_cap']
GN = ['g1 stall', 'g2 stress', 'g3 tip defl', 'g4 download', 'g14 min chord']

def bounds(p=P):
    lo = [1.5, p['AR_rng'][0], p['tcr_rng'][0], p['lam_rng'][0], p['Acap_rng'][0]]
    hi = [p['b_max'], p['AR_rng'][1], p['tcr_rng'][1], p['lam_rng'][1], p['Acap_rng'][1]]
    return np.array(lo), np.array(hi)

def solve(p=P, obj='E', nstart=25, seed=0, fixed=None):
    lo, hi = bounds(p); rng = np.random.default_rng(seed); best = None
    f = (lambda x: model(x, p)['E']) if obj == 'E' else (lambda x: model(x, p)['Pcr'])
    bnds = list(zip(lo/SC, hi/SC))
    if fixed:  # dict idx->value
        for i, v in fixed.items(): bnds[i] = (v/SC[i], v/SC[i])
    for _ in range(nstart):
        x0 = lo + rng.random(5)*(hi-lo)
        if fixed:
            for i, v in fixed.items(): x0[i] = v
        res = minimize(lambda z: f(z*SC)/100, x0/SC, method='SLSQP', bounds=bnds,
                       constraints=[{'type': 'ineq', 'fun': lambda z: -g_nl(z*SC, p)}],
                       options=dict(maxiter=300, ftol=1e-10))
        x = res.x*SC
        if res.success and g_nl(x, p).max() < 1e-6 and (best is None or f(x) < f(best)): best = x
    return best

def report(x, p=P, title=''):
    r = model(x, p); g = g_nl(x, p); lo, hi = bounds(p)
    print(f"--- {title}")
    print(f"x: b={x[0]:.3f} m  AR={x[1]:.2f}  t/c={x[2]:.3f}  lam={x[3]:.2f}  A_cap={x[4]*1e6:.1f} mm2")
    print(f"S={r['S']:.3f} m2  c_root={r['croot']:.3f} m  c_mean={r['cbar']:.3f} m  m_wing={r['m_wing']:.3f} kg  m_tot={r['mtot']:.2f} kg")
    print(f"CD0_wing={r['CD0']:.4f}  e={r['e']:.3f}  CL={r['CL']:.3f}  D={r['D']:.2f} N  L/D={r['LD']:.1f}  Pcr={r['Pcr']:.1f} W")
    print(f"Vstall={r['Vstall']:.2f} m/s  download={r['dl_frac']*100:.1f}% W  sigma={r['sig']/1e6:.0f} MPa  tip defl={r['dtip']*1e3:.0f} mm ({r['dtip']/(x[0]/2)*100:.1f}% semispan)")
    print(f"E_wing = {r['E']:.1f} Wh  (cruise+hover with mass feedback)")
    print("g:", {n: round(v, 3) for n, v in zip(GN, g)})
    act = [n for n, v in zip(GN, g) if v > -1e-3] + \
          [f"{NM[i]} at {'lower' if abs(x[i]-lo[i])/SC[i] < 1e-4 else 'upper'} bound" for i in range(5)
           if abs(x[i]-lo[i])/SC[i] < 1e-4 or abs(x[i]-hi[i])/SC[i] < 1e-4]
    print("active:", act, "feasible:", g.max() <= 1e-6)

if __name__ == '__main__':
    # 1) hand-picked feasible design point (for the report)
    x_f = np.array([2.4, 14.0, 0.12, 0.5, 25e-6])
    report(x_f, title='Feasible design point (hand-picked)')
    # 2) optimum
    xo = solve(); report(xo, title='SLSQP optimum (25 multistarts)')
    # 3) subsystem-only formulation (min Pcr, original objective) evaluated on the system metric
    xp = solve(obj='P'); report(xp, title='Min-Pcr only (original formulation)')
    # 4) sensitivities
    r0 = model(xo)['E']
    for k, d in [('m_other', 0.5), ('CLmax', 0.1), ('f_other', 0.003), ('sig_allow', 30e6), ('rho_skin', 0.05)]:
        p2 = dict(P); p2[k] = P[k] + d; x2 = solve(p2)
        print(f"sens {k:10s} +{d}: E {model(x2,p2)['E']-r0:+.2f} Wh, b={x2[0]:.2f} AR={x2[1]:.1f} t/c={x2[2]:.3f} lam={x2[3]:.2f}")
    # 5) plots: E vs AR for several spans (other variables re-optimised)
    plt.figure(figsize=(6, 4))
    for bb in [2.2, 2.6, 3.0]:
        ars = np.linspace(6, 14, 17); es = []
        for a in ars:
            xx = solve(nstart=6, fixed={0: bb, 1: a}); es.append(model(xx)['E'] if xx is not None else np.nan)
        plt.plot(ars, es, marker='o', ms=3, label=f'b = {bb} m')
    plt.xlabel('Aspect ratio'); plt.ylabel('E_wing [Wh]'); plt.legend(); plt.grid(alpha=.3)
    plt.title('Wing energy vs AR (t/c, taper, cap area re-optimised)'); plt.tight_layout()
    plt.savefig('/mnt/user-data/outputs/wing_E_vs_AR.png', dpi=160)
    # 6) interface values for teammates
    r = model(xo); print("\nINTERFACE -> Damian/Ashish:", {k: round(float(r[k]), 4) for k in ['S', 'm_wing', 'D', 'LD', 'dl_frac']}, 'CLmax', P['CLmax'])
