"""MECH 559 Sprint 3 - wing subsystem analysis (Antoine). Imports wing_model.py (Sprint 2 model).
1) algorithm behaviour: multistart spread  2) KKT / active-set check  3) parameter sensitivities
4) 3x3 aircraft-mass / non-wing-drag scenarios  5) beam discretisation check  6) interface for propulsion"""
import numpy as np
from scipy.optimize import minimize, nnls
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import wing_model as wm
from wing_model import P, SC, model, g_nl, bounds, solve, GN, NM

def run_start(x0, p):
    lo, hi = bounds(p)
    res = minimize(lambda z: model(z*SC, p)['E']/100, x0/SC, method='SLSQP', bounds=list(zip(lo/SC, hi/SC)),
                   constraints=[{'type': 'ineq', 'fun': lambda z: -g_nl(z*SC, p)}], options=dict(maxiter=300, ftol=1e-10))
    x = res.x*SC
    ok = res.success and g_nl(x, p).max() < 1e-6
    return x, ok, res.nit

def grad(f, z, h=1e-5):
    g = np.zeros_like(z)
    for i in range(len(z)):
        d = np.zeros_like(z); d[i] = h*max(1, abs(z[i]))
        g[i] = (f(z+d)-f(z-d))/(2*d[i])
    return g

def kkt(x, p=P, tol=1e-3):
    z = x/SC; lo, hi = bounds(p)
    f = lambda zz: model(zz*SC, p)['E']/100
    gf = grad(f, z)
    cols, names = [], []
    gv = g_nl(x, p)
    for j, nm in enumerate(GN):
        if gv[j] > -tol:
            cols.append(grad(lambda zz, j=j: g_nl(zz*SC, p)[j], z)); names.append(nm)
    for i in range(5):
        e = np.zeros(5); e[i] = 1.0
        if abs(x[i]-lo[i])/SC[i] < 1e-4: cols.append(-e); names.append(f'{NM[i]} >= lower')
        if abs(x[i]-hi[i])/SC[i] < 1e-4: cols.append(e); names.append(f'{NM[i]} <= upper')
    A = np.array(cols).T
    lam, rn = nnls(A, -gf)
    return names, lam, rn/np.linalg.norm(gf), np.linalg.matrix_rank(A, tol=1e-6), A.shape[1]

if __name__ == '__main__':
    # ---- 1) multistart behaviour
    rng = np.random.default_rng(1); lo, hi = bounds(P)
    runs = [run_start(lo + rng.random(5)*(hi-lo), P) for _ in range(60)]
    Es = np.array([model(x)['E'] if ok else np.nan for x, ok, _ in runs])
    nok = np.isfinite(Es).sum(); Eb = np.nanmin(Es)
    same = np.sum(np.abs(Es-Eb)/Eb < 1e-3)
    print(f"[1] multistart: {nok}/60 feasible convergences, {same}/60 within 0.1% of best E; "
          f"E range {np.nanmin(Es):.2f}-{np.nanmax(Es):.2f} Wh; mean iters {np.mean([r[2] for r in runs]):.0f}")
    xo = solve(P, nstart=40, seed=3); ro = model(xo)
    print("optimum x =", np.round(xo, 6), " E =", round(ro['E'], 2))
    wm.report(xo, title='Sprint 3 optimum')
    # ---- 2) KKT
    names, lam, relres, rk, nact = kkt(xo)
    print(f"[2] KKT: active set {names}")
    print("    multipliers (scaled f=E/100):", np.round(lam, 4), f" stationarity residual/|grad f| = {relres:.2e}",
          f" rank(A)={rk}/{nact}")
    # ---- 3) sensitivities (+10 % on each parameter, re-optimised)
    sens = {}
    for k in ['m_other', 'f_other', 'rho_skin', 'k_struct', 'E_mod', 'sig_allow', 'CLmax', 'CD_dl', 'FM']:
        p2 = dict(P); p2[k] = P[k]*1.10; x2 = solve(p2, nstart=12, seed=5)
        sens[k] = (model(x2, p2)['E']-ro['E'])/ro['E']/0.10
        print(f"[3] {k:10s} +10%: dlnE/dlnp = {sens[k]:+.3f}   x*: b={x2[0]:.2f} AR={x2[1]:.1f} t/c={x2[2]:.3f} lam={x2[3]:.2f} A={x2[4]*1e6:.1f}")
    # ---- 4) scenarios (non-wing mass x non-wing drag area)
    print("[4] scenarios  m_other[kg] f_oth[m2] -> E[Wh], L/D, m_wing[kg], b[m]")
    scen = []
    for m in [6.0, 7.0, 8.0]:
        for f in [0.010, 0.0134, 0.017]:
            p2 = dict(P); p2['m_other'] = m; p2['f_other'] = f; x2 = solve(p2, nstart=10, seed=7); r = model(x2, p2)
            scen.append(r['E']); print(f"    {m:4.1f} {f:.4f} -> {r['E']:6.1f} {r['LD']:4.1f} {r['m_wing']:.2f} {x2[0]:.2f}")
    print(f"    scenario energy range {min(scen):.1f}-{max(scen):.1f} Wh")
    # ---- 5) beam discretisation check
    from wing_model import beam
    d1 = beam(xo, ro['W'], P, 101); d2 = beam(xo, ro['W'], P, 401)
    print(f"[5] beam n=101 vs 401: tip deflection diff {abs(d1[1]-d2[1])/d2[1]*100:.3f}%, max stress diff {abs(d1[0]-d2[0])/d2[0]*100:.3f}%")
    # ---- 6) interface to propulsion/rotor
    kst = 1.2
    print(f"[6] interface: S={ro['S']:.3f} m2  m_wing={ro['m_wing']:.3f} kg  D={ro['D']:.2f} N  L/D={ro['LD']:.2f}  CLmax={P['CLmax']}  "
          f"Vstall={ro['Vstall']:.2f} m/s  Vtrans_min(k=1.2)={kst*ro['Vstall']:.2f} m/s  download={ro['dl_frac']*100:.2f}% W")
    print(f"    L/D vs propulsion placeholder 7.01: {(ro['LD']/7.01-1)*100:+.1f}%")
    # ---- plot: tornado of sensitivities
    ks = sorted(sens, key=lambda k: abs(sens[k]))
    plt.figure(figsize=(5.5, 3.3)); plt.barh(ks, [sens[k] for k in ks], color='tab:blue')
    plt.xlabel(r'$d\ln E^*_{wing}/d\ln p$ (+10% on parameter, re-optimised)'); plt.grid(axis='x', alpha=.3); plt.tight_layout()
    plt.savefig('/mnt/user-data/outputs/wing_sensitivity.png', dpi=170)
