import numpy as np
from scipy.optimize import minimize
from src.rotor.constraints import constraint_margins, scipy_constraints
from src.rotor.objective import mission_energy_wh

BOUNDS = [
    (0.15, 0.30),
    (0.015, 0.060),
    (150.0, 650.0),
    (150.0, 750.0),
    (-25.0, 10.0),
]

GUESSES = [
    [0.25,0.030,300.,400.,-10.],
    [0.20,0.040,400.,450.,-5.],
    [0.28,0.025,350.,500.,-15.],
    [0.18,0.050,500.,550.,0.],
    [0.27,0.045,450.,350.,-10.],
]

def run_one(Nb, x0, p):
    return minimize(
        mission_energy_wh, np.asarray(x0,float), args=(Nb,p), method="SLSQP",
        bounds=BOUNDS, constraints=scipy_constraints(Nb,p),
        options={"ftol":1e-9,"maxiter":500,"disp":False},
    )

def multistart(Nb, p):
    solutions = [run_one(Nb, x0, p) for x0 in GUESSES]
    feasible = [
        s for s in solutions
        if s.success and min(constraint_margins(s.x,Nb,p).values()) >= -1e-5
    ]
    if feasible:
        return min(feasible, key=lambda s:s.fun), solutions

    def violation(s):
        return sum(max(0.0,-v)**2 for v in constraint_margins(s.x,Nb,p).values())
    return min(solutions,key=violation), solutions
