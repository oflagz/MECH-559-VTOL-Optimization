import numpy as np
from scipy.optimize import minimize
from src.rotor.constraints import constraint_margins,scipy_constraints
from src.rotor.objective import mission_energy_wh

# theta75 in degrees; twist is fixed to zero until a separate twist variable is approved.
BOUNDS=[(0.15,0.30),(0.015,0.060),(150.,650.),(150.,750.),(5.,35.)]
GUESSES=[[0.25,0.030,300.,400.,20.],
         [0.20,0.040,400.,450.,15.],
         [0.28,0.025,350.,500.,25.],
         [0.18,0.050,500.,550.,20.],
         [0.27,0.045,450.,350.,30.]]

def run_one(Nb,x0,p):
    return minimize(mission_energy_wh,np.asarray(x0,float),args=(Nb,p),
                    method="SLSQP",bounds=BOUNDS,
                    constraints=scipy_constraints(Nb,p),
                    options={"ftol":1e-8,"maxiter":250,"disp":False})

def multistart(Nb,p):
    solutions=[run_one(Nb,x0,p) for x0 in GUESSES]
    feasible=[s for s in solutions if s.success and
              min(constraint_margins(s.x,Nb,p).values())>=-1e-5]
    if feasible:
        return min(feasible,key=lambda s:s.fun),solutions
    def violation(s):
        return sum(max(0.,-v)**2 for v in constraint_margins(s.x,Nb,p).values())
    return min(solutions,key=violation),solutions
