from src.rotor.parameters import Params
from src.rotor.constraints import constraint_margins
from optimization.rotor_slsqp import multistart

p = Params()
for blades in (2,3,4):
    best, _ = multistart(blades,p)
    margins = constraint_margins(best.x,blades,p)
    print(f"Nb={blades}: success={best.success}, energy={best.fun:.2f} Wh")
    print("  x =", best.x)
    print(f"  min constraint margin = {min(margins.values()):.5g}")
