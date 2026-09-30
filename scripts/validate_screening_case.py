from src.rotor.parameters import Params
from src.rotor.hover import hover_model

p = Params()
h = hover_model(0.25, 0.030, 300.0, 3, p)

print("Proposal rotor screening case")
print(f"Required thrust/rotor: {h['T_req']:.3f} N")
print(f"Solidity:              {h['sigma']:.4f}")
print(f"Ideal induced power:   {h['Pi']:.2f} W/rotor")
print(f"Profile power:         {h['P0']:.2f} W/rotor")
print(f"Shaft power:           {h['P_sh']:.2f} W/rotor")
print(f"Figure of merit:       {h['FM']:.3f}")
print(f"Total hover electrical:{p.Nr*h['P_el']:.1f} W")
