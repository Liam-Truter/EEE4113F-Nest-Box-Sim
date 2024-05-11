import numpy as np
import matplotlib.pyplot as plt
import csv

from materials import Material
from panel import Panel

surface_e = 1

surface_h_outer = 25
surface_h_inner = 25

T_init = 35
T_outer = 45
T_inner = 25
q_outer = 0

pine = Material(640,2.805, 0.15)
air = Material(1.146, 1005, 0.0262)
polystyrene = Material(12,1450,0.03936)

# Panel with air insulation
panelA = Panel()
panelA.set_shell_material(pine,0.004, surface_h_outer, surface_e)
panelA.set_insulation_material(air, 0.02)
panelA.set_core_material(pine, 0.02, surface_h_inner)
panelA.mesh(22)
panelA.set_temperature(T_init)

# Panel with polystyrene insulation
panelB = Panel()
panelB.set_shell_material(pine, 0.004, surface_h_outer, surface_e)
panelB.set_insulation_material(polystyrene, 0.02)
panelB.set_core_material(pine, 0.02, surface_h_inner)
panelB.mesh(22)
panelB.set_temperature(T_init)


# 1 hour simulation with timestep of 10 ms
dt = 0.01
t = np.arange(0,3600, dt)

for i in range(len(t)):
    q_inner_A = panelA.step_temperature(dt, T_outer, T_inner, q_outer)
    q_inner_B = panelB.step_temperature(dt, T_outer, T_inner, q_outer)
    if t[i] == 3500:
        print("pause")

# Thermal resistance of polystyrene panel
T_outer_A = panelA.T[0]
T_inner_A = panelA.T[-1]

delta_T_A = T_outer_A-T_inner_A

thermal_resistance_A = delta_T_A / q_inner_A # Asumming panel is 1 m^2

print(f"Air:\t\t {thermal_resistance_A:.3f} °C/W")

# Thermal resistance of air panel
T_outer_B = panelB.T[0]
T_inner_B = panelB.T[-1]

delta_T_B = T_outer_B-T_inner_B

thermal_resistance_B = delta_T_B / q_inner_B # Asumming panel is 1 m^2

print(f"Polystyrene:\t {thermal_resistance_B:.3f} °C/W")