import numpy as np
import matplotlib.pyplot as plt
import csv

from materials import Material
from panel import Panel

surface_e = 1

surface_h_outer = 25
surface_h_inner = 25

q_outer = 0

T_init = 35
T_outer = 45
T_inner = 25

# Available databases to test materials
databases = ['NCFS UCF Materials.csv', 'ORNL Materials.csv']
materials = []
for database in databases:
    # Open each database
    with open(database,'r') as csvfile:
        csvreader = csv.reader(csvfile)
        # Iterate through each row excluding header
        csvreader.__next__()
        for row in csvreader:
            name = row[0]       # Name is a string
            k = float(row[1])   # Conductivity is a float
            rho = int(row[2])   # Density is an integer
            cp = float(row[3])  # Heat Capacity is a float

            # Append material properties to list of materials
            materials.append([name, k, rho, cp])

# 1 hour simulation with timestep of 10 ms
dt = 0.02
t = np.arange(0,3600, dt)

material_name_len = max(len(material[0]) for material in materials)

for i in range(len(materials)):
    row = materials[i]
    name=row[0]
    k=row[1]
    rho=row[2]
    cp=row[3]
    material = Material(rho,cp,k)
    # 2 cm thick solid panel of given material
    panel = Panel()
    panel.set_shell_material(material,0.005, surface_h_outer, surface_e)
    panel.set_insulation_material(material, 0.005)
    panel.set_core_material(material,0.01,surface_h_inner)
    panel.mesh(8)
    panel.set_temperature(T_init)

    # Simulate temperature over time and get heat transfer
    for i in range(len(t)):
        q_inner = panel.step_temperature(dt,T_outer,T_inner, q_outer)
    
    # Calculate thermal resistance
    T_surface_outer = panel.T[0]
    T_surface_inner = panel.T[-1]
    delta_T = T_surface_outer-T_surface_inner
    thermal_resistance = delta_T / q_inner # Assuming panel is 1 m^2
    # Format all materials to same width
    name += ':'
    formatted_name = '{:<{width}}'.format(name, width=material_name_len+1)
    print(f"{formatted_name} {thermal_resistance:.3f} °C/W")