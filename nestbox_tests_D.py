import numpy as np
import matplotlib.pyplot as plt
from nestbox import *

start_time = 0
stop_time = 24
dt = 0.1
dt_hours = dt/3600
t = np.arange(start_time, stop_time, dt_hours)
t_save = np.arange(start_time, stop_time, 1/6)

sun = Solar(1400)
temperature_data = np.loadtxt('Temperatures DEC (increased).csv', delimiter=',', skiprows=1)
ambient_temperature_times = temperature_data[:,0]
ambient_temperature_values = temperature_data[:,1]
ambient_temperatures = AmbientTemperature(ambient_temperature_times, ambient_temperature_values)
environment = Environment(sun, ambient_temperatures, -25.6)

air = Material(1.146, 1005, 0.0262)
oak = Material(545, 2.385, 0.17)
pine = Material(640,2.805, 0.15)
hd_concrete = Material(2800,840, 1.4)

paint_e_dark = 0.9
paint_e_light = 0.06
surface_outer_h = 25
surface_inner_h = 25

T_init = 30
T_init_A = 31
T_init_B = 31
T_init_C = 31
T_init_D = 31

shaded_panel = Panel()
shaded_panel.set_shell_material(pine, 0.005, surface_outer_h, paint_e_light)
shaded_panel.set_insulation_material(pine, 0.005)
shaded_panel.set_core_material(pine, 0.01, surface_inner_h)
shaded_panel.mesh(4)
shaded_panel.set_temperature(T_init)
shaded_panel_dark = Panel()
shaded_panel_dark.set_shell_material(pine, 0.005, surface_outer_h, paint_e_dark)
shaded_panel_dark.set_insulation_material(pine, 0.005)
shaded_panel_dark.set_core_material(pine, 0.01, surface_inner_h)
shaded_panel_dark.mesh(4)
shaded_panel.set_temperature(T_init)

panel = Panel()
panel.set_shell_material(pine, 0.002, surface_outer_h, paint_e_light)
panel.set_insulation_material(air, 0.02)
panel.set_core_material(pine, 0.02, surface_inner_h)
panel.mesh(21)
panel.set_temperature(T_init)
panel_dark = Panel()
panel_dark.set_shell_material(pine, 0.002, surface_outer_h, paint_e_light)
panel_dark.set_insulation_material(air, 0.02)
panel_dark.set_core_material(pine, 0.02, surface_inner_h)
panel_dark.mesh(21)
panel_dark.set_temperature(T_init)

concrete_panel = Panel()
concrete_panel.set_shell_material(pine, 0.005, surface_outer_h, paint_e_light)
concrete_panel.set_insulation_material(hd_concrete, 0.03)
concrete_panel.set_core_material(pine, 0.005, surface_inner_h)
concrete_panel.mesh(8)
concrete_panel.set_temperature(T_init)

thick_concrete_panel = Panel()
thick_concrete_panel.set_shell_material(pine, 0.005, surface_outer_h, paint_e_light)
thick_concrete_panel.set_insulation_material(hd_concrete, 0.2)
thick_concrete_panel.set_core_material(pine, 0.005, surface_inner_h)
thick_concrete_panel.mesh(42)
thick_concrete_panel.set_temperature(T_init)

boxA = NestBox(initial_time=start_time,panel=panel, shaded_panel=shaded_panel, environment=environment,orientation=180, T_init=T_init_A, peltier=5)
boxB = NestBox(initial_time=start_time,panel=panel, shaded_panel=shaded_panel, environment=environment,orientation=180, T_init=T_init_B, peltier=15)
boxC = NestBox(initial_time=start_time,panel=panel, shaded_panel=shaded_panel, environment=environment,orientation=180, T_init=T_init_C, peltier=25)
boxD = NestBox(initial_time=start_time,panel=panel, shaded_panel=shaded_panel, environment=environment,orientation=180, T_init=T_init_D, peltier=50)

# Create an empty plot
plt.figure(figsize=(12, 8))
plt.xlabel('Time (hours)')
plt.ylabel('Temperature')
plt.title(f'Internal and External Temperature of Nest Boxes With Various Peltier Cooling Levels Throughout the Day')
plt.grid(True)

lineA, = plt.plot([], [], label="External")
lineBA, = plt.plot([], [], label="5W")
lineBB, = plt.plot([], [], label="15W")
lineBC, = plt.plot([], [], label="25W")
lineBD, = plt.plot([], [], label="50W")

plt.legend()

steps = len(t)
steps_save = len(t_save)
steps_per_save = steps // steps_save

internal_temperatures_A = np.zeros(steps_save)
internal_temperatures_B = np.zeros(steps_save)
internal_temperatures_C = np.zeros(steps_save)
internal_temperatures_D = np.zeros(steps_save)
external_temperatures = np.zeros(steps_save)

internal_temperatures_A[0] = T_init_A
internal_temperatures_B[0] = T_init_B
internal_temperatures_C[0] = T_init_C
internal_temperatures_D[0] = T_init_D
external_temperatures[0] = environment.get_ambient_temperature()
for step in range(1,steps_save):
    for substep in range(steps_per_save):
        boxA.step_temperature(dt)
        boxB.step_temperature(dt)
        boxC.step_temperature(dt)
        boxD.step_temperature(dt)
    
    internal_temperatures_A[step] = boxA.step_temperature(dt)
    internal_temperatures_B[step] = boxB.step_temperature(dt)
    internal_temperatures_C[step] = boxC.step_temperature(dt)
    internal_temperatures_D[step] = boxD.step_temperature(dt)
    external_temperatures[step] = boxA.environment.get_ambient_temperature()

    lineBA.set_data(t_save[:step], internal_temperatures_A[:step])
    lineBB.set_data(t_save[:step], internal_temperatures_B[:step])
    lineBC.set_data(t_save[:step], internal_temperatures_C[:step])
    lineBD.set_data(t_save[:step], internal_temperatures_D[:step])
    lineA.set_data(t_save[:step], external_temperatures[:step])

    lower_bound = min(np.min(internal_temperatures_A[:step]), np.min(internal_temperatures_B[:step]),\
                      np.min(internal_temperatures_C[:step]), np.min(internal_temperatures_D[:step]),\
                      np.min(external_temperatures[:step])) - 5
    upper_bound = max(np.max(internal_temperatures_A[:step]), np.max(internal_temperatures_B[:step]),\
                      np.max(internal_temperatures_C[:step]), np.max(internal_temperatures_D[:step]),\
                      np.max(external_temperatures[:step])) + 5

    plt.xlim(start_time, stop_time)
    plt.ylim(lower_bound,upper_bound)
    plt.pause(0.1)

print(f"5W:  {boxA.getTotalEnergyUse()/1000} kJ")
print(f"15W: {boxB.getTotalEnergyUse()/1000} kJ")
print(f"25W: {boxC.getTotalEnergyUse()/1000} kJ")
print(f"50W: {boxD.getTotalEnergyUse()/1000} kJ")

plt.show()