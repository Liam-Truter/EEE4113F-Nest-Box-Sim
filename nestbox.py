import numpy as np
from panel import Panel
from materials import Material
from environment import *

class NestBox:
    def __init__(self, internal_dimensions=(0.17,0.21,0.46), hole_size=0.055, panel=None, shaded_panel=None, T_init=30, environment=None, orientation=180, initial_time=0, T_thresh=35, peltier=0):
        """
        Create a nestbox environment with given parameters.

        Arguments:
        internal_dimensions - dimensions in meters of the internal volume (width, depth, height)
        hole_size - diameter of entry hole in meters
        panel - default panel for all sides exposed to direct solar radiation
        shaded_panel - default panel for all sides not exposed to direct solar radiation
        T_init - inital temperature of all internal volumes
        environment - the environment class detailing the solar and temperature profile as well as latitude of the location
        orientation - heading in degrees
        initial_time - time to start simulation in hours
        T_thresh - threshold temperature in °C to activate the Peltier
        peltier - cooling capacity of the Peltier device in W
        """
        self.internal_dimensions=internal_dimensions
        self.hole_size = hole_size
        self.air = Material(1.146, 1005, 0.0262) # Properties of air
        self.bird = Material(999.8, 4186, 0.6)   # Properties of water (good estimation for bird)
        self.bird_volume = 0.001                 # Assume birds take up 1l of space combined
        self.bird_metabolism = 15                # Assume birds collectively generate 10 W of heat
        
        # Get internal volume from internal dimensions - the volume displaced by birds
        self.air_volume = internal_dimensions[0] * internal_dimensions[1] * internal_dimensions[2] - self.bird_volume

        # Calculate the combined thermal mass of air and bird volumes
        self.thermal_mass = self.air.rho * self.air.cp * self.air_volume + self.bird.rho * self.bird.cp * self.bird_volume

        # Reasonable default material for panels is oak
        oak = Material(545, 2.385, 0.17)

        if panel == None:
            # If no panel is set, create sensible default
            air = self.air
            panel = Panel()
            # 2mm shell, 2cm air gap insulation, 2cm solid core 
            panel.set_shell_material(oak, 0.002)
            panel.set_insulation_material(air, 0.02)
            panel.set_core_material(oak, 0.2)
            panel.mesh(21)
        
        if shaded_panel == None:
            # If no shaded panel is set, create sensible defualt
            shaded_panel = Panel()

            # 2cm solid panel
            shaded_panel.set_shell_material(oak, 0.005)
            shaded_panel.set_insulation_material(oak, 0.005)
            shaded_panel.set_core_material(oak, 0.01)
            shaded_panel.mesh(4)
        
        # Initialize temperatures of panels
        panel.set_temperature(T_init)
        shaded_panel.set_temperature(T_init)
        
        # Set panels
        self.panel_front = Panel(shaded_panel)  # Should always be in shaded side anyway
        self.panel_top = Panel(panel)
        self.panel_left = Panel(panel)
        self.panel_right = Panel(panel)
        self.panel_shaded = Panel(shaded_panel) # Bottom and back panels combined for simplicity

        # Initialize internal temperature
        self.temperature = T_init

        # Create sensible default environment if none is set
        if environment == None:
            temperatures = AmbientTemperature([0, 24], [T_init, T_init]) # Constantly at initial temperature
            sun = Solar(1400)
            environment = Environment(sun, temperatures, -25.6)          # Latitude of the Kalahari

        self.environment = environment

        self.T_thresh = T_thresh
        self.peltier = peltier
        self.cooling_energy_use = 0
        self.time = initial_time
        self.environment.set_time_hours(initial_time)
        self.orientation = orientation

    def step_time(self, dt):
        """
        Increment the nest box simulation time in seconds.

        Arguments:
        dt - time increment in seconds
        """
        dt_hours = dt/3600          # Convert time increment to hours
        self.time += dt_hours       # Increment hour time of day
        self.time = self.time % 24  # If passing midnight, reset to zero
        self.environment.set_time_hours(self.time)
    
    def step_temperature(self, dt):
        """
        Perform an update on internal temperatures of the nestbox.

        Arguments:
        dt - time in seconds to simulate

        Returns:
        temperature - internal temperature of the nestbox after simulation step
        """
        # Solar flux density through each panel
        q_solar_front = self.environment.get_radiation(self.orientation, 0)
        q_solar_top = self.environment.get_radiation(self.orientation, 90)
        q_solar_left = self.environment.get_radiation(self.orientation - 90, 0)
        q_solar_right = self.environment.get_radiation(self.orientation + 90, 0)

        # Update environmental time
        self.environment.set_time_hours(self.time)
        # Get ambient outside temperature
        T_ambient = self.environment.get_ambient_temperature()

        # Simulate panel thermals, save the flux density into the inside of the nest
        q_front = self.panel_front.step_temperature(dt, T_ambient, self.temperature, q_solar_front)
        q_top = self.panel_front.step_temperature(dt, T_ambient, self.temperature, q_solar_top)
        q_left = self.panel_front.step_temperature(dt, T_ambient, self.temperature, q_solar_left)
        q_right = self.panel_front.step_temperature(dt, T_ambient, self.temperature, q_solar_right)
        q_shaded = self.panel_shaded.step_temperature(dt, T_ambient, self.temperature, 0)

        # Get internal dimensions
        width = self.internal_dimensions[0]
        length = self.internal_dimensions[1]
        height = self.internal_dimensions[2]

        # Calculate internal area of each panel
        A_frontback = width * height
        A_topbottom = width * length
        A_leftright = length * height
        A_shaded = A_frontback + A_topbottom

        # Convert heat flux density into heat flux from area
        Q_front = q_front * A_frontback
        Q_top = q_top * A_topbottom
        Q_left = q_left * A_leftright
        Q_right = q_right * A_leftright
        Q_shaded = q_shaded * A_shaded

        # Calculate area of the hole in m^2
        A_hole = np.pi * (0.5*self.hole_size)**2
        # Calculate natural convection through entry hole
        Q_convection = 20 * A_hole * (T_ambient - self.temperature)

        # Peltier off by default
        Q_cooling = 0
        # Activate peltier if above threshold (bang-bang) and record energy used
        if self.temperature > self.T_thresh:
            Q_cooling = -self.peltier
            self.cooling_energy_use += self.peltier * dt

        # Calculate net heat flux of the entire nestbox system into the internal volume
        Qnet = Q_front + Q_top + Q_left + Q_right + Q_convection + Q_shaded + self.bird_metabolism + Q_cooling
        # Calculate change in temperature due to heat fluxes
        dT = (Qnet * dt)/self.thermal_mass

        # Increment simulation time
        self.step_time(dt)
        # Update internal temperature
        self.temperature += dT * dt
        return self.temperature
    
    def getTotalEnergyUse(self):
        """
        Get the energy used by the system during simulation time.

        Returns:
        energy - energy used to cool the system in J
        """
        return self.cooling_energy_use
        