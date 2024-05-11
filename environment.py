import numpy as np
from scipy import interpolate

class Radiation:
    def __init__(self,angle,intensity):
        """
        Radiation with horizontal angle and intensity.

        Arguments:
        angle - angle with the horizontal axis
        intensity - 
        """
        self.angle = angle
        self.intensity = intensity

class AmbientTemperature:
    def __init__(self, times, temperatures):
        """
        Create an ambient temperature profile from an array of times and corresponding temperatures.

        Arguments:
        times - N-dimensional array of times
        temperatures - N-dimensional array of temperatures at corresponding times
        """
        self.times = times
        self.temperatures = temperatures
        # Interpolate temperature over time
        self.f = interpolate.splrep(times, temperatures, s=25)

    def get_temperature_at_hour(self, time):
        """
        Return interpolated temperature at time.

        Arguments:
        time - time of day in hours

        Returns:
        T - temperature in °C
        """
        return interpolate.BSpline(*self.f)(time)

class Solar:
    def __init__(self, peak_intensity):
        """
        Create a solar profile with a certain peak intensity.

        Arguments:
        peak_intensity - peak solar intensity at the equator
        """
        self.peak_intensity = peak_intensity
    
    def get_radiation_at_hour(self, hour, latitude=0):
        """
        Get the radiation (intensity and angle) of the solar profile at a given time and latitude.

        Arguments:
        hour - hour of the day (0-24)
        latitude - angle of latitude of the location to measure intensity
        """
        # Calculate the horizontal angle of the sun
        angle = hour/24 * 360 - 90

        # If daytime calculate intensity
        if angle >=0 and angle <= 180:
            # Sinusoidal intensity dependant on time and latitude
            intensity = self.peak_intensity * np.sin(2*np.pi*(hour-6)/24) * np.cos(np.radians(latitude))
        
        # Otherwise intensity = 0
        else:
            intensity = 0

        return Radiation(angle, intensity)

class Environment:
    def __init__(self, solar, temperature, latitude=0, time=0):
        """
        Create an environment with a solar and ambient temperature profile at a given latitude.

        Arguments:
        solar - solar profile to use
        temperature - ambient temperature profile to use
        latitude - location angle of latitude
        """
        self.solar = solar
        self.temperature = temperature
        self.latitude = latitude
        self.time = time
    
    def set_time_hours(self, time):
        """
        Set time for all environmental calculation in hours.

        Arguments:
        time - time in hours
        """
        self.time = time
    
    def set_time_seconds(self, time):
        """
        Set time for all environmental calculation in seconds.

        Arguments:
        time - time in seconds
        """
        self.time = time/3600
    
    def get_ambient_temperature(self):
        """
        Get ambient temperature at current time.

        Returns:
        temperature - ambient temperature at current time
        """
        return self.temperature.get_temperature_at_hour(self.time)
    
    def get_radiation(self, heading=0, attitude=90):
        """
        Get radiaion passing through a surface normal at heading and attitude.

        Arguments:
        heading - compass heading
        attitude - vertical angle with the surface

        Returns:
        absorbed - radiation absorbed by surface
        """
        radiation = self.solar.get_radiation_at_hour(self.time, self.latitude)
        
        # Convert heading and attitude to radians for trigonometry
        heading_rad = np.radians(heading)
        attitude_rad = np.radians(attitude)

        # Calculate x y and z components of the surface normal
        # x - east/west
        # y - north/south
        # z - up/down
        surface_x = np.sin(heading_rad) * np.cos(attitude_rad)
        surface_y = np.cos(heading_rad) * np.cos(attitude_rad)
        surface_z = np.sin(attitude_rad)

        # Convert latitude and radiation angle to radians for trigonometry
        latitude_rad = np.radians(self.latitude)
        angle_rad = np.radians(radiation.angle)
        
        # Calculate x y and z components of radiation
        radiation_x = radiation.intensity * np.cos(angle_rad)
        radiation_y = -radiation.intensity * np.sin(latitude_rad)*np.sin(angle_rad)
        radiation_x_squared = radiation_x**2
        radiation_y_squared = radiation_y**2
        radiation_z_squared = radiation.intensity **2 - radiation_x_squared - radiation_y_squared
        
        # Correct rounding errors
        if radiation_z_squared < 0:
            radiation_z=0
        else:
            radiation_z = np.sqrt(radiation_z_squared)
        
        # Absorbed radiation is dot product of components of radiation and surface normal
        absorbed = surface_x * radiation_x + surface_y * radiation_y + surface_z * radiation_z

        # If radiation flux is negative then no radiation is absorbed
        if absorbed < 0:
            absorbed = 0
        return absorbed