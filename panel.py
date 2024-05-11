import numpy as np
import matplotlib.pyplot as plt
from materials import *

def heat_gradient(T, dx, alpha, q_boundary_left, q_boundary_right, thermal_mass_density):
    """ 
    Returns an array of temperature gradients for an array of temperatures given a spacing dx and thermal diffusivity alpha.
    The boundaries are set as neuman boundary conditions (constant heat flux density)
    
    Arguments:
    T - Temperature array
    dx - Length between cells in m
    alpha - Thermal diffusivity (k/(rho*cp))

    Returns:
    dTdt - Heat gradient with respect to time
    """
    dTdt = np.zeros_like(T)
    N = len(T)

    if N > 2:
        # Solve heat conduction equation between adjacent cells for all inner cells.
        # Central differencing
        dTdt[1:N-1] = alpha * (T[:N-2] -2*T[1:N-1] + T[2:N]) / dx**2
    if N > 1:
        # Solve heat conduction equation between outer cell and neighbouring cell
        dTdt[0] += alpha * ((T[1] - T[0])/dx**2 ) + q_boundary_left/(thermal_mass_density*dx)
        dTdt[-1] += alpha * ((T[-2] - T[-1])/dx**2) + q_boundary_right/(thermal_mass_density*dx)
    else:
        # Single point only affected by heat flux
        dTdt = (q_boundary_left+q_boundary_right)/(thermal_mass_density*dx)
    return dTdt

class Panel:
    
    def __init__(self, panel=None):
        if panel == None:
            return
        self.set_shell_material(panel.shell_material, panel.shell_thickness, panel.shell_surface_h, panel.shell_surface_e)
        self.set_insulation_material(panel.insulation_material, panel.insulation_thickness)
        self.set_core_material(panel.core_material, panel.core_thickness, panel.core_surface_h)

        self.mesh(panel.N)
        self.set_temperature(panel.T)

    def set_shell_material(self, material, thickness, h=25, e=0.9):
        """
        Set the outer shell material of the panel with a certain thickness and surface conditions

        Arguments:
        material - Material which the shell is made from
        thickness - Thickness of the shell in m
        h - Outer surface convection coefficient
        e - Outer surface emissivity
        """
        self.shell_material = material
        self.shell_thickness = thickness
        self.shell_surface_h = h
        self.shell_surface_e = e

    def set_insulation_material(self, material, thickness):
        """
        Set the insulation material of the panel with a certain thickness

        Arguments:
        material - Material which the insulation is made from
        thickness - Thickness of the shell in m
        """
        self.insulation_material = material
        self.insulation_thickness = thickness

    def set_core_material(self, material, thickness, h=25):
        """
        Set the inner core material of the panel with a certain thickness and surface conditions

        Arguments:
        material - Material which the shell is made from
        thickness - Thickness of the shell in m
        h - Inner surface convection coefficient
        """
        self.core_material = material
        self.core_thickness = thickness
        self.core_surface_h = h

    def mesh(self, N):
        """
        Mesh the panel into a number of discrete points

        Arguments:
        N - number of discrete mesh points
        """
        self.N = N
        self.thickness = self.shell_thickness+self.insulation_thickness+self.core_thickness
        self.N_shell = int(self.shell_thickness/self.thickness * N) # Mesh shell as proportion of total thickness
        self.N_insulation = int(self.insulation_thickness/self.thickness * N) # Mesh insulation as proportion of total thickness
        self.N_core = int(self.core_thickness/self.thickness * N) # Mesh core as proportion of total thickness
        self.dx = self.thickness/N # Determine length of a single mesh cell
        self.x = np.arange(0,self.thickness,self.dx) # Create array of x-coordinates
        self.T = np.zeros(N) # Create temperature mesh

    def set_temperature(self, T):
        """
        Set the temperature of all mesh cells.

        Arguments:
        T - Temperature
        """
        self.T = np.zeros(self.N) + T

    def step_temperature(self, dt, T_infinity=45, T_internal=45, q_solar=1400):
        """
        Update the internal temperatures of points inside the panel.

        Arguments:
        dt - Time step in seconds
        T_infinity - Ambient temperature
        T_internal - Temperature inside nest
        q_solar - Heat flux density at outer surface due to radiation

        Returns:
        q_conv_bot - Heat flux density from the panel to the interior
        """
        T = self.T
        dTdt = np.zeros(self.N)
        
        # Bottom boundary temp change due to convection
        T_bottom = T[-1]
        h_bot = self.core_surface_h
        # Convection at bottom surface
        q_conv_bot = h_bot*(T_bottom-T_internal)

        # Fourth order Runge-Kutta
        k1 = self.rk_temperature_step(T, T_infinity, T_internal, q_solar)
        k2 = self.rk_temperature_step(T+dt*k1/2, T_infinity, T_internal, q_solar)
        k3 = self.rk_temperature_step(T+dt*k2/2, T_infinity, T_internal, q_solar)
        k4 = self.rk_temperature_step(T+dt*k3, T_infinity, T_internal, q_solar)

        dTdt = 1/6 * (k1 + 2*k2 + 2*k3 + k4)
        #dTdt = k1
        T += dTdt*dt 
        # Update internal temperature
        self.T = T
        # Return total heat transfer to inside of the nest
        return q_conv_bot
    
    def rk_temperature_step(self, T, T_infinity=45, T_internal=45, q_solar=1400):
        dTdt = np.zeros(self.N)
        
        # Calculate alpha values for each material
        # alpha = k/(rho*cp)
        alpha_shell = self.shell_material.k/(self.shell_material.rho*self.shell_material.cp)
        alpha_insulation = self.insulation_material.k/(self.insulation_material.rho*self.insulation_material.cp)
        alpha_core = self.core_material.k/(self.core_material.rho*self.core_material.cp)
        
        # Top boundary temp change due to radiation and convection
        T_top = T[0]
        h_top = self.shell_surface_h
        # Convection at top surface
        q_conv_top = h_top*(T_top-T_infinity)
        # Radiation at top surface
        q_rad_absorbed = q_solar*self.shell_surface_e
        q_rad_emitted = 5.67e-8 * self.shell_surface_e * ((T_top+273.15)**4 - (T_infinity+273.15)**4)
        q_rad_top = q_rad_absorbed - q_rad_emitted
        q_top = q_rad_top-q_conv_top

        # Bottom boundary temp change due to convection
        T_bottom = T[-1]
        h_bot = self.core_surface_h
        # Convection at bottom surface
        q_conv_bot = h_bot*(T_bottom-T_internal)
        q_bot = -q_conv_bot
        
        # Shell heat gradient
        k_shell_ins_boundary = 0.5 * (self.shell_material.k+self.insulation_material.k)  # Average conductivity between shell and insulation
        q_shell_ins = k_shell_ins_boundary*(T[self.N_shell-1] - T[self.N_shell])/self.dx # Heat transfer density from insulation to shell
        # Calculate heat gradient with boundary conditions due to heat flux at top and heat flux at insulation boundary
        dTdt[0:self.N_shell] += heat_gradient(T[0:self.N_shell], self.dx, alpha_shell,q_top,-q_shell_ins,self.shell_material.rho*self.shell_material.cp)
        
        # Insulation Indices
        i_ins_start = self.N_shell
        i_ins_end = i_ins_start + self.N_insulation

        # Insulation heat gradient
        k_ins_core_boundary = 0.5 * (self.core_material.k+self.insulation_material.k) # Average conductivity between insulation and core
        q_ins_core = k_ins_core_boundary*(T[i_ins_end-1]-T[i_ins_end])/self.dx        # Heat transfer density from core to insulation
        # Calculate heat gradient with boundary conditions due to heat flux at shell and core boundaries
        dTdt[i_ins_start:i_ins_end] += heat_gradient(T[i_ins_start:i_ins_end], self.dx, alpha_insulation, q_shell_ins, -q_ins_core,self.insulation_material.rho*self.insulation_material.cp)

        # Core heat gradient
        # Calculate heat gradient with boundary conditions due to heat flux at bottom surface and core boundary
        dTdt[i_ins_end:] += heat_gradient(T[i_ins_end:], self.dx, alpha_core,q_ins_core,q_bot,self.core_material.rho*self.core_material.cp)
        return dTdt