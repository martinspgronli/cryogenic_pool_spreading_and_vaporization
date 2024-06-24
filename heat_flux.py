"""
This Python script calculates the heat flux from the ground and air to the pool.
Calculation of ground heat flux partly assumes perfect thermal contact between ground surface and pool, see https://doi.org/10.1016/j.ijhydene.2020.06.131


Author: Martin Spillum Grønli (2024, SINTEF Energy Research)
"""


# Import
import numpy as np
import pandas as pd


# Heat transfer from the ground partly assuming perfect contact between ground surface and pool, https://doi.org/10.1016/j.ijhydene.2020.06.131
def heat_flux_ground_partly_perf(state, var):     
    t_w = state.problem_data['time_wet']    # Time (s) each ground grid cell has been wet
    heat_flux_ground = np.where((t_w >= 0) & (t_w < 4), var.thermal_cond_ground * (var.T_ground_inf - var.T_pool) * (1.5 - 0.25 * t_w) / np.sqrt(np.pi * var.thermal_diff_ground), 0) 
    heat_flux_ground += np.where(t_w >= 4, var.thermal_cond_ground * (var.T_ground_inf - var.T_pool) / np.sqrt(np.pi * var.thermal_diff_ground * t_w), 0)
    return heat_flux_ground[:,:]


# Heat flux from ground using 1D numerical model. Different domains of ground can have different types of substrate. Also tracks surface temperature.
def heat_flux_ground(state, var):
    if var.heat_flux_1d == False:   # Assume partly perfect thermal contact between ground surface and pool
        return heat_flux_ground_partly_perf(state, var)    
    else:                           # Use 1D model for ground heat flux
        if state.t == 0:        
            # Load ground heat flux from 1D model for each domain
            data_domain_1 = np.load(var.heat_flux_data_domain_1)
            data_domain_2 = np.load(var.heat_flux_data_domain_2)                        
            var.heat_flux_ground_domain_1 = data_domain_1['heat_flux']      # Ground heat flux (W/m^2) as a function of time
            var.heat_flux_ground_domain_2 = data_domain_2['heat_flux']      # Ground heat flux (W/m^2) as a function of time

            var.temperature_domain_1 = data_domain_1['temperature']         # Surface temperature (K) in domain 1
            var.temperature_domain_2 = data_domain_2['temperature']         # Surface temperature (K) in domain 2

            var.time_domain_1 = data_domain_1['time']         # Time (s) for each value of the ground heat flux
            var.time_domain_2 = data_domain_2['time']         # Time (s) for each value of the ground heat flux
            
        heat_flux_ground = state.q[0,:,:] * 0
        surface_temperature = state.q[0,:,:] * 0
        t_w = state.problem_data['time_wet']
        delta_x = var.delta_x
        delta_y = var.delta_y

        # Get time index in 1D model for each ground grid cell and domain
        time_indices_domain_1 = np.searchsorted(var.time_domain_1, t_w[:,:])           # Depending on how long a cell has been wet, find corresponding time (index) for heat flux in 1D model
        time_indices_domain_2 = np.searchsorted(var.time_domain_2, t_w[:,:])           # Depending on how long a cell has been wet, find corresponding time (index) for heat flux in 1D model

        # Extract heat flux from 1D model for each domain
        heat_flux_ground_1_domain_1 = var.heat_flux_ground_domain_1[time_indices_domain_1[:,:]]         # Extract from 1D model the heat flux depending on how long a cell has been wet (lower value of heat flux)
        heat_flux_ground_2_domain_1 = var.heat_flux_ground_domain_1[time_indices_domain_1[:,:] - 1]     # Extract from 1D model the heat flux depending on how long a cell has been wet (upper value of heat flux)
        heat_flux_ground_1_domain_2 = var.heat_flux_ground_domain_2[time_indices_domain_2[:,:]]         # Extract from 1D model the heat flux depending on how long a cell has been wet (lower value of heat flux)
        heat_flux_ground_2_domain_2 = var.heat_flux_ground_domain_2[time_indices_domain_2[:,:] - 1]     # Extract from 1D model the heat flux depending on how long a cell has been wet (upper value of heat flux)
        
        heat_flux_ground_1_domain_1 = np.where(t_w[:,:] == 0, 0, heat_flux_ground_1_domain_1)           # If ground grid cells are dry, heat flux is zero   
        heat_flux_ground_2_domain_1 = np.where(t_w[:,:] == 0, 0, heat_flux_ground_2_domain_1)
        heat_flux_ground_1_domain_2 = np.where(t_w[:,:] == 0, 0, heat_flux_ground_1_domain_2)
        heat_flux_ground_2_domain_2 = np.where(t_w[:,:] == 0, 0, heat_flux_ground_2_domain_2)
        
        heat_flux_ground_domain_1 = (heat_flux_ground_1_domain_1 + heat_flux_ground_2_domain_1) / 2     # Average of upper and lower value of heat flux in domain 1
        heat_flux_ground_domain_2 = (heat_flux_ground_1_domain_2 + heat_flux_ground_2_domain_2) / 2     # Average of upper and lower value of heat flux in domain 2
        
        N_i = len(heat_flux_ground[:]) - int(len(heat_flux_ground[:]) / 2 - 2 / delta_x) + 2            # Index for the boundary between domain 1 and 2

        heat_flux_ground[0:N_i, 0:] = heat_flux_ground_domain_1[0:N_i, 0:]
        heat_flux_ground[N_i:, :] = heat_flux_ground_domain_2[N_i:, :]

        # Extract surface temperature from 1D model for each domain
        surface_temperature_1_domain_1 = var.temperature_domain_1[time_indices_domain_1[:,:]]           # Extract from 1D model the surface temperature depending on how long a cell has been wet (lower value of temperature)
        surface_temperature_2_domain_1 = var.temperature_domain_1[time_indices_domain_1[:,:] - 1]       # Extract from 1D model the surface temperature depending on how long a cell has been wet (upper value of temperature)
        surface_temperature_1_domain_2 = var.temperature_domain_2[time_indices_domain_2[:,:]]           # Extract from 1D model the surface temperature depending on how long a cell has been wet (lower value of temperature)
        surface_temperature_2_domain_2 = var.temperature_domain_2[time_indices_domain_2[:,:] - 1]       # Extract from 1D model the surface temperature depending on how long a cell has been wet (upper value of temperature)

        surface_temperature_1_domain_1 = np.where(t_w[:,:] == 0, var.T_ground_inf, surface_temperature_1_domain_1)     # If ground grid cells are dry, surface temperature is equal to initial ground temperature
        surface_temperature_2_domain_1 = np.where(t_w[:,:] == 0, var.T_ground_inf, surface_temperature_2_domain_1)
        surface_temperature_1_domain_2 = np.where(t_w[:,:] == 0, var.T_ground_inf, surface_temperature_1_domain_2)
        surface_temperature_2_domain_2 = np.where(t_w[:,:] == 0, var.T_ground_inf, surface_temperature_2_domain_2)

        surface_temperature_domain_1 = (surface_temperature_1_domain_1 + surface_temperature_1_domain_1) / 2   # Average of upper and lower value of surface temperature in domain 1
        surface_temperature_domain_2 = (surface_temperature_1_domain_2 + surface_temperature_2_domain_2) / 2   # Average of upper and lower value of surface temperature in domain 2

        surface_temperature[0:N_i, 0:] = surface_temperature_domain_1[0:N_i, 0:]
        surface_temperature[N_i:, :] = surface_temperature_domain_2[N_i:, :]

        state.problem_data['surface_temperature'] = surface_temperature

        return heat_flux_ground[:,:]


# Heat flux from air (wind) assuming circular pool, "Methods for the Calculation of Physical Effects" (Bosch & Weterings, 2005)
def heat_flux_air(state, r_max, var):    # r_max is the maximal pool radius (m)
    heat_flux_air = state.q[0,:,:] * 0
    if r_max >= 0.1:
        Pr_air = 10 ** 9 / (1.1 * var.T_air ** 3 - 1200 * var.T_air ** 2 + 322000 * var.T_air + 1.393 * 10 ** 9)    # Prandtl's number for air
        Re = var.rho_air * var.wind_speed * 2 * r_max / var.dyn_viscosity_air                                       # Reynolds' number for turbulent flow
        Nu = 0.037 * Pr_air ** (1 / 3) * Re ** 0.8                                                      # Nusselt's number
        heat_flux_air_value = Nu * var.thermal_cond_air / (2 * r_max) * (var.T_air - var.T_pool)
    else:
        heat_flux_air_value = 0
    heat_flux_air += heat_flux_air_value
    return heat_flux_air[:,:]