"""
This Python file contains variables used in the pool spreading of LNH3.


Author: Martin Spillum Grønli (2024, SINTEF Energy Research)
"""
 

# Properties of LNH3 at 1 atm and 239.832 K
rho_l                   = 682.64    # LNH3 density at 1 atm (kg/m^3). NIST Webbok
rho_v                   = 0.88997   # Density of NH3 vapor at 1 atm and boiling point (kg/m^3). NIST Webbok
vap_enthalpy            = 1.38e6    # Enthalpy of vaporization of NH_3 (J/kg). NIST Webbok

# Properties for calculation of heat flux from ground in case heat_flux_1d=False and perfect thermal contact is assumed as in Eq. (16) in doi.org/10.1016/j.ijhydene.2020.06.131 
T_pool                  = 239.832       # LNH3 temperature (K). Atmospheric boiling point of ammonia. NIST Webbook
T_ground_inf            = 273.15 + 10   # Ground temperature at infinite depth
thermal_cond_ground     = 3.72          # Thermal conductivity (W/mK) of the ground
thermal_diff_ground     = 8.227996e-7   # Thermal diffusivity (m^2/s) of the ground 

# Properties for calculation of heat flux from air
T_air                   = 273.15 + 10   # Temperature (K) in the air 
rho_air                 = 1.225     # Density (kg/m^3) of air at T = 15 ℃ = 288 K
wind_speed              = 2         # Wind speed (m/s) 10 m above the pool 
dyn_viscosity_air       = 1.802e-5  # Dynamic viscosity (kg/ms) of air at T = 15 ℃ = 288 K
thermal_cond_air        = 0.02476   # Thermal conductivity (W/mK) of air at T = 15 ℃ = 288 K, https://www.engineersedge.com/physics/viscosity_of_air_dynamic_and_kinematic_14483.htm 


# Pool variables
V_0                     = 0.001     # Initial volume (m^3), spill_height is adjusted thereafter
spill_height            = 0         # Initial height (m) of spill. Later adjusted to fit with V_0 and initial dimensions of the spill to avoid V_0 being dependent on grid spacing
shape_spill             = "rectangular" # Geometric shape of spill, "circular" or "rectangular"
initial_radius_spill    = 0.75      # Initial spill radius (m) if shape_spill == "circular"
initial_length_spill    = 0.75      # Initial spill length (m) in x-direction if shape_spill == "rectangular"
initial_width_spill     = 1         # Initial spill width (m) in y-direction if shape_spill == "rectangular"
initial_speed_x         = 1         # Initial speed in x-direction for continuous spill (m/s)
initial_speed_y         = 0         # Initial speed in y-direction for continuous spill (m/s)
min_height              = 0         # Initial height (m) outside initial spill area   
m_dot_spill             = lambda t: 15 - 0.25 * t    # Continuous spill rate (kg/s) inside radius spill
radius_spill            = 0.75      # The continuous spill is spilled within this radius (m) if shape_spill == "circular"
length_spill            = 1         # The continuous spill is spilled within this length (m) in x-direction if shape_spill == "rectangular"
width_spill             = 0.75      # The continuous spill is spilled within this width (m) in y-direction if shape_spill == "rectangular"
t_spill_stop            = 38        # Time (s) when spill is stopped. Set to 0 if no spill

obstacles               = ["dikes", "box_1", "box_2"]  # Include obstacles, "dikes", "box_1", "box_2", "None". 
friction                = True      # Turn on Manning friction term
manning_coeff           = 0.015     # Manning roughness of concrete, https://www.engineeringtoolbox.com/mannings-roughness-d_799.html
friction_depth          = 99999     # Manning friction term is applied when pool depth is below this height (m)
dry_tol                 = 1e-5      # Height (m) where pool is regarded as dry when calculating radius with a cut-off

depth_tolerance         = 0         # Momentum is set to zero below this height (m). Spill with height < depth_tolerance will not be evaporated or counted into the pool radius
x_dimension             = 5         # 0.5 * length (m) of pool in x-direction 
y_dimension             = 5         # 0.5 * length (m) of pool in y-direction 
num_cells               = 150       # Number of grid cells in x-direction. y-direction will have the ~same grid spacing
tfinal                  = 10        # End time (s). Must be less than some number when heat_flux = True, since heat flux is not calculated longer.

num_output_times        = 100       # Number of times result is outputted (to plots and _output-folder). If too high, RAM can be used up


# Ground heat flux
heat_flux_1d            = True      # Turn on heat flux from 1D model (heat flux calculated until t_end). If not, heat flux is calculated assuming perfect contact between pool and ground, see Eq. (16) in doi.org/10.1016/j.ijhydene.2020.06.131 
NH3_wet_sand_variable_boiling_corr_freeze = '../heat_flux_data/NH3_wet_sand_variable_boiling_corr_freeze_t_end_1400.npz'     # Wet saturated sand with variable thermal properties and boiling correlations
NH3_dry_sand_variable_boiling_corr = '../heat_flux_data/NH3_dry_sand_variable_boiling_corr_t_end_2800.npz'                   # Dry sand with variable thermal properties and boiling correlations
NH3_dry_sand_variable_perf_contact = '../heat_flux_data/NH3_dry_sand_variable_perf_contact_t_end_1400.npz'                   # Dry sand with variable thermal properties and perfect thermal contact
NH3_dry_sand_const_boiling_corr = '../heat_flux_data/NH3_dry_sand_const_boiling_corr_t_end_2800.npz'                         # Dry sand with constant thermal properties and boiling correlations
NH3_dry_sand_const_perf_contact = '../heat_flux_data/NH3_dry_sand_const_perf_contact_t_end_2800.npz'                         # Dry sand with constant thermal properties and perfect thermal contact
NH3_sat_concrete_variable_boiling_corr = '../heat_flux_data/NH3_sat_concrete_variable_boiling_corr_t_end_2000.npz'           # Saturated concrete with variable thermal properties and boiling correlations
NH3_dry_concrete_const_boiling_corr = '../heat_flux_data/NH3_dry_concrete_const_boiling_corr_t_end_2000.npz'                 # Dry concrete with constant (almost no variation in thermal properties) thermal properties and boiling correlations

# Set different heat flux for different domains defined in heat_flux_ground() in heat_flux.py
heat_flux_data_domain_1 = NH3_dry_concrete_const_boiling_corr                    
heat_flux_data_domain_2 = NH3_wet_sand_variable_boiling_corr_freeze