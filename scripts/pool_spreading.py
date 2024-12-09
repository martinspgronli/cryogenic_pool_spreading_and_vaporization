# Import
import numpy as np
import pandas as pd
import scipy.integrate as integrate
import warnings
import csv
from clawpack import riemann
from clawpack import pyclaw
from heat_flux import *


# Initialize state with circular or rectangular spill
def qinit(state, var):
    x0 = 0.
    y0 = 0.
    X, Y = state.p_centers
    r = np.sqrt((X - x0) ** 2 + (Y - y0) ** 2)  # Distance from center 
    if var.V_0 > 0:
        if var.spill_height == 0:       
            spill_height = 1        # In case spill_height is set to 0 even though V_0 > 0
        if var.shape_spill == "circular":       
            r = np.sqrt((X - x0) ** 2 + (Y - y0) ** 2)  # Distance from center
            state.q[0,:,:] = spill_height * (r <= var.initial_radius_spill) + var.min_height * (r > var.initial_radius_spill)
            spill_height *= var.V_0 / (np.sum(state.q[0,:,:] * var.delta_x * var.delta_y))
            state.q[0,:,:] = spill_height * (r <= var.initial_radius_spill) + var.min_height * (r > var.initial_radius_spill)   # Adjust spill_height such that initial volume exactly equals V_0
        elif var.shape_spill == "rectangular":
            domain_spill = (X <= var.initial_length_spill / 2) * (Y <= var.initial_width_spill / 2) \
                          * (X <= var.initial_length_spill / 2) * (Y >= - var.initial_width_spill / 2) \
                          * (X >= - var.initial_length_spill / 2) * (Y <= var.initial_width_spill / 2) \
                          * (X >= - var.initial_length_spill / 2) * (Y >= - var.initial_width_spill / 2)    # Make rectangular domain
            state.q[0,:,:] = spill_height * domain_spill
            spill_height *= var.V_0 / (np.sum(state.q[0,:,:] * var.delta_x * var.delta_y))
            state.q[0,:,:] = spill_height * domain_spill + var.min_height * np.invert(domain_spill)     # Adjust spill_height such that initial volume exactly equals V_0  
    else:
        state.q[0,:,:] = 0
    state.q[1,:,:] = 0.     # Set initial velocity in x-direction to zero if no initial spill
    state.q[2,:,:] = 0.     # Set initial velocity in y-direction to zero if no initial spill


# Include obstacles by defining topography (called bathymetry in Clawpack)
def topography(x, y, obstacle_categories):
    if "None" in obstacle_categories or len(obstacle_categories) == 0:  # No obstacles
        return 0
    else:
        dikes = x * 0
        box_1 = x * 0
        box_2 = x * 0
        if "dikes" in obstacle_categories:      # A rectangular area surrounded by dikes
            height_dikes = 0.04                 # Height of dikes
            width_dikes = 1                     # Width of dikes
            x_up = - 2 - width_dikes            # Start of dike 1
            x_low = - x_up                      # Start of dike 2
            y_up = - 2 - width_dikes            # Start of dike 3
            y_low = - y_up                      # Start of dike 4
            dikes =  (x > x_up) * (x < x_up + width_dikes) #+ (x < x_low) * (x > x_low - width_dikes) \
                   #+ (y > y_up) * (y < y_up + width_dikes) + (y < y_low) * (y > y_low - width_dikes)   # Make rectangular domain
            dikes = dikes * height_dikes
          
        if "box_1" in obstacle_categories:  # A rectangular box
            height_box_1 = 0.04     # Height of box
            length_box_1 = 0.5      # Length of box (x-dir)
            width_box_1  = 0.5      # Width of box (y-dir)
            x_center = 2            # Center of box (x)
            y_center = - 1.5        # Center of box (y)
            box_1 =   (x > x_center - length_box_1 / 2) * (x < x_center + length_box_1 / 2) \
                    * (y > y_center - width_box_1 / 2) * (y < y_center + width_box_1 / 2)
            box_1 = box_1 * height_box_1

        if "box_2" in obstacle_categories:  # A rectangular box
            height_box_2 = 0.04     # Height of box
            length_box_2 = 0.5      # Length of box (x-dir)
            width_box_2  = 0.5      # Width of box (y-dir)
            x_center = 2            # Center of box (x)
            y_center = 1.5          # Center of box (y)
            box_2 =   (x > x_center - length_box_2 / 2) * (x < x_center + length_box_2 / 2) \
                    * (y > y_center - width_box_2 / 2) * (y < y_center + width_box_2 / 2)
            box_2 = box_2 * height_box_2
            
    return dikes + box_1 + box_2


# Source term for spill rate
def source_spill(state, dt, var):
    src_spill = np.zeros(state.q.shape)       # Source term for spill
    if var.m_dot_spill(state.t) > 0 and state.t < var.t_spill_stop:    
        x0 = 0.
        y0 = 0.
        X, Y = state.p_centers
        r = np.sqrt((X - x0) ** 2 + (Y - y0) ** 2)
        delta_x = var.delta_x
        delta_y = var.delta_y

        if var.shape_spill == "circular":       # If circular spill area
            h_dot_spill = var.m_dot_spill(state.t) / (var.rho_l * np.pi * var.radius_spill ** 2)
            src_temp = np.zeros(state.q[0].shape)
            src_temp[:,:] = np.where(r[:,:] <= var.radius_spill, h_dot_spill, 0)
            h_dot_spill *= var.m_dot_spill(state.t) / (var.rho_l * delta_x * delta_y * np.sum(src_temp))   # Adjust spill height such that exactly m_dot_spill is spilled.
            src_spill[0,:,:] += np.where(r[:,:] <= var.radius_spill, dt * h_dot_spill, 0)
            if var.initial_speed_x != 0 or var.initial_speed_y != 0:
                mass = var.rho_l * delta_x * delta_y * state.q[0,:,:]
                print('Total mass inside spill radius (kg):', np.sum(np.where(r[:,:] <= var.radius_spill, mass, 0)))
                state.q[1,:,:] = np.where(r[:,:] <= var.radius_spill, (var.initial_speed_x * dt * var.m_dot_spill(state.t) + mass[:,:] * state.q[1,:,:]) / (dt * var.m_dot_spill(state.t) + mass[:,:]), state.q[1,:,:])     # Velocity equals a weighted average of velocity of spill already present inside spill domain and velocity of new spill
                state.q[2,:,:] = np.where(r[:,:] <= var.radius_spill, (var.initial_speed_y * dt * var.m_dot_spill(state.t) + mass[:,:] * state.q[2,:,:]) / (dt * var.m_dot_spill(state.t) + mass[:,:]), state.q[2,:,:])   
            state.problem_data['mass_spilled_total'] += np.sum(np.where(r[:,:] <= var.radius_spill, dt * var.rho_l * delta_x * delta_y * h_dot_spill, 0))     # For some reason overestimated but should be OK

        elif var.shape_spill == "rectangular":  # If rectangular spill area
            h_dot_spill = var.m_dot_spill(state.t) / (var.rho_l * var.length_spill * var.width_spill)
            domain_spill =  (X <= var.length_spill / 2) * (Y <= var. width_spill / 2) \
                            * (X <= var.length_spill / 2) * (Y >= - var.width_spill / 2) \
                            * (X >= - var.length_spill / 2) * (Y <= var.width_spill / 2) \
                            * (X >= - var.length_spill / 2) * (Y >= - var.width_spill / 2)    # Make rectangular spill domain
            src_temp = np.zeros(state.q[0].shape)
            src_temp[:,:] = np.where(domain_spill, h_dot_spill, 0)
            h_dot_spill *= var.m_dot_spill(state.t) / (var.rho_l * delta_x * delta_y * np.sum(src_temp))
            src_spill[0,:,:] += dt * h_dot_spill * domain_spill
            if var.initial_speed_x != 0 or var.initial_speed_y != 0:    
                mass = var.rho_l * delta_x * delta_y * state.q[0,:,:]
                state.q[1,:,:] = np.where(domain_spill, (var.initial_speed_x * dt * var.m_dot_spill(state.t) + mass[:,:] * state.q[1,:,:]) / (dt * var.m_dot_spill(state.t) + mass[:,:]), state.q[1,:,:])   # Velocity equals a weighted average of velocity of spill already present inside spill domain and velocity of new spill
                state.q[2,:,:] = np.where(domain_spill, (var.initial_speed_y * dt * var.m_dot_spill(state.t) + mass[:,:] * state.q[2,:,:]) / (dt * var.m_dot_spill(state.t) + mass[:,:]), state.q[2,:,:])
            state.problem_data['mass_spilled_total'] += np.sum(np.where(domain_spill, dt * var.rho_l * delta_x * delta_y * h_dot_spill, 0))     # For some reason overestimated but should be OK

    return src_spill[:,:,:]


# Calculate maximal distance from center of spill to edge of pool. Maximal distance assumes center of spill in (x, y) = (0, 0) and equals radius if circular spill
def calculate_radius(state, var):
    x0 = 0.
    y0 = 0.
    X, Y = state.p_centers
    r = np.sqrt((X - x0) ** 2 + (Y - y0) ** 2)
    
    # Calculate radius of pool without a cut-off threshold 
    dist_to_dry_cells = np.ma.compressed(np.ma.masked_where(state.q[0,:,:] == var.min_height, r[:,:]))  
    try:
        radius_max = np.max(dist_to_dry_cells)
    except:     # If min_height > 0
        radius_max = 0
    state.problem_data['radius'] = radius_max
    
    # Calculate radius of pool with a cut-off threshold according to dry_tol
    dist_to_dry_cells_cut_off = np.ma.compressed(np.ma.masked_where((state.q[0,:,:] <= var.min_height + var.dry_tol), r[:,:]))   
    try:
        radius_max_cut_off = np.max(dist_to_dry_cells_cut_off)
    except:     # If min_height > 0
        radius_max_cut_off = 0
    state.problem_data['radius_cut_off'] = radius_max_cut_off
    print('Max size =', radius_max_cut_off, 'm')
    return radius_max, radius_max_cut_off


# Manning friction term, see https://doi.org/10.1016/j.advwatres.2009.02.010
def friction_manning(state, dt, var):    
    if var.friction and state.t > 0:             
        state.q[1,:,:] = np.where(state.q[0,:,:] < var.depth_tolerance, 0, state.q[1,:,:])    # If depth is below depth_tolerance, set momentum to zero
        state.q[2,:,:] = np.where(state.q[0,:,:] < var.depth_tolerance, 0, state.q[2,:,:])
                    
        gamma = np.zeros(state.q[0,:].shape)     
        gamma =  np.where((state.q[0,:,:] <= var.friction_depth) & (state.q[0,:,:] > var.depth_tolerance) & (state.q[1,:,:] ** 2 + state.q[2,:,:] ** 2 > 0), state.problem_data['grav'] * var.manning_coeff ** 2 * np.sqrt(state.q[1,:,:] ** 2 + state.q[2,:,:] ** 2) / (state.q[0,:,:] ** (7 / 3)), 0)

        dgamma = 1 + dt * gamma
        state.q[1,:,:] = state.q[1,:,:] / dgamma
        state.q[2,:,:] = state.q[2,:,:] / dgamma    


# Calculate time grid cells have been wet. Drying and then rewetting of ground is not accounted for.
def time_wet_calculation(state, r, radius_max, var):
    state.problem_data['initial_time_wet'] = np.where((r[:,:] <= radius_max) & (state.q[0,:,:] > var.depth_tolerance) & (state.problem_data['initial_time_wet'] == - 1), state.t, state.problem_data['initial_time_wet'])        # If cell has not been wet before, set initial_time_wet to state.t
    state.problem_data['initial_time_wet'] = np.where((r[:,:] <= radius_max) & (state.q[0,:,:] > var.depth_tolerance) & (state.t < state.problem_data['initial_time_wet']), state.t, state.problem_data['initial_time_wet'])     # If state.t less than previous state.t and cell is registered as wet previously (PyClaw did not accept previous time step), set initial_time_wet to state.t
    state.problem_data['time_wet'] = np.where((r[:,:] <= radius_max) & (state.q[0,:,:] > var.depth_tolerance) & (state.problem_data['initial_time_wet'] != - 1) & (state.t < state.problem_data['initial_time_wet']), 0, state.problem_data['time_wet'])        # If state.t less than previous state.t and cell is not registered as wet previosly, set initial_time_wet to 0
    state.problem_data['time_wet'] = np.where((r[:,:] <= radius_max) & (state.q[0,:,:] > var.depth_tolerance) & (state.problem_data['initial_time_wet'] != - 1) & (state.t > state.problem_data['initial_time_wet']), state.t - state.problem_data['initial_time_wet'], state.problem_data['time_wet'])  # If wet previously and state.t greater than initial_time_wet, update time it has been wet


# Vaporization of pool
def source_evaporate(state, dt, radius_max, var):
    src_evap = np.zeros(state.q.shape)      # Source term for vaporization

    x0 = 0.
    y0 = 0.
    X, Y = state.p_centers
    r = np.sqrt((X - x0) ** 2 + (Y - y0) ** 2)
    
    if radius_max > 0:                      # If pool is present
        delta_x = var.delta_x
        delta_y = var.delta_y
        h_dot_evap = np.zeros(state.q[0].shape)             # Height to be subtracted in each cell due to evaporation

        # Calculate time a cell has been wet
        time_wet_calculation(state, r, radius_max, var)     

        # Calculate heat flux (W/m^2) from air and ground
        heat_flux = heat_flux_air(state, radius_max, var) + heat_flux_ground(state, var)

        h_dot_evap =  heat_flux / (var.vap_enthalpy * var.rho_l)   # Height to be subtracted in each cell due to evaporation
        h_dot_evap[:,:] = np.where((state.problem_data['initial_time_wet'][:,:] == - 1) | (state.problem_data['time_wet'][:,:] == 0), 0, h_dot_evap[:,:])   

        # Add evaporation to source term for evaporation
        src_evap[0,:,:] += np.where((r[:,:] <= radius_max) & (state.q[0,:,:] > var.depth_tolerance) & (state.q[0,:,:] >= h_dot_evap * dt), - dt * h_dot_evap, 0)
        src_evap[0,:,:] += np.where((r[:,:] <= radius_max) & (state.q[0,:,:] > var.depth_tolerance) & (state.q[0,:,:] < h_dot_evap * dt), - dt * state.q[0,:,:], 0)    # Set source term to negative of height * dt if h < h_dot_evap * dt. I.e. make cell dry

        # Log evaporated mass
        mass_boil_inst = np.zeros(state.q[0].shape)
        mass_boil_inst[:,:] += np.where((r[:,:] <= radius_max) & (state.q[0,:,:] > var.depth_tolerance) & (state.q[0,:,:] >= h_dot_evap * dt), var.rho_l * delta_x * delta_y * h_dot_evap, 0)
        mass_boil_inst[:,:] += np.where((r[:,:] <= radius_max) & (state.q[0,:,:] > var.depth_tolerance) & (state.q[0,:,:] < h_dot_evap * dt), var.rho_l * delta_x * delta_y * state.q[0,:,:], 0)
        state.problem_data['mass_evap_inst_per_area'] = mass_boil_inst[:,:] / (delta_x * delta_y)
        state.problem_data['mass_evap_inst'] = np.sum(mass_boil_inst) 
        state.problem_data['mass_evap_total'] += state.problem_data['mass_evap_inst'] * dt
        print('m_dot_vap (kg/s) =', state.problem_data['mass_evap_inst'])

    return src_evap[:,:,:]


# Source terms for spill, evaporation and friction. source_outer passes variables to source
def source_outer(var):
    def source(solver, state, dt):
        print('t =', state.t, 's')
        
        state.problem_data['time'] = state.t    # Log current time
       
        src = np.zeros(state.q.shape)           # Source term            

        # Spill into pool by increasing height inside spill area
        src[:,:,:] += source_spill(state, dt, var)

        # Calculate radius of pool (without and with cut-off threshold)
        radius_max, radius_max_cut_off = calculate_radius(state, var)

        # Apply vaporization of pool by decreasing height where radius < radius_max (no cut-off) and pool height > depth_tolerance
        src[:,:,:] += source_evaporate(state, dt, radius_max, var)
        
        # Apply Manning friction source term by altering momentum, see https://doi.org/10.1016/j.advwatres.2009.02.010
        friction_manning(state, dt, var) 

        return src
    return source


# Solver and controller setup
def setup(var):
    # Define solver
    solver = pyclaw.SharpClawSolver2D(riemann.sw_aug_2D)
    solver.cfl_max     = 0.9        # Maximal allowed CFL number for a time step to be accepted
    solver.cfl_desired = 0.75       # Desired CFL number
    source = source_outer(var)      # Pass variables to uncompiled source function
    solver.dq_src = source          # Include source terms
    solver.lim_type = 1             # 1 corresponds to tvd reconstruction, while 2 corresponds to WENO reconstruction
    #solver.weno_order = 5          # Order of the WENO reconstruction, Default = 5
    #solver.limiters = pyclaw.limiters.tvd.minmod   # Limiter to be used, Defualt = limiters.tvd.minmod
    solver.dt_initial = var.tfinal / var.num_output_times   # For some reason the initial time step is too large. Set smaller
    #solver.dt_variable = True      # Allow for variable time step

    solver.bc_lower[0] = pyclaw.BC.extrap           # Non-reflection outflow over the boundaries. WARNING: Setting pyclaw.BC.wall does not conserve mass
    solver.bc_upper[0] = pyclaw.BC.extrap 
    solver.bc_lower[1] = pyclaw.BC.extrap  
    solver.bc_upper[1] = pyclaw.BC.extrap

    solver.aux_bc_lower[0] = pyclaw.BC.wall         # Set wall boundary conditions for the topography
    solver.aux_bc_upper[0] = pyclaw.BC.wall
    solver.aux_bc_lower[1] = pyclaw.BC.wall
    solver.aux_bc_upper[1] = pyclaw.BC.wall

    solver.fwave = True     # Riemann solver returns f-waves

    # Define domain
    x_lower = - var.x_dimension;  x_upper = var.x_dimension
    y_lower = - var.y_dimension;  y_upper = var.y_dimension
    mx = var.num_cells
    my = int(np.ceil((var.num_cells * var.y_dimension / var.x_dimension) / 2) * 2)   # Ensure equal grid spacing in x- and y-direction
    x = pyclaw.Dimension(x_lower, x_upper, mx, name = 'x')
    y = pyclaw.Dimension(y_lower, y_upper, my, name = 'y')  
    var.delta_x = x.delta                           # Grid spacing (m) in x-direction
    var.delta_y = y.delta                           # Grid spacing (m) in y-direction
    domain = pyclaw.Domain([x, y])

    # Define state and topography
    num_aux = 1         # Number of auxiliary fields
    state = pyclaw.State(domain, solver.num_eqn, num_aux)   # Define state
    state.aux[:,:,:] = topography(state.p_centers[0], state.p_centers[1], var.obstacles)    # Define topography
    qinit(state, var)   # Initialize state

    # Set data for state. state.problem_data is used to store data at each time step but should be avoided since it fills up memory
    state.problem_data['grav'] = 9.81               # Gravitational acceleration (m/s^2)
    state.problem_data['time'] = 0                  # Track time (s)             

    # Variables used to track each grid cell at each time step. mx times my matrices 
    state.problem_data['mass_evap_total'] = 0       # Total evaporated mass (kg)    
    state.problem_data['mass_spilled_total'] = 0    # Total mass (kg) spilled into pool. For some reason overestimated but should be OK
    state.problem_data['radius'] = var.initial_radius_spill             # Radius (m) without cut-off
    state.problem_data['radius_cut_off'] = var.initial_radius_spill     # Radius (m) with cut-off
    state.problem_data['mass_evap_inst'] = 0                            # Mass (kg/s) evaporated per time
    state.problem_data['mass_evap_inst_per_area'] = np.zeros(state.q[0,:,:].shape)  # Mass (kg/sm^2) evaporated per area per time
    
    # Variables used to track each grid cell at each time step. mx times my matrices
    state.problem_data['time_wet'] = np.zeros(state.q[0,:,:].shape)                 # Time (s) grid cells have been wet 
    state.problem_data['initial_time_wet'] = np.zeros(state.q[0,:,:].shape) - 1     # Time (s) at which grid cells got wet (- 1 if never wet) 
    state.problem_data['surface_temperature'] = np.zeros(state.q[0,:,:].shape)      # Surface temperature (K)
        
    # Set up controller and controller parameters
    claw = pyclaw.Controller()
    claw.tfinal = var.tfinal
    claw.t_spill_stop = var.t_spill_stop
    claw.fluid = var.fluid
    claw.solution = pyclaw.Solution(state, domain)
    #claw.solution = pyclaw.Solution(10, file_format='ascii')       # Restart simulation at a given time (frame number). Does not work because of a bug in PyClaw, but is now fixed if PyClaw is updated
    claw.solver = solver
    claw.output_format = 'ascii'
    claw.output_aux_onlyonce = True
    claw.num_output_times = var.num_output_times
    claw.keep_copy = True
    claw.output_format = None      # Turn off ouput to _output folder
    claw.spill_height = var.spill_height
    claw.mx = mx
    claw.my = my
    return claw


# Control result by checking mass conservation
def check_mass_conservation(claw, var):
    mass_initial = np.sum(claw.frames[0].q[0,:,:] * var.delta_x * var.delta_y) * var.rho_l
    mass_spilled = np.sum(claw.frames[- 1].problem_data['mass_spilled_total'])
    mass_evaporated = np.sum(claw.frames[- 1].problem_data['mass_evap_total'])
    mass_left = np.sum(claw.frames[- 1].q[0,:,:] * var.delta_x * var.delta_y) * var.rho_l
    print('\nFluid is set to', claw.fluid, '\n')
    print('Initial mass of spill (kg) =', mass_initial)
    print('Total spilled mass (kg) =', mass_spilled)
    print('Total evaporated mass (kg) =', mass_evaporated)
    print('Total mass left (kg) =', mass_left)
    print('(Mass initial + Mass spilled - Mass evaporated - Mass left) / (Mass initial + Mass spilled) =', (mass_initial + mass_spilled - mass_evaporated - mass_left) / (mass_initial + mass_spilled))
    print('Code is set to spill', integrate.quad(var.m_dot_spill, 0, min(claw.t_spill_stop, claw.tfinal))[0], 'kg after t = 0', '\n')
    if (mass_initial + mass_spilled - mass_evaporated - mass_left) / (mass_initial + mass_spilled) > 0.05:
        warnings.warn("Mass conservation may not be satisfied")
    if mass_spilled > 0:
        print('Total mass spilled is for some reason overestimated. Instead approximate mass spilled \u2248 mass left + mass evaporated - mass initial =', mass_left + mass_evaporated - mass_initial, 'kg')       