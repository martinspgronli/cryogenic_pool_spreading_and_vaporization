"""
This Python script uses PyClaw (part of Clawpack) to solve the 2D shallow water equations 
for a spill on a solid surface with a an arbitrary size and shape. It allows for both 
instantaneous and continuous time dependent spills with and without evaporation. 
Moreover, it is possible to place obstacles of arbitrary shapes which the fluid can flow 
both over and around. Flowing over is computationally more expensive. The code is used to 
compute the pool radius (m) and boil-off mass (kg/m^2s) from a pool of liquid ammonia (LNH3)
and liquid hydrogen (LH2).

Time-dependent heat flux from the ground and air into the pool is considered. All other
heat fluxes are neglected. The heat flux from the ground can be loaded from a 1D model 
that solves the 1D heat equation in the ground. To calculate the heat flux at a given 
temperature, run heat_flux/boiling_1d.py adjusted to the correct temperature. This will 
produce an npz-file that can be imported in the variable file for the chosen fluid.

Pool spreading is limited (due to friction) by a Manning friction source term. 
This can be turned off. The term is applied where the depth of the pool is below 
a certain depth_tolerance. This depth_tolerance and the Manning coefficient, 
manning_coeff, should be chosen carefully when simulating real-world hazards. 

The maximal size of the pool (called radius r_p) is calculated both with and without 
a "cut-off". A cut-off means that the ground is regarded as dry when the depth is below 
some small height. Additionally, the calculation of radius assumes the spill to be centered 
in (x, y) = (0, 0).

For a large t_final, the pool might become unphysically thin. Friction only reduces the
speed at which the pool spreads, but will never stop it completely. A minimal depth 
could therefore be set when developing the code further.

All units follow the SI system.

BEFORE RUNNING:
                -Set fluid variable to "H2" or "NH3" in main.py. 
                -All variables are set in the variable file for the chosen fluid


Author: Martin Spillum Grønli (2024, SINTEF Energy Research)
Code is based on https://github.com/clawpack/apps/blob/master/notebooks/pyclaw/beach.ipynb


####################################################################################################################################
####################################################################################################################################


TO-DO:
    -This script does not account for the reheating of ground where the pool has dried. 
     If dried ground is rewetted, the pool will experience the heat flux as if it is was always wet. 
     Implement solution where pool stops to cool when dry.
    -Maximal distance from center of spill area assumes spill to be centered in (x, y) = (0, 0). 
     Implement a more general solution.
"""


# Import 
import importlib
from pool_spreading import *
from plot_functions import *


# Import variables from fluid specific files (variables_NH3.py or variables_H2.py)
def import_variables(fluid):
    file_name_NH3 = 'variables_NH3'
    file_name_H2 = 'variables_H2'
    if fluid == "NH3":
        file_name = file_name_NH3
    elif fluid == "H2":
        file_name = file_name_H2
    else:
        raise ValueError("No data for this fluid! Fluid =", fluid)
    variables = importlib.import_module(file_name)
    variables.fluid = fluid
    return variables


# Main function
def main(var):
    # Run pool spreading and vaporization
    claw = setup(var)

    claw.run()
    check_mass_conservation(claw, var)
   
    # Plot results
    plot_radius_pool_cut_off(claw)
    plot_evaporated_mass(claw)
    plot_height_contour(claw)
    plot_evaporated_mass_per_area_contour(claw)
    
    plot_evaporated_mass_per_area_cross_section(claw)
    plot_heigt_cross_section(claw, var)    
    plot_heigt_cross_section(claw, var, zoomed = True)    
    plot_speed_cross_section(claw, var)
    plot_height_3D(claw, var)

    return claw


if __name__ == '__main__':
    fluid = "NH3"       # Set which fluid to be spilled ("H2" or "NH3")
    variables = import_variables(fluid)
    claw = main(variables)