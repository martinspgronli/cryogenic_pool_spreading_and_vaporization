"""
This Python script contains functions for plotting the results of the pool spreading problem. The functions are called from the main file.


Author: Martin Spillum Grønli (2024, SINTEF Energy Research)
"""


# Import
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation, ticker
import matplotlib.font_manager
from matplotlib.animation import FuncAnimation, PillowWriter
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pandas as pd


# Set plot parameters
plt.rcParams.update(plt.rcParamsDefault)
plt.rcParams.update({'font.size': 15})
#plt.rcParams["font.family"] = "Times New Roman"

 
# Plot radius with cut-off
def plot_radius_pool_cut_off(claw):
    #fig = plt.figure(figsize=[12, 7.5])
    fig = plt.figure(2)
    ax1 = fig.add_subplot(111)
    time = np.zeros(len(claw.frames))
    radius_cut_off = np.zeros(len(claw.frames))
    for i in range(len(claw.frames)):
        time[i] = claw.frames[i].problem_data['time']
        radius_cut_off[i] = claw.frames[i].problem_data['radius_cut_off']   
    ax1.plot(time, radius_cut_off, label = 'Radius with cut-off', linewidth = 2)
    ax1.set_xlim(0, claw.tfinal) 
    ax1.set_ylim(0,)
    ax1.set_xlabel(r'$t$ [s]')
    ax1.set_ylabel(r'$r_\mathrm{p}$ [m]')  
    #ax1.legend() 
    plt.show()


# Plot instantaneous and total evaporated mass as a function of time
def plot_evaporated_mass(claw):
    #fig = plt.figure(figsize=[12, 7.5])
    fig = plt.figure(3)
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twinx()
    time = np.zeros(len(claw.frames) - 2)
    mass_evap_inst = np.zeros(len(claw.frames) - 2)     # Instantaneous evaporated mass summed over all grid cells
    mass_evap_total = np.zeros(len(claw.frames) - 2)    # Total (accumulated) evaporated mass
    for i in range(len(mass_evap_inst)):
        time[i] = claw.frames[i + 2].problem_data['time']
        mass_evap_inst[i] = np.sum(claw.frames[i + 2].problem_data['mass_evap_inst'])
        mass_evap_total[i] = np.sum(claw.frames[i + 2].problem_data['mass_evap_total'])
    mass_evap_total = np.insert(mass_evap_total, 0, 0)
    ax1.plot(time, mass_evap_inst, color = 'blue', linewidth = 2)
    time = np.insert(time, 0, 0)
    ax2.plot(time, mass_evap_total, color = 'red', linewidth = 2)
    ax1.tick_params(axis = 'y', labelcolor = 'blue')
    ax2.tick_params(axis = 'y', labelcolor = 'red')
    ax1.set_xlim(0, claw.tfinal)
    ax1.set_ylim(0,)
    ax2.set_ylim(0,)
    ax1.set_xlabel(r'$t$ [s]')
    ax1.set_ylabel(r'$\dot{m}_\mathrm{vap}$ [kg/s]', color = 'blue')  
    ax2.set_ylabel(r'$m_\mathrm{vap,\,tot}$ [kg]', color = 'red')  
    plt.show()


# Plot evaporated mass per area at different times for y = 0
def plot_evaporated_mass_per_area_cross_section(claw):   
    #fig = plt.figure(figsize=[12, 7.5])
    fig = plt.figure(4)
    ax1 = fig.add_subplot(111)
    slice = int(claw.my / 2)
    frame = claw.frames[0]
    x, y = frame.state.grid.p_centers  
    
    for i in range(0, len(claw.frames), int(claw.num_output_times / claw.num_output_times)):
        frame = claw.frames[i]
        ax1.plot(x[:, 0], frame.problem_data['mass_evap_inst_per_area'][:, slice], label = 't = ' + str(np.round(frame.t, 1)) + ' s', linewidth = 2)
    ax1.set_xlabel(r'$x$ [m]')
    ax1.set_ylabel(r'$\dot{m}_\mathrm{vap}/A \mathrm{\,\,[kg/sm^2]}$')  
    #ax1.legend()
    plt.show()


# Contour plot of evaporated mass per area at given time
def plot_evaporated_mass_per_area_contour(claw):
    def change_neighbours_to_false(matrix):
        rows = len(matrix)
        cols = len(matrix[0])

        # Copy the matrix to avoid modifying it while iterating
        new_matrix = np.copy(matrix)
        # Iterate through the matrix
        for i in range(rows):
            for j in range(cols):
                if matrix[i][j] == False:
                    # Set neighboring elements to False
                    for x in range(max(0, i - 1), min(rows, i + 2)):
                        for y in range(max(0, j - 1), min(cols, j + 2)):
                            new_matrix[x][y] = False
        return new_matrix
    
    time = int(len(claw.frames)) - 1     # Plot for this t = time
    plt.figure()
    ax = plt.gca()

    x, y = claw.frames[time].state.grid.p_centers  
    z = claw.frames[time].problem_data['mass_evap_inst_per_area']
    cmap = plt.get_cmap('hot_r')
    cmap.set_under(color = 'white')
    vmin = 1e-40   
    vmax = np.max(z)
    
    n_contour_levels = 200
    min_contour = vmin
    max_contour = 0.16
    levels = np.linspace(min_contour, max_contour, n_contour_levels)

    b = np.copy(claw.frames[time].aux[0,:,:])
    mask_b = np.zeros_like(b, dtype = bool)
    mask_b[b == 0] = True
    mask_b[b != 0] = False

    z = np.ma.array(z, mask = np.invert(mask_b))
    
    im = ax.contourf(x, y, z, cmap = cmap, vmin = vmin, vmax = vmax, levels = levels, corner_mask = False)
    
    mask_b = change_neighbours_to_false(np.copy(mask_b))
    b = np.ma.array(b, mask = mask_b)
    b[b == 0.0] = np.max(claw.frames[time].aux[0,:,:])
        
    ax.contourf(x, y, b, hatches = ['//'], colors = 'grey', alpha = 0.2, corner_mask = False)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size = "5%", pad = 0.05)
    cbar = plt.colorbar(im, cax = cax)
    tick_locator = ticker.MaxNLocator(nbins = 5)

    cbar.locator = tick_locator
    cbar.update_ticks()
    cbar.set_label(r'$\dot{m}_\mathrm{evap} \mathrm{\,\,[kg/sm^2]}$')
    ax.set_xlabel(r'$x$ [m]')
    ax.set_ylabel(r'$y$ [m]')
    ax.set_ylim(- 3, 3)
    ax.set_yticks([-2, 0, 2])
    ax.set_xlim(- 2.5, 5)
    ax.tick_params(axis='x', width = 1, length = 7)
    ax.tick_params(axis='y', width = 1, length = 7)
    ax2 = ax.secondary_xaxis("top")
    ax2.tick_params(top = True, length = 7, labeltop = False, direction = 'in')
    ax.annotate('(b)', xy = (- 1.5, 2.2), size = 18)
    #ax.set_title('t = ' + str(np.round(claw.frames[time].t, 1)) + ' s')
    plt.show()


# Plot cross section of height
def plot_heigt_cross_section(claw, var, ylim = (0, 1.2), zoomed = False):
    #fig = plt.figure(figsize=[12, 7.5])
    fig = plt.figure(5)
    ax1 = fig.add_subplot(111)
    
    fills = []
    frame = claw.frames[0]
    b = frame.aux[0,:,:]
    h = frame.q[0,:,:]

    surface = np.maximum(b, h + b)

    x, y = frame.state.grid.p_centers    
    slice = int(claw.my / 2)

    fill = ax1.fill_between(x[:, 0], b[:, slice], surface[:, slice], facecolor = 'blue')
    fill2 = ax1.fill_between(x[:, 0], 0 * b[:, slice], b[:, slice], facecolor = 'brown')
    fills = [fill, fill2]
    ax1.set_xlim(- var.x_dimension, var.x_dimension)
    ax1.set_ylim(0.0, 0.12) 
    #if ylim: ax1.set_ylim(ylim)

    def fplot(frame_number):
        if zoomed:
            ax1.set_ylim(0.0, var.dry_tol * 1000) 

        fills[- 1].remove()
        fills[- 2].remove()
        frame = claw.frames[frame_number]
        b = frame.aux[0,:,:]
        h = frame.q[0,:,:]
        surface = np.maximum(b, h + b)
     
        #line.set_data(x[:,0],surface[:,slice])
        fill = ax1.fill_between(x[:, 0], b[:, slice], surface[:, slice], facecolor = 'blue', where = b[:, slice] < surface[:, slice])
        fill2 = ax1.fill_between(x[:, 0], 0 * b[:, slice], b[:, slice], facecolor = 'brown')
        fills.append(fill)
        fills.append(fill2)
        timestamp = str(round(frame_number * claw.tfinal / len(claw.frames), 3))
        #ax1.set_title(r'$t = $' + timestamp + ' s')
        ax1.set_xlabel(r'$x$ [m]')
        ax1.set_ylabel(r'$h$ [m]')

        fills.append(fill) 
        return fill,
    
    if zoomed:
        anim = animation.FuncAnimation(fig, fplot, frames = len(claw.frames), interval = 100, repeat = True)
        plt.close()
    else:
        anim = animation.FuncAnimation(fig, fplot, frames = len(claw.frames), interval = 100, repeat = True)
        plt.show()
        plt.close()
    return 0


# Contour plot of height at given time
def plot_height_contour(claw, ylim = (0, 1.2)):
    def change_neighbours_to_false(matrix):
        rows = len(matrix)
        cols = len(matrix[0])

        # Copy the matrix to avoid modifying it while iterating
        new_matrix = np.copy(matrix)
        # Iterate through the matrix
        for i in range(rows):
            for j in range(cols):
                if matrix[i][j] == False:
                    # Set neighboring elements to False
                    for x in range(max(0, i - 1), min(rows, i + 2)):
                        for y in range(max(0, j - 1), min(cols, j + 2)):
                            new_matrix[x][y] = False

        return new_matrix
    
    fig = plt.figure()
    ax1 = fig.add_subplot()

    time = int(len(claw.frames)) - 1     # Plot for this t = time
    frame = claw.frames[time]
    b = frame.aux[0,:,:]

    mask_b = np.zeros_like(b, dtype = bool)
    mask_b[b == 0] = True
    mask_b[b != 0] = False
    
    h = frame.q[0,:,:]
    h = np.ma.array(h, mask = np.invert(mask_b))

    x, y = frame.state.grid.p_centers    
    n_contour_levels = 6
    min_contour = 1e-10
    max_contour = 0.05
    levels = np.linspace(min_contour, max_contour, n_contour_levels)

    vmin = 1e-10
    vmax = 0
    for i in range(len(claw.frames)):
        if np.max(claw.frames[i].q[0, :, :]) > vmax:
            vmax = np.max(claw.frames[i].q[0, :, :])
    
    cmap = plt.get_cmap('viridis')
    cmap.set_under(color = 'white')
    cs = ax1.contourf(x, y, h, cmap = cmap, vmin = vmin, vmax = vmax, levels = levels, corner_mask = False)

    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("right", size = "5%", pad = 0.05)
    
    ax1.contour(cs, colors = 'k')   # Contour line color

    b = np.copy(frame.aux[0,:,:])
    mask_b = np.zeros_like(b, dtype = bool)
    mask_b[b == 0] = True
    mask_b[b != 0] = False
    mask_b = change_neighbours_to_false(mask_b)
    b = np.ma.array(b, mask = mask_b)
    b[b == 0.0] = np.max(frame.aux[0,:,:])
    ax1.contourf(x, y, b, colors = 'grey', alpha = 0.2, hatches = ['//'], corner_mask = False)
    cb = fig.colorbar(cs, cax = cax)

    ax1.set_ylim(- 3, 3)
    ax1.set_yticks([-2, 0, 2])
    ax1.set_xlim(- 2.5, 5)
  
    def fplot3D(frame_number):
        ax1.clear()

        frame = claw.frames[frame_number]
        b = frame.aux[0,:,:]

        mask_b = np.zeros_like(b, dtype = bool)
        mask_b[b == 0] = True
        mask_b[b != 0] = False
        b = np.ma.array(b, mask = mask_b)

        h = frame.q[0,:,:]
        h = np.ma.array(h, mask = np.invert(mask_b))

        cs = ax1.contourf(x, y, h, cmap = cmap, vmin = vmin, vmax = vmax, levels = levels, corner_mask = False)

        b = frame.aux[0,:,:]
        mask_b = np.zeros_like(b, dtype = bool)
        mask_b[b == 0] = True
        mask_b[b != 0] = False
        mask_b = change_neighbours_to_false(mask_b)
        b = np.ma.array(b, mask = mask_b)
        b[b == 0.0] = np.max(frame.aux[0,:,:])  

        cs_obstacles = ax1.contourf(x, y, b, colors = 'grey', alpha = 0.2, hatches = ['//'], corner_mask = False)

        timestamp = str(round(frame_number * claw.tfinal / len(claw.frames), 3))
        ax1.set_title(r'$t = $' + timestamp + ' s')
        ax1.contour(cs, colors = 'k')   # Add contours and contour color
        #ax1.set_xlabel(r'$x$ [m]')
        ax1.set_ylabel(r'$y$ [m]')
        ax1.set_ylim(- 3, 3)
        ax1.set_xlim(- 2.5, 5)
  
        return cs, cs_obstacles

    #anim = animation.FuncAnimation(fig, fplot3D, frames = len(claw.frames), interval = 100, repeat = False)
    ax1.tick_params(axis = 'x', width = 1, length = 7)
    ax1.tick_params(axis = 'y', width = 1, length = 7)
    ax1.annotate('(a)', xy = (- 1.5, 2.2), size = 18)
    ax1.tick_params(bottom = True, length = 7, labelbottom = False)
    cb.set_label(r'$h \mathrm{\,\,[m]}$')
    #plt.show()
    #plt.close()
    return 0


# Plot height of spill in 3D
def plot_height_3D(claw, var, ylim = (0, 1.2)):
    fig = plt.figure(figsize=[25, 10])
    ax1 = fig.add_subplot(1, 2, 1, projection='3d')
    
    fills = []
    frame = claw.frames[0]
    b = frame.aux[0,:,:]

    b[b == 0] = np.nan
    
    h = frame.q[0,:,:]

    #surface = np.maximum(b, h + b)
    x, y = frame.state.grid.p_centers    
    obstacle_color = np.where(b > 1, 'black', 0.0)
    
    
    def fplot3D(frame_number):
        ax1.clear()
        ax1.set_box_aspect((1, var.y_dimension / var.x_dimension, 1))
        

        timestamp = str(round(frame_number * claw.tfinal / len(claw.frames), 3))
        ax1.set_title(r'$t = $' + timestamp + ' s')
        
        frame = claw.frames[frame_number]
        b = frame.aux[0,:,:]
        h = frame.q[0,:,:]
        
        #surface = np.maximum(b, h + b)
        h[h == 0] = np.nan

        ax1.xaxis.pane.fill = False # Left pane
        ax1.yaxis.pane.fill = False # Right pane
        ax1.grid(False)
        ax1.set_xticklabels([])
        ax1.set_yticklabels([])
        ax1.set_zticklabels([])
        ax1.set_xticks([]) 
        ax1.set_yticks([]) 
        ax1.set_zticks([])

        # Transparent spines
        ax1.w_xaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
        ax1.w_yaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
        ax1.w_zaxis.line.set_color((1.0, 1.0, 1.0, 0.0))

        # Transparent panes
        ax1.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax1.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax1.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))


        ax1.set_zlim(0, 0.04)#claw.spill_height * 1.2)     
        return ax1.plot_surface(x, y, h, color = 'blue', alpha = 1), ax1.plot_surface(x, y, b, facecolors = obstacle_color, edgecolors = 'none', alpha = 0.6, linewidth = 0, shade = False),

    anim = animation.FuncAnimation(fig, fplot3D, frames = len(claw.frames), interval = 100, repeat = True)

    plt.show()
    #plt.close()
    return 0


# Plot cross section of speed
def plot_speed_cross_section(claw, var, ylim = (0, 1.2)):
    fig = plt.figure(figsize=[12 * (var.y_dimension / var.x_dimension), 12 * (var.x_dimension / var.y_dimension)])
    ax1 = fig.add_subplot(111)
    frame = claw.frames[0]
    with np.errstate(all = 'ignore'):     # Stop Python from printing RuntimeWarning 
        speed = np.abs(np.where(frame.q[0,:,:] == 0, 0, frame.q[1,:,:] / frame.q[0,:,:]))

    x, y = frame.state.grid.p_centers    
    slice = int(claw.mx / 2)    

    ax1.plot(x[:, 0], speed[:, slice])
    ax1.set_xlim(- var.x_dimension, var.x_dimension)
    ax1.set_ylim(0, 1.1 * np.max(speed[:, slice]))

    def fplot(frame_number):
        ax1.clear()
        frame = claw.frames[frame_number]
        with np.errstate(all='ignore'):
            speed = np.abs(np.where(frame.q[0,:,:] == 0, 0, frame.q[1,:,:] / frame.q[0,:,:]))
        plot = ax1.plot(x[:, 0], speed[:, slice])
        timestamp = str(round(frame_number * claw.tfinal / len(claw.frames), 3))
        ax1.set_title(r'$t = $' + timestamp + ' s')
        return plot,

    anim = animation.FuncAnimation(fig, fplot, frames = len(claw.frames), interval = 100, repeat = True)
    plt.close()
    return 0