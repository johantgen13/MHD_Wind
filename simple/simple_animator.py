#!/usr/bin/env python
#
# This python file is used to animate the checkpoint file data from
# the 1D nonrelativistic RUST simulation.
#
# Author: Brayden JoHantgen
# Last Update: 6/9/2026

# Importing
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
from matplotlib.animation import FuncAnimation
from matplotlib.animation import PillowWriter

# Functions
def plotparams(ax, labelsize=10):
    '''
    Basic plot params
        :param ax: axes to modify
        :type ax: matplotlib axes object
        :returns: modified matplotlib axes object
    '''
    ax.minorticks_on()
    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')
    ax.tick_params(direction='in', which='both', labelsize=labelsize)
    ax.tick_params('both', length=8, width=1.8, which='major')
    ax.tick_params('both', length=4, width=1, which='minor')
    for axis in ['top', 'bottom', 'left', 'right']:
        ax.spines[axis].set_linewidth(1.5)
    return ax

def line_num(num, directory = "time_step_files"):
    '''
    '''
    file_path = directory+"/"+str(num)+".txt"
    with open(file_path, 'r') as f:
        line_count = len(f.readlines())
    return line_count


def read_txt_files(num, directory = "time_step_files"):
    '''
    '''
    file_path = directory+"/"+str(num)+".txt"

    file_data = {}
    line_count = line_num(num, directory = directory)

    with open(file_path, 'r') as f:
        line_list = f.readlines()
        for i in range(line_count):
            line_dict_split = line_list[i].split(": ")
            float_list = line_dict_split[1].split(" ")
            float_data = []
            for j in range(len(float_list)-1):
                float_data.append(float(float_list[j]))
            line_dict = {line_dict_split[0]: float_data}
            file_data.update(line_dict)
    df = pd.DataFrame(list(file_data.items()), columns = ['Physics', 'Values'])
    return df

def two_dimension_list(td_list, xnum, ynum):
    '''
    '''
    two_d = []
    for i in range(xnum):
        fill = []
        for j in range(ynum):
            fill.append(td_list[(i * ynum + j)])
        two_d.append(fill)
    transpose = [list(row) for row in zip(*two_d)]
    return transpose

def animation_func_1D_FO(i, directory = "time_step_files"):
    '''
    Input:
        i: the frame number
    Output:
        animated_plot, the data of the plot of the corresponding frame 
                        number
    Description:
        This function takes the frame number and reads the corresponding
        pickle file. It then loads and plots the data to produce one frame
        of the animation.
    '''
    frame_df = read_txt_files(i, directory = directory)

    #plt.title('t='+str(round(frame_df["Values"][1][0], 3)), fontsize=20)

    cells = frame_df["Values"][0][0]
    rnames = np.linspace(0.0, 1.0, (int(cells)))
    r = 0.5 * (rnames[:-1] + rnames[1:])

    animated_plot_1.set_data(r, frame_df["Values"][3])
    animated_plot_2.set_data(r, frame_df["Values"][4])
    animated_plot_3.set_data(r, frame_df["Values"][5])
    animated_plot_4.set_data(r, frame_df["Values"][8])
    animated_plot_5.set_data(r, frame_df["Values"][2])

def animate_2d(i, directory = "time_step_files"):
    '''
    '''
    frame_df = read_txt_files(i, directory = directory)

    plt.title("t="+str(round(frame_df["Values"][2][0], 4)), fontsize=20)

    rho_list = two_dimension_list(frame_df["Values"][4], 512, 512)
    A.set_data(rho_list)


# Reading the Files
directory = "time_step_files"
lst = os.listdir(directory)
file_num = len(lst) - 1

fig, ax = plt.subplots(1,1)#, size=(5,5))
ax = plotparams(ax)

A = ax.imshow(np.zeros((512, 512)), cmap='plasma', interpolation='nearest',
                    vmax=3.0, extent = [-0.5, 0.5, -0.5, 0.5])

animation = FuncAnimation(fig=fig, 
                    func=animate_2d, 
                    frames=file_num,
                    interval=50,
                    repeat=False,
                    )

ax.set_xlabel(r'$x$', fontsize=15)
ax.set_ylabel(r'$y$', fontsize=15)

cbar = fig.colorbar(A, ax = ax)
cbar.set_label(r'$\rho$', size=15)
cbar.ax.tick_params(labelsize=10)
plt.show()

# 1D Sims
#fig = plt.figure(figsize=(7,7))
#gs = fig.add_gridspec(3, 2, wspace=0.3)
#ax = gs.subplots(sharex=False, sharey=False)

#ax[0,0] = plotparams(ax[0,0])
#ax[0,1].set_visible(False)
#ax[1,0] = plotparams(ax[1,0])
#ax[1,1] = plotparams(ax[1,1])
#ax[2,0] = plotparams(ax[2,0])
#ax[2,1] = plotparams(ax[2,1])

#animated_plot_1, = ax[0,0].plot([], [])
#animated_plot_2, = ax[1,0].plot([], [])
#animated_plot_3, = ax[1,1].plot([], [])
#animated_plot_4, = ax[2,0].plot([], [])
#animated_plot_5, = ax[2,1].plot([], [])

#ax[0,0].set_xlim(-0.1,1.1)
#ax[0,0].set_ylim(0.0,1.1)

#ax[1,0].set_xlim(-0.1,1.1)
#ax[1,0].set_ylim(-0.4,0.8)

#ax[1,1].set_xlim(-0.1,1.1)
#ax[1,1].set_ylim(-1.8,0.2)

#ax[2,0].set_xlim(-0.1,1.1)
#ax[2,0].set_ylim(-1.25,1.25)

#ax[2,1].set_xlim(-0.1,1.1)
#ax[2,1].set_ylim(-0.05,1.15)

#ax[2,0].set_xlabel(r'$x$', fontsize=15)
#ax[2,1].set_xlabel(r'$x$', fontsize=15)

#ax[0,0].set_ylabel(r'$\rho$', fontsize=15)
#ax[1,0].set_ylabel(r'$v_x$', fontsize=15)
#ax[1,1].set_ylabel(r'$v_y$', fontsize=15)
#ax[2,0].set_ylabel(r'$B_y$', fontsize=15)
#ax[2,1].set_ylabel(r'$P$', fontsize=15)

#ax[0,0].vlines(0.5, 0.0, 1.1, color="black", alpha=0.3)
#ax[1,0].vlines(0.5, -0.4, 0.8, color="black", alpha=0.3)
#ax[1,1].vlines(0.5, -1.8, 0.2, color="black", alpha=0.3)
#ax[2,0].vlines(0.5, -1.25, 1.25, color="black", alpha=0.3)
#ax[2,1].vlines(0.5, -0.05, 1.15, color="black", alpha=0.3)

#animation = FuncAnimation(fig=fig, 
#                    func=animation_func_1D_FO, 
#                    frames=file_num,
#                    interval=50,
#                    repeat=False
#                    )

#plt.show()