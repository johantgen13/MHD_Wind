#!/usr/bin/env python
#
# This python file is used to animate the checkpoint file data from
# the RUST simulation.
#
# Author: Brayden JoHantgen
# Last Update: 5/7/2026

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

def animation_func(i, directory = "time_step_files"):
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

    plt.title('t='+str(round(frame_df["Values"][1][0], 3)), fontsize=20)

    cells = frame_df["Values"][0][0]
    rnames = np.linspace(0.0, 1.0, (int(cells+1)))
    r = 0.5 * (rnames[:-1] + rnames[1:])

    animated_plot.set_data(r, frame_df["Values"][3])


# Reading the Files
directory = "time_step_files"
lst = os.listdir(directory)
file_num = len(lst) - 1

fig, ax = plt.subplots(1,1)
ax = plotparams(ax)

animated_plot, = ax.plot([], [])
ax.set_xlim(-0.05,1.05)
ax.set_ylim(-0.05,1.05)
ax.set_xlabel('x', fontsize=15)
ax.set_ylabel(r'$\rho$', fontsize=15)

animation = FuncAnimation(fig=fig, 
                    func=animation_func, 
                    frames=file_num,
                    interval=50,
                    repeat=False
                    )

plt.show()