## Copyright@2024 Shuoyu Yue
## Email: sy481@cam.ac.uk
## Created on 20/02/2025
## Description: Used to visualize the simulation results of written assignment

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os


## Global Parameters ##
case_id = 3
t = 0.25
gama = 1.4
var_id_list = [0, 1, 2, 3]  # density, velocity, pressure, internal energy

x0, x1 = 0.0, 2.0
y0, y1 = 0.0, 2.0
cut_pos = 1.0
nCellsX, nCellsY = 200, 200

dx = (x1 - x0) / nCellsX
dy = (y1 - y0) / nCellsY


## visualization preparation
fig, axes = plt.subplots(2, 2, figsize=(8, 7))
label_list = ['Density', 'Velocity', 'Pressure', 'Specific Internal Energy']


## extract files
folder_path = "res/Case_%d/" % case_id


## visualization function
def visualize_single(cur_ax, var_id, t, cut_pos):
    
    # open file
    file_name_1 = "T=%d.txt" % t
    file_path_1 = os.path.join(folder_path, file_name_1)
    file_name_2 = "T=%.1f.txt" % t
    file_path_2 = os.path.join(folder_path, file_name_2)
    file_name_3 = "T=%.2f.txt" % t
    file_path_3 = os.path.join(folder_path, file_name_3)
    
    try:
        with open(file_path_1) as file:
            contents = file.readlines()
    except:
        try:
            with open(file_path_2) as file:
                contents = file.readlines()
        except:
            with open(file_path_3) as file:
                contents = file.readlines()

    # data storage
    num_data_point = nCellsX // 2
    points = np.arange(0, num_data_point, 1)
    rho, v, p, E = np.zeros(num_data_point), np.zeros(num_data_point), np.zeros(num_data_point), np.zeros(num_data_point)

    for line in contents:

        data = [float(i) for i in line.strip().split(', ')]
        x, y, cur_rho, cur_vx, cur_vy, cur_p = data
        
        if abs(y - cut_pos) < dy and y >= cut_pos and x >= cut_pos:
            index = int(round((x - cut_pos) / dx - 0.5))
            rho[index] = cur_rho
            cur_v = pow(pow(cur_vx, 2) + pow(cur_vy, 2), 0.5)
            v[index] = cur_v
            p[index] = cur_p
            E[index] = cur_p / (gama - 1) / cur_rho

    data_list = [rho, v, p, E]
    
    # visualization
    cur_ax.cla()
    data = data_list[var_id]
    cur_ax.scatter(points, data_list[var_id], marker='+', s=20, c='steelblue')
    cur_ax.set_ylabel(label_list[var_id])
    cur_ax.grid(alpha = 0.2)

    # change coordinates
    cur_ax.set_xlabel('X')
    x_ticks = np.linspace(cut_pos, x1, 11).round(1)
    xtick_positions = np.linspace(0, num_data_point, 11)
    cur_ax.set_xticks(xtick_positions)
    cur_ax.set_xticklabels(x_ticks)


if __name__ == "__main__":

    # single plot
    for i, ax in enumerate(axes.flat):
        var_id = var_id_list[i]
        visualize_single(ax, var_id, t, cut_pos)
    
    plt.tight_layout()
    plt.show()