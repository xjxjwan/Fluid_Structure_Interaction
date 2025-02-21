## Copyright@2024 Shuoyu Yue
## Email: sy481@cam.ac.uk
## Created on 06/11/2024
## Description: Used to visualize the simulation results of written assignment

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os


## Global Parameters ##
case_id = 3
var_id = 2  # var_id: 0-Density, 1-Velocity, 2-Pressure, 3-Energy
t = 0.25
gama = 1.4

if case_id in [1, 2]:
    x0, x1 = 0.0, 1.0
    y0, y1 = 0.0, 1.0
    nCellsX, nCellsY = 200, 200

if case_id in [3]:
    x0, x1 = 0.0, 2.0
    y0, y1 = 0.0, 2.0
    nCellsX, nCellsY = 400, 400

dx = (x1 - x0) / nCellsX
dy = (y1 - y0) / nCellsY


## visualization preparation
fig = plt.figure(figsize = (8, 8))
ax = fig.add_subplot(111, projection='3d')
label_list = ['Density', 'Velocity', 'Pressure', 'Internal Energy']


## extract files
folder_path = "res/Case_%d/" % case_id


## visualization function
def visualize_single(cur_ax, var_id, t):
    
    file_name_1 = "T=%d.txt" % t
    file_path_1 = os.path.join(folder_path, file_name_1)
    file_name_2 = "T=%.1f.txt" % t
    file_path_2 = os.path.join(folder_path, file_name_2)
    file_name_3 = "T=%.2f.txt" % t
    file_path_3 = os.path.join(folder_path, file_name_3)

    ## data storage
    xs = np.zeros((nCellsY, nCellsX))
    ys = np.zeros((nCellsY, nCellsX))
    rho = np.zeros((nCellsY, nCellsX))
    v = np.zeros((nCellsY, nCellsX))
    p = np.zeros((nCellsY, nCellsX))
    E = np.zeros((nCellsY, nCellsX))
    
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

    for line in contents:

        data = [float(i) for i in line.strip().split(', ')]
        x, y, cur_rho, cur_vx, cur_vy, cur_p = data
        
        col = int(round((x - x0) / dx - 0.5))
        row = int(round((nCellsY - 1) - ((y - y0) / dy - 0.5)))

        xs[row, col] = x
        ys[row, col] = y
        rho[row, col] = cur_rho
        cur_v = pow(pow(cur_vx, 2) + pow(cur_vy, 2), 0.5)
        v[row, col] = cur_v
        p[row, col] = cur_p
        E[row, col] = cur_p / (gama - 1) / cur_rho

    data_list = [rho, v, p, E]
    
    # visualization
    data = data_list[var_id]
    cur_ax.plot_surface(xs, ys, data, cmap='viridis')
        
    cur_ax.set_zlabel(label_list[var_id])
    cur_ax.set_xlabel('X')
    cur_ax.set_ylabel('Y')

    # # change coordinates
    # # Define the actual coordinate values for ticks
    # x_ticks = np.linspace(x0, x1, 11).round(1)
    # y_ticks = np.linspace(y0, y1, 11).round(1)
    # y_ticks = np.flipud(y_ticks)

    # xtick_positions = np.linspace(0, nCellsX, 11)
    # ytick_positions = np.linspace(0, nCellsY, 11)

    # cur_ax.set_xticks(xtick_positions)
    # cur_ax.set_xticklabels(x_ticks)

    # cur_ax.set_yticks(ytick_positions)
    # cur_ax.set_yticklabels(y_ticks)

    # plt.savefig("res/Case_%d_%s_T=%.2fs.png" % (case_id, label_list[var_id], t))
    plt.show()


if __name__ == "__main__":

    # single plot
    visualize_single(ax, var_id, t)
    
