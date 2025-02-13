## Copyright@2024 Shuoyu Yue
## Email: sy481@cam.ac.uk
## Created on 06/11/2024
## Description: Used to visualize the simulation results of written assignment

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os


## Global Parameters ##
case_id = 1
nCells = 100
x0, x1 = 0.0, 1.0
y0, y1 = 0.0, 1.0
dx = (x1 - x0) / nCells
dy = (y1 - y0) / nCells


## visualization preparation
fig, ax = plt.subplots(1, 1, figsize=(7.5, 6))
manager = plt.get_current_fig_manager()
manager.window.wm_geometry("+600+60")
label_list = ['Density', 'VelocityX', 'VelocityY', 'Pressure']


## extract files
folder_path = "res/case_%d" % case_id
ites = []
for file_name in os.listdir(folder_path):
    ites.append(int(file_name[:-4].split("=")[-1]))
ites = sorted(ites)


## visualization function
def visualize_single(cur_ax, var_id, ite):
    
    file_name = "ite=%d.txt" % ite
    file_path = os.path.join(folder_path, file_name)

    ## data storage
    rho = np.zeros((nCells, nCells))
    vx = np.zeros((nCells, nCells))
    vy = np.zeros((nCells, nCells))
    p = np.zeros((nCells, nCells))
    phi = np.zeros((nCells, nCells))
    
    with open(file_path) as file:
        contents = file.readlines()

    for line in contents:

        data = [float(i) for i in line.strip().split(', ')]
        x, y, cur_rho, cur_vx, cur_vy, cur_p, cur_phi, time = data
        
        col = int(round((x - x0) / dx - 0.5))
        row = int(round((nCells - 1) - ((y - y0) / dy - 0.5)))

        rho[row, col] = cur_rho
        vx[row, col] = cur_vx
        vy[row, col] = cur_vy
        p[row, col] = cur_p
        phi[row, col] = cur_phi

    # get real material data (real material: phi > 0)
    rho = rho * (phi > 0) + 0 * (phi < 0)
    vx = vx * (phi > 0) + 0 * (phi < 0)
    vy = vy * (phi > 0) + 0 * (phi < 0)
    p = p * (phi > 0) + 0 * (phi < 0)
    data_list = [rho, vx, vy, p]

    # Define the actual coordinate values for ticks
    x_ticks = np.linspace(x0, x1, 10).round(1)
    y_ticks = np.linspace(y0, y1, 10).round(1)
    y_ticks = np.flipud(y_ticks)
    
    # visualization
    cur_ax.cla()
    data = data_list[var_id]
    sns.heatmap(data, ax=cur_ax, cbar=True)

    cur_ax.set_title(label_list[var_id] + ", Time = %.3fs" % time)
    cur_ax.set_xlabel('X')
    cur_ax.set_ylabel('Y')

    # change coordinates
    xtick_positions = np.linspace(0, nCells - 1, 10)
    ytick_positions = np.linspace(0, nCells - 1, 10)

    cur_ax.set_xticks(xtick_positions)
    cur_ax.set_xticklabels(x_ticks)

    cur_ax.set_yticks(ytick_positions)
    cur_ax.set_yticklabels(y_ticks)
    
    # plt.tight_layout()
    # plt.pause(0.01)
    plt.show()


if __name__ == "__main__":

    # var_id: 0-Density, 1-VelocityX, 2-VelocityY, 3-Pressure
    visualize_single(ax, 0, ites[-1])
    
    # animation
    # for ite in ites:
    #     visualize_single(ite)

