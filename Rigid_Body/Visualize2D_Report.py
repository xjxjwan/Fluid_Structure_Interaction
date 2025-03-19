## Copyright@2024 Shuoyu Yue
## Email: sy481@cam.ac.uk
## Created on 21/02/2025
## Description: Used to visualize the simulation results of written assignment

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import math


## Global Parameters ##
case_id = 7

if case_id in [1, 2, 3, 4]:
    x0, x1 = 0.0, 1.0
    y0, y1 = 0.0, 1.0
    nCellsX, nCellsY = 500, 500
    t_list = [0.2, 0.4]
    var_id_list = [0, 1, 2, 3]  # var_id: 0-Density, 1-Velocity, 2-Pressure, 3-Density Gradient

if case_id in [5, 6]:
    x0, x1 = 0.0, 2.0
    y0, y1 = 0.0, 1.0
    nCellsX, nCellsY = 500, 500
    t_list = [0.5, 1.0]
    var_id_list = [0, 1, 2, 3]  # var_id: 0-Density, 1-Velocity, 2-Pressure, 3-Density Gradient

if case_id in [7]:
    x0, x1 = 0.0, 1.0
    y0, y1 = 0.0, 1.0
    nCellsX, nCellsY = 500, 500
    t_list = [0.25, 0.5, 0.75, 1.0]
    var_id_list = [0, 3]  # var_id: 0-Density, 3-Density Gradient

dx = (x1 - x0) / nCellsX
dy = (y1 - y0) / nCellsY


## visualization preparation
if case_id in [1, 2, 3, 4]:
    fig, axes = plt.subplots(2, 4, figsize=(33, 13), dpi=300)
if case_id in [5, 6]:
    fig, axes = plt.subplots(2, 4, figsize=(32, 7), dpi=300, gridspec_kw={"hspace": 0.4, "wspace": 0.1})
if case_id in [7]:
    fig, axes = plt.subplots(2, 4, figsize=(33, 13), dpi=300)

label_list = ['Density', 'Velocity', 'Pressure', 'Density Gradient']


## extract files
if case_id in [1, 2, 3, 4, 5, 6]:
    folder_path = "res/Case_%d/" % case_id
if case_id in [7]:
    folder_path = "res_oscillation/Case_%d/" % case_id


## visualization function
def visualize_single(cur_ax, var_id, t):
    
    file_name_1 = "T=%d.txt" % t
    file_path_1 = os.path.join(folder_path, file_name_1)
    file_name_2 = "T=%.1f.txt" % t
    file_path_2 = os.path.join(folder_path, file_name_2)
    file_name_3 = "T=%.2f.txt" % t
    file_path_3 = os.path.join(folder_path, file_name_3)

    ## data storage
    rho = np.zeros((nCellsY, nCellsX))
    v = np.zeros((nCellsY, nCellsX))
    p = np.zeros((nCellsY, nCellsX))
    phi = np.zeros((nCellsY, nCellsX))
    rho_grad = np.zeros((nCellsY, nCellsX))
    
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
        x, y, cur_rho, cur_vx, cur_vy, cur_p, cur_phi = data
        
        col = int(round((x - x0) / dx - 0.5))
        row = int(round((nCellsY - 1) - ((y - y0) / dy - 0.5)))

        rho[row, col] = cur_rho
        v[row, col] = pow(pow(cur_vx, 2) + pow(cur_vy, 2), 0.5)
        p[row, col] = cur_p
        phi[row, col] = cur_phi

    # get real material data (real material: phi > 0)
    rho = rho * (phi > 0) + 0 * (phi < 0)
    v = v * (phi > 0) + 0 * (phi < 0)
    p = p * (phi > 0) + 0 * (phi < 0)

    # # determine interface cells
    # interface_location = np.zeros((nCellsY, nCellsX))
    # for i in range(1, nCellsY - 1):
    #     for j in range(1, nCellsX - 1):
    #         cur_phi = phi[i][j]
    #         phi_up, phi_down = phi[i][j + 1], phi[i][j - 1]
    #         phi_right, phi_left = phi[i + 1][j], phi[i - 1][j]
    #         if (cur_phi * phi_up < 0 or cur_phi * phi_down < 0 or cur_phi * phi_right < 0 or cur_phi * phi_left < 0):
    #             interface_location[i][j] = 1

    # calculate density gradient
    for row in range(1, nCellsY - 1):
        for col in range(1, nCellsX - 1):
            if rho[row, col] != 0:         
                diff_x = (rho[row, col + 1] - rho[row, col - 1]) / (2 * dx)
                diff_y = (rho[row + 1, col] - rho[row - 1, col]) / (2 * dy)
                abs_grad = pow(pow(diff_x, 2) + pow(diff_y, 2), 0.5)
                rho_grad[row, col] = math.exp(-20 * abs_grad / (1000 * rho[row, col]))
    
    # domain boundary condition
    rho_grad[0, :] = rho_grad[1, :]
    rho_grad[nCellsY - 1, :] = rho_grad[nCellsY - 2, :]
    rho_grad[:, 0] = rho_grad[:, 1]
    rho_grad[:, nCellsX - 1] = rho_grad[:, nCellsX - 2]

    # data for visualization
    data_list = [rho, v, p, rho_grad]

    # Define the actual coordinate values for ticks
    x_ticks = np.linspace(x0, x1, 6).round(1)
    y_ticks = np.linspace(y0, y1, 6).round(1)
    y_ticks = np.flipud(y_ticks)

    # visualization
    data = data_list[var_id]
    # rigid body
    mask = phi < 0
    sns.heatmap(np.ones_like(data) * -1, ax=cur_ax, mask=~mask, cmap=["#D3D3D3"], cbar=False)
    # fluid
    sns.heatmap(data, ax=cur_ax, mask=mask)
        
    cur_ax.set_title(label_list[var_id] + ", Time = %.2fs" % t, size=20)
    cur_ax.set_xlabel('X', size=20)
    cur_ax.set_ylabel('Y', size=20)

    # change coordinates
    xtick_positions = np.linspace(0, nCellsX, 6)
    ytick_positions = np.linspace(0, nCellsY, 6)

    cur_ax.set_xticks(xtick_positions)
    cur_ax.set_xticklabels(x_ticks, rotation = 0)

    cur_ax.set_yticks(ytick_positions)
    cur_ax.set_yticklabels(y_ticks)


if __name__ == "__main__":

    # single plot
    if case_id in [1, 2, 3, 4, 5, 6]:
        for t in t_list:
            cur_axes = axes[t_list.index(t)]
            for i, ax in enumerate(cur_axes.flat):
                var_id = var_id_list[i]
                visualize_single(ax, var_id, t)
    
    if case_id == 7:
        for var_id in var_id_list:
            cur_axes = axes[var_id_list.index(var_id)]
            for i, ax in enumerate(cur_axes.flat):
                t = t_list[i]
                visualize_single(ax, var_id, t)
    
    plt.savefig("D:/Study_Master/WrittenAssignment/Writing/Figure1%d_RigidBodyTest_Case%d.png" % (case_id, case_id), 
        bbox_inches='tight', pad_inches=0.1)
    plt.subplots_adjust(left=0.05, right=0.95)
    # plt.tight_layout()
    # plt.show()

