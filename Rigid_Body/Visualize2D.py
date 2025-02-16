## Copyright@2024 Shuoyu Yue
## Email: sy481@cam.ac.uk
## Created on 06/11/2024
## Description: Used to visualize the simulation results of written assignment

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os


## Global Parameters ##
case_id = 5
var_id = 0  # var_id: 0-Density, 1-VelocityX, 2-VelocityY, 3-Pressure
t = 1.0

if case_id in [1, 2, 3, 4]:
    x0, x1 = 0.0, 1.0
    y0, y1 = 0.0, 1.0
    fig_len_x = 7.5
    nCellsX, nCellsY = 500, 500

if case_id in [5, 6]:
    x0, x1 = 0.0, 2.0
    y0, y1 = 0.0, 1.0
    fig_len_x = 15
    nCellsX, nCellsY = 500, 500

if case_id in [7]:
    x0, x1 = 0.0, 1.0
    y0, y1 = 0.0, 1.0
    fig_len_x = 7.5
    nCellsX, nCellsY = 500, 500

dx = (x1 - x0) / nCellsX
dy = (y1 - y0) / nCellsY


## visualization preparation
fig, ax = plt.subplots(1, 1, figsize=(fig_len_x, 6))
manager = plt.get_current_fig_manager()
manager.window.wm_geometry("+300+60")
label_list = ['Density', 'VelocityX', 'VelocityY', 'Pressure']


## extract files
folder_path = "res/Case_%d/" % case_id
# ites = []
# for file_name in os.listdir(folder_path):
#     ites.append(int(file_name[:-4].split("=")[-1]))
# ites = sorted(ites)


## visualization function
def visualize_single(cur_ax, var_id, t, animating = False):
    
    file_name_1 = "T=%d.txt" % t
    file_path_1 = os.path.join(folder_path, file_name_1)
    file_name_2 = "T=%.1f.txt" % t
    file_path_2 = os.path.join(folder_path, file_name_2)
    file_name_3 = "T=%.2f.txt" % t
    file_path_3 = os.path.join(folder_path, file_name_3)

    ## data storage
    rho = np.zeros((nCellsY, nCellsX))
    vx = np.zeros((nCellsY, nCellsX))
    vy = np.zeros((nCellsY, nCellsX))
    p = np.zeros((nCellsY, nCellsX))
    phi = np.zeros((nCellsY, nCellsX))
    
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
    x_ticks = np.linspace(x0, x1, 11).round(1)
    y_ticks = np.linspace(y0, y1, 11).round(1)
    y_ticks = np.flipud(y_ticks)
    
    # visualization
    cur_ax.cla()
    data = data_list[var_id]
    sns.heatmap(data, ax=cur_ax, cbar=not(animating))
        
    if animating:
        fig.subplots_adjust(left=0.1, right=0.8, bottom=0.1, top=0.95)
        plt.pause(0.01)
    else:
        cur_ax.set_title(label_list[var_id] + ", Time = %.2fs" % t)
        cur_ax.set_xlabel('X')
        cur_ax.set_ylabel('Y')

        # change coordinates
        xtick_positions = np.linspace(0, nCellsX, 11)
        ytick_positions = np.linspace(0, nCellsY, 11)

        cur_ax.set_xticks(xtick_positions)
        cur_ax.set_xticklabels(x_ticks)

        cur_ax.set_yticks(ytick_positions)
        cur_ax.set_yticklabels(y_ticks)

        plt.savefig("res/Case_%d_%s_T=%.2fs.png" % (case_id, label_list[var_id], t))
        plt.show()


if __name__ == "__main__":

    # single plot
    visualize_single(ax, var_id, t)
    
    # # animation
    # for ite in ites:
    #     visualize_single(ax, var_id, t, animating = True)

