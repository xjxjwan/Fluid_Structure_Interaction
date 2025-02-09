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


## data storage
rho_1, rho_2 = np.zeros((nCells, nCells)), np.zeros((nCells, nCells))
vx_1, vx_2 = np.zeros((nCells, nCells)), np.zeros((nCells, nCells))
vy_1, vy_2 = np.zeros((nCells, nCells)), np.zeros((nCells, nCells))
p_1, p_2 = np.zeros((nCells, nCells)), np.zeros((nCells, nCells))
phi = np.zeros((nCells, nCells))


## visualization preparation
fig, axes = plt.subplots(2, 2, figsize=(9, 9))
manager = plt.get_current_fig_manager()
manager.window.wm_geometry("+1000+80")
label_list = ['Density', 'VelocityX', 'VelocityY', 'Pressure']
# 添加独立 colorbar 轴（放置在子图右侧）
cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])


## extract files
folder_path = "res/case_%d" % case_id
ites = []
for file_name in os.listdir(folder_path):
    ites.append(int(file_name[:-4].split("=")[-1]))
ites = sorted(ites)


## visualization
for ite in ites:
    
    file_name = "ite=%d.txt" % ite
    file_path = os.path.join(folder_path, file_name)
    
    with open(file_path) as file:
        contents = file.readlines()

    for line in contents:

        data = [float(i) for i in line.strip().split(', ')]
        x, y, cur_rho_1, cur_vx_1, cur_vy_1, cur_p_1, cur_rho_2, cur_vx_2, cur_vy_2, cur_p_2, cur_phi, time = data
        
        col = int(round((x - x0) / dx - 0.5))
        row = int(round((nCells - 1) - ((y - y0) / dy - 0.5)))

        rho_1[row, col], rho_2[row, col] = cur_rho_1, cur_rho_2
        vx_1[row, col], vx_2[row, col] = cur_vx_1, cur_vx_2
        vy_1[row, col], vy_2[row, col] = cur_vy_1, cur_vy_2
        p_1[row, col], p_2[row, col] = cur_p_1, cur_p_2
        phi[row, col] = cur_phi

    # get real material data (material 1: phi < 0)
    rho = rho_1 * (phi < 0) + rho_2 * (phi > 0)
    vx = vx_1 * (phi < 0) + vx_2 * (phi > 0)
    vy = vy_1 * (phi < 0) + vy_2 * (phi > 0)
    p = p_1 * (phi < 0) + p_2 * (phi > 0)

    data_list = [rho, vx, vy, p]

    # Define the actual coordinate values for ticks
    origin_x_ticks = np.linspace(0, nCells, 10, dtype=int)
    x_ticks = np.linspace(x0, x1, 10).round(1)
    origin_y_ticks = np.linspace(0, nCells, 10, dtype=int)
    y_ticks = np.linspace(y0, y1, 10).round(1)
    y_ticks = np.flipud(y_ticks)
    
    for i, cur_ax in enumerate(axes.flat):
        cur_ax.cla()

        # 获取当前数据
        data = data_list[i]

        # 绘制 heatmap，关闭自动刻度调整，确保轴不被覆盖
        sns.heatmap(
            data, ax=cur_ax, cbar=(i == 0), cbar_ax=(cbar_ax if i == 0 else None),
            xticklabels=False, yticklabels=False  # 禁用自动刻度
        )

        cur_ax.set_title(label_list[i])
        cur_ax.set_xlabel('X')
        cur_ax.set_ylabel('Y')

        # **确保所有子图都使用相同的刻度**
        xtick_positions = np.linspace(0, nCells - 1, 10)
        ytick_positions = np.linspace(0, nCells - 1, 10)

        cur_ax.set_xticks(xtick_positions)  # **强制应用网格索引刻度**
        cur_ax.set_xticklabels(x_ticks)  # **确保是物理坐标 0-1**

        cur_ax.set_yticks(ytick_positions)  # **确保 y 轴刻度与网格对齐**
        cur_ax.set_yticklabels(y_ticks)  # **确保是物理坐标 0-1**

    fig.suptitle(f'Time = {time:.5f}s  Case = 1 (Sod Test in X-direction)', fontsize=14)
    fig.subplots_adjust(left=0.1, right=0.88, bottom=0.1, top=0.9, wspace=0.3, hspace=0.3)
    plt.pause(0.05)
