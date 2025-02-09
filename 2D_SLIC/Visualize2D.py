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
rho = np.zeros((nCells, nCells))
vx = np.zeros((nCells, nCells))
vy = np.zeros((nCells, nCells))
p = np.zeros((nCells, nCells))


## visualization preparation
fig, axes = plt.subplots(2, 2, figsize=(9, 9))
manager = plt.get_current_fig_manager()
manager.window.wm_geometry("+980+60")
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
def visualize(ite):
    
    file_name = "ite=%d.txt" % ite
    file_path = os.path.join(folder_path, file_name)
    
    with open(file_path) as file:
        contents = file.readlines()

    for line in contents:

        data = [float(i) for i in line.strip().split(', ')]
        x, y, cur_rho, cur_vx, cur_vy, cur_p, time = data
        
        col = int(round((x - x0) / dx - 0.5))
        row = int(round((nCells - 1) - ((y - y0) / dy - 0.5)))

        rho[row, col] = cur_rho
        vx[row, col] = cur_vx
        vy[row, col] = cur_vy
        p[row, col] = cur_p

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

        # **让每个 heatmap 拥有自己的 colorbar**
        sns.heatmap(
            data, ax=cur_ax, cbar=True, cbar_kws={"shrink": 0.8},  # `shrink` 控制 colorbar 尺寸
            xticklabels=False, yticklabels=False  # 禁用 heatmap 默认的 tick labels
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

    fig.suptitle(f'Time = {time:.5f}s  Case = 2 (Sod Test in Y-direction)', fontsize=14)
    fig.subplots_adjust(left=0.1, right=0.88, bottom=0.1, top=0.9, wspace=0.3, hspace=0.3)
    # plt.pause(0.05)
    plt.show()


if __name__ == "__main__":

    visualize(ites[-1])
    # for ite in ites:
    #     visualize(ite)

