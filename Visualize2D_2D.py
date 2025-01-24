## Copyright@2024 Shuoyu Yue
## Email: sy481@cam.ac.uk
## Created on 06/11/2024
## Description: Used to visualize the simulation results of practicals in course MEC


import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd


## Global Parameters ##
case_id = 1
nCells = 100
x0, x1 = 0.0, 1.0
y0, y1 = -1.0, 1.0
dx = (x1 - x0) / nCells
dy = (y1 - y0) / nCells


# data extraction
file = open("res/case_%d.txt" % (case_id))
contents = file.readlines()
file.close()


# 2D plot ##
# data process
rho = np.zeros((nCells, nCells))
vx = np.zeros((nCells, nCells))
vy = np.zeros((nCells, nCells))
p = np.zeros((nCells, nCells))

for line in contents:
    data = [float(i) for i in line[:-1].split(', ')]
    i, j, cur_rho, cur_vx, cur_vy, cur_p = data
    i, j = int(i), int(j)
    rho[j][i] = cur_rho
    vx[j][i] = cur_vx
    vy[j][i] = cur_vy
    p[j][i] = cur_p

rho = pd.DataFrame(rho)
vx = pd.DataFrame(vx)
vy = pd.DataFrame(vy)
p = pd.DataFrame(p)
data_list = [rho, vx, vy, p]


# visualization
fig, axes = plt.subplots(2, 2, figsize = (9, 8))
label_list = ['Density', 'VelocityR', 'VelocityZ', 'Pressure']

# Define the actual coordinate values for ticks
origin_x_ticks = np.linspace(0, nCells, 10)
x_ticks = np.linspace(x0, x1, 10).round(1)
origin_y_ticks = np.linspace(0, nCells, 10)
y_ticks = np.linspace(y0, y1, 10).round(1)

for i, cur_ax in enumerate(axes.flat):
    
    sns.heatmap(data_list[i], ax = cur_ax)
    
    cur_ax.set_title(label_list[i])
    cur_ax.set_xlabel('R')
    cur_ax.set_ylabel('Z')
    cur_ax.set_xticks(origin_x_ticks, x_ticks)
    cur_ax.set_yticks(origin_y_ticks, y_ticks)
    cur_ax.invert_yaxis()

plt.tight_layout()
plt.savefig("res/Figures/p2t%d.png" % (case_id))
plt.show()