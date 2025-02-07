## Copyright @ Shuoyu Yue
## Email: sy481@cam.ac.uk
## Created on 17/01/2025
## Description: Used to visualize the simulation results of practicals

import matplotlib.pyplot as plt
import math
import os

case_id = 5
folder_path = "res/case_%d" % (case_id)
ites = []

fig, axes = plt.subplots(3, 1, figsize = (8, 9))
manager = plt.get_current_fig_manager()
manager.window.wm_geometry("+1150+100")

label_list = ['Density', 'Velocity', 'Pressure']
color_list = ['r', 'g', 'b']

for file_name in os.listdir(folder_path):
    ites.append(int(file_name[:-4].split("=")[-1]))
ites = sorted(ites)

for ite in ites:
    file_name = "ite=%d.txt" % ite
    file_path = os.path.join(folder_path, file_name)
    file = open(file_path)
    contents = file.readlines()
    file.close()

    # extract data from file
    xs, rho_1, v_1, p_1, rho_2, v_2, p_2, phi = [], [], [], [], [], [], [], []
    for line in contents:
        data = [float(i) for i in line[:-1].split(', ')]
        xs.append(data[0])
        rho_1.append(data[1])
        v_1.append(data[2])
        p_1.append(data[3])
        rho_2.append(data[4])
        v_2.append(data[5])
        p_2.append(data[6])
        phi.append(data[7])

    # get real material data
    rho, v, p = [], [], []
    for i in range(len(phi)):
        if phi[i] <= 0:
            rho.append(rho_1[i])
            v.append(v_1[i])
            p.append(p_1[i])
        else:
            rho.append(rho_2[i])
            v.append(v_2[i])
            p.append(p_2[i])
    
    # visualise
    ys_list = [rho, v, p]
    for i, ax in enumerate(axes.flat):
        ax.cla()
        ax.plot(xs, ys_list[i], '+', markersize = 5, color = color_list[i]) 

        ax.grid(alpha = 0.5)
        ax.set_xlabel('X')
        ax.set_ylabel(label_list[i])

    plt.tight_layout()
    # plt.show()
    # plt.savefig("res/Figures/" + file_name[:-4] + ".png")
    plt.pause(0.01)