import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# 全局参数
case_id = 1
nCells = 100
x0, x1 = 0.0, 1.0
y0, y1 = -1.0, 1.0
dx = (x1 - x0) / nCells
dy = (y1 - y0) / nCells

# 创建图形和 3D 轴，仅创建一次
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# 定义一个函数用于在同一窗口中更新 3D 表面图
def plot_surface(data, title, zlabel, cmap='viridis'):
    ax.cla()  # 清除当前轴内容
    ax.plot_surface(X, Y, data, cmap=cmap)
    ax.set_xlabel('r')
    ax.set_ylabel('z')
    ax.set_zlabel(zlabel)
    ax.set_title(title)

    # 保存图片并显示暂停
    plt.savefig(f"res/Figures/p2t{case_id}_tStop={tStop:.2f}_{title}.png")
    plt.pause(0.5)  # 停顿0.5秒以显示更新后的图形

# 时间步列表
time_steps = np.arange(0, 0.25, 0.01)

# 绘图主循环
for tStop in time_steps:
    # 从文件中读取数据
    file_path = f"res/case_{case_id}_tStop={tStop:.2f}.txt"
    with open(file_path, 'r') as file:
        contents = file.readlines()

    # 处理数据
    xs, ys = [], []
    rho, vx, vy, p = [], [], [], []

    for line in contents:
        data = [float(i) for i in line.strip().split(', ')]
        i, j, cur_rho, cur_vx, cur_vy, cur_p = data
        xs.append(x0 + dx * (i + 0.5))
        ys.append(y0 + dy * (j + 0.5))
        rho.append(cur_rho)
        vx.append(cur_vx)
        vy.append(cur_vy)
        p.append(cur_p)

    # 将 x 和 y 转换为网格
    X = np.array(xs).reshape((nCells, nCells))
    Y = np.array(ys).reshape((nCells, nCells))
    rho = np.array(rho).reshape((nCells, nCells))
    vx = np.array(vx).reshape((nCells, nCells))
    vy = np.array(vy).reshape((nCells, nCells))
    p = np.array(p).reshape((nCells, nCells))

    # 绘制密度的三维图像，并在同一窗口中更新
    plot_surface(rho, 'Density', 'Density')
    # 若需要其他变量的绘图，注释可以去除
    # plot_surface(vx, 'VelocityR', 'VelocityR')
    # plot_surface(vy, 'VelocityZ', 'VelocityZ')
    # plot_surface(p, 'Pressure', 'Pressure')
