import os
import subprocess
import re
import numpy as np
import matplotlib.pyplot as plt

def generate_plate_input(N=8, E=2e11, nu=0.3, thickness=0.01, force_val=-1.0):
    with open(f"./data/plate_{N}x{N}.dat", "w") as f:
        # 标题
        f.write(f"Plate Element Test - {N*N} Elements ({N}x{N} Grid)\n")

        # 节点编号函数
        node_id = lambda i, j: i * (N + 1) + j + 1
        total_nodes = (N + 1) * (N + 1)
        f.write(f"{total_nodes:<4d}   1   1   1\n")

        # 节点定义
        for i in range(N + 1):
            for j in range(N + 1):
                nid = node_id(i, j)
                x, y, z = j * 8.0 / N, i * 8.0 / N, 0.0
                # 边界设为 1 1 1，其余为 0 0 0
                bc = "1   1   1   1   1   1" if i == 0 or j == 0 or i == N or j == N else "1   1   0   0   0   1"
                f.write(f"{nid:<5d} {bc}   {x:.5f}   {y:.5f}   {z:.5f}\n")

        # 内部节点加载力（非边界一圈）
        load_nodes = []
        for i in range(1, N):  # i = 行
            for j in range(1, N):  # j = 列
                nid = node_id(i, j)
                load_nodes.append(nid)
        force_val = force_val * 8 * 8 / len(load_nodes)  # 平均分配力

        f.write(f"1   {len(load_nodes)}\n")
        for nid in load_nodes:
            f.write(f"{nid:<5d}   3   {force_val:.5f}\n")

        # 材料与单元信息
        f.write(f"3   {N*N}   1\n")
        f.write(f"1   {E:.1e}   {nu:.2f}   {thickness:.3f}\n")

        # 单元定义（4 节点单元）
        eid = 1
        for i in range(N):
            for j in range(N):
                n1 = node_id(i, j)
                n2 = node_id(i, j + 1)
                n3 = node_id(i + 1, j + 1)
                n4 = node_id(i + 1, j)
                f.write(f"{eid:<5d}   {n1}   {n2}   {n3}   {n4}   1\n")
                eid += 1


# 网格尺寸
N_list = [2,4,8,16,32,64]
L2_errors = []
h_list = []

for N in N_list:
    # 生成输入文件
    input_path = f"./data/plate_{N}x{N}.dat"
    output_path = f"./data/plate_{N}x{N}.out"
    generate_plate_input(N)
    # 运行仿真
    print(f"Running stap++ for {N}x{N}...")
    subprocess.run(["./stap++", input_path], check=True)
    # 解析输出文件，提取L2相对误差
    with open(output_path, "r") as f:
        content = f.read()
        match = re.search(r"L2 相对误差 \(w\): ([0-9.eE+-]+)", content)
        if match:
            err = float(match.group(1))
            L2_errors.append(err)
            h_list.append(8.0/N)  # 单元尺寸
            print(f"N={N}, h={8.0/N}, L2相对误差={err}")
        else:
            print(f"未找到L2误差: {output_path}")
            L2_errors.append(np.nan)
            h_list.append(8.0/N)

# 绘制收敛图
log_h = np.log(h_list)
log_err = np.log(L2_errors)
plt.figure()
plt.plot(log_h, log_err, 'o-', label='S4R Element')
plt.xlabel('log(h)')
plt.ylabel('log(L2 Error)')
plt.title('S4R Element L2 Error Convergence')
slope, _ = np.polyfit(log_h, log_err, 1)
plt.plot(log_h, slope*log_h, '--', label=f'Slope={slope:.2f}')
plt.legend()
plt.grid(True)
plt.show()
