import os
import subprocess
import re
import numpy as np
import matplotlib.pyplot as plt

def generate_plate_input(N=8, E=2e11, nu=0.3, thickness=0.01, force_val=-1.0):
    # 确保数据目录存在
    os.makedirs("./data", exist_ok=True)
    
    with open(f"./data/plate_{N}x{N}_Q4.dat", "w") as f:
        # 标题
        f.write(f"Plate Element Test - {N*N} Elements ({N}x{N} Grid) - Q4 Element\n")

        # 节点总数
        total_nodes = (N + 1) * (N + 1)
        f.write(f"{total_nodes:<4d}   1   1   1\n")

        # 节点定义 - 使用矩阵存储节点编号
        node_id = 0
        node_matrix = np.zeros((N+1, N+1), dtype=int)
        
        for i in range(N + 1):
            for j in range(N + 1):
                node_id += 1
                node_matrix[i,j] = node_id
                x, y, z = j * 8.0 / N, i * 8.0 / N, 0.0
                # 边界条件：固定边界，自由内部
                bc = "1   1   1   1   1   1" if i == 0 or j == 0 or i == N or j == N else "1   1   0   0   0   1"
                f.write(f"{node_id:<5d} {bc}   {x:.5f}   {y:.5f}   {z:.5f}\n")

        # 内部节点加载力
        load_nodes = []
        for i in range(1, N):
            for j in range(1, N):
                load_nodes.append(node_matrix[i,j])
        force_val = force_val * 8 * 8 / N / N

        f.write(f"1   {len(load_nodes)}\n")
        for nid in load_nodes:
            f.write(f"{nid:<5d}   3   {force_val:.5f}\n")

        # 材料与单元信息
        f.write(f"3   {N*N}   1\n")
        f.write(f"1   {E:.1e}   {nu:.2f}   {thickness:.3f}   1000\n")

        # 单元定义（4节点Q4单元）
        eid = 1
        for i in range(N):
            for j in range(N):
                n1 = node_matrix[i, j]
                n2 = node_matrix[i, j + 1]
                n3 = node_matrix[i + 1, j + 1]
                n4 = node_matrix[i + 1, j]
                f.write(f"{eid:<5d}   {n1}   {n2}   {n3}   {n4}   1\n")
                eid += 1

# 主程序
if __name__ == "__main__":
    # 网格尺寸列表
    N_list = [8, 16]
    L2_errors = []
    h_list = []

    for N in N_list:
        # 生成输入文件
        input_path = f"./data/plate_{N}x{N}_Q4.dat"
        output_path = f"./data/plate_{N}x{N}_Q4.out"
        generate_plate_input(N)
        
        # 运行仿真
        print(f"Running stap++ for Q4 {N}x{N}...")
        try:
            subprocess.run(["./build/stap++.exe", input_path], check=True)
            
            # 解析输出文件，提取L2相对误差
            with open(output_path, "r") as f:
                content = f.read()
                match = re.search(r"L2 相对误差 \(w\): ([0-9.eE+-]+)", content)
                if match:
                    err = float(match.group(1))
                    L2_errors.append(err)
                    h_list.append(8.0 / N)
                    print(f"N={N}, h={8.0/N}, L2相对误差={err}")
                else:
                    print(f"警告：未找到L2误差: {output_path}")
                    L2_errors.append(np.nan)
                    h_list.append(8.0 / N)
                    
        except subprocess.CalledProcessError as e:
            print(f"运行stap++失败: {e}")
            L2_errors.append(np.nan)
            h_list.append(8.0 / N)
            continue

    # 绘制收敛图（至少有2个有效数据点时才绘制）
    if len(L2_errors) >= 2 and not all(np.isnan(L2_errors)):
        log_h = np.log(np.array(h_list))
        log_err = np.log(np.array(L2_errors))
        
        plt.figure()
        plt.plot(log_h, log_err, 'o-', label='Q4 Element')
        plt.xlabel('log(h)')
        plt.ylabel('log(L2 Error)')
        plt.title('Q4 Element L2 Error Convergence')
        
        # 计算并绘制收敛斜率
        valid_indices = ~np.isnan(log_err)
        if sum(valid_indices) >= 2:
            slope, _ = np.polyfit(log_h[valid_indices], log_err[valid_indices], 1)
            plt.plot(log_h, slope*log_h, '--', label=f'Slope={slope:.2f}')
        
        plt.legend()
        plt.grid(True)
        plt.show()
    else:
        print("有效数据点不足，无法绘制收敛图")
