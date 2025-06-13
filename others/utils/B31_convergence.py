import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt

# ==============================================================================
# 0. 定义路径
# ==============================================================================
PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__))
EXECUTABLE_PATH = os.path.join(PROJECT_ROOT, 'build', 'debug', 'stap++')
DATA_DIR = os.path.join(PROJECT_ROOT, 'data')
os.makedirs(DATA_DIR, exist_ok=True)
if not os.path.exists(EXECUTABLE_PATH):
    if os.path.exists(EXECUTABLE_PATH + '.exe'): EXECUTABLE_PATH += '.exe'
    else:
        print(f"*** FATAL ERROR: Executable not found at '{EXECUTABLE_PATH}'"); exit(1)

# ==============================================================================
# 1. 定义测试问题参数
# ==============================================================================
L = 10.0; H = 0.2; B = 0.1; E = 2.1E11; NU = 0.3; P = -1000.0
G = E / (2 * (1 + NU)); Area = B * H; Iy = (B * H**3) / 12.0
Iz = (H * B**3) / 12.0; J_approx = Iy + Iz

# ==============================================================================
# 2. 定义辅助函数 
# ==============================================================================

def create_input_file(N):
    base_name = f"beam_N{N}"
    file_path = os.path.join(DATA_DIR, f"{base_name}.dat")
    nodes = []
    nodes.append("     1   1   1   1   1   1   1       0.0       0.0       0.0\n")
    for i in range(1, N):
        x_coord = (i / N) * L
        nodes.append(f"     {i+1}   0   0   0   0   0   0    {x_coord:.6f}       0.0       0.0\n")
    # 约束末端扭转自由度 Rx (第4个)
    nodes.append(f"     {N+1}   0   0   0   0   0   0    {L:.6f}       0.0       0.0\n")
    elements = []
    for i in range(N):
        elements.append(f"     {i+1}     {i+1}      {i+2}       1\n")
    file_content = f"""Cantilever Beam Convergence Test (N={N})
     {N+1}      1      1      1
{''.join(nodes)}
     1       1
     {N+1}      3   {P:.1f}
     4       {N}       1
     1   {E:.4E}  {G:.4E}  {Area:.4E}  {Iy:.4E}  {Iz:.4E}  {J_approx:.4E}
{''.join(elements)}
"""
    with open(file_path, 'w', encoding='utf-8') as f: f.write(file_content)
    return base_name

def run_stappp(base_name):
    try:
        process = subprocess.run([EXECUTABLE_PATH, base_name], check=True, capture_output=True, text=True, cwd=DATA_DIR, encoding='utf-8', errors='replace')
        return os.path.join(DATA_DIR, f"{base_name}.out")
    except subprocess.CalledProcessError as e:
        print("--- STAP++ Execution Failed ---"); print("STDERR:", e.stderr); return None


def parse_output_file(output_filename):
    """
    从.out文件中解析节点位移，同时获取 uz (z位移) 和 ry (绕y轴转角)。
    返回一个字典，键是节点ID，值是包含 (uz, ry) 的元组。
    """
    nodal_data = {}
    if not os.path.exists(output_filename): return None
    try:
        with open(output_filename, 'r', encoding='utf-8', errors='replace') as f:
            lines = f.readlines()
    except Exception as e:
        print(f"Error reading {output_filename}: {e}"); return None

    start_index = -1
    for i, line in enumerate(lines):
        if "D I S P L A C E M E N T S" in line:
            start_index = i + 3; break
    if start_index == -1: return None

    for line in lines[start_index:]:
        if not line.strip(): break
        parts = line.split()
        if len(parts) >= 7:
            try:
                node_id = int(parts[0])
                uz = float(parts[3])  # Z-displacement (index 3)
                ry = float(parts[5])  # Theta-Y (index 5)
                nodal_data[node_id] = {'uz': uz, 'ry': ry}
            except (ValueError, IndexError): continue
            
    return nodal_data


def hermite_interpolation(xi, v1, theta1, v2, theta2, h):
    """
    使用三次Hermite形函数在单元内插值。
    xi: 单元内自然坐标 (-1 to 1)
    v1, v2: 单元两端节点的位移
    theta1, theta2: 单元两端节点的转角
    h: 单元长度
    """
    H1 = 0.25 * (1 - xi)**2 * (2 + xi)
    H2 = 0.25 * (1 - xi)**2 * (1 + xi)
    H3 = 0.25 * (1 + xi)**2 * (2 - xi)
    H4 = 0.25 * (1 + xi)**2 * (xi - 1)
    return H1 * v1 + H2 * (h / 2.0) * theta1 + H3 * v2 + H4 * (h / 2.0) * theta2


def calculate_l2_error(N, nodal_data):
    """
    使用高斯求积精确计算L2范数误差。
    """
    total_squared_error = 0.0
    
    # 2点高斯积分点和权重
    gauss_points = [-1/np.sqrt(3), 1/np.sqrt(3)]
    gauss_weights = [1.0, 1.0]
    
    element_length = L / N

    for i in range(N):  # 遍历每个单元
        node1_id = i + 1
        node2_id = i + 2
        
        # 获取单元两端的节点自由度值
        uz1 = nodal_data[node1_id]['uz']
        ry1 = nodal_data[node1_id]['ry']
        uz2 = nodal_data[node2_id]['uz']
        ry2 = nodal_data[node2_id]['ry']
        
        # 在单元内进行高斯求积
        for j, xi in enumerate(gauss_points):
            # 物理坐标
            x_physical = (i + 0.5 * (1 + xi)) * element_length
            
            # 计算该点的数值解和精确解
            uz_fem = hermite_interpolation(xi, uz1, ry1, uz2, ry2, element_length)
            uz_exact = analytical_solution(x_physical)
            
            # 累加误差的平方
            squared_error = (uz_fem - uz_exact)**2
            
            # 乘以权重和雅可比行列式 (对于一维问题，雅可比是 h/2)
            total_squared_error += squared_error * gauss_weights[j] * (element_length / 2.0)
            
    return np.sqrt(total_squared_error)

def analytical_solution(x):
    """
    悬臂梁在末端受集中力P作用下的铁木辛柯梁精确解。
    """
    # 弯曲变形部分 (来自欧拉-伯努利理论)
    uz_bending = (P * x**2) / (6 * E * Iy) * (3 * L - x)
    
    # 剪切变形部分
    kappa = 5.0 / 6.0 # 剪切修正系数
    uz_shear = (P * x) / (kappa * G * Area)

    return uz_bending + uz_shear

# ==============================================================================
# 3. 主执行流程
# ==============================================================================

if __name__ == "__main__":
    
    element_counts = [2, 4, 8, 16, 32, 64, 128] # 增加更多数据点
    h_values = []
    error_values = []
    
    print("--- Starting Convergence Study for B31 Element (Type 4) ---")
    
    for N in element_counts:
        print(f"\nProcessing case with N = {N} elements...")
        base_name = create_input_file(N)
        output_file = run_stappp(base_name)
        if not output_file: break
        
        nodal_data = parse_output_file(output_file)
        if not nodal_data: print(f"  Failed to parse data for N={N}."); break
        print(f"  Parsed {len(nodal_data)} nodal data points.")

        h = L / N
        l2_error = calculate_l2_error(N, nodal_data)
        print(f"  h = {h:.4f}, L2 Error = {l2_error:.6E}")

        h_values.append(h)
        error_values.append(l2_error)
        
    print("\n--- Convergence Study Finished ---")

    # 绘制双对数曲线
    if h_values and error_values:
        h_values = np.array(h_values)
        error_values = np.array(error_values)
        
        # 在对数空间中进行线性拟合，计算收敛率 p
        # 我们只用后半部分的数据点来拟合，因为粗网格的结果可能不够稳定
        fit_range = slice(len(h_values) // 2, None)
        log_h = np.log(h_values[fit_range])
        log_error = np.log(error_values[fit_range])
        slope, intercept = np.polyfit(log_h, log_error, 1)
        
        plt.style.use('seaborn-v0_8-whitegrid')
        plt.figure(figsize=(10, 8))
        
        # 绘制FEM误差数据点
        plt.loglog(h_values, error_values, 'o', markersize=8, markerfacecolor='none', markeredgecolor='blue', label='FEM L2 Error')
        
        # 绘制拟合线
        fit_line = np.exp(intercept) * (h_values**slope)
        plt.loglog(h_values, fit_line, '--', color='red', label=f'Best Fit (slope p = {slope:.4f})')
        
        # 绘制理论收敛率参考线
        # O(h^2) 参考线
        h2_line = error_values[0] * (h_values / h_values[0])**2
        plt.loglog(h_values, h2_line, ':', color='gray', label='Reference slope p = 2')
        
        
        plt.title('Convergence Analysis of B31 Beam Element', fontsize=16)
        plt.xlabel('Element Size, h', fontsize=12)
        plt.ylabel('L2 Norm of Displacement Error', fontsize=12)
        plt.legend(fontsize=11)
        plt.gca().invert_xaxis()
        
        plot_path = os.path.join(PROJECT_ROOT, "convergence_plot_B31.png")
        plt.savefig(plot_path, dpi=300)
        print(f"\nConvergence plot saved to '{plot_path}'")
        plt.show()