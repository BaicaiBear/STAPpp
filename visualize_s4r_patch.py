import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import re

# 读取节点、单元、位移、应力
node_coords = {}
elements = []
displacements = {}
stresses = {}

with open('plate-example_8x8.dat', 'r') as f:
    lines = f.readlines()
    # 读取节点
    node_start = 2
    node_count = int(lines[1].split()[0])
    for i in range(node_count):
        parts = lines[node_start + i].split()
        node_coords[int(parts[0])] = np.array([float(parts[-3]), float(parts[-2]), float(parts[-1])])
    # 读取所有单元定义（每行6列且首列为单元号，后4列为节点号，第6列为材料号）
    for line in lines:
        parts = line.split()
        if len(parts) == 6 and parts[0].isdigit():
            n1, n2, n3, n4 = map(int, parts[1:5])
            if all(n in node_coords for n in [n1, n2, n3, n4]):
                elements.append([n1, n2, n3, n4])

# 读取out文件中的位移和应力
with open('plate-example_8x8.out', 'r') as f:
    lines = f.readlines()
    # 位移
    disp_start = None
    for idx, line in enumerate(lines):
        if ' D I S P L A C E M E N T S (S4R: theta_x, theta_y, w)' in line:
            disp_start = idx + 3
            break
    if disp_start:
        for i in range(len(node_coords)):
            if disp_start + i >= len(lines):
                break
            print(lines[disp_start + i])
            parts = lines[disp_start + i].split()
            if len(parts) < 4:
                continue  # 跳过空行或格式不对的行
            node = int(parts[0])
            disp = np.array([0, 0, float(parts[3])])
            displacements[node] = disp
            print(f"Node {node}: Displacement = {disp}")

# 变形放大系数
def_scale = 1e3  # 10e5 = 1e6

fig = plt.figure(figsize=(10,5))
ax1 = fig.add_subplot(121, projection='3d')
ax2 = fig.add_subplot(122, projection='3d')

print(f"节点数: {len(node_coords)}")
print(f"单元数: {len(elements)}")
print(f"节点编号: {list(node_coords.keys())}")
print(f"单元节点: {elements}")
print(f"位移节点: {list(displacements.keys())}")

for idx, elem in enumerate(elements):
    verts = [node_coords[n] for n in elem]
    elem_stress = stresses.get(idx+1, np.zeros(3))
    poly = Poly3DCollection([verts], alpha=0.6, edgecolor='k')
    ax1.add_collection3d(poly)
    for v in verts:
        ax1.scatter(*v, color='k')
ax1.set_title('Original Mesh')
ax1.set_xlabel('X')
ax1.set_ylabel('Y')
ax1.set_zlabel('Z')

# 变形后网格及应力云图
for idx, elem in enumerate(elements):
    verts = []
    for n in elem:
        if n not in displacements:
            disp = np.zeros(3)
        else:
            disp = displacements[n]
        # 变形后坐标 = 原坐标 + 放大因子 * 变形量
        verts.append(node_coords[n] + def_scale * disp)
    elem_stress = stresses.get(idx+1, np.zeros(3))
    poly = Poly3DCollection([verts], alpha=0.6, edgecolor='k')
    ax2.add_collection3d(poly)
    for v in verts:
        ax2.scatter(*v, color='k')
ax2.set_title('Deformed Mesh')
ax2.set_xlabel('X')
ax2.set_ylabel('Y')
ax2.set_zlabel('Z')


plt.tight_layout()
plt.show()
