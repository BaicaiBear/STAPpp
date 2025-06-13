import numpy as np
import matplotlib.pyplot as plt

# 读取节点坐标
node_coords = {}
displacements = {}

with open('data/plate_8x8.dat', 'r') as f:
    lines = f.readlines()
    node_start = 2
    node_count = int(lines[1].split()[0])
    for i in range(node_count):
        parts = lines[node_start + i].split()
        node_coords[int(parts[0])] = np.array([float(parts[-3]), float(parts[-2]), float(parts[-1])])

# 读取节点位移
with open('data/plate_8x8.out', 'r') as f:
    lines = f.readlines()
    disp_start = None
    for idx, line in enumerate(lines):
        if ' D I S P L A C E M E N T S' in line:
            disp_start = idx + 3
            break
    if disp_start:
        for i in range(len(node_coords)):
            if disp_start + i >= len(lines):
                break
            parts = lines[disp_start + i].split()
            if len(parts) < 4:
                continue
            node = int(parts[0])
            disp = np.array([float(parts[1]), float(parts[2]), float(parts[3])])
            displacements[node] = disp

# 变形放大系数
def_scale = 1e8  # 可根据实际情况调整

fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection='3d')

# 原始节点
orig_xyz = np.array([node_coords[n] for n in node_coords])
ax.scatter(orig_xyz[:,0], orig_xyz[:,1], orig_xyz[:,2], c='b', label='Original', s=20)

# 变形后节点
deformed_xyz = []
for n in node_coords:
    disp = displacements.get(n, np.zeros(3))
    deformed_xyz.append(node_coords[n] + def_scale * disp)
deformed_xyz = np.array(deformed_xyz)
ax.scatter(deformed_xyz[:,0], deformed_xyz[:,1], deformed_xyz[:,2], c='r', label='Deformed', s=20)

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('Original (blue) & Deformed (red) Nodes')
ax.legend()
plt.tight_layout()
plt.show()