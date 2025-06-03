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
                x, y, z = j * 1.0, i * 1.0, 0.0
                # 边界设为 1 1 1，其余为 0 0 0
                bc = "1   1   1" if i == 0 or j == 0 or i == N or j == N else "0   0   0"
                f.write(f"{nid:<5d} {bc}   {x:.1f}   {y:.1f}   {z:.1f}\n")

        # 内部节点加载力（非边界一圈）
        load_nodes = []
        for i in range(1, N):  # i = 行
            for j in range(1, N):  # j = 列
                nid = node_id(i, j)
                load_nodes.append(nid)

        f.write(f"1   {len(load_nodes)}\n")
        for nid in load_nodes:
            f.write(f"{nid:<5d}   3   {force_val:.1f}\n")

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

if __name__ == "__main__":
    generate_plate_input(N=150)
