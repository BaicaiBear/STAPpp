def generate_one_edge_fixed_plate():
    """
    生成一个4节点、4单元的非标准正方形薄板输入文件，
    左边一排节点全部固定（w, θx, θy=0），其余自由，无外载荷。
    单元为非标准正方形（可人为扰动部分节点坐标）。
    """
    # 节点坐标（人为扰动部分节点）
    nodes = [
        (1, 0.0, 0.0, 0.0),
        (2, 1.0, 0.0, 0.0),
        (3, 2.33, 0.0, 0.0),
        (4, 0.0, 1.0, 0.0),
        (5, 0.84, 0.7, 0.0),
        (6, 2.0, 1.1, 0.0),
        (7, 0.1, 2.11, 0.0),
        (8, 1.0, 2.0, 0.0),
        (9, 1.9, 2.2, 0.0),
    ]
    # 单元定义（4个非标准四边形）
    elements = [
        (1, 1, 2, 5, 4),
        (2, 2, 3, 6, 5),
        (3, 4, 5, 8, 7),
        (4, 5, 6, 9, 8),
    ]
    # 左边一排节点固定
    fixed_nodes = {1, 4}
    # 文件写入
    with open("./data/plate_1edge_fixed_4elem.dat", "w") as f:
        f.write("Plate Element Test - 4 Elements (1 edge fixed, 4 quad)\n")
        f.write("9   1   1   1\n")
        for nid, x, y, z in nodes:
            if nid in fixed_nodes:
                bc = "1 0.0   1 0.0   1 0.0   1 0.0   1 0.1   1 0.0"
            else:
                bc = "1 0.0   1 0.0   0 0.0   0 0.0   0 0.0   1 0.0"
            f.write(f"{nid:<5d} {bc}   {x:.3f}   {y:.3f}   {z:.3f}\n")
        # 无外载荷
        f.write("1   0\n")
        # 材料与单元信息
        f.write("3   4   1\n")
        f.write("1   2.0e+11   0.30   0.010   1000\n")
        for eid, n1, n2, n3, n4 in elements:
            f.write(f"{eid:<5d}   {n1}   {n2}   {n3}   {n4}   1\n")

if __name__ == "__main__":
    generate_one_edge_fixed_plate()