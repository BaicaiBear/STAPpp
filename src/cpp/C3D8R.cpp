/*
################################################
Author: Jinkun Liu
Email: liujk22@mails.tsinghua.edu.cn
################################################
*/
#include "C3D8R.h"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <Eigen/Dense>
using namespace std;
// using namespace Eigen;

//	Constructor
CC3D8R::CC3D8R()
{
	NEN_ = 8;  // Each element has 8 nodes
	nodes_ = new CNode*[NEN_];
	
	ND_ = 48; // 6*8
	LocationMatrix_ = new unsigned int[ND_];
	
	ElementMaterial_ = nullptr;
}

CC3D8R::~CC3D8R()
{
	// Destructor does not need to delete nodes_ and ElementMaterial_ as they are managed by the base class CElement
}

bool CC3D8R::Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList)
{
	unsigned int MSet;  // Material property set number
	unsigned int N1, N2, N3, N4, N5, N6, N7, N8;  // Node numbers

	Input >> N1 >> N2 >> N3 >> N4 >> N5 >> N6 >> N7 >> N8 >> MSet;
	ElementMaterial_ = dynamic_cast<CC3D8RMaterial*>(MaterialSets) + MSet - 1;

	nodes_[0] = &NodeList[N1 - 1];
	nodes_[1] = &NodeList[N2 - 1];
	nodes_[2] = &NodeList[N3 - 1];
	nodes_[3] = &NodeList[N4 - 1];
	nodes_[4] = &NodeList[N5 - 1];
	nodes_[5] = &NodeList[N6 - 1];
	nodes_[6] = &NodeList[N7 - 1];
	nodes_[7] = &NodeList[N8 - 1];

	return true;
}

void CC3D8R::Write(COutputter& output)
{
	output << setw(11) << nodes_[0]->NodeNumber
		   << setw(9) << nodes_[1]->NodeNumber
		   << setw(9) << nodes_[2]->NodeNumber
		   << setw(9) << nodes_[3]->NodeNumber
		   << setw(9) << nodes_[4]->NodeNumber
		   << setw(9) << nodes_[5]->NodeNumber
		   << setw(9) << nodes_[6]->NodeNumber
		   << setw(9) << nodes_[7]->NodeNumber
		   << setw(12) << ElementMaterial_->nset << endl;
}

//	Calculate element stiffness matrix 
//	Upper triangular matrix, stored as an array column by colum starting from the diagonal element
void CC3D8R::ElementStiffness(double* matrix)
{
    // 清零存储上三角的刚度矩阵
    clear(matrix, SizeOfStiffnessMatrix());

    // --- 1. 材料属性 ---
    CC3D8RMaterial* material_ = dynamic_cast<CC3D8RMaterial*>(ElementMaterial_);
    double E  = material_->E;    // 杨氏模量
    double nu = material_->nu;   // 泊松比
    double lambda = (E * nu) / ((1.0 + nu) * (1.0 - 2.0 * nu));  // 拉梅常数 λ
    double mu     = E / (2.0 * (1.0 + nu));                     // 剪切模量 μ

    // --- 2. 构造 6×6 弹性矩阵 Dmat ---
    double Dmat[6][6] = {0.0};
    // Dmat 前 3×3 部分先填 λ
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            Dmat[i][j] = lambda;
        }
    }
    // 再加上 2μ 到对角线 Dmat[0][0], Dmat[1][1], Dmat[2][2]
    for (int i = 0; i < 3; ++i) {
        Dmat[i][i] += 2.0 * mu;
    }
    // 剩余剪切分量 Dmat[3][3], Dmat[4][4], Dmat[5][5] = μ
    Dmat[3][3] = mu;
    Dmat[4][4] = mu;
    Dmat[5][5] = mu;

    // --- 3. 自然坐标下的形函数梯度（对 ξ, η, ζ 的偏导），中心点处都是常数 -0.125 / 0.125 ---
    // shapeFunctionGradNat[k][i]: k=0 对应 dN_i/dξ, k=1 对应 dN_i/dη, k=2 对应 dN_i/dζ
    static const double shapeFunctionGradNat[3][8] = {
        { -0.125,  0.125,  0.125, -0.125, -0.125,  0.125,  0.125, -0.125 },  // dN/dξ
        { -0.125, -0.125,  0.125,  0.125, -0.125, -0.125,  0.125,  0.125 },  // dN/dη
        { -0.125, -0.125, -0.125, -0.125,  0.125,  0.125,  0.125,  0.125 }   // dN/dζ
    };

    // --- 4. 读取物理坐标 （8 个节点） ---
    double physX[8], physY[8], physZ[8];
    for (int i = 0; i < 8; ++i) {
        physX[i] = nodes_[i]->XYZ[0];
        physY[i] = nodes_[i]->XYZ[1];
        physZ[i] = nodes_[i]->XYZ[2];
    }

    // --- 5. 计算雅可比矩阵 J（3×3）： J = shapeGradNat * physicalCoords ---
    double J[3][3] = {0.0};
    for (int row = 0; row < 3; ++row) {
        for (int col = 0; col < 3; ++col) {
            double sum = 0.0;
            for (int node = 0; node < 8; ++node) {
                double gradNat = shapeFunctionGradNat[row][node];
                if (col == 0)       sum += gradNat * physX[node];
                else if (col == 1)  sum += gradNat * physY[node];
                else                sum += gradNat * physZ[node];
            }
            J[row][col] = sum;
        }
    }

    // --- 6. 计算 det(J) ---
    double detJ =
          J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1])
        - J[0][1] * (J[1][0] * J[2][2] - J[1][2] * J[2][0])
        + J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);

    // --- 7. 计算 J 的逆矩阵 Jinv[3][3] = adj(J)/detJ ---
    double Jinv[3][3];
    // 先计算伴随矩阵 (cofactor 转置)
    Jinv[0][0] =   (J[1][1]*J[2][2] - J[1][2]*J[2][1]);
    Jinv[1][0] = - (J[1][0]*J[2][2] - J[1][2]*J[2][0]);
    Jinv[2][0] =   (J[1][0]*J[2][1] - J[1][1]*J[2][0]);

    Jinv[0][1] = - (J[0][1]*J[2][2] - J[0][2]*J[2][1]);
    Jinv[1][1] =   (J[0][0]*J[2][2] - J[0][2]*J[2][0]);
    Jinv[2][1] = - (J[0][0]*J[2][1] - J[0][1]*J[2][0]);

    Jinv[0][2] =   (J[0][1]*J[1][2] - J[0][2]*J[1][1]);
    Jinv[1][2] = - (J[0][0]*J[1][2] - J[0][2]*J[1][0]);
    Jinv[2][2] =   (J[0][0]*J[1][1] - J[0][1]*J[1][0]);

    // 除以 detJ
    double invDetJ = 1.0 / detJ;
    for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 3; ++c) {
            Jinv[r][c] *= invDetJ;
        }
    }

    // --- 8. 将形函数梯度从自然坐标映射到物理坐标： shapeGradXYZ = Jinv * shapeGradNat ---
    double shapeGradXYZ[3][8];
    for (int row = 0; row < 3; ++row) {
        for (int node = 0; node < 8; ++node) {
            double sum = 0.0;
            for (int k = 0; k < 3; ++k) {
                sum += Jinv[row][k] * shapeFunctionGradNat[k][node];
            }
            shapeGradXYZ[row][node] = sum;
        }
    }

    // --- 9. 组装 B 矩阵 (6×24) ---
    static double Bmat[6][24];
    // 先清零
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 24; ++j) {
            Bmat[i][j] = 0.0;
        }
    }
    // 对每个节点 i = 0..7，填充对应列
    for (int i = 0; i < 8; ++i) {
        int idx = 3 * i;
        double dN_dx = shapeGradXYZ[0][i];
        double dN_dy = shapeGradXYZ[1][i];
        double dN_dz = shapeGradXYZ[2][i];

        // normal strains
        Bmat[0][idx    ] = dN_dx;
        Bmat[1][idx + 1] = dN_dy;
        Bmat[2][idx + 2] = dN_dz;
        // shear strains γ_xy
        Bmat[3][idx    ] = dN_dy;
        Bmat[3][idx + 1] = dN_dx;
        // shear strains γ_yz
        Bmat[4][idx + 1] = dN_dz;
        Bmat[4][idx + 2] = dN_dy;
        // shear strains γ_zx
        Bmat[5][idx    ] = dN_dz;
        Bmat[5][idx + 2] = dN_dx;
    }

    // --- 10. 计算 Ke_local = Bᵀ * Dmat * B  （24×24）并乘以体积因子 detJ * 8 ---
    double Ke_local[24][24];
    // 首先清零
    for (int i = 0; i < 24; ++i) {
        for (int j = 0; j < 24; ++j) {
            Ke_local[i][j] = 0.0;
        }
    }
    // B·D·B 形式：Ke[i][j] = sum_{m=0..5, n=0..5} B[m][i] * Dmat[m][n] * B[n][j]
    for (int i = 0; i < 24; ++i) {
        for (int j = 0; j < 24; ++j) {
            double sum_mn = 0.0;
            for (int m = 0; m < 6; ++m) {
                double B_m_i = Bmat[m][i];
                if (B_m_i == 0.0) continue;
                for (int n = 0; n < 6; ++n) {
                    double B_n_j = Bmat[n][j];
                    if (B_n_j == 0.0) continue;
                    sum_mn += B_m_i * Dmat[m][n] * B_n_j;
                }
            }
            Ke_local[i][j] = sum_mn * (detJ * 8.0);
        }
    }

    // --- 11. Hourglass 控制项 ---
    // 11.1 定义 γ 模式矩阵（4×8）
    static const double gammaHG[4][8] = {
        { -1.0,  1.0,  1.0, -1.0, -1.0,  1.0,  1.0, -1.0 },  // γ1
        { -1.0, -1.0,  1.0,  1.0, -1.0, -1.0,  1.0,  1.0 },  // γ2
        { -1.0, -1.0, -1.0, -1.0,  1.0,  1.0,  1.0,  1.0 },  // γ3
        {  1.0, -1.0,  1.0, -1.0, -1.0,  1.0, -1.0,  1.0 }   // γ4
    };

    // 11.2 先计算 Jinv 每行的元素之和（用于 hourglass 映射）
    double sumInvJRow[3] = {
        Jinv[0][0] + Jinv[0][1] + Jinv[0][2],
        Jinv[1][0] + Jinv[1][1] + Jinv[1][2],
        Jinv[2][0] + Jinv[2][1] + Jinv[2][2]
    };

    // 11.3 Khg 初始化为 0（24×24）
    double Khg[24][24];
    for (int i = 0; i < 24; ++i) {
        for (int j = 0; j < 24; ++j) {
            Khg[i][j] = 0.0;
        }
    }

    // 11.4 体积 V = |detJ * 8|
    double V = fabs(detJ * 8.0);
    // l2 = (cuberoot(V))^2
    double cbrtV = pow(V, 1.0 / 3.0);
    double l2    = cbrtV * cbrtV;
    // hourglass 强度系数
    const double hgAlpha = 0.005;
    double gh = mu * V / l2 * hgAlpha;

    // 11.5 对每一个 hourglass 模式做累加
    for (int hgIdx = 0; hgIdx < 4; ++hgIdx) {
        // 11.5.1 先构造节点方向上的 hourglass 梯度，物理坐标下
        double B_hg_phys[3][8];
        for (int node = 0; node < 8; ++node) {
            double gammaVal = gammaHG[hgIdx][node];
            // B_hg_phys[row][node] = sumInvJRow[row] * gammaVal
            B_hg_phys[0][node] = sumInvJRow[0] * gammaVal;
            B_hg_phys[1][node] = sumInvJRow[1] * gammaVal;
            B_hg_phys[2][node] = sumInvJRow[2] * gammaVal;
        }

        // 11.5.2 转换成 24×1 向量 B_hg_vec
        double B_hg_vec[24];
        for (int i = 0; i < 8; ++i) {
            int idx = 3 * i;
            B_hg_vec[idx    ] = B_hg_phys[0][i];
            B_hg_vec[idx + 1] = B_hg_phys[1][i];
            B_hg_vec[idx + 2] = B_hg_phys[2][i];
        }

        // 11.5.3 用外积累加： Khg += gh * (B_hg_vec * B_hg_vecᵀ)
        for (int i = 0; i < 24; ++i) {
            double Bi = B_hg_vec[i];
            if (Bi == 0.0) continue;
            for (int j = 0; j < 24; ++j) {
                double Bj = B_hg_vec[j];
                if (Bj == 0.0) continue;
                Khg[i][j] += gh * (Bi * Bj);
            }
        }
    }

    // 11.6 将 Khg 累加到 Ke_local
    for (int i = 0; i < 24; ++i) {
        for (int j = 0; j < 24; ++j) {
            Ke_local[i][j] += Khg[i][j];
        }
    }

    // --- 12. 将 Ke_local 上三角 / 对称存储 到一维数组 matrix（按列，从对角线开始） ---
    // 存储顺序：对每一列 j，从 i=j, j-1, …, 0 一直到 i=0
    int idxPack = 0;
    for (int m = 0; m < 8; ++m) {
        for (int n = 0; n < 3; ++n) {
            int j = m * 3 + n, realj = 6 * m + n;
            for (int k = m; k >= 0; --k) {
                for (int l = 3; l >= 0; --l) {
                    int i = k * 3 + l, reali = 6 * k + l;
                    if (i > j) continue;  
                matrix[realj*(realj+1)/2+realj-reali] = Ke_local[i][j];
                }
            }
        }
    }
}

void CC3D8R::ElementStress(double* stress, double* Displacement)
{
    // --- 1. 获取材料属性 ---
    CC3D8RMaterial* material_ = dynamic_cast<CC3D8RMaterial*>(ElementMaterial_);
    double E  = material_->E;    // 杨氏模量
    double nu = material_->nu;   // 泊松比

    // --- 2. 构造 6×6 弹性矩阵 Dmat ---
    //    Dmat = | λ+2μ   λ       λ      0    0    0 |
    //           |   λ   λ+2μ     λ      0    0    0 |
    //           |   λ     λ   λ+2μ    0    0    0 |
    //           |   0     0     0     μ    0    0 |
    //           |   0     0     0     0    μ    0 |
    //           |   0     0     0     0    0    μ |
    double lambda = (E * nu) / ((1.0 + nu) * (1.0 - 2.0 * nu));
    double mu     = E / (2.0 * (1.0 + nu));
    double Dmat[6][6] = {0.0};
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            Dmat[i][j] = lambda;
        }
        Dmat[i][i] += 2.0 * mu;
    }
    Dmat[3][3] = mu;
    Dmat[4][4] = mu;
    Dmat[5][5] = mu;

    // --- 3. 收集 8 个节点的物理坐标 ---
    double physX[8], physY[8], physZ[8];
    for (int i = 0; i < 8; ++i) {
        physX[i] = nodes_[i]->XYZ[0];
        physY[i] = nodes_[i]->XYZ[1];
        physZ[i] = nodes_[i]->XYZ[2];
    }

    // --- 4. 自然坐标下形函数梯度 shapeFunctionGradNat[3][8] ---
    static const double shapeFunctionGradNat[3][8] = {
        { -0.125,  0.125,  0.125, -0.125, -0.125,  0.125,  0.125, -0.125 },  // dN/dξ
        { -0.125, -0.125,  0.125,  0.125, -0.125, -0.125,  0.125,  0.125 },  // dN/dη
        { -0.125, -0.125, -0.125, -0.125,  0.125,  0.125,  0.125,  0.125 }   // dN/dζ
    };

    // --- 5. 计算雅可比矩阵 J[3][3] ---
    double J[3][3] = {0.0};
    for (int row = 0; row < 3; ++row) {
        for (int col = 0; col < 3; ++col) {
            double sum = 0.0;
            for (int node = 0; node < 8; ++node) {
                double gNat = shapeFunctionGradNat[row][node];
                if (col == 0)       sum += gNat * physX[node];
                else if (col == 1)  sum += gNat * physY[node];
                else                sum += gNat * physZ[node];
            }
            J[row][col] = sum;
        }
    }

    // --- 6. 计算 det(J) ---
    double detJ =
          J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1])
        - J[0][1] * (J[1][0] * J[2][2] - J[1][2] * J[2][0])
        + J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);

    // --- 7. 计算 J 的伴随矩阵（cofactor 转置），然后除以 detJ 得到逆矩阵 Jinv ---
    double Jinv[3][3];
    // 伴随矩阵（代数余子式的转置）
    Jinv[0][0] =   (J[1][1]*J[2][2] - J[1][2]*J[2][1]);
    Jinv[1][0] = - (J[1][0]*J[2][2] - J[1][2]*J[2][0]);
    Jinv[2][0] =   (J[1][0]*J[2][1] - J[1][1]*J[2][0]);

    Jinv[0][1] = - (J[0][1]*J[2][2] - J[0][2]*J[2][1]);
    Jinv[1][1] =   (J[0][0]*J[2][2] - J[0][2]*J[2][0]);
    Jinv[2][1] = - (J[0][0]*J[2][1] - J[0][1]*J[2][0]);

    Jinv[0][2] =   (J[0][1]*J[1][2] - J[0][2]*J[1][1]);
    Jinv[1][2] = - (J[0][0]*J[1][2] - J[0][2]*J[1][0]);
    Jinv[2][2] =   (J[0][0]*J[1][1] - J[0][1]*J[1][0]);

    double invDetJ = 1.0 / detJ;
    for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 3; ++c) {
            Jinv[r][c] *= invDetJ;
        }
    }

    // --- 8. 将形函数梯度映射到物理坐标下 shapeGradXYZ[3][8] ---
    double shapeGradXYZ[3][8];
    for (int row = 0; row < 3; ++row) {
        for (int node = 0; node < 8; ++node) {
            double sum = 0.0;
            for (int k = 0; k < 3; ++k) {
                sum += Jinv[row][k] * shapeFunctionGradNat[k][node];
            }
            shapeGradXYZ[row][node] = sum;
        }
    }

    // --- 9. 从全局位移向量构造单元位移向量 u_e[24] ---
    double u_e[24];
    for (int i = 0; i < 8; ++i) {
        for (int j = 0; j < 3; ++j) {
            if (LocationMatrix_[6*i+j] > 0) {
                u_e[3*i+j] = Displacement[LocationMatrix_[6*i+j] - 1];
            } else {
                u_e[6*i+j] = 0.0;
            }
        }
    }

    // --- 10. 计算应变向量 eps[6] = B * u_e，直接展开 B·u_e ---
    double eps[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    for (int node = 0; node < 8; ++node) {
        int idx = 3 * node;
        double ux = u_e[idx];
        double uy = u_e[idx + 1];
        double uz = u_e[idx + 2];

        double dN_dx = shapeGradXYZ[0][node];
        double dN_dy = shapeGradXYZ[1][node];
        double dN_dz = shapeGradXYZ[2][node];

        // ε_xx
        eps[0] += dN_dx * ux;
        // ε_yy
        eps[1] += dN_dy * uy;
        // ε_zz
        eps[2] += dN_dz * uz;
        // γ_xy = ε_xy + ε_yx
        eps[3] += dN_dy * ux + dN_dx * uy;
        // γ_yz = ε_yz + ε_zy
        eps[4] += dN_dz * uy + dN_dy * uz;
        // γ_zx = ε_zx + ε_xz
        eps[5] += dN_dz * ux + dN_dx * uz;
    }

    // --- 11. 根据 σ = Dmat · eps 计算应力 σ[6] ---
    double sigma[6];
    for (int i = 0; i < 6; ++i) {
        double sum = 0.0;
        for (int j = 0; j < 6; ++j) {
            sum += Dmat[i][j] * eps[j];
        }
        sigma[i] = sum;
    }

    // --- 12. 将结果存入输出数组 stress[6] ---
    for (int i = 0; i < 6; ++i) {
        stress[i] = sigma[i];
    }
}
