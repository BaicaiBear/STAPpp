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
    CC3D8RMaterial* material = dynamic_cast<CC3D8RMaterial*>(ElementMaterial_);
    double E = material->E;
    double nu = material->nu;
    double lambda = (E * nu) / ((1.0 + nu) * (1.0 - 2.0 * nu));
    double mu     = E / (2.0 * (1.0 + nu));

    // --- 2. 构造 D 矩阵 ---
    double Dmat[6][6] = {0.0};
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j)
            Dmat[i][j] = lambda;
        Dmat[i][i] += 2.0 * mu;
    }
    Dmat[3][3] = mu;
    Dmat[4][4] = mu;
    Dmat[5][5] = mu;

    // --- 3. Gauss 点设置 ---
    const int ngp = 2;
    double gpLoc[2] = {-1.0 / sqrt(3.0), 1.0 / sqrt(3.0)};
    double gpW[2]   = {1.0, 1.0};

    double Ke_local[24][24] = {0.0};
    static const int sign[8][3] = {
        {-1, -1, -1}, {1, -1, -1}, {1, 1, -1}, {-1, 1, -1},
        {-1, -1,  1}, {1, -1,  1}, {1, 1,  1}, {-1, 1,  1}
    };

    // --- 4. 遍历 Gauss 点 ---
    for (int i = 0; i < ngp; ++i) {
    for (int j = 0; j < ngp; ++j) {
    for (int k = 0; k < ngp; ++k) {
        double xi   = gpLoc[i];
        double eta  = gpLoc[j];
        double zeta = gpLoc[k];
        double w    = gpW[i] * gpW[j] * gpW[k];

        // --- 4.1. 形函数在自然坐标下的梯度 ---
        double shapeGradNat[3][8];
        for (int n = 0; n < 8; ++n) {
            int xi_i   = sign[n][0];
            int eta_i  = sign[n][1];
            int zeta_i = sign[n][2];
            double c0 = 0.125 * (1 + eta * eta_i) * (1 + zeta * zeta_i);
            double c1 = 0.125 * (1 + xi  * xi_i)  * (1 + zeta * zeta_i);
            double c2 = 0.125 * (1 + xi  * xi_i)  * (1 + eta  * eta_i);
            shapeGradNat[0][n] = xi_i   * c0;
            shapeGradNat[1][n] = eta_i  * c1;
            shapeGradNat[2][n] = zeta_i * c2;
        }

        // --- 4.2. 物理坐标 ---
        double physX[8], physY[8], physZ[8];
        for (int n = 0; n < 8; ++n) {
            physX[n] = nodes_[n]->XYZ[0];
            physY[n] = nodes_[n]->XYZ[1];
            physZ[n] = nodes_[n]->XYZ[2];
        }

        // --- 4.3. 计算雅可比 J ---
        double J[3][3] = {0.0};
        for (int r = 0; r < 3; ++r) {
            for (int c = 0; c < 3; ++c) {
                double sum = 0.0;
                for (int n = 0; n < 8; ++n) {
                    double g = shapeGradNat[r][n];
                    sum += g * (c == 0 ? physX[n] : c == 1 ? physY[n] : physZ[n]);
                }
                J[r][c] = sum;
            }
        }
        double detJ =
            J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1])
          - J[0][1] * (J[1][0] * J[2][2] - J[1][2] * J[2][0])
          + J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);
        double invDetJ = 1.0 / detJ;

        // --- 4.4. 计算 Jinv ---
        double Jinv[3][3];
        Jinv[0][0] =  (J[1][1]*J[2][2] - J[1][2]*J[2][1]);
        Jinv[1][0] = -(J[1][0]*J[2][2] - J[1][2]*J[2][0]);
        Jinv[2][0] =  (J[1][0]*J[2][1] - J[1][1]*J[2][0]);
        Jinv[0][1] = -(J[0][1]*J[2][2] - J[0][2]*J[2][1]);
        Jinv[1][1] =  (J[0][0]*J[2][2] - J[0][2]*J[2][0]);
        Jinv[2][1] = -(J[0][0]*J[2][1] - J[0][1]*J[2][0]);
        Jinv[0][2] =  (J[0][1]*J[1][2] - J[0][2]*J[1][1]);
        Jinv[1][2] = -(J[0][0]*J[1][2] - J[0][2]*J[1][0]);
        Jinv[2][2] =  (J[0][0]*J[1][1] - J[0][1]*J[1][0]);
        for (int r = 0; r < 3; ++r)
            for (int c = 0; c < 3; ++c)
                Jinv[r][c] *= invDetJ;

        // --- 4.5. 物理坐标下的形函数梯度 ---
        double shapeGradXYZ[3][8];
        for (int r = 0; r < 3; ++r) {
            for (int n = 0; n < 8; ++n) {
                double sum = 0.0;
                for (int s = 0; s < 3; ++s)
                    sum += Jinv[r][s] * shapeGradNat[s][n];
                shapeGradXYZ[r][n] = sum;
            }
        }

        // --- 4.6. 组装 B 矩阵 ---
        double B[6][24] = {0.0};
        for (int n = 0; n < 8; ++n) {
            int idx = 3 * n;
            double dNdx = shapeGradXYZ[0][n];
            double dNdy = shapeGradXYZ[1][n];
            double dNdz = shapeGradXYZ[2][n];
            B[0][idx    ] = dNdx;
            B[1][idx + 1] = dNdy;
            B[2][idx + 2] = dNdz;
            B[3][idx    ] = dNdy;
            B[3][idx + 1] = dNdx;
            B[4][idx + 1] = dNdz;
            B[4][idx + 2] = dNdy;
            B[5][idx    ] = dNdz;
            B[5][idx + 2] = dNdx;
        }

        // --- 4.7. 累加 Ke_local ---
        for (int a = 0; a < 24; ++a) {
            for (int b = 0; b < 24; ++b) {
                double sum = 0.0;
                for (int m2 = 0; m2 < 6; ++m2) {
                    for (int n2 = 0; n2 < 6; ++n2) {
                        sum += B[m2][a] * Dmat[m2][n2] * B[n2][b];
                    }
                }
                Ke_local[a][b] += sum * detJ * w;
            }
        }
    }}} // end Gauss loops

    // --- 12. 将 Ke_local 上三角 / 对称存储 到一维数组 matrix（按列，从对角线开始） ---
    int idxPack = 0;
    for (int m = 0; m < 8; ++m) {
        for (int n = 0; n < 3; ++n) {
            int j = m * 3 + n, realj = 6 * m + n;
            for (int k = m; k >= 0; --k) {
                for (int l = 3; l >= 0; --l) {
                    int i = k * 3 + l, reali = 6 * k + l;
                    if (i > j) continue;
                    matrix[realj * (realj + 1) / 2 + realj - reali] = Ke_local[i][j];
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
