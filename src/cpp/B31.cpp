/*****************************************************************************/
/*  STAP++ : B31 两结点空间线性梁单元                                          */
/*  Author : Yanzi Wang                                                      */
/*****************************************************************************/

#include "B31.h"

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

//	Constructor
CB31::CB31()
{
	NEN_ = 2;	// Each element has 2 nodes
	nodes_ = new CNode*[NEN_];
    ND_ = 12; // 6*2, 2 nodes with 6 degrees of freedom each
    LocationMatrix_ = new unsigned int[ND_];
	ElementMaterial_ = nullptr;
}

//	Desconstructor
CB31::~CB31()
{
}

//	Read element data from stream Input
bool CB31::Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList)
{
	unsigned int MSet;	// Material property set number
	unsigned int N1, N2;	// Left node number and right node number

	Input >> N1 >> N2 >> MSet;
    ElementMaterial_ = dynamic_cast<CB31Material*>(MaterialSets) + MSet - 1;
	nodes_[0] = &NodeList[N1 - 1];
	nodes_[1] = &NodeList[N2 - 1];

	return true;
}

//	Write element data to stream
void CB31::Write(COutputter& output)
{
	output << setw(11) << nodes_[0]->NodeNumber
		   << setw(9) << nodes_[1]->NodeNumber << setw(12) << ElementMaterial_->nset << endl;
}

//	Calculate element stiffness matrix 
//	Upper triangular matrix, stored as an array column by colum starting from the diagonal element
// 在 B31.cpp 中
void CB31::ElementStiffness(double* Matrix)
{
    // --- 0. 初始化和获取参数 ---
    clear(Matrix, SizeOfStiffnessMatrix()); // 清零总刚度矩阵

    double T[12][12];
    double L;
    GetTransformationMatrix(T, L); // 获取坐标变换矩阵和单元长度

    CB31Material* material_ = dynamic_cast<CB31Material*>(ElementMaterial_);
    double E = material_->E;
    double G = material_->G; 
    double A = material_->Area;
    double Iy = material_->Iy; 
    double Iz = material_->Iz;
    double J = material_->J;   
    
    // 剪切修正系数，对于矩形截面通常取 5/6
    const double kappa = 5.0 / 6.0;
    
    // 局部坐标系下的单元刚度矩阵 (最终将被填充)
    double k_local[12][12] = {0};

    // --- 1. 轴向刚度和扭转刚度 (精确积分, 相当于1点高斯积分) ---
    k_local[0][0] = E * A / L;   k_local[6][6] = E * A / L;
    k_local[0][6] = -E * A / L;  k_local[6][0] = -E * A / L;

    k_local[3][3] = G * J / L;   k_local[9][9] = G * J / L;
    k_local[3][9] = -G * J / L;  k_local[9][3] = -G * J / L;

    // --- 2. 弯曲和剪切刚度 (通过数值积分) ---
    // 形函数 N(xi) 和其导数 dN/dxi (xi是-1到1的自然坐标)
    auto N1 = [](double xi){ return 0.5 * (1.0 - xi); };
    auto N2 = [](double xi){ return 0.5 * (1.0 + xi); };
    const double dN1_dxi = -0.5;
    const double dN2_dxi = 0.5;
    
    // 物理坐标导数与自然坐标导数的关系: d/dx = (2/L) * d/dxi
    const double d_dx = 2.0 / L;

    // --- 2a. 弯曲刚度部分 (2点高斯积分) ---
    const double gauss_points_2[] = {-1.0 / sqrt(3.0), 1.0 / sqrt(3.0)};
    const double gauss_weights_2[] = {1.0, 1.0};

    for (int i = 0; i < 2; ++i) {
        double xi = gauss_points_2[i];
        double w = gauss_weights_2[i];

        // 弯曲应变-位移矩阵 B_b (bending)
        // 对应自由度: [ry1, ry2] for y-bending, [rz1, rz2] for z-bending
        double B_b_y[2] = {dN1_dxi * d_dx, dN2_dxi * d_dx}; // for y-bending (around y-axis)
        double B_b_z[2] = {dN1_dxi * d_dx, dN2_dxi * d_dx}; // for z-bending (around z-axis)

        // 组装 K_bb = integral(B_b^T * D_b * B_b) dx
        // y-bending (关联到 ry1, ry2 自由度, 即 4, 10)
        double D_b_y = E * Iy;
        k_local[4][4]   += B_b_y[0] * D_b_y * B_b_y[0] * w * (L / 2.0);
        k_local[4][10]  += B_b_y[0] * D_b_y * B_b_y[1] * w * (L / 2.0);
        k_local[10][4]  += B_b_y[1] * D_b_y * B_b_y[0] * w * (L / 2.0);
        k_local[10][10] += B_b_y[1] * D_b_y * B_b_y[1] * w * (L / 2.0);

        // z-bending (关联到 rz1, rz2 自由度, 即 5, 11)
        double D_b_z = E * Iz;
        k_local[5][5]   += B_b_z[0] * D_b_z * B_b_z[0] * w * (L / 2.0);
        k_local[5][11]  += B_b_z[0] * D_b_z * B_b_z[1] * w * (L / 2.0);
        k_local[11][5]  += B_b_z[1] * D_b_z * B_b_z[0] * w * (L / 2.0);
        k_local[11][11] += B_b_z[1] * D_b_z * B_b_z[1] * w * (L / 2.0);
    }
    
    // --- 2b. 剪切刚度部分 (1点减缩积分, xi=0, w=2) ---
    double xi = 0.0;
    double w = 2.0;
    
    // 形函数在中心点的值
    double N1_c = N1(xi); double N2_c = N2(xi);
    
    // 剪切应变-位移矩阵 B_s (shear)
    // 对应自由度: [tz1, ry1, tz2, ry2] for y-shear, [ty1, rz1, ty2, rz2] for z-shear
    
    // Y-shear (平面 x-z, 抵抗 Tz, 关联 Ry)
    // gamma_xz = duz/dx - ry
    double B_s_y[4] = {dN1_dxi * d_dx, -N1_c, dN2_dxi * d_dx, -N2_c};
    int dof_y[] = {2, 4, 8, 10}; // [tz1, ry1, tz2, ry2]
    double D_s_y = kappa * G * A;
    for(int i=0; i<4; ++i) {
        for(int j=0; j<4; ++j) {
            k_local[dof_y[i]][dof_y[j]] += B_s_y[i] * D_s_y * B_s_y[j] * w * (L / 2.0);
        }
    }

    // Z-shear (平面 x-y, 抵抗 Ty, 关联 Rz)
    // gamma_xy = duy/dx + rz
    double B_s_z[4] = {dN1_dxi * d_dx, N1_c, dN2_dxi * d_dx, N2_c};
    int dof_z[] = {1, 5, 7, 11}; // [ty1, rz1, ty2, rz2]
    double D_s_z = kappa * G * A;
    for(int i=0; i<4; ++i) {
        for(int j=0; j<4; ++j) {
            k_local[dof_z[i]][dof_z[j]] += B_s_z[i] * D_s_z * B_s_z[j] * w * (L / 2.0);
        }
    }

    // --- 3. 坐标变换: K_global = T^T * K_local * T ---
    double temp[12][12] = {0};
    double k_global[12][12] = {0};
    
    // temp = K_local * T
    for(int i=0; i<12; i++)
        for(int j=0; j<12; j++)
            for(int m=0; m<12; m++)
                temp[i][j] += k_local[i][m] * T[m][j];
    
    // k_global = T^T * temp
    for(int i=0; i<12; i++)
        for(int j=0; j<12; j++)
            for(int n=0; n<12; n++)
                k_global[i][j] += T[n][i] * temp[n][j];

    // --- 4. 填充到一维数组 Matrix (使用项目特有的逆序列优先格式) ---
    for (unsigned int j = 0; j < ND_; j++)
    {
        for (unsigned int i = 0; i <= j; i++)
        {
            unsigned int index = (j + 1) * j / 2 + (j - i);
            if (index < SizeOfStiffnessMatrix())
            {
                Matrix[index] = k_global[i][j];
            }
        }
    }
}

// 计算单元的应力
// stress[0]: 轴力N, stress[1]: 扭矩T, stress[2]: 剪力Vy, stress[3]: 剪力Vz, stress[4]: 弯矩My, stress[5]: 弯矩Mz
void CB31::ElementStress(double* stress, double* Displacement)
{
    CB31Material* material_ = dynamic_cast<CB31Material*>(ElementMaterial_);
    // 获取坐标变换矩阵T和单元长度L
    double T[12][12];
	double L;
    GetTransformationMatrix(T, L);
    // 提取单元全局自由度位移向量u_global
    double u_global[12] = {0};
    unsigned int* LocationMatrix = GetLocationMatrix();
    for(int i=0;i<12;i++) {
        int loc = LocationMatrix[i];
        if(loc) u_global[i] = Displacement[loc-1];
    }
    // 计算u_local = T^T * u_global
    double u_local[12] = {0};
    for(int i=0;i<12;i++)
        for(int j=0;j<12;j++)
            u_local[i] += T[j][i] * u_global[j];
    // 节点1、2的位移和转角（局部坐标系）
    double u1[6] = {0}, u2[6] = {0};
    for(int i=0;i<6;i++) {
        u1[i] = u_local[i];
        u2[i] = u_local[i+6];
    }
    // 轴向变形
    double du = u2[0] - u1[0];
    // 轴力
    stress[0] = material_->E * du / L;
    // 扭矩
    stress[1] = material_->G * material_->J * (u2[3] - u1[3]) / L;
    // 剪力Vy（绕z轴弯曲，y方向剪力）
    stress[2] = 12 * material_->E * material_->Iz / (L*L*L) * (u1[1] - u2[1])
                + 6 * material_->E * material_->Iz / (L*L) * (u1[5] + u2[5]);
    // 剪力Vz（绕y轴弯曲，z方向剪力）
    stress[3] = 12 * material_->E * material_->Iy / (L*L*L) * (u1[2] - u2[2])
                - 6 * material_->E * material_->Iy / (L*L) * (u1[4] + u2[4]);
    // 弯矩My（绕y轴）
    stress[4] = material_->E * material_->Iy / L * (u2[4] - u1[4]);
    // 弯矩Mz（绕z轴）
    stress[5] = material_->E * material_->Iz / L * (u2[5] - u1[5]);
}

// 构造12×12坐标变换矩阵T，并返回单元长度L
void CB31::GetTransformationMatrix(double T[12][12], double& L) const
{
    double x1 = nodes_[0]->XYZ[0], y1 = nodes_[0]->XYZ[1], z1 = nodes_[0]->XYZ[2];
    double x2 = nodes_[1]->XYZ[0], y2 = nodes_[1]->XYZ[1], z2 = nodes_[1]->XYZ[2];
    double dx = x2 - x1, dy = y2 - y1, dz = z2 - z1;
    L = sqrt(dx*dx + dy*dy + dz*dz);
    double ex[3] = {dx/L, dy/L, dz/L};
    double ref[3] = {0, 0, 1};
    if (fabs(ex[0]) < 1e-6 && fabs(ex[1]) < 1e-6) // ex接近z轴
        ref[0] = 1, ref[1] = 0, ref[2] = 0;
    double ey[3] = {
        ref[1]*ex[2] - ref[2]*ex[1],
        ref[2]*ex[0] - ref[0]*ex[2],
        ref[0]*ex[1] - ref[1]*ex[0]
    };
    double ey_len = sqrt(ey[0]*ey[0] + ey[1]*ey[1] + ey[2]*ey[2]);
    for(int i=0;i<3;i++) ey[i] /= ey_len;
    double ez[3] = {
        ex[1]*ey[2] - ex[2]*ey[1],
        ex[2]*ey[0] - ex[0]*ey[2],
        ex[0]*ey[1] - ex[1]*ey[0]
    };
    double R[3][3] = {
        {ex[0], ex[1], ex[2]},
        {ey[0], ey[1], ey[2]},
        {ez[0], ez[1], ez[2]}
    };
    for(int i=0;i<12;i++) for(int j=0;j<12;j++) T[i][j]=0;
    for(int i=0;i<4;i++)
        for(int j=0;j<3;j++)
            for(int k=0;k<3;k++)
                T[i*3+j][i*3+k] = R[j][k];
}