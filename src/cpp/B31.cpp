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
    ElementMaterial_ = dynamic_cast<CBarMaterial*>(MaterialSets) + MSet - 1;
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
void CB31::ElementStiffness(double* Matrix)
{
	clear(Matrix, SizeOfStiffnessMatrix());
    
    // 获取坐标变换矩阵T和单元长度L
    double T[12][12];
    double L;
    GetTransformationMatrix(T, L);

    // 材料参数
    CB31Material* material_ = dynamic_cast<CB31Material*>(ElementMaterial_);
    double E = material_->E;
    double G = material_->G; 
    double A = material_->Area;
    double Iy = material_->Iy; 
    double Iz = material_->Iz;
    double J = material_->J;   

    // 局部刚度矩阵（12×12）
    double k[12][12] = {0};
    double EA_L = E*A/L;
    double GJ_L = G*J/L;
    double EIy_L = E*Iy/L;
    double EIz_L = E*Iz/L;
    double EIy_L3 = 2*E*Iy/L;
    double EIz_L3 = 2*E*Iz/L;
    double EIy_L2 = 6*E*Iy/(L*L);
    double EIz_L2 = 6*E*Iz/(L*L);
    double EIy_L4 = 4*E*Iy/L;
    double EIz_L4 = 4*E*Iz/L;
    // 轴向
    k[0][0] = k[6][6] = EA_L;
    k[0][6] = k[6][0] = -EA_L;
    // 扭转
    k[3][3] = k[9][9] = GJ_L;
    k[3][9] = k[9][3] = -GJ_L;
    // y方向弯曲
    k[1][1] = k[7][7] = EIz_L3;
    k[1][7] = k[7][1] = EIz_L2;
    k[1][5] = k[5][1] =  EIz_L2;
    k[1][11]= k[11][1]= -EIz_L2;
    k[5][5] = k[11][11]= EIz_L4;
    k[5][7] = k[7][5] = -EIz_L2;
    k[5][11]= k[11][5]= 2*EIz_L4;
    k[7][11]= k[11][7]= -EIz_L2;
    // z方向弯曲
    k[2][2] = k[8][8] = EIy_L3;
    k[2][8] = k[8][2] = EIy_L2;
    k[2][4] = k[4][2] = -EIy_L2;
    k[2][10]= k[10][2]=  EIy_L2;
    k[4][4] = k[10][10]= EIy_L4;
    k[4][8] = k[8][4] =  EIy_L2;
    k[4][10]= k[10][4]= 2*EIy_L4;
    k[8][10]= k[10][8]=  EIy_L2;

    // 坐标变换：K_global = T^T * K_local * T
    double temp[12][12] = {0};
    double k_global[12][12] = {0};
    // temp = K_local * T
    for(int i=0;i<12;i++)
        for(int j=0;j<12;j++)
            for(int m=0;m<12;m++)
                temp[i][j] += k[i][m] * T[m][j];
    // k_global = T^T * temp
    for(int i=0;i<12;i++)
        for(int j=0;j<12;j++)
            for(int n=0;n<12;n++)
                k_global[i][j] += T[n][i] * temp[n][j];

    // 按带状存储方式写入Matrix（上三角，列优先）
    int idx = 0;
    for(int j=0;j<12;j++)
        for(int i=0;i<=j;i++)
            Matrix[idx++] = k_global[i][j];
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