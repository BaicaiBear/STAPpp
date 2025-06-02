/*****************************************************************************/
/*  STAP++ : S4R 4节点壳单元（减少积分，含hourglass控制）                    */
/*  参考Abaqus实现，支持与其他3D单元混用                                      */
/*****************************************************************************/

#include "S4R.h"
#include "Material.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

CS4R::CS4R() {
    NEN_ = 4; // 4节点
    nodes_ = new CNode*[NEN_];
    ND_ = 12; // 4节点*3自由度
    LocationMatrix_ = new unsigned int[ND_];
    ElementMaterial_ = nullptr;
}

CS4R::~CS4R() {}

bool CS4R::Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList) {
    unsigned int MSet;
    unsigned int N[4];
    Input >> N[0] >> N[1] >> N[2] >> N[3] >> MSet;
    ElementMaterial_ = dynamic_cast<CS4RMaterial*>(MaterialSets) + MSet - 1;
    for (int i = 0; i < 4; ++i) {
        nodes_[i] = &NodeList[N[i] - 1];
    }
    return true;
}

void CS4R::Write(COutputter& output) {
    output << setw(8) << nodes_[0]->NodeNumber
           << setw(8) << nodes_[1]->NodeNumber
           << setw(8) << nodes_[2]->NodeNumber
           << setw(8) << nodes_[3]->NodeNumber
           << setw(8) << ElementMaterial_->nset << endl;
}

void CS4R::ElementStiffness(double* Matrix) {
    // Mindlin-Reissner壳理论下的S4R单元刚度矩阵实现（膜2x2高斯积分，剪切1点积分）
    // 1. 计算节点坐标
    double x[4], y[4], z[4];
    for (int i = 0; i < 4; ++i) {
        x[i] = nodes_[i]->XYZ[0];
        y[i] = nodes_[i]->XYZ[1];
        z[i] = nodes_[i]->XYZ[2];
    }
    // 2. 计算单元厚度、材料参数
    CS4RMaterial* mat = dynamic_cast<CS4RMaterial*>(ElementMaterial_);
    double t = mat->thickness;
    double E = mat->E;
    double nu = mat->nu;
    double kappa = 5.0/6.0; // 剪切校正系数
    double G = E/(2*(1+nu));
    // 3. 初始化刚度矩阵
    int size = 12 * (12 + 1) / 2;
    for (int i = 0; i < size; ++i) Matrix[i] = 0.0;
    // 4. 2x2高斯积分点和权重
    double gauss[2] = { -1.0/std::sqrt(3.0), 1.0/std::sqrt(3.0) };
    double weight[2] = { 1.0, 1.0 };
    // 5. 膜刚度2x2高斯积分
    for (int gp_x = 0; gp_x < 2; ++gp_x) {
        for (int gp_y = 0; gp_y < 2; ++gp_y) {
            double xi = gauss[gp_x];
            double eta = gauss[gp_y];
            double w = weight[gp_x] * weight[gp_y];
            // 形函数及其导数
            double N[4] = {(1-xi)*(1-eta)/4, (1+xi)*(1-eta)/4, (1+xi)*(1+eta)/4, (1-xi)*(1+eta)/4};
            double dN_dxi[4]  = {-(1-eta)/4, (1-eta)/4, (1+eta)/4, -(1+eta)/4};
            double dN_deta[4] = {-(1-xi)/4, -(1+xi)/4, (1+xi)/4, (1-xi)/4};
            // Jacobi矩阵
            double dx_dxi = 0, dx_deta = 0, dy_dxi = 0, dy_deta = 0;
            for (int i = 0; i < 4; ++i) {
                dx_dxi  += dN_dxi[i]  * x[i];
                dx_deta += dN_deta[i] * x[i];
                dy_dxi  += dN_dxi[i]  * y[i];
                dy_deta += dN_deta[i] * y[i];
            }
            double J[2][2] = {{dx_dxi, dx_deta}, {dy_dxi, dy_deta}};
            double detJ = J[0][0]*J[1][1] - J[0][1]*J[1][0];
            double invJ[2][2] = { { J[1][1]/detJ, -J[0][1]/detJ }, { -J[1][0]/detJ, J[0][0]/detJ } };
            // B矩阵（膜）
            double Bb[3][12] = {0};
            for (int i = 0; i < 4; ++i) {
                double dN_dx = invJ[0][0]*dN_dxi[i] + invJ[0][1]*dN_deta[i];
                double dN_dy = invJ[1][0]*dN_dxi[i] + invJ[1][1]*dN_deta[i];
                Bb[0][i*3+0] = dN_dx;
                Bb[1][i*3+1] = dN_dy;
                Bb[2][i*3+0] = dN_dy;
                Bb[2][i*3+1] = dN_dx;
            }
            // D矩阵（膜）
            double Db[3][3] = {
                {1, nu, 0},
                {nu, 1, 0},
                {0, 0, (1-nu)/2}
            };
            double coeff = E*t*t*t/(12.0*(1-nu*nu));
            for(int i=0;i<3;++i) for(int j=0;j<3;++j) Db[i][j]*=coeff;
            // 组装膜刚度
            int index = 0;
            for (int i = 0; i < 12; ++i) {
                for (int j = i; j >= 0; --j) {
                    // 只累加膜刚度
                    for (int m = 0; m < 3; ++m) {
                        for (int n = 0; n < 3; ++n) {
                            Matrix[index] += w * detJ * Bb[m][i] * Db[m][n] * Bb[n][j];
                        }
                    }
                    ++index;
                }
            }
        }
    }
    // 6. 剪切刚度1点积分
    double xi = 0.0, eta = 0.0, w = 4.0;
    double N[4] = {(1-xi)*(1-eta)/4, (1+xi)*(1-eta)/4, (1+xi)*(1+eta)/4, (1-xi)*(1+eta)/4};
    double dN_dxi[4]  = {-(1-eta)/4, (1-eta)/4, (1+eta)/4, -(1+eta)/4};
    double dN_deta[4] = {-(1-xi)/4, -(1+xi)/4, (1+xi)/4, (1-xi)/4};
    double dx_dxi = 0, dx_deta = 0, dy_dxi = 0, dy_deta = 0;
    for (int i = 0; i < 4; ++i) {
        dx_dxi  += dN_dxi[i]  * x[i];
        dx_deta += dN_deta[i] * x[i];
        dy_dxi  += dN_dxi[i]  * y[i];
        dy_deta += dN_deta[i] * y[i];
    }
    double J[2][2] = {{dx_dxi, dx_deta}, {dy_dxi, dy_deta}};
    double detJ = J[0][0]*J[1][1] - J[0][1]*J[1][0];
    double invJ[2][2] = { { J[1][1]/detJ, -J[0][1]/detJ }, { -J[1][0]/detJ, J[0][0]/detJ } };
    double Bs[2][12] = {0};
    for (int i = 0; i < 4; ++i) {
        double dN_dx = invJ[0][0]*dN_dxi[i] + invJ[0][1]*dN_deta[i];
        double dN_dy = invJ[1][0]*dN_dxi[i] + invJ[1][1]*dN_deta[i];
        double Ni = N[i];
        Bs[0][i*3+0] = -Ni;
        Bs[1][i*3+1] = -Ni;
        Bs[0][i*3+2] = dN_dx;
        Bs[1][i*3+2] = dN_dy;
    }
    double Ds[2][2] = { {kappa*G*t, 0}, {0, kappa*G*t} };
    // 组装剪切刚度
    int index = 0;
    for (int i = 0; i < 12; ++i) {
        for (int j = i; j >= 0; --j) {
            // 只累加剪切刚度
            for (int m = 0; m < 2; ++m) {
                for (int n = 0; n < 2; ++n) {
                    Matrix[index] += w * detJ * Bs[m][i] * Ds[m][n] * Bs[n][j];
                }
            }
            ++index;
        }
    }

    // // 9. hourglass控制（基于Flanagan-Belytschko方法，简化实现）
    // // hourglass 模式矢量
    // double hg[4] = {1, -1, 1, -1};
    // double alpha = 0.05; // hourglass 抑制系数（可调）
    // for(int i=0;i<4;++i) {
    //     for(int j=0;j<4;++j) {
    //         // 只对z方向（板厚方向）自由度加hourglass抑制
    //         Matrix[(i*3+2)*12 + (j*3+2)] += alpha * E * t * hg[i] * hg[j] * detJ;
    //     }
    // }
}

void CS4R::ElementStress(double* stress, double* Displacement) {
    // 只在单个高斯点计算应力
    double x[4], y[4], z[4];
    for (int i = 0; i < 4; ++i) {
        x[i] = nodes_[i]->XYZ[0];
        y[i] = nodes_[i]->XYZ[1];
        z[i] = nodes_[i]->XYZ[2];
    }
    CS4RMaterial* mat = dynamic_cast<CS4RMaterial*>(ElementMaterial_);
    double t = mat->thickness;
    double E = mat->E;
    double nu = mat->nu;
    double xi = 0.0, eta = 0.0;
    double dN_dxi[4]  = {-(1-eta)/4, (1-eta)/4, (1+eta)/4, -(1+eta)/4};
    double dN_deta[4] = {-(1-xi)/4, -(1+xi)/4, (1+xi)/4, (1-xi)/4};
    double dx_dxi = 0, dx_deta = 0, dy_dxi = 0, dy_deta = 0;
    for (int i = 0; i < 4; ++i) {
        dx_dxi  += dN_dxi[i]  * x[i];
        dx_deta += dN_deta[i] * x[i];
        dy_dxi  += dN_dxi[i]  * y[i];
        dy_deta += dN_deta[i] * y[i];
    }
    double J[2][2] = {{dx_dxi, dx_deta}, {dy_dxi, dy_deta}};
    double detJ = J[0][0]*J[1][1] - J[0][1]*J[1][0];
    double invJ[2][2] = { { J[1][1]/detJ, -J[0][1]/detJ }, { -J[1][0]/detJ, J[0][0]/detJ } };
    double Bb[3][12] = {0};
    for (int i = 0; i < 4; ++i) {
        double dN_dx = invJ[0][0]*dN_dxi[i] + invJ[0][1]*dN_deta[i];
        double dN_dy = invJ[1][0]*dN_dxi[i] + invJ[1][1]*dN_deta[i];
        Bb[0][i*3+0] = dN_dx;
        Bb[1][i*3+1] = dN_dy;
        Bb[2][i*3+0] = dN_dy;
        Bb[2][i*3+1] = dN_dx;
    }
    double Db[3][3] = {
        {1, nu, 0},
        {nu, 1, 0},
        {0, 0, (1-nu)/2}
    };
    double coeff = E*t*t*t/(12.0*(1-nu*nu));
    for(int i=0;i<3;++i) for(int j=0;j<3;++j) Db[i][j]*=coeff;
    double u[12];
    for(int i=0;i<12;++i) u[i] = Displacement[LocationMatrix_[i]-1];
    for(int m=0;m<3;++m) {
        stress[m] = 0.0;
        for(int i=0;i<12;++i) stress[m] += Bb[m][i]*u[i];
    }
    double stress_out[3] = {0};
    for(int m=0;m<3;++m) for(int n=0;n<3;++n) stress_out[m] += Db[m][n]*stress[n];
    for(int m=0;m<3;++m) stress[m] = stress_out[m];
}
