#include "Q4.h"
#include "Material.h"
#include "Node.h"
#include "Domain.h"

#include <cmath>
#include <iostream>

using namespace std;

CQ4::CQ4()
{
    NEN_ = 4;
    ND_ = 8;
    nodes_ = new CNode*[NEN_];
    LocationMatrix_ = new unsigned int[ND_];
    ElementMaterial_ = nullptr;
}

CQ4::~CQ4()
{
    delete[] nodes_;
    delete[] LocationMatrix_;
}

bool CQ4::Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList)
{
    unsigned int N1, N2, N3, N4, MSet;
    Input >> N1 >> N2 >> N3 >> N4 >> MSet;
    nodes_[0] = &NodeList[N1 - 1];
    nodes_[1] = &NodeList[N2 - 1];
    nodes_[2] = &NodeList[N3 - 1];
    nodes_[3] = &NodeList[N4 - 1];
    ElementMaterial_ = &MaterialSets[MSet - 1];
    return true;
}

void CQ4::Write(COutputter& output)
{
    output << setw(11) << nodes_[0]->NodeNumber
           << setw(9)  << nodes_[1]->NodeNumber
           << setw(9)  << nodes_[2]->NodeNumber
           << setw(9)  << nodes_[3]->NodeNumber
           << setw(12) << ElementMaterial_->nset << endl;
}

void CQ4::ElementStiffness(double* Matrix)
{
    CQ4Material* mat = dynamic_cast<CQ4Material*>(ElementMaterial_);
    if (!mat)
    {
        cerr << "Error: ElementMaterial is not of type CQ4Material." << endl;
        exit(-1);
    }
    double E = mat->E;
    double nu = mat->nu;
    double t = mat->thickness;
    bool plane_stress = mat->PlaneStress;
    double D[3][3] = { 0 };
    if (plane_stress)
    {
        double coeff = E / (1 - nu * nu);
        D[0][0] = D[1][1] = coeff;
        D[0][1] = D[1][0] = nu * coeff;
        D[2][2] = (1 - nu) / 2 * coeff;
    }
    else
    {
        double coeff = E * (1 - nu) / ((1 + nu) * (1 - 2 * nu));
        D[0][0] = D[1][1] = coeff;
        D[0][1] = D[1][0] = nu / (1 - nu) * coeff;
        D[2][2] = (1 - 2 * nu) / (2 * (1 - nu)) * coeff;
    }

    double x[4], y[4];
    for (int i = 0; i < 4; ++i)
    {
        x[i] = nodes_[i]->XYZ[0];
        y[i] = nodes_[i]->XYZ[1];
    }

    double gauss[2] = { -1.0 / sqrt(3.0), 1.0 / sqrt(3.0) };
    double weight[2] = { 1.0, 1.0 };
    double K[8][8] = { 0 };
    for (int i = 0; i < 2; ++i)
    {
        for (int j = 0; j < 2; ++j)
        {
            double xi = gauss[i];
            double eta = gauss[j];
            double w = weight[i] * weight[j];
            double dN_dxi[4]  = {-(1 - eta) / 4, (1 - eta) / 4, (1 + eta) / 4, -(1 + eta) / 4};
            double dN_deta[4] = {-(1 - xi) / 4, -(1 + xi) / 4, (1 + xi) / 4, (1 - xi) / 4};
            double dx_dxi = 0, dx_deta = 0, dy_dxi = 0, dy_deta = 0;
            for (int k = 0; k < 4; ++k)
            {
                dx_dxi  += dN_dxi[k] * x[k];
                dx_deta += dN_deta[k] * x[k];
                dy_dxi  += dN_dxi[k] * y[k];
                dy_deta += dN_deta[k] * y[k];
            }

            double J[2][2] = { {dx_dxi, dy_dxi}, {dx_deta, dy_deta} };
            double detJ = J[0][0] * J[1][1] - J[0][1] * J[1][0];
            double invJ[2][2] = {
                { J[1][1] / detJ, -J[0][1] / detJ },
                { -J[1][0] / detJ, J[0][0] / detJ }
            };

            double B[3][8] = { 0 };
            for (int k = 0; k < 4; ++k)
            {
                double dN_dx = invJ[0][0] * dN_dxi[k] + invJ[0][1] * dN_deta[k];
                double dN_dy = invJ[1][0] * dN_dxi[k] + invJ[1][1] * dN_deta[k];
                B[0][2 * k]     = dN_dx;
                B[1][2 * k + 1] = dN_dy;
                B[2][2 * k]     = dN_dy;
                B[2][2 * k + 1] = dN_dx;
            }

            for (int m = 0; m < 8; ++m)
            {
                for (int n = 0; n < 8; ++n)
                {
                    double sum = 0.0;
                    for (int p = 0; p < 3; ++p)
                        for (int q = 0; q < 3; ++q)
                            sum += B[p][m] * D[p][q] * B[q][n];
                    K[m][n] += sum * detJ * t * w;
                }
            }
        }
    }
    for (int i = 0; i < 8; ++i)
        for (int j = 0; j < 8; ++j)
            Matrix[i * 8 + j] = K[i][j];
}

void CQ4::ElementStress(double* stress, double* Displacement)
{
    CQ4Material* mat = dynamic_cast<CQ4Material*>(ElementMaterial_);
    if (!mat)
    {
        cerr << "Error: ElementMaterial is not of type CQ4Material." << endl;
        exit(-1);
    }

    double E = mat->E;
    double nu = mat->nu;
    bool plane_stress = mat->PlaneStress;
    double D[3][3] = { 0 };
    if (plane_stress)
    {
        double coeff = E / (1 - nu * nu);
        D[0][0] = D[1][1] = coeff;
        D[0][1] = D[1][0] = nu * coeff;
        D[2][2] = (1 - nu) / 2 * coeff;
    }
    else
    {
        double coeff = E * (1 - nu) / ((1 + nu) * (1 - 2 * nu));
        D[0][0] = D[1][1] = coeff;
        D[0][1] = D[1][0] = nu / (1 - nu) * coeff;
        D[2][2] = (1 - 2 * nu) / (2 * (1 - nu)) * coeff;
    }

    double x[4], y[4];
    for (int i = 0; i < 4; ++i)
    {
        x[i] = nodes_[i]->XYZ[0];
        y[i] = nodes_[i]->XYZ[1];
    }

    double xi = 0.0, eta = 0.0;
    double dN_dxi[4]  = {-(1 - eta) / 4, (1 - eta) / 4, (1 + eta) / 4, -(1 + eta) / 4};
    double dN_deta[4] = {-(1 - xi) / 4, -(1 + xi) / 4, (1 + xi) / 4, (1 - xi) / 4};
    double dx_dxi = 0, dx_deta = 0, dy_dxi = 0, dy_deta = 0;
    for (int i = 0; i < 4; ++i)
    {
        dx_dxi  += dN_dxi[i]  * x[i];
        dx_deta += dN_deta[i] * x[i];
        dy_dxi  += dN_dxi[i]  * y[i];
        dy_deta += dN_deta[i] * y[i];
    }

    double J[2][2] = { {dx_dxi, dy_dxi}, {dx_deta, dy_deta} };
    double detJ = J[0][0] * J[1][1] - J[0][1] * J[1][0];
    double invJ[2][2] = {
        { J[1][1] / detJ, -J[0][1] / detJ },
        { -J[1][0] / detJ, J[0][0] / detJ }
    };

    double B[3][8] = { 0 };
    for (int i = 0; i < 4; ++i)
    {
        double dN_dx = invJ[0][0] * dN_dxi[i] + invJ[0][1] * dN_deta[i];
        double dN_dy = invJ[1][0] * dN_dxi[i] + invJ[1][1] * dN_deta[i];
        B[0][2 * i]     = dN_dx;
        B[1][2 * i + 1] = dN_dy;
        B[2][2 * i]     = dN_dy;
        B[2][2 * i + 1] = dN_dx;
    }

    double u[8] = { 0.0 };
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 2; ++j)
        {
            int idx = LocationMatrix_[i * 2 + j];
            u[i * 2 + j] = (idx > 0) ? Displacement[idx - 1] : 0.0;
        }

    double strain[3] = { 0 };
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 8; ++j)
            strain[i] += B[i][j] * u[j];

    for (int i = 0; i < 3; ++i)
    {
        stress[i] = 0.0;
        for (int j = 0; j < 3; ++j)
            stress[i] += D[i][j] * strain[j];
    }
}

CNode* CQ4::GetNode(unsigned int i) const
{
    return nodes_[i];
}
