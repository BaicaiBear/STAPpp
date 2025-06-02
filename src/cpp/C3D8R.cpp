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
using namespace Eigen;

//	Constructor
CC3D8R::CC3D8R()
{
	NEN_ = 8;  // Each element has 8 nodes
	nodes_ = new CNode*[NEN_];
	
	ND_ = 24; // 3*8
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
	clear(matrix, SizeOfStiffnessMatrix());

	// TODO: Implement the actual stiffness matrix calculation based on the element's geometry and material properties
	// Get the Material properties
	CC3D8RMaterial* material_ = dynamic_cast<CC3D8RMaterial*>(ElementMaterial_); // Pointer to material of the element
	double E = material_->E; // Young's modulus
	double nu = material_->nu; // Poisson's ratio
	Matrix<double, 6, 6> D;
	double lambda = (E * nu) / ((1 + nu) * (1 - 2 * nu));  // Lame  λ
	double mu     = E / (2 * (1 + nu));                   // G = μ
	D.setZero();
	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
			D(i, j) = lambda;
	for (int i = 0; i < 3; ++i)
		D(i, i) += 2 * mu;
	D(3, 3) = mu;
	D(4, 4) = mu;
	D(5, 5) = mu;

	// Calculate the Jacobian Matrix
	// The Gradient of the shape functions in natural coordinates (ξ, η, ζ) for a C3D8R element
	Matrix<double, 3, 8> shapeFunctionGradNat;
	shapeFunctionGradNat <<
	-0.125,  0.125,  0.125, -0.125, -0.125,  0.125,  0.125, -0.125,  // dN/dxi
	-0.125, -0.125,  0.125,  0.125, -0.125, -0.125,  0.125,  0.125,  // dN/deta
	-0.125, -0.125, -0.125, -0.125,  0.125,  0.125,  0.125,  0.125;  // dN/dzeta
	// The physical coordinates of the nodes in the element
	Matrix<double, 8, 3> physicalCoords;
	for (int i = 0; i < 8; ++i) {
		physicalCoords(i, 0) = nodes_[i]->XYZ[0];  // x
		physicalCoords(i, 1) = nodes_[i]->XYZ[1];  // y
		physicalCoords(i, 2) = nodes_[i]->XYZ[2];  // z
	}
	//Jacobian J = naturalCoords * physicalCoords，3x3
	Matrix3d J = shapeFunctionGradNat * physicalCoords;  // 中心点的缩减积分权重
	Matrix3d Jinv = J.inverse();

	// Calculate the Gradient of the shape functions in physical coordinates
	Matrix<double, 3, 8> shapeFunctionGradXYZ = Jinv * shapeFunctionGradNat;

	// Calculate the B matrix
	Matrix<double, 6, 24> B;
    B.setZero();
    for (int i = 0; i < 8; ++i) {
        int idx = i * 3;
        double dN_dx = shapeFunctionGradXYZ(0, i);
        double dN_dy = shapeFunctionGradXYZ(1, i);
        double dN_dz = shapeFunctionGradXYZ(2, i);

        B(0, idx)     = dN_dx;
        B(1, idx + 1) = dN_dy;
        B(2, idx + 2) = dN_dz;
        B(3, idx)     = dN_dy;
        B(3, idx + 1) = dN_dx;
        B(4, idx + 1) = dN_dz;
        B(4, idx + 2) = dN_dy;
        B(5, idx)     = dN_dz;
        B(5, idx + 2) = dN_dx;
    }
	Matrix<double, 24, 24> Ke = B.transpose() * D * B * (J.determinant() * 8.0);  // 8.0 is the volume factor for the C3D8R element

	const double alpha = 0.01;

	const double gamma[4][8] = {
    {-1, 1, 1, -1, -1, 1, 1, -1},   // γ1
    {-1, -1, 1, 1, -1, -1, 1, 1},   // γ2
    {-1, -1, -1, -1, 1, 1, 1, 1},   // γ3
    { 1, -1, 1, -1, -1, 1, -1, 1}   // γ4
	};

	std::array<Matrix<double, 3, 8>, 4> B_hg_alpha;
	for (int alpha = 0; alpha < 4; ++alpha) {
		for (int i = 0; i < 8; ++i) {
			B_hg_alpha[alpha](0, i) = gamma[alpha][i]; // x
			B_hg_alpha[alpha](1, i) = gamma[alpha][i]; // y
			B_hg_alpha[alpha](2, i) = gamma[alpha][i]; // z
		}
		// 映射到全局坐标
		B_hg_alpha[alpha] = Jinv * B_hg_alpha[alpha];
	}

	Matrix<double, 24, 24> Khg = Matrix<double, 24, 24>::Zero();
	double V = std::abs(J.determinant() * 8.0);
	double l2 = std::pow(std::cbrt(V), 2);
	double gh = mu * V / l2 * alpha;       

	for (int alpha = 0; alpha < 4; ++alpha) {
		Vector<double, 24> B_hg_vec;
		B_hg_vec.setZero();
		for (int i = 0; i < 8; ++i) {
			int idx = i * 3;
			B_hg_vec(idx)     = B_hg_alpha[alpha](0, i);
			B_hg_vec(idx + 1) = B_hg_alpha[alpha](1, i);
			B_hg_vec(idx + 2) = B_hg_alpha[alpha](2, i);
		}
		Khg += gh * B_hg_vec * B_hg_vec.transpose();  // rank-1 outer product
	}
	Matrix<double, 24, 24> Ke_total = Ke + Khg;

	// Store the upper triangular part of the stiffness matrix in the provided array
	int index = 0;
	for (int j = 0; j < 24; ++j) {
		for (int i = 0; i <= j; ++i) {
			matrix[index++] = Ke_total(i, j);
		}
	}
}

void CC3D8R::ElementStress(double* stress, double* Displacement)
{
    // Get material properties from the associated material object
    CC3D8RMaterial* material_ = dynamic_cast<CC3D8RMaterial*>(ElementMaterial_);
    double E = material_->E;     // Young's modulus
    double nu = material_->nu;   // Poisson's ratio

    // Construct the 6x6 constitutive matrix D for isotropic elasticity
    Matrix<double, 6, 6> D;
    double lambda = (E * nu) / ((1 + nu) * (1 - 2 * nu));  // Lamé's first parameter
    double mu     = E / (2 * (1 + nu));                   // Shear modulus
    D.setZero();
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            D(i, j) = lambda;
    for (int i = 0; i < 3; ++i)
        D(i, i) += 2 * mu;
    D(3, 3) = mu;
    D(4, 4) = mu;
    D(5, 5) = mu;

    // Collect nodal coordinates (8 nodes, 3 coordinates each)
    Matrix<double, 8, 3> physicalCoords;
    for (int i = 0; i < 8; ++i) {
        physicalCoords(i, 0) = nodes_[i]->XYZ[0];
        physicalCoords(i, 1) = nodes_[i]->XYZ[1];
        physicalCoords(i, 2) = nodes_[i]->XYZ[2];
    }

    // Gradient of shape functions with respect to natural coordinates (ξ, η, ζ)
    Matrix<double, 3, 8> shapeFunctionGradNat;
    shapeFunctionGradNat <<
        -0.125,  0.125,  0.125, -0.125, -0.125,  0.125,  0.125, -0.125,
        -0.125, -0.125,  0.125,  0.125, -0.125, -0.125,  0.125,  0.125,
        -0.125, -0.125, -0.125, -0.125,  0.125,  0.125,  0.125,  0.125;

    // Compute the Jacobian matrix and its inverse
    Matrix3d J = shapeFunctionGradNat * physicalCoords;
    Matrix3d Jinv = J.inverse();

    // Transform gradient to physical coordinates (∂N/∂x, ∂N/∂y, ∂N/∂z)
    Matrix<double, 3, 8> shapeFunctionGradXYZ = Jinv * shapeFunctionGradNat;

    // Construct the 6x24 strain-displacement matrix B
    Matrix<double, 6, 24> B;
    B.setZero();
    for (int i = 0; i < 8; ++i) {
        int idx = i * 3;
        double dN_dx = shapeFunctionGradXYZ(0, i);
        double dN_dy = shapeFunctionGradXYZ(1, i);
        double dN_dz = shapeFunctionGradXYZ(2, i);

        B(0, idx)     = dN_dx;
        B(1, idx + 1) = dN_dy;
        B(2, idx + 2) = dN_dz;
        B(3, idx)     = dN_dy;
        B(3, idx + 1) = dN_dx;
        B(4, idx + 1) = dN_dz;
        B(4, idx + 2) = dN_dy;
        B(5, idx)     = dN_dz;
        B(5, idx + 2) = dN_dx;
    }

    // Assemble element displacement vector u_e from global displacement vector
    Vector<double, 24> u_e;
    for (int i = 0; i < 24; ++i) {
        u_e(i) = (LocationMatrix_[i] ? Displacement[LocationMatrix_[i] - 1] : 0.0);
    }

    // Compute stress using σ = D * B * u_e
    Vector<double, 6> sigma = D * (B * u_e);

    // Store computed stress in output array
    for (int i = 0; i < 6; ++i) {
        stress[i] = sigma(i);
    }
}
