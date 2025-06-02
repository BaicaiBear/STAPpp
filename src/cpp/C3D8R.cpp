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
	Matrix3d J = shapeFunctionGradNat * physicalCoords * (1.0 / 8.0);  // 中心点的缩减积分权重
	Matrix3d Jinv = J.inverse();

	// Calculate the Gradient of the shape functions in physical coordinates
	Matrix<double, 3, 8> shapeFunctionGradXYZ = Jinv * shapeFunctionGradNat;

	// Calculate the B matrix
	Matrix<double, 6, 24> B;
	for (int i = 0; i < 8; ++i) {
		int idx = i * 3;
		B(0, idx)     = shapeFunctionGradXYZ(0, i);  // dN/dx
		B(1, idx + 1) = shapeFunctionGradXYZ(1, i);  // dN/dy
		B(2, idx + 2) = shapeFunctionGradXYZ(2, i);  // dN/dz
		B(3, idx)     = shapeFunctionGradXYZ(1, i);  // dN/dy
		B(3, idx + 1) = shapeFunctionGradXYZ(0, i);  // dN/dx
		B(4, idx + 1) = shapeFunctionGradXYZ(2, i);  // dN/dz
		B(4, idx + 2) = shapeFunctionGradXYZ(1, i);  // dN/dy
		B(5, idx)     = shapeFunctionGradXYZ(2, i);  // dN/dz
		B(5, idx + 2) = shapeFunctionGradXYZ(0, i);  // dN/dx
	}
	Matrix<double, 24, 24> Ke = B.transpose() * D * B * (J.determinant() * 8.0);  // 8.0 is the volume factor for the C3D8R element

	// Store the upper triangular part of the stiffness matrix in the provided array
	int index = 0;
	for (int j = 0; j < 24; ++j) {
		for (int i = 0; i <= j; ++i) {
			matrix[index++] = Ke(i, j);
		}
	}
}

void CC3D8R::ElementStress(double* stress, double* Displacement)
{
	// TODO: Implement the actual stress calculation based on the displacement and material properties
	CC3D8RMaterial* material_ = dynamic_cast<CC3D8RMaterial*>(ElementMaterial_); // Pointer to material of the element
	cerr << "ElementStress not implemented yet for CC3D8R." << endl;
}