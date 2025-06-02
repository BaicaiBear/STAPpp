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

using namespace std;

//	Constructor
CC3D8R::CC3D8R()
{
	NEN_ = 8;  // Each element has 8 nodes
	nodes_ = new CNode*[NEN_];
	
	ND_ = 24; // 3*8
	LocationMatrix_ = new unsigned int[ND_];
	
	ElementMaterial_ = nullptr;
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
void CC3D8R::ElementStiffness(double* Matrix)
{
	clear(Matrix, SizeOfStiffnessMatrix());

	// TODO: Implement the actual stiffness matrix calculation based on the element's geometry and material properties
	CC3D8RMaterial* material_ = dynamic_cast<CC3D8RMaterial*>(ElementMaterial_); // Pointer to material of the element
	double E = material_->E; // Young's modulus
	double nu = material_->nu; // Poisson's ratio
	cerr << "ElementStiffness not implemented yet for CC3D8R." << endl;
}

void CC3D8R::ElementStress(double* stress, double* Displacement)
{
	// TODO: Implement the actual stress calculation based on the displacement and material properties
	CC3D8RMaterial* material_ = dynamic_cast<CC3D8RMaterial*>(ElementMaterial_); // Pointer to material of the element
	cerr << "ElementStress not implemented yet for CC3D8R." << endl;
}