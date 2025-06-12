/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Bar.h"

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

//	Constructor
CBar::CBar()
{
	NEN_ = 2;	// Each element has 2 nodes
	nodes_ = new CNode*[NEN_];
    
    ND_ = 12;
    LocationMatrix_ = new unsigned int[ND_];

	ElementMaterial_ = nullptr;
}

//	Desconstructor
CBar::~CBar()
{
}

//	Read element data from stream Input
bool CBar::Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList)
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
void CBar::Write(COutputter& output)
{
	output << setw(11) << nodes_[0]->NodeNumber
		   << setw(9) << nodes_[1]->NodeNumber << setw(12) << ElementMaterial_->nset << endl;
}

//	Calculate element stiffness matrix 
//	Upper triangular matrix, stored as an array column by colum starting from the diagonal element
void CBar::ElementStiffness(double* Matrix)
{
	clear(Matrix, SizeOfStiffnessMatrix());

//	Calculate bar length
	double DX[3];		//	dx = x2-x1, dy = y2-y1, dz = z2-z1
	for (unsigned int i = 0; i < 3; i++)
		DX[i] = nodes_[1]->XYZ[i] - nodes_[0]->XYZ[i];

	double DX2[6];	//  Quadratic polynomial (dx^2, dy^2, dz^2, dx*dy, dy*dz, dx*dz)
	DX2[0] = DX[0] * DX[0];
	DX2[1] = DX[1] * DX[1];
	DX2[2] = DX[2] * DX[2];
	DX2[3] = DX[0] * DX[1];
	DX2[4] = DX[1] * DX[2];
	DX2[5] = DX[0] * DX[2];

	double L2 = DX2[0] + DX2[1] + DX2[2];
	double L = sqrt(L2);

//	Calculate element stiffness matrix

	CBarMaterial* material_ = dynamic_cast<CBarMaterial*>(ElementMaterial_);	// Pointer to material of the element

	double k = material_->E * material_->Area / L / L2;

	double matrix[21];

	matrix[0] = k*DX2[0];
	matrix[1] = k*DX2[1];
	matrix[2] = k*DX2[3];
	matrix[3] = k*DX2[2];
	matrix[4] = k*DX2[4];
	matrix[5] = k*DX2[5];
	matrix[6] = k*DX2[0];
	matrix[7] = -k*DX2[5];
	matrix[8] = -k*DX2[3];
	matrix[9] = -k*DX2[0];
	matrix[10] = k*DX2[1];
	matrix[11] = k*DX2[3];
	matrix[12] = -k*DX2[4];
	matrix[13] = -k*DX2[1];
	matrix[14] = -k*DX2[3];
	matrix[15] = k*DX2[2];
	matrix[16] = k*DX2[4];
	matrix[17] = k*DX2[5];
	matrix[18] = -k*DX2[2];
	matrix[19] = -k*DX2[4];
	matrix[20] = -k*DX2[5];

	for (int i = 0; i < NEN_; ++i) {
        for (int j = 0; j < 3; ++j) {
            int f = i*3+j, realf = 6*i+j;
            for (int p = i; p >= 0; --p){
                for (int q = 2; q >= 0; --q) {
                    int k = p*3+q, realk = 6*p+q;
                    if (k > f) continue; 
					Matrix[realf*(realf+1)/2+realf-realk] = matrix[f*(f+1)/2+f-k];
                }
            }
        }
    }
}

//	Calculate element stress 
void CBar::ElementStress(double* stress, double* Displacement)
{
	CBarMaterial* material_ = dynamic_cast<CBarMaterial*>(ElementMaterial_);	// Pointer to material of the element

	double DX[3];	//	dx = x2-x1, dy = y2-y1, dz = z2-z1
	double L2 = 0;	//	Square of bar length (L^2)

	for (unsigned int i = 0; i < 3; i++)
	{
		DX[i] = nodes_[1]->XYZ[i] - nodes_[0]->XYZ[i];
		L2 = L2 + DX[i]*DX[i];
	}

	double S[6];
	for (unsigned int i = 0; i < 3; i++)
	{
		S[i] = -DX[i] * material_->E / L2;
		S[i+3] = -S[i];
	}
	
	*stress = 0.0;
	for (unsigned int i = 0; i < 6; i++)
	{
		if (LocationMatrix_[i])
			*stress += S[i] * Displacement[LocationMatrix_[i]-1];
	}
}
