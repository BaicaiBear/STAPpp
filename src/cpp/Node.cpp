/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include <iostream>
#include <iomanip>

#include "Node.h"

CNode::CNode(double X, double Y, double Z)
{
    XYZ[0] = X;        // Coordinates of the node
    XYZ[1] = Y;
    XYZ[2] = Z;
    for (int i = 0; i < NDF; ++i) {
        bcode[i] = 0;
        bcval[i] = 0.0;
    }
    NodeNumber = 0;    // Node number
};

//  Read element data from stream Input
bool CNode::Read(ifstream& Input)
{
    Input >> NodeNumber;
    for (int i = 0; i < NDF; ++i) {
        Input >> bcode[i] >> bcval[i];
    }
    Input >> XYZ[0] >> XYZ[1] >> XYZ[2];
    return true;
}

//  Output nodal point data to stream
void CNode::Write(COutputter& output)
{
    output << setw(9) << NodeNumber;
    for (int i = 0; i < NDF; ++i)
        output << setw(5) << bcode[i];
    for (int i = 0; i < NDF; ++i)
        output << setw(10) << bcval[i];
    output << setw(18) << XYZ[0] << setw(15) << XYZ[1] << setw(15) << XYZ[2] << endl;
}

//	Output equation numbers of nodal point to stream
void CNode::WriteEquationNo(COutputter& output)
{
	output << setw(9) << NodeNumber << "       ";

	for (unsigned int dof = 0; dof < CNode::NDF; dof++)	// Loop over for DOFs of node np
	{
		output << setw(5) << bcode[dof];
	}

	output << endl;
}

//	Write nodal displacement
void CNode::WriteNodalDisplacement(COutputter& output, double* Displacement)
{
    output << setw(5) << NodeNumber << "        ";
    for (unsigned int j = 0; j < NDF; j++)
    {
        if (bcode[j] == 0)
        {
            output << setw(18) << bcval[j];
        }
        else
        {
            output << setw(18) << Displacement[bcode[j] - 1];
        }
    }
    output << endl;
}
