/*****************************************************************************/
/*  STAP++ : B31 两结点空间线性梁单元                                          */
/*  Author : Yanzi Wang                                                      */
/*****************************************************************************/

#pragma once

#include "Element.h"

using namespace std;

//! B31 element class
class CB31 : public CElement
{
public:

//!	Constructor
	CB31();

//!	Desconstructor
	~CB31();

//!	Read element data from stream Input
	virtual bool Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList);

//!	Write element data to stream
	virtual void Write(COutputter& output);

//!	Calculate element stiffness matrix
	virtual void ElementStiffness(double* Matrix);

//!	Calculate element stress
	virtual void ElementStress(double* stress, double* Displacement);

private:
    void GetTransformationMatrix(double T[12][12], double& L) const;
};