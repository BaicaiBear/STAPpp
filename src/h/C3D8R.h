/*
################################################
Author: Jinkun Liu
Email: liujk22@mails.tsinghua.edu.cn
################################################
*/
#pragma once

#include "Element.h"

using namespace std;

// C3D8R element class
class CC3D8R : public CElement
{
public:

//!	Constructor
	CC3D8R();

//!	Desconstructor
	~CC3D8R();

//!	Read element data from stream Input
	virtual bool Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList);

//!	Write element data to stream
	virtual void Write(COutputter& output);

//!	Calculate element stiffness matrix
	virtual void ElementStiffness(double* Matrix);

//!	Calculate element stress
	virtual void ElementStress(double* stress, double* Displacement);
};
