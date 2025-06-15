/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#pragma once

#include "Outputter.h"

using namespace std;

//!	Material base class which only define one data member
/*!	All type of material classes should be derived from this base class */
class CMaterial
{
public:

	unsigned int nset;	//!< Number of set
	
	double E;  //!< Young's modulus

public:

//! Virtual deconstructor
    virtual ~CMaterial() {};

//!	Read material data from stream Input
	virtual bool Read(ifstream& Input) = 0;

//!	Write material data to Stream
    virtual void Write(COutputter& output) = 0;

};

//!	Material class for bar element
class CBarMaterial : public CMaterial
{
public:

	double Area;	//!< Sectional area of a bar element
    double density; //!< 材料密度

public:
	
//!	Read material data from stream Input
	virtual bool Read(ifstream& Input);

//!	Write material data to Stream
	virtual void Write(COutputter& output);
};

// S4R壳单元材料类
class CS4RMaterial : public CMaterial
{
public:
    double nu;    // 泊松比
    double thickness; // 单元厚度
    double density; // 材料密度

    virtual bool Read(ifstream& Input) override;
    virtual void Write(COutputter& output) override;
};

class CC3D8RMaterial : public CMaterial
{
public:
	// Young's modulus
	double E;
	// Poisson's ratio
	double nu;
	// 材料密度
	double density;

public:

//!	Read material data from stream Input
	virtual bool Read(ifstream& Input);

//!	Write material data to Stream
	virtual void Write(COutputter& output);
};

//! Material class for B31 space beam element
class CB31Material : public CMaterial
{
public:
    double G;    //!< Shear modulus
    double Area; //!< Sectional area
    double Iy;   //!< Moment of inertia about y
    double Iz;   //!< Moment of inertia about z
    double J;    //!< Torsional constant
    double density; // 材料密度

public:
    virtual bool Read(ifstream& Input);
    virtual void Write(COutputter& output);
};

// Q4 单元材料类定义
class CQ4Material : public CMaterial
{
public:
    double nu;           // 泊松比
    double thickness;    // 厚度
    double density;      // 材料密度
    bool PlaneStress;    // 是否为平面应力状态

    virtual bool Read(ifstream& Input) override;
    virtual void Write(COutputter& output) override;
};
