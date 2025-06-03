/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Material.h"

#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

//	Read material data from stream Input
bool CBarMaterial::Read(ifstream& Input)
{
	Input >> nset;	// Number of property set

	Input >> E >> Area;	// Young's modulus and section area

	return true;
}

//	Write material data to Stream
void CBarMaterial::Write(COutputter& output)
{
	output << setw(16) << E << setw(16) << Area << endl;
}

// S4R壳单元材料类实现
bool CS4RMaterial::Read(ifstream& Input) {
    Input >> nset;
    Input >> E >> nu >> thickness;
    return true;
}
void CS4RMaterial::Write(COutputter& output) {
    output << setw(16) << E << setw(16) << nu << setw(16) << thickness << setw(16) << endl;
}


//	Read material data from stream Input
bool CC3D8RMaterial::Read(ifstream& Input)
{
	Input >> nset;	// Number of property set

	Input >> E >> nu;	// Young's modulus and Poisson's ratio

	return true;
}

//	Write material data to Stream
void CC3D8RMaterial::Write(COutputter& output)
{
	output << setw(16) << E << setw(16) << nu << endl;
}