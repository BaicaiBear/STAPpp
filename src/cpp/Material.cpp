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
    double _;
	Input >> E >> _ >> density >> Area; // Young's modulus, section area, and density

	return true;
}

//	Write material data to Stream
void CBarMaterial::Write(COutputter& output)
{
	output << setw(16) << E << setw(16) << density << Area << setw(16) << endl;
}

// S4R壳单元材料类实现
bool CS4RMaterial::Read(ifstream& Input) {
    Input >> nset;
    Input >> E >> nu >> density >> thickness;
    return true;
}
void CS4RMaterial::Write(COutputter& output) {
    output << setw(16) << E << setw(16) << nu << setw(16) << thickness << setw(16) << density << endl;
}


//	Read material data from stream Input
bool CC3D8RMaterial::Read(ifstream& Input)
{
	Input >> nset;	// Number of property set

	Input >> E >> nu >> density;	// Young's modulus and Poisson's ratio

	return true;
}

//	Write material data to Stream
void CC3D8RMaterial::Write(COutputter& output)
{
	output << setw(16) << E << setw(16) << nu << setw(16) << density << endl;
}

bool CB31Material::Read(ifstream& Input)
{
    Input >> nset;
    double nu;
    double b, d, t1, t2, t3, t4;
    Input >> E >> nu >> density >> b >> d >> t1 >> t2 >> t3 >> t4;
    Area = b * d; // 矩形截面面积
    Iy = Iy = (1.0 / 12.0) * d * pow(b, 3)
          - (1.0 / 12.0) * (d - t1 - t2) * pow(b - t3 - t4, 3);
    Iz = (1.0 / 12.0) * b * pow(d, 3)
          - (1.0 / 12.0) * (b - t3 - t4) * pow(d - t1 - t2, 3);
    J = 4.0 * pow((b - 0.5 * (t3 + t4)) * (d - 0.5 * (t1 + t2)), 2)
          / ((b - t3 - t4) / t1
           + (b - t3 - t4) / t2
           + (d - t1 - t2) / t3
           + (d - t1 - t2) / t4);
    // 计算剪切模量 G
    G = E / (2 * (1 + nu)); // 使用泊松比计算剪切模量
    return true;
}

void CB31Material::Write(COutputter& output)
{
    output << setw(16) << E << setw(16) << G << setw(16) << Area
           << setw(16) << Iy << setw(16) << Iz << setw(16) << J << setw(16) << density << endl;
}

// Q4 单元材料类实现
bool CQ4Material::Read(ifstream& Input)
{
    Input >> nset;
    Input >> E >> nu >> density >> thickness;
    
    int stress_type;
    Input >> stress_type;
    PlaneStress = (stress_type == 1);

    return true;
}

void CQ4Material::Write(COutputter& output)
{
    output << setw(16) << E 
           << setw(16) << nu 
           << setw(16) << density 
           << setw(16) << thickness 
           << setw(16) << (PlaneStress ? 1 : 0) << endl;
}
