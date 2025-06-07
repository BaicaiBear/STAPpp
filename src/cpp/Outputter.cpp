/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include <ctime>
#include <cmath>
#include <set>
#include "S4R.h"

#include "Domain.h"
#include "Outputter.h"
#include "SkylineMatrix.h"

using namespace std;

//	Output current time and date
void COutputter::PrintTime(const struct tm* ptm, COutputter &output)
{
	const char* weekday[] = {"Sunday", "Monday", "Tuesday", "Wednesday",
							 "Thursday", "Friday", "Saturday"};
	const char* month[] = {"January", "February", "March", "April", "May", "June",
						   "July", "August", "September", "October", "November", "December"};

	output << "        (";
	output << ptm->tm_hour << ":" << ptm->tm_min << ":" << ptm->tm_sec << " on ";
	output << month[ptm->tm_mon] << " " << ptm->tm_mday << ", " << ptm->tm_year + 1900 << ", "
		   << weekday[ptm->tm_wday] << ")" << endl
		   << endl;
}

COutputter* COutputter::_instance = nullptr;

//	Constructor
COutputter::COutputter(string FileName)
{
	OutputFile.open(FileName);

	if (!OutputFile)
	{
		cerr << "*** Error *** File " << FileName << " does not exist !" << endl;
		exit(3);
	}
}

//	Return the single instance of the class
COutputter* COutputter::GetInstance(string FileName)
{
	if (!_instance)
		_instance = new COutputter(FileName);
    
	return _instance;
}

//	Print program logo
void COutputter::OutputHeading()
{
	CDomain* FEMData = CDomain::GetInstance();

	*this << "TITLE : " << FEMData->GetTitle() << endl;

	time_t rawtime;
	struct tm* timeinfo;

	time(&rawtime);
	timeinfo = localtime(&rawtime);

	PrintTime(timeinfo, *this);
}

//	Print nodal data
void COutputter::OutputNodeInfo()
{
	CDomain* FEMData = CDomain::GetInstance();

	CNode* NodeList = FEMData->GetNodeList();

	*this << "C O N T R O L   I N F O R M A T I O N" << endl
		  << endl;

	*this << setiosflags(ios::scientific) << setprecision(5);

	unsigned int NUMNP = FEMData->GetNUMNP();
	unsigned int NUMEG = FEMData->GetNUMEG();
	unsigned int NLCASE = FEMData->GetNLCASE();
	unsigned int MODEX = FEMData->GetMODEX();

	*this << "      NUMBER OF NODAL POINTS . . . . . . . . . . (NUMNP)  =" << setw(6) << NUMNP << endl;
	*this << "      NUMBER OF ELEMENT GROUPS . . . . . . . . . (NUMEG)  =" << setw(6) << NUMEG << endl;
	*this << "      NUMBER OF LOAD CASES . . . . . . . . . . . (NLCASE) =" << setw(6) << NLCASE << endl;
	*this << "      SOLUTION MODE  . . . . . . . . . . . . . . (MODEX)  =" << setw(6) << MODEX << endl;
	*this << "         EQ.0, DATA CHECK" << endl
		  << "         EQ.1, EXECUTION" << endl
		  << endl;

	*this << " N O D A L   P O I N T   D A T A" << endl << endl;
	*this << "    NODE       BOUNDARY                         NODAL POINT" << endl
		  << "   NUMBER  CONDITION  CODES                     COORDINATES" << endl;

	for (unsigned int np = 0; np < NUMNP; np++)
		NodeList[np].Write(*this);

	*this << endl;
}

//	Output equation numbers
void COutputter::OutputEquationNumber()
{
	CDomain* FEMData = CDomain::GetInstance();
	unsigned int NUMNP = FEMData->GetNUMNP();

	CNode* NodeList = FEMData->GetNodeList();

	*this << " EQUATION NUMBERS" << endl
		  << endl;
	*this << "   NODE NUMBER   DEGREES OF FREEDOM" << endl;
	*this << "        N           X    Y    Z" << endl;

	for (unsigned int np = 0; np < NUMNP; np++) // Loop over for all node
		NodeList[np].WriteEquationNo(*this);

	*this << endl;
}

//	Output element data
void COutputter::OutputElementInfo()
{
	//	Print element group control line

	CDomain* FEMData = CDomain::GetInstance();

	unsigned int NUMEG = FEMData->GetNUMEG();

	*this << " E L E M E N T   G R O U P   D A T A" << endl
		  << endl
		  << endl;

	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)
	{
		*this << " E L E M E N T   D E F I N I T I O N" << endl
			  << endl;

		ElementTypes ElementType = FEMData->GetEleGrpList()[EleGrp].GetElementType();
		unsigned int NUME = FEMData->GetEleGrpList()[EleGrp].GetNUME();

		*this << " ELEMENT TYPE  . . . . . . . . . . . . .( NPAR(1) ) . . =" << setw(5)
			  << ElementType << endl;
		*this << "     EQ.1, TRUSS ELEMENTS" << endl
			  << "     EQ.2, ELEMENTS CURRENTLY" << endl
			  << "     EQ.3, NOT AVAILABLE" << endl
			  << endl;

		*this << " NUMBER OF ELEMENTS. . . . . . . . . . .( NPAR(2) ) . . =" << setw(5) << NUME
			  << endl
			  << endl;

		switch (ElementType)
		{
			case ElementTypes::Bar: // Bar element
				OutputBarElements(EleGrp);
				break;
			case ElementTypes::C3D8R: // C3D8R element
				OutputC3D8RElements(EleGrp);
				break;
			case ElementTypes::S4R: // S4R element
				*this << " I don't think it is necesasry to output the element info for S4R." << endl << "You can use utils/visualize_s4r_patch.py to visualize the S4R element." << endl << endl;
			case ElementTypes::B31: // B31 element
				OutputB31Elements(EleGrp);
				break;
		    default:
		        *this << ElementType << " has not been implemented yet." << endl;
		        break;
		}
	}
}
//	Output C3D8R element data
void COutputter::OutputC3D8RElements(unsigned int EleGrp)
{
	CDomain* FEMData = CDomain::GetInstance();

	CElementGroup& ElementGroup = FEMData->GetEleGrpList()[EleGrp];
	unsigned int NUMMAT = ElementGroup.GetNUMMAT();

	*this << " M A T E R I A L   D E F I N I T I O N" << endl
		  << endl;
	*this << " NUMBER OF DIFFERENT SETS OF MATERIAL" << endl;
	*this << " AND CONSTANTS  . . . . . . . . . . .( NPAR(3) ) . . =" << setw(5) << NUMMAT
		  << endl
		  << endl;

	*this << "  SET       YOUNG'S     POISSON'S" << endl
		  << " NUMBER     MODULUS        RATIO" << endl
		  << "               E            nu" << endl;

	*this << setiosflags(ios::scientific) << setprecision(5);

	//	Loop over for all property sets
	for (unsigned int mset = 0; mset < NUMMAT; mset++)
    {
        *this << setw(5) << mset+1;
		ElementGroup.GetMaterial(mset).Write(*this);
    }

	*this << endl << endl
		  << " E L E M E N T   I N F O R M A T I O N" << endl;
    
	*this << " ELEMENT     NODE     NODE     NODE     NODE     NODE     NODE     NODE     NODE       MATERIAL" << endl
		  << " NUMBER-N      I        J        K        L        M        N        O        P       SET NUMBER" << endl;

	unsigned int NUME = ElementGroup.GetNUME();

	//	Loop over for all elements in group EleGrp
	for (unsigned int Ele = 0; Ele < NUME; Ele++)
    {
        *this << setw(5) << Ele+1;
		ElementGroup[Ele].Write(*this);
    }

	*this << endl;
}

//	Output bar element data
void COutputter::OutputBarElements(unsigned int EleGrp)
{
	CDomain* FEMData = CDomain::GetInstance();

	CElementGroup& ElementGroup = FEMData->GetEleGrpList()[EleGrp];
	unsigned int NUMMAT = ElementGroup.GetNUMMAT();

	*this << " M A T E R I A L   D E F I N I T I O N" << endl
		  << endl;
	*this << " NUMBER OF DIFFERENT SETS OF MATERIAL" << endl;
	*this << " AND CROSS-SECTIONAL  CONSTANTS  . . . .( NPAR(3) ) . . =" << setw(5) << NUMMAT
		  << endl
		  << endl;

	*this << "  SET       YOUNG'S     CROSS-SECTIONAL" << endl
		  << " NUMBER     MODULUS          AREA" << endl
		  << "               E              A" << endl;

	*this << setiosflags(ios::scientific) << setprecision(5);

	//	Loop over for all property sets
	for (unsigned int mset = 0; mset < NUMMAT; mset++)
    {
        *this << setw(5) << mset+1;
		ElementGroup.GetMaterial(mset).Write(*this);
    }

	*this << endl << endl
		  << " E L E M E N T   I N F O R M A T I O N" << endl;
    
	*this << " ELEMENT     NODE     NODE       MATERIAL" << endl
		  << " NUMBER-N      I        J       SET NUMBER" << endl;

	unsigned int NUME = ElementGroup.GetNUME();

	//	Loop over for all elements in group EleGrp
	for (unsigned int Ele = 0; Ele < NUME; Ele++)
    {
        *this << setw(5) << Ele+1;
		ElementGroup[Ele].Write(*this);
    }

	*this << endl;
}

// 输出B31单元数据
void COutputter::OutputB31Elements(unsigned int EleGrp)
{
	CDomain* FEMData = CDomain::GetInstance();
	CElementGroup& ElementGroup = FEMData->GetEleGrpList()[EleGrp];
	unsigned int NUMMAT = ElementGroup.GetNUMMAT();

	*this << " M A T E R I A L   D E F I N I T I O N (B31)" << endl << endl;
	*this << " NUMBER OF DIFFERENT SETS OF MATERIAL" << endl;
	*this << " AND SECTION CONSTANTS  . . . .( NPAR(3) ) . . =" << setw(5) << NUMMAT << endl << endl;
	*this << "  SET       E           G           A           Iy          Iz          J" << endl;
	*this << " NUMBER     (Pa)        (Pa)        (m2)        (m4)        (m4)        (m4)" << endl;
	*this << setiosflags(ios::scientific) << setprecision(5);
	for (unsigned int mset = 0; mset < NUMMAT; mset++)
	{
		*this << setw(5) << mset+1;
		ElementGroup.GetMaterial(mset).Write(*this); // 假定CB31Material::Write已实现
	}
	*this << endl << endl << " E L E M E N T   I N F O R M A T I O N (B31)" << endl;
	*this << " ELEMENT     NODE     NODE       MATERIAL" << endl
		  << " NUMBER-N      I        J       SET NUMBER" << endl;
	unsigned int NUME = ElementGroup.GetNUME();
	for (unsigned int Ele = 0; Ele < NUME; Ele++)
	{
		*this << setw(5) << Ele+1;
		ElementGroup[Ele].Write(*this); // 假定B31::Write已实现
	}
	*this << endl;
}

//	Print load data
void COutputter::OutputLoadInfo()
{
	CDomain* FEMData = CDomain::GetInstance();

	for (unsigned int lcase = 1; lcase <= FEMData->GetNLCASE(); lcase++)
	{
		CLoadCaseData* LoadData = &FEMData->GetLoadCases()[lcase - 1];

		*this << setiosflags(ios::scientific);
		*this << " L O A D   C A S E   D A T A" << endl
			  << endl;

		*this << "     LOAD CASE NUMBER . . . . . . . =" << setw(6) << lcase << endl;
		*this << "     NUMBER OF CONCENTRATED LOADS . =" << setw(6) << LoadData->nloads << endl
			  << endl;
		*this << "    NODE       DIRECTION      LOAD" << endl
			  << "   NUMBER                   MAGNITUDE" << endl;

		LoadData->Write(*this);

		*this << endl;
	}
}

//	Print nodal displacement
void COutputter::OutputNodalDisplacement()
{
    CDomain* FEMData = CDomain::GetInstance();
    CNode* NodeList = FEMData->GetNodeList();
    double* Displacement = FEMData->GetDisplacement();
    *this << setiosflags(ios::scientific);

    *this << " D I S P L A C E M E N T S" << endl << endl;
    *this << "  NODE           X-DISPLACEMENT    Y-DISPLACEMENT    Z-DISPLACEMENT    THETA-X    THETA-Y    THETA-Z" << endl;
    for (unsigned int np = 0; np < FEMData->GetNUMNP(); np++) {
		NodeList[np].WriteNodalDisplacement(*this, Displacement);
	}
    // 输出有限元解与精确解的L2误差
    OutputL2Error();
}

// 计算精确解w（板中心线挠度），可根据需要扩展为Mx等
// 这里只实现w的精确解，参数含义与Python一致
#include <vector>
#include <cmath>

double ExactSolutionW(double x, double y, double D, double q, double a, double b, double nu)
{
    // 只实现w的精确解，简化实现
    double pi = M_PI;
    double K = -4 * q * a * a / pow(pi, 3);
    int m_all[4] = {1, 3, 5, 7};
    double E[8] = {0};
    E[1] = 0.3722 * K;
    E[3] = -0.0380 * K;
    E[5] = -0.0178 * K;
    E[7] = -0.0085 * K;
    double w1 = 0.0, w2 = 0.0, w3 = 0.0;
    for (int mi = 0; mi < 4; ++mi) {
        int m = m_all[mi];
        double a_m = m * pi * b / (2 * a);
        double b_m = m * pi * a / (2 * b);
        double A1 = 4 * q * pow(a, 4) / (pow(pi, 5) * D) * pow(-1, (m-1)/2) / pow(m, 5);
        double B1 = (a_m * tanh(a_m) + 2) / (2 * cosh(a_m));
        double C1 = 1 / (2 * cosh(a_m));
        double D1 = m * pi / a;
        double A2 = -a * a / (2 * pi * pi * D) * E[m] * pow(-1, (m-1)/2) / (m * m * cosh(a_m));
        double B2 = a_m * tanh(a_m);
        double D2 = m * pi / a;
        double A3 = -b * b / (2 * pi * pi * D) * E[m] * pow(-1, (m-1)/2) / (m * m * cosh(b_m));
        double B3 = b_m * tanh(b_m);
        double D3 = m * pi / b;
        w1 += A1 * cos(D1 * x) * (1 - B1 * cosh(D1 * y) + C1 * D1 * y * sinh(D1 * y));
        w2 += A2 * cos(D2 * x) * (D2 * y * sinh(D2 * y) - B2 * cosh(D2 * y));
        w3 += A3 * cos(D3 * y) * (D3 * x * sinh(D3 * x) - B3 * cosh(D3 * x));
    }
    return w1 + w2 + w3;
}

// 输出有限元解与精确解的L2误差
void COutputter::OutputL2Error()
{
    CDomain* FEMData = CDomain::GetInstance();
    CNode* NodeList = FEMData->GetNodeList();
    double* Displacement = FEMData->GetDisplacement();
    unsigned int NUMNP = FEMData->GetNUMNP();
    double D = 1.0, nu = 0, E = 0, thickness = 0, q = 0, a = 0, b = 0;
    bool found = false;
    unsigned int NUMEG = FEMData->GetNUMEG();
    for (unsigned int eg = 0; eg < NUMEG; ++eg) {
        CElementGroup& EleGrp = FEMData->GetEleGrpList()[eg];
        if (EleGrp.GetElementType() == ElementTypes::S4R) {
            CS4RMaterial* mat = dynamic_cast<CS4RMaterial*>(&EleGrp.GetMaterial(0));
            if (mat) {
                E = mat->E;
                nu = mat->nu;
                thickness = mat->thickness;
                found = true;
                break;
            }
        }
    }
    if (!found) {
        *this << "[L2误差] 未找到S4R材料参数，无法计算精确解。" << endl;
        return;
    }
    // 板刚度D
    D = E * pow(thickness, 3) / (12.0 * (1.0 - nu * nu));
    // 板长宽a,b（节点x/y最大最小值）
    double x_min = 1e20, x_max = -1e20, y_min = 1e20, y_max = -1e20;
    for (unsigned int np = 0; np < NUMNP; ++np) {
        double x = NodeList[np].XYZ[0];
        double y = NodeList[np].XYZ[1];
        if (x < x_min) x_min = x;
        if (x > x_max) x_max = x;
        if (y < y_min) y_min = y;
        if (y > y_max) y_max = y;
    }
    a = x_max - x_min;
    b = y_max - y_min;
    // 均布载荷q（取第一个荷载工况所有z向集中荷载之和/面积）
    q = 0.0;
    if (FEMData->GetNLCASE() > 0) {
        CLoadCaseData* LoadData = &FEMData->GetLoadCases()[0];
        for (unsigned int i = 0; i < LoadData->nloads; ++i) {
            // dof==3 代表Z方向
            if (LoadData->dof[i] == 3) {
                q += LoadData->load[i];
            }
        }
        q /= (a * b);
    }
    // 计算L2误差（节点坐标中心化）
    double l2_sum = 0.0;
    double l2_exact = 0.0;
    for (unsigned int np = 0; np < NUMNP; ++np) {
        double x = NodeList[np].XYZ[0] - a/2.0 - x_min;
        double y = NodeList[np].XYZ[1] - b/2.0 - y_min;
        unsigned int eqn = NodeList[np].bcode[2];
        double w_fem = (eqn ? Displacement[eqn - 1] : 0.0);
        double w_exact = ExactSolutionW(x, y, D, q, a, b, nu);
        l2_sum += (w_fem - w_exact) * (w_fem - w_exact);
        l2_exact += w_exact * w_exact;
    }
    double l2_error = sqrt(l2_sum / NUMNP);
    double l2_rel = sqrt(l2_sum / l2_exact);
    *this << "L2 绝对误差 (w): " << l2_error << endl;
    *this << "L2 相对误差 (w): " << l2_rel << endl << endl;
}

//	Calculate stresses
void COutputter::OutputElementStress()
{
    CDomain* FEMData = CDomain::GetInstance();

    double* Displacement = FEMData->GetDisplacement();

    unsigned int NUMEG = FEMData->GetNUMEG();

    for (unsigned int EleGrpIndex = 0; EleGrpIndex < NUMEG; EleGrpIndex++)
    {
        *this << " S T R E S S  C A L C U L A T I O N S  F O R  E L E M E N T  G R O U P" << setw(5)
              << EleGrpIndex + 1 << endl
              << endl;

        CElementGroup& EleGrp = FEMData->GetEleGrpList()[EleGrpIndex];
        unsigned int NUME = EleGrp.GetNUME();
        ElementTypes ElementType = EleGrp.GetElementType();

        switch (ElementType)
        {
            case ElementTypes::Bar: // Bar element
                *this << "  ELEMENT             FORCE            STRESS" << endl
                    << "  NUMBER" << endl;

                double stress;

                for (unsigned int Ele = 0; Ele < NUME; Ele++)
                {
                    CElement& Element = EleGrp[Ele];
                    Element.ElementStress(&stress, Displacement);

                    CBarMaterial& material = *dynamic_cast<CBarMaterial*>(Element.GetElementMaterial());
                    *this << setw(5) << Ele + 1 << setw(22) << stress * material.Area << setw(18)
                        << stress << endl;
                }

                *this << endl;
				break;

			case ElementTypes::C3D8R: // C3D8R element
				*this << "  ELEMENT      SIGMA_X      SIGMA_Y      SIGMA_Z      TAU_XY      TAU_YZ      TAU_ZX     VON MISES" << endl
					<< "  NUMBER" << endl;

				double stresses[6];
				for (unsigned int Ele = 0; Ele < NUME; Ele++)
				{
					CElement& Element = EleGrp[Ele];
					Element.ElementStress(stresses, Displacement);

					// Calculate von Mises stress
					double s1 = stresses[0];
					double s2 = stresses[1];
					double s3 = stresses[2];
					double tau12 = stresses[3];
					double tau23 = stresses[4];
					double tau31 = stresses[5];
					double vonMises = sqrt(0.5 * ((s1-s2)*(s1-s2) + (s2-s3)*(s2-s3) + (s3-s1)*(s3-s1) + 6.0 * (tau12*tau12 + tau23*tau23 + tau31*tau31)));

					*this << setw(5) << Ele + 1 
						<< setw(12) << stresses[0] << setw(12) << stresses[1] << setw(12) << stresses[2]
						<< setw(12) << stresses[3] << setw(12) << stresses[4] << setw(12) << stresses[5]
						<< setw(12) << vonMises << endl;
				}

				*this << endl;
                break;
            case ElementTypes::S4R:
                *this << "  ELEMENT     Mx           My           Mxy          Qx           Qy" << endl;
                *this << "  NUMBER" << endl;
                for (unsigned int Ele = 0; Ele < NUME; Ele++) {
                    CElement& Element = EleGrp[Ele];
                    double stress[5];
                    Element.ElementStress(stress, Displacement);
                    *this << setw(5) << Ele + 1;
                    for (int i = 0; i < 5; ++i) {
                        *this << setw(13) << stress[i];
                    }
                    *this << endl;
                }
                *this << endl;
            	break;
			case ElementTypes::B31: // B31 element
				*this << "  ELEMENT     N-FORCE      Mx         My         TORSION      SIGMA" << endl;
				for (unsigned int Ele = 0; Ele < NUME; Ele++)
				{
					double stress[6] = {0}; 
					CElement& Element = EleGrp[Ele];
					Element.ElementStress(stress, Displacement); 
					*this << setw(5) << Ele+1;
					for(int i=0;i<6;++i) *this << setw(12) << stress[i];
					*this << endl;
				}
				*this << endl;
				break;
			default: // Invalid element type
				cerr << "*** Error *** Elment type " << ElementType
					<< " has not been implemented.\n\n";
        }
	}    
}

//	Print total system data
void COutputter::OutputTotalSystemData()
{
	CDomain* FEMData = CDomain::GetInstance();

	*this << "	TOTAL SYSTEM DATA" << endl
		  << endl;

	*this << "     NUMBER OF EQUATIONS . . . . . . . . . . . . . .(NEQ) = " << FEMData->GetNEQ()
		  << endl
		  << "     NUMBER OF MATRIX ELEMENTS . . . . . . . . . . .(NWK) = " << FEMData->GetStiffnessMatrix()->size()
		  << endl
		  << "     MAXIMUM HALF BANDWIDTH  . . . . . . . . . . . .(MK ) = " << FEMData->GetStiffnessMatrix()->GetMaximumHalfBandwidth()
		  << endl
		  << "     MEAN HALF BANDWIDTH . . . . . . . . . . . . . .(MM ) = " << FEMData->GetStiffnessMatrix()->size() / FEMData->GetNEQ() << endl
		  << endl
		  << endl;
}

#ifdef _DEBUG_

//	Print column heights for debuging
void COutputter::PrintColumnHeights()
{
	*this << "*** _Debug_ *** Column Heights" << endl;

	CDomain* FEMData = CDomain::GetInstance();

	unsigned int NEQ = FEMData->GetNEQ();
	CSkylineMatrix<double> *StiffnessMatrix = FEMData->GetStiffnessMatrix();
	unsigned int* ColumnHeights = StiffnessMatrix->GetColumnHeights();

	for (unsigned int col = 0; col < NEQ; col++)
	{
		if (col + 1 % 10 == 0)
		{
			*this << endl;
		}

		*this << setw(8) << ColumnHeights[col];
	}

	*this << endl
		  << endl;
}

//	Print address of diagonal elements for debuging
void COutputter::PrintDiagonalAddress()
{
	*this << "*** _Debug_ *** Address of Diagonal Element" << endl;

	CDomain* FEMData = CDomain::GetInstance();

	unsigned int NEQ = FEMData->GetNEQ();
	CSkylineMatrix<double> *StiffnessMatrix = FEMData->GetStiffnessMatrix();
	unsigned int* DiagonalAddress = StiffnessMatrix->GetDiagonalAddress();

	for (unsigned int col = 0; col <= NEQ; col++)
	{
		if (col + 1 % 10 == 0)
		{
			*this << endl;
		}

		*this << setw(8) << DiagonalAddress[col];
	}

	*this << endl
		  << endl;
}

//	Print banded and full stiffness matrix for debuging
void COutputter::PrintStiffnessMatrix()
{
	*this << "*** _Debug_ *** Banded stiffness matrix" << endl;

	CDomain* FEMData = CDomain::GetInstance();

	unsigned int NEQ = FEMData->GetNEQ();
	CSkylineMatrix<double> *StiffnessMatrix = FEMData->GetStiffnessMatrix();
	unsigned int* DiagonalAddress = StiffnessMatrix->GetDiagonalAddress();

	*this << setiosflags(ios::scientific) << setprecision(5);

	for (unsigned int i = 0; i < DiagonalAddress[NEQ] - DiagonalAddress[0]; i++)
	{
		*this << setw(14) << (*StiffnessMatrix)(i);

		if ((i + 1) % 6 == 0)
		{
			*this << endl;
		}
	}

	*this << endl
		  << endl;

	*this << "*** _Debug_ *** Full stiffness matrix" << endl;

	for (int I = 1; I <= NEQ; I++)
	{
		for (int J = 1; J <= NEQ; J++)
		{
			int J_new = (J > I) ? J : I;
			int I_new = (J > I) ? I : J;
			int H = DiagonalAddress[J_new] - DiagonalAddress[J_new - 1];
			if (J_new - I_new - H >= 0)
			{
				*this << setw(14) << 0.0;
			}
			else
			{
				*this << setw(14) << (*StiffnessMatrix)(I_new, J_new);
			}
		}

		*this << endl;
	}

	*this << endl;
}

//	Print displacement vector for debuging
void COutputter::PrintDisplacement()
{
	*this << "*** _Debug_ *** Displacement vector" << endl;

	CDomain* FEMData = CDomain::GetInstance();

	unsigned int NEQ = FEMData->GetNEQ();
	double* Force = FEMData->GetForce();

	*this << setiosflags(ios::scientific) << setprecision(5);

	for (unsigned int i = 0; i < NEQ; i++)
	{
		if ((i + 1) % 6 == 0)
		{
			*this << endl;
		}

		*this << setw(14) << Force[i];
	}

	*this << endl
		  << endl;
}

#endif
#include <fstream>
#include <vector>
#include <string>

void COutputter::OutputVTK(const std::string& filename) {
    CDomain* FEMData = CDomain::GetInstance();
    unsigned int NUMNP = FEMData->GetNUMNP();
    unsigned int NUMEG = FEMData->GetNUMEG();
    CNode* NodeList = FEMData->GetNodeList();
    double* Displacement = FEMData->GetDisplacement();

    std::ofstream vtkfile(filename);
    if (!vtkfile) return;

    // VTK头部
    vtkfile << "# vtk DataFile Version 3.0\n";
    vtkfile << "STAP++ Results\n";
    vtkfile << "ASCII\n";
    vtkfile << "DATASET UNSTRUCTURED_GRID\n";

    // 节点坐标
    vtkfile << "POINTS " << NUMNP << " float\n";
    for (unsigned int i = 0; i < NUMNP; ++i) {
        vtkfile << NodeList[i].XYZ[0] << " " << NodeList[i].XYZ[1] << " " << NodeList[i].XYZ[2] << "\n";
    }

    // 统计单元总数和单元节点总数
    unsigned int total_cells = 0, total_cell_size = 0;
    std::vector<std::pair<int, std::vector<unsigned int>>> cell_info; // {VTK类型, 节点编号}
    for (unsigned int eg = 0; eg < NUMEG; ++eg) {
        CElementGroup& EleGrp = FEMData->GetEleGrpList()[eg];
        unsigned int NUME = EleGrp.GetNUME();
        ElementTypes etype = EleGrp.GetElementType();
        for (unsigned int e = 0; e < NUME; ++e) {
            std::vector<unsigned int> conn;
            int vtk_type = 0;
            CElement& elem = EleGrp[e];
            CNode** elem_nodes = elem.GetNodes();
            if (etype == ElementTypes::Bar) {
                // Bar: 2节点，VTK_LINE=3
                conn.push_back(elem_nodes[0]->NodeNumber - 1);
                conn.push_back(elem_nodes[1]->NodeNumber - 1);
                vtk_type = 3;
            } else if (etype == ElementTypes::S4R) {
                // S4R: 4节点，VTK_QUAD=9
                for (int k = 0; k < 4; ++k)
                    conn.push_back(elem_nodes[k]->NodeNumber - 1);
                vtk_type = 9;
            } else if (etype == ElementTypes::C3D8R) {
                // C3D8R: 8节点，VTK_HEXAHEDRON=12
                for (int k = 0; k < 8; ++k)
                    conn.push_back(elem_nodes[k]->NodeNumber - 1);
                vtk_type = 12;
            }
            if (!conn.empty()) {
                cell_info.push_back({vtk_type, conn});
                total_cells++;
                total_cell_size += (1 + conn.size());
            }
        }
    }
    // 单元连接
    vtkfile << "CELLS " << total_cells << " " << total_cell_size << "\n";
    for (const auto& c : cell_info) {
        vtkfile << c.second.size();
        for (auto nid : c.second) vtkfile << " " << nid;
        vtkfile << "\n";
    }
    // 单元类型
    vtkfile << "CELL_TYPES " << total_cells << "\n";
    for (const auto& c : cell_info) vtkfile << c.first << "\n";

    // 节点数据：位移
    vtkfile << "POINT_DATA " << NUMNP << "\n";
    vtkfile << "VECTORS Displacement float\n";
    for (unsigned int i = 0; i < NUMNP; ++i) {
        double ux = 0, uy = 0, uz = 0;
        if (NodeList[i].bcode[0]) ux = Displacement[NodeList[i].bcode[0] - 1];
        if (NodeList[i].bcode[1]) uy = Displacement[NodeList[i].bcode[1] - 1];
        if (NodeList[i].bcode[2]) uz = Displacement[NodeList[i].bcode[2] - 1];
        vtkfile << ux << " " << uy << " " << uz << "\n";
    }

    // 单元数据：应力
    vtkfile << "CELL_DATA " << total_cells << "\n";
    vtkfile << "SCALARS Stress float 1\nLOOKUP_TABLE default\n";
    unsigned int cell_idx = 0;
    for (unsigned int eg = 0; eg < NUMEG; ++eg) {
        CElementGroup& EleGrp = FEMData->GetEleGrpList()[eg];
        unsigned int NUME = EleGrp.GetNUME();
        ElementTypes etype = EleGrp.GetElementType();
        for (unsigned int e = 0; e < NUME; ++e) {
            double val = 0;
            if (etype == ElementTypes::Bar) {
                double stress = 0;
                EleGrp[e].ElementStress(&stress, Displacement);
                val = stress;
            } else if (etype == ElementTypes::C3D8R) {
                double stress[6];
                EleGrp[e].ElementStress(stress, Displacement);
                // 取等效应力（von Mises）
                double s1 = stress[0], s2 = stress[1], s3 = stress[2];
                double tau12 = stress[3], tau23 = stress[4], tau31 = stress[5];
                val = sqrt(0.5 * ((s1-s2)*(s1-s2) + (s2-s3)*(s2-s3) + (s3-s1)*(s3-s1) + 6.0 * (tau12*tau12 + tau23*tau23 + tau31*tau31)));
            } else if (etype == ElementTypes::S4R) {
                double stress[5];
                EleGrp[e].ElementStress(stress, Displacement);
                // 取Mx为主（可根据需求改为其它分量）
                val = stress[0];
            }
            vtkfile << val << "\n";
            cell_idx++;
        }
    }
    vtkfile.close();
}
