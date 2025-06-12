/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Domain.h"
#include "Material.h"
#include "Element.h"

using namespace std;


CDomain* CDomain::_instance = nullptr;

//	Constructor
CDomain::CDomain()
{
	Title[0] = '0';
	MODEX = 0;

	NUMNP = 0;
	NodeList = nullptr;
	
	NUMEG = 0;
	EleGrpList = nullptr;
	
	NLCASE = 0;
	NLOAD = nullptr;
	LoadCases = nullptr;
	
	NEQ = 0;

	Force = nullptr;
	StiffnessMatrix = nullptr;
}

//	Desconstructor
CDomain::~CDomain()
{
	delete [] NodeList;

	delete [] EleGrpList;

	delete [] NLOAD;
	delete [] LoadCases;

	delete [] Force;
	delete StiffnessMatrix;
}

//	Return pointer to the instance of the Domain class
CDomain* CDomain::GetInstance()
{
	if (!_instance) 
		_instance = new CDomain();
	
	return _instance;
}

//	Read domain data from the input data file
bool CDomain::ReadData(string FileName, string OutFile)
{
	Input.open(FileName);

	if (!Input) 
	{
		cerr << "*** Error *** File " << FileName << " does not exist !" << endl;
		exit(3);
	}

	COutputter* Output = COutputter::GetInstance(OutFile);

//	Read the heading line
	Input.getline(Title, 256);
	Output->OutputHeading();

//	Read the control line
	Input >> NUMNP >> NUMEG >> NLCASE >> MODEX;

//	Read nodal point data
	if (ReadNodalPoints())
        Output->OutputNodeInfo();
    else
        return false;

//	Update equation number
	CalculateEquationNumber();
	Output->OutputEquationNumber();

//	Read load data
	if (ReadLoadCases())
        Output->OutputLoadInfo();
    else
        return false;

//	Read element data
	if (ReadElements())
        Output->OutputElementInfo();
    else {
        Output->OutputElementInfo();
		return false;
	}

	return true;
}

//	Read nodal point data
bool CDomain::ReadNodalPoints()
{

//	Read nodal point data lines
	NodeList = new CNode[NUMNP];

//	Loop over for all nodal points
	for (unsigned int np = 0; np < NUMNP; np++)
    {
		if (!NodeList[np].Read(Input))
			return false;
    
        if (NodeList[np].NodeNumber != np + 1)
        {
            cerr << "*** Error *** Nodes must be inputted in order !" << endl
            << "   Expected node number : " << np + 1 << endl
            << "   Provided node number : " << NodeList[np].NodeNumber << endl;
        
            return false;
        }
    }

	return true;
}

//	Calculate global equation numbers corresponding to every degree of freedom of each node
void CDomain::CalculateEquationNumber()
{
	NEQ = 0;
	for (unsigned int np = 0; np < NUMNP; np++)	// Loop over for all node
	{
		for (unsigned int dof = 0; dof < CNode::NDF; dof++)	// Loop over for DOFs of node np
		{
			if (NodeList[np].bcode[dof]) 
				NodeList[np].bcode[dof] = 0;
			else
			{
				NEQ++;
				NodeList[np].bcode[dof] = NEQ;
			}
		}
	}
}

//	Read load case data
bool CDomain::ReadLoadCases()
{
//	Read load data lines
	LoadCases = new CLoadCaseData[NLCASE];	// List all load cases

//	Loop over for all load cases
	for (unsigned int lcase = 0; lcase < NLCASE; lcase++)
    {
        unsigned int LL;
        Input >> LL;
        
        if (LL != lcase + 1)
        {
            cerr << "*** Error *** Load case must be inputted in order !" << endl
            << "   Expected load case : " << lcase + 1 << endl
            << "   Provided load case : " << LL << endl;
            
            return false;
        }

        LoadCases[lcase].Read(Input);
    }

	return true;
}

// Read element data
bool CDomain::ReadElements()
{
    EleGrpList = new CElementGroup[NUMEG];

//	Loop over for all element group
	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)
        if (!EleGrpList[EleGrp].Read(Input))
            return false;
    
    return true;
}

//	Calculate column heights
void CDomain::CalculateColumnHeights()
{
#ifdef _DEBUG_
    COutputter* Output = COutputter::GetInstance();
    *Output << setw(9) << "Ele = " << setw(22) << "Location Matrix" << endl;
#endif
	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)		//	Loop over for all element groups
    {
        CElementGroup& ElementGrp = EleGrpList[EleGrp];
        unsigned int NUME = ElementGrp.GetNUME();
		for (unsigned int Ele = 0; Ele < NUME; Ele++)	//	Loop over for all elements in group EleGrp
        {
            CElement& Element = ElementGrp[Ele];
            // Generate location matrix
            Element.GenerateLocationMatrix();
#ifdef _DEBUG_
            unsigned int* LocationMatrix = Element.GetLocationMatrix();
            
            *Output << setw(9) << Ele+1;
            for (int i=0; i<Element.GetND(); i++)
                *Output << setw(5) << LocationMatrix[i];
            *Output << endl;
#endif

            StiffnessMatrix->CalculateColumnHeight(Element.GetLocationMatrix(), Element.GetND());
        }
    }
    
    StiffnessMatrix->CalculateMaximumHalfBandwidth();
    
#ifdef _DEBUG_
    *Output << endl;
	Output->PrintColumnHeights();
#endif

}

//    Allocate storage for matrices Force, ColumnHeights, DiagonalAddress and StiffnessMatrix
//    and calculate the column heights and address of diagonal elements
void CDomain::AllocateMatrices()
{
    //    Allocate for global force/displacement vector
    Force = new double[NEQ];
    
    //  Create the banded stiffness matrix
    StiffnessMatrix = new CSkylineMatrix<double>(NEQ);
    
    //    Calculate column heights
    CalculateColumnHeights();
    
    //    Calculate address of diagonal elements in banded matrix
    StiffnessMatrix->CalculateDiagnoalAddress();
    
    //    Allocate for banded global stiffness matrix
    StiffnessMatrix->Allocate();
    
    COutputter* Output = COutputter::GetInstance();
    Output->OutputTotalSystemData();
}

//	Assemble the banded gloabl stiffness matrix
void CDomain::AssembleStiffnessMatrix()
{
//	Loop over for all element groups
	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)
	{
        CElementGroup& ElementGrp = EleGrpList[EleGrp];
        unsigned int NUME = ElementGrp.GetNUME();

		unsigned int size = ElementGrp[0].SizeOfStiffnessMatrix();
		double* Matrix = new double[size];

//		Loop over for all elements in group EleGrp
		for (unsigned int Ele = 0; Ele < NUME; Ele++)
        {
            CElement& Element = ElementGrp[Ele];
            Element.ElementStiffness(Matrix);
            StiffnessMatrix->Assembly(Matrix, Element.GetLocationMatrix(), Element.GetND());
        }

		delete[] Matrix;
		Matrix = nullptr;
	}

#ifdef _DEBUG_
	COutputter* Output = COutputter::GetInstance();
	Output->PrintStiffnessMatrix();
#endif

}

//	Assemble the global nodal force vector for load case LoadCase
bool CDomain::AssembleForce(unsigned int LoadCase)
{
    if (LoadCase > NLCASE) 
        return false;

    CLoadCaseData* LoadData = &LoadCases[LoadCase - 1];

    clear(Force, NEQ);

    // 1. 集中荷载
    for (unsigned int lnum = 0; lnum < LoadData->nloads; lnum++)
    {
        unsigned int dof = NodeList[LoadData->node[lnum] - 1].bcode[LoadData->dof[lnum] - 1];
        if(dof) // The DOF is activated
            Force[dof - 1] += LoadData->load[lnum];
    }

    // 2. 单元重力体力
    const double g = -10.0; // 重力加速度
    for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)
    {
        CElementGroup& ElementGrp = EleGrpList[EleGrp];
        unsigned int NUME = ElementGrp.GetNUME();
        for (unsigned int Ele = 0; Ele < NUME; Ele++)
        {
            CElement& Element = ElementGrp[Ele];
            CMaterial* mat = Element.GetElementMaterial();
            double density = 0.0, area = 1.0, thickness = 1.0;
            // 针对不同单元类型提取参数
            if (auto* matS4R = dynamic_cast<CS4RMaterial*>(mat)) {
                density = matS4R->density;
                thickness = matS4R->thickness;
            } else if (auto* matBar = dynamic_cast<CBarMaterial*>(mat)) {
                density = matBar->density;
                area = matBar->Area;
            } else if (auto* matB31 = dynamic_cast<CB31Material*>(mat)) {
                density = matB31->density;
                area = matB31->Area;
            } else if (auto* matC3D8R = dynamic_cast<CC3D8RMaterial*>(mat)) {
                density = matC3D8R->density;
            }
            // 节点数和自由度
            unsigned int NEN = Element.GetNEN();
            CNode** nodes = Element.GetNodes();
            unsigned int* LocationMatrix = Element.GetLocationMatrix();
            unsigned int ND = Element.GetND();

            // 计算单元重力节点力（以均匀分布为例）
            // S4R: 面密度*厚度*面积/4, B31/Bar: 线密度*长度/2, C3D8R: 体密度*体积/8
            if (auto* matS4R = dynamic_cast<CS4RMaterial*>(mat)) {
                // S4R壳单元
                // 计算单元面积
                double x[4], y[4];
                for (int i = 0; i < 4; ++i) {
                    x[i] = nodes[i]->XYZ[0];
                    y[i] = nodes[i]->XYZ[1];
                }
                double area4 = 0.5 * fabs(
                    (x[0]*y[1] + x[1]*y[2] + x[2]*y[3] + x[3]*y[0])
                  - (y[0]*x[1] + y[1]*x[2] + y[2]*x[3] + y[3]*x[0])
                );
                double fz = density * thickness * area4 * g / 4.0;
                for (int i = 0; i < 4; ++i) {
                    int loc = LocationMatrix[i*6+2]; // z方向自由度
                    if (loc) Force[loc-1] += fz;
                }
            } else if (auto* matBar = dynamic_cast<CBarMaterial*>(mat)) {
                // Bar单元
                double x1 = nodes[0]->XYZ[0], y1 = nodes[0]->XYZ[1], z1 = nodes[0]->XYZ[2];
                double x2 = nodes[1]->XYZ[0], y2 = nodes[1]->XYZ[1], z2 = nodes[1]->XYZ[2];
                double L = sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1));
                double fz = density * area * L * g / 2.0;
                for (int i = 0; i < 2; ++i) {
                    int loc = LocationMatrix[i*6+2]; // z方向
                    if (loc) Force[loc-1] += fz;
                }
            } else if (auto* matB31 = dynamic_cast<CB31Material*>(mat)) {
                // B31梁单元
                double x1 = nodes[0]->XYZ[0], y1 = nodes[0]->XYZ[1], z1 = nodes[0]->XYZ[2];
                double x2 = nodes[1]->XYZ[0], y2 = nodes[1]->XYZ[1], z2 = nodes[1]->XYZ[2];
                double L = sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1));
                double fz = density * area * L * g / 2.0;
                for (int i = 0; i < 2; ++i) {
                    int loc = LocationMatrix[i*6+2]; // z方向
                    if (loc) Force[loc-1] += fz;
                }
            } else if (auto* matC3D8R = dynamic_cast<CC3D8RMaterial*>(mat)) {
                // C3D8R体单元
                double x[8], y[8], z[8];
                for (int i = 0; i < 8; ++i) {
                    x[i] = nodes[i]->XYZ[0];
                    y[i] = nodes[i]->XYZ[1];
                    z[i] = nodes[i]->XYZ[2];
                }
                // 体积近似为六面体对角线分割法
                double V = fabs(
                    (x[0]*y[1]*z[2] + x[1]*y[2]*z[3] + x[2]*y[3]*z[0] + x[3]*y[0]*z[1])
                  - (z[0]*y[1]*x[2] + z[1]*y[2]*x[3] + z[2]*y[3]*x[0] + z[3]*y[0]*x[1])
                ); // 可用更精确方法
                double fz = density * V * g / 8.0;
                for (int i = 0; i < 8; ++i) {
                    int loc = LocationMatrix[i*6+2]; // z方向
                    if (loc) Force[loc-1] += fz;
                }
            }
        }
    }

    return true;
}

