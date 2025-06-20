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
#include "Bar.h"
#include "Outputter.h"
#include "Clock.h"
#include "EigenSolver.h"
#include <iostream>

using namespace std;

int main(int argc, char *argv[])
{
    cout << "**********************************************************" << endl;
    cout << "Program Launch: STAP++" << endl;
	if (argc < 2 || argc > 3) //  Print help message
	{
	    cout << "Usage: stap++ InputFileName [SolverType]\n";
        cout << "SolverType: 0 - Skyline LDLT (Original)\n";
        cout << "            1 - Eigen Direct Solver (SimplicialLDLT)\n";
        cout << "            2 - Eigen Iterative Solver (Conjugate Gradient)\n";
		exit(1);
	}

	string filename(argv[1]);
    size_t found = filename.find_last_of('.');

    // If the input file name is provided with an extension
    if (found != std::string::npos) {
        if (filename.substr(found) == ".dat")
            filename = filename.substr(0, found);
        else {
            // The input file name must has an extension of 'dat'
            cout << "*** Error *** Invalid file extension: "
                 << filename.substr(found+1) << endl;
            exit(1);
        }
    }

    string InFile = filename + ".dat";
	string OutFile = filename + ".out";

	CDomain* FEMData = CDomain::GetInstance();

    // Set solver type if provided
    if (argc == 3) {
        unsigned int solverType = atoi(argv[2]);
        if (solverType > 2) {
            cout << "*** Warning *** Invalid solver type: " << solverType << endl;
            cout << "Using default solver type: 0 (Skyline LDLT)" << endl;
            solverType = 0;
        }
        FEMData->SetSolverType(solverType);
    }

    Clock timer;
    timer.Start();

//  Read data and define the problem domain
	if (!FEMData->ReadData(InFile, OutFile))
	{
		cerr << "*** Error *** Data input failed!" << endl;
		exit(1);
	}
    double time_input = timer.ElapsedTime();

    COutputter* Output = COutputter::GetInstance();
    if (!FEMData->GetMODEX())
    {
        *Output << "Data check completed !" << endl << endl;
        return 0;
    }
//  Allocate global vectors and matrices, such as the Force, ColumnHeights,
//  DiagonalAddress and StiffnessMatrix, and calculate the column heights
//  and address of diagonal elements
	FEMData->AllocateMatrices();
    
//  Assemble the banded gloabl stiffness matrix
	FEMData->AssembleStiffnessMatrix();

    
    double time_assemble = timer.ElapsedTime();

//  Solve the linear equilibrium equations for displacements
    // Get solver type (0: Skyline LDLT, 1: Eigen Direct, 2: Eigen CG)
    unsigned int solverType = FEMData->GetSolverType();
    
    // Create the appropriate solver
    CLDLTSolver* skylineSolver = nullptr;
    CEigenDirectSolver* eigenDirectSolver = nullptr;
    CEigenCGSolver* eigenCGSolver = nullptr;
    
    cout << "Using solver type: ";
    switch (solverType) {
        case 0: // Skyline LDLT
            cout << "Skyline LDLT (Original)" << endl;
            skylineSolver = new CLDLTSolver(FEMData->GetStiffnessMatrix());
            skylineSolver->LDLT(); // Perform L*D*L(T) factorization
            break;
        case 1: // Eigen Direct
            cout << "Eigen Direct Solver (SimplicialLDLT)" << endl;
            eigenDirectSolver = new CEigenDirectSolver(FEMData->GetEigenStiffnessMatrix());
            eigenDirectSolver->Factorize(); // Perform factorization
            break;
        case 2: // Eigen CG
            cout << "Eigen Iterative Solver (Conjugate Gradient)" << endl;
            eigenCGSolver = new CEigenCGSolver(FEMData->GetEigenStiffnessMatrix());
            break;
        default:
            cout << "Skyline LDLT (Default)" << endl;
            skylineSolver = new CLDLTSolver(FEMData->GetStiffnessMatrix());
            skylineSolver->LDLT(); // Perform L*D*L(T) factorization
            break;
    }

#ifdef _DEBUG_
    Output->PrintStiffnessMatrix();
#endif
        
//  Loop over for all load cases
    for (unsigned int lcase = 0; lcase < FEMData->GetNLCASE(); lcase++)
    {
//      Assemble righ-hand-side vector (force vector)
        FEMData->AssembleForce(lcase + 1);
            
//      Reduce right-hand-side force vector and back substitute
        switch (solverType) {
            case 0: // Skyline LDLT
                skylineSolver->BackSubstitution(FEMData->GetForce());
                break;
            case 1: // Eigen Direct
                eigenDirectSolver->Solve(FEMData->GetForce());
                break;
            case 2: // Eigen CG
                eigenCGSolver->Solve(FEMData->GetForce());
                break;
            default:
                skylineSolver->BackSubstitution(FEMData->GetForce());
                break;
        }

        *Output << " LOAD CASE" << setw(5) << lcase + 1 << endl << endl << endl;

#ifdef _DEBUG_
        Output->PrintDisplacement();
#endif
            
        Output->OutputNodalDisplacement();

//      Calculate and output stresses of all elements
        Output->OutputElementStress();
    }

    double time_solution = timer.ElapsedTime();
    
    timer.Stop();
    
    *Output << "\n S O L U T I O N   T I M E   L O G   I N   S E C \n\n"
            << "     TIME FOR INPUT PHASE = " << time_input << endl
            << "     TIME FOR CALCULATION OF STIFFNESS MATRIX = " << time_assemble - time_input << endl
            << "     TIME FOR FACTORIZATION AND LOAD CASE SOLUTIONS = " << time_solution - time_assemble << endl << endl
            << "     T O T A L   S O L U T I O N   T I M E = " << time_solution << endl << endl;

    // 导出VTK文件，便于ParaView后处理
    Output->OutputVTK(filename + ".vtk");

    // Clean up solvers
    delete skylineSolver;
    delete eigenDirectSolver;
    delete eigenCGSolver;

	return 0;
}
