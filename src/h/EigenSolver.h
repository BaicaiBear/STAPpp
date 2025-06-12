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

#include "EigenMatrix.h"
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCholesky>
#include <memory>

//! Base class for all Eigen-based solvers
class CEigenSolverBase
{
protected:
    CEigenMatrix& K;
    
public:
    //! Constructor
    CEigenSolverBase(CEigenMatrix* K) : K(*K) {};
    
    //! Destructor
    virtual ~CEigenSolverBase() {};
    
    //! Solve the linear system
    virtual void Solve(double* Force) = 0;
};

//! Direct solver using Eigen's SimplicialLDLT
class CEigenDirectSolver : public CEigenSolverBase
{
private:
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver_;
    bool factorized_;
    
public:
    //! Constructor
    CEigenDirectSolver(CEigenMatrix* K) : CEigenSolverBase(K), factorized_(false) {};
    
    //! Factorize the matrix
    void Factorize() {
        // Assemble the matrix if not already done
        K.Assemble();
        
        // 替换 SimplicialLDLT 为 SparseLU
        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver_;
        solver_.analyzePattern(K.GetMatrix());
        solver_.factorize(K.GetMatrix());
        if (solver_.info() != Eigen::Success) {
            std::cerr << "*** Error *** LU Factorization failed!" << std::endl;
            exit(4);
        }
        
        factorized_ = true;
    }
    
    //! Solve the linear system
    void Solve(double* Force) override {
        if (!factorized_) {
            Factorize();
        }
        
        // Convert Force to Eigen vector
        Eigen::VectorXd b(K.dim());
        for (unsigned int i = 0; i < K.dim(); ++i) {
            b(i) = Force[i];
        }
        
        // Solve the system
        Eigen::VectorXd x = solver_.solve(b);
        
        if (solver_.info() != Eigen::Success) {
            std::cerr << "*** Error *** Solving failed!" << std::endl;
            exit(4);
        }
        
        // Copy the solution back to Force
        for (unsigned int i = 0; i < K.dim(); ++i) {
            Force[i] = x(i);
        }
    }
};

//! Iterative solver using Eigen's ConjugateGradient
class CEigenCGSolver : public CEigenSolverBase
{
private:
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> solver_;
    double tolerance_;
    int maxIterations_;
    
public:
    //! Constructor
    CEigenCGSolver(CEigenMatrix* K, double tolerance = 1e-6, int maxIterations = 1000) 
        : CEigenSolverBase(K), tolerance_(tolerance), maxIterations_(maxIterations) {
        solver_.setTolerance(tolerance_);
        solver_.setMaxIterations(maxIterations_);
    }
    
    //! Set tolerance
    void SetTolerance(double tolerance) {
        tolerance_ = tolerance;
        solver_.setTolerance(tolerance_);
    }
    
    //! Set maximum iterations
    void SetMaxIterations(int maxIterations) {
        maxIterations_ = maxIterations;
        solver_.setMaxIterations(maxIterations_);
    }
    
    //! Solve the linear system
    void Solve(double* Force) override {
        // Assemble the matrix if not already done
        K.Assemble();
        
        // Convert Force to Eigen vector
        Eigen::VectorXd b(K.dim());
        for (unsigned int i = 0; i < K.dim(); ++i) {
            b(i) = Force[i];
        }
        
        // Initialize the solver
        solver_.compute(K.GetMatrix());
        
        if (solver_.info() != Eigen::Success) {
            std::cerr << "*** Error *** Solver initialization failed!" << std::endl;
            exit(4);
        }
        
        // Solve the system
        Eigen::VectorXd x = solver_.solve(b);
        
        if (solver_.info() != Eigen::Success) {
            std::cerr << "*** Error *** Solving failed!" << std::endl;
            exit(4);
        }
        
        // Output solver statistics
        std::cout << "CG solver converged in " << solver_.iterations() << " iterations." << std::endl;
        std::cout << "CG solver error: " << solver_.error() << std::endl;
        
        // Copy the solution back to Force
        for (unsigned int i = 0; i < K.dim(); ++i) {
            Force[i] = x(i);
        }
    }
};
