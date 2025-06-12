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

#include "Eigen/Sparse"
#include "Eigen/Dense"
#include <vector>
#include <iostream>

//! CEigenMatrix class is used to store the FEM stiffness matrix using Eigen sparse matrix
class CEigenMatrix
{
private:
    //! Dimension of the stiffness matrix
    unsigned int NEQ_;

    //! Eigen sparse matrix to store the stiffness matrix
    Eigen::SparseMatrix<double> matrix_;

    //! Triplet list for assembling the sparse matrix
    std::vector<Eigen::Triplet<double>> tripletList_;

public:
    //! Constructor
    CEigenMatrix() : NEQ_(0) {}

    //! Constructor with dimension
    CEigenMatrix(unsigned int N) : NEQ_(N), matrix_(N, N) {}

    //! Destructor
    ~CEigenMatrix() {}

    //! Initialize the matrix with dimension N
    void Initialize(unsigned int N) {
        NEQ_ = N;
        matrix_.resize(N, N);
        matrix_.setZero();
        tripletList_.clear();
    }

    //! Reserve space for non-zero entries
    void Reserve(unsigned int nonZeros) {
        tripletList_.reserve(nonZeros);
    }

    //! Add an element to the triplet list
    void AddElement(unsigned int i, unsigned int j, double value) {
        // Convert from 1-based to 0-based indexing
        tripletList_.push_back(Eigen::Triplet<double>(i-1, j-1, value));
    }

    //! Assemble the sparse matrix from triplet list
    void Assemble() {
        matrix_.setFromTriplets(tripletList_.begin(), tripletList_.end());
        matrix_.makeCompressed();
        tripletList_.clear(); // Free memory
    }

    //! Assemble the element stiffness matrix to the global stiffness matrix
    void Assembly(double* Matrix, unsigned int* LocationMatrix, size_t ND) {
        // Assemble global stiffness matrix
        for (unsigned int j = 0; j < ND; j++) {
            unsigned int Lj = LocationMatrix[j];    // Global equation number corresponding to jth DOF of the element
            if (!Lj) continue;
            
            // Address of diagonal element of column j in the one dimensional element stiffness matrix
            unsigned int DiagjElement = (j+1)*j/2;
            
            for (unsigned int i = 0; i <= j; i++) {
                unsigned int Li = LocationMatrix[i];    // Global equation number corresponding to ith DOF of the element
                
                if (!Li) continue;
                
                // Add to triplet list (both i,j and j,i for symmetry)
                AddElement(Li, Lj, Matrix[DiagjElement + j - i]);
                if (i != j) {
                    AddElement(Lj, Li, Matrix[DiagjElement + j - i]);
                }
            }
        }
    }

    //! Return the dimension of the stiffness matrix
    inline unsigned int dim() const {
        return NEQ_;
    }

    //! Get the Eigen sparse matrix
    inline Eigen::SparseMatrix<double>& GetMatrix() {
        return matrix_;
    }

    //! Get a specific element (1-based indexing for compatibility)
    double operator()(unsigned int i, unsigned int j) const {
        // Convert from 1-based to 0-based indexing
        return matrix_.coeff(i-1, j-1);
    }
};
