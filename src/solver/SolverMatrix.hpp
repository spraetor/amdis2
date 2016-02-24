/** \file SolverMatrix.h */

#pragma once

// AMDiS includes
#include "DOFMatrix.hpp"

namespace AMDiS
{
  /**
   * \brief Helper class to provide complete matrix for ITL solvers and preconditioners
   */

  // Primary template
  template <typename Matrix>
  class SolverMatrix;


  // Specialization for DOFMatrix, i.e. ScalarProblem
  template <>
  class SolverMatrix<DOFMatrix>
  {
  public:
    SolverMatrix()
      : matrix(0)
    {}

    void setMatrix(DOFMatrix const& A)
    {
      matrix= &A.getBaseMatrix();
    }

    DOFMatrix::base_matrix_type const& getMatrix() const
    {
      return *matrix;
    }

  private:
    DOFMatrix::base_matrix_type const* matrix;
  };


  // Specialization for matrices of matrices, i.e. VectorProblem
  template <>
  class SolverMatrix<Matrix<DOFMatrix*>>
  {
    using MatrixType = Matrix<DOFMatrix*>;
    
  public :
    SolverMatrix()
      : originalMat(NULL)
    {}

    void setMatrix(MatrixType const& A)
    {
      originalMat = &A;
      matrix.change_dim(0,0);
    }

    DOFMatrix::base_matrix_type const& getMatrix() const
    {
      if( num_rows(matrix) == 0)
        buildMatrix();
      return matrix;
    }

    DOFMatrix::base_matrix_type const& getSubMatrix(int i, int j) const
    {
      TEST_EXIT(i>=0 && j>=0 && i < originalMat->getNumRows() 
			     && j < originalMat->getNumCols())
	("SubMatrix indices out of range!\n");

      return (*originalMat)[i][j]->getBaseMatrix();
    }

    MatrixType const* getOriginalMat() const
    {
      return originalMat;
    }

  private:
    void buildMatrix() const;

  private:
    mutable DOFMatrix::base_matrix_type matrix;

    MatrixType const* originalMat;
  };

} // end namespace AMDiS
