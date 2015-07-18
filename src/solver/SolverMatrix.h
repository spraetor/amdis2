/******************************************************************************
 *
 * AMDiS - Adaptive multidimensional simulations
 *
 * Copyright (C) 2013 Dresden University of Technology. All Rights Reserved.
 * Web: https://fusionforge.zih.tu-dresden.de/projects/amdis
 *
 * Authors: 
 * Simon Vey, Thomas Witkowski, Andreas Naumann, Simon Praetorius, et al.
 *
 * This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
 * WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 *
 *
 * This file is part of AMDiS
 *
 * See also license.opensource.txt in the distribution.
 * 
 ******************************************************************************/


/** \file SolverMatrix.h */

#ifndef AMDIS_SOLVERMATRIX_H
#define AMDIS_SOLVERMATRIX_H

#include "DOFMatrix.h"

namespace AMDiS {

  /**
   * \brief Helper class to provide complete matrix for ITL solvers and preconditioners
   */

  // General case undefined !
  template <typename Matrix>
  class SolverMatrix { };


  // DOFMatrix, i.e. ScalarProblem
  template <>
  class SolverMatrix<DOFMatrix>
  {
  public:
    SolverMatrix()
      : matrix(0)
    {}

    void setMatrix(const DOFMatrix& A)
    {
      matrix= &A.getBaseMatrix();
    }

    const DOFMatrix::base_matrix_type& getMatrix() const 
    { 
      return *matrix; 
    }

  private:
    const DOFMatrix::base_matrix_type* matrix;
  };


  // VectorProblem
  template <>
  class SolverMatrix<Matrix<DOFMatrix*> >
  {
  public :
    SolverMatrix<Matrix<DOFMatrix* > >()
      : originalMat(NULL)
    {}

    void setMatrix(const Matrix<DOFMatrix*>& A)
    {
      originalMat = &A;
      matrix.change_dim(0,0);
    }

    const DOFMatrix::base_matrix_type& getMatrix() const 
    { 
      if( num_rows(matrix) == 0)
        buildMatrix();
      return matrix; 
    }

    const DOFMatrix::base_matrix_type& getSubMatrix(int i, int j) const 
    { 
      TEST_EXIT(i>=0 && j>=0 && i<originalMat->getNumRows() && j<originalMat->getNumCols())
	("SubMatrix indices out of range!\n");
	
      return (*originalMat)[i][j]->getBaseMatrix(); 
    }

    const Matrix<DOFMatrix*>* getOriginalMat() const
    {
      return originalMat;
    }
    
  private:
    void buildMatrix() const;    
    
  private:      
    mutable DOFMatrix::base_matrix_type matrix;
    
    const Matrix<DOFMatrix*>* originalMat;
  };

} // end namespace AMDiS

#endif // AMDIS_SOLVERMATRIX_H
