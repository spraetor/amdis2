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


#ifndef AMDIS_PETSCTYPES_H
#define AMDIS_PETSCTYPES_H

#include <vector>
#include <map>
#include <string>

#include "MatrixStreams.h"

#include <petsc.h>
#include <petscmat.h> 
#include <petscvec.h>

namespace AMDiS {
  
  /// pair of nested matrices and blockmatrix
  struct PetscMatrix
  {
    PetscMatrix() 
      :	matrix(PETSC_NULL),
	isNested(false),
	assembled(false) {}
	
    ~PetscMatrix() { destroy(); }
    
    std::vector<Mat> nestMat;
    Mat matrix;
    
    void destroy();
    bool isNested;
    bool assembled;
  };
  
  struct PetscMatrixNested : PetscMatrix
  {
    PetscMatrixNested() : PetscMatrix() 
    {
      isNested = true;
    }
  };
  
  
  /// pair of nested vectors and blockvector
  struct PetscVector
  {
    PetscVector() 
      : vector(PETSC_NULL), 
	isNested(false),
	assembled(false) {}
	
    ~PetscVector() { destroy(); }
    
    std::vector<Vec> nestVec;
    Vec vector;
    
    void destroy();
    bool isNested;
    bool assembled;
  };
  
  struct PetscVectorNested : PetscVector
  {
    PetscVectorNested() : PetscVector()
    {
      isNested = true;
    }
  };
  
  /// structure for parameter mapping
  struct PetscParameters
  {
    std::map<std::string,bool> matSolverPackage;
    std::map<std::string,std::string> solverMap;
    std::map<std::string,std::string> preconMap;
    std::map<std::string,bool> emptyParam;
    
    PetscParameters();
  };

  
  // matrix-streams for Petsc Types
  // ___________________________________________________________________________
  
  /// fill Petsc Mat datastructure directly
  inline void operator<<(Mat& mat, const DOFMatrix& rhs);
  
  /// fill PetscMatrix
  template< typename Mapper >
  inline void operator<<(PetscMatrix& mat, MatMap<const SolverMatrix<Matrix<DOFMatrix* > >, Mapper >& rhs);
  
  /// fill nested PetscMatrix
  template< typename Mapper >
  inline void operator<<(PetscMatrixNested& mat, MatMap<const SolverMatrix<Matrix<DOFMatrix* > >, Mapper >& mapper);

  
  /// fill Petsc Vec datastructure directly from DOFVector
  inline void operator<<(Vec& petscVec, const DOFVector<double>& vec);
  
  /// fill Petsc Vec datastructure directly from SystemVector
  inline void operator<<(Vec& petscVec, const SystemVector& vec);
  
  inline void operator<<(PetscVector& petscVec, const SystemVector& rhs);
  
  template< typename Mapper >
  inline void operator<<(PetscVector& petscVec, VecMap<const SystemVector, Mapper>& rhs);
  
  template< typename Mapper >
  inline void operator<<(PetscVectorNested& petscVec, VecMap<const SystemVector, Mapper>& rhs);
  
  
  /// fill AMDiS DOFVector from Petsc Vec datastructur
  inline void operator>>(const Vec& petscVec, DOFVector<double>& vec);
  
  /// fill AMDiS SystemVector from Petsc Vec datastructur
  inline void operator>>(const Vec& petscVec, SystemVector& vec);
  
  /// fill AMDiS SystemVector from PetscVector
  template< typename Mapper >
  void operator>>(const PetscVector& dest, VecMap<SystemVector, Mapper>& rhs);
  
  /// fill AMDiS SystemVector from nested PetscVector
  template< typename Mapper >
  void operator>>(const PetscVectorNested& dest, VecMap<SystemVector, Mapper>& rhs);
  
} // end namespace AMDiS

#include "solver/PetscTypes.hh"

#endif // AMDIS_PETSCTYPES_H
