#pragma once

// std c++ headers
#include <map>
#include <string>
#include <vector>

// PETSc headers
#include <petsc.h>
#include <petscmat.h>
#include <petscvec.h>

// AMDiS includes
#include "MatrixStreams.hpp"

namespace AMDiS
{
  /// pair of nested matrices and blockmatrix
  struct PetscMatrix
  {
    PetscMatrix()
      :	matrix(PETSC_NULL),
        isNested(false),
        assembled(false) {}

    ~PetscMatrix()
    {
      destroy();
    }

    std::vector<Mat> nestMat;
    Mat matrix;

    void destroy();
    bool isNested;
    bool assembled;
  };

  struct PetscMatrixNested : public PetscMatrix
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

    ~PetscVector()
    {
      destroy();
    }

    std::vector<Vec> nestVec;
    Vec vector;

    void destroy();
    bool isNested;
    bool assembled;
  };

  struct PetscVectorNested : public PetscVector
  {
    PetscVectorNested() : PetscVector()
    {
      isNested = true;
    }
  };

  /// structure for parameter mapping
  struct PetscParameters
  {
    std::map<std::string, bool> matSolverPackage;
    std::map<std::string, std::string> solverMap;
    std::map<std::string, std::string> preconMap;
    std::map<std::string, bool> emptyParam;

    PetscParameters();
  };


  // matrix-streams for Petsc Types
  // ---------------------------------------------------------------------------

  /// fill Petsc Mat datastructure directly
  inline void operator<<(Mat& mat, DOFMatrix const& rhs);

  /// fill PetscMatrix
  template <class Mapper>
  inline void operator<<(PetscMatrix& mat, 
			 MatMap<const SolverMatrix<Matrix<DOFMatrix*>>, Mapper>& rhs);

  /// fill nested PetscMatrix
  template <class Mapper>
  inline void operator<<(PetscMatrixNested& mat, 
			 MatMap<const SolverMatrix<Matrix<DOFMatrix*>>, Mapper>& mapper);


  /// fill Petsc Vec datastructure directly from DOFVector
  inline void operator<<(Vec& petscVec, DOFVector<double> const& vec);

  /// fill Petsc Vec datastructure directly from SystemVector
  inline void operator<<(Vec& petscVec, SystemVector const& vec);

  inline void operator<<(PetscVector& petscVec, SystemVector const& rhs);

  template <class Mapper>
  inline void operator<<(PetscVector& petscVec, 
			 VecMap<const SystemVector, Mapper>& rhs);

  template <class Mapper>
  inline void operator<<(PetscVectorNested& petscVec, 
			 VecMap<const SystemVector, Mapper>& rhs);


  /// fill AMDiS DOFVector from Petsc Vec datastructur
  inline void operator>>(Vec const& petscVec, DOFVector<double>& vec);

  /// fill AMDiS SystemVector from Petsc Vec datastructur
  inline void operator>>(Vec const& petscVec, SystemVector& vec);

  /// fill AMDiS SystemVector from PetscVector
  template <class Mapper>
  void operator>>(PetscVector const& dest, 
		  VecMap<SystemVector, Mapper>& rhs);

  /// fill AMDiS SystemVector from nested PetscVector
  template <class Mapper>
  void operator>>(PetscVectorNested const& dest, 
		  VecMap<SystemVector, Mapper>& rhs);

} // end namespace AMDiS

#include "solver/PetscTypes.hh"
