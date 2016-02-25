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


/** \file PetscSolverSchur.h */

#include "parallel/PetscSolver.hpp"

#ifndef AMDIS_PETSC_SOLVER_SCHUR_H
#define AMDIS_PETSC_SOLVER_SCHUR_H

namespace AMDiS
{
  namespace Parallel
  {

    /** \ingroup Solver
     *
     * \brief
     * PETSc solver using fieldsplit
     */
    class PetscSolverSchur : public PetscSolver
    {
    public:
      /// Creator class
      class Creator : public LinearSolverCreator
      {
      public:
        virtual ~Creator() {}

        /// Returns a new PetscSolver object.
        LinearSolverInterface* create()
        {
          return new PetscSolverSchur(this->name);
        }
      };

      PetscSolverSchur(std::string name)
        : PetscSolver(name)
      {}

      void fillPetscMatrix(Matrix<DOFMatrix*>* mat);

      void fillPetscRhs(SystemVector* vec);

      void solvePetscMatrix(SystemVector& vec, AdaptInfo& adaptInfo);

      void destroyMatrixData()
      {}

      void destroyVectorData()
      {}

      Flag getBoundaryDofRequirement()
      {
        return
          MeshDistributor::BOUNDARY_SUBOBJ_SORTED |
          MeshDistributor::BOUNDARY_FILL_INFO_SEND_DOFS;
      }

    protected:
      void updateDofData(int nComponents);

      /// Takes a DOF matrix and sends the values to the global PETSc matrix.
      void setDofMatrix(DOFMatrix* mat, int dispMult = 1,
                        int dispAddRow = 0, int dispAddCol = 0);

      /// Takes a DOF vector and sends its values to a given PETSc vector.
      void setDofVector(Vec& petscVec, DOFVector<double>* vec,
                        int disMult = 1, int dispAdd = 0, bool rankOnly = false);

    protected:
      int nBoundaryDofs;

      int rStartBoundaryDofs;

      int nOverallBoundaryDofs;

      std::set<DegreeOfFreedom> boundaryDofs;

      std::map<DegreeOfFreedom, DegreeOfFreedom> mapGlobalBoundaryDof;

      int nInteriorDofs;

      int rStartInteriorDofs;

      int nOverallInteriorDofs;

      std::set<DegreeOfFreedom> interiorDofs;

      std::map<DegreeOfFreedom, DegreeOfFreedom> mapGlobalInteriorDof;

      Mat matA11, matA12, matA21, matA22;

      IS interiorIs, boundaryIs;

      Vec petscSolVec;
    };

  } // end namespace Parallel

} // end namespace AMDiS

#endif
