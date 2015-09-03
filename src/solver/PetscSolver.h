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


/** \file PetscSolver.h */

#ifndef AMDIS_PETSC_SOLVER_SEQ_H
#define AMDIS_PETSC_SOLVER_SEQ_H

#ifdef HAVE_SEQ_PETSC

#include "solver/LinearSolver.h"
#include "solver/PetscTypes.h"
#include "solver/MatrixStreams.h"
#include "Timer.h"
#include <vector>
#include <iostream>
#include <boost/mpl/bool.hpp>

#include <petsc.h>
#include <petscksp.h>
#include <petscmat.h>
#include <petscvec.h>
#include <petscsys.h>
#include <petscao.h>

namespace AMDiS
{

  /**
   * \ingroup Solver
   *
   * \brief Common base class for wrappers to use Petsc preconditioners in AMDiS.
   */
  template<class MatrixType, class VectorType>
  struct PetscPreconditionerInterface : public PreconditionerInterface
  {
    PetscPreconditionerInterface(std::string prefix, std::string name)
      : prefix(prefix), name(name) {}

    virtual ~PetscPreconditionerInterface() {}

    virtual void init(PC pc, const SolverMatrix<Matrix<DOFMatrix*>>& A, const MatrixType& fullMatrix)
    {
      PetscOptionsInsertString(("-" + prefix + "pc_type " + name).c_str());
      PCSetFromOptions(pc);
    }

    virtual void exit() {}

  protected:
    std::string prefix;
    std::string name;
  };


  /// PETSc preconditioners using MATSEQAIJ and VECSEQ data structures
  typedef PetscPreconditionerInterface<PetscMatrix, PetscVector> PetscPreconditioner;

  /// PETSc preconditioner using MATNEST and VECNEST data structures
  typedef PetscPreconditionerInterface<PetscMatrixNested, PetscVectorNested> PetscPreconditionerNested;


  /// Runner that can be passed to LinearSolver to redirect the solution to PETSc
  template<typename MatrixType, typename VectorType>
  class PetscRunner : public RunnerInterface
  {
  public:
    /// Constructor of standard PETSc runner. Reads ksp and pc parameters from initfile.
    PetscRunner(LinearSolverInterface* oem);

    /// Destructor, deletes preconditioner object \ref preconditioner
    ~PetscRunner()
    {
      if (preconditioner)
      {
        delete preconditioner;
        preconditioner = NULL;
      }
    }

    /// Initialize the solver \ref ksp and preconditioner \ref pc
    void init(const SolverMatrix<Matrix<DOFMatrix*>>& A, const MatrixType& fullMatrix);

    /// Solve the linear equation \f$ A\cdot x = b \f$ by applying the PETSc solver \ref ksp
    int solve(const MatrixType& A, VectorType& x, const VectorType& b);

    /// Destroy solver \ref ksp and preconditioner \ref pc
    virtual void exit()
    {
      preconditioner->exit();
      KSPDestroy(&ksp);
    }

    /// Get the PETSc solver \ref ksp
    KSP getSolver()
    {
      return ksp;
    }

    /// Get the PETSc preconditioner \ref pc
    PC getPc()
    {
      return pc;
    }

    /// Returns preconditioner object \ref preconditioner
    PreconditionerInterface* getLeftPrecon()
    {
      return preconditioner;
    }

    /// Returns preconditioner object \ref preconditioner
    PreconditionerInterface* getRightPrecon()
    {
      return preconditioner;
    }


    static void createSubSolver(KSP& ksp_, Mat m, std::string kspPrefix_);

    static void setSolver(KSP ksp_, std::string kspPrefix_,
                          KSPType kspType, PCType pcType,
                          PetscReal rtol = PETSC_DEFAULT,
                          PetscReal atol = PETSC_DEFAULT,
                          PetscInt maxIt = PETSC_DEFAULT);


  protected:
    /// Creates a preconditioner object \ref preconditioner
    void setPrecon();

  protected:
    LinearSolverInterface& oem;

    /// PETSc solver object
    KSP ksp;

    /// PETSc preconditioner object
    PC pc;

    /// KSP database prefix
    std::string kspPrefix;

    bool zeroStartVector;	///< initialize startsolution for iterative solver with 0
    bool initialized;		///< true when \ref init() was called
    bool matSolverPackage;	///< true when using a direct solver

    /// Preconditioner object, that sets parameters to \ref pc
    PetscPreconditionerInterface<MatrixType, VectorType>* preconditioner;
  };


  /**
   * \ingroup Solver
   *
   * \brief
   * Wrapper for the external PETSc solver:
   *   http://www.mcs.anl.gov/petsc/
   *
   * This is a suite of data structures and routines for the
   * scalable (parallel) solution of scientific applications.
   */
  typedef LinearSolver<PetscMatrix, PetscVector, PetscRunner<PetscMatrix, PetscVector>> PetscSolver;
  typedef LinearSolver<PetscMatrixNested, PetscVectorNested, PetscRunner<PetscMatrixNested, PetscVectorNested>> PetscSolverNested;

} // end namespace AMDiS

#endif // HAVE_SEQ_PETSC

#include "solver/PetscSolver.hh"

#endif // AMDIS_PETSC_SOLVER_SEQ_H
