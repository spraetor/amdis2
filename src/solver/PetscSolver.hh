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


#ifdef HAVE_SEQ_PETSC

namespace AMDiS {  
  
  template<typename M, typename V>
  PetscRunner<M,V>::PetscRunner(LinearSolverInterface* oem_)
    : oem(*oem_),
      kspPrefix(""),
      zeroStartVector(false),
      initialized(false),
      preconditioner(NULL)
  {
    PetscParameters params;
    matSolverPackage = false;
    
    // set the parameter prefix
    Parameters::get(oem.getName() + "->petsc prefix", kspPrefix); // important when more than one solver to configure
    
    // set the solver
    std::string solverName = "petsc";
    Parameters::get(oem.getName(), solverName);
    if (solverName == "petsc") 
      Parameters::get(oem.getName() + "->ksp_type", solverName);
    
    std::string kspSolver = params.solverMap[solverName];
    
    if (params.matSolverPackage[kspSolver]) {
      // direct solvers
      PetscOptionsInsertString(("-" + kspPrefix + "ksp_type preonly").c_str());
      PetscOptionsInsertString(("-" + kspPrefix + "pc_type lu").c_str());
      PetscOptionsInsertString(("-" + kspPrefix + "pc_factor_mat_solver_package " + (kspSolver != "direct" ? kspSolver : "umfpack")).c_str());
      oem.setMaxIterations(1);
      zeroStartVector = true;
      matSolverPackage = true;
    } else if (!params.emptyParam[kspSolver]) {    
      // other solvers
      PetscOptionsInsertString(("-" + kspPrefix + "ksp_type " + kspSolver).c_str());
    }
    
    // set the preconditioner
    setPrecon();
    
    if (oem.getInfo() >= 20)
      PetscOptionsInsertString(("-" + kspPrefix + "ksp_monitor_true_residual").c_str());
    else if (oem.getInfo() >= 10)
      PetscOptionsInsertString(("-" + kspPrefix + "ksp_monitor").c_str());
    
    // command line string
    std::string kspString = "";
    Parameters::get(oem.getName() + "->ksp", kspString);
    if (kspString != "")
      PetscOptionsInsertString(kspString.c_str());
  }

  
  template<typename MatrixType, typename VectorType>
  void PetscRunner<MatrixType, VectorType>::init(const SolverMatrix<Matrix<DOFMatrix*> >& A, const MatrixType& fullMatrix)
  {
    if (initialized)
      exit();
    
    createSubSolver(ksp, fullMatrix.matrix, kspPrefix);
    
    // set tolerance options from LinearSolverInterface-parameters
    KSPSetTolerances(ksp, oem.getRelative(), oem.getTolerance(), PETSC_DEFAULT, oem.getMaxIterations());

    // Do not delete the solution vector, use it for the initial guess.
    if (!zeroStartVector)
      KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
  
    KSPGetPC(ksp, &pc);  
    
    preconditioner->init(pc, A, fullMatrix);    
    initialized = true;
  }

  
  template<typename MatrixType, typename VectorType>
  int PetscRunner<MatrixType, VectorType>::solve(const MatrixType& A, VectorType& x, const VectorType& b) 
  {
    if (zeroStartVector)
      VecSet(x.vector, 0.0);
    
    PetscErrorCode solverError = KSPSolve(ksp, b.vector, x.vector);
    
    PetscInt nIter = 0;
    KSPGetIterationNumber(ksp, &nIter);
    PetscReal residual_norm = -1.0;
    KSPGetResidualNorm(ksp, &residual_norm);
    
    if (residual_norm <= 0.0) {
      Vec r;
      VecDuplicate(b.vector, &r);
      KSPBuildResidual(ksp, PETSC_NULL, PETSC_NULL, &r);
      VecNorm(r, NORM_2, &residual_norm);
    }
    
    oem.setErrorCode(solverError);
    oem.setIterations(nIter);
    oem.setResidual(residual_norm);
    
    return solverError;
  }

  
  template<typename M, typename V>
  void PetscRunner<M,V>::createSubSolver(KSP &ksp_, Mat m, std::string kspPrefix_)
  {
    KSPCreate(PETSC_COMM_SELF, &ksp_);
#if (PETSC_VERSION_MINOR >= 5)
    KSPSetOperators(ksp_, m, m);
#else
    KSPSetOperators(ksp_, m, m, SAME_NONZERO_PATTERN);
#endif
    KSPSetOptionsPrefix(ksp_, kspPrefix_.c_str());
    KSPSetFromOptions(ksp_);
  }


  template<typename M, typename V>
  void PetscRunner<M,V>::setSolver(KSP ksp_, std::string kspPrefix_,
				   KSPType kspType,  PCType pcType, 
				   PetscReal rtol, PetscReal atol, PetscInt maxIt)
  {
    KSPSetType(ksp_, kspType);
    KSPSetTolerances(ksp_, rtol, atol, PETSC_DEFAULT, maxIt);
    if (kspPrefix_ != "")
      KSPSetOptionsPrefix(ksp_, kspPrefix_.c_str());
    KSPSetFromOptions(ksp_);

    PC pc_;
    KSPGetPC(ksp_, &pc_);
    PCSetType(pc_, pcType);
    PCSetFromOptions(pc_);
  }
  
  
  template<typename MatrixType, typename VectorType>
  void PetscRunner<MatrixType, VectorType>::setPrecon()
  {
    std::string precon = "";
    std::string initFileStr = oem.getName() + "->pc_type";
    Parameters::get(initFileStr, precon);
    if (!precon.size() || precon == "no" || precon == "0") {
      initFileStr = oem.getName() + "->left precon";
      Parameters::get(initFileStr, precon);
    } 
    if (!precon.size() || precon == "no" || precon == "0") {
      initFileStr = oem.getName() + "->right precon";
      Parameters::get(initFileStr, precon);
    }
    // no preconditioner should be created
    if (!precon.size() || precon == "no" || precon == "0")
      return;
    
    PetscParameters params;
    if (preconditioner) {
      delete preconditioner;
      preconditioner = NULL;
    }
    
    // Creator for the left preconditioner
    CreatorInterfaceName< PetscPreconditionerInterface<MatrixType, VectorType> >* creator(NULL);

    try {
      creator = dynamic_cast<CreatorInterfaceName< PetscPreconditionerInterface<MatrixType, VectorType> >*>(
	CreatorMap<PetscPreconditionerInterface<MatrixType, VectorType> >::getCreator(precon, initFileStr) );
    } catch (...) {
      creator = NULL;
    }
    
    if (creator) {
      creator->setName(initFileStr);
      preconditioner = creator->create();
    } else if (!matSolverPackage && !params.emptyParam[precon]) {
      precon = (params.preconMap.find(precon) != params.preconMap.end() ? params.preconMap[precon] : precon);
      preconditioner = new PetscPreconditionerInterface<MatrixType, VectorType>(kspPrefix, precon);
    } else {
      ERROR_EXIT((std::string("There is no creator for the given preconditioner '") + precon + "'").c_str());
    }
  }
}

#endif // HAVE_SEQ_PETSC
