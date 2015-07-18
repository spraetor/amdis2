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


#include "parallel/ParallelProblemStat.h"
#include "parallel/ParallelSolver.h"
#include "parallel/MeshDistributor.h"
#include "parallel/MpiHelper.h"

#include "parallel/ParallelMapper.h"
#include "solver/LinearSolverInterface.h"

#ifdef HAVE_PARALLEL_MTL4
  #include "parallel/PITL_Solver.h"
  #include "solver/KrylovPreconditioner.h"
#elif defined HAVE_PARALLEL_PETSC
  #include "parallel/PetscSolverFeti.h"
  #include "parallel/PetscSolverSchur.h"
  #include "parallel/PetscSolverGlobalBlockMatrix.h"
  #include "parallel/PetscSolverGlobalMatrix.h"
  #include "parallel/PetscSolverNavierStokes.h"
#endif
#include "Global.h"


#if defined HAVE_PARALLEL_PETSC || defined HAVE_SEQ_PETSC
  #include "solver/PetscTypes.h"
#endif

namespace AMDiS { namespace Parallel {

  double ParallelProblemStat::initTimeStamp = 0.0;
  bool ParallelProblemStat::initialized = false;


  ParallelProblemStat::ParallelProblemStat(std::string nameStr,
					   ProblemIterationInterface *problemIteration)
    : ProblemStatSeq(nameStr, problemIteration),
      meshDistributor(NULL)
  {
    initTimeStamp = MPI::Wtime();
    mpi::globalMin(initTimeStamp);
    addSolvers();
  }


  void ParallelProblemStat::initialize(Flag initFlag,
					   ProblemStatSeq *adoptProblem,
					   Flag adoptFlag)
  {
    FUNCNAME("ParallelProblemStat::initialize()");
    MSG("ParallelProblemStat::initialize()\n");
    
    MSG("Initialization phase 0 needed %.5f seconds\n", 
	MPI::Wtime() - initTimeStamp);

    ProblemStatSeq::initialize(initFlag, adoptProblem, adoptFlag);

    MeshDistributor::addProblemStatGlobal(this);
    meshDistributor = MeshDistributor::globalMeshDistributor;
    meshDistributor->addInterchangeVector(getSolution());
        
    ParallelSolver *parallelSolver = dynamic_cast<ParallelSolver*>(solver);
    TEST_EXIT(parallelSolver != NULL)
      ("ParallelProblem loaded, but no ParallelSolver selected! This does not fit together.\n");

    parallelSolver->setMeshDistributor(meshDistributor, 0);
    
    // For the additional component with extra mesh and fespace, we don't need to put them
    // into the parallel system.
    std::vector<const FiniteElemSpace*> tmpFeSpaces, tmpComponentSpaces(getComponentSpaces().begin(),
							   getComponentSpaces().begin() + getNumComponents());
							   
    for (size_t i = 0; i < tmpComponentSpaces.size(); i++)
      if (std::find(tmpFeSpaces.begin(), tmpFeSpaces.end(), tmpComponentSpaces[i]) == tmpFeSpaces.end())
	tmpFeSpaces.push_back(tmpComponentSpaces[i]);
      
    parallelSolver->init(tmpComponentSpaces, tmpFeSpaces);
  }


  void ParallelProblemStat::buildAfterCoarsen(AdaptInfo *adaptInfo, Flag flag,
						  bool assembleMatrix,
						  bool assembleVector)
  {
    FUNCNAME("ParallelProblemStat::buildAfterCoarsen()");

    TEST_EXIT(MeshDistributor::globalMeshDistributor != NULL)
      ("No Meshdistributor! Should not happen!\n");

    MeshDistributor::globalMeshDistributor->checkMeshChange();
    ProblemStatSeq::buildAfterCoarsen(adaptInfo, flag, 
				      assembleMatrix, assembleVector);
  }


  void ParallelProblemStat::addPeriodicBC(BoundaryType type, int row, int col)
  {
    if (MeshDistributor::globalMeshDistributor->isInitialized())
      return;

    ProblemStatSeq::addPeriodicBC(type, row, col);
  }
  

  void ParallelProblemStat::addSolvers()
  {
    if (!initialized) {
    initialized = true;
    LinearSolverCreator *creator;

#if defined HAVE_PARALLEL_MTL4
    creator = new P_CGSolver::Creator;
    CreatorMap< LinearSolverInterface >::addCreator("p_mtl_cg", creator);

    creator = new P_CGSSolver::Creator;
    CreatorMap< LinearSolverInterface >::addCreator("p_mtl_cgs", creator);

//     creator = new P_BiCGSolver::Creator;
//     CreatorMap< LinearSolverInterface >::addCreator("p_mtl_bicg", creator);

    creator = new P_BiCGStabSolver::Creator;
    CreatorMap< LinearSolverInterface >::addCreator("p_mtl_bicgstab", creator);

    creator = new P_BiCGStab2Solver::Creator;
    CreatorMap< LinearSolverInterface >::addCreator("p_mtl_bicgstab2", creator);

    creator = new P_BiCGStabEllSolver::Creator;
    CreatorMap< LinearSolverInterface >::addCreator("p_mtl_bicgstab_ell", creator);

    creator = new P_QMRSolver::Creator;
    CreatorMap< LinearSolverInterface >::addCreator("p_mtl_qmr", creator);

    creator = new P_TFQMRSolver::Creator;
    CreatorMap< LinearSolverInterface >::addCreator("p_mtl_tfqmr", creator);

    creator = new P_GMResSolver::Creator;
    CreatorMap< LinearSolverInterface >::addCreator("p_mtl_gmres", creator);

    creator = new P_FGMResSolver::Creator;
    CreatorMap< LinearSolverInterface >::addCreator("p_mtl_fgmres", creator);

    creator = new P_IDRsSolver::Creator;
    CreatorMap< LinearSolverInterface >::addCreator("p_mtl_idr_s", creator);

    creator = new P_MinResSolver::Creator;
    CreatorMap< LinearSolverInterface >::addCreator("p_mtl_minres", creator);
    
    creator = new P_PreOnly::Creator;
    CreatorMap< LinearSolverInterface >::addCreator("p_mtl_preonly", creator);
    CreatorMap< LinearSolverInterface >::addCreator("p_mtl_richardson", creator);

#elif defined HAVE_PARALLEL_PETSC
    creator = new PetscSolverGlobalMatrix::Creator;
    CreatorMap< LinearSolverInterface >::addCreator("p_petsc_petsc", creator); // standard PETSc creator
    
    std::map<std::string,std::string>::iterator it;
    PetscParameters params;
    for (it = params.solverMap.begin();
	 it!= params.solverMap.end();
	 it++) {
      CreatorMap< LinearSolverInterface >::addCreator("p_petsc_" + it->first, creator);
    }  
    
    creator = new PetscSolverSchur::Creator;
    CreatorMap< LinearSolverInterface >::addCreator("p_petsc_petsc-schur", creator);
    
    creator = new PetscSolverGlobalBlockMatrix::Creator;
    CreatorMap< LinearSolverInterface >::addCreator("p_petsc_petsc-block", creator);
    
    creator = new PetscSolverFeti::Creator;
    CreatorMap< LinearSolverInterface >::addCreator("p_petsc_petsc-feti", creator);
      
    creator = new PetscSolverNavierStokes::Creator;
    CreatorMap< LinearSolverInterface >::addCreator("p_petsc_petsc-navierstokes", creator);
    CreatorMap< LinearSolverInterface >::addCreator("petsc-navierstokes", creator);
    
#elif defined  HAVE_BDDC_ML    
    creator = new BddcMlSolver::Creator;
    CreatorMap< LinearSolverInterface >::addCreator("bddcml", creator);
#endif
    }
  }
  
} // end namespace Parallel

  
#ifdef HAVE_PARALLEL_MTL4
  template< > 
  void CreatorMap<Parallel::ParallelPreconditioner>::addDefaultCreators() 
  {
    Parallel::ParallelPreconditionCreator *creator;
    
    creator =  new Parallel::P_DiagonalPreconditioner::Creator;
    addCreator("diag", creator);
    
    creator = new Parallel::P_ILUPreconditioner::Creator;
    addCreator("ilu", creator);
    
    creator = new Parallel::P_ICPreconditioner::Creator;
    addCreator("ic", creator);
    
    creator =  new Parallel::P_IdentityPreconditioner::Creator;
    addCreator("no", creator);

    creator =  new KrylovPreconditionerParallel::Creator;
    addCreator("krylov", creator);
    addCreator("solver", creator);
  }
#endif
  
} // end namespace AMDiS
