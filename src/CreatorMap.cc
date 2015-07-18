#include "BasisFunction.h"
#include "Bubble.h"
#include "CreatorMap.h"
#include "MTL4Types.h"
#include "solver/LinearSolverInterface.h"
#include "solver/ITL_Solver.h"
#include "solver/BITL_Solver.h"
#include "solver/ITL_Preconditioner.h"
#include "solver/UmfPackSolver.h"
#include "solver/KrylovPreconditioner.h"
#include "MatrixVector.h"
#include "SystemVector.h"
#include "Lagrange.h"
#include "LeafData.h"
#include "SurfaceRegion_ED.h"
#include "ElementRegion_ED.h"
#include "DOFMatrix.h"
#include "est/Estimator.h"
#include "est/RecoveryEstimator.h"
#include "est/ResidualEstimator.h"
#include "est/SimpleResidualEstimator.h"
#include "time/RosenbrockMethod.h"
#include "nonlin/NonLinSolver.h"

#ifdef HAVE_SEQ_PETSC
  #include "solver/PetscSolver.h"
#endif

#ifdef MTL_HAS_HYPRE
  #include "solver/HypreSolver.h"
#endif

#if defined HAVE_PARALLEL_PETSC || defined HAVE_SEQ_PETSC
  #include "solver/PetscTypes.h"
#endif

namespace AMDiS {

  template<>
  void CreatorMap<LinearSolverInterface>::addDefaultCreators()
  {
    LinearSolverCreator *creator;

    // creators for MTL4 krylov solvers
    // _________________________________________________________________________
    
    creator = new CGSolver::Creator;
    addCreator("mtl_cg", creator);

    creator = new CGSSolver::Creator;
    addCreator("mtl_cgs", creator);

    creator = new BiCGSolver::Creator;
    addCreator("mtl_bicg", creator);

    creator = new BiCGStabSolver::Creator;
    addCreator("mtl_bicgstab", creator);

    creator = new BiCGStab2Solver::Creator;
    addCreator("mtl_bicgstab2", creator);

    creator = new BiCGStabEllSolver::Creator;
    addCreator("mtl_bicgstab_ell", creator);

    creator = new QMRSolver::Creator;
    addCreator("mtl_qmr", creator);

    creator = new TFQMRSolver::Creator;
    addCreator("mtl_tfqmr", creator);

    creator = new GMResSolver::Creator;
    addCreator("mtl_gmres", creator);

    creator = new GcrSolver::Creator;
    addCreator("mtl_gcr", creator);    

    creator = new FGMResSolver::Creator;
    addCreator("mtl_fgmres", creator);   
    
    creator = new IDRsSolver::Creator;
    addCreator("mtl_idr_s", creator);

    creator = new MinResSolver::Creator;
    addCreator("mtl_minres", creator);

    creator = new PreOnly::Creator;
    addCreator("mtl_preonly", creator);
    addCreator("mtl_richardson", creator);

#ifdef HAVE_UMFPACK
    creator = new UmfPackSolver::Creator;
    addCreator("mtl_umfpack", creator);
    addCreator("mtl_direct", creator);
#endif
    
#ifdef MTL_HAS_HYPRE
    creator = new HypreSolver::Creator;
    addCreator("mtl_hypre", creator);
#endif 
    
    // creators for block krylov solvers
    // _________________________________________________________________________
    
    creator = new B_CGSolver::Creator;
    addCreator("bmtl_cg", creator);

    creator = new B_CGSSolver::Creator;
    addCreator("bmtl_cgs", creator);

//     creator = new B_BiCGSolver::Creator;
//     addCreator("mtl_bicg", creator);

    creator = new B_BiCGStabSolver::Creator;
    addCreator("bmtl_bicgstab", creator);

//     creator = new B_BiCGStab2Solver::Creator;
//     addCreator("bmtl_bicgstab2", creator);      // MTL Error

//     creator = new B_BiCGStabEllSolver::Creator;
//     addCreator("bmtl_bicgstab_ell", creator);   // MTL Error

//     creator = new QMRSolver::Creator;
//     addCreator("mtl_qmr", creator);

    creator = new B_TFQMRSolver::Creator;
    addCreator("bmtl_tfqmr", creator);

    creator = new B_GMResSolver::Creator;
    addCreator("bmtl_gmres", creator);

    creator = new B_GcrSolver::Creator;
    addCreator("bmtl_gcr", creator);    

    creator = new B_FGMResSolver::Creator;
    addCreator("bmtl_fgmres", creator);   
    
//     creator = new B_IDRsSolver::Creator;
//     addCreator("bmtl_idr_s", creator);     // MTL Error

    creator = new B_MinResSolver::Creator;
    addCreator("bmtl_minres", creator);

    creator = new B_PreOnly::Creator;
    addCreator("bmtl_preonly", creator);
    addCreator("bmtl_richardson", creator);
    
    // creators for PETSc solvers
    // _________________________________________________________________________
    
#if defined HAVE_SEQ_PETSC
    // sequential PETSc-Solver
    creator = new PetscSolver::Creator;
    addCreator("petsc_petsc", creator); // standard creator for petsc solver
    addCreator("petsc", creator);
    
    LinearSolverCreator* creator2 = new PetscSolverNested::Creator;
    addCreator("bpetsc_petsc", creator2); // standard creator for petsc solver
    addCreator("bpetsc", creator2);
        
    std::map<std::string,std::string>::iterator it;
    PetscParameters params;
    for (it = params.solverMap.begin();
	 it!= params.solverMap.end();
	 it++) {
      CreatorMap< LinearSolverInterface >::addCreator("petsc_" + it->first, creator);
      CreatorMap< LinearSolverInterface >::addCreator("bpetsc_" + it->first, creator2);
    }
#endif
  }


  template<>
  void CreatorMap<ITL_PreconditionerBase< MTLTypes::MTLMatrix, MTLTypes::MTLVector > >::addDefaultCreators()
  {
    typedef CreatorInterfaceName<ITL_PreconditionerBase< MTLTypes::MTLMatrix, MTLTypes::MTLVector > > PreconditionCreator;
    PreconditionCreator *creator;

    creator =  new DiagonalPreconditioner::Creator;
    addCreator("diag", creator);

    creator =  new MassLumpingPreconditioner::Creator;
    addCreator("lumping", creator);
    addCreator("masslumping", creator);

    creator = new ILUPreconditioner::Creator;
    addCreator("ilu", creator);

    creator = new ICPreconditioner::Creator;
    addCreator("ic", creator);

    creator =  new IdentityPreconditioner::Creator;
    addCreator("no", creator);

    creator =  new KrylovPreconditionerSeq::Creator;
    addCreator("krylov", creator);
    addCreator("solver", creator);
  }


  template<>
  void CreatorMap<ITL_PreconditionerBase<BlockMTLMatrix, MTLTypes::MTLVector> >::addDefaultCreators()
  {
    addCreator("no", new BlockIdentityPreconditioner::Creator);
    addCreator("diag", new BlockDiagonalPreconditioner::Creator);
  }


#if defined HAVE_SEQ_PETSC
  template<>
  void CreatorMap<PetscPreconditioner>::addDefaultCreators() { }
  
  template<>
  void CreatorMap<PetscPreconditionerNested>::addDefaultCreators() { }
#endif

  template<>
  void CreatorMap<Estimator>::addDefaultCreators()
  {
    EstimatorCreator *creator;

    creator = new ResidualEstimator::Creator;
    addCreator("residual", creator);

    creator = new SimpleResidualEstimator::Creator;
    addCreator("simple-residual", creator);

    creator = new RecoveryEstimator::Creator;
    addCreator("recovery", creator);
  }


  template<>
  void CreatorMap<BasisFunction>::addDefaultCreators()
  {
    BasisFunctionCreator *creator;

    // Lagrange-functions up to order 4
    for (int i = 1; i <= MAX_DEGREE; i++) {
      creator = new Lagrange::Creator(i);
      addCreator("P" + std::to_string(i), creator);
      addCreator("Lagrange" + std::to_string(i), creator);
      addCreator("CG" + std::to_string(i), creator);
    }

    // linear Lagrange functions plus element bubble function
    creator = new Bubble::Creator;
    addCreator("P1+bubble", creator);
  }


  template<>
  void CreatorMap<ElementData>::addDefaultCreators()
  {
    CreatorInterface<ElementData> *creator;

    creator = new LeafDataEstimatable::Creator;
    addCreator("LeafDataEstimatable", creator);

    creator = new LeafDataEstimatableVec::Creator;
    addCreator("LeafDataEstimatableVec", creator);

    creator = new LeafDataCoarsenable::Creator;
    addCreator("LeafDataCoarsenable", creator);

    creator = new LeafDataCoarsenableVec::Creator;
    addCreator("LeafDataCoarsenableVec", creator);

    creator = new LeafDataPeriodic::Creator;
    addCreator("LeafDataPeriodic", creator);

    creator = new SurfaceRegion_ED::Creator;
    addCreator("SurfaceRegion_ED", creator);

    creator = new ElementRegion_ED::Creator;
    addCreator("ElementRegion_ED", creator);
  }


  template<>
  void CreatorMap<RosenbrockMethod>::addDefaultCreators()
  {
    addCreator("ros2", new Ros2::Creator);
    addCreator("rowda3", new Rowda3::Creator);
    addCreator("ros3p", new Ros3p::Creator);
    addCreator("rodasp", new Rodasp::Creator);
    addCreator("rosI2P1", new ROSI2P1::Creator);
    addCreator("rosI2P2", new ROSI2P2::Creator);
    addCreator("rosI2Pw", new ROSI2Pw::Creator);
    addCreator("rosI2PW", new ROSI2PW::Creator);
    addCreator("ros3Pw", new Ros3Pw::Creator);
    addCreator("ros34PW2", new Ros34PW2::Creator);
  }


  template<>
  void CreatorMap<NonLinSolver>::addDefaultCreators()
  {
    addCreator("newton", new Newton::Creator);
    addCreator("newtonArmijo", new NewtonArmijo::Creator);
  }
  
} // end namespace AMDiS
