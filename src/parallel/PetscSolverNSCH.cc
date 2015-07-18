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



#include "parallel/PetscSolverNSCH.h"
#include "parallel/PetscHelper.h"
#include "TransformDOF.h"




namespace AMDiS { namespace Parallel {
  
  using namespace std;
  
  
  PetscErrorCode pcShell(PC pc, Vec b, Vec x) // solve Px=b
  {FUNCNAME("PCApply()");
  
  void *ctx;
  PCShellGetContext(pc, &ctx);
  CahnHilliardData2* data = static_cast<CahnHilliardData2*>(ctx);
  
  /// extract vectors
  Vec b1, b2, b34, b5, x1, x2, x34, x5;
  data->globalMatrixSolver->extractVectorComponent(b, data->dim+1, &b1);
  data->globalMatrixSolver->extractVectorComponent(b, data->dim+2, &b2);
  data->globalMatrixSolver->extractVectorComponent(b, 0, &b34, data->dim);  
  data->globalMatrixSolver->extractVectorComponent(b, data->dim, &b5);
  
  data->globalMatrixSolver->extractVectorComponent(x, data->dim+1, &x1);
  data->globalMatrixSolver->extractVectorComponent(x, data->dim+2, &x2);
  data->globalMatrixSolver->extractVectorComponent(x, 0, &x34, data->dim);  
  data->globalMatrixSolver->extractVectorComponent(x, data->dim, &x5);
    
  
  /// solve Cahn-Hilliard Preconditioner  
    Vec y1, y2;
     VecDuplicate(b1, &y1);
     VecDuplicate(b2, &y2);
    
     KSPSolve(data->kspMassCH, b1, y1); 			// M*y1 = b1    
     MatMultAdd(data->matMinusDeltaK, y1, b2, x1); 	// -> x1 := b2-delta*K*y1
     
     KSPSolve(data->kspLaplaceCH, x1, y2); 			// (M+eps*sqrt(delta))*y2 = x1
     MatMult(data->matMassCH, y2, x1); 			// x1 := M*y2
     
     KSPSolve(data->kspLaplaceCH, x1, x2);			// (M+eps*sqrt(delta))*x2 = x1
     double factor = (*data->eps)/sqrt(*data->delta);
     VecCopy(x2, x1); 					// x1 := x2
     VecAXPBYPCZ(x1, 1.0, factor, -factor, y1, y2);	// x1 = 1*y1 + factor*y2 - factor*x1
     
     VecDestroy(&y1);
     VecDestroy(&y2);     
    /**/ 
     
  /// solve Navier-Stokes Preconditioner  
  Vec tmp34, tmp5, tmp5_2, tmp34_2;
  VecDuplicate(b34, &tmp34);
  VecDuplicate(b34, &tmp34_2);
  VecDuplicate(b5, &tmp5);
  VecDuplicate(b5, &tmp5_2);
  
  KSPSolve(data->kspVelocity, b34, tmp34);
  VecScale(tmp34, -1.0);
  MatMultAdd(data->matDiv, tmp34, b5, tmp5);  
  
  /// approximierte Schur Komplement
   KSPSolve(data->kspLaplace, tmp5, x5);
   { //project out constant Null-space
      int vecSize;
      VecGetSize(x5, &vecSize);
      PetscScalar vecSum;
      VecSum(x5, &vecSum);
      vecSum = vecSum / static_cast<PetscScalar>(-1.0 * vecSize);
      VecShift(x5, vecSum); 
      //VecView(y, PETSC_VIEWER_STDOUT_WORLD);
    }
    MatMult(data->matConDif, x5, tmp5);
    KSPSolve(data->kspMass, tmp5, x5);
    VecScale(x5,-1.0);
    
  MatMultAdd(data->matGrad, x5, b34, tmp34); 
  KSPSolve(data->kspVelocity, tmp34, x34);
  VecScale(x5,-1.0);
/**/  

  VecDestroy(&tmp34);
  VecDestroy(&tmp5);
  VecDestroy(&tmp5_2);
  VecDestroy(&tmp34_2);
  VecDestroy(&b1);
  VecDestroy(&b2);
  VecDestroy(&b34);
  VecDestroy(&b5);
  VecDestroy(&x1);
  VecDestroy(&x2);
  VecDestroy(&x34);
  VecDestroy(&x5);
  
    PetscFunctionReturn(0);
  }
  
  
  
  PetscSolverNSCH::PetscSolverNSCH(string name)
  : PetscSolverGlobalMatrix(name),      
  pressureNullSpace(true),
  useOldInitialGuess(false),      
  velocitySolutionMode(0),
  massSolutionMode(0),
  laplaceSolutionMode(0),
  regularizeLaplace(0),
  massMatrixSolverCH(NULL),
  laplaceMatrixSolverCH(NULL), 
  deltaKMatrixSolver(NULL),
  massMatrixSolver(NULL),
  laplaceMatrixSolver(NULL),
  conDifMatrixSolver(NULL),
  nu(NULL),
  invTau(NULL),
  solution(NULL),
  phase(NULL)
  {    
    Parameters::get(initFileStr + "->use old initial guess", 
		    useOldInitialGuess);
    
    Parameters::get(initFileStr + "->navierstokes->velocity solver", 
		    velocitySolutionMode);
    
    Parameters::get(initFileStr + "->navierstokes->mass solver", 
		    massSolutionMode);
    
    Parameters::get(initFileStr + "->navierstokes->laplace solver", 
		    laplaceSolutionMode);
    Parameters::get(initFileStr + "->navierstokes->regularize laplace", 
		    regularizeLaplace);
    
    
  }
  
  
  void PetscSolverNSCH::solvePetscMatrix(SystemVector &vec, 
						  AdaptInfo *adaptInfo)
  {
    FUNCNAME("PetscSolverNSCH::solvePetscMatrix()");
    
       if (useOldInitialGuess) {      
          VecSet(getVecSolInterior(), 0.0);
          
          for (int i = 0; i < solution->getSize(); i++)
    	setDofVector(getVecSolInterior(), solution->getDOFVector(i), i, true);
          
          vecSolAssembly();
          KSPSetInitialGuessNonzero(kspInterior, PETSC_TRUE);
  }
  
    PetscSolverGlobalMatrix::solvePetscMatrix(vec, adaptInfo);
  }
  
  
  void PetscSolverNSCH::initSolver(KSP &ksp)
  {
    FUNCNAME("PetscSolverNSCH::initSolver()");
    
    // Create FGMRES based outer solver    
    MSG("CREATE POS 1: %p\n", &ksp);
    KSPCreate(domainComm, &ksp);
#if (PETSC_VERSION_MINOR >= 5)
    KSPSetOperators(ksp, getMatInterior(), getMatInterior());
#else
    KSPSetOperators(ksp, getMatInterior(), getMatInterior(), SAME_NONZERO_PATTERN);
#endif
    KSPMonitorSet(ksp, KSPMonitorTrueResidualNorm, PETSC_NULL, PETSC_NULL);
    petsc_helper::setSolver(ksp, "ch_", KSPFGMRES, PCSHELL, getRelative(), getTolerance(), getMaxIterations());
    setConstantNullSpace(ksp, componentSpaces[0]->getMesh()->getDim() , true);
  }
  
  
  void PetscSolverNSCH::initPreconditioner(PC pc)
  {
    FUNCNAME("PetscSolverNSCH::initPreconditioner()");
       
    MPI::COMM_WORLD.Barrier();
    double wtime = MPI::Wtime();
    int dim = componentSpaces[0]->getMesh()->getDim();     
    pressureComponent=dim;
    const FiniteElemSpace *cahnHilliardFeSpace = componentSpaces[dim+1];
//     const FiniteElemSpace *velocityFeSpace= componentSpaces[0];
    const FiniteElemSpace *pressureFeSpace = componentSpaces[pressureComponent];
        
    PCSetType(pc, PCSHELL);
    PCShellSetApply(pc, pcShell);
    PCShellSetContext(pc, &matShellContext);
    matShellContext.globalMatrixSolver = (this);
    matShellContext.mpiCommGlobal= &(meshDistributor->getMpiComm(0));
    matShellContext.dim = dim;
    
    
    /// Init Cahn-Hilliard part
    
    DOFMatrix laplaceMatrixCH(cahnHilliardFeSpace, cahnHilliardFeSpace);
    Operator laplaceOpCH(cahnHilliardFeSpace, cahnHilliardFeSpace);
    
    DOFMatrix massMatrixCH(cahnHilliardFeSpace, cahnHilliardFeSpace);
    Operator massOpCH(cahnHilliardFeSpace, cahnHilliardFeSpace);
    
    DOFMatrix deltaKMatrix(cahnHilliardFeSpace, cahnHilliardFeSpace);
    Operator laplaceOp2(cahnHilliardFeSpace, cahnHilliardFeSpace);
    
    { 
      Simple_ZOT zot;
      massOpCH.addTerm(&zot);
      laplaceOpCH.addTerm(&zot); // M
      Simple_SOT sot2((*eps)*sqrt(*delta));
      laplaceOpCH.addTerm(&sot2); // eps*sqrt(delta)*K
      Simple_SOT sot(-(*delta));    
      laplaceOp2.addTerm(&sot); // -delta*K
      
      massMatrixCH.assembleOperator(massOpCH);        
      massMatrixSolverCH = createSubSolver(dim+1, "mass_");       
      massMatrixSolverCH->fillPetscMatrix(&massMatrixCH);     
      // === matrix (M + eps*sqrt(delta)*K) ===
      laplaceMatrixCH.assembleOperator(laplaceOpCH);
      laplaceMatrixSolverCH = createSubSolver(dim+1, "laplace_");
      laplaceMatrixSolverCH->fillPetscMatrix(&laplaceMatrixCH);     
      
      // === matrix (-delta*K) ===
      deltaKMatrix.assembleOperator(laplaceOp2);
      deltaKMatrixSolver = createSubSolver(dim+1, "laplace2_");
      deltaKMatrixSolver->fillPetscMatrix(&deltaKMatrix);      
    }
    
    // === Setup solver ===
    matShellContext.kspMassCH = massMatrixSolverCH->getSolver();
    matShellContext.kspLaplaceCH = laplaceMatrixSolverCH->getSolver();
    matShellContext.matMassCH = massMatrixSolverCH->getMatInterior();    
    matShellContext.matMinusDeltaK = deltaKMatrixSolver->getMatInterior();   
    matShellContext.eps = eps;
    matShellContext.delta = delta;
    
    petsc_helper::setSolver(matShellContext.kspMassCH, "", KSPCG, PCJACOBI, 0.0, 1e-14, 2);
    petsc_helper::setSolver(matShellContext.kspLaplaceCH, "laplace_", KSPRICHARDSON, PCHYPRE, 0.0, 1e-14, 1);
  
    
    /// Init Navier-Stokes part
    
    DOFMatrix massMatrix(pressureFeSpace, pressureFeSpace);
      Operator massOp(pressureFeSpace, pressureFeSpace);
      ZeroOrderTerm *massTerm = new Simple_ZOT;
      massOp.addTerm(massTerm);
      massMatrix.assembleOperator(massOp);
      delete massTerm;
    massMatrixSolver = createSubSolver(pressureComponent, "mass_");
    massMatrixSolver->fillPetscMatrix(&massMatrix);
    matShellContext.kspMass = massMatrixSolver->getSolver();   
    
    DOFMatrix laplaceMatrix(pressureFeSpace, pressureFeSpace);    
      Operator laplaceOp(pressureFeSpace, pressureFeSpace);
      SecondOrderTerm *laplaceTerm = new Simple_SOT;
      laplaceOp.addTerm(laplaceTerm);
      laplaceMatrix.assembleOperator(laplaceOp);
      delete laplaceTerm;    
    laplaceMatrixSolver = createSubSolver(pressureComponent, string("laplace_"));
    laplaceMatrixSolver->fillPetscMatrix(&laplaceMatrix);
        
    // === Create convection-diffusion operator ===
    
    DOFVector<double> vx(pressureFeSpace, "vx");
    DOFVector<double> vy(pressureFeSpace, "vy");
    DOFVector<double> vz(pressureFeSpace, "vz");
    DOFVector<double> vp(pressureFeSpace, "vp");
    vx.interpol(solution->getDOFVector(0));
    vy.interpol(solution->getDOFVector(1));
    if (dim == 3) vz.interpol(solution->getDOFVector(2));
    DOFMatrix conDifMatrix(pressureFeSpace, pressureFeSpace);
    {
      Operator conDifOp(pressureFeSpace, pressureFeSpace);
      ZeroOrderTerm *conDif0 = NULL;
      SecondOrderTerm *conDif1 = NULL;
      FirstOrderTerm *conDif2 = NULL, *conDif3 = NULL, *conDif4 = NULL;
      vp.interpol(solution->getDOFVector(dim+2));
      
      densityFunctionTau = new LinearInterpolation(*rho1,*rho2,*invTau);
      viscosityFunction = new LinearInterpolation(*nu1,*nu2);
      densityFunction = new LinearInterpolation2(*rho1,*rho2);      
      
      conDif0 = new VecAtQP_ZOT(&vp, densityFunctionTau);
      conDifOp.addTerm(conDif0);	
      conDif1 = new VecAtQP_SOT(&vp, viscosityFunction );
      conDifOp.addTerm(conDif1);	
      conDif2 = new Vec2AtQP_FOT(&vx, &vp, densityFunction, 0);
      conDifOp.addTerm(conDif2, GRD_PHI);
      
      conDif3 = new Vec2AtQP_FOT(&vy, &vp, densityFunction, 1);
      conDifOp.addTerm(conDif3, GRD_PHI);
      
      if (dim == 3) {
	conDif4 = new Vec2AtQP_FOT(&vz, &vp, densityFunction, 2);
	conDifOp.addTerm(conDif4, GRD_PHI);
      }
         
      conDifMatrix.assembleOperator(conDifOp);
      
      delete conDif0;
      delete conDif1;
      delete conDif2;
      delete conDif3;
      if (dim == 3) delete conDif4;      
    }
    
    conDifMatrixSolver = createSubSolver(pressureComponent, "condif_");
    conDifMatrixSolver->fillPetscMatrix(&conDifMatrix);
    matShellContext.matConDif = conDifMatrixSolver->getMatInterior();    
    
    extractMatrixComponent(mat[0][0], dim, 1, 0, dim, &(matShellContext.matDiv));
    extractMatrixComponent(mat[0][0], 0, dim, dim, 1, &(matShellContext.matGrad));
    extractMatrixComponent(mat[0][0], 0, dim, 0, dim, &(matShellContext.velocityMat));
    
    ///erstelle kspVelocity
    KSPCreate((meshDistributor->getMpiComm(0)), &(matShellContext.kspVelocity));
    
#if (PETSC_VERSION_MINOR >= 5)
    KSPSetOperators(matShellContext.kspVelocity, matShellContext.velocityMat, matShellContext.velocityMat);
#else
    KSPSetOperators(matShellContext.kspVelocity, matShellContext.velocityMat, matShellContext.velocityMat, SAME_NONZERO_PATTERN);
#endif
    
    ///regularisiere LaplaceMatrix
  if (regularizeLaplace)
  {
    PetscInt rows[1];
    rows[0]=0;
    MatZeroRows(laplaceMatrixSolver->getMatInterior(), 1, rows, 0, PETSC_NULL, PETSC_NULL);
    KSPCreate((meshDistributor->getMpiComm(0)), &(matShellContext.kspLaplace));
#if (PETSC_VERSION_MINOR >= 5)
    KSPSetOperators(matShellContext.kspLaplace, laplaceMatrixSolver->getMatInterior(), laplaceMatrixSolver->getMatInterior());
#else
    KSPSetOperators(matShellContext.kspLaplace, laplaceMatrixSolver->getMatInterior(), laplaceMatrixSolver->getMatInterior(), SAME_NONZERO_PATTERN);
#endif
  }
   else
   {  matShellContext.kspLaplace=laplaceMatrixSolver->getSolver();
     setConstantNullSpace(matShellContext.kspLaplace);
   }
    
  
    
    // === Setup solver ===
    
    switch (massSolutionMode) {
      case 0:
	petsc_helper::setSolver(matShellContext.kspMass, "mass_", KSPCG, PCJACOBI, 0.0, 1e-14, 2);
	break;
      case 1:
	petsc_helper::setSolverWithLu(matShellContext.kspMass, "mass_", KSPRICHARDSON,     PCLU, MATSOLVERMUMPS, 0.0, 1e-14, 1);
	break;
      default:
	ERROR_EXIT("No mass solution mode %d available!\n", massSolutionMode);
    }
    
    switch (laplaceSolutionMode) {
      case 0:
	petsc_helper::setSolver(matShellContext.kspLaplace, "laplace_", KSPRICHARDSON, PCHYPRE, 0.0, 1e-14, 1);
	break;
      case 1:
	petsc_helper::setSolverWithLu(matShellContext.kspLaplace, "mass_",   KSPRICHARDSON, PCLU, MATSOLVERMUMPS,   0.0, 1e-14, 1);
	break;
      default:
	ERROR_EXIT("No laplace solution mode %d available!\n", laplaceSolutionMode);
    }
    
    switch (velocitySolutionMode) {
      case 0:      
	petsc_helper::setSolver(matShellContext.kspVelocity, "", KSPRICHARDSON, PCHYPRE, 0.0, 1e-14, 1);
	break;
      case 1:
	petsc_helper::setSolverWithLu(matShellContext.kspVelocity, "", KSPPREONLY,   PCLU, MATSOLVERMUMPS , 0.0, 1e-14, 1);
	break;
      default:
	ERROR_EXIT("No velocity solution mode %d available!\n", velocitySolutionMode);
    }
    
 
    PCSetFromOptions(pc);
    
    MSG("Setup of Cahn-Hilliard preconditioner needed %.5f seconds\n", 
	MPI::Wtime() - wtime);
     }
     
     
     void PetscSolverNSCH::exitPreconditioner(PC pc)
     {
       FUNCNAME("PetscSolverNSCH::exitPreconditioner()");
       
       massMatrixSolverCH->destroyMatrixData();
       massMatrixSolverCH->destroyVectorData();
       laplaceMatrixSolverCH->destroyMatrixData();
       laplaceMatrixSolverCH->destroyVectorData();    
       deltaKMatrixSolver->destroyMatrixData();
       deltaKMatrixSolver->destroyVectorData();
       
       
       massMatrixSolverCH->destroyVectorData();
       laplaceMatrixSolverCH->destroyVectorData();
       deltaKMatrixSolver->destroyVectorData();
       
       delete massMatrixSolverCH;
       massMatrixSolverCH = NULL;
       
       delete laplaceMatrixSolverCH;
       laplaceMatrixSolverCH = NULL;
       
       delete deltaKMatrixSolver;
       deltaKMatrixSolver = NULL;
       
       massMatrixSolver->destroyMatrixData();
       massMatrixSolver->destroyVectorData();
       laplaceMatrixSolver->destroyMatrixData();
       laplaceMatrixSolver->destroyVectorData();
       conDifMatrixSolver->destroyMatrixData();
       conDifMatrixSolver->destroyVectorData();
       
       massMatrixSolver->destroyVectorData();
       laplaceMatrixSolver->destroyVectorData();
       conDifMatrixSolver->destroyVectorData();
       
       
       delete massMatrixSolver;
       massMatrixSolver = NULL;
       
       delete laplaceMatrixSolver;
       laplaceMatrixSolver = NULL;
       
       delete conDifMatrixSolver;
       conDifMatrixSolver = NULL;
       
       KSPDestroy(&(matShellContext.kspVelocity));
       if (regularizeLaplace)
        KSPDestroy(&(matShellContext.kspLaplace));
       
       MatDestroy(&(matShellContext.matGrad));
       MatDestroy(&(matShellContext.matDiv));
       MatDestroy(&(matShellContext.velocityMat));
       
       delete densityFunction;
       delete densityFunctionTau;
       delete viscosityFunction;
     }
  } }
  
