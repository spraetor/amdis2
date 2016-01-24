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


/** \file PetscSolverFeti.h */

#include <map>
#include "AMDiS_fwd.h"
#include "parallel/MpiHelper.h"
#include "parallel/PetscSolver.h"
#include "parallel/PetscSolverFetiStructs.h"
#include "parallel/ParallelDofMapping.h"
#include "parallel/ParallelTypes.h"

#ifndef AMDIS_PETSC_SOLVER_FETI_H
#define AMDIS_PETSC_SOLVER_FETI_H

namespace AMDiS
{
  namespace Parallel
  {

    /** \ingroup Solver
     *
     * \brief
     * FETI-DP implementation based on PETSc.
     */
    class PetscSolverFeti : public PetscSolver
    {
    public:
      typedef PetscSolver super;

      /// Creator class
      class Creator : public LinearSolverCreator
      {
      public:
        virtual ~Creator() {}

        /// Returns a new PetscSolver object.
        LinearSolverInterface* create()
        {
          return new PetscSolverFeti(this->name);
        }
      };

      /// Constructor of FETI-DP solver class.
      PetscSolverFeti(std::string name);

      virtual void init(std::vector<const FiniteElemSpace*>& componentSpaces,
                        std::vector<const FiniteElemSpace*>& feSpaces,
                        bool createGlobalMapping = true);

      /// After mesh changes, or if the solver is called the first time, this
      /// function creates all information about primal nodes, dual nodes and
      /// lagrange constraints.
      void createFetiData();

      /// Assemble the sequentially created matrices to the global matrices
      /// required by the FETI-DP method.
      void fillPetscMatrix(Matrix<DOFMatrix*>* mat);

      /// Assembles the global rhs vectors from the sequentially created ones.
      void fillPetscRhs(SystemVector* vec);

      /// Solve the system using FETI-DP method.
      void solvePetscMatrix(SystemVector& vec, AdaptInfo& adaptInfo);

      /// Just for the super trick
      void solveGlobal(Vec& rhs, Vec& sol);

      /// Destroys all matrix data structures.
      void destroyMatrixData();

      /// Detroys all vector data structures.
      void destroyVectorData();

      /// Returns flags to denote which information of the boundary DOFs are
      /// required by the FETI-DP solver.
      Flag getBoundaryDofRequirement()
      {
        return
          MeshDistributor::BOUNDARY_SUBOBJ_SORTED |
          MeshDistributor::BOUNDARY_FILL_INFO_SEND_DOFS |
          MeshDistributor::BOUNDARY_FILL_INFO_RECV_DOFS;
      }

      /// Initialization of the FETI-DPdata structures.
      void initialize();

      int getNumberOfPrimals()
      {
        return primalDofMap.getOverallDofs();
      }

      int getNumberOfRankPrimals()
      {
        return primalDofMap.getRankDofs();
      }

      int getNumberOfDuals()
      {
        return dualDofMap.getOverallDofs();
      }

      int getNumberOfRankDuals()
      {
        return dualDofMap.getRankDofs();
      }

      int getNumberOfLagrange()
      {
        return lagrangeMap.getOverallDofs();
      }

    protected:
      ///
      void createDirichletData(Matrix<DOFMatrix*>& mat);

      /// Defines which boundary nodes are primal. Creates global index of
      /// the primal variables.
      void createPrimals(int component);

      /// Defines the set of dual variables and creates the global index of
      /// dual variables.
      void createDuals(int component);

      ///
      void createInterfaceNodes(int component);

      /// Create Lagrange multiplier variables corresponding to the dual
      /// variables.
      void createLagrange(int component);

      void createAugmentedLagrange(int component);

      /// Creates a global index of the B variables.
      void createIndexB(int component);

      /// Creates the Lagrange multiplier constraints and assembles them
      /// to \ref mat_lagrange.
      void createMatLagrange();

      std::vector<std::vector<BoundaryObject>> getCoarseEdges();

      std::vector<std::vector<BoundaryObject>> getCoarseFaces();

      void createMatAugmentedLagrange();

      bool testWirebasketEdge(BoundaryObject& edge,
                              const FiniteElemSpace* feSpace);

      ///
      void createPreconditionerMatrix(Matrix<DOFMatrix*>* mat);

      /// Creates PETSc KSP solver object for solving the Schur complement
      /// system on the primal variables, \ref ksp_schur_primal
      void createSchurPrimalKsp();

      ///
      void createMatExplicitSchurPrimal();

      ///
      void createMatExplicitAugmentedSchurPrimal();

      /// Destroys PETSc KSP solver object \ref ksp_schur_primal
      void destroySchurPrimalKsp();

      /// Creates PETSc KSP solver object for the FETI-DP operator, \ref ksp_feti
      void createFetiKsp();

      ///
      void createFetiExactKsp();

      ///
      void createFetiInexactKsp();

      ///
      void createFetiInexactReducedKsp();

      ///
      void createFetiPreconLumped(PC pc);

      ///
      void createFetiPreconDirichlet(PC pc);

      /// Destroys FETI-DP operator, \ref ksp_feti
      void destroyFetiKsp();

      ///
      void destroyFetiExactKsp();

      ///
      void destroyFetiInexactKsp();

      ///
      void destroyFetiInexactReducedKsp();

      /// Create the null space of the FETI-DP operator (if there is one) and
      /// attachets it to the corresponding matrices and KSP objects.
      void createNullSpace();

      /// In debug modes, this function runs some debug tests on the FETI
      /// matrices. In optimized mode, nothing is done here.
      void dbgMatrix(Matrix<DOFMatrix*>* mat);

      /** \brief
      * Recovers AMDiS solution vector from PETSc's solution vectors of the
      * FETI-DP system. First, the B variables can locally be copied to the
      * corresponding entries in the DOF vectors. The primal variable must
      * be communicated such that all ranks sharing a primal get a copy of
      * the corresponding value.
      *
      * \param[in]   vec_sol_b        Global PETSc vector of the solution of
      *                               the B variables.
      * \param[in]   vec_sol_primal   Global PETSc vector of the solution of
      *                               the primal variables.
      * \param[out]  vec              SystemVector containing all solution
      *                               DOF vectors.
      */
      void recoverSolution(Vec& vec_sol_b,
                           Vec& vec_sol_primal,
                           SystemVector& vec);

      ///
      void recoverInterfaceSolution(Vec& vecInterface,
                                    SystemVector& vec);

      ///
      void solveFeti(Vec& rhsInterior, Vec& rhsCoarse,
                     Vec& solInterior, Vec& solCoarse);

      ///
      void solveFetiExact(Vec& rhsInterior, Vec& rhsCoarse,
                          Vec& solInterior, Vec& solCoarse);

      ///
      void solveFetiInexact(Vec& rhsInterior, Vec& rhsCoarse,
                            Vec& solInterior, Vec& solCoarse);

      ///
      void solveFetiInexactReduced(Vec& rhsInterior, Vec& rhsCoarse,
                                   Vec& solInterior, Vec& solCoarse);

      ///
      void resetStatistics();

      ///
      void printStatistics();

      /// Checks whether a given DOF is a primal DOF in a given component.
      inline bool isPrimal(int component, DegreeOfFreedom dof)
      {
        return primalDofMap[component].isSet(dof);
      }

      /// Checks whether a given DOF is a dual DOF in a given component.
      inline bool isDual(int component, DegreeOfFreedom dof)
      {
        return dualDofMap[component].isSet(dof);
      }

      /// Checks whether a given DOF is an interface DOF in a given component.
      inline bool isInterface(int component, DegreeOfFreedom dof)
      {
        if (component == pressureComponent)
          return interfaceDofMap[component].isSet(dof);

        return false;
      }

    protected:
      /// Type of FETI-DP solver, i.e., exact or some inexact version
      FetiSolverType fetiSolverType;

      ///
      ParallelDofMapping dofMapSubDomain;

      /// Mapping from primal DOF indices to a global index of primals.
      ParallelDofMapping primalDofMap;

      /// Mapping from dual DOF indices to a global index of duals.
      ParallelDofMapping dualDofMap;

      /// Mapping from interface DOF indices to a global index of interface
      /// nodes. This is mainly used for Stokes-like solvers, where the pressure
      /// interface nodes are neither primal nor dual.
      ParallelDofMapping interfaceDofMap;

      /// Index for each non primal DOF to the global index of B variables (thus,
      /// all pure local variables).
      ParallelDofMapping localDofMap;

      /// Stores to each dual DOF index the index of the first Lagrange
      /// constraint that is assigned to this DOF.
      ParallelDofMapping lagrangeMap;

      /// Mapping of pure local DOF indices, thus no primal and no dual DOFs are
      /// in this map. Is used for the Dirichlet preconditioner only.
      ParallelDofMapping interiorDofMap;

      /// Stores to all dual boundary DOFs in each FE space the set of
      /// ranks which contain this global DOF.
      std::map<const FiniteElemSpace*, DofIndexToPartitions> boundaryDofRanks;

      /// Global PETSc matrix of Lagrange variables.
      Mat mat_lagrange;

      ///
      Mat mat_augmented_lagrange;

      /// 0: Solve the Schur complement on primal variables with iterative solver.
      /// 1: Create the Schur complement matrix explicitly and solve it with a
      ///    direct solver.
      int schurPrimalSolver;

      /// PETSc solver object to solve the Schur complement on the
      /// primal variables.
      KSP ksp_schur_primal;

      /// Matrix object that defines a matrix-free implementation for the action
      /// of the Schur complement on the primal variables.
      Mat mat_schur_primal;

      /// Data for MatMult operation in matrix \ref mat_schur_primal
      SchurPrimalData schurPrimalData;

      ///
      SchurPrimalAugmentedData schurPrimalAugmentedData;

      /// PETSc solver object to solve a system with FETI-DP.
      KSP ksp_feti;

      /// Matrix object that defines a matrix-free implementation for the action
      /// of the FETI-DP operator.
      Mat mat_feti;

      /// Data for MatMult operation in matrix \ref mat_feti
      FetiData fetiData;

      ///
      FetiInexactData fetiInexactData;

      ///
      FetiInexactPreconData fetiInexactPreconData;

      /// Defines which preconditioner should be used to solve the reduced
      /// FETI-DP system.
      FetiPreconditioner fetiPreconditioner;

      /// Preconditioner object for the reduced FETI-DP system.
      PC precon_feti;

      Mat mat_lagrange_scaled;

      FetiDirichletPreconData fetiDirichletPreconData;

      FetiLumpedPreconData fetiLumpedPreconData;

      FetiInterfaceLumpedPreconData fetiInterfaceLumpedPreconData;

      FetiKspData fetiKspData;

      /// Matrices for Dirichlet preconditioner.
      Mat mat_interior_interior, mat_duals_duals, mat_interior_duals, mat_duals_interior;

      KSP ksp_interior;

      int levelMode;

      // If true, FETI-DP subdomains are on MPI::COMM_SELF
      bool subDomainIsLocal;

      PetscSolver* subdomain;

      // Just a trick for multi level things, should be removed and generalized!
      PetscSolver* mlSubdomain;

      PetscSolver* massMatrixSolver;

      bool printTimings;

      bool augmentedLagrange;

      int nRankEdges;

      int nOverallEdges;

      /// There are two different dirichlet modes:
      ///   0: dirichlet rows are zeroed and a diagonal element is set to one.
      ///   1: dirichlet rows are removed (this mode does not work correctly, but
      ///      many function are prepered to make use of it)
      int dirichletMode;

      /// If true, the FETI-DP solver is applied to a Stokes like problem. Thus,
      /// there is a pressure variable which is not part of the coarse grid
      /// problem.
      bool stokesMode;

      /// Only used if \ref stokesMode is enabled. In this case, this variable
      /// defines the component number of the pressure variable.
      int pressureComponent;

      /// Maps from component number to set of DOFs which are Dirichlet DOfs in
      /// this component.
      std::map<int, std::set<DegreeOfFreedom>> dirichletRows;

      friend class PetscSolverFetiDebug;
    };
  } // end namespace Parallel
} // end namespace AMDiS

#endif
