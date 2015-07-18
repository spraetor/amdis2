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



/** \file ParallelCoarseSpaceSolver.h */

#ifndef AMDIS_PARALLEL_COARSE_SPACE_SOLVER_H
#define AMDIS_PARALLEL_COARSE_SPACE_SOLVER_H

#include <mpi.h>
#include <vector>
#include <map>
#include <petsc.h>
#include "AMDiS_fwd.h"
#include "parallel/ParallelDofMapping.h"
#include "parallel/MeshDistributor.h"
#include "parallel/MatrixNnzStructure.h"
#include "parallel/ParallelSolver.h"
#include "solver/LinearSolverInterface.h"

namespace AMDiS
{
  namespace Parallel
  {
    
    /** \ingroup Solver
    * 
    * \brief
    * This class implements a block structured PETSc matrix/vec which seperates
    * the discretization of the interior of subdomains and the discretization 
    * of the coarse space. Thus, we have one matrix block for the interior and
    * one matrix block for the coarse space plus the coupling blocks. Some notes:
    * - For a single level domain decomposition method (e.g. the standad
    *   FETI-DP method), the interior matrix is local to the current rank and the
    *   coarse space matrix is a globally distributed matrix.
    * - There are different coarse spaces for different components possible. In
    *   this case, there are as many blocks as there are different coarse spaces
    *   plus one block for the interior matrix.
    * - This class also manages the creation of the corresponding non zero 
    *   structure of the matrices.
    */
    class ParallelCoarseSpaceSolver : public ParallelSolver
    {      
    public:
      /// Constructor
      ParallelCoarseSpaceSolver(std::string name);
      
      /// Destructor
      virtual ~ParallelCoarseSpaceSolver() {};
      
      /// Specialization of ParallelSolver::init()
      virtual void init(std::vector<const FiniteElemSpace*> &componentSpaces,
			std::vector<const FiniteElemSpace*> &feSpaces,
			bool createGlobalMapping = true);
      
      /** \brief
      * Sets the coarse space for all or a specific component.
      * 
      * \param[in]  coarseDofs  Coarse space DOF mapping.
      * \param[in]  component   If the standard value -1 is used, the coarse
      *                         space DOF mapping is set for all components
      *                         of the equation. Otherwise, the coarse space
      *                         DOF mapping is set only for the given one.
      */
      void setCoarseSpaceDofMapping(ParallelDofMapping *coarseDofs, 
				    int component = -1);
      

      /// Create a parallel distributed PETSc vector based on this mapping.
      void createVec(ParallelDofMapping& map, Vec &vec)
      {
	VecCreateMPI(map.getMpiComm(), map.getRankDofs(), map.getOverallDofs(), &vec);
      }

      /// Create a parallel distributed PETsc vector based on this mapping but
      /// with a different (larger) global size. This is used in multi-level 
      /// method to embed a local vector into a subdomain spaned by several
      /// ranks.
      void createVec(ParallelDofMapping& map, Vec &vec, int nGlobalRows)
      {
	VecCreateMPI(map.getMpiComm(), map.getRankDofs(), nGlobalRows, &vec);
      }

      void createLocalVec(ParallelDofMapping& map, Vec &vec)
      {
	VecCreateSeq(PETSC_COMM_SELF, map.getRankDofs(), &vec);
      }
      
      /// Creates matrices and vectors with respect to the coarse space.
      void createMatVec(Matrix<DOFMatrix*>& seqMat);

      /// Run PETSc's matrix assembly routines.
      void matAssembly();

      /// Run PETSc's vector assembly routines on rhs vectors.
      void vecRhsAssembly();

      /// Run PETSc's vector assembly routines on solution vectors.
      void vecSolAssembly();

      /// Destroys PETSc matrix objects.
      void matDestroy();

      /// Destroys PETSc vector objects.
      void vecDestroy();


      /// Just for super trick
      std::vector<std::vector<Mat> >& getMat()
      {
	return mat;
      }

      /// Just for super trick
      std::vector<Vec>& getVecRhs()
      {
	return vecRhs;
      }


      /// Get interior matrix.
      inline Mat& getMatInterior()
      {
	TEST_EXIT_DBG(mat.size() > 0)("No matrix data!\n");
	return mat[0][0];
      }

      /// Get coarse space matrix.
      inline Mat& getMatCoarse(int coarseSpace0 = 0, int coarseSpace1 = 0)
      {
	TEST_EXIT_DBG(static_cast<int>(mat.size()) > coarseSpace0 + 1)("No matrix data!\n");
	TEST_EXIT_DBG(static_cast<int>(mat.size()) > coarseSpace1 + 1)("No matrix data!\n");
	return mat[coarseSpace0 + 1][coarseSpace1 + 1];
      }

      /// Get coupling matrix of the interior and some coarse space.
      inline Mat& getMatInteriorCoarse(int coarseSpace = 0)
      {
	TEST_EXIT_DBG(static_cast<int>(mat.size()) > coarseSpace + 1)("No matrix data!\n");
	return mat[0][coarseSpace + 1];
      }

      /// Get coupling of some coarse space matrix and the interior.
      inline Mat& getMatCoarseInterior(int coarseSpace = 0)
      {
	TEST_EXIT_DBG(static_cast<int>(mat.size()) > coarseSpace + 1)("No matrix data!\n");
	return mat[coarseSpace + 1][0];
      }

      /// Get the coarse space matrix of some system component.
      inline Mat& getMatCoarseByComponent(int rowComp, int colComp = -1)
      {
	int rowMatIndex = componentIthCoarseMap[rowComp] + 1;
	int colMatIndex = componentIthCoarseMap[(colComp == -1 ? rowComp : colComp)] + 1;
	return mat[rowMatIndex][colMatIndex];
      }

      /// Get coupling matrix of the interior and the coarse space of a 
      /// system component.
      inline Mat& getMatInteriorCoarseByComponent(int comp)
      {
	int matIndex = componentIthCoarseMap[comp] + 1;
	return mat[0][matIndex];
      }

      /// Get coupling matrix of the coarse space of a system component and the
      /// interior matrix.
      inline Mat& getMatCoarseInteriorByComponent(int comp)
      {
	int matIndex = componentIthCoarseMap[comp] + 1;
	return mat[matIndex][0];
      }

      /// Get the RHS vector of the interior.
      inline Vec& getVecRhsInterior()
      {
	return vecRhs[0];
      }

      /// Get the RHS vector of some coarse space.
      inline Vec& getVecRhsCoarse(int coarseSpace = 0)
      {
	return vecRhs[coarseSpace + 1];
      }

      /// Get the RHS vector of the coarse space of a system component.
      inline Vec& getVecRhsCoarseByComponent(int comp)
      {
	int vecIndex = componentIthCoarseMap[comp] + 1;
	return vecRhs[vecIndex];
      }

      /// Get the solution vector of the interior.
      inline Vec& getVecSolInterior()
      {
	return vecSol[0];
      }

      /// Get the solution vector of some coarse space.
      inline Vec& getVecSolCoarse(int coarseSpace = 0)
      {
	FUNCNAME("ParallelCoarseSpaceSolver::getVecSolCoarse()");

	TEST_EXIT_DBG(coarseSpace + 1 < static_cast<int>(vecSol.size()))
	  ("Wrong component %d, vecSol has only %d!\n", coarseSpace + 1, vecSol.size());

	return vecSol[coarseSpace + 1];
      }

      /// Get the solution vector of the coarse space of a system component.
      inline Vec& getVecSolCoarseByComponent(int comp)
      {
	int vecIndex = componentIthCoarseMap[comp] + 1;
	return vecSol[vecIndex];
      }
      
      /** \brief
      * Checks whether a given DOF index in some component is a coarse space DOF.
      * Note (TODO): The specification of both, the component number and FE 
      * space is not really necessary. Rewrite this!
      *
      * \param[in]  component  Component number of the system.
      * \param[in]  dof        DOF index
      *
      * \return     True, if the dof is a coarse space DOF in the component. 
      *             False otherwise.
      */
      inline bool isCoarseSpace(int component,
				DegreeOfFreedom dof)
      {
	FUNCNAME("ParallelCoarseSpaceSolver::isCoarseSpace()");
	
	if (coarseSpaceMap.empty())
	  return false;

	TEST_EXIT_DBG(coarseSpaceMap.count(component))
	  ("Component %d has no coarse space defined!\n", component);

	return (*(coarseSpaceMap[component]))[component].isSet(dof);
      }


      /// Returns whether the solver has a coarse grid.
      inline bool hasCoarseSpace() 
      {
	return (!coarseSpaceMap.empty());
      }
      
    protected:
      /// Prepare internal data structures. First, it create \ref uniqueCoarseMap
      /// and \ref componentIthCoarseMap . Both are used to create the correct 
      /// number of matrix and vectors in \ref mat and \ref vec.    
      virtual void prepare();
      
      /// Computes the values of \ref rStartInterior and 
      /// \ref nGlobalOverallInterior.
      void updateSubdomainData();

    private:
      /// Checks for mesh changes. Returns true if the mesh has been changed
      /// until the last matrix creation. Is used to rebuild matrix nz structure.
      bool checkMeshChange();

    protected:
      /// Matrix of PETSc matrices. mat[0][0] is the interior discretization
      /// matrix, mat[1][1] corresponds to the first coarse space and so on. 
      /// mat[i][j], with i not equal to j, are the coupling between the interior
      /// and the coarse space matrices, and between the different coarse spaces
      /// respectively.
      std::vector<std::vector<Mat> > mat;

      /// Solution and RHS vectors. vec[0] is the interior vector, vec[1] the 
      /// first coarse space vector and so on.
      std::vector<Vec> vecSol, vecRhs;

      /// Matrix of objects to control the matrix non zero structure of the 
      /// corresponding PETSc matrices stored in \ref mat.
      std::vector<std::vector<MatrixNnzStructure> > nnz;

      /** \brief
      * Stores for each system component (i.e. each PDE variable) the coarse
      * space that is used for its discretization.
      *
      * Example: We solve the Stokes equation in 2D with a different coarse 
      * space for the velocity unknowns (component 0 and 1) and the pressure
      * (component 2). Than:
      *    componentIthCoarseMap[0] = 0
      *    componentIthCoarseMap[1] = 0
      *    componentIthCoarseMap[2] = 1
      * The indices can directly be used to access the correspondig parallel
      * DOF mapping in \ref uniqueCoarseMap.
      */
      std::vector<int> componentIthCoarseMap;

      /// Stores a set of all coarse space DOF mapping. All entries are unique.
      std::vector<ParallelDofMapping*> uniqueCoarseMap;

      /// Stores the mesh change index of the mesh the nnz structure was created
      /// for. Therefore, if the mesh change index is higher than this value, we
      /// have to create a new nnz structure for PETSc matrices, because the mesh
      ///  has been changed and therefore also the assembled matrix structure.
      int lastMeshNnz;

      /// If this variable is set to true, the non-zero matrix structure is
      /// created each time from scratch by calling \ref createPetscNnzStrcuture.
      /// This can be necessary if the number of non-zeros in the matrix varies
      /// though the mesh does not change. This may happen if there are many
      /// operators using DOFVectors from old timestep containing many zeros due to
      /// some phase fields.
      bool alwaysCreateNnzStructure;

      /// Offset for the interior DOFs of the local interior with respect to the
      /// subdomain. In the case of a one-level method, each local interior
      /// is exactly one subdomain. In the case of a multi-level method, one
      /// subdomain may consists of several rank domains. This value defines than
      /// the offset ot rank's interior rows to the subdomain's interior rows.
      int rStartInterior;

      /// Number of overall rows in subdomain's interior. For one-level methods,
      /// this value is equal to the number of rows in rank's interior. See also
      /// explenation for \ref rStarInterior.
      int nGlobalOverallInterior;

      /// Parallel DOF mapping of the (optional) coarse space. Allows to define
      /// different coarse spaces for different components.
      std::map<int, ParallelDofMapping*> coarseSpaceMap;   
      
      MPI::Intracomm domainComm;

      MPI::Intracomm coarseSpaceComm;
      
      /// Prefix string for parameters in init file.
      std::string initFileStr;
    };

  } // end namespace Parallel
} // end namespace AMDiS

#endif

