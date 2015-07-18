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


/** \file ParallelSolver.h */

#ifndef AMDIS_PARALLEL_SOLVER_H
#define AMDIS_PARALLEL_SOLVER_H

#include <mpi.h>
#include <vector>
#include "solver/LinearSolverInterface.h"
#include "AMDiS_fwd.h"

namespace AMDiS
{
  namespace Parallel
  {
    
    /** \ingroup Solver
    * 
    * \brief
    * base class for all parallel solvers, e.g. \ref PetscSolver, or \ref PMTL4Solver
    */
    class ParallelSolver : public LinearSolverInterface
    {
    public: // methods
      
      /// Constructor. Parameter \param globalIndices should be true for most 
      /// \ref PetscSolver and false for \ref PMTL4Solver
      ParallelSolver(std::string name, bool globalIndices);
      
      virtual ~ParallelSolver();
      
      virtual void init(std::vector<const FiniteElemSpace*> &componentSpaces,
			std::vector<const FiniteElemSpace*> &feSpaces,
			bool createGlobalMapping = true);

      /// Set \ref ParallelDofMapping for the interior DOFs.
      void setDofMapping(ParallelDofMapping *mapping)
      {
	interiorMap = mapping;
      }

      /// Set \ref MeshDistributor object
      void setMeshDistributor(MeshDistributor *m, int level = 0)
      {
	meshDistributor = m;
	meshLevel = level;
      }

      /// Returns the \ref ParallelDofMapping for the interior DOFs.
      ParallelDofMapping* getDofMapping()
      {
	return interiorMap;
      }

      ParallelDofMapping& getDofMap()
      {
	return dofMap;
      }

      std::vector<const FiniteElemSpace*>& getComponentSpaces()
      {
	return componentSpaces;
      }
      
    protected: // members

      /// FE spaces of all components for the stationary problem the specific
      /// solver object is registered to.
      std::vector<const FiniteElemSpace*> componentSpaces;

      /// Set of unique FE spaces in \ref componentSpaces.
      std::vector<const FiniteElemSpace*> feSpaces;

      /// Parallel DOF mapping for the interior.
      ParallelDofMapping *interiorMap;
      
      /// ???
      ParallelDofMapping dofMap;

      /// If the parallel DOF mappaings of this solver are registered to the
      /// mesh distributor object, this variable is set to true to remove them
      /// in the destructor.
      bool parallelDofMappingsRegistered;
      
      /// Pointer to a mesh distributor object.
      MeshDistributor *meshDistributor;
      
      /// Level of subdomain/interior discretization. Is used for multi-level
      /// methods only.
      int meshLevel;
    };
    
  } // end namespace Parallel
} // end namespace AMDiS

#endif // AMDIS_PARALLEL_SOLVER_H

