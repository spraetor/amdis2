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
#include "parallel/ParallelSolver.hpp"
#include "parallel/ParallelDofMapping.hpp"
#include "parallel/MeshDistributor.hpp"

namespace AMDiS
{
  namespace Parallel
  {

    ParallelSolver::ParallelSolver(std::string name, bool globalIndices)
      : LinearSolverInterface(name),
        interiorMap(NULL),
        dofMap(FESPACE_WISE, globalIndices),
        parallelDofMappingsRegistered(false),
        meshDistributor(NULL),
        meshLevel(0)
    {
      setDofMapping(&dofMap);
    }


    ParallelSolver::~ParallelSolver()
    {
      if (parallelDofMappingsRegistered)
        meshDistributor->removeDofMap(dofMap);
    }

    void ParallelSolver::init(std::vector<const FiniteElemSpace*>& fe0,
                              std::vector<const FiniteElemSpace*>& fe1,
                              bool createGlobalMapping)
    {
      FUNCNAME("ParallelSolver::init()");

      TEST_EXIT(meshDistributor)("No mesh distributor object defined!\n");
      TEST_EXIT(fe0.size())("No component spaces provided!\n");
      TEST_EXIT(fe1.size())("No FE spaces provided!\n");

      componentSpaces = fe0;
      feSpaces = fe1;

      MeshLevelData& levelData = meshDistributor->getMeshLevelData();
      if (createGlobalMapping)
      {
        parallelDofMappingsRegistered = true;

        dofMap.init(componentSpaces, feSpaces);
        dofMap.setMpiComm(levelData.getMpiComm(meshLevel));

        // The meshes and fespaces of the dofmap come from the main problem.
        // The meshes and feSpaces of dofcomms may come from diff problems.
        dofMap.setDofComms(meshDistributor->getDofComms(), meshLevel);
        dofMap.clear();
        meshDistributor->registerDofMap(dofMap);
      }
    }

  }
} // end namespaces
