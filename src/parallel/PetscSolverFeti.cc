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


#include "MatrixVector.h"
#include "parallel/PetscHelper.h"
#include "parallel/PetscSolverFeti.h"
#include "parallel/PetscSolverFetiDebug.h"
#include "parallel/PetscSolverFetiMonitor.h"
#include "parallel/PetscSolverFetiStructs.h"
#include "parallel/PetscSolverFetiOperators.h"
#include "parallel/PetscSolverFetiTimings.h"
#include "parallel/StdMpi.h"
#include "parallel/MpiHelper.h"
#include "parallel/PetscSolverGlobalMatrix.h"
#include "io/VtkWriter.h"

namespace AMDiS
{
  namespace Parallel
  {

    using namespace std;

    PetscSolverFeti::PetscSolverFeti(string name)
      : PetscSolver(name),
        fetiSolverType(EXACT),
        dofMapSubDomain(FESPACE_WISE, true),
        primalDofMap(COMPONENT_WISE),
        dualDofMap(COMPONENT_WISE),
        interfaceDofMap(COMPONENT_WISE),
        localDofMap(COMPONENT_WISE),
        lagrangeMap(COMPONENT_WISE),
        interiorDofMap(COMPONENT_WISE),
        schurPrimalSolver(0),
        levelMode(1),
        subDomainIsLocal(true),
        subdomain(NULL),
        massMatrixSolver(NULL),
        printTimings(false),
        augmentedLagrange(false),
        nRankEdges(0),
        nOverallEdges(0),
        dirichletMode(0),
        stokesMode(false),
        pressureComponent(-1)
    {
      FUNCNAME("PetscSolverFeti::PetscSolverFeti()");

      string preconditionerName = "";
      Parameters::get(name + "->left precon", preconditionerName);
      if (preconditionerName == "" || preconditionerName == "no")
      {
        MSG("Create FETI-DP solver with no preconditioner!\n");
        fetiPreconditioner = FETI_NONE;
      }
      else if (preconditionerName == "dirichlet")
      {
        MSG("Create FETI-DP solver with Dirichlet preconditioner!\n");
        fetiPreconditioner = FETI_DIRICHLET;
      }
      else if (preconditionerName == "lumped")
      {
        MSG("Create FETI-DP solver with lumped preconditioner!\n");
        fetiPreconditioner = FETI_LUMPED;
      }
      else
      {
        ERROR_EXIT("Preconditioner \"%s\" not available!\n",
                   preconditionerName.c_str());
      }

      preconditionerName = "";
      Parameters::get(name + "->right precon", preconditionerName);
      if (preconditionerName != "" && preconditionerName != "no")
      {
        ERROR_EXIT("FETI-DP does not support right preconditioning! (parameter \"%s->right precon\" has value \"%s\")\n",
                   name.c_str(), preconditionerName.c_str());
      }

      Parameters::get(name + "->feti->schur primal solver", schurPrimalSolver);
      TEST_EXIT(schurPrimalSolver == 0 || schurPrimalSolver == 1)
      ("Wrong solver \"%d\"for the Schur primal complement!\n",
       schurPrimalSolver);

      Parameters::get(name + "->feti->stokes mode", stokesMode);
      if (stokesMode)
      {
        Parameters::get(name + "->feti->pressure component", pressureComponent);
        TEST_EXIT(pressureComponent >= 0)
        ("FETI-DP in Stokes mode, no pressure component defined!\n");
      }

      Parameters::get(name + "->feti->augmented lagrange", augmentedLagrange);

      Parameters::get(name + "->feti->symmetric", isSymmetric);

      {
        MSG("WARNING: CHECK THIS HERE BEFORE GOING INTO RUNNING MULTILEVEL FETI-DP!\n");
        Parameters::get("parallel->level mode", levelMode);
      }

      {
        int tmp = 0;
        Parameters::get(name + "->feti->inexact", tmp);
        if (tmp == 1)
          fetiSolverType = INEXACT;
        if (tmp == 2)
          fetiSolverType = INEXACT_REDUCED;
      }

      Parameters::get("parallel->print timings", printTimings);
    }


    void PetscSolverFeti::init(vector<const FiniteElemSpace*>& fe0,
                               vector<const FiniteElemSpace*>& fe1,
                               bool createGlobalMapping)
    {
      FUNCNAME_DBG("PetscSolverFeti::init()");

      super::init(fe0, fe1, createGlobalMapping);

      MeshLevelData& levelData = meshDistributor->getMeshLevelData();
      int nLevels = levelData.getNumberOfLevels();
      TEST_EXIT_DBG(nLevels >= 1)("nLevels < 1! Should not happen!\n");

      if (createGlobalMapping)
      {
        if (meshLevel + 1 < nLevels &&
            levelData.getMpiComm(meshLevel + 1) != MPI::COMM_SELF)
        {
          dofMapSubDomain.init(componentSpaces, feSpaces);
          dofMapSubDomain.setMpiComm(levelData.getMpiComm(meshLevel + 1));
          dofMapSubDomain.setDofComms(meshDistributor->getDofComms(), meshLevel + 1);
          dofMapSubDomain.clear();
          meshDistributor->registerDofMap(dofMapSubDomain);
        }
      }
    }


    void PetscSolverFeti::initialize()
    {
      FUNCNAME("PetscSolverFeti::initialize()");

#if (DEBUG != 0)
      MSG("Init FETI-DP on mesh level %d\n", meshLevel);
#endif

      TEST_EXIT_DBG(meshLevel + 2 <=
                    meshDistributor->getMeshLevelData().getNumberOfLevels())
      ("Mesh hierarchy does not contain at least %d levels!\n", meshLevel + 2);

      MeshLevelData& levelData = meshDistributor->getMeshLevelData();

      subDomainIsLocal = (levelData.getMpiComm(meshLevel + 1) == MPI::COMM_SELF);

      if (subdomain == NULL)
      {
        string subSolverInitStr = name + "->subsolver";
        string solverType = "petsc";
        Parameters::get(subSolverInitStr, solverType);
        solverType = "p_" + solverType;
        LinearSolverCreator* solverCreator =
          dynamic_cast<LinearSolverCreator*>(CreatorMap<LinearSolverInterface>::getCreator(solverType, name));
        TEST_EXIT(solverCreator)
        ("No valid solver type found in parameter \"%s\"\n",
         name.c_str());
        solverCreator->setName(subSolverInitStr);
        subdomain = dynamic_cast<PetscSolver*>(solverCreator->create());
        subdomain->setSymmetric(isSymmetric);
        subdomain->setHandleDirichletRows(dirichletMode == 0);
        subdomain->setMeshDistributor(meshDistributor, meshLevel + 1);
        subdomain->init(componentSpaces, feSpaces);

        delete solverCreator;
      }

      primalDofMap.init(componentSpaces, feSpaces);
      dualDofMap.init(componentSpaces, feSpaces, false);
      localDofMap.init(componentSpaces, feSpaces, !subDomainIsLocal);
      lagrangeMap.init(componentSpaces, feSpaces);

      if (stokesMode)
        interfaceDofMap.init(componentSpaces, feSpaces);

      if (fetiPreconditioner == FETI_DIRICHLET)
      {
        TEST_EXIT(levelMode == 1)
        ("Dirichlet preconditioner not yet implemented for multilevel FETI-DP\n");

        interiorDofMap.init(componentSpaces, feSpaces, false);
      }
    }


    void PetscSolverFeti::createDirichletData(Matrix<DOFMatrix*>& mat)
    {
      FUNCNAME("PetscSolverFeti::createDirichletData()");

      if (dirichletMode == 1)
      {
        int nComponents = mat.getNumRows();
        for (int component = 0; component < nComponents; component++)
        {
          DOFMatrix* dofMat = mat[component][component];
          if (!dofMat)
            continue;

          dirichletRows[component] = dofMat->getDirichletRows();
        }
      }
    }


    void PetscSolverFeti::createFetiData()
    {
      FUNCNAME("PetscSolverFeti::createFetiData()");

      double timeCounter = MPI::Wtime();

      MeshLevelData& levelData = meshDistributor->getMeshLevelData();

      primalDofMap.clear();
      dualDofMap.clear();
      lagrangeMap.clear();
      localDofMap.clear();
      if (fetiPreconditioner == FETI_DIRICHLET)
        interiorDofMap.clear();

      primalDofMap.setDofComms(meshDistributor->getDofComms(), meshLevel);
      lagrangeMap.setDofComms(meshDistributor->getDofComms(), meshLevel);

      primalDofMap.setMpiComm(levelData.getMpiComm(meshLevel));
      dualDofMap.setMpiComm(levelData.getMpiComm(meshLevel));
      lagrangeMap.setMpiComm(levelData.getMpiComm(meshLevel));
      localDofMap.setMpiComm(levelData.getMpiComm(meshLevel + 1));
      if (fetiPreconditioner == FETI_DIRICHLET)
        interiorDofMap.setMpiComm(levelData.getMpiComm(meshLevel + 1));

      localDofMap.setDofComms(meshDistributor->getDofComms(), meshLevel + 1);

      if (stokesMode)
      {
        interfaceDofMap.clear();
        interfaceDofMap.setDofComms(meshDistributor->getDofComms(), meshLevel);
        interfaceDofMap.setMpiComm(levelData.getMpiComm(0));
      }

      int nComponents = componentSpaces.size();
      for (int component = 0; component < nComponents; component++)
      {
        createPrimals(component);
        createDuals(component);
        createInterfaceNodes(component);
        createIndexB(component);
      }

      primalDofMap.update();
      dualDofMap.update();
      localDofMap.update();

      if (fetiPreconditioner == FETI_DIRICHLET)
        interiorDofMap.update();

      if (stokesMode)
        interfaceDofMap.update();

      for (int component = 0; component < nComponents; component++)
      {
        createLagrange(component);
        createAugmentedLagrange(component);
      }

      lagrangeMap.update();


      // === ===

      if (subDomainIsLocal)
      {
        MSG("WARNING: MAKE GENERAL!\n");

        rStartInterior = 0;
        int localDofs = localDofMap.getOverallDofs();

        mpi::getDofNumbering(domainComm, localDofs,
                             rStartInterior, nGlobalOverallInterior);
      }
      else
      {
        MSG("WARNING: MAKE GENERAL!\n");

        MeshLevelData& levelData = meshDistributor->getMeshLevelData();

        int groupRowsInterior = 0;
        if (levelData.getMpiComm(1).Get_rank() == 0)
          groupRowsInterior = localDofMap.getOverallDofs();

        mpi::getDofNumbering(domainComm, groupRowsInterior,
                             rStartInterior, nGlobalOverallInterior);

        int tmp = 0;
        if (levelData.getMpiComm(1).Get_rank() == 0)
          tmp = rStartInterior;

        levelData.getMpiComm(1).Allreduce(&tmp, &rStartInterior, 1,
                                          MPI_INT, MPI_SUM);
      }


      for (int i = 0; i < static_cast<int>(componentSpaces.size()); i++)
      {
        const FiniteElemSpace* feSpace = componentSpaces[i];

        MSG("FETI-DP data for %d-ith component (FE space %p):\n", i, feSpace);

        if (i == pressureComponent)
        {
          MSG("  nRankInterface = %d  nOverallInterface = %d\n",
              interfaceDofMap[i].nRankDofs,
              interfaceDofMap[i].nOverallDofs);
        }
        else
        {
          MSG("  nRankPrimals = %d   nLocalPrimals = %d  nOverallPrimals = %d\n",
              primalDofMap[i].nRankDofs,
              primalDofMap[i].nLocalDofs,
              primalDofMap[i].nOverallDofs);

          MSG("  nRankDuals = %d  nOverallDuals = %d\n",
              dualDofMap[i].nRankDofs,
              dualDofMap[i].nOverallDofs);

          MSG("  nRankLagrange = %d  nOverallLagrange = %d\n",
              lagrangeMap[i].nRankDofs,
              lagrangeMap[i].nOverallDofs);

          MSG("  nRankLocal = %d  nOverallLocal = %d\n",
              localDofMap[i].nRankDofs,
              localDofMap[i].nOverallDofs);
        }
      }

      subdomain->setDofMapping(&localDofMap);
      subdomain->setCoarseSpaceDofMapping(&primalDofMap);
      if (stokesMode)
        subdomain->setCoarseSpaceDofMapping(&interfaceDofMap, pressureComponent);

      if (printTimings)
      {
        MPI::COMM_WORLD.Barrier();
        timeCounter = MPI::Wtime() - timeCounter;
        MSG("FETI-DP timing 01: %.5f seconds (creation of basic data structures)\n",
            timeCounter);
      }


      bool writePrimals = false;
      Parameters::get("parallel->debug->write primals", writePrimals);
      if (writePrimals)
        PetscSolverFetiDebug::writePrimalFiles(*this);
    }


    void PetscSolverFeti::createPrimals(int component)
    {
      FUNCNAME("PetscSolverFeti::createPrimals()");

      if (component == pressureComponent)
        return;

      const FiniteElemSpace* feSpace = componentSpaces[component];

      // === Define all vertices on the interior boundaries of the macro mesh ===
      // === to be primal variables.                                          ===

      // Set of DOF indices that are considered to be primal variables.
      DofContainerSet& vertices =
        meshDistributor->getBoundaryDofInfo(feSpace, meshLevel).geoDofs[VERTEX];

      DofIndexSet primals;
      for (DofContainerSet::iterator it = vertices.begin();
           it != vertices.end(); ++it)
      {

        if (dirichletRows[component].count(**it))
          continue;

        if (meshLevel == 1 && not (*interiorMap)[component].isSet(**it))
          continue;

        if (subDomainIsLocal)
        {
          primals.insert(**it);
        }
        else
        {
          double e = 1e-8;
          WorldVector<double> c;
          feSpace->getMesh()->getDofIndexCoords(*it, feSpace, c);

          if ((fabs(c[0]) < e && fabs(c[1] - 12.5) < e) ||
              (fabs(c[0] - 25.0) < e && fabs(c[1] - 12.5) < e) ||
              (fabs(c[0] - 12.5) < e && fabs(c[1]) < e) ||
              (fabs(c[0] - 12.5) < e && fabs(c[1] - 25.0) < e) ||
              (fabs(c[0] - 12.5) < e && fabs(c[1] - 12.5) < e))
          {
            MSG("PRIMAL COORD %f %f\n", c[0], c[1]);
            primals.insert(**it);
          }
          else
          {
            MSG("OMMIT SOME PRIMAL!\n");
          }
        }
      }

      // === Calculate the number of primals that are owned by the rank and ===
      // === create local indices of the primals starting at zero.          ===

      for (DofIndexSet::iterator it = primals.begin(); it != primals.end(); ++it)
      {
        if (dofMap[feSpace].isRankDof(*it))
        {
          primalDofMap[component].insertRankDof(*it);
        }
        else
        {
          primalDofMap[component].insertNonRankDof(*it);
        }
      }
    }


    void PetscSolverFeti::createDuals(int component)
    {
      FUNCNAME("PetscSolverFeti::createDuals()");

      if (component == pressureComponent)
        return;

      const FiniteElemSpace* feSpace = componentSpaces[component];

      // === Create global index of the dual nodes on each rank. ===

      DofContainer allBoundaryDofs;
      meshDistributor->getAllBoundaryDofs(feSpace, meshLevel, allBoundaryDofs);

      for (DofContainer::iterator it = allBoundaryDofs.begin();
           it != allBoundaryDofs.end(); ++it)
      {
        if (dirichletRows[component].count(**it))
          continue;

        if (isPrimal(component, **it))
          continue;

        if (meshLevel == 1 && not (*interiorMap)[component].isSet(**it))
          continue;

        if (subDomainIsLocal || dofMapSubDomain[feSpace].isRankDof(**it))
          dualDofMap[component].insertRankDof(**it);
      }
    }


    void PetscSolverFeti::createInterfaceNodes(int component)
    {
      FUNCNAME("PetscSolverFeti::createInterfaceNodes()");

      if (component != pressureComponent)
        return;

      const FiniteElemSpace* feSpace = componentSpaces[component];
      DofContainer allBoundaryDofs;
      meshDistributor->getAllBoundaryDofs(feSpace, meshLevel, allBoundaryDofs);

      for (DofContainer::iterator it = allBoundaryDofs.begin();
           it != allBoundaryDofs.end(); ++it)
      {
        if (dirichletRows[component].count(**it))
          continue;

        if (dofMap[feSpace].isRankDof(**it))
          interfaceDofMap[component].insertRankDof(**it);
        else
          interfaceDofMap[component].insertNonRankDof(**it);
      }
    }


    void PetscSolverFeti::createLagrange(int component)
    {
      FUNCNAME("PetscSolverFeti::createLagrange()");

      if (component == pressureComponent)
        return;

      const FiniteElemSpace* feSpace = componentSpaces[component];
      Mesh* mesh = feSpace->getMesh();
      boundaryDofRanks[feSpace].clear();

      // Stores for all rank owned communication DOFs, if the counterpart is
      // a rank owned DOF in its subdomain. Thus, the following map stores to
      // each rank number all DOFs that fulfill this requirenment.
      map<int, std::set<DegreeOfFreedom>> subDomainRankDofs;

      if (not subDomainIsLocal)
      {
        StdMpi<vector<int>> stdMpi(domainComm);

        for (DofComm::Iterator it(meshDistributor->getDofComm(mesh, meshLevel).getRecvDofs(), feSpace);
             !it.end(); it.nextRank())
        {

          vector<int> dofs;
          dofs.reserve(it.getDofs().size());

          for (; !it.endDofIter(); it.nextDof())
          {
            if (dofMapSubDomain[feSpace].isRankDof(it.getDofIndex()))
              dofs.push_back(1);
            else
              dofs.push_back(0);
          }

          stdMpi.send(it.getRank(), dofs);
        }

        for (DofComm::Iterator it(meshDistributor->getDofComm(mesh, meshLevel).getSendDofs(), feSpace);
             !it.end(); it.nextRank())
          stdMpi.recv(it.getRank());

        stdMpi.startCommunication();

        for (DofComm::Iterator it(meshDistributor->getDofComm(mesh, meshLevel).getSendDofs(), feSpace);
             !it.end(); it.nextRank())
          for (; !it.endDofIter(); it.nextDof())
            if (!isPrimal(component, it.getDofIndex()) &&
                stdMpi.getRecvData(it.getRank())[it.getDofCounter()] == 1)
              subDomainRankDofs[it.getRank()].insert(it.getDofIndex());
      }

      if (dualDofMap[component].nLocalDofs == 0)
        return;


      // === Create for each dual node that is owned by the rank, the set ===
      // === of ranks that contain this node (denoted by W(x_j)).         ===

      int mpiRank = domainComm.Get_rank();
      for (DofComm::Iterator it(meshDistributor->getDofComm(mesh, meshLevel).getSendDofs(), feSpace);
           !it.end(); it.nextRank())
      {
        for (; !it.endDofIter(); it.nextDof())
        {
          if (!isPrimal(component, it.getDofIndex()))
          {
            boundaryDofRanks[feSpace][it.getDofIndex()].insert(mpiRank);

            // If the subdomain is local, always add the counterpart rank,
            // otherwise check if the other rank is the owner of the DOF in
            // its subdomain.
            if (subDomainIsLocal ||
                subDomainRankDofs[it.getRank()].count(it.getDofIndex()))
              boundaryDofRanks[feSpace][it.getDofIndex()].insert(it.getRank());
          }
        }
      }


      // === Communicate these sets for all rank owned dual nodes to other ===
      // === ranks that also have this node.                               ===

      StdMpi<vector<std::set<int>>> stdMpi(meshDistributor->getMpiComm(meshLevel));

      for (DofComm::Iterator it(meshDistributor->getDofComm(mesh, meshLevel).getSendDofs(), feSpace);
           !it.end(); it.nextRank())
        for (; !it.endDofIter(); it.nextDof())
          if (!isPrimal(component, it.getDofIndex()))
            if (subDomainIsLocal || subDomainRankDofs[it.getRank()].count(it.getDofIndex()))
              stdMpi.getSendData(it.getRank()).push_back(boundaryDofRanks[feSpace][it.getDofIndex()]);

      stdMpi.updateSendDataSize();

      for (DofComm::Iterator it(meshDistributor->getDofComm(mesh, meshLevel).getRecvDofs(), feSpace);
           !it.end(); it.nextRank())
      {
        bool recvFromRank = false;
        for (; !it.endDofIter(); it.nextDof())
        {
          if (!isPrimal(component, it.getDofIndex()))
          {
            if (subDomainIsLocal || dofMapSubDomain[feSpace].isRankDof(it.getDofIndex()))
            {
              recvFromRank = true;
              break;
            }
          }
        }

        if (recvFromRank)
          stdMpi.recv(it.getRank());
      }

      stdMpi.startCommunication();

      for (DofComm::Iterator it(meshDistributor->getDofComm(mesh, meshLevel).getRecvDofs(), feSpace);
           !it.end(); it.nextRank())
      {
        int i = 0;
        for (; !it.endDofIter(); it.nextDof())
        {
          if (!isPrimal(component, it.getDofIndex()))
          {
            if (subDomainIsLocal || dofMapSubDomain[feSpace].isRankDof(it.getDofIndex()))
            {
              boundaryDofRanks[feSpace][it.getDofIndex()] =
                stdMpi.getRecvData(it.getRank())[i++];
            }
            else
            {
              lagrangeMap[component].insertNonRankDof(it.getDofIndex());
            }
          }
        }
      }


      // === Reserve for each dual node, on the rank that owns this node, the ===
      // === appropriate number of Lagrange constraints.                      ===

      int nRankLagrange = 0;
      DofMap& dualMap = dualDofMap[component].getMap();
      for (DofMap::iterator it = dualMap.begin(); it != dualMap.end(); ++it)
      {
        if (dofMap[feSpace].isRankDof(it->first))
        {
          lagrangeMap[component].insertRankDof(it->first, nRankLagrange);
          int degree = boundaryDofRanks[feSpace][it->first].size();
          nRankLagrange += (degree * (degree - 1)) / 2;
        }
        else
        {
          lagrangeMap[component].insertNonRankDof(it->first);
        }
      }
      lagrangeMap[component].nRankDofs = nRankLagrange;
    }


    void PetscSolverFeti::createAugmentedLagrange(int component)
    {
      FUNCNAME("PetscSolverFeti::createAugmentedLagrange()");

      if (!augmentedLagrange)
        return;
    }


    void PetscSolverFeti::createIndexB(int component)
    {
      FUNCNAME("PetscSolverFeti::createIndexB()");

      const FiniteElemSpace* feSpace = componentSpaces[component];
      DOFAdmin* admin = feSpace->getAdmin();

      // === To ensure that all interior node on each rank are listen first in ===
      // === the global index of all B nodes, insert all interior nodes first, ===
      // === without defining a correct index.                                 ===

      int nLocalInterior = 0;
      for (int i = 0; i < admin->getUsedSize(); i++)
      {
        if (admin->isDofFree(i) ||
            isPrimal(component, i) ||
            isDual(component, i) ||
            isInterface(component, i) ||
            dirichletRows[component].count(i))
          continue;

        if (meshLevel == 1 &&  not (*interiorMap)[component].isSet(i))
          continue;

        if (subDomainIsLocal)
        {
          localDofMap[component].insertRankDof(i, nLocalInterior);

          if (fetiPreconditioner == FETI_DIRICHLET)
            interiorDofMap[component].insertRankDof(i, nLocalInterior);

          nLocalInterior++;
        }
        else
        {
          if (dofMapSubDomain[feSpace].isRankDof(i))
            localDofMap[component].insertRankDof(i);
          else
            localDofMap[component].insertNonRankDof(i);

          TEST_EXIT_DBG(fetiPreconditioner == FETI_NONE)
          ("Not yet implemnted!\n");
        }
      }

      // === And finally, add the global indicies of all dual nodes. ===

      for (DofMap::iterator it = dualDofMap[component].getMap().begin();
           it != dualDofMap[component].getMap().end(); ++it)
      {
        if (subDomainIsLocal)
        {
          localDofMap[component].insertRankDof(it->first);
        }
        else
        {
          if (dofMapSubDomain[feSpace].isRankDof(it->first))
            localDofMap[component].insertRankDof(it->first);
          else
            localDofMap[component].insertNonRankDof(it->first);
        }
      }
    }


    void PetscSolverFeti::createMatLagrange()
    {
      FUNCNAME("PetscSolverFeti::createMatLagrange()");

      double wtime = MPI::Wtime();
      int mpiRank = domainComm.Get_rank();

      // === Create distributed matrix for Lagrange constraints. ===

      MatCreateAIJ(domainComm,
                   lagrangeMap.getRankDofs(), localDofMap.getRankDofs(),
                   lagrangeMap.getOverallDofs(), nGlobalOverallInterior,
                   2, PETSC_NULL, 2, PETSC_NULL,
                   &mat_lagrange);
      MatSetOption(mat_lagrange, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);


      Vec vec_scale_lagrange;
      createVec(lagrangeMap, vec_scale_lagrange);


      // === Create for all duals the corresponding Lagrange constraints. On ===
      // === each rank we traverse all pairs (n, m) of ranks, with n < m,    ===
      // === that contain this node. If the current rank number is r, and    ===
      // === n == r, the rank sets 1.0 for the corresponding constraint, if  ===
      // === m == r, than the rank sets -1.0 for the corresponding           ===
      // === constraint.                                                     ===

      for (unsigned int component = 0; component < componentSpaces.size(); component++)
      {
        DofMap& dualMap = dualDofMap[component].getMap();

        for (DofMap::iterator it = dualMap.begin(); it != dualMap.end(); ++it)
        {
          TEST_EXIT_DBG(boundaryDofRanks[componentSpaces[component]].count(it->first))
          ("Should not happen!\n");

          // Global index of the first Lagrange constriant for this node.
          int index = lagrangeMap.getMatIndex(component, it->first);

          // Copy set of all ranks that contain this dual node.
          vector<int> W(boundaryDofRanks[componentSpaces[component]][it->first].begin(),
                        boundaryDofRanks[componentSpaces[component]][it->first].end());
          // Number of ranks that contain this dual node.
          int degree = W.size();

          TEST_EXIT_DBG(degree > 1)("Should not happen!\n");

          int counter = 0;
          for (int i = 0; i < degree; i++)
          {
            for (int j = i + 1; j < degree; j++)
            {
              if (W[i] == mpiRank || W[j] == mpiRank)
              {
                MatSetValue(mat_lagrange,
                            index + counter,
                            localDofMap.getMatIndex(component, it->first) + rStartInterior,
                            (W[i] == mpiRank ? 1.0 : -1.0),
                            INSERT_VALUES);
              }
              counter++;
            }
          }

          // === Create scaling factors for scaling the lagrange matrix, which ===
          // === is required for FETI-DP preconditioners.                      ===

          if (dofMap[componentSpaces[component]].isRankDof(it->first))
          {
            int nConstraints = (degree * (degree - 1)) / 2;
            for (int i = 0; i < nConstraints; i++)
            {
              VecSetValue(vec_scale_lagrange,
                          index + i,
                          1.0 / static_cast<double>(degree),
                          INSERT_VALUES);
            }
          }
        }
      }

      MatAssemblyBegin(mat_lagrange, MAT_FINAL_ASSEMBLY);
      MatAssemblyEnd(mat_lagrange, MAT_FINAL_ASSEMBLY);


#if (DEBUG != 0)
      {
        int nZeroRows = PetscSolverFetiDebug::testZeroRows(mat_lagrange);
        int m,n;
        MatGetSize(mat_lagrange, &m ,&n);
        mpi::globalAdd(domainComm, nZeroRows);
        MSG("Lagrange matrix has %d zero rows and global size of %d %d!\n", nZeroRows, m, n);

        TEST_EXIT(nZeroRows == 0)("Lagrange matrix has zero rows!\n");
      }
#endif


      // === If required, create \ref mat_lagrange_scaled ===

      VecAssemblyBegin(vec_scale_lagrange);
      VecAssemblyEnd(vec_scale_lagrange);

      if (fetiPreconditioner != FETI_NONE ||
          fetiSolverType == INEXACT ||
          stokesMode)
      {
        MatDuplicate(mat_lagrange, MAT_COPY_VALUES, &mat_lagrange_scaled);
        MatDiagonalScale(mat_lagrange_scaled, vec_scale_lagrange, PETSC_NULL);
      }

      VecDestroy(&vec_scale_lagrange);


      // === Print final timings. ===

      if (printTimings)
      {
        MPI::COMM_WORLD.Barrier();
        MSG("FETI-DP timing 05: %.5f seconds (creation of lagrange constraint matrix)\n",
            MPI::Wtime() - wtime);
      }
    }


    bool PetscSolverFeti::testWirebasketEdge(BoundaryObject& edge, const FiniteElemSpace* feSpace)
    {
      FUNCNAME("PetscSolverFeti::testWirebasketEdge()");

      return true;

      if (meshDistributor->getMesh()->getDim() == 2)
        return true;

      if (meshDistributor->getIntBoundary(meshLevel).getDegreeOwn(edge) != 3)
        return false;

      return true;

      Element* el = edge.el;
      int i0 = el->getVertexOfEdge(edge.ithObj, 0);
      int i1 = el->getVertexOfEdge(edge.ithObj, 1);
      DegreeOfFreedom d0 = el->getDof(i0, 0);
      DegreeOfFreedom d1 = el->getDof(i1, 0);
      WorldVector<double> c0, c1;
      el->getMesh()->getDofIndexCoords(d0, feSpace, c0);
      el->getMesh()->getDofIndexCoords(d1, feSpace, c1);
      bool xe = fabs(c0[0] - c1[0]) < 1e-8;
      bool ye = fabs(c0[1] - c1[1]) < 1e-8;
      bool ze = fabs(c0[2] - c1[2]) < 1e-8;
      int counter = static_cast<int>(xe) + static_cast<int>(ye) + static_cast<int>(ze);
      return (counter == 2);
    }


    vector<vector<BoundaryObject>> PetscSolverFeti::getCoarseEdges()
    {
      FUNCNAME("PetscSolverFeti::getAugmentedLagrange()");

      InteriorBoundary& intBound = meshDistributor->getIntBoundary(meshLevel);
      std::set<BoundaryObject> allEdges;
      for (InteriorBoundary::iterator it(intBound.getOwn()); !it.end(); ++it)
        if (it->rankObj.subObj == EDGE &&
            testWirebasketEdge(it->rankObj, feSpaces[0]) &&
            allEdges.count(it->rankObj) == 0)
          allEdges.insert(it->rankObj);

      int nEmptyEdges = 0;
      vector<vector<BoundaryObject>> data;
      for (std::set<BoundaryObject>::iterator it = allEdges.begin();
           it != allEdges.end(); ++it)
      {
        DofContainer edgeDofs;
        it->el->getAllDofs(feSpaces[0], *it, edgeDofs);
        if (edgeDofs.size() == 0)
        {
          nEmptyEdges++;
        }
        else
        {
          vector<BoundaryObject> oneBoundary;
          oneBoundary.push_back(*it);
          data.push_back(oneBoundary);
        }
      }

      int nEdges = allEdges.size();
      mpi::globalAdd(nEdges);
      mpi::globalAdd(nEmptyEdges);

      MSG("Coarse space edges: %d (empty: %d)\n", nEdges, nEmptyEdges);

      return data;
    }


    vector<vector<BoundaryObject>> PetscSolverFeti::getCoarseFaces()
    {
      FUNCNAME("PetscSolverFeti::getAugmentedLagrange()");

      InteriorBoundary& intBound = meshDistributor->getIntBoundary(meshLevel);
      map<int, std::set<BoundaryObject>> allFaces;
      for (InteriorBoundary::iterator it(intBound.getOwn()); !it.end(); ++it)
        if (it->rankObj.subObj == FACE)
          allFaces[it.getRank()].insert(it->rankObj);

#if 0
      std::set<BoundaryObject> allMyEdges;
      TraverseStack stack;
      ElInfo* elInfo = stack.traverseFirst(meshDistributor->getMesh(), 0, Mesh::CALL_EL_LEVEL | Mesh::FILL_BOUND);
      while (elInfo)
      {
        Element* el = elInfo->getElement();
        for (int i = 0; i < el->getGeo(EDGE); i++)
        {
          BoundaryObject bobj(el, elInfo->getType(), EDGE, i);
          if (intBound.getDegreeOwn(bobj) == 1 && elInfo->getBoundary(EDGE, i) == INTERIOR)
          {
            allMyEdges.insert(bobj);
          }
        }
        elInfo = stack.traverseNext(elInfo);
      }



      for (map<int, std::set<BoundaryObject>>::iterator it = allFaces.begin();
           it != allFaces.end(); ++it)
      {
        if (it->second.size() == 2)
        {
          vector<AtomicBoundary>& bs = intBound.getOwn()[it->first];
          for (int i = 0; i < static_cast<int>(bs.size()); i++)
          {
            if (bs[i].rankObj.subObj == EDGE &&
                intBound.getDegreeOwn(bs[i].rankObj) == 1 &&
                allMyEdges.count(bs[i].rankObj))
            {
              MSG("FOUND AN EDGE: %d %d %d\n",
                  bs[i].rankObj.elIndex, bs[i].rankObj.subObj, bs[i].rankObj.ithObj);
              it->second.insert(bs[i].rankObj);
            }
          }
        }
      }
#endif


      int nEmptyFaces = 0;
      vector<vector<BoundaryObject>> data;
      for (map<int, std::set<BoundaryObject>>::iterator it = allFaces.begin();
           it != allFaces.end(); ++it)
      {
        vector<BoundaryObject> oneBoundary;
        for (std::set<BoundaryObject>::iterator bIt = it->second.begin();
             bIt != it->second.end(); ++bIt)
        {
          DofContainer faceDofs;
          bIt->el->getAllDofs(feSpaces[0], *bIt, faceDofs);
          if (faceDofs.size())
            oneBoundary.push_back(*bIt);
        }
        if (oneBoundary.size())
          data.push_back(oneBoundary);
        else
          nEmptyFaces++;
      }

      int nFaces = allFaces.size();
      mpi::globalAdd(nFaces);
      mpi::globalAdd(nEmptyFaces);

      MSG("Coarse space faces: %d (empty: %d)\n", nFaces, nEmptyFaces);

      return data;
    }


    void PetscSolverFeti::createMatAugmentedLagrange()
    {
      FUNCNAME("PetscSolverFeti::createMatAugmentedLagrange()");

      if (!augmentedLagrange)
        return;

      double wtime = MPI::Wtime();

      vector<vector<BoundaryObject>> allEdges = getCoarseEdges();
      vector<vector<BoundaryObject>> allFaces = getCoarseFaces();
      allEdges.insert(allEdges.end(), allFaces.begin(), allFaces.end());

      nRankEdges = allEdges.size();
      int rStartEdges = 0;
      mpi::getDofNumbering(domainComm, nRankEdges, rStartEdges, nOverallEdges);

      MSG("nRankEdges = %d, nOverallEdges = %d\n", nRankEdges, nOverallEdges);

      nRankEdges *= componentSpaces.size();
      rStartEdges *= componentSpaces.size();
      nOverallEdges *= componentSpaces.size();

      MatCreateAIJ(domainComm,
                   nRankEdges, lagrangeMap.getRankDofs(),
                   nOverallEdges, lagrangeMap.getOverallDofs(),
                   2, PETSC_NULL, 2, PETSC_NULL,
                   &mat_augmented_lagrange);
      MatSetOption(mat_augmented_lagrange, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

      int rowCounter = rStartEdges;
      for (vector<vector<BoundaryObject>>::iterator it = allEdges.begin();
           it != allEdges.end(); ++it)
      {
        for (int component = 0; component < static_cast<int>(componentSpaces.size()); component++)
        {
          for (vector<BoundaryObject>::iterator edgeIt = it->begin();
               edgeIt != it->end(); ++edgeIt)
          {

            DofContainer edgeDofs;
            edgeIt->el->getAllDofs(componentSpaces[component], *edgeIt, edgeDofs);

            TEST_EXIT(edgeDofs.size())("Should not happen!\n");

            for (DofContainer::iterator it = edgeDofs.begin();
                 it != edgeDofs.end(); ++it)
            {
              TEST_EXIT(isPrimal(component, **it) == false)
              ("Should not be primal!\n");

              int col = lagrangeMap.getMatIndex(component, **it);
              double value = 1.0;
              MatSetValue(mat_augmented_lagrange, rowCounter, col, value, INSERT_VALUES);
            }
          }
          rowCounter++;
        }
      }

      MatAssemblyBegin(mat_augmented_lagrange, MAT_FINAL_ASSEMBLY);
      MatAssemblyEnd(mat_augmented_lagrange, MAT_FINAL_ASSEMBLY);

      int nZeroRows = PetscSolverFetiDebug::testZeroRows(mat_augmented_lagrange);
      int m,n;
      MatGetSize(mat_augmented_lagrange, &m ,&n);
      MSG("Augmented lagrange matrix has %d zero rows and global size of %d %d!\n", nZeroRows, m, n);

      if (printTimings)
      {
        MPI::COMM_WORLD.Barrier();
        MSG("FETI-DP timing 05a: %.5f seconds (creation of augmented lagrange constraint matrix)\n",
            MPI::Wtime() - wtime);
      }
    }


    void PetscSolverFeti::createSchurPrimalKsp()
    {
      FUNCNAME("PetscSolverFeti::createSchurPrimalKsp()");

      if (schurPrimalSolver == 0)
      {
        MSG("Create iterative schur primal solver on level %d!\n", meshLevel);

        if (augmentedLagrange == false)
        {
          schurPrimalData.subSolver = subdomain;

          VecCreateMPI(meshDistributor->getMeshLevelData().getMpiComm(meshLevel),
                       localDofMap.getRankDofs(),
                       nGlobalOverallInterior, &(schurPrimalData.tmp_vec_b));

          createVec(primalDofMap, schurPrimalData.tmp_vec_primal);


          MatCreateShell(domainComm,
                         primalDofMap.getRankDofs(),
                         primalDofMap.getRankDofs(),
                         primalDofMap.getOverallDofs(),
                         primalDofMap.getOverallDofs(),
                         &schurPrimalData,
                         &mat_schur_primal);
          MatShellSetOperation(mat_schur_primal, MATOP_MULT,
                               (void(*)(void))petscMultMatSchurPrimal);
        }
        else
        {
          schurPrimalAugmentedData.subSolver = subdomain;
          schurPrimalAugmentedData.nestedVec = true;

          createVec(localDofMap, schurPrimalAugmentedData.tmp_vec_b0, nGlobalOverallInterior);
          createVec(localDofMap, schurPrimalAugmentedData.tmp_vec_b1, nGlobalOverallInterior);
          createVec(primalDofMap, schurPrimalAugmentedData.tmp_vec_primal);
          createVec(lagrangeMap, schurPrimalAugmentedData.tmp_vec_lagrange);

          schurPrimalAugmentedData.mat_lagrange = &mat_lagrange;
          schurPrimalAugmentedData.mat_augmented_lagrange = &mat_augmented_lagrange;

          MatCreateShell(domainComm,
                         primalDofMap.getRankDofs() + nRankEdges,
                         primalDofMap.getRankDofs() + nRankEdges,
                         primalDofMap.getOverallDofs() + nOverallEdges,
                         primalDofMap.getOverallDofs() + nOverallEdges,
                         &schurPrimalAugmentedData,
                         &mat_schur_primal);
          MatShellSetOperation(mat_schur_primal, MATOP_MULT,
                               (void(*)(void))petscMultMatSchurPrimalAugmented);
        }

        KSPCreate(domainComm, &ksp_schur_primal);
#if (PETSC_VERSION_MINOR >= 5)
        KSPSetOperators(ksp_schur_primal, mat_schur_primal, mat_schur_primal);
#else
        KSPSetOperators(ksp_schur_primal, mat_schur_primal, mat_schur_primal, SAME_NONZERO_PATTERN);
#endif
        KSPSetOptionsPrefix(ksp_schur_primal, "schur_primal_");
        KSPSetType(ksp_schur_primal, KSPGMRES);
        KSPSetFromOptions(ksp_schur_primal);
      }
      else
      {
        MSG("Create direct schur primal solver!\n");

        double wtime = MPI::Wtime();

        // === Create explicit matrix representation of the Schur primal system. ===

        if (!augmentedLagrange)
          createMatExplicitSchurPrimal();
        else
          createMatExplicitAugmentedSchurPrimal();


        // === Create KSP solver object and set appropriate solver options. ===

        KSPCreate(domainComm, &ksp_schur_primal);
#if (PETSC_VERSION_MINOR >= 5)
        KSPSetOperators(ksp_schur_primal, mat_schur_primal, mat_schur_primal);
#else
        KSPSetOperators(ksp_schur_primal, mat_schur_primal, mat_schur_primal, SAME_NONZERO_PATTERN);
#endif
        KSPSetOptionsPrefix(ksp_schur_primal, "schur_primal_");
        KSPSetType(ksp_schur_primal, KSPPREONLY);
        PC pc_schur_primal;
        KSPGetPC(ksp_schur_primal, &pc_schur_primal);
        PCSetType(pc_schur_primal, PCLU);
        PCFactorSetMatSolverPackage(pc_schur_primal, MATSOLVERMUMPS);
        KSPSetFromOptions(ksp_schur_primal);


        // === And finally print timings, if required. ===

        if (printTimings)
        {
          MPI::COMM_WORLD.Barrier();
          MatInfo minfo;
          MatGetInfo(mat_schur_primal, MAT_GLOBAL_SUM, &minfo);
          MSG("Schur primal matrix nnz = %f\n", minfo.nz_used);

          MSG("FETI-DP timing 06: %.5f seconds (creation of schur primal matrix)\n",
              MPI::Wtime() - wtime);

          wtime = MPI::Wtime();
          KSPSetUp(ksp_schur_primal);
          KSPSetUpOnBlocks(ksp_schur_primal);
          MPI::COMM_WORLD.Barrier();
          MSG("FETI-DP timing 07: %.5f seconds (factorization of primal schur matrix).\n",
              MPI::Wtime() - wtime);
        }
      }
    }


    void PetscSolverFeti::createMatExplicitSchurPrimal()
    {
      FUNCNAME("PetscSolverFeti::createMatExplicitSchurPrimal()");

      int creationMode = 0;
      Parameters::get(name + "->feti->schur primal creation mode", creationMode);
      if (creationMode == 0)
      {
        // matK = inv(A_BB) A_BPi
        Mat matK;
        petsc_helper::blockMatMatSolve(domainComm,
                                       subdomain->getSolver(),
                                       subdomain->getMatInteriorCoarse(),
                                       matK);

        // mat_schur_primal = A_PiPi - A_PiB inv(A_BB) A_BPi
        //                  = A_PiPi - A_PiB matK
        MatMatMult(subdomain->getMatCoarseInterior(), matK, MAT_INITIAL_MATRIX,
                   PETSC_DEFAULT, &mat_schur_primal);
        MatAYPX(mat_schur_primal, -1.0, subdomain->getMatCoarse(), DIFFERENT_NONZERO_PATTERN);

        MatDestroy(&matK);
      }
      else
      {
        schurPrimalData.subSolver = subdomain;

        createVec(localDofMap, schurPrimalData.tmp_vec_b, nGlobalOverallInterior);
        createVec(primalDofMap, schurPrimalData.tmp_vec_primal);
        Mat tmp;
        MatCreateShell(domainComm,
                       primalDofMap.getRankDofs(),
                       primalDofMap.getRankDofs(),
                       primalDofMap.getOverallDofs(),
                       primalDofMap.getOverallDofs(),
                       &schurPrimalData,
                       &tmp);
        MatShellSetOperation(tmp, MATOP_MULT,
                             (void(*)(void))petscMultMatSchurPrimal);
        MatComputeExplicitOperator(tmp, &mat_schur_primal);

        MatDestroy(&tmp);
        schurPrimalData.subSolver = NULL;
        VecDestroy(&schurPrimalData.tmp_vec_b);
        VecDestroy(&schurPrimalData.tmp_vec_primal);
      }
    }


    void PetscSolverFeti::createMatExplicitAugmentedSchurPrimal()
    {
      FUNCNAME("PetscSolverFeti::createMatExplicitAugmentedSchurPrimal()");

      int creationMode = 0;
      Parameters::get("parallel->feti->schur primal creation mode", creationMode);
      if (creationMode == 0)
      {
        // qj = -Q J
        Mat qj;
        MatMatMult(mat_augmented_lagrange, mat_lagrange, MAT_INITIAL_MATRIX,
                   PETSC_DEFAULT, &qj);
        MatScale(qj, -1.0);


        // matTmp = inv(A_BB) A_BPi
        Mat matTmp;
        petsc_helper::blockMatMatSolve(domainComm,
                                       subdomain->getSolver(),
                                       subdomain->getMatInteriorCoarse(),
                                       matTmp);


        // mat00 = A_PiPi - A_PiB inv(A_BB) A_BPi
        //       = A_PiPi - A_PiB matTmp
        Mat mat00;
        MatMatMult(subdomain->getMatCoarseInterior(), matTmp, MAT_INITIAL_MATRIX,
                   PETSC_DEFAULT, &mat00);
        MatAYPX(mat00, -1.0, subdomain->getMatCoarse(), DIFFERENT_NONZERO_PATTERN);


        // mat10 = -Q J inv(A_BB) A_BPi
        //       = qj matTmp
        Mat mat10;
        MatMatMult(qj, matTmp, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &mat10);


        // matTmp = inv(A_BB) trans(J) trans(Q)
        Mat qT, jTqT;
        MatTranspose(mat_augmented_lagrange, MAT_INITIAL_MATRIX, &qT);
        MatTransposeMatMult(mat_lagrange, qT, MAT_INITIAL_MATRIX, PETSC_DEFAULT,
                            &jTqT);
        petsc_helper::blockMatMatSolve(domainComm,
                                       subdomain->getSolver(),
                                       jTqT,
                                       matTmp);
        MatDestroy(&qT);
        MatDestroy(&jTqT);

        // mat01 = -A_PiB inv(A_BB) trans(J) trans(Q)
        //       = -A_PiB matTmp
        Mat mat01;
        MatMatMult(subdomain->getMatCoarseInterior(), matTmp, MAT_INITIAL_MATRIX,
                   PETSC_DEFAULT, &mat01);
        MatScale(mat01, -1.0);

        // mat11 = -Q J inv(A_BB) trans(J) trans(Q)
        //       = qj matTmp
        Mat mat11;
        MatMatMult(qj, matTmp, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &mat11);

        MatDestroy(&matTmp);
        MatDestroy(&qj);

        Mat nestMat[4] = {mat00, mat01, mat10, mat11};
        MatCreateNest(PETSC_COMM_WORLD, 2, PETSC_NULL, 2, PETSC_NULL, nestMat, &matTmp);
        petsc_helper::matNestConvert(matTmp, mat_schur_primal);

        MatDestroy(&mat00);
        MatDestroy(&mat01);
        MatDestroy(&mat10);
        MatDestroy(&mat11);
        MatDestroy(&matTmp);
      }
      else
      {
        Mat tmp;

        schurPrimalAugmentedData.subSolver = subdomain;
        schurPrimalAugmentedData.nestedVec = false;

        createVec(localDofMap, schurPrimalAugmentedData.tmp_vec_b0, nGlobalOverallInterior);
        createVec(localDofMap, schurPrimalAugmentedData.tmp_vec_b1, nGlobalOverallInterior);
        createVec(primalDofMap, schurPrimalAugmentedData.tmp_vec_primal);
        createVec(lagrangeMap, schurPrimalAugmentedData.tmp_vec_lagrange);

        schurPrimalAugmentedData.mat_lagrange = &mat_lagrange;
        schurPrimalAugmentedData.mat_augmented_lagrange = &mat_augmented_lagrange;

        MatCreateShell(domainComm,
                       primalDofMap.getRankDofs() + nRankEdges,
                       primalDofMap.getRankDofs() + nRankEdges,
                       primalDofMap.getOverallDofs() + nOverallEdges,
                       primalDofMap.getOverallDofs() + nOverallEdges,
                       &schurPrimalAugmentedData,
                       &tmp);
        MatShellSetOperation(tmp, MATOP_MULT,
                             (void(*)(void))petscMultMatSchurPrimalAugmented);

        MatComputeExplicitOperator(tmp, &mat_schur_primal);

        MatDestroy(&tmp);
        schurPrimalAugmentedData.subSolver = NULL;
        schurPrimalAugmentedData.mat_lagrange = NULL;
        schurPrimalAugmentedData.mat_augmented_lagrange = NULL;
        VecDestroy(&schurPrimalAugmentedData.tmp_vec_b0);
        VecDestroy(&schurPrimalAugmentedData.tmp_vec_b1);
        VecDestroy(&schurPrimalAugmentedData.tmp_vec_primal);
        VecDestroy(&schurPrimalAugmentedData.tmp_vec_lagrange);
      }
    }


    void PetscSolverFeti::destroySchurPrimalKsp()
    {
      FUNCNAME("PetscSolverFeti::destroySchurPrimal()");

      if (schurPrimalSolver == 0)
      {
        if (augmentedLagrange == false)
        {
          schurPrimalData.subSolver = NULL;

          VecDestroy(&schurPrimalData.tmp_vec_b);
          VecDestroy(&schurPrimalData.tmp_vec_primal);
        }
        else
        {
          schurPrimalAugmentedData.subSolver = NULL;
          schurPrimalAugmentedData.mat_lagrange = NULL;
          schurPrimalAugmentedData.mat_augmented_lagrange = NULL;

          VecDestroy(&schurPrimalAugmentedData.tmp_vec_b0);
          VecDestroy(&schurPrimalAugmentedData.tmp_vec_b1);
          VecDestroy(&schurPrimalAugmentedData.tmp_vec_primal);
          VecDestroy(&schurPrimalAugmentedData.tmp_vec_lagrange);
        }
      }

      MatDestroy(&mat_schur_primal);
      KSPDestroy(&ksp_schur_primal);
    }


    void PetscSolverFeti::createFetiKsp()
    {
      FUNCNAME("PetscSolverFeti::createFetiKsp()");

      switch (fetiSolverType)
      {
      case EXACT:
        createFetiExactKsp();
        break;
      case INEXACT:
        createFetiInexactKsp();
        break;
      case INEXACT_REDUCED:
        createFetiInexactReducedKsp();
        break;
      default:
        ERROR_EXIT("Should not happen!\n");
      }
    }


    void PetscSolverFeti::createFetiExactKsp()
    {
      FUNCNAME("PetscSolverFeti::createFetiExactKsp()");

      // === Create FETI-DP solver object. ===

      fetiData.mat_lagrange = &mat_lagrange;
      fetiData.subSolver = subdomain;
      fetiData.ksp_schur_primal = &ksp_schur_primal;

      VecCreateMPI(meshDistributor->getMeshLevelData().getMpiComm(meshLevel),
                   localDofMap.getRankDofs(),
                   nGlobalOverallInterior, &(fetiData.tmp_vec_b0));

      createVec(lagrangeMap, fetiData.tmp_vec_lagrange);
      createVec(primalDofMap, fetiData.tmp_vec_primal0);

      if (stokesMode == false)
      {
        MatCreateShell(domainComm,
                       lagrangeMap.getRankDofs(),
                       lagrangeMap.getRankDofs(),
                       lagrangeMap.getOverallDofs(),
                       lagrangeMap.getOverallDofs(),
                       &fetiData, &mat_feti);
        if (augmentedLagrange == false)
        {
          MatShellSetOperation(mat_feti, MATOP_MULT,
                               (void(*)(void))petscMultMatFeti);
        }
        else
        {
          fetiData.mat_augmented_lagrange = &mat_augmented_lagrange;
          createVec(primalDofMap, fetiData.tmp_vec_primal1);
          MatShellSetOperation(mat_feti, MATOP_MULT,
                               (void(*)(void))petscMultMatFetiAugmented);
        }
      }
      else
      {
        TEST_EXIT_DBG(!augmentedLagrange)("Not yet implemented!\n");

        createVec(localDofMap, fetiData.tmp_vec_b1, nGlobalOverallInterior);
        createVec(primalDofMap, fetiData.tmp_vec_primal1);
        createVec(interfaceDofMap, fetiData.tmp_vec_interface);

        MatCreateShell(domainComm,
                       interfaceDofMap.getRankDofs() + lagrangeMap.getRankDofs(),
                       interfaceDofMap.getRankDofs() + lagrangeMap.getRankDofs(),
                       interfaceDofMap.getOverallDofs() + lagrangeMap.getOverallDofs(),
                       interfaceDofMap.getOverallDofs() + lagrangeMap.getOverallDofs(),
                       &fetiData, &mat_feti);
        MatShellSetOperation(mat_feti, MATOP_MULT,
                             (void(*)(void))petscMultMatFetiInterface);
      }

      KSPCreate(domainComm, &ksp_feti);
#if (PETSC_VERSION_MINOR >= 5)
      KSPSetOperators(ksp_feti, mat_feti, mat_feti);
#else
      KSPSetOperators(ksp_feti, mat_feti, mat_feti, SAME_NONZERO_PATTERN);
#endif
      KSPSetOptionsPrefix(ksp_feti, "feti_");
      KSPSetType(ksp_feti, KSPGMRES);
      KSPSetTolerances(ksp_feti, 0, 1e-8, 1e+3, 1000);
      KSPSetFromOptions(ksp_feti);


      // === Set KSP monitor. ===

      bool monitor = false;
      Parameters::get(name + "->feti->monitor", monitor);
      if (monitor)
      {
        if (stokesMode)
          KSPMonitorSet(ksp_feti, KSPMonitorFetiStokes, &fetiKspData, PETSC_NULL);
        else
          KSPMonitorSet(ksp_feti, KSPMonitorTrueResidualNorm, PETSC_NULL, PETSC_NULL);
      }


      // === Create null space objects. ===

      createNullSpace();


      switch (fetiPreconditioner)
      {
      case FETI_DIRICHLET:
        KSPGetPC(ksp_feti, &precon_feti);
        createFetiPreconDirichlet(precon_feti);
        break;

      case FETI_LUMPED:
        KSPGetPC(ksp_feti, &precon_feti);
        createFetiPreconLumped(precon_feti);
        break;

      default:
        break;
      }
    }


    void PetscSolverFeti::createFetiInexactKsp()
    {
      FUNCNAME("PetscSolverFeti::createFetiInexactKsp()");

      // === Init solver ===

      int localSize =
        localDofMap.getRankDofs() +
        primalDofMap.getRankDofs() +
        lagrangeMap.getRankDofs();

      int globalSize =
        nGlobalOverallInterior +
        primalDofMap.getOverallDofs() +
        lagrangeMap.getOverallDofs();

      fetiInexactData.matBB = &(subdomain->getMatInterior());
      fetiInexactData.matBPi = &(subdomain->getMatInteriorCoarse());
      fetiInexactData.matPiB = &(subdomain->getMatCoarseInterior());
      fetiInexactData.matPiPi = &(subdomain->getMatCoarse());
      fetiInexactData.mat_lagrange = &mat_lagrange;

      createVec(localDofMap, fetiInexactData.tmp_vec_b0);
      createVec(localDofMap, fetiInexactData.tmp_vec_b1);

      MatCreateShell(domainComm,
                     localSize, localSize, globalSize, globalSize,
                     &fetiInexactData, &mat_feti);

      MatShellSetOperation(mat_feti, MATOP_MULT,
                           (void(*)(void))petscMultMatFetiInexact);

      KSPCreate(domainComm, &ksp_feti);
#if (PETSC_VERSION_MINOR >= 5)
      KSPSetOperators(ksp_feti, mat_feti, mat_feti);
#else
      KSPSetOperators(ksp_feti, mat_feti, mat_feti, SAME_NONZERO_PATTERN);
#endif
      KSPSetOptionsPrefix(ksp_feti, "feti_");
      KSPSetType(ksp_feti, KSPGMRES);
      KSPSetTolerances(ksp_feti, 0, 1e-8, 1e+3, 1000);
      KSPSetFromOptions(ksp_feti);


      // === Init preconditioner ===

      fetiInexactPreconData.ksp_schur = ksp_schur_primal;
      fetiInexactPreconData.ksp_interior = subdomain->getSolver();
      fetiInexactPreconData.matPiB = &(subdomain->getMatCoarseInterior());
      fetiInexactPreconData.matBPi = &(subdomain->getMatInteriorCoarse());
      fetiInexactPreconData.mat_lagrange = &mat_lagrange;
      createVec(localDofMap, fetiInexactPreconData.tmp_vec_b0);

      KSPCreate(domainComm, &(fetiInexactPreconData.ksp_pc_feti));
#if (PETSC_VERSION_MINOR >= 5)
      KSPSetOperators(fetiInexactPreconData.ksp_pc_feti, mat_lagrange, mat_lagrange);
#else
      KSPSetOperators(fetiInexactPreconData.ksp_pc_feti, mat_lagrange, mat_lagrange, SAME_NONZERO_PATTERN);
#endif
      KSPGetPC(fetiInexactPreconData.ksp_pc_feti,
               &(fetiInexactPreconData.pc_feti));
      createFetiPreconLumped(fetiInexactPreconData.pc_feti);
      PCSetUp(fetiInexactPreconData.pc_feti);

      // === Setup pc object ===

      PC pc;
      KSPGetPC(ksp_feti, &pc);
      PCSetType(pc, PCSHELL);
      PCShellSetApply(pc, pcInexactFetiShell);
      PCShellSetContext(pc, &fetiInexactPreconData);
    }


    void PetscSolverFeti::createFetiInexactReducedKsp()
    {
      FUNCNAME("PetscSolverFeti::createFetiInexactReducedKsp()");

      ERROR_EXIT("Not yet implemented!\n");
    }


    void PetscSolverFeti::createFetiPreconLumped(PC pc)
    {
      FUNCNAME("PetscSolverFeti::createFetiPreconLumped()");

      FetiLumpedPreconData* lumpedData =
        (stokesMode ? &fetiInterfaceLumpedPreconData : &fetiLumpedPreconData);

      lumpedData->mat_lagrange_scaled = &mat_lagrange_scaled;
      lumpedData->mat_duals_duals = &mat_duals_duals;

      VecCreateMPI(meshDistributor->getMeshLevelData().getMpiComm(meshLevel),
                   localDofMap.getRankDofs(),
                   nGlobalOverallInterior, &(lumpedData->tmp_vec_b0));
      MatGetVecs(mat_duals_duals, PETSC_NULL,
                 &(lumpedData->tmp_vec_duals0));
      MatGetVecs(mat_duals_duals, PETSC_NULL,
                 &(lumpedData->tmp_vec_duals1));

      for (int component = 0; component < static_cast<int>(componentSpaces.size());
           component++)
      {
        if (stokesMode && component == pressureComponent)
          continue;

        DofMap& dualMap = dualDofMap[component].getMap();
        for (DofMap::iterator it = dualMap.begin(); it != dualMap.end(); ++it)
        {
          DegreeOfFreedom d = it->first;
          int matIndexLocal = localDofMap.getLocalMatIndex(component, d);
          int matIndexDual = dualDofMap.getLocalMatIndex(component, d);
          lumpedData->localToDualMap[matIndexLocal] = matIndexDual;
        }
      }

      if (stokesMode)
      {
        // === Create mass matrix solver ===

        const FiniteElemSpace* pressureFeSpace =
          componentSpaces[pressureComponent];

        // Create parallel DOF mapping in pressure space.
        ParallelDofMapping* massMapping = NULL;
        if (massMatrixSolver)
        {
          massMapping = massMatrixSolver->getDofMapping();
        }
        else
        {
          massMapping =
            new ParallelDofMapping(COMPONENT_WISE, true);
          massMapping->init(pressureFeSpace, pressureFeSpace);
          massMapping->setDofComms(meshDistributor->getDofComms(), meshLevel);
          massMapping->setMpiComm(meshDistributor->getMeshLevelData().getMpiComm(0));
        }
        (*massMapping)[0] = interfaceDofMap[pressureComponent];
        massMapping->update();

        DOFMatrix massMatrix(pressureFeSpace, pressureFeSpace);
        Operator op(pressureFeSpace, pressureFeSpace);
        Simple_ZOT zot;
        op.addTerm(&zot);
        massMatrix.assembleOperator(op);

        if (!massMatrixSolver)
        {
          massMatrixSolver = new PetscSolverGlobalMatrix("");
          massMatrixSolver->setKspPrefix("mass_");
          massMatrixSolver->setMeshDistributor(meshDistributor, meshLevel);
          massMatrixSolver->setDofMapping(massMapping);
        }

        massMatrixSolver->fillPetscMatrix(&massMatrix);

        int r, c;
        MatGetSize(massMatrixSolver->getMatInterior(), &r, &c);
        MatInfo info;
        MatGetInfo(massMatrixSolver->getMatInterior(), MAT_GLOBAL_SUM, &info);
        MSG("MASS MAT INFO:  size = %d x %d   nnz = %d\n",
            r, c, static_cast<int>(info.nz_used));

        fetiInterfaceLumpedPreconData.ksp_mass = massMatrixSolver->getSolver();


        // === Create tmp vectors ===

        createVec(localDofMap, fetiInterfaceLumpedPreconData.tmp_vec_b1);
        createVec(primalDofMap, fetiInterfaceLumpedPreconData.tmp_primal);
        fetiInterfaceLumpedPreconData.subSolver = subdomain;
      }

      // === Set PC object ===

      PCSetType(pc, PCSHELL);
      if (stokesMode)
      {
        PCShellSetContext(pc, static_cast<void*>(&fetiInterfaceLumpedPreconData));
        PCShellSetApply(pc, petscApplyFetiInterfaceLumpedPrecon);
      }
      else
      {
        PCShellSetContext(pc, static_cast<void*>(&fetiLumpedPreconData));
        PCShellSetApply(pc, petscApplyFetiLumpedPrecon);
      }
    }


    void PetscSolverFeti::createFetiPreconDirichlet(PC pc)
    {
      FUNCNAME("PetscSolverFeti::createFetiPreconDirichlet()");

      TEST_EXIT(subDomainIsLocal)
      ("Check for localDofMap.getLocalMatIndex, which should not work for multilevel FETI-DP!\n");

      TEST_EXIT(!stokesMode)
      ("Stokes mode does not yet support the Dirichlet precondition!\n");

      KSPCreate(PETSC_COMM_SELF, &ksp_interior);
#if (PETSC_VERSION_MINOR >= 5)
      KSPSetOperators(ksp_interior, mat_interior_interior, mat_interior_interior);
#else
      KSPSetOperators(ksp_interior, mat_interior_interior, mat_interior_interior, SAME_NONZERO_PATTERN);
#endif
      KSPSetOptionsPrefix(ksp_interior, "precon_interior_");
      KSPSetType(ksp_interior, KSPPREONLY);
      PC pc_interior;
      KSPGetPC(ksp_interior, &pc_interior);
      if (isSymmetric)
      {
        PCSetType(pc_interior, PCCHOLESKY);
        PCFactorSetMatSolverPackage(pc_interior, MATSOLVERMUMPS);
      }
      else
      {
        PCSetType(pc_interior, PCLU);
        PCFactorSetMatSolverPackage(pc_interior, MATSOLVERUMFPACK);
      }
      KSPSetFromOptions(ksp_interior);

      fetiDirichletPreconData.mat_lagrange_scaled = &mat_lagrange_scaled;
      fetiDirichletPreconData.mat_interior_interior = &mat_interior_interior;
      fetiDirichletPreconData.mat_duals_duals = &mat_duals_duals;
      fetiDirichletPreconData.mat_interior_duals = &mat_interior_duals;
      fetiDirichletPreconData.mat_duals_interior = &mat_duals_interior;
      fetiDirichletPreconData.ksp_interior = &ksp_interior;

      VecCreateMPI(meshDistributor->getMeshLevelData().getMpiComm(meshLevel),
                   localDofMap.getRankDofs(),
                   nGlobalOverallInterior, &(fetiDirichletPreconData.tmp_vec_b));
      MatGetVecs(mat_duals_duals, PETSC_NULL,
                 &(fetiDirichletPreconData.tmp_vec_duals0));
      MatGetVecs(mat_duals_duals, PETSC_NULL,
                 &(fetiDirichletPreconData.tmp_vec_duals1));
      MatGetVecs(mat_interior_interior, PETSC_NULL,
                 &(fetiDirichletPreconData.tmp_vec_interior));

      TEST_EXIT_DBG(subDomainIsLocal)
      ("Should not happen, check usage of localDofMap!\n");

      for (unsigned int component = 0; component < componentSpaces.size(); component++)
      {
        DofMap& dualMap = dualDofMap[component].getMap();
        for (DofMap::iterator it = dualMap.begin(); it != dualMap.end(); ++it)
        {
          DegreeOfFreedom d = it->first;
          int matIndexLocal = localDofMap.getLocalMatIndex(component, d);
          int matIndexDual = dualDofMap.getLocalMatIndex(component, d);
          fetiDirichletPreconData.localToDualMap[matIndexLocal] = matIndexDual;
        }
      }

      PCSetType(pc, PCSHELL);
      PCShellSetContext(pc, static_cast<void*>(&fetiDirichletPreconData));
      PCShellSetApply(pc, petscApplyFetiDirichletPrecon);

      // For the case, that we want to print the timings, we force the LU
      // factorization of the local problems to be done here explicitly.
      if (printTimings)
      {
        double wtime = MPI::Wtime();
        KSPSetUp(ksp_interior);
        KSPSetUpOnBlocks(ksp_interior);
        MPI::COMM_WORLD.Barrier();
        MSG("FETI-DP timing 08: %.5f seconds (factorization of Dirichlet preconditoner matrices)\n",
            MPI::Wtime() - wtime);
      }
    }


    void PetscSolverFeti::destroyFetiKsp()
    {
      FUNCNAME("PetscSolverFeti::destroyFetiKsp()");

      switch (fetiSolverType)
      {
      case EXACT:
        destroyFetiExactKsp();
        break;
      case INEXACT:
        destroyFetiInexactKsp();
        break;
      case INEXACT_REDUCED:
        destroyFetiInexactReducedKsp();
        break;
      default:
        ERROR_EXIT("Should not happen!\n");
      }
    }


    void PetscSolverFeti::destroyFetiExactKsp()
    {
      FUNCNAME("PetscSolverFeti::destroyFetiExactKsp()");

      // === Destroy FETI-DP solver object. ===

      fetiData.mat_lagrange = PETSC_NULL;
      fetiData.subSolver = NULL;
      fetiData.ksp_schur_primal = PETSC_NULL;

      VecDestroy(&fetiData.tmp_vec_b0);
      VecDestroy(&fetiData.tmp_vec_lagrange);
      VecDestroy(&fetiData.tmp_vec_primal0);

      if (augmentedLagrange)
      {
        fetiData.mat_augmented_lagrange = PETSC_NULL;
        VecDestroy(&fetiData.tmp_vec_primal1);
      }

      if (stokesMode)
      {
        VecDestroy(&fetiData.tmp_vec_b1);
        VecDestroy(&fetiData.tmp_vec_primal1);
        VecDestroy(&fetiData.tmp_vec_interface);
      }

      MatDestroy(&mat_feti);
      KSPDestroy(&ksp_feti);


      // === Destroy FETI-DP preconditioner object. ===

      switch (fetiPreconditioner)
      {
      case FETI_DIRICHLET:
        KSPDestroy(&ksp_interior);

        fetiDirichletPreconData.mat_lagrange_scaled = NULL;
        fetiDirichletPreconData.mat_interior_interior = NULL;
        fetiDirichletPreconData.mat_duals_duals = NULL;
        fetiDirichletPreconData.mat_interior_duals = NULL;
        fetiDirichletPreconData.mat_duals_interior = NULL;
        fetiDirichletPreconData.ksp_interior = NULL;

        VecDestroy(&fetiDirichletPreconData.tmp_vec_b);
        VecDestroy(&fetiDirichletPreconData.tmp_vec_duals0);
        VecDestroy(&fetiDirichletPreconData.tmp_vec_duals1);
        VecDestroy(&fetiDirichletPreconData.tmp_vec_interior);
        MatDestroy(&mat_lagrange_scaled);
        break;

      case FETI_LUMPED:
      {
        FetiLumpedPreconData& lumpedData =
          (stokesMode ? fetiInterfaceLumpedPreconData : fetiLumpedPreconData);

        lumpedData.mat_lagrange_scaled = NULL;
        lumpedData.mat_duals_duals = NULL;

        VecDestroy(&lumpedData.tmp_vec_b0);
        VecDestroy(&lumpedData.tmp_vec_duals0);
        VecDestroy(&lumpedData.tmp_vec_duals1);
      }
      break;
      default:
        break;
      }
    }


    void PetscSolverFeti::destroyFetiInexactKsp()
    {
      FUNCNAME("PetscSolverFeti::destroyFetiInexactKsp()");

      VecDestroy(&(fetiInexactData.tmp_vec_b0));
      VecDestroy(&(fetiInexactData.tmp_vec_b1));
      MatDestroy(&mat_feti);
      KSPDestroy(&ksp_feti);

      VecDestroy(&(fetiInexactPreconData.tmp_vec_b0));
      KSPDestroy(&(fetiInexactPreconData.ksp_pc_feti));
    }


    void PetscSolverFeti::destroyFetiInexactReducedKsp()
    {
      FUNCNAME("PetscSolverFeti::destroyFetiInexactReducedKsp()");
    }


    void PetscSolverFeti::createNullSpace()
    {
      FUNCNAME("PetscSolverFeti::createNullSpace()");

      if (!stokesMode)
        return;

      const FiniteElemSpace* pressureFeSpace = componentSpaces[pressureComponent];

      Vec ktest0, ktest1;
      createLocalVec(localDofMap, ktest0);
      createLocalVec(localDofMap, ktest1);
      DofMap& m = localDofMap[pressureComponent].getMap();
      for (DofMap::iterator it = m.begin(); it != m.end(); ++it)
      {
        if (dofMap[pressureFeSpace].isRankDof(it->first))
        {
          int index = localDofMap.getLocalMatIndex(pressureComponent, it->first);
          VecSetValue(ktest0, index, 1.0, INSERT_VALUES);
        }
      }
      VecAssemblyBegin(ktest0);
      VecAssemblyEnd(ktest0);
      MatMult(subdomain->getMatInterior(), ktest0, ktest1);

      PetscScalar* valarray;
      Vec ktest2, ktest3;
      VecGetArray(ktest1, &valarray);
      VecCreateMPIWithArray(PETSC_COMM_WORLD, 1,
                            localDofMap.getRankDofs(), nGlobalOverallInterior,
                            valarray, &ktest2);
      createVec(localDofMap, ktest3, nGlobalOverallInterior);

      Vec vecArray[2];
      createVec(interfaceDofMap, vecArray[0]);
      createVec(lagrangeMap, vecArray[1]);

      VecSet(vecArray[0], 1.0);
      MatMult(subdomain->getMatInteriorCoarse(1), vecArray[0], ktest3);
      VecAXPY(ktest3, 1.0, ktest2);
      MatMult(mat_lagrange_scaled, ktest3, vecArray[1]);
      VecScale(vecArray[1], -1.0);

      Vec nullSpaceBasis;
      VecCreateNest(domainComm, 2, PETSC_NULL, vecArray, &nullSpaceBasis);

#if (DEBUG != 0)
      PetscSolverFetiDebug::writeNullSpace(*this, nullSpaceBasis);
#endif

      MatNullSpace matNullSpace;
      MatNullSpaceCreate(domainComm, PETSC_FALSE, 1, &nullSpaceBasis,
                         &matNullSpace);
      MatSetNullSpace(mat_feti, matNullSpace);
      KSPSetNullSpace(ksp_feti, matNullSpace);
      MatNullSpaceDestroy(&matNullSpace);

      VecDestroy(&ktest0);
      VecDestroy(&ktest1);
      VecDestroy(&ktest2);
      VecDestroy(&ktest3);

      VecDestroy(&(vecArray[0]));
      VecDestroy(&(vecArray[1]));
      VecDestroy(&nullSpaceBasis);
    }


    void PetscSolverFeti::dbgMatrix(Matrix<DOFMatrix*>* mat)
    {
      FUNCNAME("PetscSolverFeti::dbgMatrix()");

      if (levelMode == 2 && meshLevel == 0)
      {
        MSG("WARNING: MAKE MORE GENERAL!\n");
        return;
      }

#if (DEBUG != 0)
      PetscInt nRow, nCol;
      MatGetLocalSize(subdomain->getMatInterior(), &nRow, &nCol);
      mpi::globalAdd(nRow);
      mpi::globalAdd(nCol);

      MatInfo minfo;
      MatGetInfo(subdomain->getMatInterior(), MAT_GLOBAL_SUM, &minfo);
      int nnz = static_cast<int>(minfo.nz_used);
      mpi::globalAdd(nnz);
      MSG("Interior matrices [%d x %d] nnz = %d\n", nRow, nCol, nnz);

      MatGetSize(subdomain->getMatCoarse(), &nRow, &nCol);
      MatGetInfo(subdomain->getMatCoarse(), MAT_GLOBAL_SUM, &minfo);
      MSG("Primal matrix [%d x %d] nnz = %d\n", nRow, nCol,
          static_cast<int>(minfo.nz_used));

      MatGetSize(subdomain->getMatCoarseInterior(), &nRow, &nCol);
      MatGetInfo(subdomain->getMatCoarseInterior(), MAT_GLOBAL_SUM, &minfo);
      MSG("Primal-Interior matrix [%d x %d] nnz = %d\n", nRow, nCol,
          static_cast<int>(minfo.nz_used));

      MatGetSize(subdomain->getMatInteriorCoarse(), &nRow, &nCol);
      MatGetInfo(subdomain->getMatInteriorCoarse(), MAT_GLOBAL_SUM, &minfo);
      MSG("Interior-Primal matrix [%d x %d] nnz = %d\n", nRow, nCol,
          static_cast<int>(minfo.nz_used));

      if (stokesMode)
      {
        MatGetSize(subdomain->getMatCoarse(1, 1), &nRow, &nCol);
        MatGetInfo(subdomain->getMatCoarse(1, 1), MAT_GLOBAL_SUM, &minfo);
        MSG("Interface matrix [%d x %d] nnz = %d\n", nRow, nCol,
            static_cast<int>(minfo.nz_used));

        MatGetSize(subdomain->getMatCoarseInterior(1), &nRow, &nCol);
        MatGetInfo(subdomain->getMatCoarseInterior(1), MAT_GLOBAL_SUM, &minfo);
        MSG("Interface-Interior matrix [%d x %d] nnz = %d\n", nRow, nCol,
            static_cast<int>(minfo.nz_used));

        MatGetSize(subdomain->getMatInteriorCoarse(1), &nRow, &nCol);
        MatGetInfo(subdomain->getMatInteriorCoarse(1), MAT_GLOBAL_SUM, &minfo);
        MSG("Interior-Interface matrix [%d x %d] nnz = %d\n", nRow, nCol,
            static_cast<int>(minfo.nz_used));

        MatGetSize(subdomain->getMatCoarse(1, 0), &nRow, &nCol);
        MatGetInfo(subdomain->getMatCoarse(1, 0), MAT_GLOBAL_SUM, &minfo);
        MSG("Interface-Primal matrix [%d x %d] nnz = %d\n", nRow, nCol,
            static_cast<int>(minfo.nz_used));

        MatGetSize(subdomain->getMatCoarse(0, 1), &nRow, &nCol);
        MatGetInfo(subdomain->getMatCoarse(0, 1), MAT_GLOBAL_SUM, &minfo);
        MSG("Primal-Interface matrix [%d x %d] nnz = %d\n", nRow, nCol,
            static_cast<int>(minfo.nz_used));
      }
#endif

      int writeInteriorMatrix = -1;
      Parameters::get("parallel->debug->write interior matrix",
                      writeInteriorMatrix);

      if (writeInteriorMatrix >= 0 &&
          writeInteriorMatrix == MPI::COMM_WORLD.Get_rank())
      {
        PetscViewer petscView;
        PetscViewerBinaryOpen(PETSC_COMM_SELF, "interior.mat",
                              FILE_MODE_WRITE, &petscView);
        MatView(subdomain->getMatInterior(), petscView);
        PetscViewerDestroy(&petscView);
      }


      bool checkInteriorMatrix = false;;
      Parameters::get("parallel->debug->check interior matrix",
                      checkInteriorMatrix);
      if (checkInteriorMatrix)
      {
        int nZeroRows = PetscSolverFetiDebug::testZeroRows(subdomain->getMatInterior());
        MSG("Interior matrix has %d zero rows!\n", nZeroRows);
      }

      bool printDirichlet = false;
      Parameters::get("parallel->debug->print dirichlet information",
                      printDirichlet);
      if (printDirichlet)
      {
        int nComponents = mat->getSize();
        for (int component = 0; component < nComponents; component++)
        {
          DOFMatrix* seqMat = (*mat)[component][component];
          if (!seqMat)
            continue;

          const FiniteElemSpace* feSpace = seqMat->getRowFeSpace();
          TEST_EXIT(feSpace)("Should not happen!\n");
          std::set<DegreeOfFreedom>& dirichletRows = seqMat->getDirichletRows();
          for (std::set<DegreeOfFreedom>::iterator it = dirichletRows.begin();
               it != dirichletRows.end(); ++it)
          {
            if (localDofMap[component].isSet(*it))
            {
              MSG("Dirichlet dof %d in component %d with interior mat index %d\n",
                  *it, component, localDofMap.getMatIndex(component, *it));
            }
          }
        }
      }


      int writeCoarseMatrix = 0;
      Parameters::get("parallel->debug->write coarse matrix",
                      writeCoarseMatrix);

      if (writeCoarseMatrix > 0)
      {
        PetscViewer petscView;
        PetscViewerBinaryOpen(PETSC_COMM_WORLD, "coarse.mat",
                              FILE_MODE_WRITE, &petscView);
        MatView(subdomain->getMatCoarse(), petscView);
        PetscViewerDestroy(&petscView);
      }


      int writeSchurPrimalMatrix = 0;
      Parameters::get("parallel->debug->write schur primal matrix",
                      writeSchurPrimalMatrix);
      if (writeSchurPrimalMatrix)
      {
        PetscViewer petscView;
        PetscViewerBinaryOpen(PETSC_COMM_WORLD, "schurprimal.mat",
                              FILE_MODE_WRITE, &petscView);
        MatView(mat_schur_primal, petscView);
        PetscViewerDestroy(&petscView);
      }
    }


    void PetscSolverFeti::recoverSolution(Vec& vec_sol_b,
                                          Vec& vec_sol_primal,
                                          SystemVector& vec)
    {
      FUNCNAME("PetscSolverFeti::recoverSolution()");

      // === Get local part of the solution for B variables. ===

      PetscScalar* localSolB;
      VecGetArray(vec_sol_b, &localSolB);

      PetscInt bsize;
      VecGetLocalSize(vec_sol_b, &bsize);

      // === Create scatter to get solutions of all primal nodes that are ===
      // === contained in rank's domain.                                  ===

      unsigned int nComponents = vec.getSize();

      vector<PetscInt> globalIsIndex, localIsIndex;
      globalIsIndex.reserve(primalDofMap.getLocalDofs());
      localIsIndex.reserve(primalDofMap.getLocalDofs());

      {
        int cnt = 0;
        for (unsigned int component = 0; component < nComponents; component++)
        {
          DofMap& dofMap = primalDofMap[component].getMap();
          for (DofMap::iterator it = dofMap.begin(); it != dofMap.end(); ++it)
          {
            globalIsIndex.push_back(primalDofMap.getMatIndex(component, it->first));
            localIsIndex.push_back(cnt++);
          }
        }

        TEST_EXIT_DBG(cnt == primalDofMap.getLocalDofs())
        ("Should not happen!\n");
      }

      IS globalIs, localIs;
      ISCreateGeneral(PETSC_COMM_SELF,
                      globalIsIndex.size(),
                      &(globalIsIndex[0]),
                      PETSC_USE_POINTER,
                      &globalIs);

      ISCreateGeneral(PETSC_COMM_SELF,
                      localIsIndex.size(),
                      &(localIsIndex[0]),
                      PETSC_USE_POINTER,
                      &localIs);

      Vec local_sol_primal;
      VecCreateSeq(PETSC_COMM_SELF, localIsIndex.size(), &local_sol_primal);

      VecScatter primalScatter;
      VecScatterCreate(vec_sol_primal, globalIs, local_sol_primal, localIs, &primalScatter);
      VecScatterBegin(primalScatter, vec_sol_primal, local_sol_primal,
                      INSERT_VALUES, SCATTER_FORWARD);
      VecScatterEnd(primalScatter, vec_sol_primal, local_sol_primal,
                    INSERT_VALUES, SCATTER_FORWARD);

      ISDestroy(&globalIs);
      ISDestroy(&localIs);
      VecScatterDestroy(&primalScatter);

      PetscScalar* localSolPrimal;
      VecGetArray(local_sol_primal, &localSolPrimal);

      // === And copy from PETSc local vectors to the DOF vectors. ===

      int cnt = 0;
      for (unsigned int component = 0; component < nComponents; component++)
      {
        DOFVector<double>& dofVec = *(vec.getDOFVector(component));

        for (DofMap::iterator it = localDofMap[component].getMap().begin();
             it != localDofMap[component].getMap().end(); ++it)
        {
          if (subDomainIsLocal)
          {
            int petscIndex = localDofMap.getLocalMatIndex(component, it->first);
            dofVec[it->first] = localSolB[petscIndex];
          }
          else
          {
            if (dofMapSubDomain[componentSpaces[component]].isRankDof(it->first))
            {
              int petscIndex = localDofMap.getLocalMatIndex(component, it->first);
              TEST_EXIT(petscIndex < bsize)("Should not happen!\n");
              dofVec[it->first] = localSolB[petscIndex];
            }
          }
        }

        for (DofMap::iterator it = primalDofMap[component].getMap().begin();
             it != primalDofMap[component].getMap().end(); ++it)
          dofVec[it->first] = localSolPrimal[cnt++];
      }

      VecRestoreArray(vec_sol_b, &localSolB);
      VecRestoreArray(local_sol_primal, &localSolPrimal);
      VecDestroy(&local_sol_primal);
    }


    void PetscSolverFeti::recoverInterfaceSolution(Vec& vecInterface, SystemVector& vec)
    {
      FUNCNAME("PetscSolverFeti::recoverInterfaceSolution()");

      if (!stokesMode)
        return;

      vector<PetscInt> globalIsIndex, localIsIndex;
      globalIsIndex.reserve(interfaceDofMap.getLocalDofs());
      localIsIndex.reserve(interfaceDofMap.getLocalDofs());

      int cnt = 0;
      DofMap& dofMap = interfaceDofMap[pressureComponent].getMap();
      for (DofMap::iterator it = dofMap.begin(); it != dofMap.end(); ++it)
      {
        globalIsIndex.push_back(interfaceDofMap.getMatIndex(pressureComponent,
                                it->first));
        localIsIndex.push_back(cnt++);
      }

      IS globalIs, localIs;
      ISCreateGeneral(PETSC_COMM_SELF,
                      globalIsIndex.size(),
                      &(globalIsIndex[0]),
                      PETSC_USE_POINTER,
                      &globalIs);

      ISCreateGeneral(PETSC_COMM_SELF,
                      localIsIndex.size(),
                      &(localIsIndex[0]),
                      PETSC_USE_POINTER,
                      &localIs);

      Vec local_sol_interface;
      VecCreateSeq(PETSC_COMM_SELF, localIsIndex.size(), &local_sol_interface);

      VecScatter interfaceScatter;
      VecScatterCreate(vecInterface, globalIs, local_sol_interface, localIs, &interfaceScatter);
      VecScatterBegin(interfaceScatter, vecInterface, local_sol_interface,
                      INSERT_VALUES, SCATTER_FORWARD);
      VecScatterEnd(interfaceScatter, vecInterface, local_sol_interface,
                    INSERT_VALUES, SCATTER_FORWARD);

      ISDestroy(&globalIs);
      ISDestroy(&localIs);
      VecScatterDestroy(&interfaceScatter);


      PetscScalar* localSolInterface;
      VecGetArray(local_sol_interface, &localSolInterface);

      // === And copy from PETSc local vectors to the DOF vectors. ===

      cnt = 0;
      DOFVector<double>& dofVec = *(vec.getDOFVector(pressureComponent));
      for (DofMap::iterator it = interfaceDofMap[pressureComponent].getMap().begin();
           it != interfaceDofMap[pressureComponent].getMap().end(); ++it)
      {
        dofVec[it->first] = localSolInterface[cnt++];
      }

      VecRestoreArray(local_sol_interface, &localSolInterface);
      VecDestroy(&local_sol_interface);
    }


    void PetscSolverFeti::fillPetscMatrix(Matrix<DOFMatrix*>* mat)
    {
      FUNCNAME("PetscSolverFeti::fillPetscMatrix()");

      // === Create all sets and indices. ===

      initialize();

      createDirichletData(*mat);

      createFetiData();

      // === Create matrices for the FETI-DP method. ===

      if (printTimings)
        MPI::COMM_WORLD.Barrier();
      double wtime = MPI::Wtime();

      subdomain->fillPetscMatrix(mat);


      // === SUPER TRICK ===

      if (meshLevel == 1)
      {
        MSG("START MAT TRICK!\n");
        mlSubdomain = new PetscSolverGlobalMatrix("");
        mlSubdomain->setSymmetric(isSymmetric);
        mlSubdomain->setHandleDirichletRows(dirichletMode == 0);
        mlSubdomain->setMeshDistributor(meshDistributor, meshLevel);
        mlSubdomain->init(componentSpaces, feSpaces);
        mlSubdomain->setDofMapping(interiorMap);
        mlSubdomain->setCoarseSpaceDofMapping(coarseSpaceMap[0]);
        mlSubdomain->fillPetscMatrix(mat);
        this->mat = mlSubdomain->getMat();
        MSG("END MAT TRICK!\n");
      }




      if (printTimings)
      {
        MPI::COMM_WORLD.Barrier();
        MSG("FETI-DP timing 02: %.5f seconds (creation of interior matrices)\n",
            MPI::Wtime() - wtime);

        // For the case, that we want to print the timings, we force the LU
        // factorization of the local problems to be done here explicitly.
        wtime = MPI::Wtime();
        KSPSetUp(subdomain->getSolver());
        KSPSetUpOnBlocks(subdomain->getSolver());
        MPI::COMM_WORLD.Barrier();
        MSG("FETI-DP timing 04: %.5f seconds (factorization of subdomain matrices)\n",
            MPI::Wtime() - wtime);
      }

      // === Create matrices for FETI-DP preconditioner. ===
      createPreconditionerMatrix(mat);

      // === Create and fill PETSc matrix for Lagrange constraints. ===
      createMatLagrange();

      // === ===
      createMatAugmentedLagrange();

      // === Create PETSc solver for the Schur complement on primal variables. ===
      createSchurPrimalKsp();

      // === Create PETSc solver for the FETI-DP operator. ===
      createFetiKsp();

      // === If required, run debug tests. ===
      dbgMatrix(mat);
    }


    void PetscSolverFeti::fillPetscRhs(SystemVector* vec)
    {
      FUNCNAME("PetscSolverFeti::fillPetscRhs()");

      subdomain->fillPetscRhs(vec);

      if (meshLevel == 1)
      {
        MSG("START VEC TRICK!\n");
        mlSubdomain->fillPetscRhs(vec);
        this->vecRhs = mlSubdomain->getVecRhs();
        MSG("END VEC TRICK!\n");
      }
    }


    void PetscSolverFeti::createPreconditionerMatrix(Matrix<DOFMatrix*>* mat)
    {
      FUNCNAME("PetscSolverFeti::createPreconditionerMatrix()");

      if (fetiPreconditioner == FETI_NONE && fetiSolverType == EXACT)
        return;

      double wtime = MPI::Wtime();
      int nRowsDual = dualDofMap.getRankDofs();

      MatCreateSeqAIJ(PETSC_COMM_SELF,
                      nRowsDual, nRowsDual, 100, PETSC_NULL,
                      &mat_duals_duals);
      MatSetOption(mat_duals_duals, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

      if (fetiPreconditioner == FETI_DIRICHLET)
      {
        int nRowsInterior = interiorDofMap.getRankDofs();

        MatCreateSeqAIJ(PETSC_COMM_SELF,
                        nRowsInterior, nRowsInterior, 100, PETSC_NULL,
                        &mat_interior_interior);
        MatSetOption(mat_interior_interior,
                     MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

        MatCreateSeqAIJ(PETSC_COMM_SELF,
                        nRowsInterior, nRowsDual, 100, PETSC_NULL,
                        &mat_interior_duals);
        MatSetOption(mat_interior_duals,
                     MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

        MatCreateSeqAIJ(PETSC_COMM_SELF,
                        nRowsDual, nRowsInterior, 100, PETSC_NULL,
                        &mat_duals_interior);
        MatSetOption(mat_duals_interior,
                     MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
      }

      // === Prepare traverse of sequentially created matrices. ===

      using mtl::tag::row;
      using mtl::tag::nz;
      using mtl::begin;
      using mtl::end;
      namespace traits = mtl::traits;
      typedef DOFMatrix::base_matrix_type Matrix;

      typedef traits::range_generator<row, Matrix>::type cursor_type;
      typedef traits::range_generator<nz, cursor_type>::type icursor_type;

      vector<int> colsLocal, colsLocalOther;
      vector<double> valuesLocal, valuesLocalOther;
      colsLocal.reserve(300);
      colsLocalOther.reserve(300);
      valuesLocal.reserve(300);
      valuesLocalOther.reserve(300);


      // === Traverse all sequentially created matrices and add the values to ===
      // === the global PETSc matrices.                                       ===

      int nComponents = mat->getSize();
      for (int rowComponent = 0; rowComponent < nComponents; rowComponent++)
      {
        for (int colComponent = 0; colComponent < nComponents; colComponent++)
        {
          DOFMatrix* dofMat = (*mat)[rowComponent][colComponent];

          if (!dofMat)
            continue;

          TEST_EXIT_DBG(dofMat->getRowFeSpace() == componentSpaces[rowComponent])
          ("Wrong matrix row FE space!\n");
          TEST_EXIT_DBG(dofMat->getColFeSpace() == componentSpaces[colComponent])
          ("Wrong matrix col FE space!!\n");

          if (stokesMode &&
              (rowComponent == pressureComponent ||
               colComponent == pressureComponent))
            continue;

          // 	const FiniteElemSpace *rowFeSpace = dofMat->getRowFeSpace();
          // 	const FiniteElemSpace *colFeSpace = dofMat->getColFeSpace();

          traits::col<Matrix>::type col(dofMat->getBaseMatrix());
          traits::const_value<Matrix>::type value(dofMat->getBaseMatrix());

          // Traverse all rows.
          for (cursor_type cursor = begin<row>(dofMat->getBaseMatrix()),
               cend = end<row>(dofMat->getBaseMatrix()); cursor != cend; ++cursor)
          {

            if (dirichletRows[rowComponent].count(cursor.value()))
              continue;

            if (isPrimal(rowComponent, cursor.value()))
              continue;

            if (fetiPreconditioner == FETI_DIRICHLET)
            {
              colsLocal.clear();
              colsLocalOther.clear();
              valuesLocal.clear();
              valuesLocalOther.clear();

              // Traverse all columns.
              for (icursor_type icursor = begin<nz>(cursor), icend = end<nz>(cursor);
                   icursor != icend; ++icursor)
              {

                if (dirichletRows[colComponent].count(col(*icursor)))
                  continue;

                if (isPrimal(colComponent, col(*icursor)))
                  continue;

                if (!isDual(rowComponent, cursor.value()))
                {
                  if (!isDual(colComponent, col(*icursor)))
                  {
                    int colIndex =
                      interiorDofMap.getLocalMatIndex(colComponent, col(*icursor));
                    colsLocal.push_back(colIndex);
                    valuesLocal.push_back(value(*icursor));
                  }
                  else
                  {
                    int colIndex =
                      dualDofMap.getLocalMatIndex(colComponent, col(*icursor));
                    colsLocalOther.push_back(colIndex);
                    valuesLocalOther.push_back(value(*icursor));
                  }
                }
                else
                {
                  if (!isDual(colComponent, col(*icursor)))
                  {
                    int colIndex =
                      interiorDofMap.getLocalMatIndex(colComponent, col(*icursor));
                    colsLocalOther.push_back(colIndex);
                    valuesLocalOther.push_back(value(*icursor));
                  }
                  else
                  {
                    int colIndex =
                      dualDofMap.getLocalMatIndex(colComponent, col(*icursor));
                    colsLocal.push_back(colIndex);
                    valuesLocal.push_back(value(*icursor));
                  }
                }
              }  // for each nnz in row

              if (!isDual(rowComponent, cursor.value()))
              {
                int rowIndex =
                  interiorDofMap.getLocalMatIndex(rowComponent, cursor.value());
                MatSetValues(mat_interior_interior, 1, &rowIndex, colsLocal.size(),
                             &(colsLocal[0]), &(valuesLocal[0]), INSERT_VALUES);

                if (colsLocalOther.size())
                  MatSetValues(mat_interior_duals, 1, &rowIndex, colsLocalOther.size(),
                               &(colsLocalOther[0]), &(valuesLocalOther[0]), INSERT_VALUES);
              }
              else
              {
                int rowIndex =
                  dualDofMap.getLocalMatIndex(rowComponent, cursor.value());
                MatSetValues(mat_duals_duals, 1, &rowIndex, colsLocal.size(),
                             &(colsLocal[0]), &(valuesLocal[0]), INSERT_VALUES);
                if (colsLocalOther.size())
                  MatSetValues(mat_duals_interior, 1, &rowIndex, colsLocalOther.size(),
                               &(colsLocalOther[0]), &(valuesLocalOther[0]), INSERT_VALUES);
              }
            }

            if (fetiPreconditioner == FETI_LUMPED || fetiSolverType == INEXACT)
            {
              if (!isDual(rowComponent, cursor.value()))
                continue;

              colsLocal.clear();
              valuesLocal.clear();

              // Traverse all columns.
              for (icursor_type icursor = begin<nz>(cursor), icend = end<nz>(cursor);
                   icursor != icend; ++icursor)
              {
                if (dirichletRows[colComponent].count(col(*icursor)))
                  continue;

                if (!isDual(colComponent, col(*icursor)))
                  continue;

                int colIndex =
                  dualDofMap.getLocalMatIndex(colComponent, col(*icursor));
                colsLocal.push_back(colIndex);
                valuesLocal.push_back(value(*icursor));
              }

              int rowIndex = dualDofMap.getLocalMatIndex(rowComponent, cursor.value());
              MatSetValues(mat_duals_duals, 1, &rowIndex, colsLocal.size(),
                           &(colsLocal[0]), &(valuesLocal[0]), INSERT_VALUES);
            }
          }
        }
      }

      // === Start global assembly procedure for preconditioner matrices. ===

      if (fetiPreconditioner == FETI_LUMPED || fetiSolverType == INEXACT)
      {
        MatAssemblyBegin(mat_duals_duals, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(mat_duals_duals, MAT_FINAL_ASSEMBLY);
      }

      if (fetiPreconditioner == FETI_DIRICHLET)
      {
        MatAssemblyBegin(mat_interior_interior, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(mat_interior_interior, MAT_FINAL_ASSEMBLY);

        MatAssemblyBegin(mat_duals_duals, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(mat_duals_duals, MAT_FINAL_ASSEMBLY);

        MatAssemblyBegin(mat_interior_duals, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(mat_interior_duals, MAT_FINAL_ASSEMBLY);

        MatAssemblyBegin(mat_duals_interior, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(mat_duals_interior, MAT_FINAL_ASSEMBLY);
      }

      if (printTimings)
      {
        MPI::COMM_WORLD.Barrier();
        MSG("FETI-DP timing 03: %.5f seconds (creation of preconditioner matrices)\n",
            MPI::Wtime() - wtime);
      }
    }


    void PetscSolverFeti::solveFeti(Vec& rhsInterior, Vec& rhsCoarse,
                                    Vec& solInterior, Vec& solCoarse)
    {
      FUNCNAME("PetscSolverFeti::solveFeti()");

      switch (fetiSolverType)
      {
      case EXACT:
        solveFetiExact(rhsInterior, rhsCoarse, solInterior, solCoarse);
        break;
      case INEXACT:
        solveFetiInexact(rhsInterior, rhsCoarse, solInterior, solCoarse);
        break;
      case INEXACT_REDUCED:
        solveFetiInexactReduced(rhsInterior, rhsCoarse, solInterior, solCoarse);
        break;
      default:
        ERROR_EXIT("Should not happen!\n");
      }
    }


    void PetscSolverFeti::solveFetiExact(Vec& rhsInterior, Vec& rhsCoarse,
                                         Vec& solInterior, Vec& solCoarse)
    {
      FUNCNAME("PetscSolverFeti::solveFetiExact()");

      // === Some temporary vectors. ===

      Vec tmp_b1;
      VecCreateMPI(meshDistributor->getMeshLevelData().getMpiComm(meshLevel),
                   localDofMap.getRankDofs(),
                   nGlobalOverallInterior, &tmp_b1);

      Vec tmp_primal1;
      createVec(primalDofMap, tmp_primal1);

      Vec tmp_lagrange;
      MatGetVecs(mat_lagrange, PETSC_NULL, &tmp_lagrange);

      // === Create RHS and solution vectors. ===

      Vec vecRhs, vecSol;
      Vec vecRhsLagrange, vecSolLagrange;
      MatGetVecs(mat_lagrange, PETSC_NULL, &vecRhsLagrange);
      MatGetVecs(mat_lagrange, PETSC_NULL, &vecSolLagrange);

      vecRhs = vecRhsLagrange;
      vecSol = vecSolLagrange;

      VecDuplicate(vecSol, &fetiKspData.draft);

      // === Create reduced RHS ===

      double wtime = MPI::Wtime();

      //    d = L inv(K_BB) f_B - L inv(K_BB) K_BPi inv(S_PiPi) [f_Pi - K_PiB inv(K_BB) f_B]
      // vecRhs = L * inv(K_BB) * f_B
      subdomain->solveGlobal(rhsInterior, solInterior);

      MatMult(mat_lagrange, solInterior, vecRhsLagrange);

      // solCoarse = M_PiB * inv(K_BB) * f_B
      MatMult(subdomain->getMatCoarseInterior(), solInterior, solCoarse);

      // solCoarse = f_Pi - M_PiB * inv(K_BB) * f_B
      VecAXPBY(solCoarse, 1.0, -1.0, rhsCoarse);

      double wtime2 = MPI::Wtime();
      // solCoarse = inv(S_PiPi) (f_Pi - M_PiB * inv(K_BB) * f_B)
      KSPSolve(ksp_schur_primal, solCoarse, solCoarse);
      if (printTimings)
      {
        MSG("FETI-DP timing 09a: %.5f seconds (create rhs vector)\n",
            MPI::Wtime() - wtime2);
      }

      MatMult(subdomain->getMatInteriorCoarse(), solCoarse, solInterior);
      subdomain->solveGlobal(solInterior, solInterior);
      MatMult(mat_lagrange, solInterior, tmp_lagrange);
      VecAXPY(vecRhsLagrange, -1.0, tmp_lagrange);

      if (printTimings)
      {
        MPI::COMM_WORLD.Barrier();
        MSG("FETI-DP timing 09: %.5f seconds (create rhs vector)\n",
            MPI::Wtime() - wtime);
        wtime = MPI::Wtime();
        FetiTimings::reset();
      }

      // === Optionally run some debug procedures on the FETI-DP system. ===
      PetscSolverFetiDebug::debugFeti(*this, vecRhs);

      // === Solve with FETI-DP operator. ===
      KSPSolve(ksp_feti, vecRhs, vecSol);

      if (printTimings)
      {
        MPI::COMM_WORLD.Barrier();
        MSG("FETI-DP timing 10: %.5f seconds (application of FETI-DP operator)\n",
            MPI::Wtime() - wtime);
        wtime = MPI::Wtime();

        MSG("FETI-DP timing 10a: %.5f [ %.5f %.5f ] seconds (FETI-DP KSP Solve)\n",
            FetiTimings::fetiSolve, FetiTimings::fetiSolve01, FetiTimings::fetiSolve02);
        MSG("FETI-DP timing 10b: %.5f seconds (FETI-DP KSP Solve)\n",
            FetiTimings::fetiPreconditioner);
      }


      // === Solve for u_primals. ===

      MatMultTranspose(mat_lagrange, vecSol, solInterior);
      VecAYPX(solInterior, -1.0, rhsInterior);

      subdomain->solveGlobal(solInterior, tmp_b1);
      MatMult(subdomain->getMatCoarseInterior(), tmp_b1, solCoarse);
      VecAYPX(solCoarse, -1.0, rhsCoarse);

      KSPSolve(ksp_schur_primal, solCoarse, solCoarse);

      // === Solve for u_b. ===

      MatMult(subdomain->getMatInteriorCoarse(), solCoarse, tmp_b1);
      VecAXPY(solInterior, -1.0, tmp_b1);
      subdomain->solveGlobal(solInterior, solInterior);


      // === Print timings and delete data. ===

      if (printTimings)
      {
        MPI::COMM_WORLD.Barrier();
        MSG("FETI-DP timing 11: %.5f seconds (Inner solve and solution recovery)\n",
            MPI::Wtime() - wtime);
      }

      VecDestroy(&vecRhs);
      VecDestroy(&vecSol);
      VecDestroy(&tmp_b1);
      VecDestroy(&tmp_lagrange);
      VecDestroy(&tmp_primal1);
    }


    void PetscSolverFeti::solveFetiInexact(Vec& rhsInterior, Vec& rhsCoarse,
                                           Vec& solInterior, Vec& solCoarse)
    {
      FUNCNAME("PetscSolverFeti::solveFetiInexact()");

      Vec tmpLagrange0, tmpLagrange1;
      createVec(lagrangeMap, tmpLagrange0);
      createVec(lagrangeMap, tmpLagrange1);
      VecSet(tmpLagrange0, 0.0);
      VecSet(tmpLagrange1, 0.0);

      Vec nestVecRhs[3];
      nestVecRhs[0] = rhsInterior;
      nestVecRhs[1] = rhsCoarse;
      nestVecRhs[2] = tmpLagrange0;

      Vec nestVecSol[3];
      nestVecSol[0] = solInterior;
      nestVecSol[1] = solCoarse;
      nestVecSol[2] = tmpLagrange1;

      Vec vecRhs, vecSol;
      VecCreateNest(domainComm, 3, PETSC_NULL, nestVecRhs, &vecRhs);
      VecCreateNest(domainComm, 3, PETSC_NULL, nestVecSol, &vecSol);

      KSPSolve(ksp_feti, vecRhs, vecSol);

      VecDestroy(&vecRhs);
      VecDestroy(&vecSol);
      VecDestroy(&tmpLagrange0);
      VecDestroy(&tmpLagrange1);
    }


    void PetscSolverFeti::solveFetiInexactReduced(Vec& rhsInterior, Vec& rhsCoarse,
        Vec& solInterior, Vec& solCoarse)
    {
      FUNCNAME("PetscSolverFeti::solveFetiInexactReduced()");

      ERROR_EXIT("Not yet implemented!\n");
    }


    void PetscSolverFeti::destroyMatrixData()
    {
      FUNCNAME("PetscSolverFeti::destroyMatrixData()");

      MatDestroy(&mat_lagrange);

      if (augmentedLagrange)
        MatDestroy(&mat_augmented_lagrange);

      // === Destroy preconditioner data structures. ===

      if (fetiPreconditioner != FETI_NONE)
        MatDestroy(&mat_duals_duals);

      if (fetiPreconditioner == FETI_DIRICHLET)
      {
        MatDestroy(&mat_interior_interior);
        MatDestroy(&mat_interior_duals);
        MatDestroy(&mat_duals_interior);
      }

      destroySchurPrimalKsp();

      destroyFetiKsp();

      subdomain->destroyMatrixData();
    }


    void PetscSolverFeti::destroyVectorData()
    {
      FUNCNAME("PetscSolverFeti::destroyVectorData()");

      subdomain->destroyVectorData();
    }


    void PetscSolverFeti::solvePetscMatrix(SystemVector& vec, AdaptInfo& adaptInfo)
    {
      FUNCNAME("PetscSolverFeti::solvePetscMatrix()");

      Vec solInterior;
      VecCreateMPI(meshDistributor->getMeshLevelData().getMpiComm(meshLevel),
                   localDofMap.getRankDofs(),
                   nGlobalOverallInterior,
                   &solInterior);

      Vec solCoarse;
      createVec(primalDofMap, solCoarse);

      solveFeti(subdomain->getVecRhsInterior(),
                subdomain->getVecRhsCoarse(),
                solInterior,
                solCoarse);

      // === And recover AMDiS solution vectors. ===

      recoverSolution(solInterior, solCoarse, vec);

      VecDestroy(&solInterior);
      VecDestroy(&solCoarse);

      MeshDistributor::globalMeshDistributor->synchVector(vec);
    }


    void PetscSolverFeti::solveGlobal(Vec& rhs, Vec& sol)
    {
      FUNCNAME("PetscSolverFeti::solveGlobal()");

      Vec rhsInterior, rhsCoarse, solInterior, solCoarse;
      VecCreateMPI(domainComm,
                   localDofMap.getRankDofs(),
                   nGlobalOverallInterior,
                   &rhsInterior);
      createVec(primalDofMap, rhsCoarse);
      VecDuplicate(rhsInterior, &solInterior);
      VecDuplicate(rhsCoarse, &solCoarse);

      int offset = 0;
      {
        int domainLocal = 0, nSuperLocal = 0;
        if (domainComm.Get_rank() == 0)
          domainLocal = interiorMap->getOverallDofs();

        mpi::getDofNumbering(meshDistributor->getMpiComm(meshLevel - 1),
                             domainLocal, offset, nSuperLocal);

        int tmp = 0;
        if (domainComm.Get_rank() == 0)
          tmp = offset;

        domainComm.Allreduce(&tmp, &offset, 1, MPI_INT, MPI_SUM);
      }

      vector<int> localFromRhs, coarseFromRhs;
      vector<int> rhsToLocal, rhsToCoarse;

      int nComponents = componentSpaces.size();
      for (int i = 0; i < nComponents; i++)
      {
        DofMap& dMap = localDofMap[i].getMap();

        for (DofMap::iterator it = dMap.begin(); it != dMap.end(); ++it)
        {
          int matL = localDofMap.getMatIndex(i, it->first) + rStartInterior;
          int matI = interiorMap->getMatIndex(i, it->first) + offset;

          localFromRhs.push_back(matL);
          rhsToLocal.push_back(matI);
        }
      }

      for (int i = 0; i < nComponents; i++)
      {
        DofMap& dMap = primalDofMap[i].getMap();

        for (DofMap::iterator it = dMap.begin(); it != dMap.end(); ++it)
        {
          int matL = primalDofMap.getMatIndex(i, it->first);
          int matI = interiorMap->getMatIndex(i, it->first) + offset;

          coarseFromRhs.push_back(matL);
          rhsToCoarse.push_back(matI);
        }
      }

      copyVec(rhs, rhsInterior, rhsToLocal, localFromRhs);
      copyVec(rhs, rhsCoarse, rhsToCoarse, coarseFromRhs);

      solveFeti(rhsInterior, rhsCoarse, solInterior, solCoarse);

      copyVec(solInterior, sol, localFromRhs, rhsToLocal);
      copyVec(solCoarse, sol, coarseFromRhs, rhsToCoarse);

      MPI::COMM_WORLD.Barrier();

      VecDestroy(&rhsInterior);
      VecDestroy(&rhsCoarse);
      VecDestroy(&solInterior);
      VecDestroy(&solCoarse);
    }
  }
} // end namespaces
