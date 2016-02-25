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


#include "parallel/ParallelDofMapping.hpp"
#include "parallel/StdMpi.hpp"

namespace AMDiS
{
  namespace Parallel
  {

    using namespace std;


    void DofToMatIndex::getReverse(int rowIndex, int& component, int& dofIndex) const
    {
      for (map<int, boost::container::flat_map<DegreeOfFreedom, int>>::const_iterator it0 = data.begin();
           it0 != data.end(); ++it0)
        for (boost::container::flat_map<DegreeOfFreedom, int>::const_iterator it1 = it0->second.begin();
             it1 != it0->second.end(); ++it1)
          if (it1->second == rowIndex)
          {
            component = it0->first;
            dofIndex = it1->first;
            return;
          }

      component = -1;
      dofIndex = -1;
    }


    ComponentDofMap::ComponentDofMap()
      : dofComm(NULL),
        feSpace(NULL),
        globalMapping(false)
    {
      clear();
    }


    void ComponentDofMap::clear()
    {
      dofMap.clear();

      nonRankDofs.clear();

      nRankDofs = 0;
      nLocalDofs = 0;
      nOverallDofs = 0;
      rStartDofs = 0;
    }


    void ComponentDofMap::update()
    {
      // === Compute local indices for all rank owned DOFs. ===

      for (DofMap::iterator it = dofMap.begin(); it != dofMap.end(); ++it)
        if (it->second.local == -1 && nonRankDofs.count(it->first) == 0)
          it->second.local = nRankDofs++;

      // === Compute number of local and global DOFs in the mapping. ===

      nOverallDofs = 0;
      rStartDofs = 0;
      mpi::getDofNumbering(mpiComm, nRankDofs, rStartDofs, nOverallDofs);

      // === If required, compute also the global indices. ===

      if (globalMapping)
      {
        computeGlobalMapping();
        computeNonLocalIndices();
      }
    }


    void ComponentDofMap::computeGlobalMapping()
    {
      for (DofMap::iterator it = dofMap.begin(); it != dofMap.end(); ++it)
        it->second.global = it->second.local + rStartDofs;
    }


    void ComponentDofMap::computeNonLocalIndices()
    {
      FUNCNAME_DBG("ComponentDofMap::computeNonLocalIndices()");

      TEST_EXIT_DBG(dofComm)("No DOF communicator defined!\n");

      // === Send all global indices of DOFs that are owned by the rank to all ===
      // === other ranks that also include this DOF.                           ===

      StdMpi<vector<int>> stdMpi(mpiComm);

      for (DofComm::Iterator it(dofComm->getSendDofs(), feSpace);
           !it.end(); it.nextRank())
      {
        int rank = it.getRank();

        // if (meshLevel > 0)
        //   rank = levelData->mapRank(rank, 0, meshLevel);

        for (; !it.endDofIter(); it.nextDof())
        {
          TEST_EXIT_DBG(dofMap.count(it.getDofIndex()))
          ("DOF index not in the dofmap. Something is wrong.\n");
          if (dofMap.count(it.getDofIndex()) &&
              !nonRankDofs.count(it.getDofIndex()))
            stdMpi.getSendData(rank).
            push_back(dofMap[it.getDofIndex()].global);
        }
      }

      stdMpi.updateSendDataSize();


      // === Check from which ranks this rank must receive some data. ===

      for (DofComm::Iterator it(dofComm->getRecvDofs(), feSpace);
           !it.end(); it.nextRank())
      {
        bool recvFromRank = false;
        for (; !it.endDofIter(); it.nextDof())
        {
          if (nonRankDofs.count(it.getDofIndex()))
          {
            recvFromRank = true;
            break;
          }
        }

        if (recvFromRank)
        {
          int rank = it.getRank();
          // if (meshLevel > 0)
          //   rank = levelData->mapRank(rank, 0, meshLevel);

          stdMpi.recv(rank);
        }
      }


      // === Start communication to exchange global indices. ===

      stdMpi.startCommunication();


      // === And set the global indices for all DOFs that are not owned by rank. ===

      for (DofComm::Iterator it(dofComm->getRecvDofs(), feSpace);
           !it.end(); it.nextRank())
      {
        int rank = it.getRank();
        // if (meshLevel > 0)
        //   rank = levelData->mapRank(rank, 0, meshLevel);

        int i = 0;
        for (; !it.endDofIter(); it.nextDof())
        {
          TEST_EXIT_DBG(nonRankDofs.count(it.getDofIndex()))("Something is wrong.\n");

          if (nonRankDofs.count(it.getDofIndex()))
            dofMap[it.getDofIndex()].global = stdMpi.getRecvData(rank)[i++];
        }
      }

    }


    void FeSpaceData::init(vector<const FiniteElemSpace*>& f0,
                           vector<const FiniteElemSpace*>& f1,
                           bool globalMapping)
    {
      componentSpaces = f0;
      feSpaces = f1;

      for (vector<const FiniteElemSpace*>::iterator it = feSpaces.begin();
           it != feSpaces.end(); ++it)
      {
        if (componentData.count(*it))
          componentData.find(*it)->second.clear();
        else
          componentData.insert(make_pair(*it, ComponentDofMap()));

        componentData[*it].setFeSpace(*it);
        componentData[*it].setMesh((*it)->getMesh());
        componentData[*it].setGlobalMapping(globalMapping);
      }
    }


    void ComponentData::init(vector<const FiniteElemSpace*>& f0,
                             vector<const FiniteElemSpace*>& f1,
                             bool globalMapping)
    {
      componentSpaces = f0;
      feSpaces = f1;

      for (unsigned int component = 0; component < componentSpaces.size(); component++)
      {
        if (componentData.count(component))
          componentData.find(component)->second.clear();
        else
          componentData.insert(make_pair(component, ComponentDofMap()));

        componentData[component].setFeSpace(componentSpaces[component]);
        componentData[component].setMesh(componentSpaces[component]->getMesh());
        componentData[component].setGlobalMapping(globalMapping);
      }
    }


    ParallelDofMapping::ParallelDofMapping(DofMappingMode mode,
                                           bool matIndexFromGlobal)
      :  globalMapping(true),
         needMatIndexFromGlobal(matIndexFromGlobal),
         nRankDofs(1),
         nLocalDofs(1),
         nOverallDofs(1),
         rStartDofs(1),
         mode(mode)
    {
      switch (mode)
      {
      case COMPONENT_WISE:
        data = new ComponentData();
        break;
      case FESPACE_WISE:
        data = new FeSpaceData();
        break;
      }

      nRankDofs = -1;
      nLocalDofs = -1;
      nOverallDofs = -1;
      rStartDofs = -1;
    }

    void ParallelDofMapping::init(vector<const FiniteElemSpace*>& fe,
                                  vector<const FiniteElemSpace*>& uniqueFe,
                                  bool b)
    {

      globalMapping = b;
      componentSpaces = fe;
      feSpaces = uniqueFe;

      data->init(fe, uniqueFe, globalMapping);
    }

    void ParallelDofMapping::clear(Mesh* mesh)
    {

      for (ComponentIterator& it = data->getIteratorData(); !it.end(); it.next())
        if (!mesh || it->getMesh() == mesh)
          it->clear();

      nRankDofs = -1;
      nLocalDofs = -1;
      nOverallDofs = -1;
      rStartDofs = -1;

      if (!mesh)
        dofToMatIndex.clear();
      else
        for (unsigned int component = 0; component < componentSpaces.size();
             component++)
        {
          if (mesh == (*data)[component].getMesh())
            dofToMatIndex.clear(component);
        }
    }


    void ParallelDofMapping::setMpiComm(MPI::Intracomm& m)
    {
      mpiComm = m;

      for (ComponentIterator& it = data->getIteratorData(); !it.end(); it.next())
        it->setMpiComm(m);
    }


    void ParallelDofMapping::setDofComms(std::map<Mesh*, MultiLevelDofComm>& dofComms, int level)
    {
      FUNCNAME("ParallelDofMapping::setDofComms()");

      for (ComponentIterator& it = data->getIteratorData(); !it.end(); it.next())
      {

        TEST_EXIT(dofComms.find(const_cast<Mesh*>(it->getMesh())) != dofComms.end())
        ("DofComm and ParallelDofMapping not match.\n");

        it->setDofComm(dofComms[const_cast<Mesh*>(it->getMesh())][level]);
      }
    }


    int ParallelDofMapping::computeRankDofs()
    {
      int result = 0;
      for (ComponentIterator& it = data->getIteratorComponent(); !it.end(); it.next())
        result += it->nRankDofs;

      return result;
    }


    int ParallelDofMapping::computeLocalDofs()
    {
      int result = 0;
      for (ComponentIterator& it = data->getIteratorComponent(); !it.end(); it.next())
        result += it->nLocalDofs;

      return result;
    }


    int ParallelDofMapping::computeOverallDofs()
    {
      int result = 0;
      for (ComponentIterator& it = data->getIteratorComponent(); !it.end(); it.next())
        result += it->nOverallDofs;

      return result;
    }


    int ParallelDofMapping::computeStartDofs()
    {
      int result = 0;
      for (ComponentIterator& it = data->getIteratorComponent(); !it.end(); it.next())
        result += it->rStartDofs;

      return result;
    }

    void ParallelDofMapping::update(Mesh* mesh)
    {
      // First, update all FE space DOF mappings.
      for (ComponentIterator& it = data->getIteratorData(); !it.end(); it.next())
        if (!mesh || it->getMesh() == mesh)
          it->update();

      // Compute all numbers from this mappings.
      nRankDofs = computeRankDofs();
      nLocalDofs = computeLocalDofs();
      nOverallDofs = computeOverallDofs();
      rStartDofs = computeStartDofs();
      // And finally, compute the matrix indices.
      computeMatIndex(needMatIndexFromGlobal);
    }

    void ParallelDofMapping::updateMatIndex()
    {
      computeMatIndex(needMatIndexFromGlobal);
    }


    void ParallelDofMapping::computeMatIndex(bool globalIndex)
    {
      FUNCNAME_DBG("ParallelDofMapping::computeMatIndex()");

      dofToMatIndex.clear();

      // The offset is always added to the local matrix index. The offset for the
      // DOFs in the first FE spaces is the smalled global index of a DOF that is
      // owned by the rank.
      int offset = rStartDofs;

      // === Create the matrix indices for all component FE spaces. ===

      for (unsigned int component = 0; component < componentSpaces.size();
           component++)
      {

        // Traverse all DOFs of the FE space and create for all rank owned DOFs
        // a matrix index.
        DofMap& dofMap = (*data)[component].getMap();
        for (DofMap::iterator it = dofMap.begin(); it != dofMap.end(); ++it)
        {
          if ((*data)[component].isRankDof(it->first))
          {
            int globalMatIndex = it->second.local + offset;
            if (globalIndex)
              dofToMatIndex.add(component, it->second.global, globalMatIndex);
            else
              dofToMatIndex.add(component, it->first, globalMatIndex);
          }
        }
        // Increase the offset for the next FE space by the number of DOFs owned
        // by the rank in the current FE space.
        offset += (*data)[component].nRankDofs;

        // If there are no non local DOFs, continue with the next FE space.
        if (!globalMapping)
          continue;


        DofComm* dofComm = &((*data)[component].getDofComm());
        TEST_EXIT_DBG(dofComm != NULL)("No communicator given!\n");

        // === Communicate the matrix indices for all DOFs that are on some ===
        // === interior boundaries.                                         ===

        StdMpi<vector<DegreeOfFreedom>> stdMpi(mpiComm);
        for (DofComm::Iterator it(dofComm->getSendDofs(), componentSpaces[component]);
             !it.end(); it.nextRank())
        {
          vector<DegreeOfFreedom> sendGlobalDofs;

          for (; !it.endDofIter(); it.nextDof())
          {
            TEST_EXIT_DBG((*data)[component].isRankDof(it.getDofIndex()))
            ("Send dof %d is not a rank dof\n", it.getDofIndex());

            if (dofMap.count(it.getDofIndex()))
            {
              if (globalIndex)
                sendGlobalDofs.push_back(dofToMatIndex.get(component,
                                         dofMap[it.getDofIndex()].global));
              else
                sendGlobalDofs.push_back(dofToMatIndex.get(component, it.getDofIndex()));
            }
          }

          int rank = it.getRank();
          // if (meshLevel > 0)
          //   rank = levelData->mapRank(rank, 0, meshLevel);
          stdMpi.send(rank, sendGlobalDofs);
        }

        for (DofComm::Iterator it(dofComm->getRecvDofs(), componentSpaces[component]);
             !it.end(); it.nextRank())
        {
          int rank = it.getRank();
          // if (meshLevel > 0)
          //   rank = levelData->mapRank(rank, 0, meshLevel);
          stdMpi.recv(rank);
        }

        stdMpi.startCommunication();

        for (DofComm::Iterator it(dofComm->getRecvDofs(), componentSpaces[component]);
             !it.end(); it.nextRank())
        {
          int rank = it.getRank();
          // if (meshLevel > 0)
          //   rank = levelData->mapRank(rank, 0, meshLevel);

          int counter = 0;
          for (; !it.endDofIter(); it.nextDof())
          {
            if (dofMap.count(it.getDofIndex()))
            {
              DegreeOfFreedom d = stdMpi.getRecvData(rank)[counter++];
              if (globalIndex)
                dofToMatIndex.add(component, dofMap[it.getDofIndex()].global, d);
              else
                dofToMatIndex.add(component, it.getDofIndex(), d);
            }
          }
        }

        // === Communicate DOFs on periodic boundaries. ===

        stdMpi.clear();

        for (DofComm::Iterator it(dofComm->getPeriodicDofs(), componentSpaces[component]);
             !it.end(); it.nextRank())
        {
          vector<DegreeOfFreedom> sendGlobalDofs;

          for (; !it.endDofIter(); it.nextDof())
            if (dofMap.count(it.getDofIndex()))
            {
              if (globalIndex)
              {
                sendGlobalDofs.push_back(dofMap[it.getDofIndex()].global);
                sendGlobalDofs.push_back(dofToMatIndex.get(component,
                                         dofMap[it.getDofIndex()].global));
              }
              else
              {
                sendGlobalDofs.push_back(dofToMatIndex.get(component, it.getDofIndex()));
              }
            }

          int rank = it.getRank();
          // if (meshLevel > 0)
          //   rank = levelData->mapRank(rank, 0, meshLevel);

          stdMpi.send(rank, sendGlobalDofs);
          stdMpi.recv(rank);
        }

        stdMpi.startCommunication();

        for (DofComm::Iterator it(dofComm->getPeriodicDofs(), componentSpaces[component]);
             !it.end(); it.nextRank())
        {

          // 	int rank = it.getRank();
          // if (meshLevel > 0)
          //   rank = levelData->mapRank(rank, 0, meshLevel);

          size_t counter = 0;

          for (; !it.endDofIter(); it.nextDof())
          {
            if (dofMap.count(it.getDofIndex()))
            {
              if (globalIndex)
              {
                TEST_EXIT_DBG(counter + 2 <= stdMpi.getRecvData(it.getRank()).size())
                ("Should not happen!\n");

                dofToMatIndex.add(component,
                                  stdMpi.getRecvData(it.getRank())[counter],
                                  stdMpi.getRecvData(it.getRank())[counter + 1]);
                counter += 2;
              }
              else
              {
                dofToMatIndex.add(component,
                                  it.getDofIndex(),
                                  stdMpi.getRecvData(it.getRank())[counter++]);
              }
            }
          }
        }
      }
    }


#ifndef HAVE_PARALLEL_MTL4
    void ParallelDofMapping::createIndexSet(IS& is,
                                            int firstComponent,
                                            int nComponents)
    {
      FUNCNAME_DBG("ParallelDofMapping::createIndexSet()");

      int firstRankDof = -1;
      ComponentDofMap& compMap = (*data)[firstComponent];
      DofMap& dofMap = compMap.getMap();
      for (DofMap::iterator it = dofMap.begin(); it != dofMap.end(); ++it)
      {
        if (compMap.isRankDof(it->first))
        {
          if (needMatIndexFromGlobal)
            firstRankDof = it->second.global;
          else
            firstRankDof = it->first;
          break;
        }
      }

      TEST_EXIT_DBG(firstRankDof >= 0)("No rank DOF found!\n");

      int nRankRows = 0;
      for (int i = firstComponent; i < firstComponent + nComponents; i++)
        nRankRows += (*data)[i].nRankDofs;

      int firstMatIndex = dofToMatIndex.get(firstComponent, firstRankDof);
      ISCreateStride(mpiComm, nRankRows, firstMatIndex, 1, &is);
    }
#endif


    void ParallelDofMapping::printInfo()
    {
      FUNCNAME("ParallelDofMapping::printInfo()");

      MSG("=== Parallel DOF mapping debug information ===\n");
      if (mode == COMPONENT_WISE)
      {
        MSG("  mapping is defined by component numbers!\n");
      }
      else
      {
        MSG("  mapping is defined by FE spaces!\n");
      }
      MSG("  matrix index is based on global DOF indices: %d\n",
          needMatIndexFromGlobal);

      MSG("  nRankDofs = %d   nLocalDofs = %d   nOverallDofs = %d  rStartDofs = %d\n",
          nRankDofs, nLocalDofs, nOverallDofs, rStartDofs);

      int nComponents = componentSpaces.size();
      int nFeSpaces = feSpaces.size();
      MSG("  number of components: %d   number of different FE spaces: %d\n",
          nComponents, nFeSpaces);

      for (int i = 0; i < nComponents; i++)
      {
        MSG("  component %d:\n", i);
        MSG("    dof-to-mat-index has %d mappings\n", dofToMatIndex.getSize(i));
        if (dofToMatIndex.getSize(i) > 0)
        {
          MSG("    dof-to-mat-index starts with (%d -> %d) and ends with (%d -> %d)\n",
              dofToMatIndex.getData(i).begin()->first,
              dofToMatIndex.getData(i).begin()->second,
              (dofToMatIndex.getData(i).end() - 1)->first,
              (dofToMatIndex.getData(i).end() - 1)->second);

          DofToMatIndex::MapType::iterator it;

          std::stringstream fn;
          fn << "map_indices-p" << MPI::COMM_WORLD.Get_rank() << ".dat";
          std::ofstream out(fn.str().c_str(), ios::out);
          for (it = dofToMatIndex.getData(i).begin(); it != dofToMatIndex.getData(i).end(); it++)
          {
            // 	  MSG("    dof-to-mat-index: (%d -> %d)\n", it->first, it->second);
            out << it->first << " " << it->second << "\n";
          }
          out.close();
        }
      }
    }
  }
}
