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



/** \file MeshLevelData.h */

#ifndef AMDIS_MESH_LEVEL_DATA_H
#define AMDIS_MESH_LEVEL_DATA_H


#include <iostream>
#include <set>
#include <vector>
#include <mpi.h>
#include "Global.h"

namespace AMDiS
{
  namespace Parallel
  {

    class MeshLevelData
    {
    public:
      MeshLevelData();

      void init(std::set<int>& neighbourRanks);

      void addLevel(std::set<int>& ranksInDomain, int domainId);

      void addLevelMode1();

      // Writes all data of this object to an output stream.
      void serialize(std::ostream& out);

      // Reads the object data from an input stream.
      void deserialize(std::istream& in);

      void print();

      std::set<int>& getLevelRanks(int level)
      {
        TEST_EXIT_DBG(level < nLevel)("Should not happen!\n");

        return levelRanks[level];
      }

      std::set<int>& getLevelNeighbours(int level)
      {
        TEST_EXIT_DBG(level < nLevel)("Should not happen!\n");

        return levelNeighbours[level];
      }

      inline int getNumberOfLevels()
      {
        return nLevel;
      }

      MPI::Intracomm& getMpiComm(int level)
      {
        FUNCNAME("MeshLevelData::getMpiComm()");

        TEST_EXIT_DBG(level < nLevel)
        ("Asked for level %d, but defined only for %d levels!\n", level, nLevel);

        return mpiComms[level];
      }

      MPI::Group& getMpiGroup(int level)
      {
        TEST_EXIT_DBG(level < nLevel)("Should not happen!\n");

        return mpiGroups[level];
      }

      int getLevelId(int level, int rank)
      {
        FUNCNAME("MeshLevelData::getLevelId()");

        TEST_EXIT_DBG(level < nLevel)
        ("Wrong level number: %d  nLevel = %d\n", level, nLevel);

        return mpiComms[level].Get_rank();
      }

      int mapRank(int fromRank, int fromLevel, int toLevel);

      bool rankInSubdomain(int rank, int level)
      {
        TEST_EXIT_DBG(level < nLevel)("Should not happen!\n");

        return static_cast<bool>(levelRanks[level].count(rank));
      }

      bool isLevelBoundary(int rank, int level);

      bool rankInLevel(int rank, int level);

    protected:
      int nLevel;

      std::vector<std::set<int>> levelRanks;

      std::vector<std::set<int>> levelNeighbours;

      std::vector<MPI::Intracomm> mpiComms;

      std::vector<MPI::Group> mpiGroups;
    };

  }
}

#endif
