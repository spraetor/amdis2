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


#include <boost/lexical_cast.hpp>
#include "parallel/MeshLevelData.h"
#include "Global.h"

namespace AMDiS { namespace Parallel {

  using namespace std;


  MeshLevelData::MeshLevelData()
  {
    std::set<int> neighbours;
    init(neighbours);
  }


  void MeshLevelData::init(std::set<int> &neighbourRanks)
  {
    levelRanks.resize(1);
    
    int mpiSize = MPI::COMM_WORLD.Get_size();
    for (int i = 0; i < mpiSize; i++)
      levelRanks[0].insert(i);
    nLevel = 1;
    
    levelNeighbours.resize(1);
    levelNeighbours[0] = neighbourRanks;

    mpiComms.resize(1);
    mpiComms[0] = MPI::COMM_WORLD;

    mpiGroups.resize(1);
    mpiGroups[0] = mpiComms[0].Get_group();
  }


  void MeshLevelData::addLevel(std::set<int> &ranksInDomain, int domainId)
  {
    FUNCNAME("MeshLevelData()::addLevel()");

    TEST_EXIT(nLevel >= 1)("Mesh level structure is not initialized()");

    levelRanks.push_back(ranksInDomain);
    nLevel++;

    levelNeighbours.resize(2);    
    levelNeighbours[1].clear();
    for (std::set<int>::iterator it = levelNeighbours[0].begin(); 
	 it != levelNeighbours[0].end(); ++it)
      if (levelRanks[1].count(*it) == 0)
	levelNeighbours[1].insert(*it);

    mpiComms.resize(2);
    mpiComms[1] = mpiComms[0].Split(domainId, mpiComms[0].Get_rank());
    
    mpiGroups.resize(2);
    mpiGroups[1] = mpiComms[1].Get_group();
  }


  void MeshLevelData::addLevelMode1()
  {
    nLevel++;
    mpiComms.push_back(MPI::COMM_SELF);
    mpiGroups.push_back(MPI::COMM_SELF.Get_group());

    levelRanks.resize(nLevel);
    levelNeighbours.resize(nLevel);

    levelRanks[nLevel - 1].insert(MPI::COMM_WORLD.Get_rank());
  }


  void MeshLevelData::serialize(ostream &out)
  {
    ERROR_EXIT("Not yet implemented!\n");
  }

  
  void MeshLevelData::deserialize(istream &in)
  {
    ERROR_EXIT("Not yet implemented!\n");
  }


  void MeshLevelData::print()
  {
    FUNCNAME("MeshLevelData::print()");

    using boost::lexical_cast;

    MSG("Print mesh level structure with %d levels: \n", nLevel);

    for (int i = 0; i < nLevel; i++) {
      string ranks = "ranks in level " + lexical_cast<string>(i) + ":";
      for (std::set<int>::iterator it = levelRanks[i].begin(); 
	   it != levelRanks[i].end(); ++it)
	ranks += " " + lexical_cast<string>(*it);

      string neighbours = "neighbours in level " + lexical_cast<string>(i) + ": ";
      for (std::set<int>::iterator it = levelNeighbours[i].begin(); 
	   it != levelNeighbours[i].end(); ++it)
	neighbours += " " + lexical_cast<string>(*it);

      if (ranks.length() < 250)
	MSG("  %s\n", ranks.c_str());
      else
	MSG("  ranks string to long!\n");

      if (neighbours.length() < 250)
	MSG("%s\n", neighbours.c_str());
      else
	MSG("  neighbours string to long!\n");
    }
  }


  int MeshLevelData::mapRank(int fromRank, int fromLevel, int toLevel)
  {
    FUNCNAME("MeshLevelData::mapRank()");

    if (fromLevel == toLevel)
      return fromRank;
    
    int toRank = -1;
    
    MPI::Group::Translate_ranks(mpiGroups[fromLevel], 1, &fromRank, 
				mpiGroups[toLevel], &toRank);
    
    TEST_EXIT_DBG(toRank != MPI::UNDEFINED)
      ("Cannot map rank %d from level %d to level %d\n",
       fromRank, fromLevel, toLevel);
    
    return toRank;			
  }


  bool MeshLevelData::isLevelBoundary(int rank, int level)
  {
    FUNCNAME("MeshLevelData::isLevelBoundary()");

    TEST_EXIT_DBG(level < nLevel)("Should not happen!\n");

    if (nLevel == 1)
      return true;

    if (level + 1 < nLevel && mpiComms[level + 1] == MPI::COMM_SELF)
      return true;

    return static_cast<bool>(levelRanks[level].count(rank));
  }


  bool MeshLevelData::rankInLevel(int rank, int level)
  {
    return static_cast<bool>(levelRanks[level].count(rank));
  }

} }
