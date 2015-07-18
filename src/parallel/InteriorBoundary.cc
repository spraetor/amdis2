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


#include "parallel/InteriorBoundary.h"
#include "parallel/ElementObjectDatabase.h"
#include "parallel/MeshDistributor.h"
#include "FiniteElemSpace.h"
#include "BasisFunction.h"
#include "Serializer.h"
#include "VertexVector.h"

namespace AMDiS { namespace Parallel {

  using namespace std;

  void InteriorBoundary::create(MeshLevelData &levelData,
				int level,
				ElementObjectDatabase &elObjDb)
  { 
    FUNCNAME("InteriorBoundary::create()");

    own.clear();
    other.clear();
    periodic.clear();

    macroElIndexMap = elObjDb.getElIndexMap();
    Mesh *mesh = elObjDb.getMesh();
    TEST_EXIT_DBG(mesh)("Should not happen!\n");
    TEST_EXIT_DBG(level < levelData.getNumberOfLevels())
      ("Should not happen!\n");

    MPI::Intracomm mpiComm = levelData.getMpiComm(level);

    if (mpiComm == MPI::COMM_SELF)
      return;

//     int levelMpiRank = mpiComm.Get_rank();
    int globalMpiRank = MPI::COMM_WORLD.Get_rank();
    std::set<int> levelRanks = levelData.getLevelRanks(level);

    // === Create interior boundary data structure. ===
    
    for (int geoPos = 0; geoPos < mesh->getDim(); geoPos++) {
      GeoIndex geoIndex = INDEX_OF_DIM(geoPos, mesh->getDim());

      elObjDb.resetIterator();
      while (elObjDb.iterate(geoIndex)) {
	flat_map<int, ElementObjectData>& objData = elObjDb.getIterateData();

	// Test, if this is a boundary object of this rank.
	if (!(objData.count(globalMpiRank) && objData.size() > 1))
	  continue;

	// Test, if the boundary object defines an interior boundary within the
	// ranks of the MPI group. If not, go to next element.
	bool boundaryWithinMpiGroup = false;
	if (levelData.getNumberOfLevels() == 1) {
	  boundaryWithinMpiGroup = true;
	} else {
	  TEST_EXIT_DBG(level + 1 < levelData.getNumberOfLevels())
	    ("Level = %d    but number of levels = %d\n",
	     level, levelData.getNumberOfLevels());

	  for (flat_map<int, ElementObjectData>::iterator it = objData.begin();
	       it != objData.end(); ++it) {
	    if (!levelData.rankInLevel(it->first, level + 1)) {
	      boundaryWithinMpiGroup = true;
	      break;
	    }
	  }
	}
	
	if (!boundaryWithinMpiGroup)
	  continue;

	int owner = elObjDb.getIterateOwner(level);
	ElementObjectData& rankBoundEl = objData[globalMpiRank];
	
	AtomicBoundary bound;
	bound.maxLevel = elObjDb.getIterateMaxLevel();
	bound.rankObj.elIndex = rankBoundEl.elIndex;
	bound.rankObj.elType = elObjDb.getElementType(rankBoundEl.elIndex);
	bound.rankObj.subObj = geoIndex;
	bound.rankObj.ithObj = rankBoundEl.ithObject;

	if (geoIndex == FACE) {
	  for (int edgeNo = 0; edgeNo < 3; edgeNo++) {
	    Element* el = elObjDb.getElementPtr(bound.rankObj.elIndex, mesh);
	    int edgeOfFace = 
	      el->getEdgeOfFace(bound.rankObj.ithObj, edgeNo);
	    
	    bound.rankObj.excludedSubstructures.push_back(make_pair(EDGE, edgeOfFace));
	  }
	}
	

	if (owner == globalMpiRank) {
	  for (flat_map<int, ElementObjectData>::iterator it2 = objData.begin();
	       it2 != objData.end(); ++it2) {
	    if (it2->first == globalMpiRank)
	      continue;

	    if (!levelData.rankInLevel(it2->first, level))
	      continue;

	    bound.neighObj.elIndex = it2->second.elIndex;
	    bound.neighObj.elType = elObjDb.getElementType(it2->second.elIndex);
	    bound.neighObj.subObj = geoIndex;
	    bound.neighObj.ithObj = it2->second.ithObject;
	    
	    bound.type = INTERIOR;

	    int rankInLevel = levelData.mapRank(it2->first, 0, level);
	    AtomicBoundary& b = getNewOwn(rankInLevel);
	    b = bound;
	    if (geoIndex == EDGE)
	      b.neighObj.reverseMode = 
		elObjDb.getEdgeReverseMode(rankBoundEl, it2->second);
	    if (geoIndex == FACE)
	      b.neighObj.reverseMode = 
		elObjDb.getFaceReverseMode(rankBoundEl, it2->second);
	  }
	  
	} else {
	  TEST_EXIT_DBG(objData.count(owner) == 1)
	    ("Should not happen!\n");
	  
	  ElementObjectData& ownerBoundEl = objData[owner];
  
	  bound.neighObj.elIndex = ownerBoundEl.elIndex;
	  bound.neighObj.elType = -1;
	  bound.neighObj.subObj = geoIndex;
	  bound.neighObj.ithObj = ownerBoundEl.ithObject;
	  
	  bound.type = INTERIOR;

	  int rankInLevel = levelData.mapRank(owner, 0, level);	  
	  AtomicBoundary& b = getNewOther(rankInLevel);
	  b = bound;	    
	  if (geoIndex == EDGE)
	    b.rankObj.reverseMode =
	      elObjDb.getEdgeReverseMode(rankBoundEl, ownerBoundEl);
	  if (geoIndex == FACE)
	    b.rankObj.reverseMode = 
	      elObjDb.getFaceReverseMode(rankBoundEl, ownerBoundEl);
	}
      }
    }

    // === Create periodic boundary data structure. ===

    int removePeriodicBoundary = 0;
    Parameters::get("parallel->remove periodic boundary", removePeriodicBoundary);


    for (PerBoundMap<DegreeOfFreedom>::iterator it = elObjDb.getPeriodicVertices().begin();
	 it != elObjDb.getPeriodicVertices().end(); ++it) {
      if (elObjDb.isInRank(it->first.first, globalMpiRank) == false)
	continue;

      ElementObjectData& perDofEl0 = 
	elObjDb.getElementsInRank(it->first.first)[globalMpiRank];

      for (flat_map<int, ElementObjectData>::iterator elIt = elObjDb.getElementsInRank(it->first.second).begin();
	   elIt != elObjDb.getElementsInRank(it->first.second).end(); ++elIt) {

	int otherElementRank = elIt->first;
	ElementObjectData& perDofEl1 = elIt->second;

	AtomicBoundary bound;
	bound.rankObj.elIndex = perDofEl0.elIndex;
	bound.rankObj.elType = elObjDb.getElementType(perDofEl0.elIndex);
	bound.rankObj.subObj = VERTEX;
	bound.rankObj.ithObj = perDofEl0.ithObject;

	bound.neighObj.elIndex = perDofEl1.elIndex;
	bound.neighObj.elType = elObjDb.getElementType(perDofEl1.elIndex);
	bound.neighObj.subObj = VERTEX;
	bound.neighObj.ithObj = perDofEl1.ithObject;


	if (removePeriodicBoundary) {
	  bound.type = INTERIOR;
	  
	  flat_map<int, ElementObjectData> objData = 
	    elObjDb.getElementsInRank(it->first.first);
	  objData.insert(elObjDb.getElementsInRank(it->first.second).begin(),
			 elObjDb.getElementsInRank(it->first.second).end());
	  
	  int owner = std::max(elObjDb.getOwner(it->first.first, level),
			       elObjDb.getOwner(it->first.second, level));

	  // Consider the following configuration of four elements partitioned
	  // to four ranks:
	  //
	  //  |--------|--------|(*)
	  //  |       /|       /|
	  //  |  3   / |  0   / |
	  //  |     /  |     /  |
	  //  |    /   |    /   |
	  //  |   /    |   /    |
	  //  |  /  2  |  /  1  |
	  //  | /      | /      |
	  //  |/-------|/-------|
	  //
	  // We are in now interested in vertex (*). For the first, rank 1 owns
	  // the interior boundary with rank 0 on this vertex. When we remove 
	  // periodic boundaries, which are applied on the left and right outer
	  // boundary, rank 3 becomes the owner of this vertex. Rank 0 and rank 1
	  // will have an interior boundary with rank 3 on this vertex, but rank
	  // 0 and rank 1 have no more interior boundaries on vertex (*). Thus we
	  // have to search and remove for this scenario.
	  
	  if (owner != globalMpiRank) {
	    removeOwn(bound.rankObj);
	    removeOther(bound.rankObj);
	  }


	  // And than we can add the boundary to either own or other part of the
	  // interior database.

	  if (owner == globalMpiRank) {
	    int rankInLevel = levelData.mapRank(otherElementRank, 0, level);
	    AtomicBoundary& b = getNewOwn(rankInLevel);	  
	    b = bound;	  
	  } else {
	    ElementObjectData& ownerBoundEl = objData[owner];
	    bound.neighObj.elIndex = ownerBoundEl.elIndex;
	    bound.neighObj.elType = -1;
	    bound.neighObj.ithObj = ownerBoundEl.ithObject;

	    int rankInLevel = levelData.mapRank(owner, 0, level);
	    AtomicBoundary& b = getNewOther(rankInLevel);
	    b = bound;
	  }

	} else {
	  bound.type = it->second;
	  
	  int rankInLevel = levelData.mapRank(otherElementRank, 0, level);
	  AtomicBoundary& b = getNewPeriodic(rankInLevel);
	  b = bound;	  
	}
      }
    }


    for (PerBoundMap<DofEdge>::iterator it = elObjDb.getPeriodicEdges().begin();
	 it != elObjDb.getPeriodicEdges().end(); ++it) {
      if (elObjDb.isInRank(it->first.first, globalMpiRank) == false)
	continue;

      ElementObjectData& perEdgeEl0 = 
	elObjDb.getElementsInRank(it->first.first)[globalMpiRank];

      for (flat_map<int, ElementObjectData>::iterator elIt = elObjDb.getElementsInRank(it->first.second).begin();
 	   elIt != elObjDb.getElementsInRank(it->first.second).end(); ++elIt) {
      
	int otherElementRank = elIt->first;
	ElementObjectData& perEdgeEl1 = elIt->second;

	AtomicBoundary bound;	    	    
	bound.rankObj.elIndex = perEdgeEl0.elIndex;
	bound.rankObj.elType = elObjDb.getElementType(perEdgeEl0.elIndex);
	bound.rankObj.subObj = EDGE;
	bound.rankObj.ithObj = perEdgeEl0.ithObject;
	
	bound.neighObj.elIndex = perEdgeEl1.elIndex;
	bound.neighObj.elType = elObjDb.getElementType(perEdgeEl1.elIndex);
	bound.neighObj.subObj = EDGE;
	bound.neighObj.ithObj = perEdgeEl1.ithObject;

	if (removePeriodicBoundary) {
	  bound.type = INTERIOR;

	  flat_map<int, ElementObjectData> objData = 
	    elObjDb.getElementsInRank(it->first.first);
	  objData.insert(elObjDb.getElementsInRank(it->first.second).begin(),
			 elObjDb.getElementsInRank(it->first.second).end());
	  
	  int owner = std::max(elObjDb.getOwner(it->first.first, level),
			       elObjDb.getOwner(it->first.second, level));
	  ElementObjectData& rankBoundEl = objData[globalMpiRank];

	  // See comments in the same part of code for VERTEX
	  if (owner != globalMpiRank) {
	    removeOwn(bound.rankObj);
	    removeOther(bound.rankObj);
	  }

	  if (owner == globalMpiRank) {
	    int rankInLevel = levelData.mapRank(owner, 0, level);
	    AtomicBoundary& b = getNewOwn(rankInLevel); 
	    b = bound;	  
	  } else {
	    int rankInLevel = levelData.mapRank(owner, 0, level);
	    AtomicBoundary& b = getNewOther(rankInLevel);
	    b = bound;
	    
	    ElementObjectData& ownerBoundEl = objData[owner];    
	    b.neighObj.elIndex = ownerBoundEl.elIndex;
	    b.neighObj.elType = -1;
	    b.neighObj.ithObj = ownerBoundEl.ithObject;
	    b.rankObj.reverseMode = 
	      elObjDb.getEdgeReverseMode(rankBoundEl, ownerBoundEl);
	  }
	} else {
	  bound.type = it->second;

	  int rankInLevel = levelData.mapRank(otherElementRank, 0, level);
	  AtomicBoundary& b = getNewPeriodic(rankInLevel);
	  b = bound;
	  
	  if (globalMpiRank > otherElementRank)
	    b.neighObj.reverseMode = 
	      elObjDb.getEdgeReverseMode(perEdgeEl0, perEdgeEl1);
	  else
	    b.rankObj.reverseMode = 
	      elObjDb.getEdgeReverseMode(perEdgeEl0, perEdgeEl1);
	}
      }
    }


    for (PerBoundMap<DofFace>::iterator it = elObjDb.getPeriodicFaces().begin();
	 it != elObjDb.getPeriodicFaces().end(); ++it) {
      if (elObjDb.isInRank(it->first.first, globalMpiRank) == false)
	continue;

      TEST_EXIT_DBG(elObjDb.getElements(it->first.first).size() == 1)
 	("Should not happen!\n");
      TEST_EXIT_DBG(elObjDb.getElements(it->first.second).size() == 1)
 	("Should not happen!\n");

      ElementObjectData& perFaceEl0 = 
	elObjDb.getElementsInRank(it->first.first)[globalMpiRank];

      for (flat_map<int, ElementObjectData>::iterator elIt = elObjDb.getElementsInRank(it->first.second).begin();
 	   elIt != elObjDb.getElementsInRank(it->first.second).end(); ++elIt) {
      
	int otherElementRank = elIt->first;
	ElementObjectData& perFaceEl1 = elIt->second;

	AtomicBoundary bound;	    	    
	bound.rankObj.elIndex = perFaceEl0.elIndex;
	bound.rankObj.elType = elObjDb.getElementType(perFaceEl0.elIndex);
	bound.rankObj.subObj = FACE;
	bound.rankObj.ithObj = perFaceEl0.ithObject;
	
	bound.neighObj.elIndex = perFaceEl1.elIndex;
	bound.neighObj.elType = elObjDb.getElementType(perFaceEl1.elIndex);
	bound.neighObj.subObj = FACE;
	bound.neighObj.ithObj = perFaceEl1.ithObject;
	
	if (removePeriodicBoundary) {
	  bound.type = INTERIOR;

	  flat_map<int, ElementObjectData> objData = 
	    elObjDb.getElementsInRank(it->first.first);
	  objData.insert(elObjDb.getElementsInRank(it->first.second).begin(),
			 elObjDb.getElementsInRank(it->first.second).end());
	  
	  int owner = std::max(elObjDb.getOwner(it->first.first, level),
			       elObjDb.getOwner(it->first.second, level));
	  ElementObjectData& rankBoundEl = objData[globalMpiRank];

	  if (owner == globalMpiRank) {
	    int rankInLevel = levelData.mapRank(otherElementRank, 0, level);
	    AtomicBoundary& b = getNewOwn(rankInLevel);
	    b = bound;	  
	  } else {
	    int rankInLevel = levelData.mapRank(owner, 0, level);
	    AtomicBoundary& b = getNewOther(rankInLevel);
	    b = bound;
	    
	    ElementObjectData& ownerBoundEl = objData[owner];    
	    b.neighObj.elIndex = ownerBoundEl.elIndex;
	    b.neighObj.elType = -1;
	    b.neighObj.ithObj = ownerBoundEl.ithObject;
	    b.rankObj.reverseMode = 
	      elObjDb.getFaceReverseMode(rankBoundEl, ownerBoundEl);
	  }
	} else {
	  bound.type = it->second;
	  
	  int rankInLevel = levelData.mapRank(otherElementRank, 0, level);
	  AtomicBoundary& b = getNewPeriodic(rankInLevel);
	  b = bound;
	  
	  if (globalMpiRank > otherElementRank)
	    b.neighObj.reverseMode = 
	      elObjDb.getFaceReverseMode(perFaceEl0, perFaceEl1);
	  else
	    b.rankObj.reverseMode = 
	      elObjDb.getFaceReverseMode(perFaceEl0, perFaceEl1);
	}
      }
    }


    // === Once we have this information, we must care about the order of the ===
    // === atomic bounds in the three boundary handling object. Eventually    ===
    // === all the boundaries have to be in the same order on both ranks that ===
    // === share the bounday.                                                 ===

    StdMpi<vector<AtomicBoundary> > stdMpi(mpiComm);
    stdMpi.send(own);
    stdMpi.recv(other);
    stdMpi.startCommunication();


    // === The information about all neighbouring boundaries has been         ===
    // === received. So the rank tests if its own atomic boundaries are in    ===
    // === the same order. If not, the atomic boundaries are swaped to the    ===
    // === correct order.                                                     ===

    for (RankToBoundMap::iterator rankIt = other.begin();
	 rankIt != other.end(); ++rankIt) {

      int rank = rankIt->first;

      // === We have received from "rank" the ordered list of element       ===
      // === indices. Now, we have to sort the corresponding list in this   ===
      // === rank to get the same order.                                    ===
     
      for (unsigned int j = 0; j < rankIt->second.size(); j++) {

	// If the expected object is not at place, search for it.

	BoundaryObject &receivedBound = 
	  stdMpi.getRecvData()[rank][j].rankObj;

	if ((rankIt->second)[j].neighObj != receivedBound) {
	  unsigned int k = j + 1;

	  for (; k < rankIt->second.size(); k++) {
 	    if ((rankIt->second)[k].neighObj == receivedBound)
	      break;
	  }

	  // The element must always be found, because the list is just in
	  // another order.
	  TEST_EXIT_DBG(k < rankIt->second.size())
	    ("Cannot find element %d/%d/%d received from rank %d!\n",
	     receivedBound.elIndex, 
	     receivedBound.subObj, 
	     receivedBound.ithObj,
	     rank);

	  // Swap the current with the found element.
	  AtomicBoundary tmpBound = (rankIt->second)[k];
	  (rankIt->second)[k] = (rankIt->second)[j];
	  (rankIt->second)[j] = tmpBound;	
	}
      }
    }


    // === Do the same for the periodic boundaries. ===

    if (periodic.size() > 0) {
      stdMpi.clear();

      RankToBoundMap sendBounds, recvBounds;
      for (RankToBoundMap::iterator rankIt = periodic.begin();
	   rankIt != periodic.end(); ++rankIt) {

	if (rankIt->first == globalMpiRank)
	  continue;

	if (rankIt->first < globalMpiRank)
	  sendBounds[rankIt->first] = rankIt->second;
	else
	  recvBounds[rankIt->first] = rankIt->second;	
      }

      stdMpi.send(sendBounds);
      stdMpi.recv(recvBounds);
      stdMpi.startCommunication();

      for (RankToBoundMap::iterator rankIt = periodic.begin();
	   rankIt != periodic.end(); ++rankIt) {

 	if (rankIt->first <= globalMpiRank)
 	  continue;
  
	for (unsigned int j = 0; j < rankIt->second.size(); j++) {
	  BoundaryObject &recvRankObj = 
	    stdMpi.getRecvData()[rankIt->first][j].rankObj;
	  BoundaryObject &recvNeighObj = 
	    stdMpi.getRecvData()[rankIt->first][j].neighObj;

	  if (periodic[rankIt->first][j].neighObj != recvRankObj ||
	      periodic[rankIt->first][j].rankObj != recvNeighObj) {
	    unsigned int k = j + 1;	    
	    for (; k < rankIt->second.size(); k++)
	      if (periodic[rankIt->first][k].neighObj == recvRankObj &&
		  periodic[rankIt->first][k].rankObj == recvNeighObj)
		break;
	    
	    // The element must always be found, because the list is just in 
	    // another order.
	    TEST_EXIT_DBG(k < rankIt->second.size())("Should never happen!\n");

	    // Swap the current with the found element.
	    AtomicBoundary tmpBound = (rankIt->second)[k];
	    (rankIt->second)[k] = (rankIt->second)[j];
	    (rankIt->second)[j] = tmpBound;	
	  } 
	}
      }     
    } // periodicBoundary.boundary.size() > 0


    // === Run verification procedure. ===
    
    verifyBoundary();
  }


  int InteriorBoundary::getDegreeOwn(BoundaryObject &bObj)
  {
    int counter = 0;

    for (RankToBoundMap::iterator it = own.begin(); it != own.end(); ++it) {
      for (vector<AtomicBoundary>::iterator bIt = it->second.begin(); 
	   bIt != it->second.end(); ++bIt) {
	if (bIt->rankObj == bObj) {
	  counter++;
	  break;
	}
      }
    }

    return counter;
  }


  void InteriorBoundary::serialize(ostream &out)
  {
    serialize(out, own);
    serialize(out, other);
    serialize(out, periodic);
  }


  void InteriorBoundary::serialize(ostream &out,
				   RankToBoundMap& boundary)
  {
    int mSize = boundary.size();
    SerUtil::serialize(out, mSize);
    for (RankToBoundMap::iterator it = boundary.begin(); 
	 it != boundary.end(); ++it) {
      int rank = it->first;
      int boundSize = it->second.size();
      SerUtil::serialize(out, rank);
      SerUtil::serialize(out, boundSize);
      for (int i = 0; i < boundSize; i++) {
	AtomicBoundary &bound = (it->second)[i];

	SerUtil::serialize(out, bound.rankObj.elIndex);
	SerUtil::serialize(out, bound.rankObj.elType);
	SerUtil::serialize(out, bound.rankObj.subObj);
	SerUtil::serialize(out, bound.rankObj.ithObj);
	SerUtil::serialize(out, bound.rankObj.reverseMode);
	serializeExcludeList(out, bound.rankObj.excludedSubstructures);

	SerUtil::serialize(out, bound.neighObj.elIndex);
	SerUtil::serialize(out, bound.neighObj.elType);
	SerUtil::serialize(out, bound.neighObj.subObj);
	SerUtil::serialize(out, bound.neighObj.ithObj);
	SerUtil::serialize(out, bound.neighObj.reverseMode);
	serializeExcludeList(out, bound.neighObj.excludedSubstructures);

	SerUtil::serialize(out, bound.type);
      }
    }
  }


  void InteriorBoundary::deserialize(istream &in, Mesh *mesh)				     
  {
    map<int, Element*> elIndexMap;
    mesh->getElementIndexMap(elIndexMap);

    deserialize(in, own, elIndexMap);
    deserialize(in, other, elIndexMap);
    deserialize(in, periodic, elIndexMap);
  }


  void InteriorBoundary::deserialize(istream &in, 
				     RankToBoundMap& boundary,
				     map<int, Element*> &elIndexMap)
  {
    FUNCNAME_DBG("InteriorBoundary::deserialize()");

    int mSize = 0;
    SerUtil::deserialize(in, mSize);
    for (int i = 0; i < mSize; i++) {
      int rank = 0;
      int boundSize = 0;
      SerUtil::deserialize(in, rank);
      SerUtil::deserialize(in, boundSize);

      boundary[rank].resize(boundSize);
      for (int i = 0; i < boundSize; i++) {
	AtomicBoundary &bound = boundary[rank][i];

	SerUtil::deserialize(in, bound.rankObj.elIndex);
	SerUtil::deserialize(in, bound.rankObj.elType);
	SerUtil::deserialize(in, bound.rankObj.subObj);
	SerUtil::deserialize(in, bound.rankObj.ithObj);
	SerUtil::deserialize(in, bound.rankObj.reverseMode);
	deserializeExcludeList(in, bound.rankObj.excludedSubstructures);

	SerUtil::deserialize(in, bound.neighObj.elIndex);
	SerUtil::deserialize(in, bound.neighObj.elType);
	SerUtil::deserialize(in, bound.neighObj.subObj);
	SerUtil::deserialize(in, bound.neighObj.ithObj);
	SerUtil::deserialize(in, bound.neighObj.reverseMode);
	deserializeExcludeList(in, bound.neighObj.excludedSubstructures);

	SerUtil::deserialize(in, bound.type);

	TEST_EXIT_DBG(elIndexMap.count(bound.rankObj.elIndex) == 1)
	  ("Cannot find element with index %d for deserialization!\n", 
	   bound.rankObj.elIndex);

	TEST_EXIT_DBG(elIndexMap[bound.rankObj.elIndex]->getIndex() == 
		      bound.rankObj.elIndex)("Should not happen!\n");

	bound.rankObj.el = elIndexMap[bound.rankObj.elIndex];

	// For the case of periodic interior boundaries, a rank may have an
	// boundary with itself. In this case, also the pointer to the neighbour
	// object must be set correctly.
	if (elIndexMap.count(bound.neighObj.elIndex))
	  bound.neighObj.el = elIndexMap[bound.neighObj.elIndex];
	else
	  bound.neighObj.el = NULL;
      }
    }
  }


  void InteriorBoundary::verifyBoundary()
  {
    FUNCNAME("InteriorBoundary::verifyBoundary()");

    // === Check if no other boundery rank object is also included in own ===
    // === boundary part, which would make no sence.                      ===

    std::set<BoundaryObject> allOwnBounds;
    for (map<int, vector<AtomicBoundary> >::iterator it = own.begin(); 
	 it != own.end(); ++it)
      for (vector<AtomicBoundary>::iterator bIt = it->second.begin();
	   bIt != it->second.end(); ++bIt)
	allOwnBounds.insert(bIt->rankObj);

    for (map<int, vector<AtomicBoundary> >::iterator it = other.begin(); 
	 it != other.end(); ++it)
      for (vector<AtomicBoundary>::iterator bIt = it->second.begin();
	   bIt != it->second.end(); ++bIt)
	if (allOwnBounds.count(bIt->rankObj)) {
	  ERROR_EXIT("Boundary %d/%d/%d is included in both, own and other interior boundary part!\n",
		     bIt->rankObj.elIndex,
		     bIt->rankObj.subObj,
		     bIt->rankObj.ithObj);
	}
  }


  AtomicBoundary& InteriorBoundary::getNewOwn(int rank)
  {
    int size = own[rank].size();
    own[rank].resize(size + 1);
    return own[rank][size];
  }


  AtomicBoundary& InteriorBoundary::getNewOther(int rank)
  {
    int size = other[rank].size();
    other[rank].resize(size + 1);
    return other[rank][size];
  }


  AtomicBoundary& InteriorBoundary::getNewPeriodic(int rank)
  {
    FUNCNAME("InteriorBoundary::getNewPeriodic()");

    int size = periodic[rank].size();
    periodic[rank].resize(size + 1);
    return periodic[rank][size];
  }


  bool InteriorBoundary::checkOther(AtomicBoundary& bound, int rank)
  {
    FUNCNAME("InteriorBoundary::checkOther()");

    int globalMpiRank = MPI::COMM_WORLD.Get_rank();
    TEST_EXIT(rank > globalMpiRank)
      ("Wrong MPI ranks: %d %d\n", rank, globalMpiRank);

    bool found = false;

    if (other.count(rank)) {
      vector<AtomicBoundary>& otherBounds = other[rank];
      for (int i = 0; i < static_cast<int>(otherBounds.size()); i++) {
	if (otherBounds[i].rankObj == bound.rankObj) {
	  if (otherBounds[i].neighObj != bound.neighObj) {
	    ERROR_EXIT("If this really can happen, I have to thing about this!\n");
	  } else {
	    found = true;
	    break;
	  }
	}
      }
    }

    return found;
  }


  bool InteriorBoundary::removeOwn(BoundaryObject& bound)
  {
    bool removed = false;

    for (map<int, vector<AtomicBoundary> >::iterator it = own.begin();
	 it != own.end(); ++it) {
      for (vector<AtomicBoundary>::iterator bIt = it->second.begin();
	   bIt != it->second.end(); ++bIt) {
	if (bIt->rankObj == bound) {
	  it->second.erase(bIt);
	  removed = true;
	  break;
	}
      }
    }

    return removed;
  }


  bool InteriorBoundary::removeOther(BoundaryObject& bound)
  {
    bool removed = false;

    for (map<int, vector<AtomicBoundary> >::iterator it = other.begin();
	 it != other.end(); ++it) {
      for (vector<AtomicBoundary>::iterator bIt = it->second.begin();
	   bIt != it->second.end(); ++bIt) {
	if (bIt->rankObj == bound) {
	  it->second.erase(bIt);
	  removed = true;
	  break;
	}
      }
    }

    return removed;    
  }


  void InteriorBoundary::serializeExcludeList(ostream &out, 
					      ExcludeList &list)
  {
    int size = list.size();
    SerUtil::serialize(out, size);
    for (int i = 0; i < size; i++) {
      SerUtil::serialize(out, list[i].first);
      SerUtil::serialize(out, list[i].second);
    }
  }


  void InteriorBoundary::deserializeExcludeList(istream &in, 
						ExcludeList &list)
  {
    int size = 0;
    SerUtil::deserialize(in, size);
    list.resize(0);
    list.reserve(size);

    for (int i = 0; i < size; i++) {
      GeoIndex a;
      int b;

      SerUtil::deserialize(in, a);
      SerUtil::deserialize(in, b);
      list.push_back(make_pair(a, b));
    }
  }


  void MultiLevelInteriorBoundary::create(MeshLevelData &levelData,
					  ElementObjectDatabase &elObjDb)
  {
    levelIntBound.clear();

    int nLevel = levelData.getNumberOfLevels();
    for (int level = 0; level < nLevel; level++)
      levelIntBound[level].create(levelData, level, elObjDb);
  }


  void MultiLevelInteriorBoundary::serialize(ostream &out)
  {
    FUNCNAME("MultiLevelInteriorBoundary::serialize()");
    ERROR_EXIT("Not yet implemented!\n");
  }


  void MultiLevelInteriorBoundary::deserialize(istream &in, Mesh *mesh)
  {
    FUNCNAME("MultiLevelInteriorBoundary::deserialize()");
    ERROR_EXIT("Not yet implemented!\n");
  }

} }
