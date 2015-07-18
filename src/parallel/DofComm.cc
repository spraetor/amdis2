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


#include "parallel/DofComm.h"
#include "parallel/InteriorBoundary.h"
#include "parallel/MeshLevelData.h"
#include "FiniteElemSpace.h"
#include "Debug.h"
#include "ElementDofIterator.h"
#include "DOFVector.h"

namespace AMDiS { namespace Parallel {

  using namespace std;


  void DofComm::init(vector<const FiniteElemSpace*> &fe)
  {
    feSpaces = fe;
    mesh = feSpaces[0]->getMesh();
    
    sendDofs.clear();
    recvDofs.clear();
    periodicDofs.clear();   
  }

  void DofComm::create(Mesh* mesh, InteriorBoundary &boundary)
  {
    FUNCNAME("DofComm::create()");
    TEST_EXIT(mesh == this->mesh) ("No equal mesh.\n");

    sendDofs.clear();
    recvDofs.clear();
    periodicDofs.clear();  
    
    createContainer(mesh, boundary, boundary.getOwn(), sendDofs);
    createContainer(mesh, boundary, boundary.getOther(), recvDofs);

#if (DEBUG != 0)
    {
      std::set<DegreeOfFreedom> sds;
      for (DofComm::Iterator it(sendDofs, feSpaces[0]);
         !it.end(); it.nextRank())
      for (; !it.endDofIter(); it.nextDof())
         sds.insert(it.getDofIndex());

      for (DofComm::Iterator it2(recvDofs, feSpaces[0]);
         !it2.end(); it2.nextRank())
        for (; !it2.endDofIter(); it2.nextDof())
          TEST_EXIT(!sds.count(it2.getDofIndex()))("Send and recv have same dofs.\n");
    }
#endif
  }



  void DofComm::createContainer(Mesh* mesh,
				InteriorBoundary &boundary,
                                RankToBoundMap &rankToBoundMap,
				DataType &data)
  {
    // === Fill data. ===

    for (unsigned int i = 0; i < feSpaces.size(); i++)
      for (InteriorBoundary::iterator it(rankToBoundMap); !it.end(); ++it) 
	boundary.getElementPtr(it->rankObj.elIndex, mesh)
	  ->getAllDofs(feSpaces[i], it->rankObj, data[it.getRank()][feSpaces[i]]);

    // === Remove empty data containers. ===
   
    DataIter dit = data.begin();
    while (dit != data.end()) {
      FeMapIter it = dit->second.begin();
      while (it != dit->second.end()) {
	if (it->second.size() == 0) {
	  const FiniteElemSpace *fe = it->first;
	  ++it;
	  dit->second.erase(fe);
	} else
	  ++it;
      }
      
      if (dit->second.size() == 0)
	data.erase(dit++);
      else
	++dit;
    }
  }


  void DofComm::serialize(ostream &out)
  {
    FUNCNAME("DofComm:serialize()");

    ERROR_EXIT("MUSS DAS WIRKLICH SEIN????\n");
  }

  
  int DofComm::getNumberDofs(DataType &data, 
			     const FiniteElemSpace *feSpace,
			     bool countDouble)
  {
    DofContainerSet dofSet;
    DofContainer dofVec;

    for (DataIter rankIt = data.begin(); rankIt != data.end(); ++rankIt) {
      for (FeMapIter feIt = rankIt->second.begin();
	   feIt != rankIt->second.end(); ++feIt) {
	if (feIt->first == feSpace) {
	  if (countDouble) {
	    dofVec.insert(dofVec.end(), feIt->second.begin(), feIt->second.end());
	  } else {
	    dofSet.insert(feIt->second.begin(), feIt->second.end());
	  }
	}
      }
    }

    if (countDouble)
      return static_cast<int>(dofVec.size());    
    return static_cast<int>(dofSet.size());
  }
  

  int DofComm::getDegree(const FiniteElemSpace *feSpace, 
			 const DegreeOfFreedom *dof)
  {
    int degree = 0;

    for (map<int, FeMapType>::iterator it = sendDofs.begin();
	 it != sendDofs.end(); ++it) {
      DofContainer &dc = it->second[feSpace];
      if (find(dc.begin(), dc.end(), dof) != dc.end())
	degree++;
    }

    for (map<int, FeMapType>::iterator it = recvDofs.begin();
	 it != recvDofs.end(); ++it) {
      DofContainer &dc = it->second[feSpace];
      if (find(dc.begin(), dc.end(), dof) != dc.end())
	degree++;
    }

    return degree;
  }


  bool DofComm::Iterator::setNextFeMap()
  {
    FUNCNAME_DBG("DofComm::Iterator::setNextFeMap()");

    if (dataIter != data.end()) {
      TEST_EXIT_DBG(dataIter->second.size())("Should not happen!\n");

      feMapIter = dataIter->second.begin();
      
      if (traverseFeSpace != NULL) {
	if ((dataIter->second.count(traverseFeSpace) == 0))
	  return false;
	
	while (feMapIter->first != traverseFeSpace &&
	       feMapIter != dataIter->second.end())
	  ++feMapIter;
	
	TEST_EXIT_DBG(feMapIter != dataIter->second.end() &&
		      feMapIter->first == traverseFeSpace)
	  ("Should not happen!\n");
      }
      
      dofIter = feMapIter->second.begin();      
      dofCounter = 0;
    }

    return true;
  }


  void MultiLevelDofComm::init(MeshLevelData &levelData,
			       vector<const FiniteElemSpace*> &fe)
  {
    int nLevel = levelData.getNumberOfLevels();
    for (int level = 0; level < nLevel; level++)
      levelDofComm[level].init(fe);
  }

  void MultiLevelDofComm::create(Mesh* mesh, MultiLevelInteriorBoundary &boundary)
  {
    for (map<int, DofComm>::iterator it = levelDofComm.begin();
	 it != levelDofComm.end(); ++it)
      it->second.create(mesh, boundary[it->first]);
  }

} }
