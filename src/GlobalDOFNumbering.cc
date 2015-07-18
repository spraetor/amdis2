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


#include "GlobalDOFNumbering.h"
#include "MeshStructure.h"

namespace AMDiS {

  GlobalDOFNumbering::GlobalDOFNumbering(std::vector<MeshStructure*> &meshCodes,
					 std::vector< std::vector<DegreeOfFreedom> > &dofCodes,
					 int dofsPerElement)
    : dofsPerElement_(dofsPerElement),
      globalIndexCounter_(0)
  {
    int rank, mpiSize = static_cast<int>(meshCodes.size());

    std::vector< std::vector<DegreeOfFreedom>::iterator> it(mpiSize);

    // init
    for(rank = 0; rank < mpiSize; rank++) {
      meshCodes[rank]->reset();
      it[rank] = dofCodes[rank].begin();
    }

    // start recursion
    while(meshCodes[0]->getCurrentElement() < meshCodes[0]->getNumElements()) {
      // one call for each macro element
      handleRecursive(meshCodes, it);
    }
  }

  void GlobalDOFNumbering::handleRecursive(std::vector<MeshStructure*> &meshCodes,
					   std::vector< std::vector<DegreeOfFreedom>::iterator> &it)
  {
    int rank, mpiSize = static_cast<int>(meshCodes.size());
    int dofNr;

    // === handle ===
    int firstRankUsed = -1;
    for(rank = 0; rank < mpiSize; rank++) {
      if(meshCodes[rank] != NULL) {
	firstRankUsed = rank;
	break;
      }
    }

    TEST_EXIT(firstRankUsed > -1)("no used rank in recursion\n");

    for(dofNr = 0; dofNr < dofsPerElement_; dofNr++) {
      bool newDOF = (localToGlobal_[firstRankUsed][*(it[firstRankUsed])] == 0);

      DegreeOfFreedom globalIndex = -1;
      if(newDOF) {
	globalIndex = globalIndexCounter_++;
      }

      for(rank = 0; rank < mpiSize; rank++) {
	if(meshCodes[rank]) {
	  if(newDOF) {
	    DegreeOfFreedom localIndex = *(it[rank]);
	    localToGlobal_[rank][localIndex] = globalIndex + 1;
	    globalToLocal_[rank][globalIndex] = localIndex + 1;
	  }
	  ++(it[rank]);
	}
      }
    }

    // === recursion ===
    std::vector<MeshStructure*> recursiveCodes;
    bool cont = false;
    for(rank = 0; rank < mpiSize; rank++) {
      if(meshCodes[rank]) {
	if(!meshCodes[rank]->isLeafElement()) {
	  recursiveCodes[rank] = meshCodes[rank];
	  cont = true;
	} else {
	  recursiveCodes[rank] = NULL;
	}
	meshCodes[rank]->nextElement();
      } else {
	recursiveCodes[rank] = NULL;
      }
    }

    if(cont) {
      handleRecursive(recursiveCodes, it); // left tree
      handleRecursive(recursiveCodes, it); // right tree
    }
  }

  DegreeOfFreedom GlobalDOFNumbering::getLocalIndex(int rank, DegreeOfFreedom globalIndex)
  {
    return (globalToLocal_[rank][globalIndex] - 1);
  }

  DegreeOfFreedom GlobalDOFNumbering::getGlobalIndex(int rank, DegreeOfFreedom localIndex)
  {
    return (localToGlobal_[rank][localIndex] - 1);
  }

}
