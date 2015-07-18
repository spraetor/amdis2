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


#include "DualTraverse.h"
#include "Mesh.h"
#include "ElInfo.h"

namespace AMDiS {

  bool DualTraverse::traverseFirst(Mesh *mesh1, 
				   Mesh *mesh2, 
				   int level1, 
				   int level2, 
				   Flag flag1,
				   Flag flag2,
				   ElInfo **elInfo1,
				   ElInfo **elInfo2,
				   ElInfo **elInfoSmall,
				   ElInfo **elInfoLarge)
  {
    FUNCNAME("DualTraverse::traverseFirst()");
    // replace CALL_EL_LEVEL by CALL_MG_LEVEL (covers whole domain)
    if (flag1.isSet(Mesh::CALL_EL_LEVEL)) {
      flag1 &= ~Mesh::CALL_EL_LEVEL;
      flag1 |= Mesh::CALL_MG_LEVEL;
      level1_ = level1;
    } else {
      level1_ = -1;
    }

    if (flag2.isSet(Mesh::CALL_EL_LEVEL)) {
      flag2 &= ~Mesh::CALL_EL_LEVEL;
      flag2 |= Mesh::CALL_MG_LEVEL;
      level2_ = level2;
    } else {
      level2_ = -1;
    }

    // replace CALL_LEAF_EL_LEVEL by CALL_MG_LEVEL (covers whole domain)
    if (flag1.isSet(Mesh::CALL_LEAF_EL_LEVEL)) {
      flag1 &= ~Mesh::CALL_LEAF_EL_LEVEL;
      flag1 |= Mesh::CALL_MG_LEVEL;
      level1_ = level1;
      callLeafElLevel1_ = true;
    } else {
      level1_ = -1;
      callLeafElLevel1_ = false;
    }

    if (flag2.isSet(Mesh::CALL_LEAF_EL_LEVEL)) {
      flag2 &= ~Mesh::CALL_LEAF_EL_LEVEL;
      flag2 |= Mesh::CALL_MG_LEVEL;
      level2_ = level2;
      callLeafElLevel2_ = true;
    } else {
      level2_ = -1;
      callLeafElLevel2_ = false;
    }

    // call standard traverse
    *elInfo1 = stack1.traverseFirst(mesh1, level1, flag1);
    while (*elInfo1 != NULL && skipEl1(*elInfo1)) {
      *elInfo1 = stack1.traverseNext(*elInfo1);
    }

    *elInfo2 = stack2.traverseFirst(mesh2, level2, flag2);
    while (*elInfo2 != NULL && skipEl2(*elInfo2)) {
      *elInfo2 = stack2.traverseNext(*elInfo2);
    }
 
    // finished ?
    if (*elInfo1 == NULL || *elInfo2 == NULL) {
      TEST_EXIT(*elInfo1 == *elInfo2)("invalid dual traverse\n");
      return false;
    }

    rest = 1.0;

    bool accepted = check(elInfo1, elInfo2, elInfoSmall, elInfoLarge);

    // check for non domain covering level traverse
    if (!accepted ||
	(level1_ != -1 && (*elInfo1)->getLevel() != level1_) ||
	(callLeafElLevel1_ && !(*elInfo1)->getElement()->isLeaf()) ||
	(level2_ != -1 && (*elInfo2)->getLevel() != level2_) ||
	(callLeafElLevel2_ && !(*elInfo2)->getElement()->isLeaf())) {
      return traverseNext(elInfo1, elInfo2, elInfoSmall, elInfoLarge);
    }

    fillSubElInfo(*elInfo1, *elInfo2, *elInfoSmall, *elInfoLarge);

    return true;
  }


  bool DualTraverse::traverseNext(ElInfo **elInfo1,
				  ElInfo **elInfo2,
				  ElInfo **elInfoSmall,
				  ElInfo **elInfoLarge)
  {
    FUNCNAME("DualTraverse::traverseNext()");
    // call standard traverse
    if (inc1) {
      do {
	*elInfo1 = stack1.traverseNext(*elInfo1);
      } while(*elInfo1 != NULL && skipEl1(*elInfo1));
    }
    if (inc2) {
      do {
	*elInfo2 = stack2.traverseNext(*elInfo2);
      } while (*elInfo2 != NULL && skipEl2(*elInfo2));
    }

    // finished ?
    if (*elInfo1 == NULL || *elInfo2 == NULL) {
      TEST_EXIT(*elInfo1 == *elInfo2)("invalid dual traverse\n");
      return false;
    }

    // finished ?
    if (*elInfo1 == NULL || *elInfo2 == NULL) {
      TEST_EXIT(*elInfo1 == *elInfo2)("invalid dual traverse\n");
      return false;
    }

    bool accepted = check(elInfo1, elInfo2, elInfoSmall, elInfoLarge);

    // check for non domain covering level traverse
    if (!accepted ||
	(level1_ != -1 && (*elInfo1)->getLevel() != level1_) ||
	(callLeafElLevel1_ && !(*elInfo1)->getElement()->isLeaf()) ||
	(level2_ != -1 && (*elInfo2)->getLevel() != level2_) ||
	(callLeafElLevel2_ && !(*elInfo2)->getElement()->isLeaf())) {
      return traverseNext(elInfo1, elInfo2, elInfoSmall, elInfoLarge);
    }

    fillSubElInfo(*elInfo1, *elInfo2, *elInfoSmall, *elInfoLarge);

    return true;
  }


  void DualTraverse::prepareNextStep(ElInfo **elInfo1,
				     ElInfo **elInfo2,
				     ElInfo **elInfoSmall,
				     ElInfo **elInfoLarge)
  {
    // which is the smaller element ?
    *elInfoSmall = 
      (*elInfo1)->getLevel() > (*elInfo2)->getLevel() ?
      *elInfo1 :
      *elInfo2;
    *elInfoLarge = 
      (*elInfo1)->getLevel() <= (*elInfo2)->getLevel() ?
      *elInfo1 :
      *elInfo2;

    // update rest
    rest -= 1.0 / (1 << ((*elInfoSmall)->getLevel() - (*elInfoLarge)->getLevel()));

    if (rest < 1e-32) {
      // large element finished -> increment both elements
      rest = 1.0; 
      inc1 = true;
      inc2 = true;
    } else {
      // increment only small element
      inc1 = (*elInfo1 == *elInfoSmall) ? true : false;
      inc2 = (*elInfo2 == *elInfoSmall) ? true : false;
    }
  }


  void DualTraverse::fillSubElInfo(ElInfo *elInfo1, 
				   ElInfo *elInfo2,
				   ElInfo *elInfoSmall,
				   ElInfo *elInfoLarge)
  {
    if (!fillSubElemMat)
      return;

    if (elInfo1 == elInfoSmall)
      stack1.fillRefinementPath(*elInfoSmall, *elInfo2);
    else
      stack2.fillRefinementPath(*elInfoSmall, *elInfo1);    
  }


  int DualTraverse::getFace(DualElInfo *dualElInfo, int largeFace)
  {
    FUNCNAME_DBG("DualTraverse::getFace()");

    TEST_EXIT_DBG(dualElInfo)("No dual element info object!\n");

    ElInfo *largeElInfo = dualElInfo->largeElInfo;
    ElInfo *smallElInfo = dualElInfo->smallElInfo;

    TEST_EXIT_DBG(largeElInfo->getLevel() <= smallElInfo->getLevel())
      ("Should not happen!\n");

    if (largeElInfo->getLevel() == smallElInfo->getLevel())
      return largeFace;

    TEST_EXIT_DBG(smallElInfo->getRefinementPathLength() != 0)
      ("Refinement path is empty. This should not happen!\n");

    return -1;
  }
}
