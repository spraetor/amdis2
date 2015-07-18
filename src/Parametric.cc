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


#include "Parametric.h"
#include "ElInfo.h"
#include "DOFVector.h"
#include "FixVec.h"

namespace AMDiS {

  ElInfo *ParametricFirstOrder::addParametricInfo(ElInfo *elInfo) 
  {
    elInfo->setParametric(true);
    int dow = Global::getGeo(WORLD);
    Element *element = elInfo->getElement();
    const DegreeOfFreedom **dof = element->getDof();
  
    for (int i = 0; i < elInfo->getElement()->getGeo(VERTEX); i++) {
      if (elInfo->getFillFlag().isSet(Mesh::FILL_COORDS))
	for (int j = 0; j < dow; j++)
	  elInfo->getCoord(i)[j] = (*(*dofCoords_)[j])[dof[i][0]];

      if (elInfo->getFillFlag().isSet(Mesh::FILL_OPP_COORDS)) {
	TEST_EXIT(elInfo->getFillFlag().isSet(Mesh::FILL_NEIGH))
	  ("FILL_NEIGH not set\n");

	if (elInfo->getNeighbour(i)) {
	  const DegreeOfFreedom **neighDof = elInfo->getNeighbour(i)->getDof();
	  for (int j = 0; j < dow; j++)
	    elInfo->getOppCoord(i)[j] = 
	      (*(*dofCoords_)[j])[neighDof[elInfo->getOppVertex(i)][0]];	  
	}
      }
    }

    return elInfo;
  }

  ElInfo *ParametricFirstOrder::removeParametricInfo(ElInfo *elInfo) 
  {
    elInfo->setParametric(false);
    return elInfo;
  }

}
