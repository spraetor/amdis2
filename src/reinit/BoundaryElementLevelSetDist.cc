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


#include "BoundaryElementLevelSetDist.h"

namespace reinit {

using namespace AMDiS;
int 
BoundaryElementLevelSetDist::calcDistOnBoundaryElement(
                             ElInfo *elInfo,
			     FixVec<double, VERTEX> &dVec)
{
  // Get intersection information.
  int  elStatus = elLS->createElementLevelSet(elInfo);
  if (elStatus != ElementLevelSet::LEVEL_SET_BOUNDARY)
    return elStatus;
  
  const double  *elVertLevelSetVec = elLS->getElVertLevelSetVec();
  
  // Set distance to values of level set function in element vertices.
  for (int i=0; i<=dim; ++i) {
    dVec[i] = fabs(elVertLevelSetVec[i]);
  }
  
  return elStatus;
}

}
