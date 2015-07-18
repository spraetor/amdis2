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


#include "BoundaryElementEdgeDist.h"

namespace reinit {

using namespace AMDiS;
int 
BoundaryElementEdgeDist::calcDistOnBoundaryElement(
			 ElInfo *elInfo,
			 FixVec<double, VERTEX> &dVec)
{
  WorldVector<double> vert1, vert2;
  int vert1Ind, vert2Ind;
  double length;
  double val;
  
  
  // Get intersection information.
  int  elStatus = elLS->createElementLevelSet(elInfo);
  if (elStatus != ElementLevelSet::LEVEL_SET_BOUNDARY)
    return elStatus;
  
  VectorOfFixVecs<DimVec<double> > *elIntersecPoints = 
    elLS->getElIntersecPoints();
  int numIntersecPoints = elLS->getNumElIntersecPoints();
  
  // For all vertices: calculate distance to intersection points.
  for (int i = 0; i < numIntersecPoints; ++i) {
    vert1Ind = -1; 
    vert2Ind = -1;
    length = 0;
    
    // Get vertices on edge with intersection point.
    for (int j=0; j<=dim; ++j) {
      if ((*elIntersecPoints)[i][j] > 1.e-15) {
	if (vert1Ind == -1) vert1Ind = j;
	else {
	  vert2Ind = j;
	  break;
	}
      }
    }
    
    // Get length of edge with intersection point.
    vert1 = elInfo->getCoord(vert1Ind);
    vert2 = elInfo->getCoord(vert2Ind);
    for (int j=0; j<dim; ++j) {
      length += (vert1[j] - vert2[j])*(vert1[j]-vert2[j]);
    }
    length = sqrt(length);
    
    // Get distance of vert1 and vert2 to intersection point.
    val = (*elIntersecPoints)[i][vert2Ind] * length;
    dVec[vert1Ind] = val;
    dVec[vert2Ind] = length - val;
  }
  
  return elStatus;
}

}
