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



#ifndef BOUNDARYELEMENTNORMALDIST_H
#define BOUNDARYELEMENTNORMALDIST_H

#include "ElInfo.h"
#include "FixVec.h"
#include "ElementLevelSet.h"
#include "BoundaryElementDist.h"

namespace reinit {
  
  using namespace AMDiS;

class BoundaryElementNormalDist : public BoundaryElementDist
{
public:

  BoundaryElementNormalDist(ElementLevelSet *elLS_, int dim_)
    : BoundaryElementDist(elLS_, dim_)
  {}

  /**
   * Calculates distance from the interface for all vertices of a boundary 
   * element.
   * Distance is here the normal distance.
   *
   * Return value: Status of element elInfo.
   */
  int calcDistOnBoundaryElement(ElInfo *elInfo,
				FixVec<double, VERTEX> &dVec);

 protected:
  /**
   * Calculates the distance of vertex vert in element elInfo
   * from the plane given by the normal vector of the plane normalVec
   * and one point lying in the plane planeVec.
   */
  double calculateDistLevelSetWithNormal(ElInfo *elInfo, 
					 int vert, 
					 WorldVector<double> &normalVec,
					 WorldVector<double> &planeVec);
};

}

using reinit::BoundaryElementNormalDist;

#endif  // BOUNDARYELEMENTNORMALDIST_H
