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


#include "BoundaryElementNormalDist.h"

namespace reinit
{

  using namespace AMDiS;
  int
  BoundaryElementNormalDist::calcDistOnBoundaryElement(
    ElInfo* elInfo,
    FixVec<double, VERTEX>& dVec)
  {

    // Get intersection information.
    int  elStatus = elLS->createElementLevelSet(elInfo);
    if (elStatus != ElementLevelSet::LEVEL_SET_BOUNDARY)
      return elStatus;

    VectorOfFixVecs<DimVec<double>>* elIntersecPoints =
                                   elLS->getElIntersecPoints();

    // Calculate (unit) normal to intersection plane.
    FixVec<WorldVector<double>, DIMEN> planeVecs(dim);
    WorldVector<double> normalVec;
    for (int i=0; i<dim; ++i)
    {
      elInfo->coordToWorld((*elIntersecPoints)[i], (planeVecs[i]));
    }
    calcNormal(planeVecs, normalVec);

    // Calculate normal distance for all vertices.
    for (int i=0; i<=dim; ++i)
    {

      dVec[i] = calculateDistLevelSetWithNormal(elInfo,
                i,
                normalVec,
                planeVecs[0]);
    }

    return elStatus;
  }

  double
  BoundaryElementNormalDist::calculateDistLevelSetWithNormal(
    ElInfo* elInfo,
    int vert,
    WorldVector<double>& normalVec,
    WorldVector<double>& planeVec)
  {
    // Get world coordinates of vertex.
    const WorldVector<double>& vertex = elInfo->getCoord(vert);

    // Calculate distance.
    double dist = 0.0;
    for (int i=0; i<dim; ++i)
    {
      dist += (vertex[i]-planeVec[i]) * normalVec[i];
    }

    return fabs(dist);
  }

}
