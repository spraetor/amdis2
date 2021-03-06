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




#ifndef BOUNDARYELEMENTLEVELSETDIST_H
#define BOUNDARYELEMENTLEVELSETDIST_H

#include "ElInfo.h"
#include "FixVec.h"
#include "ElementLevelSet.h"
#include "BoundaryElementDist.h"

namespace reinit
{

  using namespace AMDiS;

  class BoundaryElementLevelSetDist : public BoundaryElementDist
  {
  public:

    BoundaryElementLevelSetDist(ElementLevelSet* elLS_, int dim_)
      : BoundaryElementDist(elLS_, dim_)
    {}

    /**
     * Calculates distance from the interface for all vertices of a boundary
     * element.
     * Distance is here the value of the level set function at boundary
     * vertices..
     *
     * Return value: Status of element elInfo.
     */
    int calcDistOnBoundaryElement(ElInfo* elInfo,
                                  FixVec<double, VERTEX>& dVec);
  };

}

using reinit::BoundaryElementLevelSetDist;

#endif  // BOUNDARYELEMENTLEVELSETDIST_H
