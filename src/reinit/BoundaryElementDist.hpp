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



#ifndef BOUNDARYELEMENTDIST_H
#define BOUNDARYELEMENTDIST_H

#include "ElInfo.h"
#include "FixVec.h"

#include "ElementLevelSet.h"

namespace reinit
{

  using namespace AMDiS;

  class BoundaryElementDist
  {
  public:

    BoundaryElementDist(ElementLevelSet* elLS_, int dim_)
      : dim(dim_),
        elLS(elLS_)
    {
      FUNCNAME("BoundaryElementDist::BoundaryElementDist()");

      TEST_EXIT(dim == 2 || dim == 3)
      ("function only works for dimension 2 !\n");
    }

    /**
     * Virtual destructor.
     */
    virtual ~BoundaryElementDist() {};

    /**
     * Calculates distance from the interface for all vertices of a boundary
     * element.
     *
     * Pure virtual: Realizations in
     *               BoundaryElementLevelSetDist - boundary value initialization
     *                                             by level set function
     *               BoundaryElementTopDist      - topological distance
     *                                             (distance within element)
     *               BoundaryElementEdgeDist     - distance along edges
     *               BoundaryElementNormalDist   - normal distance to
     *                                             intersection plane
     */
    virtual int calcDistOnBoundaryElement(ElInfo* elInfo,
                                          FixVec<double, VERTEX>& dVec) = 0;

  protected:
    /**
     * Calculate the (unit) normal for the plane defined by the world vectors in
     * vecs.
     */
    void calcNormal(const FixVec<WorldVector<double>, DIMEN>& vecs,
                    WorldVector<double>& normalVec);

    /**
     * Normal calculation: 2-dimensional case.
     */
    void calcNormal_2d(const FixVec<WorldVector<double>, DIMEN>& vecs,
                       WorldVector<double>& normalVec);

    /**
     * Normal calculation: 3-dimensional case.
     */
    void calcNormal_3d(const FixVec<WorldVector<double>, DIMEN>& vecs,
                       WorldVector<double>& normalVec);

  protected:
    /**
     * Dimension.
     */
    int dim;

    /**
     * Holds level set function and functionalities for intersection point
     * calculation.
     */
    ElementLevelSet* elLS;
  };

}

using reinit::BoundaryElementDist;

#endif  // BOUNDARYELEMENTDIST_H
