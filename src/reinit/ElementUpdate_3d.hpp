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




#ifndef ELEMENTUPDATE_3D_H
#define ELEMENTUPDATE_3D_H

#include "FixVec.h"
#include "ElementUpdate.h"
#include "ElementUpdate_2d.h"
#include "VelocityExt.h"

namespace reinit
{

  using namespace AMDiS;

  class ElementUpdate_3d : public ElementUpdate
  {
  public:
    ElementUpdate_3d(VelocityExt* velExt_ = NULL)
      : ElementUpdate(velExt_)
    {
      elUpdate2d = new ElementUpdate_2d(velExt_);
    }

    ~ElementUpdate_3d()
    {
      delete elUpdate2d;
    }

    /**
     * Realization of ElementUpdate::calcElementUpdate.
     * Calculates the Bornemann update for an element. The update is
     * calculated for the last vertex of \ref vert, denoted by xh.
     */
    double calcElementUpdate(const FixVec<WorldVector<double> *, VERTEX>& vert_,
                             FixVec<double, VERTEX>& uhVal);

  protected:
    /**
     * Sorts vertices of the face opposite to xh and translates the element
     * to the point of origin.
     *
     * Denote A the vertex \ref vert[0], B the vertex \ref vert[1], ...
     * after the sort, and uhVal(A), uhVal(B), ... the corresponding function
     * values.
     * Then, after the sort, uhVal(A) <= uhVal(B) <= uhVal(C).
     * Finally, translation of element moves A to zero, uhVal(A) is set to
     * zero, uhVal(B), uhVal(C) correspondingly.
     */
    void sortAndTranslateElement(FixVec<WorldVector<double> *, VERTEX>& vert,
                                 FixVec<double, VERTEX>& uhVal,
                                 double& uhVal0_orig);

    /**
     * Calculates rotation of element.
     *
     * Denote A1, B1, ... the vertices of element after rotation.
     * Then distance of B1 to x-y-plane is uhVal(B) and
     * distance of C1 to x-y-plane is uhVal(C).
     *
     * Return value: true - rotation exists
     *               false - rotation does not exist
     *
     * Assumptions: --> A == zero vector (then AB == B, ...)
     *              --> uhVal(A) == 0
     */
    bool rotateElement(const FixVec<WorldVector<double> *, VERTEX>& vert,
                       const FixVec<double, VERTEX>& uhVal,
                       WorldVector<double>& B1,
                       WorldVector<double>& C1,
                       WorldVector<double>& D1);

    /**
     * Determines position of D2 with respect to the triangle AB2C2.
     * A is the point of origin.
     *
     * Note: Limited to triangles with angles not greater than Pi/2
     *       (acute triangles).
     */
    int calcPosition(const WorldVector<double>& B2,
                     const WorldVector<double>& C2,
                     const WorldVector<double>& D2);

    /**
     * Calculate 2-dimensional Bornemann update for face of element. The
     * update is calculated for vertex with coordinates \ref *C2d.
     */
    double calcFaceUpdate(WorldVector<double>* A2d,
                          WorldVector<double>* B2d,
                          WorldVector<double>* C2d,
                          double& uhValA2d,
                          double& uhValB2d);

  protected:
    /**
     * Used for Bornemann update on faces of element.
     */
    ElementUpdate_2d* elUpdate2d;

    /**
     * Flags to decribe position of D2 with respect to AB2C2.
     */
    static const int IN_ABC   = 0;
    static const int VERT_A   = 1;
    static const int VERT_B   = 2;
    static const int VERT_C   = 3;
    static const int EDGE_AB  = 4;
    static const int EDGE_AC  = 5;
    static const int EDGE_BC  = 6;
  };

}

using reinit::ElementUpdate_3d;

#endif // ELEMENTUPDATE_3D_H
