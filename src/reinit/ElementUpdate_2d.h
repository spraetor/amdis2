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




#ifndef ELEMENTUPDATE_2D_H
#define ELEMENTUPDATE_2D_H

#include "FixVec.h"
#include "ElementUpdate.h"
#include "VelocityExt.h"

namespace reinit
{

  using namespace AMDiS;

  class ElementUpdate_2d : public ElementUpdate
  {
  public:
    ElementUpdate_2d(VelocityExt* velExt_ = NULL)
      : ElementUpdate(velExt_)
    {}

    /** \brief
     * Realization of ElementUpdate::calcElementUpdate. Calculates the Bornemann
     * update for an element.
     */
    double calcElementUpdate(const FixVec<WorldVector<double> *, VERTEX>& vert,
                             FixVec<double, VERTEX>& uhVal);

  private:
    WorldVector<double> xhminusYh, zhminusYh, xhminusZh;
  };

}

using reinit::ElementUpdate_2d;

#endif // ELEMENTUPDATE_2D_H

