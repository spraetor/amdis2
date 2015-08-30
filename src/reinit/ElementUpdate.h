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




#ifndef ELEMENTUPDATE_H
#define ELEMENTUPDATE_H

#include "FixVec.h"
#include "VelocityExt.h"

namespace reinit
{

  using namespace AMDiS;

  class ElementUpdate
  {
  public:
    ElementUpdate(VelocityExt* velExt_)
      : velExt(velExt_)
    {}

    virtual ~ElementUpdate() {}

    /**
     * Pure virtual function.
     * Calculate Bornemann update on an element. Realization is implemented for
     * dimensions 2 and 3.
     */
    virtual double calcElementUpdate(const FixVec<WorldVector<double> *, VERTEX>& vert,
                                     FixVec<double, VERTEX>& uhVal) = 0;

  protected:
    /**
     * Object needed to extrapolate velocity from the interface.
     */
    VelocityExt* velExt;
  };

}

using reinit::ElementUpdate;

#endif // ELEMENTUPDATE_H
