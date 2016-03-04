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



/** \file TranslateLsFct.h */

#ifndef AMDIS_TRANSLATELSFCT_H
#define AMDIS_TRANSLATELSFCT_H

#include "ElementFunction.h"
#include "FixVec.h"

namespace compositeFEM
{

  using namespace AMDiS;

  template<typename T>
  class TranslateLsFct : public ElementFunction<T>
  {
  public:
    /// Constructor
    TranslateLsFct(ElementFunction<T>* f_, double c_)
      : ElementFunction<T>(),
        f(f_),
        c(c_)
    {}

    /**
     * evaluation at given coordinates.
     *
     * Defines new level set function. Zero level set of this level set
     * function is level set c of level set function f.
     */
    const T& operator()(const DimVec<double>& bary) const
    {
      f->setElInfo(this->elInfo);
      return ((*f)(bary) - c);
    };

  protected:
    /// Level set function which is translated.
    ElementFunction<T>* f;

    /// Translation parameter.
    double c;
  };

}

using compositeFEM::TranslateLsFct;

#endif  // AMDIS_TRANSLATELSFCT_H
