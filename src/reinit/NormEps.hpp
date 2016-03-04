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




#ifndef NORMEPS_H
#define NORMEPS_H

#include "FixVec.h"
#include "Global.h"
#include "Initfile.h"

namespace reinit
{

  using namespace AMDiS;

  class NormEps
  {
  public:

    /// Calculates regularized norm of WorldVector x.
    static double calcNormEps(const WorldVector<double>& x);

    /// Calculates square of regularized norm of WorldVector x.
    static double calcNormEps2(const WorldVector<double>& x);

    /// Set regularization epsilon eps.
    static inline void setEps()
    {
      FUNCNAME("NormEps::setEps()");

      eps = 0.0;
      Parameters::get("NormEps->epsilon for norm regularization", eps);

      TEST_EXIT(eps > 1.e-15)("illegal eps for norm regularization !\n");
    }

  protected:
    /// Epsilon for regularization.
    static double eps;
  };

}

using reinit::NormEps;

#endif // NORMEPS_H
