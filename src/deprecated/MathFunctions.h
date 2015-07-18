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



/** \file MathFunctions.h */

#ifndef AMDIS_MATHFUNCTIONS_H
#define AMDIS_MATHFUNCTIONS_H

#include "Global.h"
#include <limits>

namespace AMDiS 
{
  //converts signed distance to phasefield
  inline double Phi1(double r, double eps) { return 0.5 * (1 - tanh(3 * r / eps)); }
  inline double Phi2(double r, double eps) { return 0.5 * (1 + tanh(3 * r / eps)); }
  
  //levelset: positive (1) in the set, negative (-1) outside, zero on the boundary
  inline double LevelSet(double r)
  {
    if (r < 0)
      return 1;
    if (r > 0)
      return -1;
    return 0;
  }

  inline double Phi1ToR(double p1, double eps) 
  {
    double x = std::max(-1.0 + std::numeric_limits<double>::epsilon(), 
		    std::min(1.0 - std::numeric_limits<double>::epsilon(), p1));
    return eps / 3.0 * log((1 + x) / (1 - x)) * 0.5;
  }

  inline double Phi2ToR(double p2, double eps) 
  {
    double x = std::max(-1.0 + std::numeric_limits<double>::epsilon(), 
		    std::min(1.0 - std::numeric_limits<double>::epsilon(), 1 + 2 * p2));
    return eps / 3.0 * log( (1 + x) / (1 - x) );
  }

}
#endif
