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



/** \file SubElInfo.h */

#ifndef AMDIS_SUBELINFO_H
#define AMDIS_SUBELINFO_H

#include "AMDiS_fwd.h"
#include "FixVec.h"

namespace AMDiS {
  
  /** \brief
   * The class \ref SubElInfo holds all information on a subelement (element to
   * which it belongs, its vertices in barycentric coordinates with respect to 
   * the element, corresponding determinant, ...).
   *
   * Main routines:
   * SubElInfo() - Creates a \ref SubElInfo for a subelement of the element 
   *               \ref elInfo. The vertices of the subelement given in 
   *               barycentric coordinates with respect to element are handed over 
   *               to this routine. The determinant corresponding to subelement
   *               is calculated.
   */
  class SubElInfo
  {
  public:
    /// Constructor
    SubElInfo(VectorOfFixVecs<DimVec<double> > *lambda, const ElInfo *elInfo); 

    /// Destructor
    ~SubElInfo() 
    {
      if (lambda)
	delete lambda;
    }

    /// Get b-th coordinate of the a-th vertex of subelement (barycentriccoordinates). 
    double getLambda(int a, int b) const 
    {
      if (lambda) 
	return (*lambda)[a][b]; 
      else 
	return 0.0;
    }

    /// Get coordinates of a-th vertex of subelement (barycentric coordinates).
    const DimVec<double>& getLambda(int a) const 
    { 
      return (*lambda)[a]; 
    }

    /// Get coordinates of all vertices of subelement (barycentric coordinates).
    VectorOfFixVecs<DimVec<double> > *getLambda() const 
    { 
      return lambda; 
    }

    /// Get determinant corresponding to subelement.
    double getDet() const 
    {
      return det;
    }

  protected:

    /// Contains elInfo of the element that contains subelement.
    const ElInfo *elInfo;

    /// Barycentrc coordinates of the vertices of subelement.
    VectorOfFixVecs<DimVec<double> > *lambda;

    /// Determinant corresponding to the (World-)subelement.
    double det;
  };

}
#endif  // AMDIS_SUBELINFO_H
