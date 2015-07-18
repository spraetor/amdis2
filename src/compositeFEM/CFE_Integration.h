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



/** \file CFE_Integration.h */

#ifndef AMDIS_CFE_INTEGRATION_H
#define AMDIS_CFE_INTEGRATION_H

#include "ElementFunction.h"
#include "Quadrature.h"
#include "ElementLevelSet.h"


namespace compositeFEM {
  
  using namespace AMDiS;

  class CFE_Integration
  {
  public:
    /// Calculates integral of function f on domain where level set function is negative.
    static double integrate_onNegLs(ElementFunction<double> *f, 
				    ElementLevelSet *elLS,
				    int deg = 1, 
				    Quadrature *q = NULL);

    /**
     * Calculates surface integral of function f on the zero level set.
     *
     * Note: Quadrature q is a quadrature formula for dimension dim-1.
     */
    static double integrate_onZeroLs(ElementFunction<double> *f, 
				     ElementLevelSet *elLS,
				     int deg = 1, 
				     Quadrature *q = NULL);
  protected:
    /// Calculates determinant for surface given through surfVert.
    static double calcSurfaceDet(ElInfo *loc_elInfo,
				 VectorOfFixVecs<DimVec<double> > &surfVert);
  };

}

#endif  // AMDIS_CFE_INTEGRATION_H
