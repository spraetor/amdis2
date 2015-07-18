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



/** \file SurfaceQuadrature.h */

#ifndef AMDIS_SURFACEQAUDRATURE_H
#define AMDIS_SURFACEQAUDRATURE_H

#include "FixVec.h"
#include "Quadrature.h"

namespace AMDiS {

  /** 
   * \ingroup Integration
   *
   * \brief
   * Quadrature for element surfaces. Uses the functionality of the standard
   * Quadrature class for dim-1 but calculates new barycentric coordinates
   * for the element sides.
   */
  class SurfaceQuadrature : public Quadrature
  {
  public:
    /// Constructs a SurfaceQuadrature based on a standard Quadrature of dim-1.
    SurfaceQuadrature(Quadrature *quad, VectorOfFixVecs<DimVec<double> >& coords);

    /// Destructor.
    ~SurfaceQuadrature()
    {}

    /// Adapts SurfaceQuadrature to \ref coords.
    void scaleSurfaceQuadrature(VectorOfFixVecs<DimVec<double> > &coords);

  protected:
    /// Pointer to the original quadrature
    Quadrature *quad;

    VectorOfFixVecs<DimVec<double> > coords;
  };

}

#endif // AMDIS_SURFACEQAUDRATURE_H
