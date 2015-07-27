/** \file SurfaceQuadrature.h */

#pragma once

#include "FixVec.h"
#include "Quadrature.h"

namespace AMDiS 
{
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

    /// Adapts SurfaceQuadrature to \ref coords.
    void scaleSurfaceQuadrature(VectorOfFixVecs<DimVec<double> > &coords);

  protected:
    /// Pointer to the original quadrature
    Quadrature *quad;

    VectorOfFixVecs<DimVec<double> > coords;
  };

} // end namespace AMDiS
