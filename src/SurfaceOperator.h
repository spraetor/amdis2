/** \file SurfaceOperator.h */

#pragma once

#include "AMDiS_fwd.h"
#include "Operator.h"
#include "FixVec.h"

namespace AMDiS
{
  /**
   * \ingroup Integration
   *
   * \brief
   * A SurfaceOperator is a Operator used for surface integration. Instead of a
   * normal Quadrature it uses a SurfaceQuadrature with integration points at
   * one element side instead of integration points in the inner of the element.
   * The determinant in ElInfo is replaced by the surface determinant, so the
   * assemblage can be done with the standard Assembler classes.
   * The SurfaceQuadrature is used for the implementation of Neumann and Robin
   * boundary conditions.
   */
  class SurfaceOperator : public Operator
  {
  public:
    /// Creates a SurfaceOperator conforming to operat for the given \ref coords.
    SurfaceOperator(Operator const& operat,
                    VectorOfFixVecs<DimVec<double>>& coords);

    /// Adapt surface quadratures to \ref coords.
    void adaptSurfaceOperator(VectorOfFixVecs<DimVec<double>>& coords);

    /** \brief
     * Implementation of \ref Operator::getElementMatrix(). Repalces the
     * determinant by the surface determinant and deligates the call to
     * the base class function.
     */
    virtual void getElementMatrix(ElInfo const* elInfo,
                                  ElementMatrix& userMat,
                                  double factor = 1.0) override;

    /** \brief
     * Implementation of \ref Operator::getElementVector(). Repalces the
     * determinant by the surface determinant and deligates the call to
     * the base class function.
     */
    virtual void getElementVector(ElInfo const* elInfo,
                                  DenseVector<double>& userVec,
                                  double factor = 1.0) override;

  protected:
    VectorOfFixVecs<DimVec<double>> coords_;

    /// Surface quadratures
    SurfaceQuadrature* quad2;
    SurfaceQuadrature* quad1GrdPsi;
    SurfaceQuadrature* quad1GrdPhi;
    SurfaceQuadrature* quad0;
  };

} // end namespace
