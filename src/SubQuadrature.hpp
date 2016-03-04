#pragma once

#include "AMDiS_fwd.hpp"
#include "Quadrature.hpp"

namespace AMDiS
{
  class SubQuadrature : public Quadrature
  {
  public:
    SubQuadrature(Quadrature* quad, int dim_);

    void scaleQuadrature(VectorOfFixVecs<DimVec<double>>& coords);

    int getSubDim() const
    {
      return subDim_;
    }

  protected:
    Quadrature* quad_;

    int subDim_;
  };

} // end namespace AMDiS
