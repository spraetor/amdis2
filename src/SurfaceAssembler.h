/** \file SurfaceAssembler.h */

#pragma once

#include "AMDiS_fwd.h"
#include "Assembler.h"
#include "FixVec.h"

namespace AMDiS
{
  /**
   * \ingroup Integration
   *
   * \brief
   */
  class SurfaceAssembler : public Assembler
  {
  public:
    /// Creates a SurfaceAssembler conforming to operate for the given \ref coords.
    SurfaceAssembler(Operator* operat,
                     FiniteElemSpace const* rowFeSpace,
                     FiniteElemSpace const* colFeSpace,
                     VectorOfFixVecs<DimVec<double>>& coords);

    /// Destructor
    ~SurfaceAssembler();

    /// Adapt surface quadratures to \ref coords.
    void adaptSurfaceAssembler(VectorOfFixVecs<DimVec<double>>& coords);

    ///
    bool initElementVector(ElInfo const* elInfo);

  protected:
    VectorOfFixVecs<DimVec<double>> coords_;
  };

} // end namespace AMDiS
