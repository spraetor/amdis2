#pragma once

#include "AMDiS_fwd.hpp"

namespace AMDiS
{
  class SubAssembler;

  struct LazyOperatorTermBase
  {
    template <class List>
    void insertFeSpaces(List& /*feSpaces*/) const {}

    constexpr int getDegree() const
    {
      return 0;
    }

    void initElement(ElInfo const*          /*elInfo*/,
                     SubAssembler*          /*subAssembler*/,
                     Quadrature*            /*quad*/,
                     BasisFunction const*   /*basisFct*/ = NULL) {}
  };
  
} // end namespace AMDiS
