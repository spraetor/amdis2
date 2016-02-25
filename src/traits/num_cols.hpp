/** \file num_cols.hpp */

#pragma once

// AMDiS headers
#include "AMDiS_fwd.hpp"
#include "traits/traits_fwd.hpp"
#include "traits/scalar_types.hpp"

namespace AMDiS
{
  template <class T, class = Requires_t<concepts::Arithmetic<T>>>
  inline size_t num_cols(T const& /*v*/)
  {
    return 1;
  }

  /// size implementation for AMDiS::VectorOfFixVecs
  template <class FixVecType>
  inline size_t num_cols(VectorOfFixVecs<FixVecType> const& v)
  {
    return v.getSize() == 0 ? 0 : v[0].getSize();
  }

  /// AMDiS::MatrixOfFixVecs
  template <class FixVecType>
  inline size_t num_cols(MatrixOfFixVecs<FixVecType> const& v)
  {
    return v.getNumberOfColumns();
  }

} // end namespace AMDiS
