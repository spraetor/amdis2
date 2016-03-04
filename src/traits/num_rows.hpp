/** \file num_rows.hpp */

#pragma once

// AMDiS headers
#include "AMDiS_fwd.hpp"
#include "traits/traits_fwd.hpp"
#include "traits/scalar_types.hpp"

namespace AMDiS
{
  template <class T, class = Requires_t<concepts::Arithmetic<T>>>
  inline size_t num_rows(T const& /*v*/)
  {
    return 1;
  }

  /// size implementation for AMDiS::VectorOfFixVecs
  template <class FixVecType>
  inline size_t num_rows(VectorOfFixVecs<FixVecType> const& v)
  {
    return v.getSize();
  }

  /// AMDiS::MatrixOfFixVecs
  template <class FixVecType>
  inline size_t num_rows(MatrixOfFixVecs<FixVecType> const& v)
  {
    return v.getNumberOfRows();
  }

} // end namespace AMDiS
