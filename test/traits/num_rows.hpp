/** \file num_rows.hpp */

#pragma once

// AMDiS headers
#include "traits/traits_fwd.hpp"
#include "traits/scalar_types.hpp"

namespace AMDiS
{
  template <class T, class = Requires_t<concepts::Arithmetic<T>>>
  inline size_t num_rows(T const& /*v*/)
  {
    return 1;
  }

} // end namespace AMDiS
