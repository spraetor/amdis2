/** \file is_symmetric.hpp */

#pragma once

namespace AMDiS
{
  template <class Matrix>
  inline bool is_symmetric(const Matrix& matrix, double tol = 1.e-5)
  {
    using mtl::tag::major;
    using mtl::tag::nz;
    using mtl::begin;
    using mtl::end;
    namespace traits= mtl::traits;

    typename traits::row<Matrix>::type                                 row(matrix);
    typename traits::col<Matrix>::type                                 col(matrix);
    typename traits::const_value<Matrix>::type                         value(matrix);

    typedef typename traits::range_generator<major, Matrix>::type      cursor_type;
    typedef typename traits::range_generator<nz, cursor_type>::type    icursor_type;

    for (cursor_type cursor = begin<major>(matrix), cend = end<major>(matrix); cursor != cend; ++cursor)
      for (icursor_type icursor = begin<nz>(cursor), icend = end<nz>(cursor); icursor != icend; ++icursor)
        // Compare each non-zero entry with its transposed
        if (abs(value(*icursor) - matrix[col(*icursor)][row(*icursor)]) > tol)
          return false;
    return true;
  }

} // end namspace AMDiS
