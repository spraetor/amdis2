/** \file resize.hpp */

#pragma once

#include "AMDiS_fwd.h"
#include <traits/traits_fwd.hpp>

namespace AMDiS
{
  namespace traits
  {
    template <typename T>
    struct resize<T, typename enable_if<is_scalar<T>>::type>
    {
      void operator()(T& v, size_t r) const
      {
        TEST_EXIT_DBG(r == 1)
        ("Skalars can not be resized!\n");
      }
    };

    // == vectors ===

    template <typename T>
    struct resize<T, typename enable_if<and_<is_vector<T>, is_mtl<T>>>::type>
    {
      void operator()(T& v, size_t r) const
      {
        v.change_dim(r);
      }
    };


    /// change_dim implementation for STL vectors
    template <typename Value>
    struct resize<std::vector<Value>>
    {
      void operator()(std::vector<Value>& v, size_t r)
      {
        v.resize(r);
      }
    };


    // === matrices ===

    template <typename T>
    struct resize<T, typename enable_if<and_<is_matrix<T>, is_mtl<T>>>::type>
    {
      void operator()(T& m, size_t r, size_t c)
      {
        return m.change_dim(r,c);
      }
    };

  } // end namespace traits


  /// resize function for vectors
  template <typename Collection>
  typename boost::enable_if<typename traits::is_vector<Collection>::type>::type
  resize(Collection& c, size_t rows)
  {
    traits::resize<Collection>()(c, rows);
  }

  /// resize function for matrices
  template <typename Collection>
  typename boost::enable_if<typename traits::is_matrix<Collection>::type>::type
  resize(Collection& c, size_t rows, size_t cols)
  {
    traits::resize<Collection>()(c, rows, cols);
  }

} // end namespace AMDiS

#endif // AMDIS_TYPE_TRAITS_RESIZE_HPP
