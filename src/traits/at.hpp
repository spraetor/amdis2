/** \file at.hpp */

#pragma once

#include <traits/basic.hpp>
#include <traits/traits_fwd.hpp>
#include <traits/scalar_types.hpp>
#include <traits/concepts_base.hpp>

namespace AMDiS
{
  namespace traits
  {
    // == scalars ===
    template <class T>
    struct at<T, Requires_t<concepts::Arithmetic<T>> >
    {
      using value_type = T;
      value_type& operator()(T& v, size_t r)
      {
        TEST_EXIT_DBG(r == 0)("Scalars have only 1 component!\n");
        return v;
      }
      value_type const& operator()(T const& v, size_t r) const
      {
        TEST_EXIT_DBG(r == 0)("Scalars have only 1 component!\n");
        return v;
      }
      value_type& operator()(T& v, size_t r, size_t c)
      {
        TEST_EXIT_DBG(r == 0 && c == 0)("Scalars have only 1 row/column!\n");
        return v;
      }
      value_type const& operator()(T const& v, size_t r, size_t c) const
      {
        TEST_EXIT_DBG(r == 0 && c == 0)("Scalars have only 1 row/column!\n");
        return v;
      }
    };


    // == vectors ===
    template <class T>
    struct at<T, Requires_t<concepts::Vector<T>> >
    {
      using value_type = Value_t<T>;
      value_type& operator()(T& v, size_t r)
      {
        return v[r];
      }
      value_type const& operator()(T const& v, size_t r) const
      {
        return v[r];
      }
    };


    // === matrices ===
    template <typename T>
    struct at<T, Requires_t<concepts::Matrix<T>> >
    {
      using value_type = Value_t<T>;
      value_type& operator()(T& v, size_t r, size_t c)
      {
        return v(r,c);
      }
      value_type const& operator()(T const& v, size_t r, size_t c) const
      {
        return v(r,c);
      }
    };

  } // end namespace traits


  /// access components of vectors
  template <class Collection>
  auto at(Collection&& c, size_t rows) RETURNS
  (
    traits::at<Collection>()(std::forward<Collection>(c), rows)
  )

  /// access components of matrices
  template <class Collection>
  auto at(Collection&& c, size_t rows, size_t cols) RETURNS
  (
    traits::at<Collection>()(std::forward<Collection>(c), rows, cols)
  )

} // end namespace AMDiS
