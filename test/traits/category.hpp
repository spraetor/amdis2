/** \file category.hpp */

#pragma once

// std c++ headers
#include <type_traits>

#include "traits/traits_fwd.hpp"
#include "traits/tag.hpp"

namespace AMDiS
{
  namespace traits
  {
    // scalars
    template <class T>
    struct category<T, Requires_t<concepts::Arithmetic<T>> >
    {
      typedef tag::scalar tag;
      typedef T           value_type;
      //       typedef size_t      size_type;
    };

    template <class T>
    struct category<WorldVector<WorldVector<T>>>
    {
      typedef tag::matrix  tag;
      typedef T            value_type;
      typedef int          size_type;
    };

    template <class T>
    struct category<const T> : category<T> {};

    // -------------------------------------------------------------------------

    template <class T, class S, class = void>
    struct IsConvertible : std::is_convertible<T, S> {};

    template <class T, class S>
    struct IsConvertible<T, S, Requires_t<and_<is_vector<T>, is_vector<S>>>>
      : std::is_convertible<Value_t<T>, Value_t<S>> {};

    template <class T, class S>
    struct IsConvertible<T, S, Requires_t<and_<is_matrix<T>, is_matrix<S>>>>
      : std::is_convertible<Value_t<T>, Value_t<S>> {};

  } // end namespace traits

} // end namespace AMDiS
