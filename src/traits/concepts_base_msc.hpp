/** \file concepts_base.hpp */

#pragma once

// std c++ headers
#include <utility>
#include <type_traits>

// AMDiS headers
#include "traits/basic.hpp"
#include "traits/meta_basic.hpp"

// macro to generate concept-checks
#define HAS_MEMBER_GENERATE(name) \
  template <class, class, class = void> struct has_member_ ## name : false_ {}; \
  template <class F, class Return, class Arg0>                             \
  struct has_member_ ## name <F, Return(Arg0),                             \
    Requires_t< std::is_convertible<                                       \
      decltype(std::declval<F>().name (std::declval<Arg0>())),             \
      Return >> >                                                          \
  : true_ {};                                                              \
  template <class F, class Return, class Arg0, class Arg1>                 \
  struct has_member_ ## name <F, Return(Arg0, Arg1),                       \
    Requires_t< std::is_convertible<                                       \
      decltype(std::declval<F>().name (std::declval<Arg0>(), std::declval<Arg1>())),  \
      Return >> >                                                          \
  : true_ {};                                                              \
  template <class F, class Return, class Arg0, class Arg1, class Arg2>     \
  struct has_member_ ## name <F, Return(Arg0, Arg1, Arg2),                 \
    Requires_t< std::is_convertible<                                       \
      decltype(std::declval<F>().name (std::declval<Arg0>(), std::declval<Arg1>(), std::declval<Arg2>())), \
      Return >> >                                                          \
  : true_ {};

#define HAS_MEMBER(name) has_member_ ## name

namespace AMDiS
{
  namespace traits
  {
    // Some type traits to test for type-attributes of a class
    // -------------------------------------------------------------------------

    namespace detail
    {
      template <class T, class = typename T::value_type >
      true_ HasValueTypeImpl(int);

      template <class T> false_ HasValueTypeImpl(...);


      template <class T, class = typename T::size_type >
      true_ HasSizeTypeImpl(int);

      template <class T> false_ HasSizeTypeImpl(...);


      template <class T, class = typename T::result_type >
      true_ HasResultTypeImpl(int);

      template <class T> false_ HasResultTypeImpl(...);

    } // end namespace detail


    template <class T>
    using HasValueType = decltype(detail::HasValueTypeImpl<T>(int{}));

    template <class T>
    using HasSizeType = decltype(detail::HasSizeTypeImpl<T>(int{}));

    template <class T>
    using HasResultType = decltype(detail::HasResultTypeImpl<T>(int{}));


    // Type traits to test whether a class is a functor, i.e. has a operator()
    // -------------------------------------------------------------------------

    // source: http://stackoverflow.com/questions/25603240/checking-callable-template-parameter-types
    template <class, class, class = void>
    struct IsFunctor : false_ {};

    template <class F, class Return, class Arg0>
    struct IsFunctor<F, Return(Arg0),
        Requires_t< std::is_same< decltype(std::declval<F>()(std::declval<Arg0>())),
                                  Return >> >
      : true_ {};

    template <class F, class Return, class Arg0, class Arg1>
    struct IsFunctor<F, Return(Arg0, Arg1),
        Requires_t< std::is_same< decltype(std::declval<F>()(std::declval<Arg0>(), std::declval<Arg1>())),
                                  Return >> >
      : true_ {};

    template <class F, class Return, class Arg0, class Arg1, class Arg2>
    struct IsFunctor<F, Return(Arg0, Arg1, Arg2),
        Requires_t< std::is_same< decltype(std::declval<F>()(std::declval<Arg0>(), std::declval<Arg1>(), std::declval<Arg2>())),
                                  Return >> >
      : true_ {};



    template <class, class, class = void>
    struct IsFunctorWeak : false_ {};

    template <class F, class Return, class Arg0>
    struct IsFunctorWeak<F, Return(Arg0),
        Requires_t< std::is_convertible< decltype(std::declval<F>()(std::declval<Arg0>())),
                                         Return >> >
      : true_ {};

    template <class F, class Return, class Arg0, class Arg1>
    struct IsFunctorWeak<F, Return(Arg0, Arg1),
        Requires_t< std::is_convertible< decltype(std::declval<F>()(std::declval<Arg0>(), std::declval<Arg1>())),
                                         Return >> >
      : true_ {};

    template <class F, class Return, class Arg0, class Arg1, class Arg2>
    struct IsFunctorWeak<F, Return(Arg0, Arg1, Arg2),
        Requires_t< std::is_convertible< decltype(std::declval<F>()(std::declval<Arg0>(), std::declval<Arg1>(), std::declval<Arg2>())),
                                         Return >> >
      : true_ {};

  } // end namespace traits


  /// Namespace that implements basic concept checks
  namespace concepts
  {
    /// @defgroup concepts Concepts definitions, using decltype and other c++11 features.

    // -------------------------------------------------------------------------
#if 0
    template <class, class, class = void>
    struct check_vector : false_ {};

    template <class F, class Return, class... Args>
    struct check_vector<F, Return(Args...),
        Requires_t< std::is_same< decltype(std::declval<F>()[std::declval<Args>()...]),
                                  Return >> >
      : true_ {};

    template <class, class, class = void>
    struct check_vector_weak : false_ {};

    template <class F, class Return, class... Args>
    struct check_vector_weak<F, Return(Args...),
        Requires_t< std::is_convertible< decltype(std::declval<F>()[std::declval<Args>()...]),
                                         Return >> >
      : true_ {};
#endif


    namespace detail
    {
      template <class T, class = decltype(&T::operator[])>
      true_ VectorImpl(int);

      template <class T> false_ VectorImpl(...);

    } // end namespace detail

    template <class T>
    using Vector = decltype(detail::VectorImpl<T>(int{}));

    template <class T>
    using Matrix = traits::IsFunctorWeak<T, Value_t<T>(size_t, size_t)>;

  } // end namespace concepts


  /// Namespace that implements basic concepts-checks using enable_if and
  /// type-traits from namespace \ref concepts
  namespace requires {}

} // end namespace AMDiS
